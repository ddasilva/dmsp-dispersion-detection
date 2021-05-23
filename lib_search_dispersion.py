"""Searches DMSP dataset for dispersion events, as defined by the dispersion
event detection algorithm developed by da Silva in 2020-2021.
"""
from datetime import datetime, timedelta
import h5py
from intervaltree import IntervalTree
from numba import jit
import numpy as np
import pandas as pd
import progressbar
import pytz
from spacepy import pycdf

from lib_util import find_moving_average


INTERVAL_LENGTH = 30              # seconds
INTEGRAL_THRESHOLD = 0.8          # units of Log(eV)
MIN_POS_FRAC = .8                 # fraction

MIN_AVG_IFLUX_SHEATH = 1e4        # units of diff en flux 
MIN_AVG_EFLUX_SHEATH = 1e6        # units of diff en flux 
MIN_PEAK_EFLUX_SHEATH = 10**6.5   # units of diff en flux

MIN_MLAT = 55                     # degrees
MIN_ION_VALID_ENERGY = 50         # eV; workaround for noise at low energies
MAX_SHEATH_ENERGY = 3.1e3         # eV


MAX_EIC_ENERGY = 3.15e3           # eV
MIN_BZ_STRENGTH = 3.0             # nT, must be >=0 (sign will be assigned)
MIN_IFLUX_AT_EIC = 10**6          # units of diff en flux

OMNIWEB_FILL_VALUE = 9999         # fill value for msising omniweb data


def read_omniweb_files(omniweb_files, silent=False):
    """Read OMNIWeb files into a single dictionary.

    Args
      omniweb_files: string path to cdf files
    Returns
      dictionary mapping parameters to file
    """
    # Read OMNIWeb Data
    # ------------------------------------------------------------------------------------
    t_items = []
    Bx_items = []
    By_items = []
    Bz_items = []
    
    for omniweb_file in sorted(omniweb_files):
        # Open file
        if not silent:
            print(f'Loading {omniweb_file}')
        omniweb_cdf = pycdf.CDF(omniweb_file)

        # Read the data
        t_items.append(np.array([time.replace(tzinfo=pytz.utc)
                           for time in omniweb_cdf['Epoch'][:]]))

        Bx_items.append(omniweb_cdf['BX_GSE'][:])
        By_items.append(omniweb_cdf['BY_GSM'][:])
        Bz_items.append(omniweb_cdf['BZ_GSM'][:])

        # Close file
        omniweb_cdf.close()

    # Merge arrays list of items
    omniweb_fh = {}
    omniweb_fh['t'] = np.concatenate(t_items)
    omniweb_fh['Bx'] = np.concatenate(Bx_items)
    omniweb_fh['By'] = np.concatenate(By_items)
    omniweb_fh['Bz'] = np.concatenate(Bz_items)

    return omniweb_fh


def read_dmsp_flux_file(dmsp_flux_filename):
    """Read DMSP Flux hdf or cdf file and load variables into a dictionary.
    
    Args
      dmsp_flux_filename: string path to hdf or cdf file
    Returns
      dictionary mapping parameters to file
    """
    # Read DMSP Data in HDF5 Format
    # ------------------------------------------------------------------------------------
    if dmsp_flux_filename.endswith(('.hdf5', '.hdf', '.h5')):
        # Open file
        hdf = h5py.File(dmsp_flux_filename, 'r')
        
        # Populate file handle dictionary 
        dmsp_flux_fh = {}    
        dmsp_flux_fh['t'] = np.array(
            [datetime(1970, 1, 1, tzinfo=pytz.utc) + timedelta(seconds=i)
             for i in hdf['Data']['Array Layout']['timestamps'][:]]
        )
        dmsp_flux_fh['ch_energy'] = hdf['Data']['Array Layout']['ch_energy'][:]
        dmsp_flux_fh['mlat'] = hdf['Data']['Array Layout']['1D Parameters']['mlat'][:]
        dmsp_flux_fh['mlt'] = hdf['Data']['Array Layout']['1D Parameters']['mlt'][:]
        dmsp_flux_fh['el_d_ener'] = hdf['Data']['Array Layout']['2D Parameters']['el_d_ener'][:]
        dmsp_flux_fh['ion_d_ener'] = hdf['Data']['Array Layout']['2D Parameters']['ion_d_ener'][:]
        
        # Close file
        hdf.close()
        
    # Read DMSP Data in CDF Format
    # ------------------------------------------------------------------------------------
    elif dmsp_flux_filename.endswith('cdf'):

        # Open file
        cdf = pycdf.CDF(dmsp_flux_filename)

        dmsp_flux_fh = {}    
        dmsp_flux_fh['t'] = np.array([t.replace(tzinfo=pytz.utc) for t in cdf['Epoch'][:]])
        dmsp_flux_fh['ch_energy'] = cdf['CHANNEL_ENERGIES'][:][::-1]
        dmsp_flux_fh['mlat'] = cdf['SC_AACGM_LAT'][:]
        dmsp_flux_fh['mlt'] = cdf['SC_AACGM_LTIME'][:]
        dmsp_flux_fh['el_d_ener'] = cdf['ELE_DIFF_ENERGY_FLUX'][:].T
        dmsp_flux_fh['ion_d_ener'] = cdf['ION_DIFF_ENERGY_FLUX'][:].T

        # Close file
        cdf.close()
        
    # Crash for any other format
    # ------------------------------------------------------------------------------------
    else:
        raise RuntimeError(f'Invalid file format for {dmsp_flux_filename}')
        
    # Compute (simple) derived variables
    # ------------------------------------------------------------------------------------
    # Ion and Electron Average Flux
    ch_i = dmsp_flux_fh['ch_energy'].searchsorted(MIN_ION_VALID_ENERGY)
    ch_j = dmsp_flux_fh['ch_energy'].searchsorted(MAX_SHEATH_ENERGY)
    
    dmsp_flux_fh['iflux_avg_sheath'] = np.mean(
        dmsp_flux_fh['ion_d_ener'][ch_i:ch_j, :], axis=0
    )
    dmsp_flux_fh['eflux_avg_sheath'] = np.mean(
        dmsp_flux_fh['el_d_ener'][:ch_j, :], axis=0
    )
    # Ion and Electron Peak Flux in Sheath
    dmsp_flux_fh['iflux_peak_sheath'] = np.max(
        dmsp_flux_fh['ion_d_ener'][ch_i:ch_j, :], axis=0
    )
    dmsp_flux_fh['eflux_peak_sheath'] = np.max(
        dmsp_flux_fh['el_d_ener'][:ch_j, :], axis=0
    )

    
    return dmsp_flux_fh


def estimate_Eic(dmsp_flux_fh, i, j, frac=.1):    
    """Calculates the Eic parameter-- energy level with 10% of peak flux.

    Starting at the energy level of peak flux, works downwards and selects the
    energy level at approximately 10% of the peak flux.

    See Also: Lockwood 1992, Journal of Geophysical Research
    
    Args
      dmsp_flux_fh: file handle (as returned by read_dmsp_flux_file)
      i: start index to limit search
      j: end index to limit search
      frac: algorithm parameter; fraction of peak flux to use for energy search
    Returns
      Eic: floating point numpy array
    """
    ch_bot = dmsp_flux_fh['ch_energy'].searchsorted(MIN_ION_VALID_ENERGY)
    flux_max = dmsp_flux_fh['ion_d_ener'][ch_bot:, i:j].max(axis=0)
    flux_max_ind = dmsp_flux_fh['ion_d_ener'][ch_bot:, i:j].argmax(axis=0) + ch_bot # over time

    # holds channels with flux above frac * max flux
    threshold_match_mask = (dmsp_flux_fh['ion_d_ener'][:, i:j] > frac * flux_max)
    fill_mask = np.zeros(dmsp_flux_fh['t'][i:j].shape, dtype=bool)
    
    for jj, kk in enumerate(flux_max_ind):
        threshold_match_mask[kk:, jj] = 0

        # if all points under max flux energy under threshold
        if np.all(threshold_match_mask[:kk, jj] == False):
            fill_mask[jj] = True

    # select last true value under max flux energy
    ind = len(dmsp_flux_fh['ch_energy']) - 1 - \
        threshold_match_mask[::-1, :].argmax(axis=0)

    # if no points > frac * max flux under max flux energy, then just use
    # max flux energy
    Eic = dmsp_flux_fh['ch_energy'][ind].copy()
    Eic[fill_mask] = dmsp_flux_fh['ch_energy'][flux_max_ind[fill_mask]]     
    
    return Eic


def estimate_log_Eic_smooth_derivative(
        dmsp_flux_fh, eic_window_size=11, eic_deriv_window_size=5
):
    """Calculate the smoothed derivative of the smoothed Log10(Eic) parameter.
    
    Args
      dmsp_flux_fh: file handle returned by read_files()
    Returns
      dLogEicdt_smooth: smoothed derivative corresponding to times found in dmsp_flux_fh['t']
      Eic_smooth: smoothed Eic cooresponding to times found in dmsp_flux_fh['t']
    """
    # Calcualte the Eic parameter. Throughout this function, the Eic is in
    # log-space.
    Eic = np.log10(estimate_Eic(dmsp_flux_fh, i=0, j=dmsp_flux_fh['t'].size))

    en_inds = np.log10(dmsp_flux_fh['ch_energy']).searchsorted(Eic)
    en_inds = [min(i, 18) for i in en_inds]
    flux_at_Eic = np.array(
        [dmsp_flux_fh['ion_d_ener'][en_ind, j] for (j, en_ind)
         in enumerate(en_inds)]
    )
    mask = (flux_at_Eic > MIN_IFLUX_AT_EIC)

    #Eic_smooth = Eic.copy()
    #Eic_smooth = smooth_Eic(Eic, mask, 1)
    Eic_smooth = smooth_Eic(Eic, mask, eic_window_size)
    
    # Find the smoothed derivative of the smoothed log10(Eic) function. For sake
    # of simplicity, the derivative is estimated with a forward difference.
    dLogEic = Eic_smooth.copy()
    dLogEic[:-1] = np.diff(Eic_smooth)
    dLogEic[-1] = dLogEic[-2]                       # copy last value to retain shape
    
    dt = [delta.total_seconds() for delta in np.diff(dmsp_flux_fh['t'])]
    dt.append(dt[-1])                         # copy last value to retain shape)
    dt = np.array(dt)
    
    dLogEicdt_smooth = dLogEic / dt

    #dLogEicdt_smooth = find_moving_average(dLogEic/dt, eic_deriv_window_size)
    #dLogEicdt_smooth[np.isnan(dLogEicdt_smooth)] = 0
        
    return dLogEicdt_smooth, Eic_smooth


def walk_and_integrate(dmsp_flux_fh, omniweb_fh, dLogEicdt_smooth, Eic_smooth,
                       interval_length, reverse_effect=False,
                       return_integrand=False):
    """Walk through windows in the file and test for matching intervals with
    integration of the metric function.
    
    Args
      dmsp_flux_fh: file handle returned by read_files()
      omniweb_fh: file handle returned by read_files()
      dLogEicdt_smooth: smoothed derivative corresponding to times found in
        dmsp_flux_fh['t']
      Eic_smooth: smoothed Eic cooresponding to times found in dmsp_flux_fh['t']
      interval_length: length of interval
      reverse_effect: Search for effects in the opposite direction with a magnetic
        field set to the opposite of the coded threshold.
      return_integrand: return array of integrand values
    """
    # Make interval length into timedelta
    interval_length = timedelta(seconds=interval_length)
    bar = progressbar.ProgressBar()
    matching_intervals = IntervalTree()
    integrand_save = np.zeros(dmsp_flux_fh['t'].size)
    integral_save = np.nan * np.zeros(dmsp_flux_fh['t'].size)
    
    # Walk through timeseries 
    for start_time_idx, start_time in bar(list(enumerate(dmsp_flux_fh['t']))):
        # First, check that the minutely Bz measurement from OMNIWeb associated
        # with the midpoint of this interval is less than the threshold.
        B_test_time = start_time + timedelta(seconds=INTERVAL_LENGTH/2)
        B_i = omniweb_fh['t'].searchsorted(B_test_time)
        if B_i == omniweb_fh['Bz'].size:
            continue
        Bx = omniweb_fh['Bx'][B_i]
        By = omniweb_fh['By'][B_i]
        Bz = omniweb_fh['Bz'][B_i]
        
        if Bz > OMNIWEB_FILL_VALUE:
            continue
        elif (not reverse_effect) and Bz > -MIN_BZ_STRENGTH:
            continue
        elif reverse_effect and Bz < MIN_BZ_STRENGTH:
            continue

        # Second, check that the MLT associated with the interval is in the day-
        # side region.
        mlt_test_time = start_time + timedelta(seconds=INTERVAL_LENGTH/2)
        mlt_i = dmsp_flux_fh['t'].searchsorted(mlt_test_time)
        if mlt_i == dmsp_flux_fh['mlt'].size:
            continue
        mlt = dmsp_flux_fh['mlt'][mlt_i]

        if not (6 < mlt < 18):
            continue
        
        # Determine end time of interval. If less than `interval_length` from
        # the end of the file, the interval may be less than `interval_length`.
        end_time_idx = dmsp_flux_fh['t'].searchsorted(start_time + interval_length)

        if end_time_idx > dmsp_flux_fh['t'].size:
            end_time_idx = dmsp_flux_fh['t'].size

        end_time = dmsp_flux_fh['t'][end_time_idx - 1]

        # Only check the interval if the magnetic latitude range
        # |mlat| > 60 deg.
        mlat = dmsp_flux_fh['mlat'][start_time_idx:end_time_idx]

        if not np.all(np.abs(mlat) > MIN_MLAT):
            continue
        
        # Setup integrand for integration. This contains multiple
        # multiplicative terms to control different aspects of the value.
        mlat_direction = -np.sign(np.diff(np.abs(mlat)))

        if reverse_effect:
            mlat_direction *= -1
                
        iflux_avg_sheath_mask = (
            dmsp_flux_fh['iflux_avg_sheath'][start_time_idx:end_time_idx]
            > MIN_AVG_IFLUX_SHEATH).astype(int)
        eflux_avg_sheath_mask = (
            dmsp_flux_fh['eflux_avg_sheath'][start_time_idx:end_time_idx]
            > MIN_AVG_EFLUX_SHEATH).astype(int)        
        eflux_peak_sheath_mask = (
            dmsp_flux_fh['eflux_peak_sheath'][start_time_idx:end_time_idx]
            > MIN_PEAK_EFLUX_SHEATH).astype(int)
        Eic_in_range = (
            Eic_smooth[start_time_idx:end_time_idx]
            < np.log10(MAX_EIC_ENERGY)).astype(int)

        en_inds = np.log10(dmsp_flux_fh['ch_energy']).searchsorted(Eic_smooth[start_time_idx:end_time_idx])
        en_inds = [min(i, 18) for i in en_inds]
        flux_at_Eic = dmsp_flux_fh['ion_d_ener'][en_inds, np.arange(start_time_idx, end_time_idx)]
        flux_at_Eic[np.isnan(Eic_smooth[start_time_idx:end_time_idx])] = np.nan        
        with np.errstate(invalid='ignore'):
            flux_at_Eic_mask = (
                np.isfinite(flux_at_Eic) & (flux_at_Eic > MIN_IFLUX_AT_EIC)
            ).astype(int)
        
        integrand = (
            mlat_direction *
            iflux_avg_sheath_mask[:-1] *
            eflux_avg_sheath_mask[:-1] *
            eflux_peak_sheath_mask[:-1] *
            Eic_in_range[:-1] *
            flux_at_Eic_mask[:-1] *
            dLogEicdt_smooth[start_time_idx:end_time_idx-1]
        )

        t = dmsp_flux_fh['t'][start_time_idx:end_time_idx]
        dt = [delta.total_seconds() for delta in np.diff(t)]

        if integrand.size > 0:
            integrand_save[start_time_idx] = integrand[0]
        else:
            continue
        
        # Compute integral of the integrand over increments of dt using
        # rectangular integraion
        t = dmsp_flux_fh['t'][start_time_idx:end_time_idx]
        dt = np.array([delta.total_seconds() for delta in np.diff(t)])

        integral = np.sum(integrand * dt)        
        pos_frac = integral / np.sum(np.abs(integrand)*dt)
        integral_save[start_time_idx] = integral
                
        # The test/accept condition on the integral value and the fraction of
        # values above zero.
        if integral > INTEGRAL_THRESHOLD and pos_frac > MIN_POS_FRAC:
            matching_intervals[start_time:end_time] = {
                'integral': integral,
                'Bx': Bx,
                'By': By,
                'Bz': Bz,
            }

    # Merge overlapping intervals into common intervals. Retain the 
    # metadata attached to each.
    def reducer(current, new_data):
        for key in new_data:
            if key not in current:
                current[key] = set()
            current[key].add(new_data[key])

        return current

    matching_intervals.merge_overlaps(
        data_reducer=reducer, data_initializer={}
    )

    # Convert to pandas dataframe for easy output to terminal of matching
    # intervals and associated metadata.
    df_match_rows = []
    
    for interval in sorted(matching_intervals):
        df_match_rows.append([
            interval.begin, interval.end,
            np.mean(list(interval.data['Bx'])),
            np.mean(list(interval.data['By'])),
            np.mean(list(interval.data['Bz'])),
            np.mean(list(interval.data['integral'])),
        ])

    df_match = pd.DataFrame(df_match_rows, columns=[
        'start_time', 'end_time',
        'Bx_mean', 'By_mean', 'Bz_mean',
        'integral_mean'
    ])

    if return_integrand:
        return df_match, integrand_save, integral_save, pos_frac
    else:
        return df_match


@jit(nopython=True)
def smooth_Eic(Eic, keep_mask, window_size):
    """Smooth Eic with a mask of points to include in moving average.
    
    Arguments
      Eic: floating point numpy array
      keep_mask: points to include in moving average
      window_size: integer, must be odd
    Returns
      Smoothed Eic array
    """
    assert window_size % 2 == 1, 'Window size must be odd'

    Eic_smooth = Eic.copy()
    Eic_smooth[:] = np.nan
    
    for i in range(Eic.size):
        if not keep_mask[i]:
            Eic_smooth[i] = np.nan
        else:
            total = 0.0
            count = 0
            for di in range(-window_size//2, window_size//2 + 1):
                if i + di > 0 and i + di < Eic.size and keep_mask[i + di]:
                    total += Eic[i + di]
                    count += 1
            if count > 0:  # else left as nan
                Eic_smooth[i] = total / count

    return Eic_smooth

