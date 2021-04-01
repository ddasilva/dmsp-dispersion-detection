#!/usr/bin/env python
"""Searches DMSP dataset for dispersion events, as defined by the dispersion
event detection algorithm developed by da Silva in 2020-2021.
"""
from datetime import datetime, timedelta
import h5py
from intervaltree import IntervalTree
import numpy as np
import pandas as pd
import progressbar
import pytz
from spacepy import pycdf


INTERVAL_LENGTH = 60              # seconds
INTEGRAL_THRESHOLD = 0.9

MIN_ENERGY_EL_PDENS = 10**1.5     # Electron partial density lower energy, eV
MAX_ENERGY_EL_PDENS = 10**3.5     # Electron partial density upper energy, eV
MIN_ION_DENS = 10**7.25           # cm^-3
MIN_EL_PDENS = 10**5              # cm^-3
MIN_MLAT = 55                     # degrees
MIN_ENERGY_ANALYZED = 50          # eV
MAX_ENERGY_ANALYZED = 10**3.5     # eV
MIN_BZ_STRENGTH = 3.0             # nT, must be >=0 (sign will be assigned)
MIN_FLUX_AT_EIC = 10**6.0         # spectrogram units
MIN_POS_INTEGRAL_FRAC = 0.8       # fraction
OMNIWEB_FILL_VALUE = 9999         # Bz fill value for msising omniweb data
    

def read_omniweb_files(omniweb_files, silent=False):
    """Read OMNIWeb files into a single dictionary.

    Args
      filenames: string path to cdf files
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


def read_dmsp_file(dmsp_filename):
    """Read DMSP hdf or cdf file and load variables into a dictionary.
    
    Args
      filename: string path to hdf or cdf file
    Returns
      dictionary mapping parameters to file
    """
    # Read DMSP Data in HDF5 Format
    # ------------------------------------------------------------------------------------
    if dmsp_filename.endswith(('.hdf5', '.hdf', '.h5')):
        # Open file
        hdf = h5py.File(dmsp_filename, 'r')
        
        # Populate file handle dictionary 
        dmsp_fh = {}    
        dmsp_fh['t'] = np.array(
            [datetime(1970, 1, 1, tzinfo=pytz.utc) + timedelta(seconds=i)
             for i in hdf['Data']['Array Layout']['timestamps'][:]]
        )
        dmsp_fh['ch_energy'] = hdf['Data']['Array Layout']['ch_energy'][:]
        dmsp_fh['mlat'] = hdf['Data']['Array Layout']['1D Parameters']['mlat'][:]
        dmsp_fh['el_d_ener'] = hdf['Data']['Array Layout']['2D Parameters']['el_d_ener'][:]
        dmsp_fh['ion_d_ener'] = hdf['Data']['Array Layout']['2D Parameters']['ion_d_ener'][:]
        
        # Close file
        hdf.close()
        
    # Read DMSP Data in CDF Format
    # ------------------------------------------------------------------------------------
    elif dmsp_filename.endswith('cdf'):

        # Open file
        cdf = pycdf.CDF(dmsp_filename)
 
        dmsp_fh = {}    
        dmsp_fh['t'] = np.array([t.replace(tzinfo=pytz.utc) for t in cdf['Epoch'][:]])
        dmsp_fh['ch_energy'] = cdf['CHANNEL_ENERGIES'][:][::-1]
        dmsp_fh['mlat'] = cdf['SC_AACGM_LAT'][:]
        dmsp_fh['el_d_ener'] = cdf['ELE_DIFF_ENERGY_FLUX'][:].T
        dmsp_fh['ion_d_ener'] = cdf['ION_DIFF_ENERGY_FLUX'][:].T

        # Close file
        cdf.close()
        
    # Crash for any other format
    # ------------------------------------------------------------------------------------
    else:
        raise RuntimeError(f'Invalid file format for {dmsp_filename}')
        
    # Compute (simple) derived variables
    # ------------------------------------------------------------------------------------
    # Ion Density
    dmsp_fh['ion_dens'] = 4 * np.pi * np.trapz(
        dmsp_fh['ion_d_ener'].T / dmsp_fh['ch_energy'],
        dmsp_fh['ch_energy'], axis=1
    )
    
    # Electron partial density
    ch_i = dmsp_fh['ch_energy'].searchsorted(MIN_ENERGY_EL_PDENS)
    ch_j = dmsp_fh['ch_energy'].searchsorted(MAX_ENERGY_EL_PDENS)    
    dmsp_fh['el_pdens'] = 4 * np.pi * np.trapz(
        dmsp_fh['el_d_ener'][ch_i:ch_j].T / dmsp_fh['ch_energy'][ch_i:ch_j],
        dmsp_fh['ch_energy'][ch_i:ch_j], axis=1
    )
    
    return dmsp_fh


def estimate_Eic(dmsp_fh, i, j, frac=.1):    
    """Calculates the Eic parameter-- energy level with 10% of peak flux.

    Starting at the energy level of peak flux, works downwards and selects the
    energy level at approximately 10% of the peak flux.

    See Also: Lockwood 1992, Journal of Geophysical Research
    
    Args
      dmsp_fh: file handle (as returned by read_files)
      i: start index to limit search
      j: end index to limit search
      frac: algorithm parameter; fraction of peak flux to use for energy search
    Returns
      Eic: floating point numpy array
    """
    ch_bot = dmsp_fh['ch_energy'].searchsorted(MIN_ENERGY_ANALYZED)
    ch_top = dmsp_fh['ch_energy'].searchsorted(MAX_ENERGY_ANALYZED)
    flux_max = dmsp_fh['ion_d_ener'][ch_bot:ch_top, i:j].max(axis=0)
    flux_max_ind = dmsp_fh['ion_d_ener'][ch_bot:ch_top, i:j].argmax(axis=0) + ch_bot # over time

    # holds channels with flux above frac * max flux
    threshold_match_mask = (dmsp_fh['ion_d_ener'][:ch_top, i:j] > frac * flux_max)
    fill_mask = np.zeros(dmsp_fh['t'][i:j].shape, dtype=bool)
    
    for jj, kk in enumerate(flux_max_ind):
        threshold_match_mask[kk:, jj] = 0

        # if all points under max flux energy under threshold
        if np.all(threshold_match_mask[:kk, jj] == False):
            fill_mask[jj] = True

    # select last true value under max flux energy
    ind = ch_top - 1 - threshold_match_mask[::-1, :].argmax(axis=0)

    # if no points > frac * max flux under max flux energy, then just use
    # max flux energy
    Eic = dmsp_fh['ch_energy'][ind].copy()
    Eic[fill_mask] = dmsp_fh['ch_energy'][flux_max_ind[fill_mask]] 

    return Eic


def find_moving_average(a, window_size) :
    """Smooth array `a` using moving average with window size `n`.

    This code was adapted from an example found on stack overflow.

    Args
      a: numeric numpy array
      window_size: integer window size
    Returns
      smoothed array using moving average.
    """
    return np.convolve(a, np.ones((window_size,))/window_size, mode='same')


def estimate_log_Eic_smooth_derivative(
        dmsp_fh, eic_window_size=11, eic_deriv_window_size=5
):
    """Calculate the smoothed derivative of the smoothed Log10(Eic) parameter.
    
    Args
      dmsp_fh: file handle returned by read_files()
    Returns
      dEicdt_smooth: smoothed derivative corresponding to times found in dmsp_fh['t']
      Eic_smooth: smoothed Eic cooresponding to times found in dmsp_fh['t']
    """
    # Calcualte the Eic parameter. Throughout this function, the Eic is in
    # log-space.
    Eic = np.log10(estimate_Eic(dmsp_fh, i=0, j=dmsp_fh['t'].size))
    Eic_smooth = find_moving_average(Eic, eic_window_size)
    
    # Find the smoothed derivative of the smoothed log10(Eic) function. For sake
    # of simplicity, the derivative is estimated with a forward difference.
    dEic = Eic_smooth.copy()
    dEic[:-1] = np.diff(Eic_smooth)
    dEic[-1] = dEic[-2]                       # copy last value to retain shape
    
    dt = [delta.total_seconds() for delta in np.diff(dmsp_fh['t'])]
    dt.append(dt[-1])                         # copy last value to retain shape)
    
    dEicdt_smooth = find_moving_average(dEic/dt, eic_deriv_window_size)
    dEicdt_smooth[np.isnan(dEicdt_smooth)] = 0
        
    return dEicdt_smooth, Eic_smooth


def walk_and_integrate(dmsp_fh, omniweb_fh, dEicdt_smooth, Eic_smooth, interval_length,
                       reverse_effect=False, return_integrand=False):
    """Walk through windows in the file and test for matching intervals with
    integration of the metric function.
    
    Args
      dmsp_fh: file handle returned by read_files()
      omniweb_fh: file handle returned by read_files()
      dEicdt_smooth: smoothed derivative corresponding to times found in dmsp_fh['t']
      Eic_smooth: smoothed Eic cooresponding to times found in dmsp_fh['t']
      interval_length: length of interval
      reverse_effect: Search for effects in the opposite direction with a magnetic
        field set to the opposite of the coded threshold.
      return_integrand: return array of integrand values
    """
    # Make interval length into timedelta
    interval_length = timedelta(seconds=interval_length)
    bar = progressbar.ProgressBar()
    matching_intervals = IntervalTree()
    integrand_save = np.zeros(dmsp_fh['t'].size)
    integral_save = np.nan * np.zeros(dmsp_fh['t'].size)
    pos_integral_frac_save = np.nan * integrand_save.copy() 
    
    # Walk through timeseries 
    for start_time_idx, start_time in bar(list(enumerate(dmsp_fh['t']))):
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

        # Determine end time of interval. If less than `interval_length` from
        # the end of the file, the interval may be less than `interval_length`.
        end_time_idx = dmsp_fh['t'].searchsorted(start_time + interval_length)

        if end_time_idx > dmsp_fh['t'].size:
            end_time_idx = dmsp_fh['t'].size

        end_time = dmsp_fh['t'][end_time_idx - 1]

        # Only check the interval if the magnetic latitude range
        # |mlat| > 60 deg.
        mlat = dmsp_fh['mlat'][start_time_idx:end_time_idx]

        if not np.all(np.abs(mlat) > MIN_MLAT):
            continue
        
        # Setup integrand for integration. This contains multiple
        # multiplicative terms to control different aspects of the value.
        mlat_direction = -np.sign(np.diff(np.abs(mlat)))

        if reverse_effect:
            mlat_direction *= -1
        
        with np.errstate(divide='ignore'):
            ion_dens_log = np.log10(dmsp_fh['ion_dens'][start_time_idx:end_time_idx] + .01)
            el_pdens_log = np.log10(dmsp_fh['el_pdens'][start_time_idx:end_time_idx] + .01)
        
        ion_dens_mask = (ion_dens_log > np.log10(MIN_ION_DENS)).astype(int)
        el_pdens_mask = (el_pdens_log > np.log10(MIN_EL_PDENS)).astype(int)        
        
        en_inds = np.log10(dmsp_fh['ch_energy']).searchsorted(Eic_smooth[start_time_idx:end_time_idx])
        en_inds = [min(i, 18) for i in en_inds]
        flux_at_Eic = dmsp_fh['ion_d_ener'][en_inds, np.arange(start_time_idx, end_time_idx)]
        flux_at_Eic[np.isnan(Eic_smooth[start_time_idx:end_time_idx])] = np.nan

        with np.errstate(invalid='ignore'):
            flux_at_Eic_mask = (
                np.isfinite(flux_at_Eic) & (flux_at_Eic > MIN_FLUX_AT_EIC)
            ).astype(int)
        
        integrand = (
            mlat_direction *
            ion_dens_mask[:-1] *
            el_pdens_mask[:-1] *
            flux_at_Eic_mask[:-1] *
            dEicdt_smooth[start_time_idx:end_time_idx-1]
        )

        t = dmsp_fh['t'][start_time_idx:end_time_idx]
        dt = [delta.total_seconds() for delta in np.diff(t)]

        if integrand.size > 0:
            integrand_save[start_time_idx] = integrand[0]
        else:
            continue
        
        # Compute integral of the integrand over increments of dt using
        # rectangular integraion
        t = dmsp_fh['t'][start_time_idx:end_time_idx]
        dt = np.array([delta.total_seconds() for delta in np.diff(t)])

        integral = np.sum(integrand * dt)
                
        upper_integral = np.sum(integrand[integrand > 0] * dt[integrand > 0])
        abs_integral =  np.sum(np.abs(integrand) * dt)
        
        if abs_integral > 0:
            pos_integral_frac = upper_integral / abs_integral
        else:
            pos_integral_frac = 0
        
        pos_integral_frac_save[start_time_idx] = pos_integral_frac
        integral_save[start_time_idx] = integral
                
        # The test/accept condition on the integral value and the fraction of
        # values above zero.
        if (integral > INTEGRAL_THRESHOLD and pos_integral_frac > MIN_POS_INTEGRAL_FRAC):
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
        return df_match, integrand_save, integral_save, pos_integral_frac_save
    else:
        return df_match

