#!/usr/bin/env python
"""Searches DMSP dataset for dispersion events, as defined
by the dispersion event detection algorithm developed by da
Silve in 2020.
"""
import argparse
from datetime import datetime, timedelta
import h5py
from intervaltree import IntervalTree
import numpy as np
import os
import pandas as pd
import progressbar
import pytz
import sys


DEFAULT_INTERVAL_LENGTH = 60

DEFAULT_MAX_FLUX_THRESHOLD = 6.0  # eV

DEFAULT_INTEGRAL_THRESHOLD = 0.7  # Log eV

DEFAULT_FRAC_ABOVE_ZERO_THRESHOLD = .7


def main():
    """Main function of the program."""
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('dmsp_file')
    parser.add_argument('-I', '--interval-length', default=DEFAULT_INTERVAL_LENGTH,
                        type=float,
                        help='Interval length (seconds)')
    parser.add_argument('-F', '--max-flux-threshold', default=DEFAULT_MAX_FLUX_THRESHOLD,
                        type=float,
                        help='Max Flux Threshold')
    parser.add_argument('-S', '--integral-threshold', default=DEFAULT_INTEGRAL_THRESHOLD,
                        type=float,
                        help='Integral threshold')
    parser.add_argument('-f', '--frac-above-zero-threshold',
                        default=DEFAULT_FRAC_ABOVE_ZERO_THRESHOLD,
                        type=float,
                        help='Fraction above zero threshold')

    args = parser.parse_args()

    # Check file exists
    if not os.path.exists(args.dmsp_file):
        print(f'File {dmsp_file} does not exist', file=sys.stderr)
        sys.exit(1)

    # Calculate  smoothed dEic/dt in log-space
    fh = read_file(args.dmsp_file)

    dEicdt_smooth = estimate_log_Eic_smooth_derivative(fh)

    df_match = walk_and_integrate(
        fh, dEicdt_smooth, args.interval_length, args.max_flux_threshold,
        args.integral_threshold, args.frac_above_zero_threshold
    )

    # Output intervals to terminal
    print(df_match.to_string(index=False))
    

def read_file(filename):
    """Read DMSP hdf file and read variables into a dictionary.
    
    Args
      filename: string path to hdf file
    Returns
      dictionary mapping parameters to file
    """
    # Open file
    hdf = h5py.File(filename, 'r')

    # Populate file handle dictionary 
    fh = {}    
    fh['t'] = np.array(
        [datetime(1970, 1, 1, tzinfo=pytz.utc) + timedelta(seconds=i)
         for i in hdf['Data']['Array Layout']['timestamps'][:]]
    )
    fh['ch_energy'] = hdf['Data']['Array Layout']['ch_energy'][:]
    fh['mlat'] = hdf['Data']['Array Layout']['1D Parameters']['mlat'][:]
    fh['mlt'] = hdf['Data']['Array Layout']['1D Parameters']['mlt'][:]
    fh['ion_d_flux'] = hdf['Data']['Array Layout']['2D Parameters']['ion_d_flux'][:]
    fh['ion_d_ener'] = hdf['Data']['Array Layout']['2D Parameters']['ion_d_ener'][:]

    # Close file
    hdf.close()

    return fh


def estimate_Eic(fh, i, j, frac=.1):    
    """Calculates the Eic parameter-- energy level with 10% of peak flux.

    Starting at the energy level of peak flux, works downwards and selects the
    energy level at approximately 10% of the peak flux.

    See Also: Lockwood 1992, Journal of Geophysical Research
    
    Args
      fh: file handle (as returned by read_file)
      i: start index to limit search
      j: end index to limit search
      frac: algorithm parameter; fraction of peak flux to use for energy search
    Returns
      Eic: floating point numpy array
    """
    flux_max = fh['ion_d_ener'][:, i:j].max(axis=0)
    flux_max_ind = fh['ion_d_ener'][:, i:j].argmax(axis=0)

    # holds channels with flux above frac * max flux
    threshold_match_mask = (fh['ion_d_ener'][:, i:j] > frac * flux_max) 
    fill_mask = np.zeros(fh['t'][i:j].shape, dtype=bool)
    
    for jj, kk in enumerate(flux_max_ind):
        threshold_match_mask[kk:, jj] = 0

        # if all points under max flux energy under threshold
        if np.all(threshold_match_mask[:kk, jj] == False):
            fill_mask[jj] = True

    # select last true value under max flux energy
    n_channels = len(fh['ch_energy'])
    ind = n_channels - 1 - threshold_match_mask[::-1, :].argmax(axis=0)

    # if no points > frac * max flux under max flux energy, then just use
    # max flux energy
    Eic = fh['ch_energy'][ind].copy()
    Eic[fill_mask] = fh['ch_energy'][flux_max_ind[fill_mask]] 
                                                           
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


def estimate_log_Eic_smooth_derivative(fh, eic_window_size=11,
                                       eic_deriv_window_size=5,
                                       return_Eic=False):
    """Calculate the smoothed derivative of the smoothed Log10(Eic) parameter.
    
    Args
      fh: file handle returned by read_file()
    Returns
      dEicdt_smooth: smoothed derivative corresponding to times found in fh['t']
    """
    # Calcualte the Eic parameter. Throughout this function, the Eic is in
    # log-space.
    Eic = np.log10(estimate_Eic(fh, i=0, j=fh['t'].size))
    Eic_smooth = find_moving_average(Eic, eic_window_size)

    # Find the smoothed derivative of the smoothed log10(Eic) function. For sake
    # of simplicity, the derivative is estimated with a forward difference.
    dEic = Eic_smooth.copy()
    dEic[:-1] = np.diff(Eic_smooth)
    dEic[-1] = dEic[-2]                       # copy last value to retain shape
    
    dt = [delta.total_seconds() for delta in np.diff(fh['t'])]
    dt.append(dt[-1])                         # copy last value to retain shape)
    
    dEicdt_smooth = find_moving_average(dEic/dt, eic_deriv_window_size)

    if return_Eic:
        return dEicdt_smooth, Eic_smooth
    else:
        return dEicdt_smooth


def walk_and_integrate(fh, dEicdt_smooth, interval_length, max_flux_threshold,
                       integral_threshold, frac_above_zero_threshold,
                       return_integrand=False):
    """Walk through windows in the file and test for matching intervals with
    integration of the metric function.
    
    Args
      fh: file handle returned by read_file()
      dEicdt_smooth: smoothed derivative corresponding to times found in fh['t']
      interval_length: length of interval
      return_integrand: return array of integrand values

      == thresholds ==
      max_flux_threshold: threshold for maximum threshold condition of algorithm
      integral_threshold: threshold for accept/reject condition with integral
          value
      frac_above_zero_treshold: threshold for accept/reject condition with 
          fraction of ingrand values above zero
    """
    # Make interval length into timedelta
    interval_length = timedelta(seconds=interval_length)
    bar = progressbar.ProgressBar()
    matching_intervals = IntervalTree()
    integrals = []
    integrand_save = np.zeros(fh['t'].size)
    
    # Walk through timeseries 
    for start_time_idx, start_time in bar(list(enumerate(fh['t']))):
        # Determine end time of interval. If less than `interval_length` from
        # the end of the file, the interval may be less than `interval_length`.
        end_time_idx = fh['t'].searchsorted(start_time + interval_length)

        if end_time_idx > fh['t'].size:
            end_time_idx = fh['t'].size

        end_time = fh['t'][end_time_idx - 1]

        # Only check the interval if in the magnetic latitude range
        # |mlat| > 60 deg.
        mlat = fh['mlat'][start_time_idx:end_time_idx]

        if not np.all(np.abs(mlat) > 60):
            continue
        
        # Setup integrand for integration. This contains multiple
        # multiplicative terms to control different aspects of the value.
        mlat_direction = -np.sign(np.diff(np.abs(mlat)))

        ion_d_ener = fh['ion_d_ener'][:, start_time_idx:end_time_idx]
        ion_d_flux = fh['ion_d_flux'][:, start_time_idx:end_time_idx]

        num_dens = integrate_number_density(ion_d_flux, fh['ch_energy'])
        
        with np.errstate(divide='ignore'):
            max_flux = np.log10(ion_d_ener).max(axis=0)
        max_flux_mask = (max_flux > max_flux_threshold).astype(int)        

        #num_dens_mask = (num_dens > 1e10).astype(int)
        num_dens_mask = np.ones_like(num_dens)
        
        integrand = (
            mlat_direction *
            max_flux_mask[:-1] *
            num_dens_mask[:-1] *
            dEicdt_smooth[start_time_idx:end_time_idx-1]
        )

        t = fh['t'][start_time_idx:end_time_idx]
        dt = [delta.total_seconds() for delta in np.diff(t)]

        if integrand.size > 0:
            integrand_save[start_time_idx:end_time_idx-1] = integrand
        else:
            continue
        
        # Compute integral of the integrand over increments of dt using
        # rectangular integraion
        t = fh['t'][start_time_idx:end_time_idx]
        dt = [delta.total_seconds() for delta in np.diff(t)]

        integral = np.sum(integrand * dt)
        
        # Collect metadata used supplemented with interval if the interval is
        # accepted
        frac_above_zero = np.sum(integrand > 0) / integrand.size
        integrand_min = np.min(integrand)
        integrand_mean = np.mean(integrand)
        integrand_max = np.max(integrand)
        
        # The test/accept condition on the integral value and the fraction of
        # values above zero.
        if (integral > integral_threshold and
            frac_above_zero > frac_above_zero_threshold):
            matching_intervals[start_time:end_time] = {
                'integral': integral,
                'frac_above_zero': frac_above_zero,
                'integrand_min': integrand_min,
                'integrand_mean': integrand_mean,
                'integrand_max': integrand_max,
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
            min(interval.data['integral']),
            np.mean(list(interval.data['integral'])),
            max(interval.data['integral']),
            min(interval.data['frac_above_zero']),
            np.mean(list(interval.data['frac_above_zero'])),
            max(interval.data['frac_above_zero']),
            min(interval.data['integrand_min']),
            np.mean(list(interval.data['integrand_mean'])),
            max(interval.data['integrand_max']),
        ])

    df_match = pd.DataFrame(df_match_rows, columns=[
        'start_time', 'end_time',
        'integral_min', 'integral_mean', 'integral_max',
        'frac_above_zero_min', 'frac_above_zero_mean', 'frac_above_zero_max',
        'integrand_min', 'integrand_mean', 'integrand_max',
    ])

    if return_integrand:
        return df_match, integrand_save
    else:
        return df_match


def integrate_number_density(ion_d_flux, ch_energy):
    """Integrates the ion_d_flux flux (number flux) over energy using
    trapezoidal integration.

    Args
      ion_d_flux: number flux per energy
      ch_energy: channel energy
    Returns
      number density
    """
    y = ion_d_flux.copy()
    x = np.transpose([ch_energy]*y.shape[1])
    num_dens = np.trapz(y, x, axis=0)

    return num_dens


if __name__ == '__main__':
    main()
