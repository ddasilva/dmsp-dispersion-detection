"""Implements an unpublished double dispersion automatic identification
method developed by da Silva. 
"""
from datetime import timedelta
from intervaltree import IntervalTree
import numpy as np
import pandas as pd
from scipy.signal import find_peaks, savgol_filter
from progressbar import ProgressBar

DEFAULT_INTEGRAL_THRESHOLD = 0.1 # Default threshold for scores
INTERVAL_LENGTH = 60             # seconds
MIN_IFLUX_AT_PEAK = 10**(6.5)    # Minimum ion flux at peak energy
MAX_PEAK_ENERGY = 10**(3.5)      # Anayze energies under this level
MIN_MLAT = 50                    # Minimum magnetic latitude (degrees)
EP_SMOOTH_WINDOW_SIZE = 11 
OMNIWEB_FILL_VALUE = 9999        # fill value for msising omniweb data


def walk_and_integrate(dmsp_flux_fh, omniweb_fh, reverse_effect, threshold):
    """Walk through windows in the file, and estimate score for each window.
    
    Args
      dmsp_flux_fh: file handle returned by read_files()
      omniweb_fh: file handle returned by read_files()
      reverse_effect: Search for effects in the opposite direction with a
        magnetic field set to the opposite of the coded threshold.
      threshold: scalar score threshold
    Returns
      df_match: Pandas dataframe holding results in each row
    """
    # Calculate coefficients -------------------------------------------------    
    mlt_coeff = (dmsp_flux_fh['mlt'] > 6) & (dmsp_flux_fh['mlt'] <= 18)
    mlt_coeff = mlt_coeff.astype(float)

    mlat_coeff = np.abs(dmsp_flux_fh['mlat']) > MIN_MLAT
    mlat_coeff = mlat_coeff.astype(float)

    mlat_dir_coeff = np.zeros(dmsp_flux_fh['t'].size)
    mlat_dir_coeff[:-1] = -np.sign(np.diff(np.abs(dmsp_flux_fh['mlat'])))
    mlat_dir_coeff[-1] = mlat_dir_coeff[-2]

    if reverse_effect:
        mlat_dir_coeff *= -1
    
    # Iterate through time
    # ------------------------------------------------------------------------
    matching_intervals = IntervalTree()
    bar = ProgressBar()
    prog_bar_iter = bar(list(enumerate(dmsp_flux_fh['t'][:-1])))

    for start_time_idx, start_time in prog_bar_iter:
        # Set end time index ------------------------------------------------
        end_time = start_time + timedelta(seconds=INTERVAL_LENGTH)
        end_time_idx = dmsp_flux_fh['t'].searchsorted(end_time)

        if end_time_idx > dmsp_flux_fh['t'].size:
            end_time_idx = dmsp_flux_fh['t'].size

        # Skip condition -----------------------------------------------------
        i, j = start_time_idx, end_time_idx

        if j - i < 8:
            continue
        
        if np.all(mlt_coeff[i:j] * mlt_coeff[i:j] == 0):
            continue
            
        # Do peakfinding -----------------------------------------------------
        t, lower_Ep, upper_Ep = calculate_dual_Ep(
            dmsp_flux_fh, start_i=i, stop_j=j
        )
        
        # Calculate Derivative of lower_Ep and upper_Ep ----------------------
        # Convert to log space 
        lower_Ep = np.log10(lower_Ep)
        upper_Ep = np.log10(upper_Ep)

        # smooth
        mask = np.ones(t.size, dtype=bool)
        lower_Ep_smooth = clean_Ep(lower_Ep, mask, EP_SMOOTH_WINDOW_SIZE)
        upper_Ep_smooth = clean_Ep(upper_Ep, mask, EP_SMOOTH_WINDOW_SIZE)
        
        # Calculate difference in times (should be uniformly spaced)
        dt = [delta.total_seconds() for delta in np.diff(t)]
        dt.append(dt[-1])                       # copy last value to retain shape
        dt = np.array(dt)

        # Calculate difference in Ep's
        lower_Ep_deriv = np.zeros(t.size)
        upper_Ep_deriv = np.zeros(t.size)
        
        lower_Ep_deriv[:-1] = np.diff(lower_Ep_smooth) / dt[:-1]
        upper_Ep_deriv[:-1] = np.diff(upper_Ep_smooth) / dt[:-1]
        
        lower_Ep_deriv[-1] = lower_Ep_deriv[-2]
        upper_Ep_deriv[-1] = upper_Ep_deriv[-2]
        
        # Set derivative = 0 at nulls
        lower_Ep_deriv[np.isnan(lower_Ep_deriv)] = 0
        upper_Ep_deriv[np.isnan(upper_Ep_deriv)] = 0

        # Look up magnetic field value at center of the interval -------------
        B_test_time = start_time + timedelta(seconds=INTERVAL_LENGTH/2)
        B_i = omniweb_fh['t'].searchsorted(B_test_time)
        if B_i == omniweb_fh['Bz'].size:
            continue
        Bx = omniweb_fh['Bx'][B_i]
        By = omniweb_fh['By'][B_i]
        Bz = omniweb_fh['Bz'][B_i]

        if Bz > OMNIWEB_FILL_VALUE:
            continue
        #elif reverse_effect and Bz < 0:
        #    continue
        elif not reverse_effect and Bz > 0:
            continue
         
        # Calculate integrand ----------------------------------------------------
        lower_integrand = (
            mlt_coeff[i:j] *
            mlat_coeff[i:j] *
            mlat_dir_coeff[i:j] *
            lower_Ep_deriv
        )
        upper_integrand = (
            mlt_coeff[i:j] *
            mlat_coeff[i:j] *
            mlat_dir_coeff[i:j] *
            upper_Ep_deriv
        )
        
        # Calculate total score for upper and lower curves
        lower_score = np.sum(lower_integrand * dt)
        upper_score = np.sum(upper_integrand * dt)
        
        if lower_score > threshold and upper_score > threshold:
            matching_intervals[start_time:end_time] = {
                'Bx': Bx,
                'By': By,
                'Bz': Bz,
                't': t,
                'lower_integrand': lower_integrand,
                'upper_integrand': upper_integrand,
                'lower_Ep': lower_Ep,
                'upper_Ep': upper_Ep,
            }

    # Merge overlapping intervals into common intervals, retaining context
    # associated with each.
    # ------------------------------------------------------------------------
    def reducer(current, new_data):
        for key in new_data:
            if key not in current:
                current[key] = []
            current[key].append(new_data[key])

        return current

    matching_intervals.merge_overlaps(
        data_reducer=reducer, data_initializer={}
    )

    # Convert to pandas dataframe for easy output to terminal of matching
    # intervals and associated metadata.
    # ------------------------------------------------------------------------
    df_match_rows = []
    
    for interval in sorted(matching_intervals):
        df_match_rows.append([
            interval.begin, interval.end,
            np.mean(list(interval.data['Bx'])),
            np.mean(list(interval.data['By'])),
            np.mean(list(interval.data['Bz'])),
            interval.data['t'][0],
            interval.data['lower_integrand'][0],
            interval.data['upper_integrand'][0],
            interval.data['lower_Ep'][0],
            interval.data['upper_Ep'][0],
        ])

    df_match = pd.DataFrame(df_match_rows, columns=[
        'start_time', 'end_time', 'Bx_mean', 'By_mean', 'Bz_mean',
        't', 'lower_integrand', 'upper_integrand', 'lower_Ep', 'upper_Ep', 
    ])

    return df_match

    
def calculate_dual_Ep(dmsp_flux_fh, start_i=None, stop_j=None, _cache={}):
    """Calculate two time series associated with lower and upper peak
    energies.

    Args
      dmsp_flux_fh: Dicitonary mapping variable names to arrays, associated with
        DMSP precipitating flux file.
      start_i: Start index of the file to analyze, defaults to start of file
      stop_j: Stop index (exclusive) to analyze, defaults to end of file
    Returns
      t: time axis associated with other return values, array of datetimes
      lower_Ep: Lower curve of peak energies
      upper_Ep: Upper curve of peak energies
    """
    # Iinitialize default keyword arguments to full size of the time series
    # available
    # ------------------------------------------------------------------------
    if start_i is None:
        start_i = 0
    if stop_j is None:
        stop_j = dmsp_flux_fh['t'].size

    # Initialize vairables that will be used during iteration
    # ------------------------------------------------------------------------
    t = dmsp_flux_fh['t']
    E = dmsp_flux_fh['ch_energy']
    lower_Ep = np.zeros(stop_j - start_i, dtype=float)
    upper_Ep = np.zeros(stop_j - start_i, dtype=float)
    lower_Flux = np.zeros(stop_j - start_i, dtype=float)
    upper_Flux = np.zeros(stop_j - start_i, dtype=float)
    
    last_filled = 'lower'

    # Iterate through the timesteps. Performs peakfinding and branches logic
    # based on number of peaks.
    # ------------------------------------------------------------------------
    for ii, i in enumerate(range(start_i, stop_j)):
        data = dmsp_flux_fh['ion_d_ener'][:, i]
        
        key = (id(dmsp_flux_fh['ion_d_ener']), i)

        if key in _cache:
            smoothed = _cache[key]
        else:
            smoothed = savgol_filter(data, 5, 2)
            _cache[key] = smoothed

        #smoothed = savgol_filter(data, 5, 2)
        #smoothed = np.convolve(data, kernel, mode='same')
        #smoothed = data

        peaks = find_peaks(smoothed)[0]
        peaks = np.array([peak for peak in peaks if
                          data[peak] > MIN_IFLUX_AT_PEAK and
                          E[peak] < MAX_PEAK_ENERGY])
        
        if len(peaks) == 0:
            lower_Ep[ii] = np.nan
            lower_Flux[ii] = np.nan
            upper_Ep[ii] = np.nan
            upper_Flux[ii] = np.nan
            last_filled = 'none'
            
        elif len(peaks) == 1:
            Enext = E[peaks[0]]
            Fluxnext = data[peaks[0]]

            cond0 = (ii == 0)
            cond1 = (np.isnan(lower_Ep[ii-1]) and np.isnan(upper_Ep[ii-1]))
            cond2 = (np.isfinite(lower_Ep[ii-1]) and lower_Ep[ii-1] >= Enext)
            cond3 = (np.isnan(upper_Ep[ii-1]))
            cond4 = (last_filled in ('none', 'lower', 'both')) 

            if ((cond0 or cond1 or cond2 or cond3) and cond4):
                lower_Ep[ii] = Enext
                lower_Flux[ii] = Fluxnext
                upper_Ep[ii] = np.nan
                upper_Flux[ii] = np.nan
                last_filled = 'lower'                
            else:
                lower_Ep[ii] = np.nan
                lower_Flux[ii] = np.nan
                upper_Ep[ii] = Enext
                upper_Flux[ii] = Fluxnext
                last_filled = 'upper'

        elif len(peaks) >= 2:
            # Pick peaks with top-2 fluxes
            I = np.argsort(smoothed[peaks])[-2:]
            low_E, hi_E = sorted(E[peaks[I]])
            low_Flux, hi_Flux = data[peaks][np.argsort(E[peaks[I]])]
            
            # if high is closer to last lower line, then probably a false
            # signal
            if np.abs(hi_E - lower_Ep[ii-1]) < np.abs(low_E - lower_Ep[ii-1]):
                lower_Ep[ii] = hi_E
                upper_Ep[ii] = np.nan
                last_filled = 'lower'
            else:
                lower_Ep[ii], upper_Ep[ii] = low_E, hi_E
                lower_Flux[ii], upper_Flux[ii] = low_Flux, hi_Flux
                last_filled = 'both'
            
    #lower_Ep = np.log10(lower_Ep)
    #upper_Ep = np.log10(upper_Ep)

    # Filter out lone points -------------------------------------------------    
    # for i in range(1, lower_Ep.size - 1):
    #     isfin1 = np.isfinite(lower_Ep[i-1:i+2])
    #     isfin2 = np.isfinite(upper_Ep[i-1:i+2])
    
    #     if np.all(isfin1 == [False, True, False]):
    #         lower_Ep[i] = np.nan
    #     if np.all(isfin2 == [False, True, False]):
    #         upper_Ep[i] = np.nan

    # Return
    return t[start_i:stop_j], lower_Ep, upper_Ep



#@jit(nopython=True)
def clean_Ep(Ep, keep_mask, window_size):
    """Smooth Ep with a mask of points to include in moving average.
    
    Arguments
      Ep: floating point numpy array
      keep_mask: points to include in moving average
      window_size: integer, must be odd
    Returns
      Smoothed Ep array
    """    
    assert (window_size is None) or (window_size % 2 == 1),\
        'Window size must be odd'

    Ep_clean = Ep.copy()
    
    for i in range(Ep.size):
        if not keep_mask[i]:
            Ep_clean[i] = np.nan
        elif (window_size is None) and keep_mask[i]:
            Ep_clean[i] = Ep[i]
        elif (window_size is not None):
            total = 0.0
            count = 0
            
            for di in range(-window_size//2, window_size//2 + 1):
                if i + di > 0 and i + di < Ep.size and keep_mask[i + di]:
                    total += Ep[i + di]
                    count += 1
            
            if count > 0:  # else left as nan
                Ep_clean[i] = total / count
                
    return Ep_clean

