"""Utility functions used by this codebase."""


from datetime import datetime, timedelta
import h5py
import numpy as np
import pytz
from spacepy import pycdf

import lib_dasilva2022


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
    n_items = []
    
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
        n_items.append(omniweb_cdf['proton_density'][:])
        
        # Close file
        omniweb_cdf.close()

    # Merge arrays list of items
    omniweb_fh = {}
    omniweb_fh['t'] = np.concatenate(t_items)
    omniweb_fh['Bx'] = np.concatenate(Bx_items)
    omniweb_fh['By'] = np.concatenate(By_items)
    omniweb_fh['Bz'] = np.concatenate(Bz_items)
    omniweb_fh['n'] = np.concatenate(n_items)
    
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
    ch_i = dmsp_flux_fh['ch_energy'].searchsorted(lib_dasilva2022.MIN_ION_VALID_ENERGY)
    ch_j = dmsp_flux_fh['ch_energy'].searchsorted(lib_dasilva2022.MAX_SHEATH_ENERGY)
    
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

