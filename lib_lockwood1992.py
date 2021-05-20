"""Tools for calculating the reconnection rate using ionospheric plasma/field
measurements and solar wind data.

Uses method adapted from Lockwood 1992.
"""
from astropy import constants, units
from datetime import datetime, timedelta
import h5py
import math
from matplotlib.dates import date2num
import numpy as np
import pytz

from lib_util import find_moving_average

# DMSP Orbital Parameters, taken from:
#   https://directory.eoportal.org/web/eoportal/-/dmsp
DMSP_ALTITUDE = 811 * units.km
DMSP_PERIOD = 101.6 * units.min

# Length of path the particle precipation travels down along the cusp field
# line. In Lockwood 1992, it is suggested the true length is between 15-30 Re.
# This is also represented by the the d' parameter.
PRECIP_TRAVEL_PATHS = np.arange(15, 35, 5)


def read_dmsp_magn_files(dmsp_magn_files, silent=False):
    """Read DMSP magnetometer files into a single dictionary.

    Args
      dmsp_magn_files: string path to hdf files
    Returns
      dictionary mapping parameters to file
    """
    t_items = []
    Btotal_items = []
    
    for dmsp_magn_file in sorted(dmsp_magn_files):
        # Open file
        if not silent:
            print(f'Loading {dmsp_magn_file}')

        dmsp_magn_hdf = h5py.File(dmsp_magn_file, 'r')
        struct_array = dmsp_magn_hdf['Data']['Table Layout'][:]
        dmsp_magn_hdf.close()
        
        # Read the data and do simple transform step
        t = [
            datetime(1970, 1, 1, tzinfo=pytz.utc) + timedelta(seconds=i)
            for i in struct_array['ut1_unix']
        ]
        
        Btotal = (
            struct_array['bd']**2 +
            struct_array['b_forward']**2 +
            struct_array['b_perp']**2
        )

        # Append to accumulated lists
        t_items.append(t)
        Btotal_items.append(Btotal)

    # Merge arrays list of items
    omniweb_fh = {}
    omniweb_fh['t'] = np.concatenate(t_items)
    omniweb_fh['Btotal'] = np.concatenate(Btotal_items)

    return omniweb_fh

def estimate_reconn_rate(dmsp_flux_fh, Eic_smooth, i=None, j=None):
    """Estimate reconnection rate using in-situ ionospheric measurements and
    solar wind data. 

    Args
      dmsp_flux_fh: file handle (as returned by read_dmsp_flux_file)
      Eic_smooth: smoothed Eic cooresponding to times found in dmsp_flux_fh['t'].
        Note: In units of Log10(eV)
    Returns
      Ey_final: Estimated reconnection rate, in units of mV/m.
    """
    # Calculate required parameters. Interpolate everything to the time axis
    # of Eic.
    # ------------------------------------------------------------------------
    Eic = (10**Eic_smooth) * units.eV    
    t = dmsp_flux_fh['t']
    
    if i is not None and j is not None:
        Eic = Eic[i:j]
        t = t[i:j]
    
    # dEic/dt
    dEic = np.diff(Eic.value).tolist()
    dEic.append(dEic[-1])
    dEic = np.array(dEic) * units.eV

    dt = [delta.total_seconds() for delta in np.diff(t)]
    dt.append(dt[-1])
    dt = np.array(dt) * units.s
    
    dEicdt = find_moving_average(dEic/dt, 5)
    
    # Magnetc field at MP
    Bmp = 50 * units.nT

    # Magnetic field at ionosphere
    Bi = 5e-5 * units.T

    # Satellite velocity
    Vs = (
        2
        * math.pi
        * (constants.R_earth + DMSP_ALTITUDE)
        / DMSP_PERIOD
    )

    # Misc
    alpha = 0
    m = constants.m_p
    
    # Calculate the Reconnection Rate using equations derived in the paper
    # ------------------------------------------------------------------------
    Ey_all = []

    for precip_travel_path in PRECIP_TRAVEL_PATHS:
        d = precip_travel_path * constants.R_earth 
        Ey = (
            (Bi * Vs * np.cos(alpha))
            /
            (1 + (d/2) * np.sqrt(m/2) * (Eic**(-3/2)) * np.abs(dEicdt))
        )
        dy = np.sqrt(Bi / Bmp)
        Ey_final = Ey / dy       
        Ey_all.append(Ey_final.to(units.mV/units.m))

    Ey_all = np.array(Ey_all)

    return PRECIP_TRAVEL_PATHS, Ey_all