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


# Length of path the particle precipation travels down along the cusp field
# line. In Lockwood 1992, it is suggested the true length is between 15-30 Re.
# This is also represented by the the d' parameter.
PRECIP_TRAVEL_PATHS = np.array([10., 20., 30.])


def estimate_reconn_rate(dmsp_flux_fh, Eic, i=None, j=None):
    """Estimate reconnection rate using in-situ ionospheric measurements and
    solar wind data. 

    Args
      dmsp_flux_fh: file handle (as returned by read_dmsp_flux_file)
      Eic: Eic cooresponding to times found in dmsp_flux_fh['t'].
        Note: In units of Log10(eV)
    Returns
      Ey_final: Estimated reconnection rate, in units of mV/m.
    """
    # Calculate required parameters. Interpolate everything to the time axis
    # of Eic.
    # ------------------------------------------------------------------------
    Eic = (10**Eic) * units.eV    
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
    
    dEicdt = dEic/dt
    
    # Magnetc field at MP
    Bmp = 50 * units.nT

    # Magnetic field at ionosphere
    Bi = 5e-5 * units.T

    # Satellite velocity
    Vs = 7.8 * units.km / units.s

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
            (1 + (d/2) * np.sqrt(m/2) * Eic**(-3/2) * np.abs(dEicdt))
        )
        dy = np.sqrt(Bi / Bmp)
        Ey_final = Ey / dy       
        Ey_all.append(Ey_final.to(units.mV/units.m))

    Ey_all = np.array(Ey_all)

    return PRECIP_TRAVEL_PATHS, Ey_all
