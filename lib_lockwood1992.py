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
PRECIP_TRAVEL_PATHS = np.array([10.])


def estimate_reconn_rate(dmsp_flux_fh, Eic, i=None, j=None):
    """Estimate reconnection rate using in-situ ionospheric measurements and
    solar wind data. 

    Args
      dmsp_flux_fh: file handle (as returned by read_dmsp_flux_file)
      Eic: Eic cooresponding to times found in dmsp_flux_fh['t'].
        Note: In units of Log10(eV)
    Returns
      precip_travel_paths: precipitation travel paths (d' in paper)
      Ey_iono: reconnection rate at ionosphere (Ey in paper)
      Ey_mpause: reconnetion rate at magnetopause (Ey' in paper)
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

    # Estimate error in numerical derivative from Taylor series remainder
    ddEicddt = np.gradient(np.gradient(Eic, t, edge_order=2), t, edge_order=2)
    dEicdt_errest = dt*np.abs(ddEicddt)/2
    
    # Magnetc field at MP
    Bmp = 50 * units.nT

    # Magnetic field at ionosphere
    Bi = 50e3 * units.nT

    # Satellite velocity
    Vs = 7.8 * units.km / units.s

    # Misc
    alpha = 0
    m = constants.m_p
    
    # Calculate the Reconnection Rate using equations derived in the paper
    # ------------------------------------------------------------------------
    Ey_iono = []         # reconnection rate at ionosphere
    Ey_mpause = []       # reconnection rate at magnetopause

    for precip_travel_path in PRECIP_TRAVEL_PATHS:
        d = precip_travel_path * constants.R_earth 
        Ey = (
            (Bi * Vs * np.cos(alpha))
            /
            (1 + (d/2) * np.sqrt(m/2) * Eic**(-3/2) * np.abs(dEicdt))
        )

        Ey_errest = np.abs(
            (Bi * Vs * np.cos(alpha)*(d/2) * np.sqrt(m/2) * Eic**(-3/2))
            /
            (1 + (d/2) * np.sqrt(m/2) * Eic**(-3/2) * np.abs(dEicdt))**2.0
        )*dEicdt_errest

        dy = np.sqrt(Bi / Bmp)
        
        Ey_final = Ey / dy
        Ey_final_errest = Ey_errest / dy

        Ey_iono.append(Ey.to(units.mV/units.m))
        Ey_iono_errest.append(Ey_errest.to(units.mV/units.m))
        Ey_mpause.append(Ey_final.to(units.mV/units.m))

    Ey_iono = np.array(Ey_iono)
    Ey_iono_errest = np.array(Ey_iono_errest)
    Ey_mpause = np.array(Ey_mpause)
    Ey_mpause_errest = np.array(Ey_mpause_errest)

    # Remove points where DeltaEic=0
    for i in range(PRECIP_TRAVEL_PATHS.size):
        Ey_iono[i, :-1][np.diff(Eic)==0] = np.nan
        Ey_mpause[i, :-1][np.diff(Eic)==0] = np.nan
    
    return PRECIP_TRAVEL_PATHS, Ey_iono, Ey_mpause, Ey_errest, Ey_mpause_errest
