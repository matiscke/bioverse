""" Define new functions for planet simulation here. Function arguments should be provided a default value."""

# Add import statements as necessary
import numpy as np

import os
from astropy.coordinates import SkyCoord
import astropy.units as u


# Bioverse modules and constants
from .classes import Table
from . import util
from .util import CATALOG
from .constants import CONST, ROOT_DIR

def read_stars_Gaia(d, filename='gcns_catalog.dat', d_max=120., M_st_min=0.075, M_st_max=2.0, R_st_min=0.095,
                    R_st_max=2.15, T_min=0., T_max=10., inc_binary=0):  # , mult=0):
    """ Reads a list of stellar properties from the Gaia nearby stars catalog.

    Parameters
    ----------
    d : Table
        An empty Table object.
    filename : str, optional
        Filename containing the Gaia target catalog.
    d_max : float, optional
        Maximum distance to which to simulate stars, in parsecs.
    M_st_min : float, optional
        Minimum stellar mass, in solar units.
    M_st_max : float, optional
        Maximum stellar mass, in solar units.
    R_st_min : float, optional
        Minimum stellar radius, in solar units.
    R_st_max : float, optional
        Maximum stellar radius, in solar units.
    T_min : float, optional
        Minimum stellar age, in Gyr.
    T_max : float, optional
        Maximum stellar age, in Gyr.
    inc_binary : bool, optional
        Include binary stars? Default = False.
    mult : float, optional
        Multiple on the total number of stars simulated. If > 1, duplicates some entries from the LUVOIR catalog.

    Returns
    -------
    d : Table
        Table containing the sample of real stars.
    """
    # Read the catalog with column names
    path = filename if os.path.exists(filename) else ROOT_DIR + '/' + filename
    catalog = np.genfromtxt(path, unpack=False, names=True, dtype=None, encoding=None)
    for name in catalog.dtype.names:
        d[name.strip()] = list(catalog[name])  # *int(mult)

    # Missing values (TODO: this part is ironically incomplete)
    if 'd' not in d.keys():
        d['d'] = np.cbrt(np.random.uniform(0, d_max ** 3, len(d)))
    if 'x' not in d.keys():
        cost, phi = np.random.uniform(-1, 1, len(d)), np.random.uniform(0, 2 * np.pi, len(d))
        r = d['d'] * np.sin(np.arccos(cost))
        d['x'], d['y'], d['z'] = r * np.cos(phi), r * np.sin(phi), d['d'] * cost
    if 'age' not in d.keys():
        d['age'] = np.random.uniform(T_min, T_max, size=len(d))
    if 'logL' not in d.keys():
        d['logL'] = np.log10(d['L_st'])
    if 'star_name' not in d.keys():
        d['star_name'] = np.char.array(np.full(len(d), 'REAL-')) + np.char.array(np.arange(len(d)).astype(str))
    if 'RV' not in d.keys():
        d['RV'] = np.random.uniform(-200, 200, size=len(d))

    # Enforce a maximum distance
    d = d[d['d'] < d_max]

    # Enforce a min/max mass & radius
    d = d[(d['M_st'] < M_st_max) & (d['M_st'] > M_st_min)]
    d = d[(d['R_st'] < R_st_max) & (d['R_st'] > R_st_min)]

    # Include/exclude stars in binary systems, default is (0/False) exclude binaries
    if inc_binary == 0:
        d = d[(d['binary'] == False)]
        # d.reset_index(inplace=True,drop=True)

    # Assign stellar IDs and names
    d['starID'] = np.arange(len(d), dtype=int)
    d['simulated'] = np.zeros(len(d), dtype=bool)

    # Draw a random age for each system
    d['age'] = np.random.uniform(T_min, T_max, len(d))

    # Add ecliptic coordinates
    ra = np.array(d['ra'])
    dec = np.array(d['dec'])
    dist = np.array(d['d'])
    c = SkyCoord(ra * u.degree, dec * u.degree, distance=dist * u.pc, frame='icrs', equinox='J2016.0')
    c_ec = c.transform_to('heliocentrictrueecliptic')
    d['helio_ecl_lon'] = np.array(c_ec.lon.value)
    d['helio_ecl_lat'] = np.array(c_ec.lat.value)

    return d
import os
from astropy.coordinates import SkyCoord
import astropy.units as u


# Bioverse modules and constants
from .classes import Table
from . import util
from .util import CATALOG
from .constants import CONST, ROOT_DIR

def read_stars_Gaia(d, filename='gcns_catalog.dat', d_max=120., M_st_min=0.075, M_st_max=2.0, R_st_min=0.095,
                    R_st_max=2.15, T_min=0., T_max=10., inc_binary=0, seed=42):  # , mult=0):
    """ Reads a list of stellar properties from the Gaia nearby stars catalog.

    Parameters
    ----------
    d : Table
        An empty Table object.
    filename : str, optional
        Filename containing the Gaia target catalog.
    d_max : float, optional
        Maximum distance to which to simulate stars, in parsecs.
    M_st_min : float, optional
        Minimum stellar mass, in solar units.
    M_st_max : float, optional
        Maximum stellar mass, in solar units.
    R_st_min : float, optional
        Minimum stellar radius, in solar units.
    R_st_max : float, optional
        Maximum stellar radius, in solar units.
    T_min : float, optional
        Minimum stellar age, in Gyr.
    T_max : float, optional
        Maximum stellar age, in Gyr.
    inc_binary : bool, optional
        Include binary stars? Default = False.
    seed : int, optional
        seed for the random number generators.
    mult : float, optional
        Multiple on the total number of stars simulated. If > 1, duplicates some entries from the LUVOIR catalog.

    Returns
    -------
    d : Table
        Table containing the sample of real stars.
    """

    np.random.seed(seed)

    # Read the catalog with column names
    path = filename if os.path.exists(filename) else ROOT_DIR + '/' + filename
    catalog = np.genfromtxt(path, unpack=False, names=True, dtype=None, encoding=None)
    for name in catalog.dtype.names:
        d[name.strip()] = list(catalog[name])  # *int(mult)

    # Missing values (TODO: this part is ironically incomplete)
    if 'd' not in d.keys():
        d['d'] = np.cbrt(np.random.uniform(0, d_max ** 3, len(d)))
    if 'x' not in d.keys():
        cost, phi = np.random.uniform(-1, 1, len(d)), np.random.uniform(0, 2 * np.pi, len(d))
        r = d['d'] * np.sin(np.arccos(cost))
        d['x'], d['y'], d['z'] = r * np.cos(phi), r * np.sin(phi), d['d'] * cost
    if 'age' not in d.keys():
        d['age'] = np.random.uniform(T_min, T_max, size=len(d))
    if 'logL' not in d.keys():
        d['logL'] = np.log10(d['L_st'])
    if 'star_name' not in d.keys():
        d['star_name'] = np.char.array(np.full(len(d), 'REAL-')) + np.char.array(np.arange(len(d)).astype(str))
    if 'RV' not in d.keys():
        d['RV'] = np.random.uniform(-200, 200, size=len(d))

    # Enforce a maximum distance
    d = d[d['d'] < d_max]

    # Enforce a min/max mass & radius
    d = d[(d['M_st'] < M_st_max) & (d['M_st'] > M_st_min)]
    d = d[(d['R_st'] < R_st_max) & (d['R_st'] > R_st_min)]

    # Include/exclude stars in binary systems, default is (0/False) exclude binaries
    if inc_binary == 0:
        d = d[(d['binary'] == False)]
        # d.reset_index(inplace=True,drop=True)

    # Assign stellar IDs and names
    d['starID'] = np.arange(len(d), dtype=int)
    d['simulated'] = np.zeros(len(d), dtype=bool)

    # Draw a random age for each system
    d['age'] = np.random.uniform(T_min, T_max, len(d))

    # Add ecliptic coordinates
    ra = np.array(d['ra'])
    dec = np.array(d['dec'])
    dist = np.array(d['d'])
    c = SkyCoord(ra * u.degree, dec * u.degree, distance=dist * u.pc, frame='icrs', equinox='J2016.0')
    c_ec = c.transform_to('heliocentrictrueecliptic')
    d['helio_ecl_lon'] = np.array(c_ec.lon.value)
    d['helio_ecl_lat'] = np.array(c_ec.lat.value)

    return d

def M_R(observed):
    """Just for this test, add a silly mass measurement"""
    observed['M'] = observed['R']*2
    return observed


def label_lateM(table, Mst_threshold=0.4):
    # is_lateM=False for all planets
    table['is_late'] = np.zeros(len(table))

    table['is_late'] = table['M_st'] < Mst_threshold
    return table


def occurrence_hypo(theta, X):
    """Define the hypothesis that the planet occurrence is lower around late Ms than around early Ms.
    """
    return None


def magma_ocean(d, funform='exp_decay', f_magma=0.1, lambda_a=2., a_cut=0.1, radius_reduction=0.2):
    """Assign a fraction of planets global magma oceans that change the planet's radius.

    Parameters:
    -----------
    d: Table
        The population of planets.
    funform: str
        The functional form of the magma ocean probability as a function of effective semi-major axis (sma).
        Can be 'exp_decay' or 'step'.
    f_magma: float
        The fraction of planets that have global magma oceans.
    lambda_a: float
        Decay parameter for the semi-major axis dependence of having a global magma ocean.
    a_cut: float
        cutoff effective sma for magma oceans. Defines position of the exponential decay or step.
    radius_reduction: float
        The fraction by which a planet's radius is reduced due to a global magma ocean.

    Returns
    -------
    d : Table
        Table containing the sample of simulated planets with new columns 'has_magmaocean' and 'is_small'.

    """
    # randomly assign planets to have a magma ocean, depending on the semi-major axes
    if funform == 'exp_decay':
        # assign magma ocean with a probability following a exponential decay function
        P_magma = f_magma * np.exp(-(d['a_eff']/a_cut)**lambda_a) # HAS TO BE REPLACED WITH MODEL OUTPUT FOR MAGMA OCEAN PLANETS
    elif funform == 'step':
        # assign a magma ocean with probability f_magma, but only to planets with a_eff < a_cut
        P_magma = f_magma * (d['a_eff'] < a_cut)
    else:
        raise ValueError('funform must be either "exp_decay" or "step"')

    d['has_magmaocean'] = np.random.uniform(0, 1, len(d)) < P_magma

    # reduce the radius of the planets with magma oceans
    mask = d['has_magmaocean']
    try:
        d.loc[mask,'R'] *= (1 - radius_reduction) # HAS TO BE REPLACED WITH MODEL OUTPUT FOR MAGMA OCEAN PLANETS
    except AttributeError:
        d['R'][mask] *= (1 - radius_reduction) # HAS TO BE REPLACED WITH MODEL OUTPUT FOR MAGMA OCEAN PLANETS

    # define planets with smaller radius than expected. THIS SHOULD BE REPLACED WITH SOMETHING MORE REALISTIC
    d['is_small'] = d['R'] < np.mean(d['R'])

    return d
