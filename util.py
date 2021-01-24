# Python imports
from astropy.modeling.blackbody import blackbody_lambda
from astropy import units as u
from ast import literal_eval as ast_literal_eval
import copy
import gc
import glob
import importlib.util
import inspect
import matplotlib
#matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
import numpy as np
import os
from pandas import DataFrame
import pickle
import scipy
from scipy import stats
import shutil
import sys
import time
from warnings import warn

# Top-level code directory
ROOT_DIR = os.path.dirname(os.path.realpath(__file__))
UI_DIR = ROOT_DIR+'/UI/'
ATMOSPHERE_TEMPLATES_DIR = ROOT_DIR+'/Templates/Atmospheres/'
SURVEYS_DIR = ROOT_DIR+'/Surveys/'
MODELS_DIR = ROOT_DIR+'/Objects/Models/'
GENERATORS_DIR = ROOT_DIR+'/Generators/'
INSTRUMENTS_DIR = ROOT_DIR+'/Instruments/'
OBJECTS_DIR = ROOT_DIR+'/Objects/'
PLOTS_DIR = ROOT_DIR+'/Plots/'
RESULTS_DIR = ROOT_DIR+'/Results/'
CATALOG_FILE = ROOT_DIR+'/Gaia.csv'

# Program version
VERSION = "0.0.1"

# HELA install directory (this will be moved to a config file later)
# HELA_DIR = ROOT_DIR+'/../HELA/'
# sys.path.append(HELA_DIR)
# import rfretrieval

# Physical constants (in cgs where applicable)
CONST = {}
CONST['T_eff_sol'] = 5777.
CONST['yr_to_day'] = 365.2422
CONST['AU_to_solRad'] = 215.03215567054767
CONST['rho_Earth'] = 5.51
CONST['g_Earth'] = 980.7
CONST['amu_to_kg'] = 1.66054e-27
CONST['R_Earth'] = 6.371e8
CONST['R_Sun'] = 6.9634e10
CONST['h_Earth'] = 8.5e5

# Data types
ARRAY_TYPES = (np.ndarray,list,tuple)
LIST_TYPES = ARRAY_TYPES
FLOAT_TYPES = (float,np.float,np.float_,np.float64)
INT_TYPES = (int,np.int_,np.int64,np.integer,np.int8)
STR_TYPES = (str,np.str_)
BOOL_TYPES = (bool,np.bool_)

# Load the Gaia stellar target catalog into memory for fast access
try:
    CATALOG = np.genfromtxt(CATALOG_FILE, delimiter=',', names=True)
except FileNotFoundError:
    warn("could not load {:s} - try running util.update_stellar_catalog")

# Progress bar if tqdm is installed, else a dummy function
try:
    from tqdm import tqdm
except ImportError:
    tqdm = None
    
def bar(arg,do_bar=True):
    if tqdm is not None and do_bar:
        return tqdm(arg)
    else:
        return arg

# Function to get the input value types (try: int, float, bool, str)
def get_type(x):
    try:
        int(x)
        return int
    except ValueError:
        try:
            float(x)
            return float
        except ValueError:
            if x.strip().upper() in ['TRUE','FALSE']:
                return bool
            else:
                return str

# Imports a function given the filename and the name of the function
def import_function_from_file(function_name,filename):
        # Get the module name from the filename (assume they are the same minus .py)
        module_name = '.'.join(filename.strip('/').split('/')[-1].split('.')[:-1])

        # Import the module
        spec = importlib.util.spec_from_file_location(module_name,filename)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)

        # Return the function
        return mod.__dict__[function_name]

def get_description_from_function(function):
    """ Gets a function description from its docstring. Reads up until the first double
    line break. For example, this second line will be read.

    This third line will not be read.

    Parameters
    ----------
    function : function
        Function to pull the description from.

    Returns
    -------
    description : str
        Returns the first paragraph of the function's docstring, or '' if there isn't one.
    """
    docstring = function.__doc__
    if docstring is None:
        return ''
    
    # Get the lines corresponding to the first paragraph
    description = ''
    for line in docstring.strip().split('\n'):
        if line == '':
            break
        description += line.replace('  ',' ').replace('  ',' ').replace('  ',' ')
    return description

# Creates a new "program" using the list of functions from a CSV file (e.g. program.dat)
def create_program(filename=ROOT_DIR+'/program.dat'):
    program = []
    with open(filename,'r') as f:
        for line in f.readlines():
            if line[0].strip() == '#': continue
            function_name,path = line.split('#')[0].split(',')[:2]
            description = ','.join(line.split('#')[0].split(',')[2:])

            # Load the function and gets its default arguments
            func = import_function_from_file(function_name,path)
            items = list(inspect.signature(func).parameters.items())[1:] # Chops off the first two arguments (s,p)
            args = {}
            for k,v in items:
                # Default argument description is "None"
                args[k] = [v.default,None]

            program.append([function_name,path,description,args])
    return np.array(program)

# Load the program from a .pkl file
def import_program(filename=ROOT_DIR+'/program.pkl'):
    return np.array(pickle.load(open(filename,'rb')))

# Saves the program to a .pkl file, including modified arguments
def export_program(program,filename=ROOT_DIR+'/program.pkl'):
    pickle.dump(program,open(filename,'wb'))

# Adjusts the position of an item in the program
def program_adjust(program,idx,shift):
    shift = int(shift)
    s2 = np.copy(program)
    if shift == 1 and (idx < len(program)-1):
        s2[[idx,idx+1]] = s2[[idx+1,idx]]
    elif shift == -1 and (idx > 0):
        s2[[idx-1,idx]] = s2[[idx,idx-1]]
    return s2

# Returns the "colors" of a planet based on its class and orbit
def get_planet_colors(d):
    c_or = 'blue' if d['class1'] == 'cold' else 'green' if d['class1'] == 'warm' else 'red'
    c_pl = 'red' if d['class2'] == 'rocky' else 'green' if d['class2'] == 'super-Earth' else 'blue'
    return c_pl,c_or

# Cycles through the values in a list
def cycle_index(vals,val,delta):
    N = len(vals)
    idx = list(vals).index(val)
    idx_new = idx + delta
    if idx_new < 0: idx_new = N + idx_new
    if idx_new >= N: idx_new = idx_new - N
    return vals[idx_new]

# Fills an array with the appropriate "nan" value for its type
def nan_fill(a,dtype=None):
    return np.full(a.shape,np.nan,dtype=a.dtype if dtype is None else dtype)

# Translates an array of object counts (e.g. [1,2,2,1,3]) into a longer list
# of the individual objects' order of creation (e.g., [1,1,2,1,2,1,1,2,3])
# Uses a numpy trick from: https://codereview.stackexchange.com/questions/83018/vectorized-numpy-version-of-arange-with-multiple-start-stop
def get_order(N):
    start = np.repeat(np.zeros(len(N)),N)
    counter = np.arange(N.sum())-np.repeat(N.cumsum()-N,N)
    return (start+counter).astype(int)

# Translates the 'subset' specified for a model into a mask
def mask_from_model_subset(pl,subset):
    mask = True
    for key in subset.keys():
        val = subset[key]
        if type(val) in LIST_TYPES:
            mask = mask & ((pl[key]>val[0])&(pl[key]<val[1]))
        else:
            mask = mask & (pl[key]==val)
    return mask

def compute_bin_centers(bins):
    """ Given a set of N bin edges, returns N-1 bin centers and half-widths. """
    return (bins[...,1:]+bins[...,:-1])/2., (bins[...,1:]-bins[...,:-1])/2.

def compute_eta_Earth(d, by_type=True):
    """ Computes the value of eta Earth for a simulated sample of planets. Note this could be inaccurate if
    there are stars without planets which are usually not listed in the simulated sample, although the algorithm
    does attempt to correct for this.

    Parameters
    ----------
    d : Table
        Simulated sample of planets.
    by_type : bool, optional
        If True, calculate eta Earth separately for each spectral type.
    """
    # Determine the fraction of stars with planets
    f_pl = len(np.unique(d['starID']))/max(d['starID'])

    # Calculate eta Earth for all stars
    N_st = max(d['starID'])
    eta = d['EEC'].sum()/N_st
    print("eta Earth = {:.2f} per star".format(eta))

    # Calculate eta Earth for each spectral type
    if by_type:
        for SpT in ['F','G','K','M']:
            # Determine the number of stars for this spectral type, accounting for stars without planets
            mask = d['SpT'] == SpT
            if mask.sum() == 0:
                continue
            N_st = len(np.unique(d['starID'][mask]))/f_pl
            eta = d['EEC'][mask].sum()/N_st
            print("          = {:.3f} per {:s} star (+- {:.5f})".format(eta, SpT, d['EEC'][mask].sum()**0.5/N_st))

        # Also calculate for M1-M6 and M7-L0 stars
        mask = (d['T_eff_st'] > 3050) & (d['T_eff_st'] <= 3900)
        N_st = len(np.unique(d['starID'][mask]))/f_pl
        eta = d['EEC'][mask].sum()/N_st
        print("          = {:.3f} per M1-M4 star (+- {:.5f})".format(eta, d['EEC'][mask].sum()**0.5/N_st))

        mask = (d['T_eff_st'] > 1950) & (d['T_eff_st'] <= 3050)
        N_st = len(np.unique(d['starID'][mask]))/f_pl
        eta = d['EEC'][mask].sum()/N_st
        print("          = {:.3f} per M5-L2 star (+- {:.5f})".format(eta, d['EEC'][mask].sum()**0.5/N_st))

def compute_occurrence_multiplier(optimistic=False, optimistic_factor=3, N_pts=30):
    """ Determines the multiplier for occurrence rates and planet periods as a function of stellar mass. """
    # Scaling factors for spectral type (Mulders+2015, Table 1 and Figure 4)
    M_st_0 = [0.17, 0.69, 0.97, 1.43]
    f_N = 1/np.array([0.35, 0.55, 0.75, 1.0])
    f_a = 1/np.array([1.6, 1.4, 1.2, 1.0])


    # Edit the M star point to be 3.5x solar
    f_N[0] = 3.5 * f_N[2]

    # In the "optimistic" case, assume TRAPPIST-1 analogs (M = 0.08) have 5x as many planets as the typical Kepler star
    if optimistic:
        # Add a 0.8 solar mass point with 10x as many planets as solar type stars
        M_st_0 = np.append([0.08], M_st_0)
        f_N = np.append([optimistic_factor*f_N[0]], f_N)
        f_a = np.append([f_a[0]], f_a)

    # Fit a second-order polynomial to f_N and f_a
    p_N = np.polyfit(M_st_0, f_N, deg=2)
    p_a = np.polyfit(M_st_0, f_a, deg=2)

    # Evaluate over a grid of masses and convert semi-major axis to period
    x = np.logspace(np.log10(0.05), np.log10(2), N_pts)
    #y_N = np.polyval(p_N, x)
    #y_P = np.polyval(p_a, x)**1.5
    y_N = np.interp(x, M_st_0, f_N)
    y_P = np.interp(x, M_st_0, f_a)**1.5

    # Normalize these factors so that y = 1 for the typical Kepler planet host
    # (M ~ 0.965 per Kopparapu+2018)
    M_st_typ = 0.965
    #y_N /= np.polyval(p_N, M_st_typ)
    #y_P /= np.polyval(p_a, M_st_typ)**1.5
    y_N /= np.interp(M_st_typ, M_st_0, f_N)
    y_P /= np.interp(M_st_typ, M_st_0, f_a)**1.5

    return x, y_N, y_P

def update_stellar_catalog(d_max=100, filename=CATALOG_FILE):
    """ Updates the catalog of nearby sources from Gaia DR2 and saves it to a file. Requires astroquery. """
    from astroquery.gaia import Gaia
    
    Nmax = 100000
    query = "SELECT TOP {:d} source_id, parallax, ra, dec, teff_val, phot_g_mean_mag".format(Nmax)+\
            " FROM gaiadr2.gaia_source"+\
            " WHERE (parallax>={:f}) AND (teff_val>=4000)".format(1000/d_max)+\
            " ORDER BY parallax DESC"

    job = Gaia.launch_job(query)
    table = job.get_results()
    table.write(filename, overwrite=True)

    if len(table) == Nmax:
        print("Warning! Gaia query exceeds row limit!")
    print("Catalog updated. Don't forget to restart!")

def get_xyz(pl,t=0,M=None,n=3):
    # Computes the x/y/z separation in AU at time(s) t, assuming the mean longitude at t=0 is M0 and the unit is days 
    # n determines the number of iterations of Newton's method for solving Kepler's equation
    
    # Determine the mean longitude at time(s) t
    if M is None:
        M = pl['M0']+(2*np.pi*t)/pl['P']
    
    # Eccentric orbit
    if np.any(pl['e'] != 0):
        # Increment M by 2pi (so the solver doesn't break near M = 0)
        M += 2 * np.pi
        
        # Eccentric anomaly (Kepler's equation solved w/ Newton's method with one iteration)
        E = M
        for i in range(n):
            E = E - (E-pl['e']*np.sin(E)-M)/(1-pl['e']*np.cos(E))
            
        # Check that the equation is properly solved
        sinE = np.sin(E)
        if (np.abs(M-(E-pl['e']*sinE)) > (0.002*np.pi)).any():
            print("Kepler's equation failed to solve! (e = {:.5f})".format(pl['e']))
        
        # Distance
        cosE = np.cos(E)
        r = pl['a']*(1-pl['e']*cosE)
        
        # True anomaly
        cos_nu = (cosE-pl['e'])/(1-pl['e']*cosE)
        sin_nu = ((1-pl['e']**2)**0.5*sinE)/(1-pl['e']*cosE)
        
    # Circular orbit
    else:
        nu = M
        r = pl['a']
        sin_nu = np.sin(M)
        cos_nu = np.cos(M)
    
    # Compute some intermediate terms
    cos_w,sin_w = np.cos(pl['w_AP']),np.sin(pl['w_AP'])
    cos_w_nu = cos_nu*cos_w-sin_nu*sin_w
    sin_w_nu = sin_nu*cos_w+sin_w*cos_nu
        
    # Compute x,y,z. z is the direction towards/from the observer.
    x = r*(np.cos(pl['w_LAN'])*cos_w_nu-np.sin(pl['w_LAN'])*sin_w_nu*pl['cos(i)'])
    y = r*(np.sin(pl['w_LAN'])*cos_w_nu+np.cos(pl['w_LAN'])*sin_w_nu*pl['cos(i)'])
    z = r*(np.sin(np.arccos(pl['cos(i)']))*sin_w_nu)
    
    return x,y,z

