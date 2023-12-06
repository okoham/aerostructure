import os
import inspect

import numpy as np
# import matplotlib.pyplot as plt 
import pandas as pd
import scipy.interpolate as spi


"""
Flexible circular frames supported by a shell. Moments, forces and displacements 
due to concentrated loads and couples (Data item 83043)

Flexible circular frames supporting a shell. The effect of adjacent frames and the 
longitudinal flexibility of the shell (data item 03.06.17)

"""


# Read coefficients CMM, CMT, ... from table xxx

# ESDU data contain some errors:
# - Error in table for d=2000, theta=15: CTN should be 0.98 instead of 0.098 (corrected)
# - Error in table for d=300, theta=45: CMM should be 0.018 in stead of -0.018 (corrected)
# - Error in table for d=4000, theta=15, CNM: 1.04 instead of 1.4 (corrected)
# - CQT, d=4000, theta=50 (also CηT): -0.09 instead of -0.9 (corrected)
# - CQT, d=4000, theta=45 (also CηT): -0.02 instead of -0.2 (corrected) 

# get the path of the directory where this module and the csv file live
_path = inspect.getfile(inspect.currentframe())
_dirname = os.path.dirname(_path)
coeff_table = pd.read_csv(os.path.join(_dirname, "esdu83043.csv"))

# values at which coefficient values are given
dvals = [0, 10, 50, 100, 300, 1000, 2000, 4000, 8000, 25000]
theta_deg = np.linspace(0, 180, num=37)
theta_rad = np.radians(theta_deg)

QUANTITIES = ["M", "T", "N", "δ", "ψ", "η", "Q"]
NQ = len(QUANTITIES)



# Set up bivariate interpolation. 
# The interpolator functions take two values: theta (usually array), and d (usually scalar)
# Note that the interpolators return 2D arrays, and that input values need to be monotonically
# increasing
_interpolators = {}
curve_ids = [f"C{x}{y}" for x in 'MTNQδψ' for y in 'MTN']
for curve_id in curve_ids:
    values = coeff_table.pivot(index="theta", columns=['d'], values=curve_id).values
    _interpolators[curve_id] = spi.RectBivariateSpline(theta_rad, dvals, values)
# the coefficients for eta are the same as for Q:
for k in "MTN":
    _interpolators["Cη" + k] = _interpolators["CQ" + k]


# Coefficients are defined for 0 ... 180 deg.
# if theta > 180: return sign * func(2pi - theta)
_signs = {
    'CMM': -1, 'CTM': -1, 'CNM': +1, 'CδM': -1, 'CψM': +1, 'CQM': +1, 'CηM': +1, 
    'CMT': -1, 'CTT': -1, 'CNT': +1, 'CδT': -1, 'CψT': +1, 'CQT': +1, 'CηT': +1,
    'CMN': +1, 'CTN': +1, 'CNN': -1, 'CδN': +1, 'CψN': -1, 'CQN': -1, 'CηN': -1
}


def coefficients(theta, d, curve_id):
    """Return interpolated coefficients, as a 1d numpy array. 

    Parameters:
    theta    ... angle (in radians) between applied load and section of interest.
    d        ... effective frame to shell stiffness parameter
    curve_id ... one of CMM, CMT, ..."""

    # the interpolation requires arrays with monotonically inceasing values. 
    # The inputs are normally not lilke that. Therefore, don't use vectorised 
    # computation, but compute each value separately. Maybe not so fast, but 
    # works reliably. 
    # Then, make it an array, so that we can iterate over it
    theta = np.mod(theta, 2*np.pi)
    theta = np.atleast_1d(theta)

    func = _interpolators[curve_id]
    sign = _signs[curve_id]

    c = []    # empy list of coefficients
    for th in theta: 
        if th <= np.pi:
            c.append(func(th, d))
        else:
            c.append(sign*func(2*np.pi - th, d))
    
    return np.squeeze(c)


# given theta, M, N, T, return all the stuff for one applies load.
# for different applied loads, can easily add up stuff.
# alpha is the point where load is applied
# phi is always the same array, -180...180. phi = 0 is top centre line
# theta is the angle of section under consideration wrt. the load application point A
# alpha is position of load application point

def compute_single_load(quantity, phi, alpha, Ma, Ta, Na, R, K, d=0):
    """
    Return quantities (one of M, T, N, Q, δ, ψ, η) for a single
    applied load, as a 1D array having the same length as phi.

    Parameters: 
    quantity  ... selects which quantity is computed. One of "M", "T", "N", "Q", 'δ', 'ψ', 'η'.
    phi   ... positions for computation (in radians, phi=0 is on top/at 12h, +ve CW). 
              float or list of float or 1d-array.
              In many cases, if forces/deformations for the complete frame is desired, 
              provide a complete array from 0 to 2pi, like np.linspace(-np.pi, np.pi, num=181)
    alpha ... angle where load is applied (in radians, from 12h position, +ve CW). 
    Ma    ... applied moment at position alpha (Nm, float). Positive moment creates tension in 
              inner flange.
    Ta    ... applied tangential force at position alpha (N, float). Positive force tends to 
              rotate frame CCW.
    Na    ... applied radial force at position alpha (N, float). Positive force outward.
    R     ... radius of frame (m, float)
    K     ... shell skin resisting force per unit tangential deflection (N/m, float).
              K = G t R / L'
    d     ... effective frame to shell stiffness parameter (float). 
              d = G t R^4 / (E I L') 
              d = 0 corresponds to a rigid frame, resulting in a sine shaped skin shear flow distribution.
    """

    assert quantity in QUANTITIES
    
    # theta ia angle between applied load (at alpha) and section of interest (at phi).
    theta = np.array(phi) - alpha

    match quantity:
        case 'M': 
            return coefficients(theta, d, 'CMM')*Ma          + coefficients(theta, d, 'CMT')*Ta*R     + coefficients(theta, d, 'CMN')*Na*R
        case 'T': 
            return coefficients(theta, d, 'CTM')*Ma/R        + coefficients(theta, d, 'CTT')*Ta       + coefficients(theta, d, 'CTN')*Na
        case 'N': 
            return coefficients(theta, d, 'CNM')*Ma/R        + coefficients(theta, d, 'CNT')*Ta       + coefficients(theta, d, 'CNN')*Na
        case 'δ': 
            return coefficients(theta, d, 'CδM')*Ma/(K*R)    + coefficients(theta, d, 'CδT')*Ta/K     + coefficients(theta, d, 'CδN')*Na/K
        case 'ψ': 
            return coefficients(theta, d, 'CψM')*Ma/(K*R**2) + coefficients(theta, d, 'CψT')*Ta/(K*R) + coefficients(theta, d, 'CψN')*Na/(K*R)
        case 'η': 
            return coefficients(theta, d, 'CηM')*Ma/(K*R)    + coefficients(theta, d, 'CηT')*Ta/K     + coefficients(theta, d, 'CηN')*Na/K
        case 'Q': 
            return coefficients(theta, d, 'CQM')*Ma/R**2     + coefficients(theta, d, 'CQT')*Ta/R     + coefficients(theta, d, 'CQN')*Na/R

    
def compute_combined_load(quantity, phi, applied_loads, R, K, d=0):
    """
    Return quantities (one of M, T, N, Q, δ, ψ, η) for a set of loads, 
    as a 1D array having the same length as phi.

    Parameters: 
    quantity  ... selects which quantity is computed. One of "M", "T", "N", "Q", 'δ', 'ψ', 'η'.
    phi   ... positions for computation (in radians, phi=0 is on top/at 12h, +ve CW). 
              float or list of float or 1d-array.
              In many cases, if forces/deformations for the complete frame is desired, 
              provide a complete array from 0 to 2pi, like np.linspace(-np.pi, np.pi, num=181)
    applied_loads ... list of dictionaries, one dictionary for each applied load. Dictionary keys are
              'alpha', 'Ma', 'Ta', 'Na'. See function `compute_single_load` for description.
    R     ... radius of frame (float)
    K     ... stiffness parameter (float)
    d     ... stiffness parameter (float). d=0 corresponds to a rigid frame, resulting in a 
              sine shaped skin shear flow distribution.
    """

    assert quantity in QUANTITIES
    
    total = np.zeros_like(phi)
    for f in applied_loads:
        # print("alpha, Ma, Ta, Na", alpha, Ma, Ta, Na)
        total += compute_single_load(quantity, phi, f['alpha'], f['Ma'], f['Ta'], f['Na'], R, K, d)   
    return total     


def compute_all(applied_loads, R, K, d=0, phi=np.linspace(-np.pi, np.pi, 361)):
    """Return a dataframe with all result quantities."""

    totals = pd.DataFrame(0, index=pd.Index(phi, name="phi_rad"), 
                          columns=['M', 'T', 'N', 'δ', 'ψ', 'η', 'Q'])
    # totals = pd.DataFrame(phi, columns=["phi_rad"])
    for q in QUANTITIES:
        totals[q] = compute_combined_load(q, phi, applied_loads, R, K, d) 

    return totals


def plot_loads(totals):
    """Return a matplotlib figure, containing all result plots.
    
    Parameter:
    totals: pandas Dataframe as returned from compute_all"""

    import matplotlib.pyplot as plt

    # TODO: plot applied loads, too!
    # TODO: make this customisable and a bit nicer.
    # TODO: what's load the reference point (on the section) for applied and section loads? Section centroid? Or OML? 
    
    #try:
    #    phi_rad = totals["phi_rad"]
    #except KeyError:
    phi_rad = totals.index
    phi_deg = np.degrees(phi_rad)
    
    plot_data = [
        ("M", "Bending moment M"),
        ("T", "Tangential force T"),
        ("N", "Shear force N"),
        ("Q", "Shear flow to skin Q"),
        ("δ", "radial derformation δ"),
        ("ψ", "ψ"),
        ("η", "η")
    ]

    fig, axes = plt.subplots(NQ, 1, figsize=(16,16), sharex=True)
    for i in range(NQ):
        quantity, ylabel = plot_data[i]
        axes[i].plot(phi_deg, totals[quantity])
        axes[i].set_ylabel(ylabel)
        if i == NQ-1: 
            axes[i].set_xlabel("phi (deg)")
        axes[i].grid()

    return fig