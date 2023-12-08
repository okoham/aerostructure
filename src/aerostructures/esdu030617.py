import numpy as np
import scipy.interpolate as spi

# we do interpolation in log space.
# Intermediate values (2e2, 2e3, ...) have been omitted.

# key: parameter p = R^6t'/IL, values: log10x = log10(GtR^4/EIL'), log10y = log10(GtR^4/EIL) 
CURVES = {
    1e2: (np.array([1.301, 1.712, 2.2152, 2.7875, 3.4827, 4.0358, 4.5351, 5.2957, 5.941]),
          np.array([1.3935, 1.6947, 1.8685, 1.9573, 1.9921, 1.9959, 1.9959, 1.9959, 1.9959])),
    1e3: (np.array([1.301, 1.7773, 2.1038, 2.6531, 3.1754, 3.8284, 4.3277, 4.7733, 5.4416, 5.987]),
          np.array([1.3935, 1.8183, 2.0345, 2.2971, 2.4245, 2.4863, 2.4979, 2.4979, 2.5018, 2.5056])),
    1e4: (np.array([1.301 , 1.5891, 1.9809, 2.3343, 2.6953, 3.1178, 3.5134, 3.9398, 4.3815, 4.8732, 5.4109, 5.987 ]),
          np.array([1.3935, 1.687 , 2.0732, 2.3898, 2.6524, 2.8802, 3.0269, 3.1273, 3.1698, 3.1891, 3.1968, 3.1968])),
    1e5: (np.array([1.301 , 1.785 , 2.1384, 2.5071, 2.8835, 3.2753, 3.7746, 4.2317, 4.7925, 5.384 , 5.9909]),
          np.array([1.3974, 1.8801, 2.2315, 2.5983, 2.9381, 3.2239, 3.4671, 3.5984, 3.6756, 3.7027, 3.7181])),
    1e6: (np.array([1.301 , 1.639 , 1.9847, 2.3458, 2.799 , 3.1908, 3.6978, 4.1318, 4.6427, 5.0806, 5.4954, 5.9947]),
          np.array([1.3974, 1.7372, 2.0809, 2.4361, 2.8918, 3.2625, 3.6486, 3.8957, 4.0695, 4.1545, 4.2008, 4.2355])),
    2e6: (np.array([1.301 , 1.7082, 2.1422, 2.4841, 2.9142, 3.4635, 3.9897, 4.4468, 4.8924, 5.3994, 5.9909]),
          np.array([1.3935, 1.8029, 2.2353, 2.5751, 3.0076, 3.5096, 3.8649, 4.0965, 4.2433, 4.3398, 4.3784]))
}

P = np.array([1e2, 1e3, 1e4, 1e5, 1e6, 2e6])
LOGP = np.log10(P)

# set up spline interpolators in log space
interpolators = {}
for p, (log10x, log10y) in CURVES.items():
    interpolators[p] = spi.CubicSpline(log10x, log10y)


def deff(d0, p):
    """Return the effective frame to shell relative stiffness parameter, defined by 
    deff = GtR^4/EIL' (see Data Item No. Struct. 03.06.17)

    G ... shear modulus for skin material
    E ... Young's modulus for frame material
    I ... effective moment of inertia of cross section of frame including an
          allowance for effect of attached skin
    R ... radius measured to centroid of frame section
    L ... distance between adjacent frames
    t ... thickness of shell skin effective in shear
    t' ... equivalent thickness of shell skin effective in carrying longitudinal stress

    Non-uniformity of frame spacing: L = 2 L1 L2 / (L1 + L2), where L1, L2 are the 
    distances to the nearest frames.

    See data item: There is something on non-uniformity of frame stiffness...

    Parameters:
    d0   = GtR^4/EIL
    p    = R^6t'/IL
    
    Returns: 
    deff ... effective skin-frame stiffness ratio, 
    
    """

    # compute one function value at d0 for each curve
    logdeff_i = np.array([func(np.log10(d0)) for func in interpolators.values()])
    # then do linear interpolation
    logdeff = np.interp(np.log10(p), LOGP, logdeff_i)
    return 10**logdeff


if __name__ == '__main__':
    pass