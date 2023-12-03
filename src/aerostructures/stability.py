import math
import numpy as np
import scipy.optimize as spo

def sigma_irb_hsb(pitch: float, thickness: float, ec: float, rc02: float, C=1.5) -> float:
    """
    Return inter-rivet buckling stress (in MPa) according to HSB 45131-01, Issue C, 2004.

    Parameters:
    pitch     ... fastener pitch in loading direction (s, mm)
    thickness ... sheet thickness (mm)
    C         ... clamping factor
                  C = 1.5: solid rivet, flush head or countersunk normal head
                  C = 3.0: solid rivet, normal head
                  C = 2.0: lockbolt (rivet, close tolerance), flush head
                  C = 4.0: lockbolt (rivet, close tolerance), protruding head
                  C = 1.0: blind rivet, flush head or countersunk normal head
                  C = 2.0: blind rivet, normal head
    Ec        ... Compressive modulus (MPa)
    Rc02      ... compressive yield strength (MPa)
    """

    assert pitch > 0
    assert thickness > 0
    assert ec > 0
    assert rc02 > 0
    assert 1 <= C <= 4 

    psi = 1.1027 * (pitch/thickness) * math.sqrt(rc02/(C * ec))
    if psi <= 1.5275: # inelastic
        return (1 - 0.3027*psi**1.5) * rc02
    else: # elastic regime
        return rc02/psi**2


def weff_bruhn(thickness, ec, fst, condition):
    """
    Return effective width (im mm) of postbuckled skin, attached to a stiffener 
    (stringer or frame). 

    Parameters: 

    """
    pass


def weff_simple(thickness):
    """20 t, on either side of rivet. 
    For aluminium, this can be derived from the Bruhn formula
    2 weff = 1.92 t \sqrt(Ec/Fcc), with E=70000 and Fcc =350 (assumung Fcc approx. Fcy ) 
    """
    return 20*thickness


def strength_flange_buckling(ec, thickness, width):
    """Return (elastic) compressive buckling stress of a free flange, according to 
    Bruhn, C6.1
    """
    return 0.388 * ec * (thickness/width)**2

# ----------------------------------

class Plate_SS_Ti(object):

    """
    Elastic buckling of simply supported rectangular plate, under combined bending and shear, 
    using the simple method of TI C 2512.
    """
    
    def __init__(self, E, nu, t, a, b):
        self.E = E
        self.nu = nu
        self.t = t
        self.a = a
        self.b = b
        self.alpha = ar = a/b

        # Euler reference stress, # TI, eq. 19
        self._sige = math.pi**2 * self.E * (self.t/self.b)**2 / (12. * (1 - self.nu**2))

        # compression in x dir, eq. 20
        if ar <= 1.:
            self.ksx = (ar + 1./ar)**2
        else:
            self.ksx = 4.
        # buckling stress, pure cpmpression
        self.sigc_cr_0 = self.ksx * self._sige

        # cisaillement, eq. 22
        if ar <= 1.:
            self.ktau = 4.0 + 5.34/ar**2
        else:
            self.ktau = 5.34 + 4.0/ar**2   
        # buckling stress, pure shear
        self.tau_cr_0 = self.ktau * self._sige
        
        # Flexion pure, eq. 21
        ar = self.alpha
        if ar <= 2./3.:
            self.kbx = 15.87 + 8.6*ar**2 + 1.87/ar**2
        else:
            self.kbx = 23.9      
        # buckling stress, pure bending
        self.sigb_cr_0 = self.kbx * self._sige


    def __str__(self):
        return f"Plate_SS_Ti: sigc_cr_0 = {self.sigc_cr_0:.1f}, sigb_cr_0 = {self.sigb_cr_0:.1f}, tau_cr_0 = {self.tau_cr_0:.1f}"


    def sig_cr_0(self):
        """Return vector of critical stresses, if each component is acting alone: 
        [sigc_cr_0, sigb_cr_0, tau_cr_0]"""
        return np.array([self.sigc_cr_0, self.sigb_cr_0, self.tau_cr_0])

    
    @staticmethod
    def _residual(x, sc, sb, tau, sccr0, sbcr0, tcr0):
        """
        Return the residual for solving the interaction equation. See TI C2512, Eq. 25 and 27.
    
        Parameters:
        x                  ... ratio of citical stress to applied stress (RF combined).
        sc, sb, tau        ... applied stresses in compression, bending, shear. NOTE: compression is positive!
        sccr0, sbcr0, tcr0 ... critical stresses if each component is applied alone
        """
        sccr, sbcr, tcr = x*sc, x*sb, x*tau
        return sccr/sccr0 + (sbcr/sbcr0)**2 + (tcr/tcr0)**2 - 1


    def rf(self, sig_applied, debug=False):
        """
        Return buckling RF, ratio of critical stresses to applied stresses, for 
        applied stressses acting together. 
        Returns np.inf if a pure tensile load is given.

        Parameters:
        sig_applied: [sigc, sigb, tau] ... applied stresses (NOTE sigc: COMPRESSION POSITIVE)
        """

        # unpack applied stresses
        sc, sb, tau = sig_applied

        # if sc negative, and sb = tau = 0, then the function raises an error. This is ok, because
        # it's pure tension and there is not buckling in that case. For that case, np.inf is returned.
        TOL = 1e-3
        if (sc <= 0) and (abs(sb) <= TOL) and (abs(tau) <= TOL):
            return np.inf

        # TODO: find a good bracket. 
        res = spo.root_scalar(Plate_SS_Ti._residual, method='bisect', x0=0, bracket=[0, 1000],
                       args=(sc, sb, tau, self.sigc_cr_0, self.sigb_cr_0, self.tau_cr_0))
        
        if debug: 
            print(res)
        return res.root