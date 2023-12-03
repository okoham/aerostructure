import math

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