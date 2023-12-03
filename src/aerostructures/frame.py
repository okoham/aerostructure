import numpy as np
import pandas as pd


class FrameI(object):
    """
    Stress analyis of a (machined) I section frame, with doubler and skin. 
    The frame does not have stringer cutouts. 
    Frame, skin, doubler may be of different materials.
    The doubler may be used to represent a distinct doubler part, or it may represent a skin padup.
    

    Parameters: 
    w_sk ... skin width, or frame pitch
    t_sk ... skin thickness
    w_do ... doubler thickness
    t_do ... doubler thickness
    w_af ... width of attached flange
    t_af ... thickness of attached flange
    t_wb ... thickness of web
    w_ff ... width of free flange
    t_ff ... thickness of free flange
    h    ... total height of frame and skin
    mat_sk ... skin material (a dictionary)
    mat_do ... doubler material
    mat_fr ... frame material
    
    500, 't_sk': 2.0, 'mat_sk': alu, 
      'w_do': 50, 't_do': 1.0, 'mat_do': ti, 
      'w_af': 45, 't_af': 4.0, 
      't_wb': 2.0, 
      'w_ff': 45, 't_ff': 2.0,
      "h": 100, 'mat_fr': alu
    """ 
    
    def __init__(self, spec):

        self.spec = spec
        z0 = 0
        z1 = spec['t_sk']
        z2 = spec['t_sk'] + spec['t_do']
        z3 = spec['t_sk'] + spec['t_do'] + spec['t_af']
        z4 = spec['h'] - spec['t_ff']
        z5 = spec['h']

        materials = self.materials = {"sk": spec['mat_sk'], "fr": spec['mat_fr'], "do": spec['mat_do'], 
                                      "af": spec['mat_fr'], "wb": spec['mat_fr'], "ff": spec['mat_fr']}

        # a systematic list of stuff
        # columns: label, w, h, mat_name, Et, z00, z02
        geo = [
            ("sk", spec['w_sk'], spec['t_sk'], materials['sk']['name'], materials['sk']['et'], z0, z1),
            ("do", spec['w_do'], spec['t_do'], materials['do']['name'], materials['do']['et'], z1, z2),
            ("af", spec['w_af'], spec['t_af'], materials['fr']['name'], materials['fr']['et'], z2, z3),
            ("wb", spec['t_wb'], z4-z3,        materials['fr']['name'], materials['fr']['et'], z3, z4),
            ("ff", spec['w_ff'], spec['t_ff'], materials['fr']['name'], materials['fr']['et'], z4, z5),
        ]
        self.df = df = pd.DataFrame(geo, columns=['element', 'w', 'h', 'mat', 'et', 'zo0', 'zo2'])

        # TODO: drop doubler if not present??
        df['eid'] = [0, 1, 2, 3, 4]
        # df['elems'] = ["sk", "do", "af", "wb", "ff"]

        # all this is generic; can share this with other classes

        # EA, centroid
        df['zo1'] = (self.df["zo0"] + self.df["zo2"])/2
        df["A"] = self.df.w * self.df.h
        df["EA"] = self.df.et * self.df.A

        self.A = df.A.sum()
        self.EA = df.EA.sum()

        
        self.zcg = (self.df.zo1 * self.df.EA).sum() / self.EA
        print(f"zo_cg = {self.zcg:.1f}")

        self.df["zs0"] = self.df.zo0 - self.zcg
        self.df["zs1"] = self.df.zo1 - self.zcg
        self.df["zs2"] = self.df.zo2 - self.zcg

        # EI wrt 0
        df['eio'] = eio = df.et * (df.w * df.h**3/12 + df.zo1**2 * df.A)
        self.EI = eio.sum() - self.EA * self.zcg**2
        
        print(f"A = {df.A.sum():.1f}")
        print(f"EA = {self.EA:.1f}")
        print(f"EI = {self.EI:.1f}")

        # a dataframe for stresses
        # TODO: maybe an explicit python loop is easier to understand than this pandas stuff
        cols = ['element', 'et', 'eid', 'zs0', 'zs1', 'zs2']
        self.dfs = df[cols].melt(id_vars=["eid", "element", "et"], value_name="zs", var_name="pos").sort_values(["zs", 'eid']).reset_index(drop=True)
        self.dfs['zo'] = self.dfs['zs'] + self.zcg


    def stiffness_ratio(self):
        # stiffening ratio: frame stiffess (w/o doubler) / total stiffness
        ea_frame = 0 
        for _, row in self.df.iterrows(): 
            if row.element in ["ff", "wb", "af"]:
                ea_frame += row.EA
        return ea_frame / self.EA


    def mass(self):
        # mass per unit length
        mass = 0
        for i, row in self.df.iterrows():
            elem = row.element
            mass += self.materials[elem]['rho'] * row.A
        return mass

    
    def stresses(self, Fx, Fz, My, lcase=''):
        df = self.dfs.copy()
        et = df.et
        zs = df.zs
        # normal stresses
        df['sx'] = et * (Fx/self.EA + zs * My / self.EI)
        aweb = self.df.loc[self.df.element == "wb", "A"].values[0]
        # shear stress
        df['tau'] = 0.0
        df.loc[df.element == "wb", 'tau'] = Fz/aweb
        df['lcase'] = lcase
        return df
        

def strength_pt(ftu, reduction_factor=1):
    """Return allowable stress in tension, for flanges, doubler, skin.
    The reduction factor (a_net/a_gross) accounts for bolt holes etc."""
    return ftu * reduction_factor


def strength_pc(fcy, reduction_factor=1):
    """Return allowable plain compression stress, for flanges, doubler, skin.
    The reduction factor (a_net/a_gross) accounts for open holes etc."""
    return fcy * reduction_factor


def strength_crippling(fcy, fb):
    """Return crippling strength for free flange, using the estimation given 
    in ESDU 78020, eq. 2.1: f_c = \sqrt (c2 fb), where c2 is the compressive 
    yield strength and fb is the elastic flange buckling stress. 
    
    Parameters:
    mat ... Material (dict)
    fb ... flangre buckling stress
    """
    return np.sqrt(fcy * fb)

# def sigma_vm(sx, sy, sxy):
#     return np.sqrt(sx**2 + sy**3 )