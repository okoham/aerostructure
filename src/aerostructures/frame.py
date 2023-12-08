import numpy as np
import pandas as pd
from stability import weff_simple, sigma_irb_hsb, strength_flange_buckling, Plate_SS_Ti


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

        # TODO: compute different cross section properties
        # 1. gross section (for mass, ...)
        # 2. net section with skin fully effective: strength analysis, skin in tension
        # 3. net section with reduced effective skin (strength in compression) 

        # EA, centroid
        df['zo1'] = (self.df["zo0"] + self.df["zo2"])/2
        df["A"] = self.df.w * self.df.h
        df["EA"] = self.df.et * self.df.A

        self.A = df.A.sum()
        self.EA = df.EA.sum()

        
        self.zcg = (self.df.zo1 * self.df.EA).sum() / self.EA
        # print(f"zo_cg = {self.zcg:.1f}")

        self.df["zs0"] = self.df.zo0 - self.zcg
        self.df["zs1"] = self.df.zo1 - self.zcg
        self.df["zs2"] = self.df.zo2 - self.zcg

        # EI wrt 0
        df['eio'] = eio = df.et * (df.w * df.h**3/12 + df.zo1**2 * df.A)
        self.EI = eio.sum() - self.EA * self.zcg**2
        
        # print(f"A = {df.A.sum():.1f}")
        # print(f"EA = {self.EA:.1f}")
        # print(f"EI = {self.EI:.1f}")

        # a dataframe for stresses
        # TODO: maybe an explicit python loop is easier to understand than this pandas stuff???
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


    def mass_total(self):
        # mass per unit length, including skin and doubler
        mass = 0
        for i, row in self.df.iterrows():
            elem = row.element
            mass += self.materials[elem]['rho'] * row.A
        return mass

    
    def mass_frame(self):
        # mass per unit length of frame alone, without skin/doubler
        mass = 0
        for i, row in self.df.iterrows():
            elem = row.element
            if elem in ('ff', 'wb', 'af'):
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


    def static_strength_analysis(self, M, T, N, Q):
        """ 
        M ... bending moment in frame at section under consideration, positive if it
              produces tension on inside of frame
        T ... direct force (tangential) in frame, positive if it produces tension
        N ... shear force (radial) in frame, positive if acting outward on clockwise
              side of section under consideration
        Q ... shear flow from shell to frame, positive if tending to rotate frame
              clockwise
        """

        res = {}
        # res.update(self.spec)

        res['unit_mass'] = self.mass()

        stresses_t = self.stresses(T, N, M)

        # sort out if skin is tension or compression
        sx_skin = stresses_t.query("element == 'sk' & pos == 'zs1'")['sx'].values[0]
        if sx_skin < 0: # compression, use effective skin width
            compr_spec = {k: v for k, v in self.spec.items()}
            # weff: that's per side of the frame! 
            # TODO: add rivet pitch  
            compr_spec['w_sk'] = weff = 2*weff_simple(self.spec['t_sk'])
            section = FrameI(compr_spec)
            res['weff'] = weff
            res['sk_compression'] = True
            stresses_c = section.stresses(T, N, M)
            stresses = stresses_c
        else: 
            res['weff'] = self.spec['w_sk']
            res['sk_compression'] = False
            section = self
            stresses = stresses_t

        # TODO: add stresses to result

        # skin plain tension/compression
        mat_sk = self.materials['sk']
        t_sk = self.spec['t_sk']
        sigma_sk = stresses_t.query("element == 'sk' & pos == 'zs1'")['sx'].values[0]
        #anet_agross = 1
        res['sigma_sk'] = sigma_sk
        res["R_sk_pt"] = R_sk_pt = strength_pt(mat_sk['ftu'])
        res["R_sk_pc"] = R_sk_pc = strength_pc(mat_sk['fcy'])
        # TODO: use actual fastener pitch and fixity
        res["R_sk_irb"] = R_sk_irb = sigma_irb_hsb(25, t_sk, mat_sk['ec'], mat_sk['fcy'], C=2.0)
        if sigma_sk > 0:
            res["rf_sk_pt"] = R_sk_pt / sigma_sk
            res["rf_sk_pc"] = np.nan 
        else: 
            res["rf_sk_pt"] = np.nan
            res["rf_sk_pc"] = R_sk_pc / abs(sigma_sk)
            res["rf_sk_irb"] = R_sk_irb / abs(sigma_sk)

        # doubler plain tension/compression
        # we don't do IRB here because it's clamped.
        mat_do = self.materials['do']
        t_do = self.spec['t_do']
        # TODO: Abzug für niete
        # anet_agross = 0.9
        res['sigma_do'] = sigma_do = stresses_t.query("element == 'do' & pos == 'zs1'")['sx'].values[0]
        res["R_do_pt"] = R_do_pt = strength_pt(mat_do['ftu'])
        res["R_do_pc"] = R_do_pc = strength_pc(mat_do['fcy'])
        if sigma_do > 0:
            res["rf_do_pt"] = R_do_pt / sigma_do
            res["rf_do_pc"] = np.nan 
        else: 
            res["rf_do_pt"] = np.nan
            res["rf_do_pc"] = R_do_pc / abs(sigma_do)

        # attached flange
        mat_fr = self.materials['fr']
        t_af = self.spec['t_af']
        w_af = self.spec['w_af']
        t_wb = self.spec['t_wb']
        res['sigma_af'] = sigma_af = stresses_t.query("element == 'af' & pos == 'zs1'")['sx'].values[0]
        # TODO: reduction factor
        res["R_af_pt"] = R_af_pt = strength_pt(mat_fr['ftu'])
        res["R_af_pc"] = R_af_pc = strength_pc(mat_fr['fcy'])
        res["R_af_lb"] = R_af_lb = strength_flange_buckling(mat_fr['ec'], t_af, (w_af - t_wb)/2)
        res["R_af_crp"] = R_af_crp = strength_crippling(mat_fr['fcy'], R_af_lb)
        if sigma_af > 0:
            res["rf_af_pt"] = R_af_pt / sigma_af
            res["rf_af_pc"] = np.nan 
            #res["rf_af_lb"] = np.nan 
            res["rf_af_crp"] = np.nan 
        else: 
            res["rf_af_pt"] = np.nan
            res["rf_af_pc"] = R_af_pc / abs(sigma_af)
            #res["rf_ff_lb"] = R_ff_lb / abs(sigma_ff)
            res["rf_af_crp"] = R_af_crp / abs(sigma_af) 

        
        # free flange
        mat_fr = self.materials['fr']
        t_ff = self.spec['t_ff']
        w_ff = self.spec['w_ff']
        res['sigma_ff'] = sigma_ff = stresses_t.query("element == 'ff' & pos == 'zs1'")['sx'].values[0]
        res["R_ff_pt"] = R_ff_pt = strength_pt(mat_fr['ftu'], reduction_factor=1)
        res["R_ff_pc"] = R_ff_pc = strength_pc(mat_fr['fcy'])
        res["R_ff_lb"] = R_ff_lb = strength_flange_buckling(mat_fr['ec'], t_ff, (w_ff - t_wb)/2)
        res["R_ff_crp"] = R_ff_crp = strength_crippling(mat_fr['fcy'], R_ff_lb)
        if sigma_ff > 0:
            res["rf_ff_pt"] = R_ff_pt / sigma_ff
            res["rf_ff_pc"] = np.nan 
            res["rf_ff_lb"] = np.nan 
            res["rf_ff_crp"] = np.nan 
        else: 
            res["rf_ff_pt"] = np.nan
            res["rf_ff_pc"] = R_ff_pc / abs(sigma_ff)
            res["rf_ff_lb"] = R_ff_lb / abs(sigma_ff)
            res["rf_ff_crp"] = R_ff_crp / abs(sigma_ff) 

        
        # web buckling
        # TODO: all edges simply supported is maybe a bit conservative?
        res["h_wb"] = h_wb = self.df.query("element == 'wb'")['h'].values[0]
        # TODO: give web stiffener distance as input. For the time being, we conservatively use a = 5b, 
        # which is basically not web stiffeners
        web_buckling = Plate_SS_Ti(mat_fr['ec'], mat_fr['nu'], t_wb, 5*h_wb, h_wb)
        res["wb_sc_cr_0"] = web_buckling.sigc_cr_0
        res["wb_sb_cr_0"] = web_buckling.sigb_cr_0
        res["wb_tau_cr_0"] = web_buckling.tau_cr_0
        smax = stresses.query("element == 'wb'").sx.max()
        smin = stresses.query("element == 'wb'").sx.min()
        res["tau_wb"] = tau = stresses.query("element == 'wb'").tau.abs().max() # they are all the same
        sx = (smax + smin)/2
        res["sc_wb"] = sc = - sx   # NOTE: compression +ve!
        res["sb_wb"] = sb = smax - sx
        res["rf_lb_wb"] = web_buckling.rf([sc, sb, tau])

        # web strength
        # TODO: reduction factor
        for _, row in stresses.query("element == 'wb'").iterrows():
            res[f'sig_wb_vm_{row.pos}'] = sig_vm = sigma_vm(row.sx, 0, row.tau)
            res[f'rf_wb_vm_{row.pos}'] = strength_pt(mat_fr['ftu']) / sig_vm

        # TODO: lateral stability

        # TODO: bolting to skin: flange failure

        # TODO: bolting to skin, skin failure

        # TODO: bolting to skin, bolt failure 

        
        # that's it!
        return res
        

def sigma_vm(sx, sy, tau): 
    """Return von Mises stress for plain stress state."""
    return np.sqrt(sx**2 + sy**2 - sx*sy + 3*tau**2)


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