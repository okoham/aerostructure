# -*- coding: utf-8 -*-
"""
Created on Sat Jan 12 22:09:08 2013

@author: oliver
"""

# 2D bolt group

import numpy as np
import json

dt_bolt = np.dtype([('boltid', np.str_, 10),
                    ('x', float),
                    ('y', float),
                    ('Fallow', float),
                    ('description', np.str_, 80)])


dt_boltload = np.dtype([('lcid', np.uint), # oder sollten wir einen str nehmen?
                        ('Fx', float),
                        ('Fy', float),
                        ('Mz', float),
                        ('lcname', np.str_, 80)])


dt_res = np.dtype([('boltid', np.str_, 10),
                   ('lcid', np.uint),
                   ('Fx', float),
                   ('Fy', float),
                   ('Fres', float),
                   ('Fallow', float),
                   ('RF', float)])



class BoltGroup2DResultObject(object):

    def __init__(self, geo, loads, result, label=''):
        self._geo = geo
        self._apploads = loads
        self._result = result
        self.label = label


    @property
    def results(self):
        return self._result


    @property
    def cg(self):
        return self._geo.schwerpunkt()


    @property
    def bolts(self):
        return self._geo.bolts


    @property
    def lap(self):
        return self._geo.u


    @property
    def applied_loads(self):
        return self._apploads


    def __len__(self):
        return len(self._result)


    def nth_critical_lc(self, n):
        pass


    def lowest_rfs(self, n=3):
        return np.sort(self._result, order='RF')[:n]


    def plot_geometry(self):
        import matplotlib.pyplot as plt
        # new axes
        f = plt.axes()
        # overall config
        f.set_title('Bolt Group Geometry')
        f.set_xlabel('x')
        f.set_ylabel('y')
        # plot bolts
        b = self._geo.bolts
        f.scatter(b['x'], b['y'], s=50, c='blue', marker='o', facecolors='none')
        for x, y, s, in zip(b['x'], b['y'], b['boltid']):
            f.text(x, y+2, s)
        # plot cg
        xs, ys = self._geo.schwerpunkt()
        f.scatter(xs, ys, s=50, c='blue', marker='+')
        f.text(xs, ys+2, 'CG')
        # plot load application point
        ux, uy = self._geo.u
        f.scatter(ux, uy, s=50, c='red', marker='s', facecolors='none')
        f.text(ux, uy+2, 'U')
        return f


    def plot_forces(self, lcase):
        import matplotlib.pyplot as plt
        results = self._result

        f = plt.axes()
        # overall config
        f.set_title('Bolt Forces, Case {0}'.format(lcase))
        f.set_xlabel('x')
        f.set_ylabel('y')

        d = np.sort(self._geo.bolts, order='boltid')
        f1 = results[results['lcid'] == lcase]
        f1 = np.sort(f1, order='boltid')

        f.scatter(d['x'], d['y'], s=50, c='blue', marker='o', facecolors='none')
        f.quiver(d['x'], d['y'], f1['Fx'], f1['Fy'])
        return f


    def html_results(self):
        return "<h4>BoltGroup2D Results</h4>Object: {0}<br>Calculation: {1}".format(self._geo.label, self.label)


    def _repr_html_(self):
        return self.html_results()


def geo2string(arr):
    """convert bolt definition array to string"""
    header = "    boltid          x          y     Fallow descrp              \n"
    underl = "---------- ---------- ---------- ---------- --------------------\n"
    t = "{0:>10} {1:>10.2f} {2:>10.2f} {3:>10.1f} {4:<20}\n"
    s = '\n' + header + underl
    for l in arr.flat:
        s += t.format(l['boltid'], l['x'], l['x'], l['Fallow'], l['description'])
    return s


def lcase2string(arr):
    """convert load case array to string"""
    header = "      lcid         Fx         Fy       Mz,u descrp              \n"
    underl = "---------- ---------- ---------- ---------- --------------------\n"
    t = "{0:>10} {1:>10.1f} {2:>10} {3:>10.1f} {4:<20}\n"
    s = '\n' + header + underl
    for l in arr.flat:
        s += t.format(l['lcid'], l['Fx'], l['Fy'], l['Mz'],l['lcname'])
    return s


def result2string(arr):
    """convert array of results to string"""
    # boltid, lcid, Fx, Fy, Fres, RF
    header = "    boltid     Fallow       lcid         Fx         Fy          F         RF\n"
    underl = "---------- ---------- ---------- ---------- ---------- ---------- ----------\n"
    t = "{0:>10} {1:>10.1f} {2:>10} {3:>10.1f} {4:>10.1f} {5:>10.1f} {6:>10.2f}\n"
    s = '\n' + header + underl
    for l in arr.flat:
        s += t.format(l['boltid'], l['Fallow'], l['lcid'], l['Fx'],l['Fy'], l['Fres'], l['RF'])
    return s


class BoltGroup2D(object):
    def __init__(self, bolts, lap, label=''):
        # feed me with a record array
        # each has properties x, y, Pall, name
        self.bolts = np.array(bolts, dtype=dt_bolt)
        self.u = np.array(lap, dtype=float) # load application point
        self.label = label


    def save_json(self, fpath):
        di = {}
        di['_type'] = 'BoltGroup2D'
        di['bolts'] = self.bolts.tolist()
        di['lap'] = self.u.tolist()
        di['label'] = self.label
        json.dump(di)
        

    # def plot3d(self):
    #     import visual
        
    #     scene1 = visual.display(title='Boltgroup', exit=False)
    #     scene1.background = visual.color.white
    #     scene1.forground = visual.color.blue
    #     scene1.select()
        
    #     bolt_length = 10.
    #     bolt_diameter = 5.
        
    #     straight = [(0, 0, 0), (0, 0, bolt_length)]
        
    #     for b in self.bolts:
    #         cb = visual.shapes.circle(pos=(b['x'], b['y']), radius=bolt_diameter/2.)
    #         visual.extrusion(pos=straight, shape=cb, color=visual.color.yellow)
    

    @property
    def nbolts(self):
        return len(self.bolts)


    def schwerpunkt(self):
        nenner = np.sum(self.bolts['Fallow'])
        xs = np.dot(self.bolts['x'], self.bolts['Fallow'])/ nenner
        ys = np.dot(self.bolts['y'], self.bolts['Fallow'])/ nenner
        return (xs, ys)


    def loads_wrt_sp(self, px, py, mzu):
        # moment auf SP umrechnen
        xu, yu = self.u
        xs, ys = self.schwerpunkt()
        mzs = mzu + (xu-xs)*py - (yu-ys)*px
        return mzs


    def analyse(self, apploads, label=''):
        apploads = np.array(apploads, dtype=dt_boltload)
        # weights for Fx, Fy
        sfall = np.sum(self.bolts['Fallow'])
        # weights for Mz
        xs, ys = self.schwerpunkt()
        nenner = (self.bolts['x'] - xs)**2 + (self.bolts['y'] - ys)**2
        nenner *= self.bolts['Fallow']
        nenner = np.sum(nenner)

        # moment auf SP umrechnen
        xu, yu = self.u
        px = apploads['Fx']
        py = apploads['Fy']
        mzs = self.loads_wrt_sp(px, py, apploads['Mz'])

        res = np.zeros((self.nbolts, len(apploads)), dtype=dt_res)

        for i in range(self.nbolts):
            # TODO: for schleife rausnehmen
            bolt = self.bolts[i]
            alpha = bolt['Fallow'] / sfall
            beta_x = bolt['Fallow'] * (bolt['y'] - ys) / nenner
            beta_y = bolt['Fallow'] * (bolt['x'] - xs) / nenner

            res[i,:]['boltid'] = bolt['boltid']
            res[i,:]['lcid'] = apploads['lcid']

            fx = px*alpha - mzs*beta_x
            fy = py*alpha + mzs*beta_y

            res[i,:]['Fx'] = fx
            res[i,:]['Fy'] = fy

            res[i,:]['Fres'] = np.sqrt(fx**2 + fy**2)
            res[i,:]['Fallow'] = bolt['Fallow']
        # finally, reserve factors
        res['RF'] = res['Fallow']/res['Fres']
        return BoltGroup2DResultObject(self, apploads, res.flatten(), label=label)


    def html_description(self):
        return "<h4>BoltGroup2D: {0}</h4>".format(self.label)


    def _repr_html_(self):
        return self.html_description()


def boltgroup_from_json(fpath):
    """Return bg2d object, from json file"""
    # should delete comments first
    d = json.load(fpath)


if __name__ == '__main__':
    pass