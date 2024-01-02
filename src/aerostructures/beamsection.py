# -*- coding: utf-8 -*-
'''
Created on 11.11.2013

@author: Oliver
'''

# from __future__ import print_function, division

#import geometry
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from numbers import Real
#import itertools
import pygeometry
from pathlib import Path

"""
simple polygon 

notation and definitions mostly in line with Wiedemann I

coordinate systems:
y0/z0 ... arbitrary
y/z ... origin at SP, axes parallel to y0/z0
ybar/zbar ... origin at SP, ybar is major, zbar minor principal axis

"""

TOL = 1.e-5
# TXT_TEMPLATE = 'c:/Users/Oliver/Projects/Bending/src/template_txt.txt'
# HTML_TEMPLATE = 'c:/Users/Oliver/Projects/Bending/src/template_html.txt'
TXT_TEMPLATE = Path(__file__).parent/'template_txt.txt'
HTML_TEMPLATE = Path(__file__).parent/'template_html.txt'

def transform_inertia(iy, iz, iyz, phi):
    # ix, iy, ixy are moments of inertia, wrt. section centroid
    # return moments of inertia, wrt. other axis, rotated phi ccw

    # wiedemann i, eq. (3.1-14)
    cp2 = math.cos(phi)**2
    sp2 = math.sin(phi)**2
    s2p = math.sin(2.*phi)
    c2p = math.cos(2.*phi)
    
    m = np.matrix([[cp2, sp2, -s2p], [sp2, cp2, s2p], [0.5*s2p, 0.5*s2p, c2p]])
    v = np.matrix([iy, iz, iyz]).T
    
    return m*v


def fig2base64(fig):
    """Return a png image, encoded as base64 string, of a matplotlib
    Figure instance"""
     
    from io import StringIO
    import base64
    
    imgdata = StringIO.StringIO()
    fig.savefig(imgdata, format='png')
    imgdata.seek(0)
    image_binary = imgdata.read()
    s64 = base64.b64encode(image_binary)
    imgdata.close()
    return s64



class PolygonSection(object):
    # notation: wiedemann I, fig. 3.1/1
    
    def __init__(self, vertices, title=""):
        
        self.vertices = np.array(vertices, dtype=float)
        assert len(self.vertices) >= 3, 'must have >= 3 vertices'
        
        self.title = str(title)
        
        # calculate area
        # check orientation
        self.area = pygeometry.polygon_area_2d(self.vertices)
        # if area is negative, polygon is CW
        assert self.area > 0., 'polygon has zero or negative area. either CW orientation or degenerate'
        
        # check self-intersection
        assert pygeometry.polygon_is_simple_2d(self.vertices), 'polygon has self-intersecting edges'

        # calculate centroid (relative to axes y0, z0)
        self.y0n, self.z0n = pygeometry.polygon_centroid_2d(self.vertices)
        
        # inertia xi, eta
        self.iy0 = pygeometry.polygon_yy_2d(self.vertices)
        self.iz0 = pygeometry.polygon_xx_2d(self.vertices)
        self.iyz0 = pygeometry.polygon_xy_2d(self.vertices)
        
        # inertia xi, eta
        self.iy = self.iy0 - self.area*self.z0n**2 
        self.iz = self.iz0 - self.area*self.y0n**2 
        self.iyz = self.iyz0 - self.area*self.y0n*self.z0n
        
        # inertia 1, 2
        t1 = 0.5*(self.iy + self.iz)
        t3 = math.sqrt(0.25*(self.iy - self.iz)**2 + self.iyz**2)
        self.i1 = t1 + t3
        self.i2 = t1 - t3
        
        # principal axis (1) angle in radians
        phi = 0.5 * math.atan2(2.0*self.iyz, self.iz - self.iy)
        ii1, dummy, dummy = transform_inertia(self.iy, self.iz, self.iyz, phi)
        if abs(ii1 - self.i1) < TOL:
            # phi belongs to i1
            self.phi1 = phi
            self.phi2 = phi + math.radians(90)
        elif abs(ii1 - self.i2) < TOL:
            # phi belongs to i2
            self.phi2 = phi
            self.phi1 = phi - math.radians(90)
        else:
            raise ValueError('Could not find angle for principal moment of inertia')
        
        # radii of gyration
        self.r1 = math.sqrt(self.i1/self.area)
        self.r2 = math.sqrt(self.i2/self.area)
        

    def section_plot(self):
        """Return a matplotlib.Figure instance showing the cross section."""
        # initialise figure
        fig, ax = plt.subplots()
        
        # plot the polygon
        poly = mpl.patches.Polygon(self.vertices, fc='none', ec='k')
        ax.add_patch(poly)
        
        # plot centroid
        ax.plot(self.y0n, self.z0n, 'ro', mec='r', mfc='none')
        
        # plot principal axes
        for r, phi, style in [(self.r1, self.phi1, 'b:'), (self.r2, self.phi2, 'g:')]:
            dy = r*math.cos(phi)
            dz = r*math.sin(phi)
            ax.plot([self.y0n-dy, self.y0n+dy], [self.z0n-dz, self.z0n+dz], style)
        
        # scaling etc.
        fig.set_figwidth(8)
        ax.set_ymargin(0.05)
        ax.set_xmargin(0.05)
        ax.set_aspect('equal')
        ax.autoscale_view()        
        return fig
    

    def vertex_coords_principal(self):
        """Return vertex coodinates in principal axes coordinates.
        numpy.array, shape (nvertices, 2)"""
        
        # vertex coordinates wrt. centroid (y, z):
        yz = self.vertex_coords_centroid()
        y = yz[:,0]
        z = yz[:,1]
        # vertex coords in principal axes system
        yzbar = np.zeros_like(yz, dtype=float)
        yzbar[:,0] =  y*np.cos(self.phi1) + z*np.sin(self.phi1)
        yzbar[:,1] = -y*np.sin(self.phi1) + z*np.cos(self.phi1)
        return yzbar
        
    
    def vertex_coords_centroid(self):
        """Return vertex coordinates relative to centroid, axes (y, z), 
           parallel to (y0, z0)."""
        return self.vertices - np.array([self.y0n, self.z0n])
    
    
    def report_html(self):
        """Generate an HTML report with cross section data. """
        from jinja2 import Template
        with open(HTML_TEMPLATE) as fin:
            s = fin.read()
        tmpl = Template(s)
        fig = self.section_plot()
        pic64 = fig2base64(fig)
        return tmpl.render(section=self, picture=pic64)
    
    
    def report_txt(self):
        from jinja2 import Template
        with open(TXT_TEMPLATE) as fin:
            s = fin.read()
        tmpl = Template(s)
        return tmpl.render(section=self)
            
    
    def report_json(self):
        d = {}
        for attr in ['area', 'y0n', 'z0n', 'iy0', 'iz0', 'iyz0', 'iy', 'iz', 'iyz',
                     'i1', 'i2', 'phi1', 'phi2', 'r1', 'r2', 'title']:
            d[attr] = getattr(self, attr)
        d['vertices'] = self.vertices.tolist()
        d['vertices_centroid'] = self.vertex_coords_centroid().tolist()
        d['vertices_principal'] = self.vertex_coords_principal().tolist()
        return d 
    
    
    def stresses(self, N, My0, Mz0, lap='centroid', lcase=''):
        """Normal stresses, due to loads applied at LAP. 
        LAP given in y0/z0 system.
        if LAP is not given, it is assumed at section centroid.
        loads given in y0/z0 coordinate system
        Other option: LAP at y0/z0 origin"""
        
        # myn, mzn are moments wrt. centroid n
        # case 1: loads given at section centroid
        
        assert isinstance(N, Real)
        assert isinstance(My0, Real)
        assert isinstance(Mz0, Real)
        
        if lap == "centroid":
            y0lap, z0lap = (self.y0n, self.z0n)
        elif lap == "origin":
            y0lap, z0lap = (0., 0.)
        elif isinstance(lap, (list, tuple, np.ndarray)):
            try:
                alap = np.array(lap, dtype=float)
            except:
                raise
            assert alap.shape == (2,), 'LAP must be shape (2,)'
            y0lap, z0lap = alap

        else:
            raise ValueError('incorrect format of lap. ')
        
        # moments wrt. centroid
        myn = My0 - N*(self.z0n - z0lap)
        mzn = Mz0 + N*(self.y0n - y0lap)
        
        co = math.cos(self.phi1)
        si = math.sin(self.phi1)

        # wiedemann i, eq. (3.1-4)
        mybar =  myn*co + mzn*si 
        mzbar = -myn*si + mzn*co
         
        # wiedemann i, eq. (3.1-21)
        sigma = []
        for i, yzbar in enumerate(self.vertex_coords_principal()):
            ybar, zbar = yzbar
            y, z = self.vertices[i]
            sigi = N/self.area
            sigi += zbar*mybar/self.i1
            sigi -= ybar*mzbar/self.i2
            sigma.append((i, y, z, sigi))
            
        dt = np.dtype([('vid', int), ('y0', float), ('z0', float), ('sig', float)])
        
        return {'N': N, 'My': My0, 'Mz': Mz0, 'LAP': (y0lap, z0lap), 
                'myn': myn, 'mzn': mzn, 'sigma': np.array(sigma, dtype=dt), 'lcase': lcase}
        

if __name__ == '__main__':
    pass