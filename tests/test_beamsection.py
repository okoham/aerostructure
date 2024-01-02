# -*- coding: utf-8 -*-

'''
Created on 11.11.2013

@author: Oliver
'''
# from __future__ import print_function, division
import unittest
#from unittest.TestCase import assertAlmostEqual
import beamsection
import math
import numpy as np



class TestRectStress(unittest.TestCase):

    w = 2.
    h = 1. 
    
    flaeche = w*h
    xs = w/2.
    ys = h/2. 
    iy = w*h**3/12.
    wy = w*h**2/6.
    wz = h*w**2/6.
    
    rec = [(0, 0), (w, 0), (w, h), (0, h)]
    sec = beamsection.PolygonSection(rec)    
    
    def testStressNCentroid(self):
        force = 100.
        expect = np.ones(4) * force/(self.w*self.h)
        res = self.sec.stresses(force, 0, 0, lap='centroid')
        np.testing.assert_almost_equal(res['sigma']['sig'], expect)
    
    def testStressNOrigin(self):
        force = 100.
        expect = np.ones(4) * force/(self.w*self.h)
        res = self.sec.stresses(force, 0, 0, lap='centroid')
        np.testing.assert_almost_equal(res['sigma']['sig'], expect)
        
    def testStressMyCentroid(self):
        moment = 100.
        expect = np.array([-1, -1, 1, 1]) * moment/self.wy
        res = self.sec.stresses(0, moment, 0, lap='centroid')
        np.testing.assert_almost_equal(res['sigma']['sig'], expect)

    def testStressMzCentroid(self):
        moment = 100.
        expect = np.array([1, -1, -1, 1]) * moment/self.wz
        res = self.sec.stresses(0, 0, moment, lap='centroid')
        np.testing.assert_almost_equal(res['sigma']['sig'], expect)

class TestSimpleRect(unittest.TestCase):

    w = 200.
    h = 1. 
    
    flaeche = w*h
    xs = w/2.
    ys = h/2. 
    
    rec = [(0, 0), (w, 0), (w, h), (0, h)]
    sec = beamsection.PolygonSection(rec)

    def testArea(self):
        self.assertAlmostEqual(self.sec.area, self.flaeche)

    def testCentroidY(self):
        self.assertAlmostEqual(self.sec.y0n, self.xs)
        
    def testCentroidZ(self):
        self.assertAlmostEqual(self.sec.z0n, self.ys)

    def testInertiaOriginY(self):
        self.assertAlmostEqual(self.sec.iy0, self.w*self.h**3/12. + self.w*self.h*self.h**2/4.)
        
    def testInertiaOriginZ(self):
        self.assertAlmostEqual(self.sec.iz0, self.h*self.w**3/12. + self.w*self.h*self.w**2/4.)
        
    def testInertiaOriginYZ(self):
        self.assertAlmostEqual(self.sec.iyz0, (self.w*self.h)**2/4.)
        
    def testInertiaCentroidY(self):
        self.assertAlmostEqual(self.sec.iy, self.w*self.h**3/12.)
        
    def testInertiaCentroidZ(self):
        self.assertAlmostEqual(self.sec.iz, self.h*self.w**3/12.)
        
    def testInertiaCentroidYZ(self):
        self.assertAlmostEqual(self.sec.iyz, 0.)
        
    def testInertiaPrincipalMax(self):
        self.assertAlmostEqual(self.sec.i1, self.h*self.w**3/12.)

    def testInertiaPrincipalMin(self):
        self.assertAlmostEqual(self.sec.i2, self.w*self.h**3/12.)

    def testInertiaPrincipalAngle(self):
        self.assertAlmostEqual(self.sec.phi1, -math.pi/2.)


class TestTriangle(unittest.TestCase):

    # Roark's formulas for stress and strain, App. A, case 11
    a = 1.95
    b = 8.45 
    d = 6.1 

    rec = [(0, 0), (b, 0), (b-a, d)]
    sec = beamsection.PolygonSection(rec)

    def testArea(self):
        self.assertAlmostEqual(self.sec.area, 0.5*self.b*self.d)

    def testCentroidY(self):
        self.assertAlmostEqual(self.sec.y0n, (2*self.b - self.a)/3.)
        
    def testCentroidZ(self):
        self.assertAlmostEqual(self.sec.z0n, self.d - self.d*2./3.)
  
    def testInertiaCentroidY(self):
        self.assertAlmostEqual(self.sec.iy, self.b*self.d**3/36.)
        
    def testInertiaCentroidZ(self):
        expect = self.b**2 - self.a*self.b + self.a**2
        expect *= self.b*self.d/36.
        self.assertAlmostEqual(self.sec.iz, expect)
        
    def testInertiaCentroidYZ(self):
        expect = self.b - 2.*self.a
        expect *= self.b*self.d**2/72.
        self.assertAlmostEqual(self.sec.iyz, expect)

    def testInertiaPrincipalAngle(self):
        z = self.d*(self.b - 2.*self.a)
        n = self.b**2 - self.a*self.b + self.a**2 - self.d**2
        expect = 0.5*math.atan2(z, n)
        self.assertAlmostEqual(self.sec.phi2, expect)

class TestAngle(unittest.TestCase):

    # Roark's formulas for stress and strain, App. A, case 8
    b = 20.
    d = 37. 
    t = 2.6 
    flaeche = t*(b + d - t)
    
    xc = b**2 + d*t - t**2
    xc /= 2.*(b + d - t)
    
    yc = d**2 + b*t - t**2
    yc /= 2.*(b + d - t)
    
    ix = (b*d**3 - (b-t)*(d-t)**3)/3. - flaeche*(d-yc)**2
    iy = (d*b**3 - (d-t)*(b-t)**3)/3. - flaeche*(b-xc)**2
    ixy = (b**2*d**2 - (b-t)**2*(d-t)**2)/4. - flaeche*(b-xc)*(d-yc)

    rec = [(0, 0), (b, 0), (b, t), (t, t), (t, d), (0, d)]
    sec = beamsection.PolygonSection(rec)


    def testArea(self):
        self.assertAlmostEqual(self.sec.area, self.flaeche)

    def testCentroidY(self):
        self.assertAlmostEqual(self.sec.y0n, self.xc)
        
    def testCentroidZ(self):
        self.assertAlmostEqual(self.sec.z0n, self.yc)
        
    def testInertiaCentroidY(self):
        self.assertAlmostEqual(self.sec.iy, self.ix)
        
    def testInertiaCentroidZ(self):
        self.assertAlmostEqual(self.sec.iz, self.iy)
        
    def testInertiaCentroidYZ(self):
        self.assertAlmostEqual(self.sec.iyz, self.ixy)
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()