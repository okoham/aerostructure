'''
Created on 14.11.2013

@author: Oliver
'''
import unittest
import pygeometry
import numpy as np
#import math


class TestPolygonArea2D(unittest.TestCase):


    def test_ccw(self):
        poly = np.array([(0,0), (1,0), (1,1), (0,1)], dtype=float)
        res = pygeometry.polygon_area_2d(poly)
        self.assertAlmostEqual(res, 1.)

    def test_cw(self):
        poly = np.array([(0,0), (0,1), (1,1), (1,0)], dtype=float)
        res = pygeometry.polygon_area_2d(poly)
        self.assertAlmostEqual(res, -1.)
        
    def test_2_vertices(self):
        poly = np.array([(0,0), (0,1)], dtype=float)
        self.assertRaises(AssertionError, pygeometry.polygon_area_2d, poly)
        
    def test_3_coincident_vertices(self):
        poly = np.array([(0,0), (0,1), (0,1)], dtype=float)
        res = pygeometry.polygon_area_2d(poly)
        self.assertAlmostEqual(res, 0.)
        #self.assertRaises(AssertionError, pygeometry.polygon_area_2d, poly)        
        
    def test_nonsimple(self):
        # this is a polygon with self-intersecting edges.
        # function just takes it as it is...
        poly = np.array([(0,0), (1,0), (0,1), (1,1)], dtype=float)
        res = pygeometry.polygon_area_2d(poly)
        self.assertAlmostEqual(res, 0.)    


class TestSegmentContainsPoint2D(unittest.TestCase):

    p00 = np.array([0.0, 0.0])
    p20 = np.array([2.0, 0.0])
    
    p10 = np.array([1.0, 0.0])
    p15 = np.array([0.5, 0.5])
    p_30 = np.array([-3.0, 0.0])
    p30 = np.array([3.0, 0.0])
    p31 = np.array([2.0, 1.0])

    
    def test_inside(self):
        u, v = pygeometry.segment_contains_point_2d(self.p00, self.p20, self.p10)
        self.assertTrue((0.0 <= u <= 1.0) and (v == 0.0), '?')
        
    def test_outside_1(self):
        u, v = pygeometry.segment_contains_point_2d(self.p00, self.p20, self.p_30)
        self.assertTrue((u < 0.) and (v == 0.0), 'u={}, v={}'.format(u, v))
   
    def test_outside_2(self):
        u, v = pygeometry.segment_contains_point_2d(self.p00, self.p20, self.p30)
        self.assertTrue((u > 1.) and (v == 0.0), 'u={}, v={}'.format(u, v))
        
    def test_left_end(self):
        u, v = pygeometry.segment_contains_point_2d(self.p00, self.p20, self.p00)
        self.assertTrue((u == 0.0) and (v == 0.0), 'u={}, v={}'.format(u, v))                

    def test_right_end(self):
        u, v = pygeometry.segment_contains_point_2d(self.p00, self.p20, self.p20)
        self.assertTrue((u == 1.0) and (v == 0.0), 'u={}, v={}'.format(u, v)) 
        
    def test_offset_1(self):
        u, v = pygeometry.segment_contains_point_2d(self.p00, self.p20, self.p15)
        self.assertTrue(v != 0.0, 'u={}, v={}'.format(u, v))
        
    def test_coincident(self):
        u, v = pygeometry.segment_contains_point_2d(self.p00, self.p00, self.p15)
        self.assertTrue(v == np.inf, 'u={}, v={}'.format(u, v))
        
    def test_coincident_2(self):
        u, v = pygeometry.segment_contains_point_2d(self.p00, self.p00, self.p00)
        self.assertTrue((0.0 <= u <= 1.0) and (v == 0.0), 'u={}, v={}'.format(u, v))
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()