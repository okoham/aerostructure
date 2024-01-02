# -*- coding: utf-8 -*-
'''
Created on 14.11.2013

@author: Oliver
'''

# from __future__ import print_function, division
import numpy as np
import itertools
import math



# polygon is called simple if it is not self-intersecting. 
# this means that each edge e_i does not intersect any other edge, except for 
# the endpoints it shares with its adjacent edges.


def polygon_has_coincident_nodes_2d(verts):
    pass

def polygon_is_simple_2d(verts):
    # check if edges do not intersect except adjacent edges at common vertex
    # das funktioniert nicht, wenn linien aufeinander liegen
    
    is_simple = True
    tol = 1.e-5
        
    n = len(verts)
    ind = np.arange(n)
    ip1 = np.roll(ind, -1)
    im1 = np.roll(ind, 1)

    edges = list(zip(ind, ip1))
    edges = np.array(edges)
    
    dt = np.dtype([('e1', int), ('e2', int), ('x', float), ('y', float)])
    expect = np.zeros(n, dtype=dt)
    
    for i, j in zip(ind, im1):
        x, y = verts[i]
        iso, jso = sorted([i, j])
        expect[i] = (iso, jso, x, y)
    
    for ie1, ie2 in itertools.combinations(ind, 2):
        v1, v2 = edges[ie1]
        p1 = verts[v1]
        p2 = verts[v2]
        
        w1, w2 = edges[ie2]
        q1 = verts[w1]
        q2 = verts[w2]

        flag, r = segments_int_2d(p1, p2, q1, q2)

        if flag == 1:
            ieso1, ieso2 = sorted([ie1, ie2])
            msk  = expect['e1'] == ieso1
            msk &= expect['e2'] == ieso2
            #msk &= expect['x'] == r[0]
            #msk &= expect['y'] == r[1]
            msk &= np.abs(expect['x'] - r[0]) < tol
            msk &= np.abs(expect['y'] - r[1]) < tol
            is_expected = np.any(msk)
            
            if not is_expected: 
                print('intersecting edges:', ieso1, ieso2, 'at:', r)
                
            is_simple &= is_expected
            
    return is_simple


def segments_int_2d(p1, p2, q1, q2):

    #     !*****************************************************************************80
    #     !
    #     !! SEGMENTS_INT_2D computes the intersection of two line segments in 2D.
    #     !
    #     !  Discussion:
    #     !
    #     !    A line segment is the finite portion of a line that lies between
    #     !    two points P1 and P2.
    #     !
    #     !    In 2D, two line segments might not intersect, even though the
    #     !    lines, of which they are portions, intersect.
    #     !
    #     !  Licensing:
    #     !
    #     !    This code is distributed under the GNU LGPL license. 
    #     !
    #     !  Modified:
    #     !
    #     !    05 January 2005
    #     !
    #     !  Author:
    #     !
    #     !    John Burkardt
    #     !
    #     !  Parameters:
    #     !
    #     !    Input, real ( kind = 8 ) P1(2), P2(2), the endpoints of the first
    #     !    segment.
    #     !
    #     !    Input, real ( kind = 8 ) Q1(2), Q2(2), the endpoints of the second
    #     !    segment.
    #     !
    #     !    Output, integer ( kind = 4 ) FLAG, records the results.
    #     !    0, the line segments do not intersect.
    #     !    1, the line segments intersect.
    #     !
    #     !    Output, real ( kind = 8 ) R(2), an intersection point, if there is one.
    #     !

    tol = 0.001
    
    nointersect = np.nan * np.ones(2)
      
    # Find the intersection of the two lines.
    ival, r = lines_exp_int_2d(p1, p2, q1, q2)
    
    # ival
    # 0, no intersection, the lines may be parallel or degenerate.
    # 1, one intersection point, returned in P.
    # 2, infinitely many intersections, the lines are identical.    
    if ival == 0:
        flag = 0 # no intersection
        return flag, nointersect #r

    # Is the intersection point part of the first line segment?
    u, v = segment_contains_point_2d(p1, p2, r)
    
    if (u < 0.0) or (1.0 < u) or (tol < v):
        # outside
        flag = 0
        return flag, nointersect 

    # Is the intersection point part of the second line segment?
    u, v = segment_contains_point_2d(q1, q2, r)
    
    if (u < 0.0) or (1.0 < u) or (tol < v):
        # outside
        flag = 0
        return flag, nointersect #r
    
    # otherwise: point is on line
    flag = 1
    return flag, r


def segment_contains_point_2d(p1, p2, p):
    # !*****************************************************************************80
    # !
    # !! SEGMENT_CONTAINS_POINT_2D reports if a line segment contains a point in 2D.
    # !
    # !  Discussion:
    # !
    # !    A line segment is the finite portion of a line that lies between
    # !    two points P1 and P2.
    # !
    # !    In exact arithmetic, point P is on the line segment between
    # !    P1 and P2 if and only if 0 <= U <= 1 and V = 0.
    # !
    # !  Licensing:
    # !
    # !    This code is distributed under the GNU LGPL license. 
    # !
    # !  Modified:
    # !
    # !    17 August 2005
    # !
    # !  Author:
    # !
    # !    John Burkardt
    # !
    # !  Parameters:
    # !
    # !    Input, real ( kind = 8 ) P1(2), P2(2), the endpoints of the line segment.
    # !
    # !    Input, real ( kind = 8 ) P(2), a point to be tested.
    # !
    # !    Output, real ( kind = 8 ) U(2), the components of P, with the first
    # !    component measured along the axis with origin at P1 and unit at P2, 
    # !    and second component the magnitude of the off-axis portion of the
    # !    vector P-P1, measured in units of (P2-P1).
    # !

    p1 = np.array(p1, dtype=float)
    p2 = np.array(p2, dtype=float)
    p = np.array(p, dtype=float)
    
    assert p1.shape == (2,)
    assert p2.shape == (2,)
    assert p.shape == (2,)
    

    normsq = np.sum((p2-p1)**2)
    tol = 1.e-5 

    if normsq <= tol:
        if np.allclose(p, p1):
            u, v = 0.5, 0.0
        else:
            u, v = 0.5, np.inf
    else:
        #u = np.zeros(2, dtype=float)
        u = np.dot((p - p1), (p2 - p1)) / normsq
        v = np.sqrt( ((u - 1.)*p1[0] - u*p2[0] + p[0])**2 +
                     ((u - 1.)*p1[1] - u*p2[1] + p[1])**2) / np.sqrt(normsq)

    return np.array([u, v])


def lines_exp_int_2d(p1, p2, q1, q2):
    # !*****************************************************************************80
    # !
    # !! LINES_EXP_INT_2D determines where two explicit lines intersect in 2D.
    # !
    # !  Discussion:
    # !
    # !    The explicit form of a line in 2D is:
    # !
    # !      the line through the points P1 and P2.
    # !
    # !  Licensing:
    # !
    # !    This code is distributed under the GNU LGPL license. 
    # !
    # !  Modified:
    # !
    # !    02 January 2005
    # !
    # !  Author:
    # !
    # !    John Burkardt
    # !
    # !  Parameters:
    # !
    # !    Input, real ( kind = 8 ) P1(2), P2(2), two points on the first line.
    # !
    # !    Input, real ( kind = 8 ) Q1(2), Q2(2), two points on the second line.
    # !
    # !    Output, integer ( kind = 4 ) IVAL, reports on the intersection:
    # !    0, no intersection, the lines may be parallel or degenerate.
    # !    1, one intersection point, returned in P.
    # !    2, infinitely many intersections, the lines are identical.
    # !
    # !    Output, real ( kind = 8 ) P(2), if IVAL = 1, P is
    # !    the intersection point.  Otherwise, P = 0.
    # !

    tol = 1.e-5

    ival = 0
    p = np.zeros(2, dtype=float)

    # Check whether either line is a point.
    point_1 = np.allclose(p1, p2)
    point_2 = np.allclose(q1, q2)

    # Convert the lines to ABC format.
    if not point_1: 
        a1, b1, c1 = line_exp2imp_2d(p1, p2)

    if not point_2: 
        a2, b2, c2 = line_exp2imp_2d(q1, q2)

    # Search for intersection of the lines.
    if point_1 and point_2:
        if np.allclose(p1, q1):
            ival = 1
            p = p1
    elif point_1:
        if abs(a2 * p1[0] + b2 * p1[1] - c2) <= tol:
            ival = 1
            p = p1
    elif point_2:
        if abs(a1 * q1[0] + b1 * q1[1] - c1) <= tol:
            ival = 1
            p = q1
    else:
        ival, p = lines_imp_int_2d(a1, b1, c1, a2, b2, c2)

    return ival, p


def lines_imp_int_2d(a1, b1, c1, a2, b2, c2):
    # !*****************************************************************************80
    # !
    # !! LINES_IMP_INT_2D determines where two implicit lines intersect in 2D.
    # !
    # !  Discussion:
    # !
    # !    The implicit form of a line in 2D is:
    # !
    # !      A * X + B * Y + C = 0
    # !
    # !  Licensing:
    # !
    # !    This code is distributed under the GNU LGPL license. 
    # !
    # !  Modified:
    # !
    # !    25 February 2005
    # !
    # !  Author:
    # !
    # !    John Burkardt
    # !
    # !  Parameters:
    # !
    # !    Input, real ( kind = 8 ) A1, B1, C1, define the first line.
    # !    At least one of A1 and B1 must be nonzero.
    # !
    # !    Input, real ( kind = 8 ) A2, B2, C2, define the second line.
    # !    At least one of A2 and B2 must be nonzero.
    # !
    # !    Output, integer ( kind = 4 ) IVAL, reports on the intersection.
    # !
    # !    -1, both A1 and B1 were zero.
    # !    -2, both A2 and B2 were zero.
    # !     0, no intersection, the lines are parallel.
    # !     1, one intersection point, returned in P.
    # !     2, infinitely many intersections, the lines are identical.
    # !
    # !    Output, real ( kind = 8 ) P(2), if IVAL = 1, then P is
    # !    the intersection point.  Otherwise, P = 0.
    # !

    tol = 1.e-5

    p = np.zeros(2, dtype=float)

    # Refuse to handle degenerate lines.
    if line_imp_is_degenerate_2d(a1, b1, c1):
        ival = -1
        return ival, p

    if line_imp_is_degenerate_2d(a2, b2, c2):
        ival = -2
        return ival, p
    
    # Set up and solve a linear system.
    a = np.array([[a1, b1], [a2, b2]])
    rhs = np.array([-c1, -c2])
      
    try:
        # If the inverse exists, then the lines intersect at the solution point.
        p = np.linalg.solve(a, rhs)
        ival = 1
    except np.linalg.LinAlgError:
        # If the inverse does not exist, then the lines are parallel
        # or coincident.  Check for parallelism by seeing if the
        # C entries are in the same ratio as the A or B entries.
        ival = 0

        if abs(a1) <= tol:
            if abs(b2*c1 - c2*b1) <= tol:
                ival = 2
        else:
            if abs(a2*c1 - c2*a1) <= tol:
                ival = 2


    return ival, p


def line_imp_is_degenerate_2d(a, b, c):
    # !*****************************************************************************80
    # !
    # !! LINE_IMP_IS_DEGENERATE_2D finds if an implicit point is degenerate in 2D.
    # !
    # !  Discussion:
    # !
    # !    The implicit form of a line in 2D is:
    # !
    # !      A * X + B * Y + C = 0
    # !
    # !  Licensing:
    # !
    # !    This code is distributed under the GNU LGPL license. 
    # !
    # !  Modified:
    # !
    # !    06 May 2005
    # !
    # !  Author:
    # !
    # !    John Burkardt
    # !
    # !  Parameters:
    # !
    # !    Input, real ( kind = 8 ) A, B, C, the implicit line parameters.
    # !
    # !    Output, logical LINE_IMP_IS_DEGENERATE_2D, is true if the
    # !    line is degenerate.
    # !

    #return a*a + b*b == 0.0
    tol = 1.e-5 
    return a*a + b*b <= tol

def line_exp_point_dist_signed_2d(p1, p2, p):
    # !*****************************************************************************80
    # !
    # !! LINE_EXP_POINT_DIST_SIGNED_2D: signed distance ( exp line, point ) in 2D.
    # !
    # !  Discussion:
    # !
    # !    The explicit form of a line in 2D is:
    # !
    # !      the line through the points P1 and P2.
    # !
    # !    The signed distance has two interesting properties:
    # !
    # !    *  The absolute value of the signed distance is the
    # !        usual (Euclidean) distance.
    # !
    # !    *  Points with signed distance 0 lie on the line,
    # !       points with a negative signed distance lie on one side
    # !         of the line,
    # !       points with a positive signed distance lie on the
    # !         other side of the line.
    # !
    # !    Assuming that C is nonnegative, then if a point is a positive
    # !    distance away from the line, it is on the same side of the
    # !    line as the point (0,0), and if it is a negative distance
    # !    from the line, it is on the opposite side from (0,0).
    # !
    # !  Licensing:
    # !
    # !    This code is distributed under the GNU LGPL license. 
    # !
    # !  Modified:
    # !
    # !    06 May 2005
    # !
    # !  Author:
    # !
    # !    John Burkardt
    # !
    # !  Parameters:
    # !
    # !    Input, real ( kind = 8 ) P1(2), P2(2), two points on the line.
    # !
    # !    Input, real ( kind = 8 ) P(2), the point whose signed distance is desired.
    # !
    # !    Output, real ( kind = 8 ) DIST_SIGNED, the signed distance from the
    # !    point to the line.
    # !
    # point is left of line: negative
    assert isinstance(p1, np.ndarray)
    assert isinstance(p2, np.ndarray)
    assert isinstance(p, np.ndarray)
    assert p1.shape == (2,)
    assert p2.shape == (2,)
    assert p.shape == (2,)    
    
    # If the explicit line degenerates to a point, the computation is easy.
    if line_exp_is_degenerate_nd(p1, p2):
        dist_signed = np.linalg.norm(p1 - p) 
    # Convert the explicit line to the implicit form A * P(1) + B * P(2) + C = 0.
    # This makes the computation of the signed distance to (X,Y) easy.
    else:
        a = p2[1] - p1[1] 
        b = p1[0] - p2[0]
        c = p2[0] * p1[1] - p1[0] * p2[1]
        dist_signed = (a*p[0] + b*p[1] + c) / math.sqrt(a*a + b*b)

    return dist_signed



def line_exp2imp_2d(p1, p2):
    # !*****************************************************************************80
    # !
    # !! LINE_EXP2IMP_2D converts an explicit line to implicit form in 2D.
    # !
    # !  Discussion:
    # !
    # !    The explicit form of a line in 2D is:
    # !
    # !      the line through the points P1 and P2.
    # !
    # !    The implicit form of a line in 2D is:
    # !
    # !      A * X + B * Y + C = 0
    # !
    # !  Licensing:
    # !
    # !    This code is distributed under the GNU LGPL license. 
    # !
    # !  Modified:
    # !
    # !    06 May 2005
    # !
    # !  Author:
    # !
    # !    John Burkardt
    # !
    # !  Parameters:
    # !
    # !    Input, real ( kind = 8 ) P1(2), P2(2), two points on the line.
    # !
    # !    Output, real ( kind = 8 ) A, B, C, the implicit form of the line.
    # !

    
    # !  Take care of degenerate cases.
    if line_exp_is_degenerate_nd(p1, p2):
        print(' ') 
        print('LINE_EXP2IMP_2D - Warning!')
        print('  The line is degenerate.')

    a = p2[1] - p1[1]
    b = p1[0] - p2[0]
    c = p2[0] * p1[1] - p1[0] * p2[1]

    norm = a*a + b*b + c*c

    if 0.0 < norm:
        a /= norm
        b /= norm
        c /= norm

    if a < 0.0:
        a *= -1
        b *= -1
        c *= -1

    return a, b, c


def line_exp_is_degenerate_nd(p1, p2):
    # !*****************************************************************************80
    # !
    # !! LINE_EXP_IS_DEGENERATE_ND finds if an explicit line is degenerate in ND.
    # !
    # !  Discussion:
    # !
    # !    The explicit form of a line in ND is:
    # !
    # !      the line through the points P1 and P2.
    # !
    # !    An explicit line is degenerate if the two defining points are equal.
    # !
    # !  Licensing:
    # !
    # !    This code is distributed under the GNU LGPL license. 
    # !
    # !  Modified:
    # !
    # !    06 May 2005
    # !
    # !  Author:
    # !
    # !    John Burkardt
    # !
    # !  Parameters:
    # !
    # !    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
    # !
    # !    Input, real ( kind = 8 ) P1(DIM_NUM), P2(DIM_NUM), two points on the line.
    # !
    # !    Output, logical LINE_EXP_IS_DEGENERATE_ND, is TRUE if the line
    # !    is degenerate.
    # !

    assert p1.shape == p2.shape
    return np.allclose(p1, p2)



def polygon_area_2d(verts):

    """
    !*****************************************************************************80
    !
    !! POLYGON_AREA_2D computes the area of a polygon in 2D.
    !
    !  Discussion:
    !
    !    AREA = 1/2 * abs ( sum ( 1 <= I <= N ) X(I) * ( Y(I+1) - Y(I-1) ) )
    !    where Y(0) should be replaced by Y(N), and Y(N+1) by Y(1).
    !
    !    If the vertices are given in counter clockwise order, the area
    !    will be positive.  If the vertices are given in clockwise order,
    !    the area will be negative.
    !
    """
    assert isinstance(verts, np.ndarray), 'verts must be numpy.ndarray'
    assert verts.dtype == float, 'verts must be float'
    assert verts.ndim == 2, 'verts must be 2D'
    assert verts.shape[1] == 2, 'second axis must be dim 2'
           
    n = len(verts)
    assert n >= 3, 'Ploygon must have at least 3 vertices'
    
    i = np.arange(n)         # 0, 1, 2, 3
    im1 = np.roll(i, 1)      # 3, 0, 1, 2
    ip1 = np.roll(i, -1)     # 1, 2, 3, 0

    x = verts[:, 0]
    y = verts[:, 1]     

    dai = x[i] * (y[ip1] - y[im1])
    
    return 0.5*np.sum(dai)


def polygon_centroid_2d(verts):
    # !*****************************************************************************80
    # !
    # !! POLYGON_CENTROID_2D computes the centroid of a polygon in 2D.
    # !
    # !  Discussion:
    # !
    # !    Denoting the centroid coordinates by CENTROID, then
    # !
    # !      CENTROID(1) = Integral ( Polygon interior ) x dx dy / Area ( Polygon )
    # !      CENTROID(2) = Integral ( Polygon interior ) y dx dy / Area ( Polygon ).
    # !
    # !    Green's theorem states that for continuously differentiable functions
    # !    M(x,y) and N(x,y),
    # !
    # !      Integral ( Polygon boundary ) ( M dx + N dy ) =
    # !      Integral ( Polygon interior ) ( dN/dx - dM/dy ) dx dy.
    # !
    # !    Using M(x,y) = 0 and N(x,y) = x*x/2, we get:
    # !
    # !      CENTROID(1) = 0.5 * Integral ( Polygon boundary ) x*x dy 
    # !                  / Area ( Polygon ),
    # !
    # !    which becomes
    # !
    # !      CENTROID(1) = 1/6 sum ( 1 <= I <= N )
    # !        ( X(I+1) + X(I) ) * ( X(I) * Y(I+1) - X(I+1) * Y(I))
    # !        / Area ( Polygon )
    # !
    # !    where, when I = N, the index "I+1" is replaced by 1.
    # !
    # !    A similar calculation gives us a formula for CENTROID(2).
    # !
    # !  Licensing:
    # !
    # !    This code is distributed under the GNU LGPL license. 
    # !
    # !  Modified:
    # !
    # !    12 July 2003
    # !
    # !  Author:
    # !
    # !    John Burkardt
    # !
    # !  Reference:
    # !
    # !    Gerard Bashein, Paul Detmer,
    # !    Centroid of a Polygon,
    # !    in Graphics Gems IV, 
    # !    edited by Paul Heckbert,
    # !    AP Professional, 1994,
    # !    T385.G6974.
    # !
    # !  Parameters:
    # !
    # !    Input, integer ( kind = 4 ) N, the number of sides of the polygon.
    # !
    # !    Input, real ( kind = 8 ) V(2,N), the coordinates of the vertices.
    # !
    # !    Output, real ( kind = 8 ) CENTROID(2), the coordinates of the centroid.
    # !

    area = polygon_area_2d(verts)
    sx = polygon_y_2d(verts)
    sy = polygon_x_2d(verts)
    
    return np.array([sy, sx]) / area
    

def polygon_x_2d(verts):
    # !*****************************************************************************80
    # !
    # !! POLYGON_X_2D integrates the function X over a polygon in 2D.
    # !
    # !  Discussion:
    # !
    # !    The polygon is bounded by the points (X(1:N), Y(1:N)).
    # !
    # !    INTEGRAL = (1/6) * sum ( 1 <= I <= N )
    # !      ( X(I)*X(I) + X(I) * X(I-1) + X(I-1)*X(I-1) ) * ( Y(I) - Y(I-1) )
    # !
    # !    where X(0) and Y(0) should be replaced by X(N) and Y(N).
    # !
    # !  Licensing:
    # !
    # !    This code is distributed under the GNU LGPL license. 
    # !
    # !  Modified:
    # !
    # !    10 July 2001
    # !
    # !  Author:
    # !
    # !    John Burkardt
    # !
    # !  Reference:
    # !
    # !    SF Bockman,
    # !    Generalizing the Formula for Areas of Polygons to Moments,
    # !    American Mathematical Society Monthly,
    # !    1989, pages 131-132.
    # !
    # !  Parameters:
    # !
    # !    Input, integer ( kind = 4 ) N, the number of vertices of the polygon.
    # !    N should be at least 3 for a nonzero result.
    # !
    # !    Input, real ( kind = 8 ) V(2,N), the coordinates of the vertices
    # !    of the polygon.  These vertices should be given in counter clockwise order.
    # !
    # !    Output, real ( kind = 8 ) RESULT, the value of the integral.
    # !
    
    n = len(verts)
    assert n >= 3, 'POLYGON_X_2D - The number of vertices must be at least 3. is {}'.format(n)

    i = np.arange(n)
    im1 = np.roll(i, 1)
    
    x = verts[:, 0]
    y = verts[:, 1]    
        
    dfi = x[i]**2 + x[i]*x[im1] + x[im1]**2
    dfi *= y[i] - y[im1]

    return np.sum(dfi)/6.


def polygon_y_2d(verts):
    # !*****************************************************************************80
    # !
    # !! POLYGON_Y_2D integrates the function Y over a polygon in 2D.
    # !
    # !  Discussion:
    # !
    # !    The polygon is bounded by the points (X(1:N), Y(1:N)).
    # !
    # !    INTEGRAL = (1/6) * sum ( 1 <= I <= N )
    # !      - ( Y(I)^2 + Y(I) * Y(I-1) + Y(I-1)^2 ) * ( X(I) - X(I-1) )
    # !
    # !    where X(0) and Y(0) should be replaced by X(N) and Y(N).
    # !
    # !  Licensing:
    # !
    # !    This code is distributed under the GNU LGPL license. 
    # !
    # !  Modified:
    # !
    # !    10 July 2001
    # !
    # !  Author:
    # !
    # !    John Burkardt
    # !
    # !  Reference:
    # !
    # !    SF Bockman,
    # !    Generalizing the Formula for Areas of Polygons to Moments,
    # !    American Mathematical Society Monthly,
    # !    1989, pages 131-132.
    # !
    # !  Parameters:
    # !
    # !    Input, integer ( kind = 4 ) N, the number of vertices of the polygon.
    # !    N should be at least 3 for a nonzero result.
    # !
    # !    Input, real ( kind = 8 ) V(2,N), the coordinates of the vertices
    # !    of the polygon.  These vertices should be given in
    # !    counter clockwise order.
    # !
    # !    Output, real ( kind = 8 ) RESULT, the value of the integral.
    # !

    n = len(verts)
    assert n >= 3, 'POLYGON_Y_2D - The number of vertices must be at least 3. is {}'.format(n)

    i = np.arange(n)
    im1 = np.roll(i, 1)
    
    x = verts[:, 0]
    y = verts[:, 1]    

    dfi = y[i]**2 + y[i]*y[im1] + y[im1]**2
    dfi *= x[i] - x[im1]
    dfi *= -1
        
    return np.sum(dfi) / 6.0


def polygon_xx_2d(verts):
    # !*****************************************************************************80
    # !
    # !! POLYGON_XX_2D integrates the function X*X over a polygon in 2D.
    # !
    # !  Discussion:
    # !
    # !    The polygon is bounded by the points (X(1:N), Y(1:N)).
    # !
    # !    INTEGRAL = (1/12) * sum ( 1 <= I <= N )
    # !      ( X(I)^3 + X(I)^2 * X(I-1) + X(I) * X(I-1)^2 + X(I-1)^3 )
    # !      * ( Y(I) - Y(I-1) )
    # !
    # !    where X(0) and Y(0) should be replaced by X(N) and Y(N).
    # !
    # !  Licensing:
    # !
    # !    This code is distributed under the GNU LGPL license. 
    # !
    # !  Modified:
    # !
    # !    10 July 2001
    # !
    # !  Author:
    # !
    # !    John Burkardt
    # !
    # !  Reference:
    # !
    # !    SF Bockman,
    # !    Generalizing the Formula for Areas of Polygons to Moments,
    # !    American Mathematical Society Monthly,
    # !    1989, pages 131-132.
    # !
    # !  Parameters:
    # !
    # !    Input, integer ( kind = 4 ) N, the number of vertices of the polygon.
    # !    N should be at least 3 for a nonzero result.
    # !
    # !    Input, real ( kind = 8 ) V(2,N), the coordinates of the vertices
    # !    of the polygon.  These vertices should be given in
    # !    counter clockwise order.
    # !
    # !    Output, real ( kind = 8 ) RESULT, the value of the integral.
    # !


    n = len(verts)
    assert n >= 3, 'POLYGON_XX_2D - The number of vertices must be at least 3. is {}'.format(n)

    i = np.arange(n)
    im1 = np.roll(i, 1)
    
    x = verts[:, 0]
    y = verts[:, 1]
        
    dfi = x[i]**3 + x[i]**2*x[im1] + x[i]*x[im1]**2 + x[im1]**3
    dfi *= y[i] - y[im1]

    return np.sum(dfi)/12.


def polygon_xy_2d(verts):
    # !*****************************************************************************80
    # !
    # !! POLYGON_XY_2D integrates the function X*Y over a polygon in 2D.
    # !
    # !  Discussion:
    # !
    # !    The polygon is bounded by the points (X(1:N), Y(1:N)).
    # !
    # !    INTEGRAL = (1/24) * sum ( 1 <= I <= N )
    # !      ( Y(I)   * ( 3 * X(I)^2 + 2 * X(I) * X(I-1) +     X(I-1)^2 )
    # !      + Y(I-1) * (     X(I)^2 + 2 * X(I) * X(I-1) + 3 * X(I-1)^2 ) )
    # !      * ( Y(I) - Y(I-1) )
    # !
    # !    where X(0) and Y(0) should be replaced by X(N) and Y(N).
    # !
    # !  Licensing:
    # !
    # !    This code is distributed under the GNU LGPL license. 
    # !
    # !  Modified:
    # !
    # !    10 July 2001
    # !
    # !  Author:
    # !
    # !    John Burkardt
    # !
    # !  Reference:
    # !
    # !    SF Bockman,
    # !    Generalizing the Formula for Areas of Polygons to Moments,
    # !    American Mathematical Society Monthly,
    # !    1989, pages 131-132.
    # !
    # !  Parameters:
    # !
    # !    Input, integer ( kind = 4 ) N, the number of vertices of the polygon.
    # !    N should be at least 3 for a nonzero result.
    # !
    # !    Input, real ( kind = 8 ) V(2,N), the coordinates of the vertices
    # !    of the polygon.  These vertices should be given in
    # !    counter clockwise order.
    # !
    # !    Output, real ( kind = 8 ) RESULT, the value of the integral.
    # !
    
    # verts must have shape (n,2)
    assert isinstance(verts, np.ndarray), 'verts must be numpy.ndarray'
    assert verts.dtype == float, 'verts must be float'
    assert verts.ndim == 2, 'verts must be 2D'
    assert verts.shape[1] == 2, 'second axis must be dim 2'
    
    n = len(verts)
    assert n >= 3, 'POLYGON_xy_2D - The number of vertices must be at least 3. is {}'.format(n)

    i = np.arange(n)
    im1 = np.roll(i, 1)
    
    x = verts[:, 0]
    y = verts[:, 1]
    
    dfi =  y[i]   * (3*x[i]**2 + 2*x[i]*x[im1] +   x[im1]**2)
    dfi += y[im1] * (  x[i]**2 + 2*x[i]*x[im1] + 3*x[im1]**2)
    dfi *= y[i] - y[im1]
    
    return np.sum(dfi) / 24.0
    

def polygon_yy_2d(verts):
    # !*****************************************************************************80
    # !
    # !! POLYGON_YY_2D integrates the function Y*Y over a polygon in 2D.
    # !
    # !  Discussion:
    # !
    # !    The polygon is bounded by the points (X(1:N), Y(1:N)).
    # !
    # !    INTEGRAL = (1/12) * sum ( 1 <= I <= N )
    # !      - ( Y(I)^3 + Y(I)^2 * Y(I-1) + Y(I) * Y(I-1)^2 + Y(I-1)^3 )
    # !      * ( X(I) - X(I-1) )
    # !
    # !    where X(0) and Y(0) should be replaced by X(N) and Y(N).
    # !
    # !  Licensing:
    # !
    # !    This code is distributed under the GNU LGPL license. 
    # !
    # !  Modified:
    # !
    # !    10 July 2001
    # !
    # !  Author:
    # !
    # !    John Burkardt
    # !
    # !  Reference:
    # !
    # !    SF Bockman,
    # !    Generalizing the Formula for Areas of Polygons to Moments,
    # !    American Mathematical Society Monthly,
    # !    1989, pages 131-132.
    # !
    # !  Parameters:
    # !
    # !    Input, integer ( kind = 4 ) N, the number of vertices of the polygon.
    # !    N should be at least 3 for a nonzero result.
    # !
    # !    Input, real ( kind = 8 ) V(2,N), the coordinates of the vertices
    # !    of the polygon.  These vertices should be given in
    # !    counter clockwise order.
    # !
    # !    Output, real ( kind = 8 ) RESULT, the value of the integral.
    # !
    
    n = len(verts)
    assert n >= 3, 'POLYGON_YY_2D - The number of vertices must be at least 3. is {}'.format(n)

    i = np.arange(n)
    im1 = np.roll(i, 1)
    
    x = verts[:, 0]
    y = verts[:, 1]
   
    dfi = y[i]**3 + y[i]**2*y[im1] + y[i]*y[im1]**2 + y[im1]**3
    dfi *= -1
    dfi *= x[i] - x[im1]
    
    return np.sum(dfi)/12.


if __name__ == '__main__':
    pass