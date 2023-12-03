import math

"""
## Buckling strength of Metal Columns HSB 41100-01 E/2006

effective length factor C

Boundary conditions     C
--------------------------
cl - cl               0.5
cl - ss               0.7
cl - guided           1.0
ss - ss               1.0
cl - free             2.0

"""

class Column(object):
    def __init__(self, ec, length, c, imin, rc02, nroc, area):
        self.ec = ec
        self.length = length
        self.c = c
        self.imin = imin
        self.rc02 = rc02
        self.nroc = nroc
        self.area = area
        
        self.gyration_radius = math.sqrt(self.imin/self.area)
        self.leff = self.c*self.length
        
    def sigma_euler(self):
        # elastic buckling stress 
        return math.pi**2 * self.ec / (self.leff/self.gyration_radius)**2
        
    def sigma_cr(self):
        # inelastic buckling stress
        xi, yi = self.transition_point()
        x = self._x()
        if x < xi:
            # eq. 12
            y = 1. - (1.-yi)*(x/xi)
        else:
            tol = 1.e-5
            n = self.nroc
            f3 = self._f3()
            yold = 1.
            ynew = 0.
            while abs(ynew-yold) > tol:
                yold = ynew
                ynew = 1./( x**2 * (1. + n * (f3*yold)**(n-1.) ))
            y = ynew
        return y*self.sigma_0()
            
            
    
    def sigma_c_fail(self):
        # failure stress, min of inelastic buckling and rc02
        # eq. 15, 16
        sc = self.sigma_cr()
        if sc < self.rc02:
            return sc
        else:
            return self.rc02
            
    def fc_fail(self):
        return self.sigma_c_fail() * self.area
    

    def sigma_c07(self):
        # eq. 3
        return self.rc02 * (3*self.rc02/(7*0.002*self.ec))**(1./(self.nroc-1.))
    
    
    def _f1(self):
        n = self.nroc
        # eq. 5
        f1 = (n-4)*n + math.sqrt((n-4.)**2 * n**2 - 3.*(n+2.)*n)
        f1 /= (n+2.)*n**2  
        return f1
    
    def _f2(self):
        n = self.nroc
        # eq. 6
        f2 = 3. + n*(n+2.)*self._f1()
        f2 /= 1. + n**2 * self._f1()
        return f2
       
    def _f3(self):
        # eq. 8
        return self._f2() * self._f1()**(1./(self.nroc-1.))
        
    def s0(self):
        #
        if self.nroc >= 6.:
            return self._f3() / (3./7.)**(1./(self.nroc-1.))
        else:
            # eq. 8
            return 1.255
        
    def sigma_0(self):
        # eq. 4
        return self.s0() * self.sigma_c07()
    
    def _x(self):
        return math.sqrt(self.sigma_0()/self.ec) * self.leff / (self.gyration_radius*math.pi)
    
    def transition_point(self):
        # eq. 10
        si = ((3./7.)*self._f1())**(1./(1-self.nroc))
        ri = 1./math.sqrt(si + (3./7.)*self.nroc*si**self.nroc)
        s0v = self.s0()
        xi = ri/math.sqrt(s0v)
        yi = si/s0v
        return xi, yi
    

if __name__ == '__main__':

    c1 = Column(206900, 560, 1, 16760, 960, 22, 143.7)

    print(c1.gyration_radius)
    print(c1.leff)
    print(c1._x())
    print(c1.sigma_c07())
    #print c1.transition_point()
    #c1.sigma_cr()
    print(c1.sigma_cr())
    print(c1.sigma_c_fail(), c1.sigma_c_fail()/c1.rc02)
    print(c1.sigma_euler())
    print(c1.fc_fail())