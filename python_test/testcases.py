from trefftzngs import *
import math
math.pi
from ngsolve import *



def simplesin(D, wavespeed):
    t = CoordCF(D)
    k = 1
    if(D==1):
        sol = CoefficientFunction((
            sin(math.pi*x)*sin(math.pi*t*wavespeed),
            cos(math.pi*x)*sin(math.pi*t*wavespeed)*math.pi,
            sin(math.pi*x)*cos(math.pi*t*wavespeed)*wavespeed*math.pi
            ))
        # sol[0] = sin( k*(wavespeed*t + x) );
        # sol[1] = k*cos(k*(wavespeed*t+x));
        # sol[2] = wavespeed*k*cos(k*(wavespeed*t+x));
        return sol
    elif(D==2):
        sq = sqrt(0.5);
        sol = CoefficientFunction((
            sin( wavespeed*t+sq*(x+y) ),
            sq*cos(wavespeed*t+sq*(x+y)),
            sq*cos(wavespeed*t+sq*(x+y)),
            wavespeed*cos(wavespeed*t+sq*(x+y))
            ))
        return sol
    elif(D==3):
        sq = sqrt(1.0/3.0);
        sol = CoefficientFunction((
            sin( wavespeed*t+sq*(x+y+z) ),
            sq*cos(wavespeed*t+sq*(x+y+z)),
            sq*cos(wavespeed*t+sq*(x+y+z)),
            sq*cos(wavespeed*t+sq*(x+y+z)),
            wavespeed*cos(wavespeed*t+sq*(x+y+z))
            ))
        return sol


def gausspw(D, wavespeed, x0 = [0.5,0.5,0.5]):
    t = CoordCF(D)
    k = 1
    delta = 1000
    if(D==1):
        g = exp( -delta*((x-x0[0])*(x-x0[0])) ),
        sol = CoefficientFunction((
            g,
            -2*delta * (x-x0[0]) * g,
            0
            ))
        return sol
    elif(D==2):
        g = exp(-delta*((x-x0[0])*(x-x0[0])+(y-x0[1])*(y-x0[1])) ),
        sol = CoefficientFunction((
            g,
            -2*delta * (x-x0[0]) * g,
            -2*delta * (y-x0[1]) * g,
            0
            ))
        return sol
    elif(D==3):
        g = exp(-delta*((x-x0[0])*(x-x0[0])+(y-x0[1])*(y-x0[1])+(z-x0[2])*(z-x0[2])) ),
        sol = CoefficientFunction((
            g,
            -2*delta * (x-x0[0]) * g,
            -2*delta * (y-x0[1]) * g,
            -2*delta * (z-x0[2]) * g,
            0
            ))
        return sol

def vertgausspw(D, wavespeed, x0=0.25):
    t = CoordCF(D)
    k = 1
    delta = 0.07
    if(D==1):
        sol = CoefficientFunction((
            exp( -(1/(delta*delta))*((x-x0)*(x-x0)) ),
            -2*(1/(delta*delta)) * (x-x0) * exp( -(1/(delta*delta))*((x-x0)*(x-x0)) ),
            0
            ))
        return sol
    elif(D==2):
        sol = CoefficientFunction((
            exp(-(1/(delta*delta))*((x-x0)*(x-x0)) ),
            -2*(1/(delta*delta))* (x-x0) * exp(-(1/(delta*delta))*((x-x0)*(x-x0)) ),
            0,
            0
            ))
        return sol
    elif(D==3):
        sol = CoefficientFunction((
            exp(-(1/(delta*delta))*((x-x0)*(x-x0)) ),
            -2*(1/(delta*delta))* (x-x0) * exp(-(1/(delta*delta))*((x-x0)*(x-x0)) ),
            0,
            0,
            0
            ))
        return sol

def standingwave(D, wavespeed):
    t = CoordCF(D)
    k = 1
    sq = sqrt(D)
    if(D==1):
        sol = CoefficientFunction((
            cos(math.pi*x)*sin(math.pi*t*wavespeed*sq)/(sq*math.pi),
            -sin(math.pi*x)*sin(math.pi*t*wavespeed*sq)/sq,
            cos(math.pi*x)*cos(math.pi*t*wavespeed*sq)*wavespeed
            ))
        return sol
    elif(D==2):
        sol = CoefficientFunction((
            cos(math.pi*x)*cos(math.pi*y)*sin(math.pi*t*wavespeed*sq)/(sq*math.pi),
            -sin(math.pi*x)*cos(math.pi*y)*sin(math.pi*t*wavespeed*sq)/sq,
            -cos(math.pi*x)*sin(math.pi*y)*sin(math.pi*t*wavespeed*sq)/sq,
            cos(math.pi*x)*cos(math.pi*y)*cos(math.pi*t*wavespeed*sq)*wavespeed
            ))
        return sol
    elif(D==3):
        sol = CoefficientFunction((
            cos(math.pi*x)*cos(math.pi*y)*cos(math.pi*z)*sin(math.pi*t*wavespeed*sq)/(sq*math.pi),
            -sin(math.pi*x)*cos(math.pi*y)*cos(math.pi*z)*sin(math.pi*t*wavespeed*sq)/sq,
            -cos(math.pi*x)*sin(math.pi*y)*cos(math.pi*z)*sin(math.pi*t*wavespeed*sq)/sq,
            -cos(math.pi*x)*cos(math.pi*y)*sin(math.pi*z)*sin(math.pi*t*wavespeed*sq)/sq,
            cos(math.pi*x)*cos(math.pi*y)*cos(math.pi*z)*cos(math.pi*t*wavespeed*sq)*wavespeed
            ))
        return sol

import ngsolve.special_functions


def singular(D, wavespeed):
    t = CoordCF(D)
    r2 = ((x)*(x)+(y)*(y))
    r = sqrt((x)*(x)+(y)*(y))
    at2 = atan2(-y,-x)+math.pi #2*atan((y)/(r+(x)))
    at2x=-(-y/r2)
    at2y=-(x/r2)
    alp = 10 #2.4048
    nrb = 2.0/3.0

    # sol = CoefficientFunction((
        # cos(alp*t)*bessel(alp*r),
        # cos(alp*t)*besselp(alp*r)*alp*x/r,
        # cos(alp*t)*besselp(alp*r)*alp*y/r,
        # -alp*sin(alp*t)*bessel(alp*r)
        # ))
    sol = CoefficientFunction((
        cos(alp*t)*sin(nrb*at2)*bessel(alp*r),
        cos(alp*t)*(cos(nrb*at2)*nrb*at2x*bessel(alp*r) + sin(nrb*at2)*besselp(alp*r)*alp*x/r),
        cos(alp*t)*(cos(nrb*at2)*nrb*at2y*bessel(alp*r) + sin(nrb*at2)*besselp(alp*r)*alp*y/r),
        -alp*sin(alp*t)*sin(nrb*at2)*bessel(alp*r)
        ))
    return sol
# cos(2*math.pi*t) * sin((2./3)*theta) * r2**(1./3),
    # cos(2*math.pi*t) * cos((2./3)*theta)*(2./3)*(1./(1+((y)/(x))**2))*((x)**(-1)) * (2./3)*r2**(-2./3)*(x),
    # cos(2*math.pi*t) * cos((2./3)*theta)*(2./3)*(1./(1+((y)/(x))**2))*(y)*(-(x)**(-2)) * (2./3)*r2**(-2./3)*(y),
    # -2*math.pi*sin(2*math.pi*t) * sin((2/3)*theta) * r2**(1./3)
