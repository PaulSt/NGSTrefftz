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


def gausspw(D, wavespeed):
    t = CoordCF(D)
    k = 1
    delta = 8000
    if(D==1):
        sol = CoefficientFunction((
            exp( -delta*((x-0.5)*(x-0.5)) ),
            -2*delta * (x-0.5) * exp( -delta*((x-0.5)*(x-0.5)) ),
            0
            ))
        return sol
    elif(D==2):
        sol = CoefficientFunction((
            exp(-delta*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)) ),
            -2*delta * (x-0.5) * exp(-delta*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)) ),
            -2*delta * (y-0.5) * exp(-delta*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)) ),
            0
            ))
        return sol
    elif(D==3):
        sol = CoefficientFunction((
            exp(-delta*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5)) ),
            -2*delta * (x-0.5) * exp(-delta*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5)) ),
            -2*delta * (y-0.5) * exp(-delta*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5)) ),
            -2*delta * (z-0.5) * exp(-delta*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5)) ),
            0
            ))
        return sol

def vertgausspw(D, wavespeed):
    t = CoordCF(D)
    k = 1
    delta = 80000
    if(D==1):
        sol = CoefficientFunction((
            exp( -delta*((x-0.25)*(x-0.25)) ),
            -2*delta * (x-0.25) * exp( -delta*((x-0.25)*(x-0.25)) ),
            0
            ))
        return sol
    elif(D==2):
        sol = CoefficientFunction((
            exp(-delta*((x-0.25)*(x-0.25)) ),
            -2*delta * (x-0.25) * exp(-delta*((x-0.25)*(x-0.25)) ),
            0,
            0
            ))
        return sol
    elif(D==3):
        sol = CoefficientFunction((
            exp(-delta*((x-0.25)*(x-0.25)) ),
            -2*delta * (x-0.25) * exp(-delta*((x-0.25)*(x-0.25)) ),
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
