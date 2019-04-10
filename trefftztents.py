from ngsolve import *
from netgen.csg import unit_cube
import time
import math
from trefftzngs import *

order = 4

c = 1
t_start = 0
t_step = 2/sqrt(3)

initmesh = Mesh(unit_cube.GenerateMesh(maxh = 0.2))
# for i in range(0,len(initmesh.GetBoundaries())):
   # initmesh.ngmesh.SetBCName(i,"neumann")

D = initmesh.dim
t = CoordCF(D)

sq = sqrt(3.0);
bdd = CoefficientFunction((
    sin(math.pi*x)*sin(math.pi*y)*sin(math.pi*z)*sin(math.pi*t*c*sq)/(sq*math.pi),
    cos(math.pi*x)*sin(math.pi*y)*sin(math.pi*z)*sin(math.pi*t*c*sq)/sq,
    sin(math.pi*x)*cos(math.pi*y)*sin(math.pi*z)*sin(math.pi*t*c*sq)/sq,
    sin(math.pi*x)*sin(math.pi*y)*cos(math.pi*z)*sin(math.pi*t*c*sq)/sq,
    sin(math.pi*x)*sin(math.pi*y)*sin(math.pi*z)*cos(math.pi*t*c*sq)*c
    ))
# sq = sqrt(1.0/3.0);
# bdd = CoefficientFunction((
    # sin( c*t+sq*(x+y+z) ),
    # sq*cos(c*t+sq*(x+y+z)),
    # sq*cos(c*t+sq*(x+y+z)),
    # sq*cos(c*t+sq*(x+y+z)),
    # c*cos(c*t+sq*(x+y+z))
    # ))
# bdd = CoefficientFunction((
    # sin(sq*x)*sin(sq*y)*sin(sq*z)*sin(t*c),
    # sq*cos(sq*x)*sin(sq*y)*sin(sq*z)*sin(t*c),
    # sq*sin(sq*x)*cos(sq*y)*sin(sq*z)*sin(t*c),
    # sq*sin(sq*x)*sin(sq*y)*cos(sq*z)*sin(t*c),
    # sin(sq*x)*sin(sq*y)*sin(sq*z)*cos(t*c)*c
    # ))
# bdd = CoefficientFunction((
# cos(math.pi*x)*cos(math.pi*y)*cos(math.pi*z)*sin(math.pi*t*c*sq)/(sq*math.pi),
# -sin(math.pi*x)*cos(math.pi*y)*cos(math.pi*z)*sin(math.pi*t*c*sq)/sq,
# -cos(math.pi*x)*sin(math.pi*y)*cos(math.pi*z)*sin(math.pi*t*c*sq)/sq,
# -cos(math.pi*x)*cos(math.pi*y)*sin(math.pi*z)*sin(math.pi*t*c*sq)/sq,
# cos(math.pi*x)*cos(math.pi*y)*cos(math.pi*z)*cos(math.pi*t*c*sq)*c
# ))
wavefront = EvolveTentsMakeWavefront(order,initmesh,t_start,bdd)

start = time.time()
# SetNumThreads(1)
# with TaskManager(pajetrace=100*1000*1000):
with TaskManager():
    wavefront = EvolveTents(order,initmesh,c,t_step,wavefront,t_start,bdd)
print("time ",time.time()-start)
error = EvolveTentsError(order,initmesh,wavefront,EvolveTentsMakeWavefront(order,initmesh,t_start + t_step,bdd))
print("error ", error)
