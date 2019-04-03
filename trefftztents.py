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

D = initmesh.dim

c = 1
sq = sqrt(3.0);
t = CoordCF(D)
bdd = CoefficientFunction((
    sin(math.pi*x)*sin(math.pi*y)*sin(math.pi*z)*sin(math.pi*t*c*sq)/(sq*math.pi),
    cos(math.pi*x)*sin(math.pi*y)*sin(math.pi*z)*sin(math.pi*t*c*sq)/sq,
    sin(math.pi*x)*cos(math.pi*y)*sin(math.pi*z)*sin(math.pi*t*c*sq)/sq,
    sin(math.pi*x)*sin(math.pi*y)*cos(math.pi*z)*sin(math.pi*t*c*sq)/sq,
    sin(math.pi*x)*sin(math.pi*y)*sin(math.pi*z)*cos(math.pi*t*c*sq)*c
    ))
wavefront = EvolveTentsMakeWavefront(order,initmesh,t_start,bdd)

start = time.time()
with TaskManager():
# SetNumThreads(1)
# with TaskManager(pajetrace=100*1000*1000):
    wavefront = EvolveTents(order,initmesh,c,t_step,wavefront,t_start,bdd)
print("time ",time.time()-start)
error = EvolveTentsError(order,initmesh,wavefront,EvolveTentsMakeWavefront(order,initmesh,t_start + t_step,bdd))
print("error ", error)
