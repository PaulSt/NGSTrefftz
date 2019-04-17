from ngsolve import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
import time
import math
from trefftzngs import *
# SetHeapSize(1000*1000*1000)

order = 3
SetNumThreads(12)

c = 1
t_start = 0
t_step = 2/sqrt(3)

h1error = []
adiam = []
timer = []
# ms = [0.4,0.3,0.2,0.1]
ms = [0,1,2,3,4]

initmesh = Mesh(unit_cube.GenerateMesh(maxh = 0.5))
for maxh in ms:
    print("RUN: ", maxh)
    # initmesh = Mesh(unit_cube.GenerateMesh(maxh = maxh))
    # initmesh = Mesh(unit_square.GenerateMesh(maxh = maxh))

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
    wavefront = EvolveTentsMakeWavefront(order,initmesh,t_start,bdd)

    fes = H1(initmesh, order=order)
    u,v = fes.TnT()
    gfu = GridFunction(fes)
    a = BilinearForm(fes)
    a += SymbolicBFI(u*v)
    a.Assemble()
    eltyp = ET.TET
    intrule = IntegrationRule(eltyp,2*order)
    irsize = len(intrule.points)

    start = time.time()
    with TaskManager():
        wavefront = EvolveTents(order,initmesh,c,t_step,wavefront,t_start,bdd)
    timer.append(time.time()-start)
    print("time ",time.time()-start)
    error = EvolveTentsError(order,initmesh,wavefront,EvolveTentsMakeWavefront(order,initmesh,t_step,bdd))
    h1error.append(error)
    print("error ", error)
    adiam.append(EvolveTentsAdiam(initmesh,1,t_step))

    if maxh != ms[-1]:
        initmesh.Refine()

print("rate")
from math import *
for i in range(len(ms)-1):
    print(log(h1error[i]/h1error[i+1])/log(2))
    # print(log(h1error[i]/h1error[i+1])/log(ms[i]/ms[i+1]))
for i in range(len(ms)):
    print("adiam ",adiam[i]," maxh ", ms[i]," error ", h1error[i], " time ", timer[i])
for i in range(len(ms)-1):
    print(log(h1error[i]/h1error[i+1])/log(adiam[i]/adiam[i+1]))
