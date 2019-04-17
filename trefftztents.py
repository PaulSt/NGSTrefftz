from ngsolve import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
import time
import math
from trefftzngs import *
SetHeapSize(1000*1000*1000)

order = 4

c = 1
t_start = 0
t_step = 2/sqrt(3)

h1error = []
ms = [0.25,0.2,0.15,0.13]

for maxh in ms:
    initmesh = Mesh(unit_cube.GenerateMesh(maxh = maxh))
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
    SetNumThreads(8)
# with TaskManager(pajetrace=100*1000*1000):
    with TaskManager():
        wavefront = EvolveTents(order,initmesh,c,t_step,wavefront,t_start,bdd)
    print("time ",time.time()-start)
    error = EvolveTentsError(order,initmesh,wavefront,EvolveTentsMakeWavefront(order,initmesh,t_step,bdd))
    h1error.append(error)
    print("error ", error)

    ipfct=IntegrationPointFunction(initmesh,intrule,wavefront)
    f = LinearForm(fes)
    f += SymbolicLFI(ipfct*v, intrule=intrule)
    f.Assemble()
    gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs()) * f.vec

    # print(sqrt(Integrate(( sin(math.pi*x)*sin(math.pi*y)*sin(math.pi*z)*sin(math.pi*t_step*c*sq)/(sq*math.pi)-gfu)*( sin(math.pi*x)*sin(math.pi*y)*sin(math.pi*z)*sin(math.pi*t_step*c*sq)/(sq*math.pi)-gfu),initmesh)))

    # corsol= CoefficientFunction((cos(math.pi*x)*sin(math.pi*y)*sin(math.pi*z)*sin(math.pi*t_step*c*sq)/sq,
        # sin(math.pi*x)*cos(math.pi*y)*sin(math.pi*z)*sin(math.pi*t_step*c*sq)/sq,
        # sin(math.pi*x)*sin(math.pi*y)*cos(math.pi*z)*sin(math.pi*t_step*c*sq)/sq))
    # print(sqrt(Integrate(( corsol-grad(gfu))*( corsol-grad(gfu)),initmesh)))
    adiam.append(EvolveTentsAdiam(initmesh,1,t_step))

    initmesh.Refine()

print("rate")
from math import *
for i in range(len(ms)-1):
    print(log(h1error[i]/h1error[i+1])/log(1/2))
    # print(log(h1error[i]/h1error[i+1])/log(ms[i]/ms[i+1]))
for i in range(len(ms)):
    print("adiam ",adiam[i]," maxh ", ms[i]," error ", h1error[i])
for i in range(len(ms)-1):
    print(log(h1error[i]/h1error[i+1])/log(adiam[i]/adiam[i+1]))
