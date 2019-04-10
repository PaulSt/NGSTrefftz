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
with TaskManager():
# with TaskManager(pajetrace=100*1000*1000):
    wavefront = EvolveTents(order,initmesh,c,t_step,wavefront,t_start,bdd)
print("time ",time.time()-start)
error = EvolveTentsError(order,initmesh,wavefront,EvolveTentsMakeWavefront(order,initmesh,t_start + t_step,bdd))
print("error ", error)

ipfct=IntegrationPointFunction(initmesh,intrule,wavefront)
f = LinearForm(fes)
f += SymbolicLFI(ipfct*v, intrule=intrule)
f.Assemble()
gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs()) * f.vec

print(sqrt(Integrate(( sin(math.pi*x)*sin(math.pi*y)*sin(math.pi*z)*sin(math.pi*t_step*c*sq)/(sq*math.pi)-gfu)*( sin(math.pi*x)*sin(math.pi*y)*sin(math.pi*z)*sin(math.pi*t_step*c*sq)/(sq*math.pi)-gfu),initmesh)))
