from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from trefftzngs import *
import netgen.gui
from ngsolve import *
from prodmesh import *
from ngsolve.solve import Tcl_Eval # for snapshots
from testcases import *

from ngsolve import *
from netgen.geom2d import unit_square
from netgen.meshing import MeshingParameters

maxh=0.01
minh=0.0005
mp = MeshingParameters (maxh = maxh)
refpoints = 500
for i in range(0, refpoints+1):
    for j in range(0, refpoints+1):
        xk = i/refpoints
        yk = j/refpoints
        mp.RestrictH (x=xk, y=yk, z=0, h=max(minh,sqrt(0.005*((xk-0.5)*(xk-0.5)+(yk-0.5)*(yk-0.5)))))

initmesh = Mesh( LshapeMesh(maxh,mp) )
# RefineAround([0.5,0.5,0],0.1,initmesh)
# RefineAround([0.5,0.5,0],0.02,initmesh)
Draw(initmesh)
input()
for i in range(0,len(initmesh.GetBoundaries())):
   initmesh.ngmesh.SetBCName(i,"neumann")

order = 2
c = 1
t_start = 0
t_step = 0.05

D = initmesh.dim
if D==3: eltyp = ET.TET
elif D==2: eltyp = ET.TRIG
elif D==1: eltyp = ET.SEGM
intrule = IntegrationRule(eltyp,2*order)
irsize = len(intrule.points)

fes = H1(initmesh, order=order)
u,v = fes.TnT()
gfu = GridFunction(fes)
a = BilinearForm(fes)
a += SymbolicBFI(u*v)
a.Assemble()
Draw(gfu,initmesh,'sol')
# Draw(gfu,initmesh,'sol',autoscale=False,min=-1,max=1)
bdd = vertgausspw(D,c)
wavefront = EvolveTentsMakeWavefront(order,initmesh,t_start,bdd )
# Draw(bdd,initmesh,'sol',autoscale=False,min=-0.05,max=0.1)
# input()
# filename = "results/mov/sol"+str(0000)
# Tcl_Eval("Ng_SnapShot .ndraw {};\n".format(filename))
# input()
# Draw(gfu,initmesh,'sol',autoscale=False,min=-0.05,max=0.1)

for t in range(0,100):
    wavefront = EvolveTents(order,initmesh,c,t_step,wavefront,t_start, bdd )

    ipfct=IntegrationPointFunction(initmesh,intrule,wavefront)
    f = LinearForm(fes)
    f += SymbolicLFI(ipfct*v, intrule=intrule)
    f.Assemble()
    gfu.vec.data = a.mat.Inverse() * f.vec
    Redraw()

    t_start += t_step
    print("time: " + str(t_start))
    filename = "results/mov/sol"+str(t).zfill(3)
    Tcl_Eval("Ng_SnapShot .ndraw {};\n".format(filename))
