from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from trefftzngs import *
import netgen.gui
from ngsolve import *
from prodmesh import *
from ngsolve.solve import Tcl_Eval # for snapshots

from ngsolve import *
from netgen.geom2d import unit_square
from netgen.meshing import MeshingParameters

maxh=0.01
# mp = MeshingParameters (maxh = maxh)
# for i in range(0, 101):
    # for j in range(0, 101):
        # xk = i/100
        # yk = j/100
        # mp.RestrictH (x=xk, y=yk, z=0, h=sqrt((xk-0.5)*(xk-0.5)+(yk-0.5)*(yk-0.5)))

order = 3
c = 1
t_start = 0
t_step = 0.02
testcase = "vertgausspw"

initmesh = Mesh( LshapeMesh(maxh) )
RefineAround([0.5,0.5,0],0.1,initmesh)
# RefineAround([0.5,0.5,0],0.02,initmesh)
Draw(initmesh)
input()
for i in range(0,len(initmesh.GetBoundaries())):
   initmesh.ngmesh.SetBCName(i,"neumann")

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
# Draw(gfu,initmesh,'sol')
# Draw(gfu,initmesh,'sol',autoscale=False,min=-1,max=1)
Draw(gfu,initmesh,'sol',autoscale=False,min=-0.15,max=0.3)
wavefront = EvolveTentsMakeWavefront(order,initmesh,c,t_start,testcase )

for t in range(0,200):
    wavefront = EvolveTents(order,initmesh,c,t_step,wavefront,t_start )

    ipfct=IntegrationPointFunction(initmesh,intrule,wavefront)
    f = LinearForm(fes)
    f += SymbolicLFI(ipfct*v, intrule=intrule)
    f.Assemble()
    gfu.vec.data = a.mat.Inverse() * f.vec
    Redraw()

    t_start += t_step
    print("time: " + str(t_start))
    # filename = "results/mov/sol"+str(t).zfill(3) +".jpg"
    # Tcl_Eval("Ng_SnapShot .ndraw {};\n".format(filename))
