from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from trefftzngs import *
import netgen.gui
from ngsolve.TensorProductTools import *
from ngsolve import *
from prodmesh import *
from ngsolve.solve import Tcl_Eval # for snapshots
from testcases import *

order = 4
c = 1
t_start = 0
t_step = 0.1

ngmesh = SegMesh(4,0,1)
# ngmesh = QadSegMesh(9,0,1)
initmesh = Mesh(ngmesh)
initmesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
# initmesh = Mesh(unit_cube.GenerateMesh(maxh = 0.4))
# initmesh = Mesh( CircleMesh(0.2) )
# for i in range(0,len(initmesh.GetBoundaries())):
   # initmesh.ngmesh.SetBCName(i,"neumann")

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
# Draw(gfu,initmesh,'sol',autoscale=False,min=-0.1,max=0.1)
# Draw(gfu,initmesh,'sol',autoscale=False,min=-0.01,max=0.01)
bdd = standingwave(D,c)
wavefront = EvolveTentsMakeWavefront(order,initmesh,c,t_start,bdd)

for t in range(0,200):
# with TaskManager():
    wavefront = EvolveTents(order,initmesh,c,t_step,wavefront,t_start,bdd)
    print("Error: " + str(EvolveTentsError(order,initmesh,wavefront,EvolveTentsMakeWavefront(order,initmesh,c,t_start + t_step,bdd))))

    ipfct=IntegrationPointFunction(initmesh,intrule,wavefront)
    f = LinearForm(fes)
    f += SymbolicLFI(ipfct*v, intrule=intrule)
    f.Assemble()
    gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs()) * f.vec
    Redraw()

    t_start += t_step
    print("time: " + str(t_start))
# filename = "results/mov/sol"+str(t).zfill(3) +".jpg"
# Tcl_Eval("Ng_SnapShot .ndraw {};\n".format(filename))
