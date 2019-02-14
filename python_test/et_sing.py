
import netgen.gui
from ngsolve import *
from ngsolve.solve import Tcl_Eval # for snapshots
from trefftzngs import *
from prodmesh import *
from testcases import *
import time

from ngsolve import *
from netgen.meshing import MeshingParameters

maxh=0.01
minh=0.005
mp = MeshingParameters (maxh = maxh)
refpoints = 500
for i in range(0, refpoints+1):
    for j in range(0, refpoints+1):
        xk = i/refpoints-0.5
        yk = j/refpoints-0.5
        mp.RestrictH (x=xk, y=yk, z=0, h=max(minh,sqrt(maxh*((xk)*(xk)+(yk)*(yk)))))

initmesh = Mesh( oLshapeMesh(maxh,mp) )
# RefineAround([0.5,0.5,0],0.1,initmesh)
# RefineAround([0.5,0.5,0],0.02,initmesh)
Draw(initmesh)

for i in range(0,len(initmesh.GetBoundaries())):
   initmesh.ngmesh.SetBCName(i,"neumann")

order = 3
c = 1
t_start = 0
t_step = 0.3

D = 2
eltyp = ET.TRIG
intrule = IntegrationRule(eltyp,2*order)
irsize = len(intrule.points)

bdd = singular(D,c)
Draw(bdd,initmesh,'u')
fes = H1(initmesh, order=order)
u,v = fes.TnT()
gfu = GridFunction(fes)
a = BilinearForm(fes)
a += SymbolicBFI(u*v)
a.Assemble()
wavefront = EvolveTentsMakeWavefront(order,initmesh,t_start,bdd )
ipfct=IntegrationPointFunction(initmesh,intrule,wavefront)
f = LinearForm(fes)
f += SymbolicLFI(ipfct*v, intrule=intrule)
f.Assemble()
gfu.vec.data = a.mat.Inverse() * f.vec
Draw(gfu,initmesh,'sol')
input()
# Draw(gfu,initmesh,'sol',autoscale=False,min=-0.1,max=0.35)

runtime = time.time()
with TaskManager():
    for t in range(0,1):
        wavefront = EvolveTents(order,initmesh,c,t_step,wavefront,t_start, bdd )

        ipfct=IntegrationPointFunction(initmesh,intrule,wavefront)
        f = LinearForm(fes)
        f += SymbolicLFI(ipfct*v, intrule=intrule)
        f.Assemble()
        gfu.vec.data = a.mat.Inverse() * f.vec
        Redraw(blocking=True)

        t_start += t_step
        print("time: " + str(t_start))
        # filename = "results/mov/sol"+str(t).zfill(3)+".jpg"
        # Tcl_Eval("Ng_SnapShot .ndraw {};\n".format(filename))
print("\n",time.time()-runtime)
