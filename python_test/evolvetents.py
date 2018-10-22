from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from trefftzngs import *
import netgen.gui
import scipy as sp
import scipy.sparse.linalg
import scipy.linalg
import time
from scipy.io import savemat
from scipy.io import loadmat
from ngsolve import *

def GetFESTrefftz(mesh,c=1):
    return FESpace("trefftzfespace", mesh, order = 4, wavespeed = c, dgjumps=True, basistype=0)

order = 4
c = 1
t_start = 0
t_step = 0.1

# ngmesh = SegMesh(4,0,1)
# ngmesh = QadSegMesh(4,0,1)
# initmesh = Mesh(ngmesh)
initmesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
# initmesh = Mesh(unit_cube.GenerateMesh(maxh = 0.5))
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
Draw(gfu,initmesh,'sol',autoscale=True,min=-1,max=1)
wavefront = EvolveTentsMakeWavefront(order,initmesh,c,t_start)

for t in range(0,500):
    wavefront = EvolveTents(order,initmesh,c,t_step,wavefront,t_start)
    print(EvolveTentsPostProcess(order,initmesh,wavefront,EvolveTentsMakeWavefront(order,initmesh,c,t_start + t_step)))

    ipfct=IntegrationPointFunction(initmesh,intrule,wavefront[0:initmesh.ne*irsize])
    f = LinearForm(fes)
    f += SymbolicLFI(ipfct*v, intrule=intrule)
    f.Assemble()
    gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs()) * f.vec
    Redraw()

    t_start += t_step
    # filename = "results/mov/sol"+str(t).zfill(3) +".jpg"
    # Tcl_Eval("Ng_SnapShot .ndraw {};\n".format(filename))

# A = mat.NumPy()[:,1:12]
# scipy.io.savemat('arrdata.mat', mdict={'arr': arr})

# fes = FESpace("trefftzfespace", initmesh,order=order,wavespeed=c,dgjumps=True,useshift=0)
# gfu = GridFunction(fes)
# gfu.vec.FV()[:] = vec
# Draw(gfu,initmesh,'bla')

#  mesh = NgsTPmesh(initmesh,c,1)
#  Draw(mesh)
