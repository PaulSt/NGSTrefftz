from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.TensorProductTools import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from trefftzngs import *
# import netgen.gui
import scipy as sp
import scipy.sparse.linalg
import scipy.linalg
import time
from scipy.io import savemat
from scipy.io import loadmat
from ngsolve.bla import VectorD
from ngsolve import *

def GetFESTrefftz(mesh,c=1):
    return FESpace("trefftzfespace", mesh, order = 4, wavespeed = c, dgjumps=True, basistype=0)

order = 4
c = 2
t_start = 0
t_step = 0.1

# ngmesh = SegMesh(4,0,1)
# ngmesh = QadSegMesh(4,0,1)
# initmesh = Mesh(ngmesh)
initmesh = Mesh(unit_square.GenerateMesh(maxh=0.4))
# initmesh = Mesh(unit_cube.GenerateMesh(maxh = 0.5))
D = initmesh.dim
if D==3: eltyp = ET.TET
elif D==2: eltyp = ET.TRIG
elif D==1: eltyp = ET.SEGM
intrule = IntegrationRule(eltyp,2*order)


fes = H1(initmesh, order=order)
u = fes.TrialFunction()  # symbolic object
v = fes.TestFunction()   # symbolic object
gfu = GridFunction(fes)  # solution
a = BilinearForm(fes)
a += SymbolicBFI(u*v)
a.Assemble()
Draw(gfu,initmesh,'sol',autoscale=False,min=-1,max=1)
# VTKOutput object
vtk = VTKOutput(ma=initmesh, coefs=[gfu], names = ["sol"], filename="result", subdivision=3)

for t in range(0,50):
    wavefront = EvolveTentsMakeWavefront(order,initmesh,c,t_start)
    wavefront = EvolveTents(order,initmesh,c,t_step,wavefront,t_start)
    print(EvolveTentsPostProcess(order,initmesh,wavefront,EvolveTentsMakeWavefront(order,initmesh,c,t_start + t_step)))

    wavefront_nograd = VectorD(int(wavefront.NumPy().size/(D+2)))
    irsize = len(intrule.points)
    for n in range(0,initmesh.ne):
        for i in range(0,irsize):
            wavefront_nograd[n*irsize + i] = wavefront[n*irsize*(D+2) + i*(D+2)]


    test=IntegrationPointFunction(initmesh,intrule,wavefront_nograd)


    f = LinearForm(fes)
    f += SymbolicLFI(test*v, intrule=intrule)
    f.Assemble()
    gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs()) * f.vec
    Redraw()
    # Exporting the results:
    # vtk.Do()

    t_start += t_step

# A = mat.NumPy()[:,1:12]
# scipy.io.savemat('arrdata.mat', mdict={'arr': arr})

# fes = FESpace("trefftzfespace", initmesh,order=order,wavespeed=c,dgjumps=True,useshift=0)
# gfu = GridFunction(fes)
# gfu.vec.FV()[:] = vec
# Draw(gfu,initmesh,'bla')

#  mesh = NgsTPmesh(initmesh,c,1)
#  Draw(mesh)
