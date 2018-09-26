from ngsolve import *
from ngsolve.comp import *
from ngsolve.TensorProductTools import *
import netgen.meshing as ngm
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from trefftzngs import *
from DGeq import *
# import netgen.gui
import scipy as sp
import scipy.sparse.linalg
import scipy.linalg
import time
from scipy.io import savemat
from scipy.io import loadmat
from netgen.csg import unit_cube
from prodmesh import *



def GetFESTrefftz(mesh,c=1):
    return FESpace("trefftzfespace", mesh, order = 4, wavespeed = c, dgjumps=True, basistype=0)

order = 9
c = 3
t_end = 1
# ngmesh = SegMesh(4,0,1)
ngmesh = QadSegMesh(4,0,1)
initmesh = Mesh(ngmesh)
# initmesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
# initmesh = Mesh(unit_cube.GenerateMesh(maxh = 0.6))
wavefront = EvolveTentsMakeWavefront(order,initmesh,c,0)
mat = EvolveTents(order,initmesh,c,t_end,wavefront,0)
print(EvolveTentsPostProcess(order,initmesh,mat,EvolveTentsMakeWavefront(order,initmesh,c,t_end)))
# A = mat.NumPy()[:,1:12]
# scipy.io.savemat('arrdata.mat', mdict={'arr': arr})

# fes = FESpace("trefftzfespace", initmesh,order=order,wavespeed=c,dgjumps=True,useshift=0)
# gfu = GridFunction(fes)
# gfu.vec.FV()[:] = vec
# Draw(gfu,initmesh,'bla')

#  mesh = NgsTPmesh(initmesh,c,1)
#  Draw(mesh)
