from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from prodmesh import *
from ngsolve import *
import netgen.gui
from trefftzngs import *
from DGeq import *
import time
from ngsolve.TensorProductTools import *
from ngsolve.comp import *
#  basemeshsize = 0.3
#  ngmeshbase = unit_square.GenerateMesh(maxh = basemeshsize)
#  mesh = PeriodicProdMesh(ngmeshbase,t_stepsize)

mesh = Mesh(SegMesh(4,0,0.5))
tpmesh = ngs_tpmesh(mesh,1)
Draw(tpmesh)
# mesh2 = Mesh(unit_square.GenerateMesh(maxh=0.4))
# mesh = Mesh(MakeTensorProductMesh(mesh2,mesh1))
