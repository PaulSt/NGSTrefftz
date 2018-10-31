#########################################################################################################################################
t_stepsize = 1/5
delta_x = 0.2
#########################################################################################################################################
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from prodmesh import *
from ngsolve import *
import netgen.gui
from trefftzngs import *
from DGeq import *
import time
# mesh = Mesh(unit_cube.GenerateMesh(maxh = 0.2,quad_dominated=True))
ngmeshbase = unit_square.GenerateMesh(maxh = delta_x)
mesh = PeriodicProdMesh(ngmeshbase,t_stepsize)

from ngsolve.TensorProductTools import *
from ngsolve.comp import *
mesh1 = Mesh(SegMesh(4,0,0.5))
mesh2 = Mesh(unit_square.GenerateMesh(maxh=0.4))
# mesh = Mesh(MakeTensorProductMesh(mesh2,mesh1))
basemesh = Mesh(SegMesh(4,0,1))
mesh = NgsTPmesh(basemesh,1,1)

mesh = Mesh( LshapeMesh(0.4) )

#########################################################################################################################################

Draw(mesh)
input()

fes = H1(mesh, order=4)
gfu = GridFunction(fes)
n = specialcf.normal(2)
print(mesh.GetBoundaries())
for i in range(2):
    gfu.Set(n[i],BND,definedon=mesh.Boundaries("inflow"))
    Draw(gfu,mesh,"gfu")
    input()
for i in range(2):
    gfu.Set(n[i],BND,definedon=mesh.Boundaries("outflow"))
    Draw(gfu,mesh,"gfu")
    input()
for i in range(2):
    gfu.Set(n[i],BND,definedon=mesh.Boundaries("dirichlet"))
    Draw(gfu,mesh,"gfu")
    input()
