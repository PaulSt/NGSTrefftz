#########################################################################################################################################
c = 1
t_stepsize = 1/5
nt_steps = 3
order = 4
k = 1
#########################################################################################################################################
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from prodmesh import ProdMesh
from ngsolve import *
import netgen.gui
from trefftzngs import *
from DGeq import *
import time
# mesh = Mesh(unit_cube.GenerateMesh(maxh = 0.2,quad_dominated=True))
ngmeshbase = unit_square.GenerateMesh(maxh = t_stepsize)
mesh = ProdMesh(ngmeshbase,t_stepsize)
#########################################################################################################################################
# fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True, basistype=0)
fes = H1(mesh, order=order)

n = specialcf.normal(3)
basemesh = Mesh(ngmeshbase)
print(mesh.GetBoundaries())
gfu = GridFunction(fes)
gfu.Set(CoefficientFunction(n[1]),BND,definedon=mesh.Boundaries("dirichlet"))
Draw(gfu,mesh,"gfu")
Draw(n,mesh,"n")
