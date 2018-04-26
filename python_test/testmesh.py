#########################################################################################################################################
c = 1
t_stepsize = 1/5
nt_steps = 3
order = 4
k = 1
#########################################################################################################################################
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from prodmesh import PeriodicProdMesh
from ngsolve import *
# import netgen.gui
from trefftzngs import *
from DGeq import *
import time
# mesh = Mesh(unit_cube.GenerateMesh(maxh = 0.2,quad_dominated=True))
ngmeshbase = unit_square.GenerateMesh(maxh = t_stepsize)
mesh = PeriodicProdMesh(ngmeshbase,t_stepsize)
#########################################################################################################################################
# fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True, basistype=0)
fes = H1(mesh, order=order)

n = specialcf.normal(3)
basemesh = Mesh(ngmeshbase)
print(mesh.GetBoundaries())
for i in range(3):
    gfu = GridFunction(fes)
    gfu.Set(n[i],BND,definedon=mesh.Boundaries("dirichlet"))
    Draw(gfu,mesh,"gfu")
    input()

D=2
u = fes.TrialFunction()
v = fes.TestFunction()
v = gU[D]
sig = CoefficientFunction(tuple([-gU[i] for i in  range(D)]))
w = gV[D]
tau = CoefficientFunction(tuple([-gV[i] for i in  range(D)]))

n = specialcf.normal(D+1)
n_t = n[D]/Norm(n)
n_x = CoefficientFunction( tuple([n[i]/Norm(n) for i in  range(D)]) )


a = LinearForm(fes)
a += SymbolicLFI(tau*tau-w*w)
a += SymbolicLFI(-v*(tau*n_x-w*n_t), BND)

Draw(n,mesh,"n")
