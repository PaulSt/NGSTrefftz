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
#########################################################################################################################################
fes = H1(mesh, order=4)

basemesh = Mesh(ngmeshbase)
n = specialcf.normal(3)
print(mesh.GetBoundaries())
for i in range(3):
    gfu = GridFunction(fes)
    gfu.Set(n[i],BND,definedon=mesh.Boundaries("default"))
    Draw(gfu,mesh,"gfu")
    input()

Draw(n,mesh,"n")

# fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True)
# D=2
# U = fes.TrialFunction()
# V = fes.TestFunction()
# gU = grad(U)
# gV = grad(V)
# v = gU[D]
# sig = CoefficientFunction(tuple([-gU[i] for i in  range(D)]))
# w = gV[D]
# tau = CoefficientFunction(tuple([-gV[i] for i in  range(D)]))

# n = specialcf.normal(D+1)
# n_t = n[D]
# n_x = CoefficientFunction( tuple([n[i] for i in  range(D)]) )

# a = LinearForm(fes)
# a += SymbolicLFI(tau*tau-w*w)
# a += SymbolicLFI(V*(-tau*n[0]+w*n[1]), element_boundary=True)
# a.Assemble()
# print(a.vec)
# import scipy as sp
# import scipy.sparse
# rows,cols,vals = a.mat.COO()
# A = sp.sparse.csr_matrix((vals,(rows,cols)))
# import numpy as np
# import numpy.linalg
# cond = np.linalg.cond(A.todense())

# print( Integrate(tau*tau-w*w + CoefficientFunction(V)*(-tau*n[0]+w*n[1]), basemesh) )
