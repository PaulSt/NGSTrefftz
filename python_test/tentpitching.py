from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from prodmesh import *
from ngsolve import *
import netgen.gui
from trefftzngs import *
from DGeq import *
from SolveTentSlab import *
import time
from ngsolve.TensorProductTools import *
from ngsolve.comp import *

c = 3
order =12
k = 15

initmesh = Mesh(SegMesh(9,0,1))
mesh = NgsTPmesh(initmesh,c,0.4)
Draw(mesh)


fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True, basistype=0)

truesol = exp(-1.3*k*((x-0.15)-c*y)*((x-0.15)-c*y)) + exp(-k*((x-0.9)+c*y)*((x-0.9)+c*y)) #sin( k*(c*y + x) )#

U0 = GridFunction(fes)
U0.Set(truesol)
Draw(U0,mesh,'U0')
# input()

v0 = grad(U0)[1]#c*k*cos(k*(c*y+x))#
sig0 = -grad(U0)[0] #-k*cos(k*(c*y+x))#-grad(U0)[1]
gD = v0


[a,f] = SolveTentSlab(fes,U0,gD,c)


gfu = GridFunction(fes, name="uDG")
rows,cols,vals = a.mat.COO()
A = sp.sparse.csr_matrix((vals,(rows,cols)))
gfu.vec.FV().NumPy()[:] = sp.sparse.linalg.spsolve(A,f.vec.FV())


L2error = Integrate((truesol - gfu)*(truesol - gfu), mesh)
sH1error = Integrate((grad(U0) - grad(gfu))*(grad(U0) - grad(gfu)), mesh)
print("error=", sqrt(L2error))
print("grad-error=", sqrt(sH1error))
Draw(gfu,mesh,'sol')
Draw(grad(gfu),mesh,'gradsol')
