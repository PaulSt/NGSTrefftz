#########################################################################################################################################
N = 4
c=2
t_steps = c*N
order = 5
k = 5
#########################################################################################################################################
from trefftzngs import *
import numpy as np
from prodmesh import CartSquare
from ngsolve import *
import netgen.gui
from DGeq import *
import time

mesh = CartSquare(N,t_steps)

#########################################################################################################################################

truesol =  sin( k*(c*y + x) )#exp(-pow(c*x+y,2)))#
v0 = c*k*cos(k*(c*y+x))#grad(U0)[0]
sig0 = -k*cos(k*(c*y+x))#-grad(U0)[1]

# for order in range(3,order):
fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True, basistype=0)
# fes = L2(mesh, order=order, flags = { "dgjumps" : True })

U0 = GridFunction(fes)
U0.Set(truesol)
Draw(U0,mesh,'U0')
# input()

start = time.clock()
[a,f] = DGeqsys(fes,truesol,v0,sig0,c,v0)
print("DGsys: ", str(time.clock()-start))
# gfu2= GridFunction(fes, name="uDG")
# gfu2.vec.data = a.mat.Inverse() * f.vec

[gfu,cond] = DGsolve(fes,a,f)

L2error = Integrate((truesol - gfu)*(truesol - gfu), mesh)
sH1error = Integrate((grad(U0) - grad(gfu))*(grad(U0) - grad(gfu)), mesh)
print("error=", L2error)
print("grad-error=", sH1error)
Draw(gfu,mesh,'sol')
Draw(grad(gfu),mesh,'gradsol')









# nmatclean = nmat[nmat.any(axis=0),:]
# nmatclean = nmatclean[:,nmat.any(axis=1)]
# nvecclean = f.vec.FV().NumPy()[nmat.any(axis=1)] #nvec
# solclean = np.linalg.solve(nmat,nvec)
# sol = np.zeros(a.mat.height)
# sol[nmat.any(axis=1)] = solclean
# print(nmat.any(axis=1))
