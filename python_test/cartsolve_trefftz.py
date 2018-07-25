#########################################################################################################################################
N = 3
c=2
t_steps = c*N
order = 3
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


# for order in range(3,order):
fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True, basistype=0)
# fes = L2(mesh, order=order, flags = { "dgjumps" : True })

[truesol,U0,sig0,v0,gD] = TestSolution(fes,c)

start = time.clock()
[a,f] = DGeqsys(fes,truesol,v0,sig0,c,gD)
print("DGsys: ", str(time.clock()-start))

start = time.clock()
[gfu,cond] = DGsolve(fes,a,f)
print("DGsolve: ", str(time.clock()-start))

L2error = Integrate((truesol - gfu)*(truesol - gfu), mesh)
sH1error = Integrate((grad(U0) - grad(gfu))*(grad(U0) - grad(gfu)), mesh)
print("error=", sqrt(L2error))
print("grad-error=", sqrt(sH1error))
Draw(gfu,mesh,'sol')
Draw(grad(gfu),mesh,'gradsol')
