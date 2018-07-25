#########################################################################################################################################
N = 3
c=2
t_steps = c*N
order = 4
#########################################################################################################################################
from trefftzngs import *
import numpy as np
from prodmesh import CartSquare
from ngsolve import *
import netgen.gui
from DGeq import *

mesh = CartSquare(N,t_steps)
#########################################################################################################################################

fes = L2(mesh, order=order, dgjumps=True)

[truesol,U0,sig0,v0,gD] = TestSolution(fes,c)

[a,f] = DGeqsys(fes,truesol,v0,sig0,c,gD,True)

[gfu,cond] = DGsolve(fes,a,f)

L2error = Integrate((truesol - gfu)*(truesol - gfu), mesh)
sH1error = Integrate((grad(U0) - grad(gfu))*(grad(U0) - grad(gfu)), mesh)
print("error=", sqrt(L2error))
print("grad-error=", sqrt(sH1error))
Draw(gfu,mesh,'sol')
Draw(grad(gfu),mesh,'gradsol')
