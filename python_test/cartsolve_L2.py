#########################################################################################################################################
N = 3
c=2
t_steps = c*N
order = 4
k = 1
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

truesol =  sin( k*(c*y + x) )#exp(-pow(c*x+y,2)))#
truesol = exp(-3*k*((x-0.5)-c*y)*((x-0.5)-c*y)) 
U0 = GridFunction(fes)
U0.Set(truesol)
Draw(U0,mesh,'U0')
v0 = grad(U0)[1]# c*k*cos(k*(c*y+x))#
sig0 = -grad(U0)[0]#-k*cos(k*(c*y+x))#

U0 = GridFunction(fes)
U0.Set(truesol)
Draw(U0,mesh,'U0')
# testhes = U0.Operator("hesse");
# Draw(testhes,mesh,'h')
#input()

[a,f] = DGeqsys(fes,truesol,v0,sig0,c,v0,True)

[gfu,cond] = DGsolve(fes,a,f)

L2error = Integrate((truesol - gfu)*(truesol - gfu), mesh)
sH1error = Integrate((grad(U0) - grad(gfu))*(grad(U0) - grad(gfu)), mesh)
print("error=", sqrt(L2error))
print("grad-error=", sqrt(sH1error))
Draw(gfu,mesh,'sol')
Draw(grad(gfu),mesh,'gradsol')
