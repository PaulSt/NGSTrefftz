#########################################################################################################################################
N = 3
c=2
t_steps = c*N
order = 4
k = 1
#########################################################################################################################################
import sys
sys.path.append("..")
from trefftzngs import *
import numpy as np
from prodmesh import CartSquare
from ngsolve import *
#import netgen.gui
from DGeq import *

mesh = CartSquare(N,t_steps)
#########################################################################################################################################

# fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True, basistype=1)
X = L2(mesh, order=order)
fes = FESpace([X,X], flags = { "dgjumps" : True })

truesol =  sin( k*(c*y + x) )#exp(-pow(c*x+y,2)))#
v0 = c*k*cos(k*(c*y+x))#grad(U0)[0]
sig0 = -k*cos(k*(c*y+x))#-grad(U0)[1]

U0 = GridFunction(fes)
U0.components[1].Set(v0)
U0.components[0].Set(sig0)
# Draw(U0.components[0],mesh,'U0')
# input()

sig,v = fes.TrialFunction()
tau,w = fes.TestFunction()

vo = v.Other()
sigo = sig.Other()
wo = w.Other()
tauo = tau.Other()

h = specialcf.mesh_size
n = specialcf.normal(2)
n_t = n[1]/Norm(n)
n_x = n[0]/Norm(n)

mean_v = 0.5*(v+vo)
mean_w = 0.5*(w+wo)
mean_sig = 0.5*(sig+sigo)
mean_tau = 0.5*(tau+tauo)

jump_vx = ( v - vo ) * n_x
jump_wx = ( w - wo ) * n_x
jump_sigx = ( sig - sigo ) * n_x
jump_taux = ( tau - tauo ) * n_x

jump_vt = ( v - vo ) * n_t
jump_wt = ( w - wo ) * n_t
jump_sigt = ( sig - sigo ) * n_t
jump_taut = ( tau - tauo ) * n_t

#jump_Ut = (U - U.Other()) * n_t

timelike = n_x**2 #IfPos(n_t,0,IfPos(-n_t,0,1)) # n_t=0
spacelike = n_t**2 #IfPos(n_x,0,IfPos(-n_x,0,1)) # n_x=0

alpha = 0.5 #pow(10,5)
beta = 0.5 #pow(10,5)
gamma = 1

a = BilinearForm(fes)
a += SymbolicBFI(  -v*(grad(tau)[0]+pow(c,-2)*grad(w)[1]) - sig*(grad(tau)[1]+grad(w)[0])  )
a += SymbolicBFI( spacelike * ( IfPos(n_t,v,vo)*(pow(c,-2)*jump_wt) + IfPos(n_t,sig,sigo)*(jump_taut) ) ,VOL,  skeleton=True ) #space like faces, no jump in x since horizontal
a += SymbolicBFI( timelike 	* ( mean_v*jump_taux + mean_sig*jump_wx + alpha*jump_vx*jump_wx + beta*jump_sigx*jump_taux ) ,VOL, skeleton=True ) #time like faces
a += SymbolicBFI( spacelike * IfPos(n_t,1,0) * ( pow(c,-2)*v*w + sig*tau ), BND, skeleton=True) #t=T (or *x)
a += SymbolicBFI( timelike 	* ( sig*n_x*w + alpha*v*w ), BND, skeleton=True) #dirichlet boundary 'timelike'
a.Assemble()

f = LinearForm(fes)
f += SymbolicLFI( spacelike * IfPos(-n_t,1,0) *  ( pow(c,-2)*v0*w + sig0*tau ), BND, skeleton=True) #t=0 (or *(1-x))
f += SymbolicLFI( timelike 	* ( v0 * (alpha*w - tau*n_x) ), BND, skeleton=True) #dirichlet boundary 'timelike'
f.Assemble()


[gfu,cond] = DGsolve(fes,a,f)

err0=U0.components[0] - gfu.components[0]
err1=U0.components[1] - gfu.components[1]
sH1error = sqrt(Integrate( err0*err0 + err1*err1, mesh))
print("grad-error=", sH1error)

Draw( gfu.components[0], mesh, "gradu" )
