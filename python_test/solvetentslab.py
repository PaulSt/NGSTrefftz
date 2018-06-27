from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from prodmesh import *
from ngsolve import *
import netgen.gui
from trefftzngs import *
from DGeq import *
import time
from ngsolve.TensorProductTools import *
from ngsolve.comp import *


def SolveTentSlab(fes,U0,gD=0,c=1):
    v0 = grad(U0)[1]
    sig0 = -grad(U0)[0]

    D = fes.mesh.dim - 1
    U = fes.TrialFunction()
    V = fes.TestFunction()
    gU = grad(U)
    gV = grad(V)

    v = gU[D]
    sig = CoefficientFunction(tuple([-gU[i] for i in  range(D)]))
    w = gV[D]
    tau = CoefficientFunction(tuple([-gV[i] for i in  range(D)]))

    vo = gU.Other()[D]
    sigo = CoefficientFunction(tuple([-gU.Other()[i] for i in  range(D)]))
    wo = gV.Other()[D]
    tauo = CoefficientFunction(tuple([-gV.Other()[i] for i in  range(D)]))

    h = specialcf.mesh_size
    n = specialcf.normal(D+1)
    n_t = n[D]/Norm(n)
    n_x = CoefficientFunction( tuple([n[i]/Norm(n) for i in  range(D)]) )

    mean_v = 0.5*(v+vo)
    mean_w = 0.5*(w+wo)
    mean_sig = 0.5*(sig+sigo)
    mean_tau = 0.5*(tau+tauo)

    jump_vx = ( v - vo ) * n_x
    jump_wx = ( w - wo ) * n_x
    jump_sigx = (( sig - sigo ) * n_x)
    jump_taux = (( tau - tauo ) * n_x)

    jump_vt = ( v - vo ) * n_t
    jump_wt = ( w - wo ) * n_t
    jump_sigt = ( sig - sigo ) * n_t
    jump_taut = ( tau - tauo ) * n_t

    jump_Ut = (U - U.Other()) * n_t
    #gfu.vec.data = a.mat.Inverse() * f.vec

    timelike = n_x*n_x # n_t=0
    spacelike = n_t**2 # n_x=0

    alpha = 0.5 #pow(10,5)
    beta = 0.5 #pow(10,5)
    gamma = 1

    a = BilinearForm(fes)
    #  a += SymbolicBFI( ( IfPos(n_t,v,vo)*(pow(c,-2)*jump_wt) + IfPos(n_t,sig,sigo)*(jump_taut) ), VOL, skeleton=True ) #space like faces, w/o x jump
    a += SymbolicBFI( ( IfPos(n_t,v,vo)*(pow(c,-2)*jump_wt+jump_taux) + IfPos(n_t,sig,sigo)*(jump_wx+jump_taut) ), VOL, skeleton=True ) #splike f
    #  a += SymbolicBFI( timelike 	* ( mean_v*jump_taux + mean_sig*jump_wx + alpha*jump_vx*jump_wx + beta*jump_sigx*jump_taux ), VOL, skeleton=True ) #time like faces
    a += SymbolicBFI( ( pow(c,-2)*v*w + sig*tau ), BND, definedon=fes.mesh.Boundaries("outflow"), skeleton=True) #t=T (or *x)
    a += SymbolicBFI( ( sig*n_x*w + alpha*v*w ), BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True) #dirichlet boundary 'timelike'
    a += SymbolicBFI( ( gamma * (-jump_Ut)*IfPos(n_t,V.Other(),V) ), VOL, skeleton=True ) #correction term to recover sol of second order system
    a += SymbolicBFI( ( gamma * U*V ), BND, definedon=fes.mesh.Boundaries("inflow"), skeleton=True ) #BND correction term to recover sol of second order system
    a.Assemble()
    f = LinearForm(fes)
    f += SymbolicLFI( ( pow(c,-2)*v0*w + sig0*tau ), BND, definedon=fes.mesh.Boundaries("inflow"), skeleton=True) #t=0 (or *(1-x))
    f += SymbolicLFI( ( gD * (alpha*w - tau*n_x) ), BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True) #dirichlet boundary 'timelike'
    f += SymbolicLFI( gamma * ( (U0)*V ), BND, definedon=fes.mesh.Boundaries("inflow"),  skeleton=True ) #rhs correction term to recover sol of second order system
    f.Assemble()

    return [a,f] 


c =2 
order = 9
k = 5

initmesh = Mesh(SegMesh(4,0,1))
mesh = NgsTPmesh(initmesh,c,1)
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





