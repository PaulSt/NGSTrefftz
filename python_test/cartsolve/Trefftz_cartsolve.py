#########################################################################################################################################
N = 4
c=2
t_steps = c*N
order = 5
k = 5
#########################################################################################################################################
from trefftzngs import *
import numpy as np

import netgen.meshing as ngm
from netgen.geom2d import unit_square
from ngsolve import *

ngmesh = ngm.Mesh()
ngmesh.SetGeometry(unit_square)
ngmesh.dim = 2
pnums = []
for j in range(t_steps + 1):
	for i in range(N + 1):
		pnums.append(ngmesh.Add(ngm.MeshPoint(ngm.Pnt(i / N, j / t_steps, 0))))

foo = ngm.FaceDescriptor(surfnr=1,domin=1,bc=1)
ngmesh.Add (foo)
ngmesh.SetMaterial(1, "mat")
for j in range(t_steps):
	for i in range(N):
		ngmesh.Add(ngm.Element2D(1, [pnums[i + j * (N + 1)],
									pnums[i + 1 + j * (N + 1)],
									pnums[i + 1 + (j + 1) * (N + 1)],
									pnums[i + (j + 1) * (N + 1)]]))
for i in range(t_steps):
	ngmesh.Add(ngm.Element1D([pnums[N + i * (N + 1)], pnums[N + (i + 1) * (N + 1)]], index=1))
	ngmesh.Add(ngm.Element1D([pnums[0 + i * (N + 1)], pnums[0 + (i + 1) * (N + 1)]], index=1))
for i in range(N):
	ngmesh.Add(ngm.Element1D([pnums[i], pnums[i + 1]], index=2))
	ngmesh.Add(ngm.Element1D([pnums[i + t_steps * (N + 1)], pnums[i + 1 + t_steps * (N + 1)]], index=2))

mesh = Mesh(ngmesh)

Draw(mesh)
# print("boundaries" + str(mesh.GetBoundaries()))
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
input()

U = fes.TrialFunction()
V = fes.TestFunction()

v = grad(U)[1]
sig = -grad(U)[0]
w = grad(V)[1]
tau = -grad(V)[0]

vo = grad(U.Other())[1]
sigo = -grad(U.Other())[0]
wo = grad(V.Other())[1]
tauo = -grad(V.Other())[0]

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

jump_Ut = (U - U.Other()) * n_t
#gfu.vec.data = a.mat.Inverse() * f.vec

timelike = n_x**2 #IfPos(n_t,0,IfPos(-n_t,0,1)) # n_t=0
spacelike = n_t**2 #IfPos(n_x,0,IfPos(-n_x,0,1)) # n_x=0

alpha = 0.5 #pow(10,5)
beta = 0.5 #pow(10,5)
gamma = 1

a = BilinearForm(fes)
# a += SymbolicBFI( spacelike * ( IfPos(n_t,v,vo)*(pow(c,-2)*jump_wt+jump_taux) + IfPos(n_t,sig,sigo)*(jump_wx+jump_taut) ) ,VOL,  skeleton=True ) #space like faces
a += SymbolicBFI( spacelike * ( IfPos(n_t,v,vo)*(pow(c,-2)*jump_wt) + IfPos(n_t,sig,sigo)*(jump_taut) ) ,VOL,  skeleton=True ) #space like faces, no jump in x since horizontal
a += SymbolicBFI( timelike 	* ( mean_v*jump_taux + mean_sig*jump_wx + alpha*jump_vx*jump_wx + beta*jump_sigx*jump_taux ) ,VOL, skeleton=True ) #time like faces
a += SymbolicBFI( spacelike * IfPos(n_t,1,0) * ( pow(c,-2)*v*w + sig*tau ), BND, skeleton=True) #t=T (or *x)
a += SymbolicBFI( timelike 	* ( sig*n_x*w + alpha*v*w ), BND, skeleton=True) #dirichlet boundary 'timelike'
a += SymbolicBFI( spacelike * ( gamma * (-jump_Ut)*IfPos(n_t,V.Other(),V) ) ,VOL,  skeleton=True ) #correction term to recover sol of second order system
a += SymbolicBFI( spacelike * ( gamma * IfPos(-n_t,1,0) * U*V ) ,BND,  skeleton=True ) #BND correction term to recover sol of second order system
a.Assemble()

f = LinearForm(fes)
f += SymbolicLFI( spacelike * IfPos(-n_t,1,0) *  ( pow(c,-2)*v0*w + sig0*tau ), BND, skeleton=True) #t=0 (or *(1-x))
f += SymbolicLFI( timelike 	* ( v0 * (alpha*w - tau*n_x) ), BND, skeleton=True) #dirichlet boundary 'timelike'
f += SymbolicLFI( spacelike * gamma * IfPos(-n_t,1,0) *  ( (truesol)*V ) ,BND,  skeleton=True ) #rhs correction term to recover sol of second order system
f.Assemble()

# gfu2= GridFunction(fes, name="uDG")
# gfu2.vec.data = a.mat.Inverse() * f.vec


gfu = GridFunction(fes, name="uDG")

nmat = np.zeros((a.mat.height,a.mat.width))
nvec = np.zeros(a.mat.width)

for i in range(a.mat.width):#gfu.vec.data = a.mat.Inverse() * f.vec
	for j in range(a.mat.height):
		nmat[j,i] = a.mat[j,i]
nvec = f.vec.FV().NumPy() #nvec

sol = np.linalg.solve(nmat,nvec)
for i in range(a.mat.height):
	gfu.vec[i] = sol[i]
print("cond nmat: ", np.linalg.cond(nmat))
# print(nmat)


U0 = GridFunction(fes)
U0.Set(truesol)
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
