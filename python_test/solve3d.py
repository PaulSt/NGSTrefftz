#########################################################################################################################################
N = 4
c=2
t_steps = 1/c
order = 5
k = 5
#########################################################################################################################################
import netgen.meshing as ngm
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve import *

# mesh = Mesh(unit_cube.GenerateMesh(maxh = 0.2,quad_dominated=True))
ngmeshbase = unit_square.GenerateMesh(maxh = 0.2)

ngmesh = ngm.Mesh()
ngmesh.dim=3
pnums = []

for p in ngmeshbase.Points():
	x,y,z = p.p
	ngmesh.Add(ngm.MeshPoint(ngm.Pnt(x,y,0)))
for p in ngmeshbase.Points():
	x,y,z = p.p
	pnums.append(ngmesh.Add(ngm.MeshPoint(ngm.Pnt(x,y,t_steps))))

El1d = ngmeshbase.Elements1D()
El2d = ngmeshbase.Elements2D()

for el in El2d:
	ngmesh.Add(ngm.Element3D(1, [el.points[0].nr,
								el.points[1].nr,
								el.points[2].nr,
								pnums[el.points[0].nr-1],
								pnums[el.points[1].nr-1],
								pnums[el.points[2].nr-1]]))

ngmesh.Add(ngm.FaceDescriptor(surfnr=1,domin=1,bc=1))
for el in El2d:
	ngmesh.Add(ngm.Element2D(1, [pnums[el.points[0].nr-1],
								pnums[el.points[1].nr-1],
								pnums[el.points[2].nr-1]]))
	ngmesh.Add(ngm.Element2D(1, [el.points[0].nr,
								el.points[1].nr,
								el.points[2].nr]))
for el in El1d:
	ngmesh.Add(ngm.Element2D(1, [el.points[0].nr,
								el.points[1].nr,
								pnums[el.points[1].nr-1],
								pnums[el.points[0].nr-1]]))


El3d = ngmesh.Elements3D()


mesh = Mesh(ngmesh)
Draw(mesh)

#
# fes = H1(mesh, order=4)
# u = GridFunction(fes)
# u.Set(x+y+z)
# Draw(u)

#########################################################################################################################################
#
from ngsolve import *
from trefftzngs import *
import numpy as np

truesol =  sin( k*(c*x + y + z) )#exp(-pow(c*x+y,2)))#
v0 = c*k*cos(k*(c*x+y+z))#grad(U0)[0]
sig0 = -k*cos(k*(c*x+y+z))#-grad(U0)[1]

# for order in range(3,order):
fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True, basistype=0)
# fes = L2(mesh, order=order, flags = { "dgjumps" : True })

U0 = GridFunction(fes)
U0.Set(truesol)
Draw(U0)
input()
#
#
#
# U = fes.TrialFunction()
# V = fes.TestFunction()
# gU = grad(U)
# gV = grad(V)
#
# v = gU[0]
# sig = CoefficientFunction(-gU[1],-gU[2])
# w = gV[0]
# tau = CoefficientFunction(-gV[1],-gV[2])
#
# vo = gU.Other()[0]
# sigo = CoefficientFunction(-gU.Other()[1],-gU.Other()[2])
# wo = gV.Other()[0]
# tauo = CoefficientFunction(-gV.Other()[1],-gV.Other()[2])
#
# h = specialcf.mesh_size
# n = specialcf.normal(3)
# n_t = n[0]/Norm(n)
# n_x = CoefficientFunction(n[1],n[2]) /Norm(n)
#
# mean_v = 0.5*(v+vo)
# mean_w = 0.5*(w+wo)
# mean_sig = 0.5*(sig+sigo)
# mean_tau = 0.5*(tau+tauo)
#
# jump_vx = ( v - vo ) * n_x
# jump_wx = ( w - wo ) * n_x
# jump_sigx = (( sig - sigo ), n_x)
# jump_taux = (( tau - tauo ) * n_x)
#
# jump_vt = ( v - vo ) * n_t
# jump_wt = ( w - wo ) * n_t
# jump_sigt = ( sig - sigo ) * n_t
# jump_taut = ( tau - tauo ) * n_t
#
# jump_Ut = (U - U.Other()) * n_t
# #gfu.vec.data = a.mat.Inverse() * f.vec
#
# timelike = n_x*n_x#n_x**2 #IfPos(n_t,0,IfPos(-n_t,0,1)) # n_t=0
# spacelike = n_t**2 #IfPos(n_x,0,IfPos(-n_x,0,1)) # n_x=0
#
# alpha = 0 #pow(10,5)
# beta = 0 #pow(10,5)
# gamma = 1
#
# a = BilinearForm(fes)
# # a += SymbolicBFI( spacelike * ( IfPos(n_t,v,vo)*(pow(c,-2)*jump_wt+jump_taux) + IfPos(n_t,sig,sigo)*(jump_wx+jump_taut) ) ,VOL,  skeleton=True ) #space like faces
# a += SymbolicBFI( spacelike * ( IfPos(n_t,v,vo)*(pow(c,-2)*jump_wt) + IfPos(n_t,sig,sigo)*(jump_taut) ) ,VOL,  skeleton=True ) #space like faces, no jump in x since horizontal
#
# #a += SymbolicBFI( timelike 	* ( mean_v*jump_taux + mean_sig*jump_wx + alpha*jump_vx*jump_wx + beta*jump_sigx*jump_taux ) ,VOL, skeleton=True ) #time like faces
# a += SymbolicBFI( timelike 	* ( mean_v*jump_taux + mean_sig*jump_wx ) ,VOL, skeleton=True ) #time like faces
#
# a += SymbolicBFI( spacelike * IfPos(n_t,1,0) * ( pow(c,-2)*v*w + sig*tau ), BND, skeleton=True) #t=T (or *x)
# a += SymbolicBFI( timelike 	* ( sig*n_x*w + alpha*v*w ), BND, skeleton=True) #dirichlet boundary 'timelike'
# #a += SymbolicBFI( spacelike * ( gamma * (jump_Ut)*IfPos(n_t,V.Other(),V) ) ,VOL,  skeleton=True ) #correction term to recover sol of second order system
# #a += SymbolicBFI( spacelike * ( gamma * IfPos(-n_t,1,0) * U*V ) ,BND,  skeleton=True ) #BND correction term to recover sol of second order system
# a.Assemble()
#
# f = LinearForm(fes)
# f += SymbolicLFI( spacelike * IfPos(-n_t,1,0) *  ( pow(c,-2)*v0*w + sig0*tau ), BND, skeleton=True) #t=0 (or *(1-x))
# f += SymbolicLFI( timelike 	* ( v0 * (alpha*w - tau*n_x) ), BND, skeleton=True) #dirichlet boundary 'timelike'
# #f += SymbolicLFI( spacelike * gamma * IfPos(-n_t,1,0) *  ( (truesol)*V ) ,BND,  skeleton=True ) #rhs correction term to recover sol of second order system
# f.Assemble()

# gfu2= GridFunction(fes, name="uDG")
# gfu2.vec.data = a.mat.Inverse() * f.vec

#
# gfu = GridFunction(fes, name="uDG")
#
# nmat = np.zeros((a.mat.height,a.mat.width))
# nvec = np.zeros(a.mat.width)
#
# for i in range(a.mat.width):#gfu.vec.data = a.mat.Inverse() * f.vec
# 	for j in range(a.mat.height):
# 		nmat[j,i] = a.mat[j,i]
# nvec = f.vec.FV().NumPy() #nvec
#
# sol = np.linalg.solve(nmat,nvec)
# for i in range(a.mat.height):
# 	gfu.vec[i] = sol[i]
# print("cond nmat: ", np.linalg.cond(nmat))
# print(nmat)
#
#
# U0 = GridFunction(fes)
# U0.Set(truesol)
# L2error = Integrate((truesol - gfu)*(truesol - gfu), mesh)
# sH1error = Integrate((grad(U0) - grad(gfu))*(grad(U0) - grad(gfu)), mesh)
# print("error=", L2error)
# print("grad-error=", sH1error)
# Draw(gfu,mesh,'sol')
# Draw(grad(gfu),mesh,'gradsol')
#
#
#
#
#
#
#
#
#
# # nmatclean = nmat[nmat.any(axis=0),:]
# # nmatclean = nmatclean[:,nmat.any(axis=1)]
# # nvecclean = f.vec.FV().NumPy()[nmat.any(axis=1)] #nvec
# # solclean = np.linalg.solve(nmat,nvec)
# # sol = np.zeros(a.mat.height)
# # sol[nmat.any(axis=1)] = solclean
# # print(nmat.any(axis=1))
