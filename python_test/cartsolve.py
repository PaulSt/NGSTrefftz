#########################################################################################################################################
N = 9
c=4
order = 6
#########################################################################################################################################
import netgen.meshing as ngm
from netgen.geom2d import unit_square
from ngsolve import *
# # Working with meshes
# ## Two-dimensional meshes
# As example we mesh a unit square [0,1]x[0,1] using quadrilaterals.
#
# We create an empty mesh and inititalize the geometry and the dimension
ngmesh = ngm.Mesh()
ngmesh.SetGeometry(unit_square)
ngmesh.dim = 2
# and add all the `MeshPoint`'s we will need for the final mesh. Similar to the one-dimensional mesh we store the `PointId`'s in the `pnums` array.
pnums = []
for i in range(N + 1):
    for j in range(N + 1):
        pnums.append(ngmesh.Add(ngm.MeshPoint(ngm.Pnt(i / N, j / N, 0))))
# Next, we add the quadrilaterals to the mesh. Before that we have so add a `FaceDescriptor` to the mesh.
# help(ngm.FaceDescriptor)
foo = ngm.FaceDescriptor(surfnr=1,domin=1,bc=1)
ngmesh.Add (foo)
ngmesh.SetMaterial(1, "mat")
for j in range(N):
    for i in range(N):
        ngmesh.Add(ngm.Element2D(1, [pnums[i + j * (N + 1)],
                                 pnums[i + (j + 1) * (N + 1)],
                                 pnums[i + 1 + (j + 1) * (N + 1)],
                                 pnums[i + 1 + j * (N + 1)]]))

# Finally we have to add boundary elements and set boundary conditions.
# horizontal boundaries
for i in range(N):
   ngmesh.Add(ngm.Element1D([pnums[N + i * (N + 1)], pnums[N + (i + 1) * (N + 1)]], index=1))
   ngmesh.Add(ngm.Element1D([pnums[0 + i * (N + 1)], pnums[0 + (i + 1) * (N + 1)]], index=1))
# vertical boundaries
for i in range(N):
   ngmesh.Add(ngm.Element1D([pnums[i], pnums[i + 1]], index=2))
   ngmesh.Add(ngm.Element1D([pnums[i + N * (N + 1)], pnums[i + 1 + N * (N + 1)]], index=2))


mesh = Mesh(ngmesh)

# mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
Draw(mesh)
# print("boundaries" + str(mesh.GetBoundaries()))
#########################################################################################################################################
from ngsolve import *
from trefftzngs import *
import numpy as np
k = 5
fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True)
# fes = L2(mesh, order=order, flags = { "dgjumps" : True })

truesol =  sin( k*(c*x + y) )#exp(-pow(c*x+y,2)))#
U0 = GridFunction(fes)
U0.Set(truesol)
v0 = c*k*cos(k*(c*x+y))#grad(U0)[0]
sig0 = -k*cos(k*(c*x+y))#-grad(U0)[1]
Draw(U0,mesh,'U0')
input()
# Draw(sig0,mesh,'sig0')
# input()
# Draw(v0,mesh,'v0')
# input()

U = fes.TrialFunction()
V = fes.TestFunction()

v = grad(U)[0]
sig = -grad(U)[1]
w = grad(V)[0]
tau = -grad(V)[1]

vo = grad(U.Other())[0]
sigo = -grad(U.Other())[1]
wo = grad(V.Other())[0]
tauo = -grad(V.Other())[1]

h = specialcf.mesh_size
n = specialcf.normal(2)
n_t = n[0]/Norm(n)
n_x = n[1]/Norm(n)

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
a += SymbolicBFI( spacelike * ( gamma * (jump_Ut)*IfPos(n_t,V.Other(),V) ) ,VOL,  skeleton=True ) #correction term to recover sol of second order system
a += SymbolicBFI( spacelike * ( gamma * IfPos(-n_t,1,0) * U*V ) ,BND,  skeleton=True ) #BND correction term to recover sol of second order system
a.Assemble()

f = LinearForm(fes)
f += SymbolicLFI( spacelike * IfPos(-n_t,1,0) *  ( pow(c,-2)*v0*w + sig0*tau ), BND, skeleton=True) #t=0 (or *(1-x))
f += SymbolicLFI( timelike 	* ( v0 * (alpha*w - tau*n_x) ), BND, skeleton=True) #dirichlet boundary 'timelike'
f += SymbolicLFI( spacelike * gamma * IfPos(-n_t,1,0) *  ( (U0)*V ) ,BND,  skeleton=True ) #rhs correction term to recover sol of second order system
f.Assemble()

gfu = GridFunction(fes, name="uDG")
# gfu.vec.data = a.mat.Inverse(inverse="pardiso") * f.vec

nmat = np.zeros((a.mat.height,a.mat.width))
nvec = np.zeros(a.mat.width)

for i in range(a.mat.width):
	#nvec[i] = f.vec[i]
	for j in range(a.mat.height):
		nmat[j,i] = a.mat[j,i]

# sol = np.linalg.solve(nmat,nvec)
# for i in range(a.mat.height):
# 	gfu.vec[i] = sol[i]

nmatclean = nmat[nmat.any(axis=0),:]
nmatclean = nmatclean[:,nmat.any(axis=1)]
nvecclean = f.vec.FV().NumPy()[nmat.any(axis=1)] #nvec

print("cond nmat: ", np.linalg.cond(nmatclean))
# print(nmat)
# print(nvec)

solclean = np.linalg.solve(nmatclean,nvecclean)
sol = np.zeros(a.mat.height)
sol[nmat.any(axis=1)] = solclean
# print(nmat.any(axis=1))


for i in range(a.mat.height):
	gfu.vec[i] = sol[i]
#gfu.vec.data = a.mat.Inverse() * f.vec
# gradu=grad(U0)
# Draw(sig0,mesh,'fun')

print("error=", Integrate((truesol - gfu)*(truesol - gfu), mesh))
print("grad-error=", Integrate((grad(U0) - grad(gfu))*(grad(U0) - grad(gfu)), mesh))
Draw(gfu,mesh,'sol')

Draw(grad(gfu),mesh,'gradsol')

# u = GridFunction(fes,"shapes")
# Draw(u)
# for i in range(2624): #126 #fes.GetNDof()):
#     print("Draw basis function ", i)
#     u.vec[:] = 0
#     u.vec[i] = 1
#     Redraw()
#     input("press key to draw next shape function")
