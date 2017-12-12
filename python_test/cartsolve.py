
import netgen.meshing as ngm
from netgen.geom2d import unit_square
from ngsolve import *
from trefftzngs import *
import numpy as np



# # Working with meshes
# ## Two-dimensional meshes
# As example we mesh a unit square [0,1]x[0,1] using quadrilaterals.
#
# We create an empty mesh and inititalize the geometry and the dimension
ngmesh = ngm.Mesh()
ngmesh.SetGeometry(unit_square)
ngmesh.dim = 2
# and add all the `MeshPoint`'s we will need for the final mesh. Similar to the one-dimensional mesh we store the `PointId`'s in the `pnums` array.
N = 2
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
Draw(mesh)

print("boundaries" + str(mesh.GetBoundaries()))


c=1
order = 3

# fes = L2(mesh, order = order, dgjumps=True)#  dirichlet="default "FESpace("l22", mesh, order = order, dgjumps = True) #
fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True)


incond = sin(x+y*np.pi)
U0 = GridFunction(fes)
U0.Set(incond)
v0 = grad(U0)[0]
sig0 = -grad(U0)[1]
# Draw(U0,mesh,'U0')
# Draw(sig0,mesh,'sig0')
# Draw(v0,mesh,'v0')

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
n_t = n[0]
n_x = n[1]

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


a = BilinearForm(fes)
# a += SymbolicBFI( (n_t!=0) * ( pow(c,-2)*IfPos(n_t,v,vo)*(jump_wt+jump_taux) + IfPos(n_t,sig,sigo)*(jump_wt+jump_taux) ) , skeleton=True ) #space like faces
# a += SymbolicBFI( (n_x!=0) * ( mean_v*jump_taux + mean_sig*jump_wx + 0.5*jump_vx*jump_wx + 0.5*jump_sigx*jump_taux ) , skeleton=True ) #time like faces
a += SymbolicBFI( (n_x!=0) * (sig + 0.5*v*n_x)*(w*n_x+tau*n_t), BND, skeleton=True) #dirichlet boundary 'timelike' new
# a += SymbolicBFI( (n_t==1) * ( pow(c,-2)*v*w + sig*tau ), BND, skeleton=True) #t=T
# a += SymbolicBFI( (n_x!=0) * ( sig*n_x*w + 0.5*v*w ), BND, skeleton=True) #dirichlet boundary 'timelike'
a.Assemble()

f = LinearForm(fes)
f += SymbolicLFI( IfPos(-n_t, 1, 0) * ( pow(c,-2)*v0*w + sig0*tau ), BND, skeleton=True) #t=0
f.Assemble()



offset = 0
nmat = np.zeros((a.mat.height-offset,a.mat.width-offset))
nvec = np.zeros(a.mat.width-offset)
for i in range(a.mat.width-offset):
	nvec[i] = f.vec[i+offset]
	for j in range(a.mat.height-offset):
		nmat[j,i] = a.mat[j+offset,i+offset]

nmatclean = nmat[:,nmat.any(axis=1)]
nmatclean = nmatclean[nmat.any(axis=1),:]
nvecclean = nvec[nmat.any(axis=1)]
solclean = np.linalg.solve(nmatclean,nvecclean)
sol = np.zeros(a.mat.height)
sol[nmat.any(axis=1)] = solclean



gfu = GridFunction(fes, name="uDG")
for i in range(a.mat.height-offset):
	gfu.vec[i+offset] = sol[i]
# gfu.vec.data = a.mat.Inverse() * f.vec
# gradu=grad(U0)
# Draw(sig0,mesh,'fun')
Draw(gfu)



# u = GridFunction(fes,"shapes")
# Draw(u)
# for i in range(2624): #126 #fes.GetNDof()):
#     print("Draw basis function ", i)
#     u.vec[:] = 0
#     u.vec[i] = 1
#     Redraw()
#     input("press key to draw next shape function")
