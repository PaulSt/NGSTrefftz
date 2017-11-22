
import netgen.meshing as ngm
from netgen.geom2d import unit_square
from ngsolve import *
from trefftzngs import *


# # Working with meshes
# ## Two-dimensional meshes
# As example we mesh a unit square [0,1]x[0,1] using quadrilaterals.
#
# We create an empty mesh and inititalize the geometry and the dimension
ngmesh = ngm.Mesh()
ngmesh.SetGeometry(unit_square)
ngmesh.dim = 2
# and add all the `MeshPoint`'s we will need for the final mesh. Similar to the one-dimensional mesh we store the `PointId`'s in the `pnums` array.
N = 5
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


c=10
V = FESpace("trefftzfespace", mesh, order = 3, wavespeed = c, dirichlet="default")#H1(mesh, dirichlet="default")
print(V.FreeDofs())
gfu = GridFunction(V)
g = x
gfu.Set(g,BND)
Draw(gfu)

# a = BilinearForm ( V );
# a += SymbolicBFI ( bn*IfPos(bn, u, u.Other()) * (v-v.Other()), VOL, skeleton=True)
# a += SymbolicBFI ( bn*IfPos(bn, u, ubnd) * v, BND, skeleton=True)
