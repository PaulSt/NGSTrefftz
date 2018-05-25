from netgen.meshing import *
from netgen.geom2d import unit_square
import scipy as sp
import scipy.sparse.linalg

ngmesh = Mesh()
ngmesh.SetGeometry(unit_square)
ngmesh.dim = 2


quads = True
Nx = 10
Ny = 1
pnums = []
for i in range(Nx + 1):
    for j in range(Ny + 1):
        pnums.append(ngmesh.Add(MeshPoint(Pnt(i / Nx, j / Ny, 0))))
        if(j==Ny):
            ngmesh.AddPointIdentification(pnums[-1-Ny],pnums[-1],identnr=1,type=2)
        # if(i==Nx):
        #     ngmesh.AddPointIdentification(pnums[j],pnums[-1],identnr=2,type=2)

ngmesh.Add (FaceDescriptor(surfnr=1,domin=1,bc=1))
ngmesh.SetMaterial(1, "mat")

for j in range(Nx):
    for i in range(Ny):
        if quads:
            ngmesh.Add(Element2D(1, [pnums[i + j * (Ny + 1)],
                                     pnums[i + (j + 1) * (Ny + 1)],
                                     pnums[i + 1 + (j + 1) * (Ny + 1)],
                                     pnums[i + 1 + j * (Ny + 1)]]))
        else:
            ngmesh.Add(Element2D(1, [pnums[i + j * (Ny + 1)],
                                     pnums[i + (j + 1) * (Ny + 1)],
                                     pnums[i + 1 + j * (Ny + 1)]]))
            ngmesh.Add(Element2D(1, [pnums[i + (j + 1) * (Ny + 1)],
                                     pnums[i + 1 + (j + 1) * (Ny + 1)],
                                     pnums[i + 1 + j * (Ny + 1)]]))

# horizontal boundaries
for i in range(Nx):
   ngmesh.Add(Element1D([pnums[Ny + i * (Ny + 1)], pnums[Ny + (i + 1) * (Ny + 1)]], index=1))
   ngmesh.Add(Element1D([pnums[0 + i * (Ny + 1)], pnums[0 + (i + 1) * (Ny + 1)]], index=1))

# vertical boundaries
for i in range(Ny):
   ngmesh.Add(Element1D([pnums[i], pnums[i + 1]], index=2))
   ngmesh.Add(Element1D([pnums[i + Nx * (Ny + 1)], pnums[i + 1 + Nx * (Ny + 1)]], index=2))

from ngsolve.solve import Draw
from ngsolve import *

mesh = Mesh(ngmesh)
Draw(mesh)

fes = L2(mesh, order=5)# , flags = { "dgjumps" : True } )

u = fes.TrialFunction()
v = fes.TestFunction()

# b = CoefficientFunction( (0.2,0.2) )
# b = CoefficientFunction( (y-0.5,0.5-x) )
b = CoefficientFunction( (0,0.2) )
bn = CoefficientFunction(specialcf.normal(2)[1])

a = BilinearForm(fes)
a += SymbolicBFI ( (-u * b*grad(v) ) )
a += SymbolicBFI ( (bn*IfPos(bn, u, u.Other()) * (v-v.Other())), VOL, skeleton=True)
# a += SymbolicBFI ( (bn*IfPos(bn, u, ubnd) * v), BND, skeleton=True)

B=BilinearForm(fes)
B+=SymbolicBFI(u*v)
B.Assemble()
rows,cols,vals = B.mat.COO()
B = sp.sparse.csr_matrix((vals,(rows,cols)))

u = GridFunction(fes)
def peak(x0,y0):
    return exp (-40 * ( (x-x0)*(x-x0) + (y-y0)*(y-y0) ))

u.Set(peak(0.5,0.4))

w = u.vec.CreateVector()

Draw (u, autoscale=False, sd=5)

t = 0
tau = 1e-3
tend = 10

input('start')
# with TaskManager():
while t < tend:
    a.Apply (u.vec, w)
    # fes.SolveM (rho=CoefficientFunction(1), vec=w)
    w.FV().NumPy()[:] = sp.sparse.linalg.spsolve(B,w.FV())
    u.vec.data -= tau * w
    t += tau
    Redraw()
