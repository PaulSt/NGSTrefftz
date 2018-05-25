from netgen.meshing import *
from netgen.geom2d import unit_square
import scipy as sp
import scipy.sparse.linalg
from ngsolve.solve import Draw
from ngsolve import *

#########################################################################################################################################
from prodmesh import *
import netgen.gui
from DGeq import *

ngmeshbase = unit_square.GenerateMesh(maxh = 1/4)
mesh = PeriodicProdMesh(ngmeshbase,0.5)
Draw(mesh)

fes = L2(mesh, order=3)# , flags = { "dgjumps" : True } )

u = fes.TrialFunction()
v = fes.TestFunction()

# b = CoefficientFunction( (0.2,0.2) )
# b = CoefficientFunction( (y-0.5,0.5-x) )
b = CoefficientFunction( (0,0,0.2) )
bn = CoefficientFunction(specialcf.normal(3)*b)

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
def peak(x0,y0,z0):
    return exp (-40 * ( (x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0) ))

u.Set(peak(0.5,0.4,0.5))

w = u.vec.CreateVector()

Draw (u, autoscale=False, sd=5)

t = 0
tau = 1e-3
tend = 10

input('start')
# with TaskManager():
while t < tend:
    a.Apply (u.vec, w)
    fes.SolveM (rho=CoefficientFunction(1), vec=w)
    # w.FV().NumPy()[:] = sp.sparse.linalg.spsolve(B,w.FV())
    u.vec.data -= tau * w
    t += tau
    Redraw()
