from ngsolve import *
from trefftzngs import *
from netgen.geom2d import unit_square
import netgen.gui
mesh=Mesh(unit_square.GenerateMesh(maxh=0.4))
intrule = IntegrationRule(TRIG,2)
test=IntegrationPointFunction(mesh,intrule)

fes = H1(mesh, order=2)
u = fes.TrialFunction()  # symbolic object
v = fes.TestFunction()   # symbolic object
gfu = GridFunction(fes)  # solution

a = BilinearForm(fes)
a += SymbolicBFI(u*v)
a.Assemble()

f = LinearForm(fes)
f += SymbolicLFI(test*v, intrule=intrule)
f.Assemble()
gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs()) * f.vec 
Draw(gfu)

