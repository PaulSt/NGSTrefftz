
#from netgen.geom2d import unit_square
from ngsolve import *
from trefftzngs import *
from myngspy import *
from netgen.csg import unit_cube
from netgen.geom2d import unit_square

mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
#mesh = Mesh(unit_cube.GenerateMesh(maxh = 0.4))
#Draw(mesh)

c = 1;
fes = FESpace("trefftzfespace", mesh, dirichlet="top|bottom|right|left", order = 3, wavespeed = c, dgjumps = True)#mesh, order = 3, wavespeed = c)

print ("freedofs: ", fes.FreeDofs())

u = fes.TrialFunction()
v = fes.TestFunction()

a = BilinearForm(fes)
a += SymbolicBFI(grad(u)*grad(v))

f = LinearForm(fes)
f += SymbolicLFI(x*v)

a.Assemble()
f.Assemble()

u = GridFunction(fes)

print ("solve")
u.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

Draw(u)
