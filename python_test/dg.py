from netgen.geom2d import unit_square
from netgen.csg import unit_cube

from ngsolve import *
from trefftzngs import *

mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
# mesh = Mesh(unit_cube.GenerateMesh(maxh = 0.4))

c = 10
order = 3
# fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps = True)
fes = FESpace("l2", mesh, order = order, dgjumps = True) #L2(mesh, order=order, flags = { "dgjumps" : True })
u = fes.TrialFunction()
v = fes.TestFunction()

jump_u = u-u.Other()
jump_v = v-v.Other()
n = specialcf.normal(2)
mean_dudn = 0.5*n * (grad(u)+grad(u.Other()))
mean_dvdn = 0.5*n * (grad(v)+grad(v.Other()))


alpha = 4
h = specialcf.mesh_size
a = BilinearForm(fes)
a += SymbolicBFI(grad(u)*grad(v))
a += SymbolicBFI(alpha*order**2/h*jump_u*jump_v, skeleton=True)
a += SymbolicBFI(alpha*order**2/h*u*v, BND, skeleton=True)
# a += SymbolicBFI(-mean_dudn*jump_v -mean_dvdn*jump_u, skeleton=True)
# a += SymbolicBFI(-n*grad(u)*v-n*grad(v)*u, BND, skeleton=True)
a.Assemble()

f = LinearForm(fes)
f += SymbolicLFI((x+y)*v)
f.Assemble()

gfu = GridFunction(fes, name="uDG")
gfu.vec.data = a.mat.Inverse() * f.vec
Draw (gfu)
#
# funk = GridFunction(fes)
# funk.Set(exp(-x*y))
# Draw(funk)

fes.FreeDofs()
