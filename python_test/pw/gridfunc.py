
#from netgen.geom2d import unit_square
from ngsolve import *
from trefftzngs import *
from netgen.csg import unit_cube
from netgen.geom2d import unit_square

mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
#mesh = Mesh(unit_cube.GenerateMesh(maxh = 0.4))
#Draw(mesh)

c = 10;
fes = FESpace("trefftzpw", mesh, order = 3, wavespeed = c)
gfu = GridFunction(fes)

uex = exp(c*(x+y))
Draw(uex,mesh,"func")

gfu.Set(uex)
gradu = grad(gfu)
Draw(gfu, mesh, 'gfu')
Draw(gradu, mesh, 'fun')
