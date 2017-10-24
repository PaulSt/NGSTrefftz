
#from netgen.geom2d import unit_square
from ngsolve import *
from trefftzngs import *
from netgen.csg import unit_cube
from netgen.geom2d import unit_square

mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
#mesh = Mesh(unit_cube.GenerateMesh(maxh = 0.4))
#Draw(mesh)

c = 1;
fes = FESpace("trefftzfespace", mesh, order = 3, wavespeed = c)
gfu = GridFunction(fes)
#kx
#uex = sin(kx*x+ky*y - c*t)
uex = sin(c*x+y)
gfu.Set(uex)
gradu = grad(gfu)
#Draw(gfu)
Draw(gradu, mesh, 'fun')
