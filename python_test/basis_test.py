
#from netgen.geom2d import unit_square
from ngsolve import *
from myngspy import *
from trefftzngs import *
from netgen.csg import unit_cube

mesh = Mesh(unit_cube.GenerateMesh(maxh = 0.2))
#Draw(mesh)
fes = FESpace("trefftzfespace", mesh, order = 3) #, flags = {"wavespeed":17})
gfu = GridFunction(fes)
#kx
#uex = sin(kx*x+ky*y - c*t)
uex = x
gfu.Set(uex)
Draw(uex)
