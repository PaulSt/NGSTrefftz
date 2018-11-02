# Test the new exported coefficientfunction
import netgen.gui
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve import *
from trefftzngs import *

# mesh = Mesh(unit_square.GenerateMesh(maxh=0.4))
mesh = Mesh(unit_cube.GenerateMesh(maxh = 0.41))

order = 4
fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = 1, dgjumps=True)
# fes = FESpace("l2", mesh, order = order, dgjumps = True) #L2(mesh, order=order, flags = { "dgjumps" : True })

u = GridFunction(fes,"shapes")
u.vec[:]=0
Draw(u)

for i in range(u.vec.size):
    print("Draw basis function ", i)
    u.vec[:] = 0
    u.vec[i] = 1
    Redraw()
    input("press key to draw next shape function")
