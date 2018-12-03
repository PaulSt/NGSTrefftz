# Test the new exported coefficientfunction
import netgen.gui
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from trefftzngs import *
from ngsolve.TensorProductTools import *
from ngsolve import *
import sys
sys.path.append("..")
from prodmesh import *

import netgen.gui
# mesh = Mesh(unit_square.GenerateMesh(maxh=0.4))
# mesh = Mesh(unit_cube.GenerateMesh(maxh = 0.41))
ngmesh = SegMesh(2,0,1)
# ngmesh = QadSegMesh(2,0,1)
initmesh = Mesh(ngmesh)

mesh = NgsTPmesh(initmesh,2,0.6)
Draw(mesh)

order = 9
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
