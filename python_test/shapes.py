# Test the new exported coefficientfunction

from netgen.geom2d import unit_square
from ngsolve import *
from trefftzngs import *

from netgen.csg import unit_cube

#mesh = Mesh(unit_square.GenerateMesh(maxh=0.4))
mesh = Mesh(unit_cube.GenerateMesh(maxh = 0.4))

#fes = FESpace("trefftzfespace", mesh, flags = { "dgjumps" : True, "order" : 3 })
fes = FESpace("trefftzfespace", mesh, order= 3)

u = GridFunction(fes,"shapes")

Draw(u)

#def printshape(i):
#    print("Draw basis function ", i)
#    u.vec[:] = 0
#    u.vec[i] = 1
#    Redraw()

#print("ndof = ", fes.ndof)
#for i in range(len(u.vec)):
#    printshape(int(input("enter number of shapefunction to print:")))


# we can use the additionally exported function here
for i in range(2624): #126 #fes.GetNDof()):
    print("Draw basis function ", i)
    u.vec[:] = 0
    u.vec[i] = 1
    Redraw()
    input("press key to draw next shape function")
