# Test the new exported coefficientfunction

from netgen.geom2d import unit_square
from ngsolve import *
from trefftzngs import *

from netgen.csg import unit_cube

mesh = Mesh(unit_cube.GenerateMesh(maxh = 0.2))

#i=0
#while (i<16):
#	print("drawing..")
#	cf = TrefftzCoefficient(i)
#	Draw(cf,mesh,"trefftzCoefficient")
#	i=i+1
#	input("..finished -  Press Enter")

fes = FESpace("trefftzfespace", mesh, order=3)

u = GridFunction(fes,"shapes")

Draw(u)

def printshape(i):
    print("Draw basis function ", i)
    u.vec[:] = 0
    u.vec[i] = 1
    Redraw()

print("ndof = ", fes.ndof)
for i in range(len(u.vec)):
    printshape(int(input("enter number of shapefunction to print:")))
