# Test the new exported coefficientfunction

from netgen.geom2d import unit_square
from ngsolve import *
from trefftzngs import *

from netgen.csg import unit_cube

mesh = Mesh(unit_cube.GenerateMesh(maxh = 0.2))

#mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
cf = TrefftzCoefficient(0)
Draw(cf,mesh,"trefftzCoefficient")
i=1
while (i<16):
	print("drawing..")
	cf = TrefftzCoefficient(i)
	Redraw()
	i=i+1
	input("..finished -  Press Enter")
