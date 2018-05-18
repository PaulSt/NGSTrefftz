#from netgen.geom2d import unit_square
from ngsolve import *
from trefftzngs import *
from netgen.csg import unit_cube
from netgen.geom2d import unit_square
import netgen.gui

# mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
# mesh = Mesh("cone/cone.vol.gz")
mesh = Mesh(unit_cube.GenerateMesh(maxh = 1))
Draw(mesh)

c = 1
order = 3
k = 1
fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = c)
# fes = L2(mesh, order=order, dgjumps = True)

gfu = GridFunction(fes)


# u = fes.TestFunction()
# v = fes.TrialFunction()

# bfi = SymbolicBFI(u*v)

# gfu.vec[:] = 0
# for id in mesh.Elements(VOL):
# 	fe = fes.GetFE(id)
# 	trafo = mesh.GetTrafo(id)
# 	mat = bfi.CalcElementMatrix(fe, trafo)
# 	for dof in fes.GetDofNrs(id):
# 		gfu.vec[dof] += 1
# 	inv = Matrix(mat.Height(), mat.Width())
# 	try:
# 		mat.Inverse(inv)
# 	except:
# 		print(id.nr, mat)


	#print(inv)
	#print(Norm(inv))
	# print(dir(bfi))
#print(gfu.vec)
#kx
#uex = sin(kx*x+ky*y - c*t)
uex = sin(k*(c*z+y+x))
# #gfu.Set(CoefficientFunction(0))
gfu.Set(uex)
#
# #print(gfu.vec)
gradu = grad(gfu)
# Draw(uex,mesh,'uex')
# #Draw(gfu, mesh, 'gfu')
print(Integrate((gfu-sin(k*(c*z+y+x)))*(gfu-sin(k*(c*z+y+x))), mesh))
print(Integrate((gradu-CoefficientFunction( tuple( (k*cos(k*(c*z+y+x)), k*cos(k*(c*z+y+x)), c*k*cos(k*(c*z+y+x)))) )) * (gradu-CoefficientFunction( tuple( (k*cos(k*(c*z+y+x)), k*cos(k*(c*z+y+x)), c*k*cos(k*(c*z+y+x))) ))), mesh))

Draw(gfu, mesh, 'gfu',draw_surf=False)
Draw(gradu, mesh, 'gradu',draw_surf=False)

