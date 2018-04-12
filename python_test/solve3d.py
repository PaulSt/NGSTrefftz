#########################################################################################################################################
c = 1
t_stepsize = 1/5
nt_steps = 3
order = 4
k = 1
#########################################################################################################################################
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from prodmesh import ProdMesh
from ngsolve import *
import netgen.gui
from trefftzngs import *
from DGeq import *
import time
# mesh = Mesh(unit_cube.GenerateMesh(maxh = 0.2,quad_dominated=True))
ngmeshbase = unit_square.GenerateMesh(maxh = t_stepsize)
mesh = ProdMesh(ngmeshbase,t_stepsize)
#Draw(mesh)
#########################################################################################################################################
truesol =  sin( k*(c*z+x+y) )#exp(-pow(c*x+y,2)))#
v0 = c*k*cos(k*(c*z+x+y))#grad(U0)[0]
sig0 = CoefficientFunction( (-k*cos(k*(c*z+x+y)),-k*cos(k*(c*z+x+y))) )#-grad(U0)[1]

fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True, basistype=0)
# fes = L2(mesh, order=order, flags = { "dgjumps" : True })

basemesh = Mesh(ngmeshbase)
Draw(basemesh)
gfu = GridFunction(fes)
gfu.Set(truesol)
#input("startcomp")
for t in range(nt_steps):
	truesol = sin( k*( c*(z+t*t_stepsize) +x+y) )
	gD = c*k* cos( k*( c*(z+t*t_stepsize) +x+y) )
	# v0 = c*k*cos(k*(c*(z+t*t_stepsize)+x+y))
	# sig0 = CoefficientFunction( (-k*cos(k*(c*(z+t*t_stepsize)+x+y)),-k*cos(k*(c*(z+t*t_stepsize)+x+y))) )

	start = time.clock()
	[a,f] = DGeqsys(fes,gfu,v0,sig0,c,gD)
	print("DGsys: ", str(time.clock()-start))

	start = time.clock()
	[gfu, cond] = DGsolve(fes,a,f)
	print("DGsolve: ", str(time.clock()-start))

	v0 = ClipCoefficientFunction(grad(gfu)[2], 2, t_stepsize)
	sig0 = CoefficientFunction( (ClipCoefficientFunction((-grad(gfu)[0]), 2, t_stepsize), ClipCoefficientFunction((-grad(gfu)[1]), 2, t_stepsize)) )

	L2error = sqrt(Integrate((truesol - gfu)*(truesol - gfu), mesh))
	gradtruesol = CoefficientFunction(( k*cos(k*(c*z+x+y)), k*cos(k*(c*z+x+y)), c*k*cos(k*(c*z+x+y)) ))
	sH1error = sqrt(Integrate((gradtruesol - grad(gfu))*(gradtruesol - grad(gfu)), mesh))
	print("L2error=", L2error)
	print("grad-error=", sH1error)

	# Draw(grad(gfu)[2],mesh,'grad')
	# Draw(gD,mesh,'gD')
	# input()
	# Draw(v0,basemesh,'gradtime')
	gfu = ClipCoefficientFunction(gfu,2,t_stepsize)
	Draw(gfu,basemesh,'sol')
	input(t)

# Draw(gfu,mesh,'sol')
# Draw(grad(gfu),mesh,'gradsol')

#
#Draw(gfu,basemesh,'sol')

#for time in np.arange(0.0, t_stepsize, 0.1):
#l2fes = L2(basemesh,order=order)
#foo = GridFunction(l2fes)

#for p in ngmeshbase.Points():
#	gfu(p.p)



#
#
# # nmatclean = nmat[nmat.any(axis=0),:]
# # nmatclean = nmatclean[:,nmat.any(axis=1)]
# # nvecclean = f.vec.FV().NumPy()[nmat.any(axis=1)] #nvec
# # solclean = np.linalg.solve(nmat,nvec)
# # sol = np.zeros(a.mat.height)
# # sol[nmat.any(axis=1)] = solclean
# # print(nmat.any(axis=1))
