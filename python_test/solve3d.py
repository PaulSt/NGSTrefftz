#########################################################################################################################################
c = 1
t_steps = 1/5
nt_steps = 1
order = 5
k = 1
#########################################################################################################################################
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from prodmesh import ProdMesh
from ngsolve import *
# import netgen.gui
from trefftzngs import *
from DGeq import *
import time
# mesh = Mesh(unit_cube.GenerateMesh(maxh = 0.2,quad_dominated=True))
ngmeshbase = unit_square.GenerateMesh(maxh = 0.3)
mesh = ProdMesh(ngmeshbase,t_steps)
#Draw(mesh)
#########################################################################################################################################
truesol =  sin( k*(c*z+x+y) )#exp(-pow(c*x+y,2)))#
v0 = c*k*cos(k*(c*z+x+y))#grad(U0)[0]
sig0 = CoefficientFunction( (-k*cos(k*(c*z+x+y)),-k*cos(k*(c*z+x+y))) )#-grad(U0)[1]
# for order in range(3,order):
fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True, basistype=0)
# fes = L2(mesh, order=order, flags = { "dgjumps" : True })
# U0 = GridFunction(fes)
# U0.Set(truesol)basemesh = Mesh(ngmeshbase)
# Draw(U0)
# input()






basemesh = Mesh(ngmeshbase)
Draw(basemesh)
gfu = GridFunction(fes)
gfu.Set(truesol)
#input("startcomp")
for t_step in range(nt_steps):
	truesol = sin( k*(c* (z+t_step*t_steps) +x+y) )
	gD = c*k* cos( k*(c* (z+t_step*t_steps) +x+y) )
	v0 = c*k*cos(k*(c*(z+t_step*t_steps)+x+y))
	sig0 = CoefficientFunction( (-k*cos(k*(c*(z+t_step*t_steps)+x+y)),-k*cos(k*(c*(z+t_step*t_steps)+x+y))) )
	start = time.clock()
	[a,f] = DGeqsys(fes,truesol,v0,sig0,c,gD)
	print("DGsys: ", str(time.clock()-start))

	start = time.clock()
	[gfu, cond] = DGsolve(fes,a,f)
	print("DGsolve: ", str(time.clock()-start))

	#sig0 = CoefficientFunction((-grad(gfu)[0],-grad(gfu)[1]))
	#v0 = grad(gfu)[2]
	L2error = Integrate((truesol - gfu)*(truesol - gfu), mesh)
	gradtruesol = CoefficientFunction(( k*cos(k*(c*z+x+y)), k*cos(k*(c*z+x+y)), c*k*cos(k*(c*z+x+y)) ))
	sH1error = sqrt(Integrate((gradtruesol - grad(gfu))*(gradtruesol - grad(gfu)), mesh))
	print("L2error=", L2error)
	print("grad-error=", sH1error)

	#Draw(gfu,basemesh,'sol')
	#input(t_step)

# Draw(gfu,mesh,'sol')
# Draw(grad(gfu),mesh,'gradsol')

#
#Draw(gfu,basemesh,'sol')

#for time in np.arange(0.0, t_steps, 0.1):
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
