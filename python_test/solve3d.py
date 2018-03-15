#########################################################################################################################################
N = 4
c=2
t_steps = 1/5
order = 5
k = 1
#########################################################################################################################################
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from prodmesh import ProdMesh
from ngsolve import *
import netgen.gui
from trefftzngs import *
from DGeq import *

# mesh = Mesh(unit_cube.GenerateMesh(maxh = 0.2,quad_dominated=True))
ngmeshbase = unit_square.GenerateMesh(maxh = 0.4)
mesh = ProdMesh(ngmeshbase,t_steps)
Draw(mesh)
#########################################################################################################################################
truesol =  sin( k*(c*z+x+y) )#exp(-pow(c*x+y,2)))#
v0 = c*k*cos(k*(c*z+x+y))#grad(U0)[0]
sig0 = CoefficientFunction( (-k*cos(k*(c*z+x+y)),-k*cos(k*(c*z+x+y))) )#-grad(U0)[1]

# for order in range(3,order):
fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True, basistype=0)
# fes = L2(mesh, order=order, flags = { "dgjumps" : True })
U0 = GridFunction(fes)
U0.Set(truesol)
# Draw(U0)
# input()

[a,f] = DGeqsys(fes,truesol,v0,sig0,c)

gfu = DGsolve(fes,a,f)

L2error = Integrate((truesol - gfu)*(truesol - gfu), mesh)
sH1error = Integrate((grad(U0) - grad(gfu))*(grad(U0) - grad(gfu)), mesh)
print("error=", L2error)
print("grad-error=", sH1error)
# Draw(gfu,mesh,'sol')
# Draw(grad(gfu),mesh,'gradsol')

#basemesh = Mesh(ngmeshbase)
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
