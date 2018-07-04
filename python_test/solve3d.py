#########################################################################################################################################
c = 2
basemeshsize = 1/5
t_stepsize = 1/3
order = 4
#########################################################################################################################################
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from prodmesh import *
from ngsolve import *
import netgen.gui
from trefftzngs import *
from DGeq import *
import time
from ngsolve.TensorProductTools import *
from ngsolve.comp import *
ngmeshbase = unit_square.GenerateMesh(maxh = basemeshsize)
mesh = PeriodicProdMesh(ngmeshbase,t_stepsize)

# mesh1 = Mesh(SegMesh(4,0,0.5))
# mesh2 = Mesh(unit_square.GenerateMesh(maxh=0.4))
# mesh = Mesh(MakeTensorProductMesh(mesh2,mesh1))
Draw(mesh)

#########################################################################################################################################
sq = sqrt(0.5)
truesol = sin( c*z+sq*(x+y) )
v0 = c*cos(c*z+sq*(x+y))
sig0 = CoefficientFunction(( -sq*cos(c*z+sq*(x+y)),-sq*cos(c*z+sq*(x+y)) ))
# truesol = exp(-100*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)) )
# v0 = 0
# sig0 = CoefficientFunction( (200 * (x-0.5) * truesol, 200 * (y-0.5) * truesol))
# truesol = -(sqrt(0.5)*x+sqrt(0.5)*y-z)*(sqrt(0.5)*x+sqrt(0.5)*y-z)
# v0 = 2*(sqrt(0.5)*x+sqrt(0.5)*y-z)
# sig0 = CoefficientFunction( (2*(sqrt(0.5)*x+sqrt(0.5)*y-z)*sqrt(0.5),2*sqrt(0.5)*(sqrt(0.5)*x+sqrt(0.5)*y-z)) )
gradtruesol = CoefficientFunction( (-sig0[0],-sig0[1], v0) )

# fes = Periodic(FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True, basistype=0))
# fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True)
fes = L2(mesh, order=order, dgjumps = True )

gfu = GridFunction(fes)
gfu.Set(truesol) # gD = c*k* cos( k*(c*z+x+y) )
gD = v0
# gD = grad(gfu)[2]

for t in range(1,5):
    gD = c*cos(c*(z+t*t_stepsize)+sq*(x+y))
    start = time.clock()
    [a,f] = DGeqsys(fes,gfu,0,0,c,gD,True,True)
    print("DGsys: ", str(time.clock()-start))

    start = time.clock()
    [gfu, cond] = DGsolve(fes,a,f)
    print("DGsolve: ", str(time.clock()-start))

    truesol = sin(c*(z+t*t_stepsize)+sq*(x+y))
    gradtruesol = CoefficientFunction( (sq*cos(c*(z+t*t_stepsize)+sq*(x+y)),sq*cos(c*(z+t*t_stepsize)+sq*(x+y)),c*cos(c*(z+t*t_stepsize)+sq*(x+y)) ) )
    L2error = sqrt(Integrate((truesol - gfu)*(truesol - gfu), mesh))
    sH1error = sqrt(Integrate((gradtruesol - grad(gfu))*(gradtruesol - grad(gfu)), mesh))
    print("L2error=", L2error)
    print("grad-error=", sH1error)

    Draw(gfu,mesh,'sol',draw_surf=False)
    # Draw(truesol, mesh, 'truesol', draw_surf=False)
    # Draw(grad(gfu),mesh,'gsol',draw_surf=False)
    input(t)
