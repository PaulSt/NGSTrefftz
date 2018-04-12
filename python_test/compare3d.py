#########################################################################################################################################
c = 1
t_steps = 1/4
nt_steps = 1
order = 8
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
SetHeapSize(100*1000*1000)
# mesh = Mesh(unit_cube.GenerateMesh(maxh = 0.2,quad_dominated=True))
ngmeshbase = unit_square.GenerateMesh(maxh = 0.5)
mesh = ProdMesh(ngmeshbase,t_steps)

#Draw(mesh)
#########################################################################################################################################
truesol =  sin( k*(c*z+x+y) )#exp(-pow(c*x+y,2)))#
v0 = c*k*cos(k*(c*z+x+y))#grad(U0)[0]
sig0 = CoefficientFunction( (-k*cos(k*(c*z+x+y)),-k*cos(k*(c*z+x+y))) )#-grad(U0)[1]

def SolveWaveeq(fes,fullsys,mesh):
    [a,f] = DGeqsys(fes,truesol,v0,sig0,c,v0,fullsys)
    [gfu, cond] = DGsolve(fes,a,f)
    L2error = sqrt(Integrate((truesol - gfu)*(truesol - gfu), mesh))
    gradtruesol = CoefficientFunction((-sig0,v0))
    sH1error = sqrt(Integrate((gradtruesol - grad(gfu))*(gradtruesol - grad(gfu)), mesh))
    dof=fes.ndof/mesh.ne
    return [dof,cond,L2error,sH1error]


truesol =  sin( k*(c*z+x+y) )#exp(-pow(c*x+y,2)))#
v0 = c*k*cos(k*(c*z+x+y))#grad(U0)[0]
sig0 = CoefficientFunction( (-k*cos(k*(c*z+x+y)),-k*cos(k*(c*z+x+y))) )#-grad(U0)[1]

ngmeshbase = unit_square.GenerateMesh(maxh = 0.5)
mesh = ProdMesh(ngmeshbase,t_steps)

solution = []

for ordr in range(3,order):
    print("run order: " + str(ordr))

    # btype = 0
    # fes = L2(mesh, order=ordr, dgjumps=True)
    # [dof,cond,L2error,sH1error] = SolveWaveeq(fes, True, mesh)
    # solution.append([btype, ordr, dof,cond,L2error,sH1error])
    # print("btype: " + 'L2' + " dof: " + str(dof) + " cond: " + str(cond) + " L2error: " + str(L2error) + " H1error: "+ str(sH1error))

    btype = 1
    fes = FESpace("trefftzfespace", mesh, order = ordr, wavespeed = c, dgjumps=True, basistype = 0)
    [dof,cond,L2error,sH1error] = SolveWaveeq(fes, False, mesh)
    solution.append([btype,ordr,dof,cond,L2error,sH1error])
    print("btype: " + 'Trefftz' + " dof: " + str(dof) + " cond: " + str(cond) + " L2error: " + str(L2error) + " H1error: "+ str(sH1error) + " run: " + str(ordr))
