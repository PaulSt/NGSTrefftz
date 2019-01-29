from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from trefftzngs import *
import netgen.gui
from ngsolve.TensorProductTools import *
from ngsolve import *
from prodmesh import *
from ngsolve.solve import Tcl_Eval # for snapshots
from testcases import *

order = 3
t_start = 0
t_step = 0.05

geo = SplineGeometry()
# geo.AddRectangle((0,0), (1,1),
                 # bcs=["b","r","t","l"],
                 # leftdomain=1)
# geo.AddRectangle((0.5,0.5), (0.75,0.75),
                 # bcs=["b2","r2","t2","l2"],
                 # leftdomain=2, rightdomain=1)

p1,p2,p3,p4 = [ geo.AppendPoint(x,y) for x,y in [(0,0), (1,0), (1,1), (0,1)] ]
p5,p6 =  [ geo.AppendPoint(x,y) for x,y in [(2,0), (2,1)] ]
geo.Append (["line", p1, p2], leftdomain=1, rightdomain=0)
geo.Append (["line", p3, p2], leftdomain=2, rightdomain=1) #got to be careful with the normal vector along interfaces
geo.Append (["line", p3, p4], leftdomain=1, rightdomain=0)
geo.Append (["line", p4, p1], leftdomain=1, rightdomain=0)
geo.Append (["line", p2, p5], leftdomain=2, rightdomain=0)
geo.Append (["line", p5, p6], leftdomain=2, rightdomain=0)
geo.Append (["line", p6, p3], leftdomain=2, rightdomain=0)

geo.SetMaterial (1, "outer")
geo.SetMaterial (2, "inner")
initmesh = Mesh(geo.GenerateMesh(maxh=0.025))
for i in range(0,len(initmesh.GetBoundaries())):
   initmesh.ngmesh.SetBCName(i,"neumann")
Draw(initmesh)

D = initmesh.dim
t = CoordCF(D)
c = Vector(2)
c[0] = 1
c[1] = 0.5

cf = CoefficientFunction(list(c))
sq = sqrt(0.5);
bdd = CoefficientFunction((
    sin( cf*t+sq*(x+y) ),
    sq*cos(cf*t+sq*(x+y)),
    sq*cos(cf*t+sq*(x+y)),
    cf*cos(cf*t+sq*(x+y))
    ))
bdd = vertgausspw(D,1)
Draw(bdd, initmesh, 'piecewise')
wavefront = EvolveTentsMakeWavefront(order,initmesh,t_start,bdd)

D = initmesh.dim
if D==3: eltyp = ET.TET
elif D==2: eltyp = ET.TRIG
elif D==1: eltyp = ET.SEGM
intrule = IntegrationRule(eltyp,2*order)
irsize = len(intrule.points)

fes = H1(initmesh, order=order)
u,v = fes.TnT()
gfu = GridFunction(fes)
a = BilinearForm(fes)
a += SymbolicBFI(u*v)
a.Assemble()

ipfct=IntegrationPointFunction(initmesh,intrule,wavefront)
f = LinearForm(fes)
f += SymbolicLFI(ipfct*v, intrule=intrule)
f.Assemble()
gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs()) * f.vec
# Draw(gfu,initmesh,'sol')
Draw(gfu,initmesh,'sol',autoscale=False,min=-0.1,max=0.35)
input()
filename = "results/mov/sol_mat0.jpg"
# Tcl_Eval("Ng_SnapShot .ndraw {};\n".format(filename))

with TaskManager():
    for t in range(0,50):
        wavefront = EvolveTents(order,initmesh,c,t_step,wavefront,t_start,bdd)
        print("Error: " + str(EvolveTentsError(order,initmesh,wavefront,EvolveTentsMakeWavefront(order,initmesh,t_start + t_step,bdd))))

        ipfct=IntegrationPointFunction(initmesh,intrule,wavefront)
        f = LinearForm(fes)
        f += SymbolicLFI(ipfct*v, intrule=intrule)
        f.Assemble()
        gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs()) * f.vec
        Redraw(blocking=True)

        t_start += t_step
        print("time: " + str(t_start))
        # filename = "results/mov/sol_mat"+str(t).zfill(3) +".jpg"
        # Tcl_Eval("Ng_SnapShot .ndraw {};\n".format(filename))
