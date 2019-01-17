from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from trefftzngs import *
import netgen.gui
from ngsolve.TensorProductTools import *
from ngsolve import *
from prodmesh import *
from ngsolve.solve import Tcl_Eval # for snapshots
from testcases import *

order = 2
t_start = 0
t_step = 0.1

# for i in range(0,len(initmesh.GetBoundaries())):
   # initmesh.ngmesh.SetBCName(i,"neumann")
geo = SplineGeometry()
geo.AddRectangle((0,0), (2,2),
                 bcs=["b","r","t","l"],
                 leftdomain=1)
geo.AddRectangle((1,1), (1.5,1.5),
                 bcs=["b2","r2","t2","l2"],
                 leftdomain=2, rightdomain=1)
geo.SetMaterial (1, "outer")
geo.SetMaterial (2, "inner")
initmesh = Mesh(geo.GenerateMesh(maxh=0.1))
Draw(initmesh)

D = initmesh.dim
t = CoordCF(D)
c = Vector(2)
c[0] = 1
c[1] = 4

cf = CoefficientFunction(list(c))
sq = sqrt(0.5);
bdd = CoefficientFunction((
    sin( cf*t+sq*(x+y) ),
    sq*cos(cf*t+sq*(x+y)),
    sq*cos(cf*t+sq*(x+y)),
    cf*cos(cf*t+sq*(x+y))
    ))
Draw(bdd, initmesh, 'piecewise')
# input()

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
Draw(gfu,initmesh,'sol')
# Draw(gfu,initmesh,'sol',autoscale=False,min=-0.1,max=0.1)
# Draw(gfu,initmesh,'sol',autoscale=False,min=-0.01,max=0.01)

wavefront = EvolveTentsMakeWavefront(order,initmesh,t_start,bdd)

for t in range(0,50):
    # with TaskManager():
    wavefront = EvolveTents(order,initmesh,c,t_step,wavefront,t_start,bdd)

    ipfct=IntegrationPointFunction(initmesh,intrule,wavefront)
    f = LinearForm(fes)
    f += SymbolicLFI(ipfct*v, intrule=intrule)
    f.Assemble()
    gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs()) * f.vec
    Redraw(blocking=True)

    t_start += t_step
    print("time: " + str(t_start))
# filename = "results/mov/sol"+str(t).zfill(3) +".jpg"
# Tcl_Eval("Ng_SnapShot .ndraw {};\n".format(filename))
