import netgen.gui
from ngsolve import *
from ngsolve.solve import Tcl_Eval # for snapshots
from trefftzngs import *
from prodmesh import *
from testcases import *
import time

from ngsolve import *
from netgen.meshing import MeshingParameters
from netgen.geom2d import unit_square

maxH = [0.1,0.08,0.06,0.04] #,0.08,0.06,0.04]
minH = [0.1,0.08,0.06,0.04] #,0.01,0.0075,0.005]

error = []
runtime = []
dof = []

for minh,maxh in zip(minH,maxH):
    mp = MeshingParameters (maxh = maxh)
    refpoints = 500
    for i in range(0, refpoints+1):
        for j in range(0, refpoints+1):
            xk = i/refpoints
            yk = j/refpoints
            r = sqrt(((xk)*(xk)+(yk)*(yk)))
            mp.RestrictH (x=xk, y=yk, z=0, h=max(minh, maxh*sqrt(sqrt(r))) )

    initmesh = Mesh( oLshapeMesh(maxh,mp) )
    # initmesh = Mesh( unit_square.GenerateMesh(mp=mp))

    # for i in range(0,len(initmesh.GetBoundaries())):
       # initmesh.ngmesh.SetBCName(i,"neumann")

    order = 3
    c = 1
    t_start = 0
    t_step = 0.5

    D = 2
    eltyp = ET.TRIG
    intrule = IntegrationRule(eltyp,2*order)
    irsize = len(intrule.points)

    bdd = singular(D,c)
    Draw(bdd,initmesh,'u')
    fes = H1(initmesh, order=order)
    u,v = fes.TnT()
    gfu = GridFunction(fes)
    a = BilinearForm(fes)
    a += SymbolicBFI(u*v)
    a.Assemble()
    wavefront = EvolveTentsMakeWavefront(order,initmesh,t_start,bdd )
    ipfct=IntegrationPointFunction(initmesh,intrule,wavefront)
    f = LinearForm(fes)
    f += SymbolicLFI(ipfct*v, intrule=intrule)
    f.Assemble()
    gfu.vec.data = a.mat.Inverse() * f.vec
    Draw(gfu,initmesh,'sol')
    input()

    runt = time.time()
    with TaskManager():
        # for t in range(0,1):
        wavefront = EvolveTents(order,initmesh,c,t_step,wavefront,t_start, bdd )
        # wavefront = EvolveTentsMakeWavefront(order,initmesh,t_start + t_step,bdd)

        ipfct=IntegrationPointFunction(initmesh,intrule,wavefront)
        f = LinearForm(fes)
        f += SymbolicLFI(ipfct*v, intrule=intrule)
        f.Assemble()
        gfu.vec.data = a.mat.Inverse() * f.vec
        # Redraw(blocking=True)

        t_start += t_step
        # print("time: " + str(t_start))
        # filename = "results/mov/sol"+str(t).zfill(3)+".jpg"
        # Tcl_Eval("Ng_SnapShot .ndraw {};\n".format(filename))
    runtime.append(time.time()-runt)
    error.append(EvolveTentsError(order,initmesh,wavefront,EvolveTentsMakeWavefront(order,initmesh,t_start,bdd)))

for i in range(len(error)):
    print("maxh", maxH[i])
    print("minh", minH[i])
    print("runtime", runtime[i])
    print("error", error[i])
    if i>0:
        print("rate",log(error[i-1]/error[i])/log(maxH[i-1]/maxH[i]))
    print()
