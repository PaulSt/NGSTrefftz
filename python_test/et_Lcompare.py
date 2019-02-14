from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from trefftzngs import *
# import netgen.gui
from ngsolve import *
from prodmesh import *
from ngsolve.solve import Tcl_Eval # for snapshots
from testcases import *
import time

from ngsolve import *
from netgen.geom2d import unit_square
from netgen.meshing import MeshingParameters

order = 2
c = 1
t_step = 0.1


maxh=[0.04,0.03,0.02,0.01,0.005001,0.005]
minh=[0,0,0,0,0,0]
# maxh=[0.01001,0.01]
# minh=[0,0]
eltyp = ET.TRIG
intrule = IntegrationRule(eltyp,2*order)
irsize = len(intrule.points)

wave=[]
runtime = []
countrun = 0

for m,k in zip(maxh, minh):
    countrun += 1
    t_start = 0
    mp = MeshingParameters (maxh = m)
    refpoints = 100
    for i in range(0, refpoints+1):
        for j in range(0, refpoints+1):
            xk = i/refpoints
            yk = j/refpoints
            r = sqrt(((xk-0.5)*(xk-0.5)+(yk-0.5)*(yk-0.5)))
            mp.RestrictH (x=xk, y=yk, z=0, h=max(k, m*sqrt(sqrt(r))) )

    if k==0:
        initmesh = Mesh( LshapeMesh(m,0) )
    else:
        initmesh = Mesh( LshapeMesh(maxh,mp) )
    # initmesh = Mesh( unit_square.GenerateMesh(maxh=m) )
    Draw(initmesh)
    for i in range(len(initmesh.GetBoundaries())):
       initmesh.ngmesh.SetBCName(i,"neumann")

    fes = L2(initmesh, order=order)
    u,v = fes.TnT()
    gfu = GridFunction(fes)
    a = BilinearForm(fes)
    a += SymbolicBFI(u*v)
    a.Assemble()
    inv = a.mat.Inverse()
    bdd = vertgausspw(2,c)
    # bdd = simplesin(2,c)
    wavefront = EvolveTentsMakeWavefront(order,initmesh,t_start,bdd )
    # Draw(gfu,initmesh,'sol')

    rt = time.time()
    with TaskManager():
        for t in range(4):
            wavefront = EvolveTents(order,initmesh,c,t_step,wavefront,t_start, bdd )

            ipfct=IntegrationPointFunction(initmesh,intrule,wavefront)
            f = LinearForm(fes)
            f += SymbolicLFI(ipfct*v, intrule=intrule)
            f.Assemble()
            gfu.vec.data = inv * f.vec
            # Redraw(blocking=True)

            t_start += t_step
            print("time: " + str(t_start))
    runtime.append(time.time()-rt)
    wave.append(gfu)

with open('lcompare'+'finehom'+".txt", "w") as text_file:
    print(runtime,file=text_file)
    print('maxh ',maxh,file=text_file)
    print('minh ',minh,file=text_file)
    for i in range(len(wave)-1):
        print('maxh ',maxh[i],' minh ',minh[i],file=text_file)
        errpre = sqrt(Integrate((wave[i-1]-wave[-1])*(wave[i-1]-wave[-1]),initmesh))
        err = sqrt(Integrate((wave[i]-wave[-1])*(wave[i]-wave[-1]),initmesh))
        print('error ', err,file=text_file)
        print('reler ', sqrt(Integrate((wave[i]-wave[-1])*(wave[i]-wave[-1]),initmesh))
            /sqrt(Integrate((wave[-1]*wave[-1]),initmesh)),file=text_file)
        print('rate ', log(errpre/err)/log(maxh[i-2]/maxh[i]),file=text_file)
    # Draw((wave[i]-wave[-1])*(wave[i]-wave[-1])/(wave[-1]*wave[-1]),initmesh,'err'+str(i))
