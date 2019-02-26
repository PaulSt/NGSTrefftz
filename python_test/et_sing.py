# import netgen.gui
from ngsolve import *
from ngsolve.solve import Tcl_Eval # for snapshots
from trefftzngs import *
from prodmesh import *
from testcases import *
import time
import scipy.special

from ngsolve import *
from netgen.meshing import MeshingParameters
from netgen.geom2d import unit_square

maxH = [0.08,0.07,0.06,0.05,0.04,  0.08,0.07,0.06,0.05,0.04]
minH = [0.08,0.07,0.06,0.05,0.04,  0.01,0.00875,0.0075,0.00625,0.005]

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
            mp.RestrictH (x=xk, y=yk, z=0, h=max(minh, maxh*sqrt(sqrt(sqrt(r)))) )

    initmesh = Mesh( oLshapeMesh(maxh,mp) )
    # initmesh = Mesh( unit_square.GenerateMesh(mp=mp))

    # for i in range(0,len(initmesh.GetBoundaries())):
       # initmesh.ngmesh.SetBCName(i,"neumann")

    order = 5
    c = 1
    t_start = 0
    t_step = 0.5

    D = 2
    eltyp = ET.TRIG
    intrule = IntegrationRule(eltyp,2*order)
    irsize = len(intrule.points)

    bdd = singular(D,c)
    t = CoordCF(D)
    r2 = ((x)*(x)+(y)*(y))
    r = sqrt((x)*(x)+(y)*(y))
    at2 = atan2(-y,-x)+math.pi #2*atan((y)/(r+(x)))
    at2x=(-y/r2)
    at2y=(x/r2)
    alp = 10 #2.4048
    nrb = 2.0/3.0
    # Draw(bdd ,initmesh,'u')
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
    # Draw(gfu,initmesh,'sol')

    runt = time.time()
    with TaskManager():
        for t in range(0,1):
            wavefront = EvolveTents(order,initmesh,c,t_step,wavefront,t_start, bdd )
            # wavefront = EvolveTentsMakeWavefront(order,initmesh,t_start + t_step,bdd)

            ipfct=IntegrationPointFunction(initmesh,intrule,wavefront)
            f = LinearForm(fes)
            f += SymbolicLFI(ipfct*v, intrule=intrule)
            f.Assemble()
            gfu.vec.data = a.mat.Inverse() * f.vec
            # Redraw(blocking=True)
            # input()

            t_start += t_step
            # print("time: " + str(t_start))
            # filename = "results/mov/sol"+str(t).zfill(3)+".jpg"
            # Tcl_Eval("Ng_SnapShot .ndraw {};\n".format(filename))
    runtime.append(time.time()-runt)
    error.append(Integrate((gfu - cos(alp*t_start)*sin(nrb*at2)*bessel(alp*r))*(gfu - cos(alp*t_start)*sin(nrb*at2)*bessel(alp*r)), initmesh))
    dof.append(EvolveTentsDofs(order,initmesh,c,t_step)) #initmesh.ne*(scipy.special.binom(3-1 + order, order) + scipy.special.binom(3-1 + order-1, order-1)))
    # error.append(EvolveTentsError(order,initmesh,wavefront,EvolveTentsMakeWavefront(order,initmesh,t_start,bdd)))

for i in range(len(error)):
    print("maxh", maxH[i])
    print("minh", minH[i])
    print("runtime", runtime[i])
    print("error", error[i])
    print("dof",dof[i])
    if i>0:
        print("rate",log(error[i-1]/error[i])/log(maxH[i-1]/maxH[i]))
        print("dofrate",log(error[i-1]/error[i])/log(sqrt(dof[i]/dof[i-1])))
    print()

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mpltools import annotation
fig = plt.figure()
ax = fig.gca()
for i in range(2):
    plt.loglog( [1/sqrt(dof[j]) for j in range(i*5,(i+1)*5)], error[i*5:(i+1)*5], '-*', label="ord="+str(1))
    annotation.slope_marker((1/sqrt(dof[i*5],error[i*5]), order+0.5-1,invert=True, ax=ax,size_frac=0.1)
#plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')
# for label in ax.xaxis.get_ticklabels()[::2]:
    # label.set_visible(False)
    #plt.title('Wavespeed: ' + str(c) + " dim: " + str(initmesh.dim) + "+1")
# plt.title(title)
plt.gca().invert_xaxis()
plt.legend()
plt.ylabel("error")
plt.xlabel("maxh")
plt.show()
#plt.savefig("results/pvtv3_"+label[yaxis]+label[xaxis]+".png")
