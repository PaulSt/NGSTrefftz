from ngsolve import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
import time
import math
from trefftzngs import *
# SetHeapSize(1000*1000*1000)


def SolveTrefftzTents(mesh, order, finaltime):
    D = initmesh.dim
    t = CoordCF(D)

    sq = sqrt(3.0);
    bdd = CoefficientFunction((
        sin(math.pi*x)*sin(math.pi*y)*sin(math.pi*z)*sin(math.pi*t*c*sq)/(sq*math.pi),
        cos(math.pi*x)*sin(math.pi*y)*sin(math.pi*z)*sin(math.pi*t*c*sq)/sq,
        sin(math.pi*x)*cos(math.pi*y)*sin(math.pi*z)*sin(math.pi*t*c*sq)/sq,
        sin(math.pi*x)*sin(math.pi*y)*cos(math.pi*z)*sin(math.pi*t*c*sq)/sq,
        sin(math.pi*x)*sin(math.pi*y)*sin(math.pi*z)*cos(math.pi*t*c*sq)*c
        ))

    TT=TrefftzTent(order,initmesh,c,bdd)
    wavefront = TT.MakeWavefront(bdd,t_start)

    start = time.time()
    # with TaskManager():
    TT.EvolveTents(t_step,wavefront)

    timing = (time.time()-start)
    print("time ",time.time()-start)
    error = TT.Error(wavefront,TT.MakeWavefront(bdd,t_step))
    print("error ", error)
    # adiam = EvolveTentsAdiam(initmesh,1,t_step)
    adiam=1

    return [error, timing, adiam]


if __name__ == '__main__':
    order = 3
    SetNumThreads(10)

    c = 1
    t_start = 0
    t_step = 2/sqrt(3)


    h1error = []
    adiam = []
    timer = []
# ms = [0.4,0.3,0.2,0.1]
    ms = [0,1]

    initmesh = Mesh(unit_cube.GenerateMesh(maxh = 1))
    # initmesh = Mesh(unit_square.GenerateMesh(maxh = 1))
    for maxh in ms:
        print("RUN: ", maxh)
        # initmesh = Mesh(unit_cube.GenerateMesh(maxh = maxh))
        [error,timing,tentdiam] =  SolveTrefftzTents(initmesh, order, t_step)
        h1error.append(error)
        adiam.append(tentdiam)
        timer.append(timing)

        if maxh != ms[-1]:
            initmesh.Refine()

    print("rate")
    from math import *
    for i in range(len(ms)-1):
        print(log(h1error[i]/h1error[i+1])/log(2))
        # print(log(h1error[i]/h1error[i+1])/log(ms[i]/ms[i+1]))
    for i in range(len(ms)):
        print("adiam ",adiam[i]," maxh ", ms[i]," error ", h1error[i], " time ", timer[i])
    for i in range(len(ms)-1):
        print(log(h1error[i]/h1error[i+1])/log(adiam[i]/adiam[i+1]))
