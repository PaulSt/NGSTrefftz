# import sys, os
# sys.path.append(os.path.join(os.path.dirname(sys.path[0])))
from ngstrefftz import *
# from ngstents import TentSlab
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve.TensorProductTools import *
from ngsolve import *
from dg import *
import time

# USE tenthight = wavespeed + 3

def SolveWaveTents(initmesh, order, c, t_step):
    """
    Solve using tent pitching
    >>> order = 4
    >>> SetNumThreads(4)
    >>> c = 1
    >>> t_step = 0.5
    >>> meshes=[ Mesh(SegMesh(4,0,math.pi)), Mesh(unit_square.GenerateMesh(maxh = 0.4)) , Mesh(unit_cube.GenerateMesh(maxh = 1))]
    >>> for initmesh in meshes:
    ...    for maxh in range(3):
    ...        SolveWaveTents(initmesh, order, c, t_step) # doctest:+ELLIPSIS
    ...        if initmesh.dim!=1:
    ...            initmesh.Refine()
    ...        else:
    ...            initmesh=Mesh(SegMesh(initmesh.ne*2,0,1))
    0.12...
    ...e-05
    ...e-06
    0.0...
    0.00...
    0.0001...
    0.1...
    0.0...
    0.003...

    same example with Neumann boundary conditions
    >>> meshes=[ Mesh(SegMesh(4,0,math.pi)), Mesh(unit_square.GenerateMesh(maxh = 0.4)) , Mesh(unit_cube.GenerateMesh(maxh = 1))]
    >>> for initmesh in meshes:
    ...    for i in range(0,len(initmesh.GetBoundaries())):
    ...        initmesh.ngmesh.SetBCName(i,"neumann")
    ...    for maxh in range(3):
    ...        SolveWaveTents(initmesh, order, c, t_step) # doctest:+ELLIPSIS
    ...        if initmesh.dim!=1:
    ...            initmesh.Refine()
    ...        else:
    ...            initmesh=Mesh(SegMesh(initmesh.ne*2,0,1))
    0.12...
    ...e-05
    ...e-06
    0.0...
    0.00...
    0.0001...
    0.1...
    0.0...
    0.003...
    """

    D = initmesh.dim
    t = CoordCF(D)
    t_start = 0

    if D==3:
        sq = sqrt(3.0);
        bdd = CoefficientFunction((
            sin(math.pi*x)*sin(math.pi*y)*sin(math.pi*z)*sin(math.pi*t*c*sq)/(sq*math.pi),
            cos(math.pi*x)*sin(math.pi*y)*sin(math.pi*z)*sin(math.pi*t*c*sq)/sq,
            sin(math.pi*x)*cos(math.pi*y)*sin(math.pi*z)*sin(math.pi*t*c*sq)/sq,
            sin(math.pi*x)*sin(math.pi*y)*cos(math.pi*z)*sin(math.pi*t*c*sq)/sq,
            sin(math.pi*x)*sin(math.pi*y)*sin(math.pi*z)*cos(math.pi*t*c*sq)*c
            ))
    elif D==2:
        sq = sqrt(2.0);
        bdd = CoefficientFunction((
            sin(math.pi*x)*sin(math.pi*y)*sin(math.pi*t*c*sq)/(sq*math.pi),
            cos(math.pi*x)*sin(math.pi*y)*sin(math.pi*t*c*sq)/sq,
            sin(math.pi*x)*cos(math.pi*y)*sin(math.pi*t*c*sq)/sq,
            sin(math.pi*x)*sin(math.pi*y)*cos(math.pi*t*c*sq)*c
            ))
    elif D==1:
        sq = sqrt(1.0);
        bdd = CoefficientFunction((
            sin(math.pi*x)*sin(math.pi*t*c*sq)/(sq*math.pi),
            cos(math.pi*x)*sin(math.pi*t*c*sq)/sq,
            sin(math.pi*x)*cos(math.pi*t*c*sq)*c
            ))

    local_ctau = True
    global_ctau = 2/3
    ts = TentSlab(initmesh, method="edge", heapsize=10*1000*1000)
    ts.SetMaxWavespeed(c)
    ts.PitchTents(dt=t_step, local_ct=local_ctau, global_ct=global_ctau)
    TT=TWave(order,ts,CoefficientFunction(c))
    TT.SetInitial(bdd)
    TT.SetBoundaryCF(bdd[D+1])
    if initmesh.ngmesh.GetBCName(0) == "neumann": TT.SetBoundaryCF(bdd[1:D+1])

    start = time.time()
    with TaskManager():
        TT.Propagate()
    timing = (time.time()-start)

    error = TT.Error(TT.GetWavefront(),TT.MakeWavefront(bdd,t_step))

    # return [error, timing, adiam]
    return error



def TestQTrefftz(order, initmesh, t_step,qtrefftz=1):
    """
    Solve using tent pitching and quasi-Trefftz basis functions
    >>> order = 4
    >>> SetNumThreads(4)
    >>> t_step = 1
    >>> for h in [4,8,16,32]:
    ...        initmesh = Mesh(SegMesh(h,0,math.pi))
    ...        TestQTrefftz(order,initmesh,t_step) # doctest:+ELLIPSIS
    0.02...
    0.001...
    ...e-05
    ...e-06
    >>> initmesh = Mesh(unit_square.GenerateMesh(maxh = 0.5))
    >>> for h in range(3):
    ...        TestQTrefftz(order,initmesh,t_step) # doctest:+ELLIPSIS
    ...        initmesh.Refine()
    0.004...
    0.001...
    ...e-05

    Compare to standard Trefftz basis
    >>> for h in [4,8,16,32]:
    ...        initmesh = Mesh(SegMesh(h,0,math.pi))
    ...        TestQTrefftz(order,initmesh,t_step,None) # doctest:+ELLIPSIS
    0.2...
    0.05...
    0.01...
    0.003...
    >>> initmesh = Mesh(unit_square.GenerateMesh(maxh = 0.5))
    >>> for h in range(3):
    ...        TestQTrefftz(order,initmesh,t_step,None) # doctest:+ELLIPSIS
    ...        initmesh.Refine()
    0.1...
    0.02...
    0.00...
    """

    # for i in range(0,len(initmesh.GetBoundaries())):
       # initmesh.ngmesh.SetBCName(i,"neumann")
    D = initmesh.dim
    t = CoordCF(D)
    t_start = 0

    c=1
    if D==1:
        ca=2.5
        bdd = CoefficientFunction((
                (x+1)**ca * exp(-sqrt(ca*(ca-1))*y),
                ca*(x+1)**(ca-1) * exp(-sqrt(ca*(ca-1))*y),
                -sqrt(ca*(ca-1)) * (x+1)**ca * exp(-sqrt(ca*(ca-1))*y)
            ))
        wavespeed=CoefficientFunction((x+1))

    else:
        ca=2.5
        bdd = CoefficientFunction((
                (x+y+1)**ca * exp(-sqrt(2*ca*(ca-1))*z),
                ca*(x+y+1)**(ca-1) * exp(-sqrt(2*ca*(ca-1))*z),
                ca*(x+y+1)**(ca-1) * exp(-sqrt(2*ca*(ca-1))*z),
                -sqrt(2*ca*(ca-1))*(x+y+1)**ca * exp(-sqrt(2*ca*(ca-1))*z)
            ))
        wavespeed=CoefficientFunction((x+y+1))


    local_ctau = True
    global_ctau = 2/3
    ts = TentSlab(initmesh, method="edge", heapsize=10*1000*1000)
    ts.SetMaxWavespeed(wavespeed)
    ts.PitchTents(dt=t_step, local_ct=local_ctau, global_ct=global_ctau)
    TT=TWave(order,ts,wavespeed,qtrefftz)
    TT.SetInitial(bdd)
    TT.SetBoundaryCF(bdd[D+1])

    start = time.time()
    with TaskManager():
        TT.Propagate()
    timing = (time.time()-start)

    error = TT.Error(TT.GetWavefront(),TT.MakeWavefront(bdd,t_step))
    # adiam = TT.MaxAdiam()

    # for t in Timers():
        # print(t)
    # return [error, timing, adiam]
    return error


def SolveWaveTentsFO(initmesh, order, c, t_step):
    """
    Solve using tent pitching first order system
    >>> order = 4
    >>> SetNumThreads(4)
    >>> c = 1
    >>> t_step = 2/sqrt(3)
    >>> initmesh = Mesh(unit_square.GenerateMesh(maxh = 0.4))
    >>> for maxh in range(3):
    ...     SolveWaveTentsFO(initmesh, order, c, t_step) # doctest:+ELLIPSIS
    ...     initmesh.Refine()
    0.0...
    0.00...
    0.0001...
    """

    D = initmesh.dim
    t = CoordCF(D)
    t_start = 0

    sq = sqrt(2.0);
    bdd = CoefficientFunction((
        cos(math.pi*x)*sin(math.pi*y)*sin(math.pi*t*c*sq)/sq,
        sin(math.pi*x)*cos(math.pi*y)*sin(math.pi*t*c*sq)/sq,
        sin(math.pi*x)*sin(math.pi*y)*cos(math.pi*t*c*sq)*c
        ))

    local_ctau = True
    global_ctau = 2/3
    ts = TentSlab(initmesh, method="edge", heapsize=10*1000*1000)
    ts.SetMaxWavespeed(c)
    ts.PitchTents(dt=t_step, local_ct=local_ctau, global_ct=global_ctau)
    TT=TWave(order,ts,CoefficientFunction(c))
    TT.SetInitial(bdd)
    TT.SetBoundaryCF(bdd[D])

    start = time.time()
    with TaskManager():
        TT.Propagate()
    timing = (time.time()-start)

    error = TT.Error(TT.GetWavefront(),TT.MakeWavefront(bdd,t_step))
    # V = L2(initmesh, order=order, dim=initmesh.dim+1)
    # u = GridFunction(V,"u")
    # TT.GetWave()

    return error


if __name__ == "__main__":
    # order = 4
    # SetNumThreads(1)
    # c = 1
    # t_step = 2/sqrt(3)
    # initmesh=Mesh(unit_cube.GenerateMesh(maxh = 0.25))
    # start = time.time()
    # print("Error",SolveWaveTents(initmesh, order, c, t_step))
    # print("PYTIME:", time.time()-start)
    # for t in Timers():
        # if 'tent' in t['name']:
            # print(t)
    # input()
    # # Error 0.0029433017038692647
    # # PYTIME: 42.57842946052551
    # # {'name': 'pitch tents', 'time': 0.011452735514933324, 'counts': 1, 'flops': 0.0, 'Gflop/s': 0.0}
    # # {'name': 'tent top bilinearform', 'time': 9.498263475069015, 'counts': 47284, 'flops': 73233459200.0, 'Gflop/s': 7.7101945415836015}
    # # {'name': 'tent top AAt', 'time': 4.70538792283988, 'counts': 47284, 'flops': 0.0, 'Gflop/s': 0.0}
    # # {'name': 'tent top calcshape', 'time': 5.122740580736957, 'counts': 47284, 'flops': 0.0, 'Gflop/s': 0.0}
    # # on new dell
    # # Error 0.0032559029714194806
    # # PYTIME: 22.690917253494263
    # # {'name': 'tent top bilinearform', 'time': 3.264730902699223, 'counts': 43518, 'flops': 67400678400.0, 'Gflop/s': 20.64509462151208}
    # # {'name': 'tent top AAt', 'time': 3.556020358299543, 'counts': 43518, 'flops': 0.0, 'Gflop/s': 0.0}
    # # {'name': 'tent top calcshape', 'time': 3.095868819707971, 'counts': 43518, 'flops': 0.0, 'Gflop/s': 0.0}

    # order = 7
    # SetNumThreads(1)
    # t_step = 1
    # h=32
    # mesh = CartSquare(h,h,xshift=3)
    # start = time.time()
    # print("error ",TestBessel(order,mesh,t_step)) # doctest:+ELLIPSIS
    # # e-13
    # print("time to solve",time.time()-start)
    # # 3.04s
    # for t in Timers():
        # if 'QTrefftz' in t['name']:
            # print(t)
    # input()

    # order = 4
    # initmesh = Mesh(unit_square.GenerateMesh(maxh = 0.1))
    # for i in range(3):
        # SetNumThreads(2**i)
        # start = time.time()
        # print(TestQTrefftz(order,initmesh,1)) # doctest:+ELLIPSIS
        # print("time to solve on ",2**i," threads ",time.time()-start)
        # input()

    import doctest
    doctest.testmod()
