from trefftzngs import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve.TensorProductTools import *
from ngsolve import *
import time

# USE tenthight = wavespeed + 3

def SolveWaveTents(initmesh, order, c, t_step):
    """
    Solve using tent pitching
    >>> order = 4
    >>> SetNumThreads(2)
    >>> c = 1
    >>> t_step = 2/sqrt(3)
    >>> meshes=[ Mesh(SegMesh(4,0,math.pi)), Mesh(unit_square.GenerateMesh(maxh = 0.4)) , Mesh(unit_cube.GenerateMesh(maxh = 1))]
    >>> for initmesh in meshes:
    ...    for maxh in range(3):
    ...        SolveWaveTents(initmesh, order, c, t_step) # doctest:+ELLIPSIS
    ...        if initmesh.dim!=1:
    ...            initmesh.Refine()
    ...        else:
    ...            initmesh=Mesh(SegMesh(initmesh.ne*2,0,1))
    0.164...
    ...e-05
    ...e-06
    0.01...
    0.001...
    0.0001...
    0.160...
    0.066...
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
    0.16...
    ...e-05
    ...e-06
    0.01...
    0.001...
    0.0001...
    0.1...
    0.057...
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

    TT=WaveTents(order,initmesh,CoefficientFunction(c),bdd)
    TT.SetWavefront(bdd,t_start)

    start = time.time()
    with TaskManager():
        TT.EvolveTents(t_step)
    timing = (time.time()-start)

    error = TT.Error(TT.GetWavefront(),TT.MakeWavefront(bdd,t_step))
    adiam = TT.MaxAdiam(t_step)

    # return [error, timing, adiam]
    return error


def TestAiry(order, mesh, t_step,qtrefftz=1):
    """
    Solve with quasi-Trefftz basis functions
    >>> order = 4
    >>> SetNumThreads(1)
    >>> t_step = 1
    >>> for h in [4,8,16,32]:
    ...        mesh = CartSquare(h,h)
    ...        TestAiry(order,mesh,t_step) # doctest:+ELLIPSIS
    5...e-05
    ...e-06
    ...e-07
    ...e-08
    """
    c=1
    bdd = CoefficientFunction((
        airy(-x-c)*cos(y),
        -airyp(-x-c)*cos(y),
        -airy(-x-c)*sin(y)
        ))
    wavespeed=CoefficientFunction(1/sqrt(c+x))

    U0=bdd[0]
    gD=bdd[2]
    v0=bdd[2]
    sig0=-bdd[1]

    fes = trefftzfespace(mesh, order=order, dgjumps=True, basistype=0, useshift=True, useqt=True)
    fes.SetWavespeed(wavespeed)
    [a,f] = DGwaveeqsys(fes,U0,v0,sig0,wavespeed,gD,True,False,alpha=0.5,beta=0.5,gamma=1,mu=0.5)
    gfu = DGsolve(fes,a,f)
    gradsol = CoefficientFunction((
        -airyp(-x-c)*cos(y),
        -airy(-x-c)*sin(y)
        ))
    dgerror = DGnormerror(fes,gfu,gradsol,wavespeed,alpha=0.5,beta=0.5)

    return dgerror


def TestAiryTent(order, initmesh, t_step,qtrefftz=1):
    """
    Solve using tent pitching and quasi-Trefftz basis functions
    >>> order = 4
    >>> SetNumThreads(2)
    >>> t_step = 1
    >>> for h in [4,8,16,32]:
    ...        initmesh = Mesh(SegMesh(h,0,math.pi))
    ...        TestAiryTent(order,initmesh,t_step) # doctest:+ELLIPSIS
    0.013...
    0.0007...
    ...e-05
    ...e-06
    >>> initmesh = Mesh(unit_square.GenerateMesh(maxh = 0.5))
    >>> for h in range(4):
    ...        TestAiryTent(order,initmesh,t_step) # doctest:+ELLIPSIS
    ...        initmesh.Refine()
    0.004...
    0.0005...
    ...e-05
    ...e-06

    Compare to standard Trefftz basis
    >>> initmesh = Mesh(unit_square.GenerateMesh(maxh = 0.5))
    >>> for h in range(4):
    ...        TestAiryTent(order,initmesh,t_step,None) # doctest:+ELLIPSIS
    ...        initmesh.Refine()
    0.013...
    0.001...
    0.0005...
    0.0001...
    """

    # for i in range(0,len(initmesh.GetBoundaries())):
       # initmesh.ngmesh.SetBCName(i,"neumann")
    D = initmesh.dim
    t = CoordCF(D)
    t_start = 0

    c=1
    if D is 1:
        bdd = CoefficientFunction((
            airy(-x-c)*cos(y),
            -airyp(-x-c)*cos(y),
            -airy(-x-c)*sin(y)
            ))
        wavespeed=CoefficientFunction(1/sqrt(c+x))
    else:
        bdd = CoefficientFunction((
            airy(-x-y-c)*cos(sqrt(2)*z),
            -airyp(-x-y-c)*cos(sqrt(2)*z),
            -airyp(-x-y-c)*cos(sqrt(2)*z),
            -airy(-x-y-c)*sin(sqrt(2)*z)*sqrt(2)
            ))
        wavespeed=CoefficientFunction(1/sqrt(c+x+y))


    TT=WaveTents(order,initmesh,wavespeed,bdd,qtrefftz)
    TT.SetWavefront(bdd,t_start)

    start = time.time()
    with TaskManager():
        TT.EvolveTents(t_step)
    timing = (time.time()-start)

    error = TT.Error(TT.GetWavefront(),TT.MakeWavefront(bdd,t_step))
    adiam = TT.MaxAdiam(t_step)

    # for t in Timers():
        # print(t)
    # return [error, timing, adiam]
    return error


def TestBessel(order, initmesh, t_step):
    """
    Solve using tent pitching and quasi-Trefftz basis functions
    >>> order = 4
    >>> SetNumThreads(2)
    >>> t_step = 1
    >>> for h in [4,8,16,32]:
    ...        initmesh = Mesh(SegMesh(h,3,4)) # need x>2
    ...        TestBessel(order,initmesh,t_step) # doctest:+ELLIPSIS

    0.013...
    0.0007...
    ...e-05
    ...e-06
    """

    D = initmesh.dim
    t = CoordCF(D)
    t_start = 0

    c=1
    if D is 1:
        bdd = CoefficientFunction((
            (sin(x)-x*cos(x))/x**2 * cos(t),
            (-2*x**(-3)*sin(x)+2*x**(-2)*cos(x)+x**(-1)*sin(x)) * cos(t),
            -(sin(x)-x*cos(x))/x**2 * sin(t)
            ))
        wavespeed=x**2-2
        BB = x**2

    TT=WaveTents(order,initmesh,wavespeed,bdd,)
    TT.SetWavefront(bdd,t_start)

    start = time.time()
    with TaskManager():
        TT.EvolveTents(t_step)
    timing = (time.time()-start)

    error = TT.Error(TT.GetWavefront(),TT.MakeWavefront(bdd,t_step))
    adiam = TT.MaxAdiam(t_step)


def TestSolution2D(fes,c,timeoffset=0):
    k = 3
    truesol = sin( k*(c*y + x) )
    v0 = c*k*cos(k*(c*y+x))
    sig0 = -k*cos(k*(c*y+x))
    gD = v0
    U0 = GridFunction(fes)
    U0.Set(truesol)
    return [truesol,U0,sig0,v0,gD]


def PostProcess(fes, truesol, sol):
    mesh = fes.mesh
    U = GridFunction(fes)
    U.Set(truesol)
    L2error = sqrt(Integrate((truesol - sol)*(truesol - sol), mesh))
    sH1error = sqrt(Integrate((grad(U) - grad(sol))*(grad(U) - grad(sol)), mesh))
    return [L2error,sH1error]


def Cartsolve2D(fes,c,fullsys=False,inputsol=None):
    """
    We can solve on a simple rectangle grid
    >>> N = 4
    >>> c = 2
    >>> order = 8
    >>> mesh = CartSquare(N,c*N)

    using Trefftz basis
    >>> fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True, basistype=0)
    >>> Cartsolve2D(fes,c) # doctest:+ELLIPSIS
    [17.0, ..., ...e-09, ...e-08]

    or normal L2 basis, requiring the full system
    >>> fes = L2(mesh, order=order, dgjumps=True)
    >>> Cartsolve2D(fes,c,True) # doctest:+ELLIPSIS
    [81.0, ..., ...e-12, ...e-10]
    """
    if inputsol is None:
        inputsol = TestSolution2D(fes,c)
    [truesol,U0,sig0,v0,gD] = inputsol

    start = time.time()
    [a,f] = DGwaveeqsys(fes,U0,v0,sig0,c,gD,fullsys,False,0.5,0.5,1)
    # print("DGsys: ", str(time.clock()-start))

    start = time.time()
    gfu = DGsolve(fes,a,f)
    cond = 0
    # print("DGsolve: ", str(time.clock()-start))

    [L2error, sH1error] = PostProcess(fes,truesol,gfu)

    dof=fes.ndof/fes.mesh.ne

    return [dof,cond,L2error,sH1error]


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
    # Error 0.0029433017038692647
    # PYTIME: 42.57842946052551
    # {'name': 'pitch tents', 'time': 0.011452735514933324, 'counts': 1, 'flops': 0.0, 'Gflop/s': 0.0}
    # {'name': 'tent top bilinearform', 'time': 9.498263475069015, 'counts': 47284, 'flops': 73233459200.0, 'Gflop/s': 7.7101945415836015}
    # {'name': 'tent top AAt', 'time': 4.70538792283988, 'counts': 47284, 'flops': 0.0, 'Gflop/s': 0.0}
    # {'name': 'tent top calcshape', 'time': 5.122740580736957, 'counts': 47284, 'flops': 0.0, 'Gflop/s': 0.0}
    # on new dell
    # Error 0.0032559029714194806
    # PYTIME: 22.690917253494263
    # {'name': 'tent top bilinearform', 'time': 3.264730902699223, 'counts': 43518, 'flops': 67400678400.0, 'Gflop/s': 20.64509462151208}
    # {'name': 'tent top AAt', 'time': 3.556020358299543, 'counts': 43518, 'flops': 0.0, 'Gflop/s': 0.0}
    # {'name': 'tent top calcshape', 'time': 3.095868819707971, 'counts': 43518, 'flops': 0.0, 'Gflop/s': 0.0}

    import doctest
    doctest.testmod()
