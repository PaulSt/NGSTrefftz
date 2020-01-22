from trefftzngs import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve.TensorProductTools import *
from ngsolve import *
import time


def SolveWaveTents(initmesh, order, c, t_step):
    """
    Solve using tent pitching
    >>> order = 4
    >>> SetNumThreads(2)
    >>> c = 1
    >>> t_step = 2/sqrt(3)
    >>> ms = [0,1,2]
    >>> meshes=[ Mesh(SegMesh(4,0,math.pi)), Mesh(unit_square.GenerateMesh(maxh = 0.4)) , Mesh(unit_cube.GenerateMesh(maxh = 1))]

    >>> for initmesh in meshes:
    ...    for maxh in ms:
    ...        SolveWaveTents(initmesh, order, c, t_step) # doctest:+ELLIPSIS
    ...        if maxh != ms[-1] and initmesh.dim!=1:
    ...            initmesh.Refine()
    ...        elif initmesh.dim==1:
    ...            initmesh=Mesh(SegMesh(initmesh.ne*2,0,1))
    0.206...
    ...e-05
    ...e-06
    0.014...
    0.001...
    0.000...
    0.155...
    0.064...
    0.003...

    if i ever feel like checking the times, here is example full output, 2 cores on my laptop:
    [0.20661893329359113, 0.0006449222564697266, 1.5707963267948966]
    [5.466177093375104e-05, 0.0021283626556396484, 0.25]
    [3.5582731954669597e-06, 0.007055044174194336, 0.125]
    [0.014848896883671806, 0.021959304809570312, 1.0302183658236426]
    [0.0017425232736635377, 0.1491250991821289, 0.6291240085324618]
    [0.00014440246076325523, 1.1876091957092285, 0.38677266732559074]
    [0.15532255816903992, 0.2532687187194824, 1.8027756377319946]
    [0.06472569782925318, 3.0337188243865967, 1.469563021303764]
    [0.0037084177550868944, 47.547746419906616, 0.7499999999996284]
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


def TestAiry(order, h, t_step):
    """
    Solve using tent pitching
    >>> order = 4
    >>> SetNumThreads(2)
    >>> t_step = 2/sqrt(3)
    >>> for h in [4,8,16,32]:
    ...        TestAiry(order,h,t_step) # doctest:+ELLIPSIS
    0.0112...
    0.0027...
    0.0004...
    ...e-05
    """

    initmesh = Mesh(SegMesh(h,0,math.pi))
    D = initmesh.dim
    t = CoordCF(D)
    t_start = 0

    c=1
    bdd = CoefficientFunction((
            airy(-x-c)*cos(t),
            -airyp(-x-c)*cos(t),
            -airy(-x-c)*sin(t)
        ))

    TT=WaveTents(order,initmesh,CoefficientFunction(1/sqrt(c+x)),bdd)
    TT.SetWavefront(bdd,t_start)

    start = time.time()
    with TaskManager():
        TT.EvolveTents(t_step)
    timing = (time.time()-start)

    error = TT.Error(TT.GetWavefront(),TT.MakeWavefront(bdd,t_step))
    adiam = TT.MaxAdiam(t_step)

    # return [error, timing, adiam]
    return error


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
    [81.0, ..., ...e-10, ...e-09]
    """
    if inputsol is None:
        inputsol = TestSolution2D(fes,c)
    [truesol,U0,sig0,v0,gD] = inputsol

    start = time.time()
    [a,f] = DGeqsys(fes,U0,v0,sig0,c,gD,fullsys,False,0.5,0.5,1)
    # print("DGsys: ", str(time.clock()-start))

    start = time.time()
    [gfu,cond] = DGsolve(fes,a,f)
    # print("DGsolve: ", str(time.clock()-start))

    [L2error, sH1error] = PostProcess(fes,truesol,gfu)

    dof=fes.ndof/fes.mesh.ne

    return [dof,cond,L2error,sH1error]


if __name__ == "__main__":
    import doctest
    doctest.testmod()
