from trefftzngs import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve.TensorProductTools import *
from ngsolve import *
import time


def TestSolution2D():
    truesol = x**2+2*y
    truesol = exp(-y)*cos(x)
    # truesol = exp(y)*exp(x)
    v0 = truesol
    sig0 = 0
    gD = 2*x
    gD = -exp(-y)*sin(x)
    U0 = truesol
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
    >>> N = 5
    >>> c = 1
    >>> order = 8
    >>> mesh = CartSquare(N,c*N)

    using Trefftz basis
    >>> fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True, basistype=0,heat=1)
    >>> Cartsolve2D(fes,c) # doctest:+ELLIPSIS
    [17.0, ..., ...e-09, ...e-08]

    or normal L2 basis, requiring the full system
    >>> fes = monomialfespace(mesh, order=order, dgjumps=True)
    >>> Cartsolve2D(fes,c,True) # doctest:+ELLIPSIS
    [81.0, ..., ...e-12, ...e-10]
    """
    if inputsol is None:
        inputsol = TestSolution2D()
    [truesol,U0,sig0,v0,gD] = inputsol

    start = time.time()
    [a,f] = DGheateqsys(fes,U0,v0,sig0,c,gD,fullsys,False,0.5,0.5,1)
    # print("DGsys: ", str(time.clock()-start))

    start = time.time()
    [gfu,cond] = DGsolve(fes,a,f)
    # print("DGsolve: ", str(time.clock()-start))

    [L2error, sH1error] = PostProcess(fes,truesol,gfu)

    dof=fes.ndof/fes.mesh.ne

    return [dof,cond,L2error,sH1error]


if __name__ == "__main__":
    # import doctest
    # doctest.testmod()

    N = 5
    c = 1
    order = 3
    mesh = CartSquare(N,c*N)
    [truesol,U0,sig0,v0,gD] = TestSolution2D()

    fes = monomialfespace(mesh, order=order, dgjumps=True)
    U = GridFunction(fes)
    U.Set(truesol)
    print("CHECK",sqrt(Integrate((truesol-U)*(truesol-U),mesh)))

    [a,f] = DGheateqsysnew(fes,fes,U0,v0,sig0,c,gD,False,False,0.5,0.5,1)
    [gfu,cond] = DGsolve(fes,a,f)
    [L2error, sH1error] = PostProcess(fes,truesol,gfu)
    print(L2error, sH1error, fes.ndof/fes.mesh.ne)

    fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True, basistype=0,heat=1)
    fes2 = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True, basistype=0,heat=1,heattest=1)
    U = GridFunction(fes)
    U.Set(truesol)
    print("CHECK",sqrt(Integrate((truesol-U)*(truesol-U),mesh)))

    [a,f] = DGheateqsysnew(fes,fes2,U0,v0,sig0,c,gD,False,False,0.5,0.5,1)
    [gfu,cond] = DGsolve(fes,a,f)
    [L2error, sH1error] = PostProcess(fes,truesol,gfu)
    print(L2error, sH1error, fes.ndof/fes.mesh.ne)
