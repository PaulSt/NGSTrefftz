from trefftzngs import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve.TensorProductTools import *
from ngsolve import *
import time
import pandas as pd


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
    gradtruesol = CoefficientFunction((-exp(-y)*sin(x),-exp(-y)*cos(x)))
    L2error = sqrt(Integrate((truesol - sol)*(truesol - sol), mesh))
    sH1error = sqrt(Integrate((gradtruesol - grad(sol))*(gradtruesol - grad(sol)), mesh))
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


    df = pd.DataFrame()
    for order in [2,3,4]:
        for N in range(5):
            mesh = CartSquare(2**N,2**N)
            fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True, basistype=0,heat=1)
            fes2 = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True, basistype=0,heat=1,heattest=1)
            U = GridFunction(fes)
            U.Set(truesol)
            print("CHECK",sqrt(Integrate((truesol-U)*(truesol-U),mesh)))

            [a,f] = DGheateqsysnew(fes,fes2,U0,v0,sig0,c,gD,False,False,0.5,0.5,1)
            [gfu,cond] = DGsolve(fes,a,f)
            [L2error, sH1error] = PostProcess(fes,truesol,gfu)
            rate = 0 if N==0 else log(errorold/L2error)/log(2)
            errorold=L2error
            print(L2error, sH1error, fes.ndof/fes.mesh.ne, rate)
            df = df.append({'hnr':N,'h':1/2**N,'p':order,'ndof':fes.ndof,'error':L2error}, ignore_index=True)
            df['p'] = df['p'].astype(int)
            df['hnr'] = df['hnr'].astype(int)

    df.to_csv('./heat.csv',index=False)
