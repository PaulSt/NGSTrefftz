from trefftzngs import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve.TensorProductTools import *
from ngsolve import *
import time


def TestSolution(c=1,timeoffset=0):
    k = 3
    truesol = sin( k*(c*y + x) )
    sig0 = -k*cos(k*(c*y+x))
    v0 = c*k*cos(k*(c*y+x))
    gD = v0
    U0 = truesol
    return [truesol,U0,sig0,v0,gD]

def TestSolutionSmoothWavespeed(timeoffset=0):
    ca=2.5
    truesol = (x+1)**ca * exp(-sqrt(ca*(ca-1))*y)
    sig0 = -ca*(x+1)**(ca-1) * exp(-sqrt(ca*(ca-1))*y)
    v0 = -sqrt(ca*(ca-1))*(x+1)**ca * exp(-sqrt(ca*(ca-1))*y)
    gD = v0
    U0 = truesol
    return [truesol,U0,sig0,v0,gD]


def PostProcess(fes, truesol, sol):
    mesh = fes.mesh
    U = GridFunction(fes)
    U.Set(truesol)
    L2error = sqrt(Integrate((truesol - sol)*(truesol - sol), mesh))
    sH1error = sqrt(Integrate((grad(U) - grad(sol))*(grad(U) - grad(sol)), mesh))
    return [L2error,sH1error]


def Cartsolve2D(fes,c,inputsol,fullsys=False):
    """
    We can solve on a simple rectangle grid
    >>> N = 4
    >>> c = 2
    >>> order = 8
    >>> mesh = CartSquare(N,c*N)

    using Trefftz basis
    >>> fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True, basistype=0)
    >>> Cartsolve2D(fes,c,TestSolution(c)) # doctest:+ELLIPSIS
    [17.0, ..., ...e-09, ...e-08]

    or normal L2 basis, requiring the full system
    >>> fes = L2(mesh, order=order, dgjumps=True)
    >>> Cartsolve2D(fes,c,TestSolution(c),True) # doctest:+ELLIPSIS
    [81.0, ..., ...e-12, ...e-10]

    or quasi-Trefftz basis for smooth wavespeeds, also requiring the full system
    >>> c = CoefficientFunction((x+1))
    >>> fes = trefftzfespace(mesh, order=order, useqt=True)
    >>> fes.SetWavespeed(c)
    >>> Cartsolve2D(fes,c,TestSolutionSmoothWavespeed(),True) # doctest:+ELLIPSIS
    [17.0, ..., ...e-11, ...e-11]
    """
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
