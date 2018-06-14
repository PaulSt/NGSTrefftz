from ngsolve import *
from trefftzngs import *

def GetFESTrefftz(mesh,c=1):
    return FESpace("trefftzfespace", mesh, order = 4, wavespeed = c, dgjumps=True, basistype=0)

def EvolveTent(fes,U0,gD=0,c=1):
    D = fes.mesh.dim - 1
    U = fes.TrialFunction()
    V = fes.TestFunction()
    gU = grad(U)
    gV = grad(V)

    v = gU[D]
    sig = CoefficientFunction(tuple([-gU[i] for i in  range(D)]))
    w = gV[D]
    tau = CoefficientFunction(tuple([-gV[i] for i in  range(D)]))

    vo = gU.Other()[D]
    sigo = CoefficientFunction(tuple([-gU.Other()[i] for i in  range(D)]))
    wo = gV.Other()[D]
    tauo = CoefficientFunction(tuple([-gV.Other()[i] for i in  range(D)]))

    h = specialcf.mesh_size
    n = specialcf.normal(D+1)
    n_t = n[D]/Norm(n)
    n_x = CoefficientFunction( tuple([n[i]/Norm(n) for i in  range(D)]) )

    mean_v = 0.5*(v+vo)
    mean_w = 0.5*(w+wo)
    mean_sig = 0.5*(sig+sigo)
    mean_tau = 0.5*(tau+tauo)

    jump_vx = ( v - vo ) * n_x
    jump_wx = ( w - wo ) * n_x
    jump_sigx = (( sig - sigo ) * n_x)
    jump_taux = (( tau - tauo ) * n_x)

    jump_vt = ( v - vo ) * n_t
    jump_wt = ( w - wo ) * n_t
    jump_sigt = ( sig - sigo ) * n_t
    jump_taut = ( tau - tauo ) * n_t

    jump_Ut = (U - U.Other()) * n_t
    #gfu.vec.data = a.mat.Inverse() * f.vec

    timelike = n_x*n_x # n_t=0
    spacelike = n_t**2 # n_x=0

    alpha = 0.5 #pow(10,5)
    beta = 0.5 #pow(10,5)
    gamma = 1

    a = BilinearForm(fes)
    a += SymbolicBFI( timelike * ( mean_v*jump_taux + mean_sig*jump_wx + alpha*jump_vx*jump_wx + beta*jump_sigx*jump_taux ), VOL, skeleton=True ) #time like faces
    a += SymbolicBFI( IfPos(n_t,1,0) * spacelike * ( pow(c,-2)*v*w + sig*tau ), element_boundary=True) #t=T (or *x)
    a += SymbolicBFI( IfPos(n_t,0,1) * spacelike * ( gamma * U*V ), element_boundary=True ) #BND correction term to recover sol of second order system
    a += SymbolicBFI( ( sig*n_x*w + alpha*v*w ), BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True) #dirichlet boundary 'timelike'
    a.Assemble()
    f = GridFunction(fes)
    rhs = BilinearForm(fes)
    rhs += SymbolicBFI( IfPos(n_t,0,1) * spacelike * ( pow(c,-2)*vo*w + sigo*tau ), element_boundary=True ) #space like faces, w/o x jump
    rhs += SymbolicBFI( IfPos(n_t,0,1) * spacelike * ( U.Other()*V ), element_boundary=True )
    rhs += SymbolicBFI( ( gD * (alpha*w - tau*n_x) ), BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True) #dirichlet boundary 'timelike'
    rhs.Apply(U0.vec,f.vec)

    gfu = GridFunction(fes, name="uDG")
    rows,cols,vals = a.mat.COO()
    A = sp.sparse.csr_matrix((vals,(rows,cols)))
    gfu.vec.FV().NumPy()[:] = sp.sparse.linalg.spsolve(A,f.vec.FV())
    return gfu









