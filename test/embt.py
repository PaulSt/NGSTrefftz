# import sys, os
# sys.path.append(os.path.join(os.path.dirname(sys.path[0])))
from ngsolve import *
from ngstrefftz import *
from dg import *
import time
SetNumThreads(3)

########################################################################
# PySVDTrefftz
########################################################################
try:
    import scipy as sp
    import numpy as np
except ImportError as e:
    pass

# import matplotlib.pylab as plt
# def trunc(values, decs=0):
    # return np.round(values*10**decs)/(10**decs)
# def spspy(sparsemat):
    # rows,cols,vals = sparsemat.COO()
    # A = sp.sparse.csr_matrix((vals,(rows,cols)))
    # # A = A.todense()
    # plt.spy(A); plt.plot(); plt.show()


def LocalP(A,eps=10**-9,printP=0):
    U,s,V = np.linalg.svd(A)
    nz = (s<eps).sum()
    nnz = (s>=eps).sum()
    S=np.diag(s)
    P = U[:,nnz:]
    # if printP==True:
        # print(trunc(A,decs=0))
        # print(trunc(U,decs=3))
        # print(trunc(P,decs=3))
    return P

def PySVDTrefftz(op,fes,eps):
    """
    # >>> fes = L2(mesh2d, order=5,  dgjumps=True)#,all_dofs_together=True)
    # >>> u,v = fes.TnT()
    # >>> uh = u.Operator("hesse")
    # >>> vh = v.Operator("hesse")
    # >>> op = (uh[0,0]+uh[1,1])*(vh[0,0]+vh[1,1])*dx
    # >>> P = PySVDTrefftz(op,fes,eps)
    # >>> a,f = dgell(fes,exactlap)
    # >>> rows,cols,vals = a.mat.COO()
    # >>> A = sp.sparse.csr_matrix((vals,(rows,cols)))
    # >>> A = A.todense()
    # >>> TA = P.transpose()@A@P

    # >>> tsol = np.linalg.solve(TA,P.transpose()@f.vec.FV())
    # >>> tpgfu = GridFunction(fes)
    # >>> tpgfu.vec.data = P@tsol
    # >>> sqrt(Integrate((tpgfu-exactlap)**2, mesh2d)) # doctest:+ELLIPSIS
    # 1...e-08
    """
    mesh = fes.mesh
    opbf = BilinearForm(fes)
    opbf += op
    opbf.Assemble()
    rows,cols,vals = opbf.mat.COO()
    A = sp.sparse.csr_matrix((vals,(rows,cols)))
    A = A.todense()
    localndof = int(fes.ndof/mesh.ne)
    trefftzndof = 2*order+1
    P = np.zeros([A.shape[1],(trefftzndof)*mesh.ne])
    for i in range(mesh.ne):
        # U,s,V = sp.linalg.svd(A[i*localndof:(i+1)*localndof,i*localndof:(i+1)*localndof])
        elmat = A[i*localndof:(i+1)*localndof,i*localndof:(i+1)*localndof]
        LP = LocalP(elmat,eps) #,i+1==mesh.ne)
        P[i*localndof:(i+1)*localndof,i*trefftzndof:(i+1)*trefftzndof] = LP
    # plt.spy(P); plt.plot(); plt.show()
    return P





########################################################################
# EmbTrefftz
########################################################################
def testembtrefftz(fes):
    """
    >>> fes = L2(mesh2d, order=5,  dgjumps=True)#,all_dofs_together=True)
    >>> testembtrefftz(fes) # doctest:+ELLIPSIS
    8...e-09
    """
    start = time.time()
    mesh = fes.mesh
    u,v = fes.TnT()
    uh = u.Operator("hesse")
    vh = v.Operator("hesse")
    op = (uh[0,0]+uh[1,1])*(vh[0,0]+vh[1,1])*dx
    # op = grad(u)*grad(v) * dx
    # op = InnerProduct(u.Operator("hesse"),v.Operator("hesse"))*dx
    with TaskManager():
        PP = TrefftzEmbedding(op,fes,eps)
    # spspy(PP)
    PPT = PP.CreateTranspose()
    a,f = dgell(fes,exactlap)
    TA = PPT@a.mat@PP
    TU = TA.Inverse()*(PPT*f.vec)
    tpgfu = GridFunction(fes)
    tpgfu.vec.data = PP*TU
    # print("trefftz time: ", time.time()-start)
    # for t in Timers():
        # if 'svdt' in t['name']:
            # print(t)
    return sqrt(Integrate((tpgfu-exactlap)**2, mesh))


def testembtrefftz_mixed(fes):
    """
    >>> fes = L2(mesh2d, order=5,  dgjumps=True)#,all_dofs_together=True)
    >>> testembtrefftz_mixed(fes) # doctest:+ELLIPSIS
    8...e-09
    """
    mesh = fes.mesh
    test_fes = L2(mesh, order=fes.globalorder-2,  dgjumps=True)#,all_dofs_together=True)
    u=fes.TrialFunction()
    v=test_fes.TestFunction()
    op = Lap(u)*(v)*dx
    startsvd = time.time()
    with TaskManager():
        PP = TrefftzEmbedding(op,fes,test_fes=test_fes)
    PPT = PP.CreateTranspose()
    a,f = dgell(fes,exactlap)
    TA = PPT@a.mat@PP
    TU = TA.Inverse()*(PPT*f.vec)
    tpgfu = GridFunction(fes)
    tpgfu.vec.data = PP*TU
    return sqrt(Integrate((tpgfu-exactlap)**2, mesh))


def testembtrefftznonsym(fes):
    """
    >>> fes = L2(mesh2d, order=5,  dgjumps=True)#,all_dofs_together=True)
    >>> testembtrefftznonsym(fes) # doctest:+ELLIPSIS
    8...e-09
    """
    start = time.time()
    mesh = fes.mesh
    u,v = fes.TnT()
    uh = u.Operator("hesse")
    vh = v.Operator("hesse")
    op = (uh[0,0]+uh[1,1])*v*dx
    with TaskManager():
        PP = TrefftzEmbedding(op,fes,eps)
    PPT = PP.CreateTranspose()
    a,f = dgell(fes,exactlap)
    TA = PPT@a.mat@PP
    TU = TA.Inverse()*(PPT*f.vec)
    tpgfu = GridFunction(fes)
    tpgfu.vec.data = PP*TU
    return sqrt(Integrate((tpgfu-exactlap)**2, mesh))

def testembtrefftzpoi(fes):
    """
    >>> fes = L2(mesh2d, order=5,  dgjumps=True)#,all_dofs_together=True)
    >>> testembtrefftzpoi(fes) # doctest:+ELLIPSIS
    3...e-09
    """
    mesh = fes.mesh
    start = time.time()
    u,v = fes.TnT()
    uh = u.Operator("hesse")
    vh = v.Operator("hesse")
    op = (uh[0,0]+uh[1,1])*(vh[0,0]+vh[1,1])*dx
    rhs = -exactpoi.Diff(x).Diff(x)-exactpoi.Diff(y).Diff(y)
    lop = -rhs*(vh[0,0]+vh[1,1])*dx
    with TaskManager():
        PP,ufv = TrefftzEmbedding(op,fes,lop,eps)
    PPT = PP.CreateTranspose()
    a,f = dgell(fes,exactpoi,rhs)
    TA = PPT@a.mat@PP
    TU = TA.Inverse()*(PPT*(f.vec-a.mat*ufv))
    tpgfu = GridFunction(fes)
    tpgfu.vec.data = PP*TU+ufv
    return sqrt(Integrate((tpgfu-exactpoi)**2, mesh))


def testembtrefftzpoi_mixed(fes):
    """
    >>> fes = L2(mesh2d, order=5,  dgjumps=True)#,all_dofs_together=True)
    >>> testembtrefftzpoi_mixed(fes) # doctest:+ELLIPSIS
    3...e-09
    """
    mesh = fes.mesh
    test_fes = L2(mesh, order=fes.globalorder-2,  dgjumps=True)#,all_dofs_together=True)
    u=fes.TrialFunction()
    v=test_fes.TestFunction()
    op = Lap(u)*(v)*dx
    rhs = -exactpoi.Diff(x).Diff(x)-exactpoi.Diff(y).Diff(y)
    lop = -rhs*v*dx
    with TaskManager():
        PP,ufv = TrefftzEmbedding(op,fes,lop,test_fes=test_fes)
    PPT = PP.CreateTranspose()
    a,f = dgell(fes,exactpoi,rhs)
    TA = PPT@a.mat@PP
    TU = TA.Inverse()*(PPT*(f.vec-a.mat*ufv))
    tpgfu = GridFunction(fes)
    tpgfu.vec.data = PP*TU+ufv
    #import netgen.gui
    #Draw(tpgfu)
    #input("")
    return sqrt(Integrate((tpgfu-exactpoi)**2, mesh))


def embtell(mesh,order):
    """
    >>> [embtell(mesh2d,5)] # doctest:+ELLIPSIS
    [...e-06]
    """

    if mesh.dim == 2:
        exact = sin(pi*(x+y))
        A = 1+x+y
        B = CF((x,-y))
        C = 3/(1+x+y)
        rhs = -sum( (A*CF((exact.Diff(x),exact.Diff(y))))[i].Diff(var) for i,var in enumerate([x,y])) + B*CF((exact.Diff(x),exact.Diff(y))) + C*exact
    elif mesh.dim == 3:
        exact = sin(pi*(x+y+z))
        A = 1+x+y+z
        B = CF((x,y,-2*z))
        C = 4/(1+x+y+z)
        rhs = -sum( (A*CF((exact.Diff(x),exact.Diff(y),exact.Diff(z))))[i].Diff(var) for i,var in enumerate([x,y,z])) + B*CF((exact.Diff(x),exact.Diff(y),exact.Diff(z))) + C*exact

    fes = L2(mesh, order=order, dgjumps=True)
    test_fes = L2(mesh, order=order-2, dgjumps=True)
    # etfes = EmbeddedTrefftzFES(fes)

    u = fes.TrialFunction()
    v = test_fes.TestFunction()
    if mesh.dim == 2:
        op = (-A*Lap(u)*v - CF((A.Diff(x),A.Diff(y)))*grad(u)*v + B*grad(u)*v + C*u*v)*dx
    elif mesh.dim == 3:
        op = -A*Lap(u)*v*dx - CF((A.Diff(x),A.Diff(y),A.Diff(z)))*grad(u)*v*dx + B*grad(u)*v*dx + C*u*v*dx
    lop = rhs*v*dx

    with TaskManager():
        PP,ufv = TrefftzEmbedding(op,fes,lop,test_fes=test_fes)
    PPT = PP.CreateTranspose()
    a,f = dgell(fes,Dbndc=exact,A=A,B=B,C=C,rhs=rhs)
    TA = PPT@a.mat@PP
    TU = TA.Inverse()*(PPT*(f.vec-a.mat*ufv))
    gfu = GridFunction(fes)
    gfu.vec.data = PP*TU+ufv

    error = sqrt(Integrate((gfu-exact)**2, mesh))
    return error

########################################################################
# Helmholtz
########################################################################

def testembtrefftzhelm(fes):
    """
    >>> fes = L2(mesh2d, order=5,  dgjumps=True, complex=True)
    >>> testembtrefftzhelm(fes) # doctest:+ELLIPSIS
    8...e-10
    """
    mesh = fes.mesh
    omega=1
    exact = exp(1j*sqrt(0.5)*(x+y))
    n = specialcf.normal(mesh.dim)
    bndc = CoefficientFunction((sqrt(0.5)*1j*exact, sqrt(0.5)*1j*exact))*n + 1j*omega*exact
    # test_fes = L2(mesh2d, order=order-2,  dgjumps=True, complex=True)
    u=fes.TrialFunction()
    v=fes.TestFunction()
    # v=test_fes.TestFunction()
    op = (-Lap(u))*(Lap(v))*dx#+1j*(-Lap(u)-u)*(Lap(v))*dx
    op += (-u)*(Lap(v))*dx#+1j*(-Lap(u)-u)*(Lap(v))*dx
    # op = (-Lap(u)-u)*((v))*dx+1j*(-Lap(u)-u)*((v))*dx
    with TaskManager():
        PP = TrefftzEmbedding(op,fes,eps=10**-8)
        # PP = TrefftzEmbedding(op,fes,test_fes=test_fes)
    PPT = PP.CreateTranspose()
    a,f = dghelm(fes,None,bndc,omega)
    TA = PPT@a.mat@PP
    TU = TA.Inverse()*(PPT*f.vec)
    tpgfu = GridFunction(fes)
    tpgfu.vec.data = PP*TU
    return sqrt(Integrate(InnerProduct(tpgfu-exact,tpgfu-exact).real, mesh))


########################################################################
# Stokes
########################################################################

def SolveStokes(mesh, k, nu, coeff, trefftz=True, ubnd=None, bndname="inflow"):
    """
    >>> nu = 1.0
    >>> zeta = cos(pi*x*(1-x)*y*(1-y))
    >>> pexact = sin(pi*(x+y))
    >>> uexact = CF((zeta.Diff(y), - zeta.Diff(x)))
    >>> graduexact = CF((uexact.Diff(x),uexact.Diff(y)),dims=(2,2)).trans
    >>> f1 = - nu*uexact[0].Diff(x).Diff(x) - nu*uexact[0].Diff(y).Diff(y) + pexact.Diff(x)
    >>> f2 = - nu*uexact[1].Diff(x).Diff(x) - nu*uexact[1].Diff(y).Diff(y) + pexact.Diff(y)
    >>> f = CF((f1,f2))

    >>> mesh = Mesh(unit_square.GenerateMesh(maxh=0.3) )
    >>> k = 5
    >>> uh, ph, ndof = SolveStokes(mesh, k, nu, f, trefftz=True) 
    >>> sqrt(Integrate(InnerProduct(uexact-uh,uexact-uh),mesh)) # doctest:+ELLIPSIS
    1...e-05
    >>> sqrt(Integrate(InnerProduct(pexact-ph,pexact-ph),mesh)) # doctest:+ELLIPSIS
    0.001...
    >>> ndof
    529
    """

    V = VectorL2(mesh, order=k, dgjumps=True)
    Q = L2(mesh, order=k - 1, dgjumps=True)
    Z = NumberSpace(mesh)
    fes = V * Q * Z
    u, v = fes.TrialFunction()[0], fes.TestFunction()[0]
    p, q = fes.TrialFunction()[1], fes.TestFunction()[1]
    lam, mu = fes.TrialFunction()[-1], fes.TestFunction()[-1]

    alpha = 20  # interior penalty param
    stab = 1e-9

    n = specialcf.normal(mesh.dim)
    h = specialcf.mesh_size

    jump_u = u - u.Other()
    jump_v = v - v.Other()
    mean_dudn = 0.5 * (grad(u) + grad(u.Other())) * n
    mean_dvdn = 0.5 * (grad(v) + grad(v.Other())) * n
    mean_q = 0.5 * n * (q + q.Other())
    mean_p = 0.5 * n * (p + p.Other())

    a = BilinearForm(fes)
    a += nu * InnerProduct(grad(u), grad(v)) * dx
    a += nu * alpha * k**2 / h * jump_u * jump_v * dx(skeleton=True)
    a += nu * (-mean_dudn * jump_v - mean_dvdn * jump_u) * dx(skeleton=True)
    a += nu * alpha * k**2 / h * u * v * ds(skeleton=True)
    a += nu * (-grad(u) * n * v - grad(v) * n * u) * ds(skeleton=True)
    a += (mean_p * jump_v + mean_q * jump_u) * dx(skeleton=True)
    a += (p * v * n + q * u * n) * ds(skeleton=True)
    if not trefftz:
        a += (-div(u) * q - div(v) * p) * dx
    a += (p * mu + q * lam) * dx
    a.Assemble()

    c = BilinearForm(fes)
    for i in a.integrators:
        c += i
    c += -stab * p * q * dx
    c += stab * lam * mu * dx
    c.Assemble()

    f = LinearForm(fes)
    f += coeff * v * dx(bonus_intorder=5)
    if ubnd:
        f += nu * alpha * k**2 / h * ubnd * v * ds(skeleton=True, definedon=mesh.Boundaries(bndname))
        f += nu * (- grad(v) * n * ubnd) * ds(skeleton=True, definedon=mesh.Boundaries(bndname))
    f.Assemble()

    gfu = GridFunction(fes)
    if trefftz:
        for i in range(V.ndof + Q.ndof, fes.ndof):
            fes.SetCouplingType(i, COUPLING_TYPE.HIDDEN_DOF)
        # ignoredofs = BitArray(fes.ndof)
        # ignoredofs[:] = False
        # for i in range(V.ndof + Q.ndof, fes.ndof):
            # ignoredofs[i] = True
        Vs = VectorL2(mesh, order=k - 2)
        Qs = L2(mesh, order=k - 1)
        test_fes = Vs * Qs
        wu, wp = test_fes.TestFunction()[0:2]

        op = nu*InnerProduct( grad(u),grad(wu) ) * dx \
            - nu*InnerProduct(grad(u)*n,wu) * dx(element_boundary=True) \
            + InnerProduct(grad(p),wu) * dx + div(u)*wp*dx

        timer = time.time()
        lop = coeff * wu * dx(bonus_intorder=10)
        PP, uf = TrefftzEmbedding(op,
                                  fes,
                                  lop,
                                  test_fes=test_fes)#, ignoredofs=ignoredofs)
        PPT = PP.CreateTranspose()
        # print(f"Trefftz embedding setup {time.time() - timer:5f} seconds")

        timer = time.time()
        TA = PPT @ a.mat @ PP
        TC = PPT @ c.mat @ PP
        Tgfu = CGSolver(TA, TC.Inverse()) * (PPT * (f.vec - a.mat * uf))
        gfu.vec.data = PP * Tgfu + uf
        # print(f"Trefftz embedding solve {time.time() - timer:5f} seconds")
        ndof = PP.shape[1]
    else:
        gfu.vec.data = CGSolver(a.mat, c.mat.Inverse()) * f.vec
        ndof = fes.ndof

    uh, ph = gfu.components[0:2]
    return uh, ph, ndof

########################################################################
# EmbTrefftzFESpace
########################################################################


def testembtrefftzfes(mesh,order):
    """
    >>> testembtrefftzfes(mesh2d,5) # doctest:+ELLIPSIS
    3...e-09
    """
    # exact = sin(10*x)*sin(10*y)
    # rhs = 200*sin(10*x)*sin(10*y)
    rhs = -exactpoi.Diff(x).Diff(x)-exactpoi.Diff(y).Diff(y)
    eps = 10**-7

    fes = L2(mesh, order=order, dgjumps=True)
    etfes = EmbeddedTrefftzFES(fes)

    u,v = fes.TnT()
    op = Lap(u)*Lap(v)*dx
    lop = -rhs*Lap(v)*dx

    uf = GridFunction(fes)
    uf.vec.data = etfes.SetOp(op,lf=lop,eps=eps)

    # etfes = Compress(etfes, etfes.FreeDofs())
    a,f = dgell(etfes,exactpoi,rhs,uf)
    # import pdb; pdb.set_trace()

    inv = a.mat.Inverse(inverse="sparsecholesky")
    gfu = GridFunction(etfes)
    gfu.vec.data = inv * f.vec
    nsol = gfu + uf

    return sqrt(Integrate((nsol-exactpoi)**2, mesh))



if __name__ == "__main__":
    import doctest
    doctest.testmod()
    #fes = L2(mesh2d, order=5,  dgjumps=True)#,all_dofs_together=True)
    #testembtrefftzpoi_mixed(fes) # doctest:+ELLIPSIS
