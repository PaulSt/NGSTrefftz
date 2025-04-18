# import sys, os
# sys.path.append(os.path.join(os.path.dirname(sys.path[0])))
from ngsolve import *
from ngstrefftz import *
from dg import *
import time
SetNumThreads(1)

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

def PySVDTrefftz(top,fes,eps):
    """
    # >>> fes = L2(mesh2d, order=5,  dgjumps=True)#,all_dofs_together=True)
    # >>> u,v = fes.TnT()
    # >>> uh = u.Operator("hesse")
    # >>> vh = v.Operator("hesse")
    # >>> top = (uh[0,0]+uh[1,1])*(vh[0,0]+vh[1,1])*dx
    # >>> P = PySVDTrefftz(top,fes,eps)
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
    opbf += top
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
    top = (uh[0,0]+uh[1,1])*(vh[0,0]+vh[1,1])*dx
    # top = grad(u)*grad(v) * dx
    # top = InnerProduct(u.Operator("hesse"),v.Operator("hesse"))*dx
    with TaskManager():
        emb = TrefftzEmbedding(top=top,eps=eps)
    PP = emb.GetEmbedding()
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
    fes_test = L2(mesh, order=fes.globalorder-2,  dgjumps=True)#,all_dofs_together=True)
    u=fes.TrialFunction()
    v=fes_test.TestFunction()
    top = Lap(u)*(v)*dx
    startsvd = time.time()
    with TaskManager():
        emb = TrefftzEmbedding(top=top)
    PP = emb.GetEmbedding()
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
    top = (uh[0,0]+uh[1,1])*v*dx
    with TaskManager():
        emb = TrefftzEmbedding(top=top,eps=eps)
    PP = emb.GetEmbedding()
    PPT = PP.CreateTranspose()
    a,f = dgell(fes,exactlap)
    TA = PPT@a.mat@PP
    TU = TA.Inverse()*(PPT*f.vec)
    tpgfu = GridFunction(fes)
    tpgfu.vec.data = PP*TU
    return sqrt(Integrate((tpgfu-exactlap)**2, mesh))

def testembtrefftzpoi(fes,test_gpsrhs = False):
    """
    >>> fes = L2(mesh2d, order=5,  dgjumps=True)#,all_dofs_together=True)
    >>> testembtrefftzpoi(fes) # doctest:+ELLIPSIS
    3...e-09
    >>> testembtrefftzpoi(fes,True) # doctest:+ELLIPSIS
    3...e-09
    """
    mesh = fes.mesh
    start = time.time()
    u,v = fes.TnT()
    uh = u.Operator("hesse")
    vh = v.Operator("hesse")
    top = (uh[0,0]+uh[1,1])*(vh[0,0]+vh[1,1])*dx
    rhs = -exactpoi.Diff(x).Diff(x)-exactpoi.Diff(y).Diff(y)
    trhs = -rhs*(vh[0,0]+vh[1,1])*dx
    with TaskManager():
        emb = TrefftzEmbedding(top=top,trhs=trhs,eps=eps)
    PP = emb.GetEmbedding()
    PPT = PP.CreateTranspose()
    if test_gpsrhs:
        uf = emb.GetParticularSolution(trhs)
    else:
        uf = emb.GetParticularSolution()
    a,f = dgell(fes,exactpoi,rhs)
    TA = PPT@a.mat@PP
    TU = TA.Inverse()*(PPT*(f.vec-a.mat*uf))
    tpgfu = GridFunction(fes)
    tpgfu.vec.data = PP*TU+uf
    return sqrt(Integrate((tpgfu-exactpoi)**2, mesh))


def testembtrefftzpoi_mixed(fes):
    """
    >>> fes = L2(mesh2d, order=5,  dgjumps=True)#,all_dofs_together=True)
    >>> testembtrefftzpoi_mixed(fes) # doctest:+ELLIPSIS
    3...e-09
    """
    mesh = fes.mesh
    fes_test = L2(mesh, order=fes.globalorder-2,  dgjumps=True)#,all_dofs_together=True)
    u=fes.TrialFunction()
    v=fes_test.TestFunction()
    top = Lap(u)*(v)*dx
    rhs = -exactpoi.Diff(x).Diff(x)-exactpoi.Diff(y).Diff(y)
    trhs = -rhs*v*dx
    with TaskManager():
        emb = TrefftzEmbedding(top=top,trhs=trhs)
    PP = emb.GetEmbedding()
    PPT = PP.CreateTranspose()
    a,f = dgell(fes,exactpoi,rhs)
    TA = PPT@a.mat@PP
    uf = emb.GetParticularSolution()
    TU = TA.Inverse()*(PPT*(f.vec-a.mat*uf))
    tpgfu = GridFunction(fes)
    tpgfu.vec.data = PP*TU+uf
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
    fes_test = L2(mesh, order=order-2, dgjumps=True)
    # etfes = EmbeddedTrefftzFES(fes)

    u = fes.TrialFunction()
    v = fes_test.TestFunction()
    if mesh.dim == 2:
        top = (-A*Lap(u)*v - CF((A.Diff(x),A.Diff(y)))*grad(u)*v + B*grad(u)*v + C*u*v)*dx
    elif mesh.dim == 3:
        top = -A*Lap(u)*v*dx - CF((A.Diff(x),A.Diff(y),A.Diff(z)))*grad(u)*v*dx + B*grad(u)*v*dx + C*u*v*dx
    trhs = rhs*v*dx

    with TaskManager():
        emb = TrefftzEmbedding(top=top,fes=fes,trhs=trhs,fes_test=fes_test)
    PP = emb.GetEmbedding()
    PPT = PP.CreateTranspose()
    uf = emb.GetParticularSolution()
    a,f = dgell(fes,Dbndc=exact,A=A,B=B,C=C,rhs=rhs)
    TA = PPT@a.mat@PP
    TU = TA.Inverse()*(PPT*(f.vec-a.mat*uf))
    gfu = GridFunction(fes)
    gfu.vec.data = PP*TU+uf

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
    # fes_test = L2(mesh2d, order=order-2,  dgjumps=True, complex=True)
    u=fes.TrialFunction()
    v=fes.TestFunction()
    # v=fes_test.TestFunction()
    top = (-Lap(u))*(Lap(v))*dx#+1j*(-Lap(u)-u)*(Lap(v))*dx
    top += (-u)*(Lap(v))*dx#+1j*(-Lap(u)-u)*(Lap(v))*dx
    # top = (-Lap(u)-u)*((v))*dx+1j*(-Lap(u)-u)*((v))*dx
    with TaskManager():
        emb = TrefftzEmbedding(top=top,fes=fes,eps=10**-8)
    PP = emb.GetEmbedding()
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

def SolveStokesTDG(mesh, k):
    """ 
    >>> SolveStokesTDG(mesh2d, 5) # doctest:+ELLIPSIS
    [1...e-05, 0.001..., 529]
    """
    nu = 1.0
    zeta = cos(pi*x*(1-x)*y*(1-y))
    pexact = sin(pi*(x+y))
    uexact = CF((zeta.Diff(y), - zeta.Diff(x)))
    graduexact = CF((uexact.Diff(x),uexact.Diff(y)),dims=(2,2)).trans
    f1 = - nu*uexact[0].Diff(x).Diff(x) - nu*uexact[0].Diff(y).Diff(y) + pexact.Diff(x)
    f2 = - nu*uexact[1].Diff(x).Diff(x) - nu*uexact[1].Diff(y).Diff(y) + pexact.Diff(y)
    rhs = CF((f1,f2))

    V = VectorL2(mesh, order=k, dgjumps=True)
    Q = L2(mesh, order=k - 1, dgjumps=True)
    Z = NumberSpace(mesh)
    fes = V * Q * Z
    u, p, z = fes.TrialFunction()
    n = specialcf.normal(mesh.dim)

    ah,ch,fh = StokesDG(fes, nu, rhs)
    a = BilinearForm(ah)
    a.Assemble()
    c = BilinearForm(ch)
    c.Assemble()
    f = LinearForm(fh)
    f.Assemble()

    gfu = GridFunction(fes)
    # for i in range(V.ndof + Q.ndof, fes.ndof):
        # fes.SetCouplingType(i, COUPLING_TYPE.HIDDEN_DOF)
    ignoredofs = BitArray(fes.ndof)
    ignoredofs[:] = False
    for i in range(V.ndof + Q.ndof, fes.ndof):
        ignoredofs[i] = True
    Vs = VectorL2(mesh, order=k - 2)
    Qs = L2(mesh, order=k - 1)
    fes_test = Vs * Qs
    wu, wp = fes_test.TestFunction()[0:2]

    top = nu*InnerProduct( grad(u),grad(wu) ) * dx \
        - nu*InnerProduct(grad(u)*n,wu) * dx(element_boundary=True) \
        + InnerProduct(grad(p),wu) * dx + div(u)*wp*dx

    timer = time.time()
    trhs = rhs * wu * dx(bonus_intorder=10)
    emb = TrefftzEmbedding(top=top, trhs=trhs, ignoredofs=ignoredofs)
    PP = emb.GetEmbedding()
    PPT = PP.CreateTranspose()
    uf = emb.GetParticularSolution()
    # print(f"Trefftz embedding setup {time.time() - timer:5f} seconds")

    TA = PPT @ a.mat @ PP
    TC = PPT @ c.mat @ PP
    Tgfu = CGSolver(TA, TC.Inverse()) * (PPT * (f.vec - a.mat * uf))
    gfu.vec.data = PP * Tgfu + uf
    # print(f"Trefftz embedding solve {time.time() - timer:5f} seconds")
    ndof = PP.shape[1]

    uh, ph = gfu.components[0:2]
    uerror = sqrt(Integrate(InnerProduct(uexact-uh,uexact-uh),mesh))
    perror = sqrt(Integrate(InnerProduct(pexact-ph,pexact-ph),mesh))
    return [uerror, perror, ndof]

def SolveStokesTDGEmbFES(mesh, k):
    """ 
    >>> SolveStokesTDGEmbFES(mesh2d, 5) # doctest:+ELLIPSIS
    [1...e-05, 0.001..., 529]
    """
    nu = 1.0
    zeta = cos(pi*x*(1-x)*y*(1-y))
    pexact = sin(pi*(x+y))
    uexact = CF((zeta.Diff(y), - zeta.Diff(x)))
    graduexact = CF((uexact.Diff(x),uexact.Diff(y)),dims=(2,2)).trans
    f1 = - nu*uexact[0].Diff(x).Diff(x) - nu*uexact[0].Diff(y).Diff(y) + pexact.Diff(x)
    f2 = - nu*uexact[1].Diff(x).Diff(x) - nu*uexact[1].Diff(y).Diff(y) + pexact.Diff(y)
    rhs = CF((f1,f2))

    V = VectorL2(mesh, order=k, dgjumps=True)
    Q = L2(mesh, order=k - 1, dgjumps=True)
    Z = NumberSpace(mesh)

    basefes = V * Q * Z
    u, p, z = basefes.TrialFunction()
    n = specialcf.normal(mesh.dim)

    ignoredofs = BitArray(basefes.ndof)
    ignoredofs[:] = False
    for i in range(V.ndof + Q.ndof, basefes.ndof):
        ignoredofs[i] = True
    Vs = VectorL2(mesh, order=k - 2)
    Qs = L2(mesh, order=k - 1)
    fes_test = Vs * Qs
    wu, wp = fes_test.TestFunction()[0:2]

    top = nu*InnerProduct( grad(u),grad(wu) ) * dx \
        - nu*InnerProduct(grad(u)*n,wu) * dx(element_boundary=True) \
        + InnerProduct(grad(p),wu) * dx + div(u)*wp*dx

    trhs = rhs * wu * dx(bonus_intorder=10)
    emb = TrefftzEmbedding(top=top, trhs=trhs, ignoredofs=ignoredofs)
    uf = emb.GetParticularSolution()
    fes = EmbeddedTrefftzFES(emb)

    ah,ch,fh = StokesDG(fes, nu, rhs)
    a = BilinearForm(fes)
    a += ah
    a.Assemble()
    # print(a.mat)
    # print((basefes.ndof-1)/mesh.ne, "+1")
    c = BilinearForm(fes)
    c += ch
    c.Assemble()
    f = LinearForm(fes)
    f += fh
    f.Assemble()

    res = f.vec.CreateVector()
    af = BilinearForm(basefes,fes)
    af += ah
    af.Apply(uf,res)
    f.vec.data -= res

    tfu = CGSolver(a.mat, c.mat.Inverse()) * f.vec
    gfu = GridFunction(basefes)
    gfu.vec.data = emb.Embed(tfu)
    gfu.vec.data += uf

    uh, ph = gfu.components[0:2]
    uerror = sqrt(Integrate(InnerProduct(uexact-uh,uexact-uh),mesh))
    perror = sqrt(Integrate(InnerProduct(pexact-ph,pexact-ph),mesh))
    return [uerror,perror,fes.ndof]

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

    u,v = fes.TnT()
    top = Lap(u)*Lap(v)*dx
    trhs = -rhs*Lap(v)*dx

    emb = TrefftzEmbedding(top=top,trhs=trhs,eps=eps)
    etfes = EmbeddedTrefftzFES(emb)
    uf = GridFunction(fes)
    uf.vec.data = emb.GetParticularSolution()
    # etfes = Compress(etfes, etfes.FreeDofs())
    a,f = dgell(etfes,exactpoi,rhs,uf)

    inv = a.mat.Inverse(inverse="sparsecholesky")
    gfu = GridFunction(etfes)
    gfu.vec.data = inv * f.vec
    nsol = gfu + uf

    return sqrt(Integrate((nsol-exactpoi)**2, mesh))



if __name__ == "__main__":
    import doctest
    doctest.testmod()
