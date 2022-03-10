# import sys, os
# sys.path.append(os.path.join(os.path.dirname(sys.path[0])))
from ngsolve import *
from ngstrefftz import *
import time
from netgen.geom2d import unit_square
from netgen.csg import unit_cube

order = 5
exactlap = exp(x)*sin(y)
exactpoi = sin(x)*sin(y)
eps = 10**-8
mesh2d = Mesh(unit_square.GenerateMesh(maxh=0.3))
mesh3d = Mesh(unit_cube.GenerateMesh(maxh = 1))
SetNumThreads(3)

Lap = lambda u : sum(Trace(u.Operator('hesse')))

########################################################################
# L2
########################################################################
def dglap(fes,bndc,rhs=0):
    """
    >>> fes = L2(mesh2d, order=order,  dgjumps=True)#,all_dofs_together=True)
    >>> a,f = dglap(fes,exactlap)
    >>> gfu = GridFunction(fes)
    >>> gfu.vec.data = a.mat.Inverse() * f.vec
    >>> sqrt(Integrate((gfu-exactlap)**2, mesh2d)) # doctest:+ELLIPSIS
    3...e-09
    """
    mesh = fes.mesh
    order = fes.globalorder
    alpha = 4
    n = specialcf.normal(mesh.dim)
    h = specialcf.mesh_size
    u = fes.TrialFunction()
    v = fes.TestFunction()

    jump_u = u-u.Other()
    jump_v = v-v.Other()
    mean_dudn = 0.5*n * (grad(u)+grad(u.Other()))
    mean_dvdn = 0.5*n * (grad(v)+grad(v.Other()))

    a = BilinearForm(fes)
    a += grad(u)*grad(v) * dx \
        +alpha*order**2/h*jump_u*jump_v * dx(skeleton=True) \
        +(-mean_dudn*jump_v-mean_dvdn*jump_u) * dx(skeleton=True) \
        +alpha*order**2/h*u*v * ds(skeleton=True) \
        +(-n*grad(u)*v-n*grad(v)*u)* ds(skeleton=True)
    a.Assemble()

    f = LinearForm(fes)
    f += alpha*order**2/h*bndc*v * ds(skeleton=True) \
         +(-n*grad(v)*bndc)* ds(skeleton=True) \
         +rhs*v*dx
    f.Assemble()
    return a,f


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
    # >>> fes = L2(mesh2d, order=order,  dgjumps=True)#,all_dofs_together=True)
    # >>> u,v = fes.TnT()
    # >>> uh = u.Operator("hesse")
    # >>> vh = v.Operator("hesse")
    # >>> op = (uh[0,0]+uh[1,1])*(vh[0,0]+vh[1,1])*dx
    # >>> P = PySVDTrefftz(op,fes,eps)
    # >>> a,f = dglap(fes,exactlap)
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
    >>> fes = L2(mesh2d, order=order,  dgjumps=True)#,all_dofs_together=True)
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
    a,f = dglap(fes,exactlap)
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
    >>> fes = L2(mesh2d, order=order,  dgjumps=True)#,all_dofs_together=True)
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
    a,f = dglap(fes,exactlap)
    TA = PPT@a.mat@PP
    TU = TA.Inverse()*(PPT*f.vec)
    tpgfu = GridFunction(fes)
    tpgfu.vec.data = PP*TU
    return sqrt(Integrate((tpgfu-exactlap)**2, mesh))


def testembtrefftznonsym(fes):
    """
    >>> fes = L2(mesh2d, order=order,  dgjumps=True)#,all_dofs_together=True)
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
    a,f = dglap(fes,exactlap)
    TA = PPT@a.mat@PP
    TU = TA.Inverse()*(PPT*f.vec)
    tpgfu = GridFunction(fes)
    tpgfu.vec.data = PP*TU
    return sqrt(Integrate((tpgfu-exactlap)**2, mesh))

def testembtrefftzpoi(fes):
    """
    >>> fes = L2(mesh2d, order=order,  dgjumps=True)#,all_dofs_together=True)
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
    a,f = dglap(fes,exactpoi,rhs)
    TA = PPT@a.mat@PP
    TU = TA.Inverse()*(PPT*(f.vec-a.mat*ufv))
    tpgfu = GridFunction(fes)
    tpgfu.vec.data = PP*TU+ufv
    return sqrt(Integrate((tpgfu-exactpoi)**2, mesh))


def testembtrefftzpoi_mixed(fes):
    """
    >>> fes = L2(mesh2d, order=order,  dgjumps=True)#,all_dofs_together=True)
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
    a,f = dglap(fes,exactpoi,rhs)
    TA = PPT@a.mat@PP
    TU = TA.Inverse()*(PPT*(f.vec-a.mat*ufv))
    tpgfu = GridFunction(fes)
    tpgfu.vec.data = PP*TU+ufv
    return sqrt(Integrate((tpgfu-exactpoi)**2, mesh))



if __name__ == "__main__":
    import doctest
    doctest.testmod()
