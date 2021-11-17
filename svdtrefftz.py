from ngsolve import *
from tngs import *
import time
from netgen.geom2d import unit_square
import scipy as sp
import numpy as np
import matplotlib.pylab as plt

def trunc(values, decs=0):
    return np.trunc(values*10**decs)/(10**decs)
def spspy(sparsemat):
    rows,cols,vals = sparsemat.COO()
    A = sp.sparse.csr_matrix((vals,(rows,cols)))
    # A = A.todense()
    plt.spy(A); plt.plot(); plt.show()

def LocalP(A,eps=10**-9):
    global printer
    U,s,V = sp.linalg.svd(A)
    nz = (s<eps).sum()
    nnz = (s>=eps).sum()
    S=np.diag(s)
    P = U[:,nnz:]
    return P

def PySVDTrefftz(op,fes,eps=10**-10):
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
        LP = LocalP(elmat)
        P[i*localndof:(i+1)*localndof,i*trefftzndof:(i+1)*trefftzndof] = LP
    return P


order = 6
mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
n = specialcf.normal(mesh.dim)
exact = exp(x)*sin(y)
alpha = 4
n = specialcf.normal(2)
h = specialcf.mesh_size

fes = L2(mesh, order=order,  dgjumps=True)#,all_dofs_together=True)
u = fes.TrialFunction()
v = fes.TestFunction()

f = LinearForm(fes)
f += alpha*order**2/h*exact*v * ds(skeleton=True) \
     +(-n*grad(v)*exact)* ds(skeleton=True)

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
f.Assemble()

uh = u.Operator("hesse")
vh = v.Operator("hesse")
op = (uh[0,0]+uh[1,1])*(vh[0,0]+vh[1,1])*dx
# op = grad(u)*grad(v) * dx
# op = InnerProduct(u.Operator("hesse"),v.Operator("hesse"))*dx



# start = time.time()
# gfu = GridFunction(fes)
# # gfu.vec.data = np.linalg.solve(A,f.vec.FV())
# gfu.vec.data = a.mat.Inverse() * f.vec
# # print("time: ", time.time()-start)
# print ("L2-error:", sqrt(Integrate((gfu-exact)**2, mesh)))

# a.DeleteMatrix()
# a.Assemble()
start = time.time()
# SetNumThreads(2)
with TaskManager():
    PP = SVDTrefftz(op,fes,10**-9)
# spspy(PP)
PPT = PP.CreateTranspose()
na = PPT@a.mat@PP
tu = na.Inverse()*(PPT*f.vec)
tpgfu = GridFunction(fes)
tpgfu.vec.data = PP*tu
# print("trefftz time: ", time.time()-start)
print ("vs Trefftz L2-error:", sqrt(Integrate((tpgfu-exact)**2, mesh)))
for t in Timers():
    if 'svdt' in t['name']:
        print(t)


# P = PySVDTrefftz(op,fes)
# plt.spy(P); plt.plot(); plt.show()
# rows,cols,vals = a.mat.COO()
# A = sp.sparse.csr_matrix((vals,(rows,cols)))
# A = A.todense()
# TA = P.transpose()@A@P

# tsol = np.linalg.solve(TA,P.transpose()@f.vec.FV())
# tpgfu = GridFunction(fes)
# tpgfu.vec.data = P@tsol
# print ("vs PythonTrefftz L2-error:", sqrt(Integrate((tpgfu-exact)**2, mesh)))
# print(P.shape,"vs",PP.shape)


fes = FESpace("trefftzfespace",mesh,order=order,eq="laplace")
u = fes.TrialFunction()
v = fes.TestFunction()

f = LinearForm(fes)
f += alpha*order**2/h*exact*v * ds(skeleton=True) \
     +(-n*grad(v)*exact)* ds(skeleton=True)

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
f.Assemble()
gfu = GridFunction(fes)
gfu.vec.data = a.mat.Inverse() * f.vec
print ("L2-error:", sqrt(Integrate((gfu-exact)**2, mesh)))
