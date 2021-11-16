from ngsolve import *
from tngs import *
import time
from netgen.geom2d import unit_square
import scipy as sp
import numpy as np
import matplotlib.pylab as plt

printer=0
def trunc(values, decs=0):
    return np.trunc(values*10**decs)/(10**decs)

def LocalP(A,eps=10**-9):
    global printer 
    U,s,V = sp.linalg.svd(A)
    nz = (s<eps).sum()
    nnz = (s>=eps).sum()
    # print(s)
    # print(nnz,nz, " vs ", 2*3+1)
    S=np.diag(s)
    # V = V[nnz:,nnz:]
    # U = U[:,nnz:]
    # P=U@S@V
    P = U[:,nnz:]
    if(printer==0):
        # print(A)
        # print(trunc(U@np.diag(s)@V,decs=0))
        # print(s<eps)
        # print(s)
        # print((s<eps)[nnz:])
        print("U \n",trunc(U,2))
        print("V \n",trunc(V,2))
        # print(trunc(U@S@V,2))
        # print("P",trunc(P))
        # print(P.shape)
        printer=1
    return P


order = 2
mesh = Mesh(unit_square.GenerateMesh(maxh=0.5))
fes = L2(mesh, order=order,  dgjumps=True)#,all_dofs_together=True)
n = specialcf.normal(mesh.dim)
exact = exp(x)*sin(y)
alpha = 4
n = specialcf.normal(2)
h = specialcf.mesh_size

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

opbf = BilinearForm(fes)
# opbf += grad(u)*grad(v) * dx
# opbf += InnerProduct(u.Operator("hesse"),v.Operator("hesse"))*dx
uh = u.Operator("hesse")
vh = v.Operator("hesse")
op = (uh[0,0]+uh[1,1])*(vh[0,0]+vh[1,1])*dx
opbf += op
opbf.Assemble()
rows,cols,vals = opbf.mat.COO()
A = sp.sparse.csr_matrix((vals,(rows,cols)))
# plt.spy(A); plt.plot(); plt.show()
A = A.todense()

localndof = int(fes.ndof/mesh.ne)
trefftzndof = 2*order+1
P = np.zeros([A.shape[1],(trefftzndof)*mesh.ne])
for i in range(mesh.ne):
    # U,s,V = sp.linalg.svd(A[i*localndof:(i+1)*localndof,i*localndof:(i+1)*localndof])
    elmat = A[i*localndof:(i+1)*localndof,i*localndof:(i+1)*localndof]
    # print(i," elmat=",trunc(elmat))
    LP = LocalP(elmat)
    P[i*localndof:(i+1)*localndof,i*trefftzndof:(i+1)*trefftzndof] = LP
    # input(i)

rows,cols,vals = a.mat.COO()
A = sp.sparse.csr_matrix((vals,(rows,cols)))
A = A.todense()
AA=P.transpose()@A@P

start = time.time()
tsol = np.linalg.solve(AA,P.transpose()@f.vec.FV())
tpgfu = GridFunction(fes)
tpgfu.vec.data = P@tsol
print("trefftz time: ", time.time()-start)

start = time.time()
gfu = GridFunction(fes)
gfu.vec.data = np.linalg.solve(A,f.vec.FV())#a.mat.Inverse() * f.vec
print("time: ", time.time()-start)


print ("L2-error:", sqrt(Integrate((gfu-exact)**2, mesh)))
print ("vs Trefftz L2-error:", sqrt(Integrate((tpgfu-exact)**2, mesh)))

plt.spy(P); plt.plot(); plt.show()
PP = SVDTrefftz(op,fes)
# rows,cols,vals = PP.COO()
# PP = sp.sparse.csr_matrix((vals,(rows,cols)))
# plt.spy(PP); plt.plot(); plt.show()
PPT = PP.CreateTranspose()
na = PPT@a.mat@PP
# tu = PP.CreateRowVector()
tu = na.Inverse()*(PPT*f.vec)
tpgfu = GridFunction(fes)
tpgfu.vec.data = PP*tu
print ("vs Trefftz L2-error:", sqrt(Integrate((tpgfu-exact)**2, mesh)))

# rows,cols,vals = na.COO()
# naa = sp.sparse.csr_matrix((vals,(rows,cols)))
# plt.spy(naa); plt.plot(); plt.show()
