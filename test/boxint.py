# %%
from ngsolve import *
from ngstrefftz import *
from ngsolve.meshes import MakeStructured2DMesh
from dg import *
SetNumThreads(3)

# ngsglobals.msg_level = 7
mesh = MakeStructured2DMesh(nx=1,ny=1, quads=False, mapping = lambda x,y: (x,y))
fes = H1(mesh, order=1, dirichlet=".*")
u,v = fes.TnT()
try:
    a = BilinearForm(fes)
    a += u*v*dbox
except Exception as e:
    print(e)
try:
    #integral = Integrate(CF(1.0)*dbox(element_boundary = True, reference_box_length=sqrt(2)/3), mesh, element_wise=True)
    # pcf = PrintCF('test.csv')
    integral = Integrate(CF(1.0)*dbox(element_boundary = False, reference_box_length=1/3), mesh, element_wise=True)
    # print("integral: ", integral)
except Exception as e:
    print(e)


def embtbox(mesh,order):
    """
    >>> [embtbox(mesh2d,5)] # doctest:+ELLIPSIS
    [...e-06]
    >>> embtbox(mesh3d,5) # doctest:+ELLIPSIS
    0.002...
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
    testfes = L2(mesh, order=order-2, dgjumps=True)
    etfes = EmbeddedTrefftzFES(fes)

    u = fes.TrialFunction()
    v = fes.TestFunction()
    if mesh.dim == 2:
        db = dbox(reference_box_length=1/3)
        op = -A*Lap(u)*v*db - CF((A.Diff(x),A.Diff(y)))*grad(u)*v*db + B*grad(u)*v*db + C*u*v*db
    elif mesh.dim == 3:
        db = dbox(reference_box_length=1/3)
        op = -A*Lap(u)*v*db - CF((A.Diff(x),A.Diff(y),A.Diff(z)))*grad(u)*v*db + B*grad(u)*v*db + C*u*v*db
    lop = rhs*v*db

    gfu = GridFunction(fes)
    gfu.vec.data = etfes.SetOp(op,lop,test_fes=testfes)

    a,f = dgell(etfes,Dbndc=exact,A=A,B=B,C=C,rhs=rhs,uf=gfu)

    tgfu = GridFunction(etfes)
    tgfu.vec.data = a.mat.Inverse() * f.vec
    gfu.vec.data += etfes.Embed(tgfu).vec 

    error = sqrt(Integrate((gfu-exact)**2, mesh))
    return error


if __name__ == "__main__":
    import doctest
    doctest.testmod()
