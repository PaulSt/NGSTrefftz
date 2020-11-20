from ngsolve import *
from trefftzngs import *
import numpy as np
import scipy as sp
import scipy.sparse.linalg
import scipy.linalg


def DGwaveeqsys(fes,U0,v0,sig0,c,gD,fullsys=False, applyrhs = False,alpha=0.5,beta=0.5,gamma=1):
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
    jump_Vt = (V - V.Other()) * n_t

    timelike = n_x*n_x # n_t=0
    spacelike = n_t**2 # n_x=0

    if(applyrhs == False):
        a = BilinearForm(fes)
        if(fullsys==True):
            HV = V.Operator("hesse")
            # a += SymbolicBFI(  -v*(-HV[0]+pow(c,-2)*HV[3]) ) #- sig*(-HV[1]+HV[2])  )
            a += SymbolicBFI(  -v*(-sum([HV[i*(D+2)] for i in range(D)]) + pow(c,-2)*HV[(D+1)*(D+1)-1]) )
        #space like faces, w/o x jump
        a += SymbolicBFI( spacelike * ( IfPos(n_t,v,vo)*(pow(c,-2)*jump_wt) + IfPos(n_t,sig,sigo)*(jump_taut) ), VOL, skeleton=True )
        # a += SymbolicBFI( spacelike * ( IfPos(n_t,v,vo)*(pow(c,-2)*jump_wt+jump_taux) + IfPos(n_t,sig,sigo)*(jump_wx+jump_taut) ), VOL, skeleton=True )
        #time like faces
        a += SymbolicBFI( timelike 	* ( mean_v*jump_taux + mean_sig*jump_wx + alpha*jump_vx*jump_wx + beta*jump_sigx*jump_taux ), VOL, skeleton=True )
        #t=T (or *x)
        a += SymbolicBFI( ( pow(c,-2)*v*w + sig*tau ), BND, definedon=fes.mesh.Boundaries("outflow"), skeleton=True)
        #dirichlet boundary 'timelike'
        a += SymbolicBFI( ( sig*n_x*w + alpha*v*w ), BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True)
        #correction term to recover sol of second order system
        a += SymbolicBFI( spacelike * ( gamma * jump_Ut*jump_Vt ), VOL, skeleton=True )
        # a += SymbolicBFI( spacelike * ( gamma * (n_t*jump_Ut)*IfPos(n_t,V.Other(),V) ), VOL, skeleton=True )
        # a += SymbolicBFI( spacelike * ( gamma * (-n_t*jump_Ut)*IfPos(n_t,V.Other(),V) ), VOL, skeleton=True )
        # a += SymbolicBFI( spacelike * ( gamma * (-jump_Ut)*IfPos(n_t,V.Other(),V) ), VOL, skeleton=True )
        a += SymbolicBFI( ( gamma * U*V ), BND, definedon=fes.mesh.Boundaries("inflow"), skeleton=True )
        a.Assemble()

        f = LinearForm(fes)
        f += SymbolicLFI( ( pow(c,-2)*v0*w + sig0*tau ), BND, definedon=fes.mesh.Boundaries("inflow"), skeleton=True) #t=0 (or *(1-x))
        f += SymbolicLFI( ( gD * (alpha*w - tau*n_x) ), BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True) #dirichlet boundary 'timelike'
        f += SymbolicLFI( gamma * ( (U0)*V ), BND, definedon=fes.mesh.Boundaries("inflow"),  skeleton=True ) #rhs correction term to recover sol of second order system
        f.Assemble()
    else:
        a = BilinearForm(fes)
        if(fullsys==True):
            HV = V.Operator("hesse")
            a += SymbolicBFI(  -v*(-sum([HV[i*(D+2)] for i in range(D)]) + pow(c,-2)*HV[(D+1)*(D+1)-1]) )
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

    return [a,f]


def DG1heateqsys(fes,U0,v0,sig0,c,gD,fullsys=False, applyrhs = False,alpha=0.5,beta=0.5,gamma=1):
    D = fes.mesh.dim - 1
    U = fes.TrialFunction()
    V = fes.TestFunction()
    gU = grad(U)
    gV = grad(V)

    v = U
    sig = CoefficientFunction(tuple([-gU[i] for i in  range(D)]))
    w = V
    tau = CoefficientFunction(tuple([-gV[i] for i in  range(D)]))

    vo = U.Other()
    sigo = CoefficientFunction(tuple([-gU.Other()[i] for i in  range(D)]))
    wo = V.Other()
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
    jump_Vt = (V - V.Other()) * n_t

    timelike = n_x*n_x # n_t=0
    spacelike = n_t**2 # n_x=0

    a = BilinearForm(fes)
    a += SymbolicBFI( 2*sig*tau )
    # a += SymbolicBFI(  (mean_v + beta*jump_sigx)*(jump_wt + jump_taux) + (mean_sig+alpha*jump_vx)*jump_wx , VOL, skeleton=True )
    # a += SymbolicBFI(  (mean_v + beta*(jump_sigx+jump_vt))*(jump_wt + jump_taux) + (mean_sig+alpha*(jump_vx+jump_sigt))*jump_wx , VOL, skeleton=True )

#split time space faces
    a += SymbolicBFI(  timelike * ((mean_v + beta*jump_sigx)*(jump_wt + jump_taux) + (mean_sig+alpha*jump_vx)*jump_wx) , VOL, skeleton=True )
    a += SymbolicBFI(  spacelike * (IfPos(n_t,v,vo)*(jump_wt + jump_taux) + (IfPos(n_t,sig,sigo))*jump_wx) , VOL, skeleton=True )
# end split time space faces

    a += SymbolicBFI( v*w, BND, definedon=fes.mesh.Boundaries("outflow"), skeleton=True)
    #dirichlet boundary 'timelike'
    a += SymbolicBFI( sig*n_x*w + alpha*v*w, BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True)
    #correction term to recover sol of second order system
    # a += SymbolicBFI( ( gamma * jump_Ut*jump_Vt ), VOL, skeleton=True )
    # a += SymbolicBFI( ( gamma * U*V ), BND, definedon=fes.mesh.Boundaries("inflow"), skeleton=True )
    a.Assemble()

    f = LinearForm(fes)
    f += SymbolicLFI( ( v0*w ), BND, definedon=fes.mesh.Boundaries("inflow"), skeleton=True) #t=0 (or *(1-x))
    f += SymbolicLFI( ( gD * (alpha*w - tau*n_x) ), BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True) #dirichlet boundary 'timelike'
    # f += SymbolicLFI( gamma * ( U0*V ), BND, definedon=fes.mesh.Boundaries("inflow"),  skeleton=True ) #rhs correction term to recover sol of second order system
    f.Assemble()

    return [a,f]



def DGheateqsys(fes,U0,v0,sig0,c,gD,fullsys=False, applyrhs = False,alpha=0.5,beta=0.5,gamma=1):
    D = fes.mesh.dim - 1
    u = fes.TrialFunction()
    v = fes.TestFunction()
    gradu = grad(u)
    gradv = grad(v)

    dtu = gradu[D]
    dtv = gradv[D]
    gu = CoefficientFunction(tuple([gradu[i] for i in  range(D)]))
    gv = CoefficientFunction(tuple([gradv[i] for i in  range(D)]))

    uo = u.Other()
    vo = v.Other()
    guo = CoefficientFunction(tuple([gradu.Other()[i] for i in  range(D)]))
    gvo = CoefficientFunction(tuple([gradv.Other()[i] for i in  range(D)]))

    h = specialcf.mesh_size
    n = specialcf.normal(D+1)
    n_t = n[D]/Norm(n)
    n_x = CoefficientFunction( tuple([n[i]/Norm(n) for i in  range(D)]) )

    mean_u = 0.5*(u+uo)
    mean_v = 0.5*(v+vo)
    mean_gu = 0.5*(gu+guo)
    mean_gv = 0.5*(gv+gvo)

    jump_ux = ( u - uo ) * n_x
    jump_vx = ( v - vo ) * n_x
    jump_gux = ( gu - guo ) * n_x
    jump_gvx = ( gv - gvo ) * n_x

    jump_ut = ( u - uo ) * n_t
    jump_vt = ( v - vo ) * n_t

    timelike = n_x*n_x # n_t=0
    spacelike = n_t**2 # n_x=0

##georgoulis
    # a = BilinearForm(fes)
    # a += SymbolicBFI( dtu*v + gu*gv )

    # a += SymbolicBFI(  -timelike * ( mean_gu*jump_vx + mean_gv*jump_ux + alpha*jump_ux*jump_vx) , VOL, skeleton=True )
    # a += SymbolicBFI(  spacelike * ( jump_ut*IfPos(n_t,v,vo) ) , VOL, skeleton=True )

    # a += SymbolicBFI( -u*v, BND, definedon=fes.mesh.Boundaries("outflow"), skeleton=True)
    # a += SymbolicBFI( -u*n_x*gv + alpha*u*v, BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True)
    # a += SymbolicBFI(  u*v , BND, definedon=fes.mesh.Boundaries("inflow"), skeleton=True )
    # a.Assemble()

    # f = LinearForm(fes)
    # f += SymbolicLFI( ( U0*v ), BND, definedon=fes.mesh.Boundaries("inflow"), skeleton=True) #t=0 (or *(1-x))
    # f += SymbolicLFI( ( gD * (alpha*v - gv*n_x) ), BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True)
    # f.Assemble()



##presi
    a = BilinearForm(fes)
    a += SymbolicBFI( - u*dtv + gu*gv )

    a += SymbolicBFI( ( -mean_gu*jump_vx ) , VOL, skeleton=True )
    a += SymbolicBFI( ( -mean_gv*jump_ux ) , VOL, skeleton=True )
    # a += SymbolicBFI( ( -0.5*jump_gux*jump_vx ) , VOL, skeleton=True )
    # a += SymbolicBFI( ( alpha/h*jump_ux*jump_vx) , VOL, skeleton=True )
    a += SymbolicBFI( ( IfPos(n_t,u,uo)*jump_vt) , VOL, skeleton=True )

    a += SymbolicBFI( u*v, BND, definedon=fes.mesh.Boundaries("outflow"), skeleton=True)
    # a += SymbolicBFI( u*v, BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True)
    # a += SymbolicBFI( u*v , BND, definedon=fes.mesh.Boundaries("inflow"), skeleton=True )
    a.Assemble()

    f = LinearForm(fes)
    f += SymbolicLFI( ( U0*v ), BND, definedon=fes.mesh.Boundaries("inflow"), skeleton=True) #t=0 (or *(1-x))
    f += SymbolicLFI( ( gD*n_x *v ), BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True)
    f.Assemble()

##new
    Hv = v.Operator("hesse")
    a = BilinearForm(fes)
    a += SymbolicBFI( - u*(dtv + Hv[0,0] ) )

    a += SymbolicBFI( ( mean_u*jump_gvx ) , VOL, skeleton=True )
    a += SymbolicBFI( ( mean_v*jump_gux ) , VOL, skeleton=True )
    # a += SymbolicBFI( ( 0.5*jump_ux*jump_gvx ) , VOL, skeleton=True )
    # a += SymbolicBFI( ( 0.5*jump_vx*jump_gux ) , VOL, skeleton=True )

    a += SymbolicBFI( ( -mean_gu*jump_vx ) , VOL, skeleton=True )
    a += SymbolicBFI( ( -mean_gv*jump_ux ) , VOL, skeleton=True )
    # a += SymbolicBFI( ( -jump_gux*jump_vx ) , VOL, skeleton=True )
    # a += SymbolicBFI( ( alpha/h*jump_ux*jump_vx) , VOL, skeleton=True )
    a += SymbolicBFI( ( IfPos(n_t,u,uo)*jump_vt) , VOL, skeleton=True )

    a += SymbolicBFI( u*v, BND, definedon=fes.mesh.Boundaries("outflow"), skeleton=True)
    # a += SymbolicBFI( u*v, BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True)
    # a += SymbolicBFI( u*v , BND, definedon=fes.mesh.Boundaries("inflow"), skeleton=True )
    # f += SymbolicLFI( ( u*n_x*gv ), BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True)
    a.Assemble()

    f = LinearForm(fes)
    f += SymbolicLFI( ( U0*v ), BND, definedon=fes.mesh.Boundaries("inflow"), skeleton=True) #t=0 (or *(1-x))
    f += SymbolicLFI( ( gD*n_x *v ), BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True)
    f.Assemble()

# new first order
    sig = CoefficientFunction(tuple([-gu[i] for i in  range(D)]))
    tau = CoefficientFunction(tuple([-gv[i] for i in  range(D)]))
    sigo = CoefficientFunction(tuple([-gradu.Other()[i] for i in  range(D)]))
    tauo = CoefficientFunction(tuple([-gradv.Other()[i] for i in  range(D)]))
    mean_sig = 0.5*(sig+sigo)
    mean_tau = 0.5*(tau+tauo)
    jump_sigx = (( sig - sigo ) * n_x)
    jump_taux = (( tau - tauo ) * n_x)


    Hv = v.Operator("hesse")
    a = BilinearForm(fes)
    # a += SymbolicBFI( - u*(dtv - Hv[0,0] ) + sig*(tau-gv) )
    a += SymbolicBFI( - u*(dtv + Hv[0,0] ) )

    a += SymbolicBFI( ( IfPos(n_t,u,uo)*jump_vt) , VOL, skeleton=True )

    a += SymbolicBFI( ( -mean_gu*jump_vx ) , VOL, skeleton=True )
    # a += SymbolicBFI( ( -mean_gv*jump_ux ) , VOL, skeleton=True )
    a += SymbolicBFI( ( mean_u*jump_gvx ) , VOL, skeleton=True )
    # a += SymbolicBFI( ( mean_v*jump_gux ) , VOL, skeleton=True )

    a += SymbolicBFI( u*v, BND, definedon=fes.mesh.Boundaries("outflow"), skeleton=True)
    a += SymbolicBFI( ( u*n_x*gv ), BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True)
    a.Assemble()

    f = LinearForm(fes)
    f += SymbolicLFI( ( U0*v ), BND, definedon=fes.mesh.Boundaries("inflow"), skeleton=True) #t=0 (or *(1-x))
    f += SymbolicLFI( ( gD*n_x*v ), BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True)
    f.Assemble()


    return [a,f]



def DGheateqsysnew(fes,fes2,U0,v0,sig0,c,gD,fullsys=False, applyrhs = False,alpha=0.5,beta=0.5,gamma=1):
    D = fes.mesh.dim - 1
    u = fes.TrialFunction()
    v = fes2.TestFunction()
    gradu = grad(u)
    gradv = grad(v)

    dtu = gradu[D]
    dtv = gradv[D]
    gu = CoefficientFunction(tuple([gradu[i] for i in  range(D)]))
    gv = CoefficientFunction(tuple([gradv[i] for i in  range(D)]))

    uo = u.Other()
    vo = v.Other()
    guo = CoefficientFunction(tuple([gradu.Other()[i] for i in  range(D)]))
    gvo = CoefficientFunction(tuple([gradv.Other()[i] for i in  range(D)]))

    h = specialcf.mesh_size
    n = specialcf.normal(D+1)
    n_t = n[D]/Norm(n)
    n_x = CoefficientFunction( tuple([n[i]/Norm(n) for i in  range(D)]) )

    mean_u = 0.5*(u+uo)
    mean_v = 0.5*(v+vo)
    mean_gu = 0.5*(gu+guo)
    mean_gv = 0.5*(gv+gvo)

    jump_ux = ( u - uo ) * n_x
    jump_vx = ( v - vo ) * n_x
    jump_gux = ( gu - guo ) * n_x
    jump_gvx = ( gv - gvo ) * n_x

    jump_ut = ( u - uo ) * n_t
    jump_vt = ( v - vo ) * n_t

    timelike = n_x*n_x # n_t=0
    spacelike = n_t**2 # n_x=0

    Hv = v.Operator("hesse")
    a = BilinearForm(trialspace=fes,testspace=fes2)
    # a += SymbolicBFI( - u*(dtv - Hv[0,0] ) + sig*(tau-gv) )
    # a += SymbolicBFI( - u*(dtv + Hv[0,0] ) )

    a += SymbolicBFI( ( IfPos(n_t,u,uo)*jump_vt) , VOL, skeleton=True )

    a += SymbolicBFI( ( -mean_gu*jump_vx ) , VOL, skeleton=True )
    a += SymbolicBFI( ( alpha*jump_ux*jump_vx ) , VOL, skeleton=True )
    # a += SymbolicBFI( ( -mean_gv*jump_ux ) , VOL, skeleton=True )
    a += SymbolicBFI( ( mean_u*jump_gvx ) , VOL, skeleton=True )
    a += SymbolicBFI( ( -alpha*jump_gux*jump_gvx ) , VOL, skeleton=True )
    # a += SymbolicBFI( ( mean_v*jump_gux ) , VOL, skeleton=True )

    a += SymbolicBFI( u*v, BND, definedon=fes.mesh.Boundaries("outflow"), skeleton=True)
    a += SymbolicBFI( ( u*n_x*gv ), BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True)
    a.Assemble()

    f = LinearForm(fes2)
    f += SymbolicLFI( ( U0*v ), BND, definedon=fes.mesh.Boundaries("inflow"), skeleton=True) #t=0 (or *(1-x))
    f += SymbolicLFI( ( gD*n_x*v ), BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True)
    f.Assemble()

    return [a,f]



def DGsolve(fes,a,f):
    gfu = GridFunction(fes, name="uDG")

    # tmp1 = f.vec.CreateVector()
    # tmp2 = f.vec.CreateVector()
    # def matvec(v):
            # tmp1.FV().NumPy()[:] = v
            # tmp2.data = a.mat * tmp1
            # return tmp2.FV().NumPy()
    # A = sp.sparse.linalg.LinearOperator( (a.mat.height,a.mat.width), matvec)
    # gfu.vec.FV().NumPy()[:], succ = sp.sparse.linalg.gmres(A, f.vec.FV().NumPy())

    rows,cols,vals = a.mat.COO()
    A = sp.sparse.csr_matrix((vals,(rows,cols)))
    # gfu.vec.FV().NumPy()[:] = sp.sparse.linalg.spsolve(A,f.vec.FV())

    # A = A.todense()
    # nmatclean = A[~(A==0).all(1).A1,:]
    # nmatclean = nmatclean[:,~(A==0).all(1).A1]
    # nvecclean = f.vec.FV().NumPy()[~(A==0).all(1).A1]
    # solclean = np.linalg.solve(nmatclean,nvecclean)
    # sol = np.zeros(a.mat.height)
    # sol[~(A==0).all(1).A1] = solclean
    # gfu.vec.FV().NumPy()[:] = sol

    cond = np.linalg.cond(A.todense())
    # cond = 0

    # print(cond)
    # print(a.mat.height)
    # print(a.mat.width)
    # nmat = np.zeros((a.mat.height,a.mat.width))
    # nvec = np.zeros(a.mat.width)

    # for i in range(a.mat.width):
            # nvec[i] = f.vec[i]/sqrt(a.mat[i,i])
            # for j in range(a.mat.height):
                    # nmat[j,i] = a.mat[j,i]/sqrt(a.mat[i,i]*a.mat[j,j])
    # sol = scipy.linalg.solve(nmat,nvec)

    # for i in range(a.mat.height):
            # gfu.vec[i] = sol[i]/sqrt(a.mat[i,i])
    # cond = np.linalg.cond(nmat)

    #print(min(np.linalg.eigvalsh(0.5*(nmat + nmat.transpose()))))
    #nmatinv = np.linalg.inv(nmat)
    #print(min(np.linalg.eigvalsh(0.5*(nmatinv + nmatinv.transpose()))))

    gfu.vec.data = a.mat.Inverse()*f.vec

    return [gfu,cond]
