import ngsolve
from ngsolve.fem import ET,IntegrationRule
from ngsolve.comp import L2,BilinearForm,LinearForm,SymbolicBFI,SymbolicLFI,FESpace
from ngsolve.fem import CoordCF
from .ngstrefftz import *

def GetWave(self,U):
    order=self.GetOrder()
    D = self.GetSpaceDim()
    initmesh = self.GetInitmesh()
    if D==3: eltyp = ET.TET
    elif D==2: eltyp = ET.TRIG
    elif D==1: eltyp = ET.SEGM
    intrule = IntegrationRule(eltyp,2*order)
    irsize = len(intrule.points)
    if U.dim==1:
        fes = L2(initmesh, order=order)
        u,v = fes.TnT()
        a = BilinearForm(fes)
        a += SymbolicBFI(u*v)
        a.Assemble()
        wavefront = self.GetWavefront()
        ipfct=IntegrationPointFunction(initmesh,intrule,wavefront)
        f = LinearForm(fes)
        f += SymbolicLFI(ipfct*v).SetIntegrationRule(eltyp,intrule)
        f.Assemble()
        U.vec.data = a.mat.Inverse(inverse='sparsecholesky') * f.vec
        # return gfu
    else:
        fes = L2(initmesh, order=order-1)**(D+1)
        u,v = fes.TnT()
        wavefront = self.GetWavefront()
        a = BilinearForm(fes)
        a += SymbolicBFI(u*v)
        a.Assemble()

        sirsize = int(wavefront.Width()/(D+1))
        f = LinearForm(fes)
        for d in range(D+1):
            ipfct=IntegrationPointFunction(initmesh,intrule,wavefront[:,d*sirsize:(d+1)*sirsize])
            f += SymbolicLFI(ipfct*v[d]).SetIntegrationRule(eltyp,intrule)
        f.Assemble()
        U.vec.data = a.mat.Inverse() * f.vec


dbox = BoxDifferentialSymbol()
dball = dbox(boxtype=BOXTYPE.BALL)

setattr(TWaveTents1, 'GetWave', GetWave)
setattr(TWaveTents2, 'GetWave', GetWave)
setattr(TWaveTents3, 'GetWave', GetWave)
setattr(QTWaveTents1, 'GetWave', GetWave)
setattr(QTWaveTents2, 'GetWave', GetWave)
