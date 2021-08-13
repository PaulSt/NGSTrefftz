from trefftzngs.trefftzngs import *
from trefftzngs.DGeq import *
from trefftzngs.prodmesh import *
from ngsolve import *
from ngsolve.fem import CoordCF
from types import MethodType

def GetWave(self):
    order=self.GetOrder()
    D = self.GetSpaceDim()
    initmesh = self.GetInitmesh()
    if D==3: eltyp = ET.TET
    elif D==2: eltyp = ET.TRIG
    elif D==1: eltyp = ET.SEGM
    intrule = IntegrationRule(eltyp,2*order)
    irsize = len(intrule.points)


    fes = L2(initmesh, order=order)
    u,v = fes.TnT()
    a = BilinearForm(fes)
    a += SymbolicBFI(u*v)
    a.Assemble()

    wavefront = self.GetWavefront()
    ipfct=IntegrationPointFunction(initmesh,intrule,wavefront)
    f = LinearForm(fes)
    f += SymbolicLFI(ipfct*v, intrule=intrule)
    f.Assemble()
    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse() * f.vec
    return gfu

setattr(WaveTents1, 'GetWave', GetWave)
setattr(WaveTents2, 'GetWave', GetWave)
setattr(WaveTents3, 'GetWave', GetWave)
setattr(QTWaveTents1, 'GetWave', GetWave)
setattr(QTWaveTents2, 'GetWave', GetWave)
# = MethodType(GetWave, None, TrefftzTent)
# TrefftzTent.GetWave = classmethod(GetWave)