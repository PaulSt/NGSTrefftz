from ngsolve import *
from ngsolve.comp import *
from ngsolve.TensorProductTools import *
import netgen.meshing as ngm
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from trefftzngs import *
from DGeq import *
#  import netgen.gui
import scipy as sp
import scipy.sparse.linalg
import scipy.linalg
import time


def GetFESTrefftz(mesh,c=1):
    return FESpace("trefftzfespace", mesh, order = 4, wavespeed = c, dgjumps=True, basistype=0)


c=2
initmesh = Mesh(SegMesh(4,0,1))
mesh = EvolveTents(initmesh,c,1)
