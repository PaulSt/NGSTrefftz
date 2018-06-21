from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from prodmesh import *
from ngsolve import *
import netgen.gui
from trefftzngs import *
from DGeq import *
import time
from ngsolve.TensorProductTools import *
from ngsolve.comp import *

mesh = Mesh(SegMesh(4,0,1))
tpmesh = NgsTPmesh(mesh,1)
Draw(tpmesh)




