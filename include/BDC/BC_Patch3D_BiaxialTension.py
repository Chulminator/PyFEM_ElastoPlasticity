import numpy as np
import math
from FemBulk import *
from MeshHandler import *

def ApplyBC(Model):
    Dim    = Model.Dim
    for ind, node in enumerate(Model.Node):
        if math.fabs(node.Coord[0] - 0.0 ) < 1e-05:
            node.BC_E[0] = 0.0

        if math.fabs(node.Coord[1] - 0.0 ) < 1e-05:
            node.BC_E[1] = 0.0

        if math.fabs(node.Coord[2] - 0.0 ) < 1e-05:
            node.BC_E[2] = 0.0

        if math.fabs(node.Coord[1] - 1.0 ) < 1e-05:
            node.BC_E[1] = 0.003/Model.totalstep

        if math.fabs(node.Coord[0] - 1.0 ) < 1e-05:
            node.BC_E[0] = 0.003/Model.totalstep
        ##############################################
    return
