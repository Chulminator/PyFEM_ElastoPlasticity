import numpy as np
import math
from FemBulk import *
from MeshHandler import *

def ApplyBC(Model):
    Dim    = Model.Dim
    for Node in Model.Node:
        if math.fabs(Node.Coord[1] - 1.75 ) < 1e-05:
            Node.BC_E[1] = 0.02/Model.totalstep

        if math.fabs(Node.Coord[1] + 1.75 ) < 1e-05:
            Node.BC_E[0] = 0.0

        if math.fabs(Node.Coord[1] + 1.75 ) < 1e-05:
            Node.BC_E[1] = 0.0

        if math.fabs(Node.Coord[1] + 1.75 ) < 1e-05:
            Node.BC_E[2] = 0.0
    return
