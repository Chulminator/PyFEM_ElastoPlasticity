import numpy as np
import math
from FemBulk import *
from MeshHandler import *


# Geometry

    # ^^^^^^^^^^^
    # |||||||||||
    # oooooooooooo
    # ____________
    #o|           |
    #o|           |
    #o|           |
    #o|___________|
    #^ooooooooooo


def ApplyBC(Model):
    Dim    = Model.Dim
    ModelNode = Model.Node
    for ind, Node in enumerate(ModelNode):
        if math.fabs(Node.Coord[1] - 0.0 ) < 1e-05:
            Node.BC_E[1] = 0.0

        if math.fabs(Node.Coord[0] - 0.0 ) < 1e-05:
            Node.BC_E[0] = 0.0

        if math.fabs(Node.Coord[1] - 1.0 ) < 1e-05:
            Node.BC_E[1] = 0.003 / Model.totalstep
    return
