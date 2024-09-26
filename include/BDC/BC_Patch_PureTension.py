import numpy as np
import math
from FemBulk import *
from MeshHandler import *

def ApplyBC(Model):
    Dim    = Model.Dim
    for ind, Node in enumerate(Model.Node):
        ## uniaxial tension ###########################
        if math.fabs(Node.Coord[1] - 0.0 ) < 1e-05:
            #Node.BC_E[0] = 0.0
            Node.BC_E[1] = 0.0

        if math.fabs(Node.Coord[0] - 0.0 ) < 1e-05:
            Node.BC_E[0] = 0.0

        if math.fabs(Node.Coord[0] - 1.0 ) < 1e-05:
            Node.BC_E[0] = 0.0

        if math.fabs(Node.Coord[1] - 1.0 ) < 1e-05:
            Node.BC_E[1] = 0.05/Model.totalstep
        ##############################################
    return
