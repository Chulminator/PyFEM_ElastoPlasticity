import numpy as np
import math
from FemBulk import *

#def ApplyBC_PlateWithHole(Fem, Node, Element):
def ApplyBC(Fem, Node, Element):
    Dim    = Fem.Dimension
    for ind, x in enumerate(Node.Coord):
        Ndof = NodeDof(Node, [Node.Id[ind]])

        if math.fabs(x[2] - 0.0 ) < 1e-05:      # For node with z = 0
            Node.BC_E[Ndof[0]*3+2] = 0.0        # z direction displacement is clamped
            if math.fabs(x[0] - 0.0 ) < 1e-05:  # For node with z = 0 and x = 0
                Node.BC_E[Ndof[0]*3] = 0.0      # x direction displacement is clamped
            if math.fabs(x[1] - 0.0 ) < 1e-05:  # For node with z = 0 and y = 0
                Node.BC_E[Ndof[0]*3+1] = 0.0    # y direction displacement is clamped

        if math.fabs(x[2] - 30.0 ) < 1e-05:     # For node with z = 30.0
            Node.BC_E[Ndof[0]*3+2] = 0.01/Fem.totalstep # z direction displacement is imposed by 0.01
    return
