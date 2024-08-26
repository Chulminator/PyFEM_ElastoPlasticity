import numpy as np
import math
from FemBulk import *
from MeshHandle import *


# Geometry

      #oooooooooooo
      #____________
    #o|           |->    traction 100
    #o|           |->    traction 100
    #o|           |->    traction 100
    #o|___________|->    traction 100
      #oooooooooooo

def ApplyBC(Fem, Node, Element):
    Dim    = Fem.Dimension
    for ind, x in enumerate(Node.Coord):
        Ndof = NodeDof(Node, [Node.Id[ind]])
        if math.fabs(x[0] - 0.0 ) < 1e-05:  # Left boundary
            Node.BC_E[Ndof[0]*Dim] = 0.0

        if math.fabs(x[1] - 0.0 ) < 1e-05:  # Bottom boundary
            Node.BC_E[Ndof[0]*Dim+1] = 0.0

        if math.fabs(x[1] - 1.0 ) < 1e-05:  # Top boundary
            Node.BC_E[Ndof[0]*Dim+1] = 0.003 /Fem.totalstep

    Facet = GenerateFacet(Fem, Element)
    GP1d = [[-np.sqrt(1./3), 1],
            [np.sqrt(1./3), 1]]


    return

