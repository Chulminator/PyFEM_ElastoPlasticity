import numpy as np
import math
from FemBulk import *

def ApplyBC(Fem, Node, Element):
    #EBC =[]
    Dim    = Fem.Dimension
    for ind,x in enumerate(Node.Coord):
        Ndof = NodeDof(Node, [Node.Id[ind]])
        if math.fabs(x[0] - 0.0 ) < 1e-05:
            Node.BC_E[Ndof[0]*Dim] = 0

        if math.fabs(x[1] - 0.0 ) < 1e-05:
            Node.BC_E[Ndof[0]*Dim+1] = 0.0

        if math.fabs(x[2] - 0.0 ) < 1e-05:
            Node.BC_E[Ndof[0]*Dim+2] = 0.0

        if math.fabs(x[2] - 1.0 ) < 1e-05:
            Node.BC_E[Ndof[0]*Dim+2] = 0.0

        if math.fabs(x[0] - 10.0 ) < 1e-05:
            Node.BC_E[Ndof[0]*Dim] = 0.02/Fem.totalstep
    return
