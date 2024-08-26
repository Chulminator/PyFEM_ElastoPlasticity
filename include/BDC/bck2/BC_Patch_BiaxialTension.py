import numpy as np
import math
from FemBulk import *
from MeshHandle import *


# Geometry

      #^^^^^^^^^^^
      #|||||||||||
      #oooooooooooo
      #____________
    #o|           | o ->
    #o|           | o ->
    #o|           | o ->
    #o|___________| o ->
    #^ooooooooooo


def ApplyBC(Fem, Node, Element):
    Dim    = Fem.Dimension
    for ind, x in enumerate(Node.Coord):
        Ndof = NodeDof(Node, [Node.Id[ind]])
        if math.fabs(x[1] - 0.0 ) < 1e-05:
            Node.BC_E[Ndof[0]*Dim+1] = 0.0

        if math.fabs(x[0] - 0.0 ) < 1e-05:
            Node.BC_E[Ndof[0]*Dim] = 0.0

        if math.fabs(x[0] - 1.0 ) < 1e-05:
            Node.BC_E[Ndof[0]*Dim] = 0.003/Fem.totalstep

        if math.fabs(x[1] - 1.0 ) < 1e-05:
            Node.BC_E[Ndof[0]*Dim+1] = 0.003/Fem.totalstep
    return
