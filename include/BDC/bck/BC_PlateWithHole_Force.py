import numpy as np
import math
from FemBulk import *

#def ApplyBC_PlateWithHole(Fem, Node, Element):
def ApplyBC(Fem, Node, Element):
    EBC =[]
    Dim    = Fem.Dimension
    for ind,x in enumerate(Node.Coord):
        Ndof = NodeDof(Node, [Node.Id[ind]])
        if math.fabs(x[0] - 0.0 ) < 1e-05:
            Node.BC_E[Ndof[0]*Dim] = 0
        if math.fabs(x[1] - 0.0 ) < 1e-05:
            Node.BC_E[Ndof[0]*Dim+1] = 0.0
        if math.fabs(x[0] - 10.0 ) < 1e-05:
            #if math.fabs(x[0] - 0.0 ) < 1e-05:
                #Node.BC_N[Ndof[0]*Dim] = 100.0/Fem.totalstep/2.0
            #if math.fabs(x[1] - 0.0 ) < 1e-05:
                #Node.BC_N[Ndof[0]*Dim] = 100.0/Fem.totalstep/2.0
            #else:
            Node.BC_N[Ndof[0]*Dim] = 100.0/Fem.totalstep
    return
