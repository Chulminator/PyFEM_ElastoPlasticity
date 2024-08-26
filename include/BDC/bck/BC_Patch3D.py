import numpy as np
import math
from FemBulk import *

#def ApplyBC_PlateWithHole(Fem, Node, Element):
def ApplyBC(Fem, Node, Element):
    #print(Node.NNode)1
    #exit(1)
    Dim    = Fem.Dimension
    for ind, x in enumerate(Node.Coord):
        Ndof = NodeDof(Node, [Node.Id[ind]])
        ## uniaxial tension ###########################
        if math.fabs(x[2] - 0.0 ) < 1e-05:
            Node.BC_E[Ndof[0]*3+2] = 0.0
            #print("\t\t\t\t",Node.Id[ind])
            if math.fabs(x[0] - 0.0 ) < 1e-05:
                Node.BC_E[Ndof[0]*3] = 0.0
                #print("\t\t",Node.Id[ind])
            if math.fabs(x[1] - 0.0 ) < 1e-05:
                Node.BC_E[Ndof[0]*3+1] = 0.0
                #print("\t\t\t",Node.Id[ind])
        if math.fabs(x[2] - 1.0 ) < 1e-05:
            Node.BC_E[Ndof[0]*3+2] = 0.03/Fem.totalstep
            #print("\t",Node.Id[ind])
        ##############################################
    #exit(1)

        #if math.fabs(x[1] - 0.0 ) < 1e-05:
            #Node.BC_E[Ndof[0]*2+1] = 0.0
        #if math.fabs(x[0] - 0.0 ) < 1e-05:
            #Node.BC_E[Ndof[0]*2] = 0.0
        #if math.fabs(x[1] - 1.0 ) < 1e-05:
            #Node.BC_E[Ndof[0]*2] = 0.03/Fem.totalstep
    return
