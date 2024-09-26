import numpy as np
import math
from FemBulk import *
from MeshHandler import *

def ApplyBC(Model):
    Dim      = Model.Dim
    for ind, Node in enumerate(Model.Node):
        if math.fabs(Node.Coord[0] - 0.0 ) < 1e-05:  # Left boundary
            Node.BC_E[0] = 0.0
        if math.fabs(Node.Coord[1] - 0.0 ) < 1e-05:  # Bottom boundary
            Node.BC_E[1] = 0.0

    if Model.Facet == None:
        ModelFacet    = GenerateFacet(Model)
        Model.NFacet  = len(ModelFacet)
        Model.Facet   = ModelFacet

        GP1d = [[-np.sqrt(1./3.), 1],
                [np.sqrt(1./3.), 1]]

    if not (Model.Element[0].ElmentType == 'q4' or Model.Element[0].ElmentType == 't3'):
        assert False, "Boundary condition is not ready for the "+Fem.ElmentType+"Mesh"


    ApplyingTraction = 100.
    for facet in Model.Facet:
        FacetNode = facet.AdjNode
        x1 = FacetNode[0].Coord[0]
        y1 = FacetNode[0].Coord[1]
        x2 = FacetNode[1].Coord[0]
        y2 = FacetNode[1].Coord[1]

        if math.fabs(x1 - 10.0 ) < 1e-05 and math.fabs(x2 - 10.0 ) < 1e-05:
            length = np.sqrt((x1-x2)**2+(y1-y2)**2)
            for ii, gp in enumerate(GP1d):
                s      = gp[0]
                weight = gp[1]
                N = np.zeros(2)
                N[0] = (1.0 - s)*0.5
                N[1] = (1.0 + s)*0.5
                FacetNode[0].BC_N[1] += N[0] * ApplyingTraction * weight * Model.width * (length*0.5) / Model.totalstep
                FacetNode[1].BC_N[1] += N[1] * ApplyingTraction * weight * Model.width * (length*0.5) / Model.totalstep
    #print(str(ind)+"="*40)
    return

