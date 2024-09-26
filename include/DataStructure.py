# ----------------------------------------------------------------
# Written by: CK in promechanics.org at Columbia University
# ----------------------------------------------------------------
import numpy as np
import sys
import time as ct
import os.path


class NodeAttrib:
    def __init__(self, x):
        self.Coord  = x
        self.Dim    = len(x)
        self.NElem  = 0

    #def __del__(self):
        #print("Node ",self.Id," is deleted")

    def SetId(self, Id):
        self.Id = Id

    def InitAttrib(self):
        # https://www.continuummechanics.org/principalstressesandstrains.html
        # https://en.wikipedia.org/wiki/Von_Mises_yield_criterion
        self.u          = np.zeros(self.Dim)  # total displacement
        self.u1         = np.zeros(self.Dim)  # total Displacement1
        self.du         = np.zeros(self.Dim)  # displacement difference
        self.u_e        = np.zeros(self.Dim)  # elastic displacement
        self.u_p        = np.zeros(self.Dim)  # plastic displacement
        self.BC_E       = np.empty(self.Dim)  # Essential boundary condtion
        self.BC_E[:]    = np.NaN
        self.BC_N       = np.zeros(self.Dim)  # Natural boundary condtion
        self.BC_N_init  = np.zeros(self.Dim)  # Initla natural boundary condtion
        self.F_int      = np.zeros(self.Dim)  # Internal force
        self.F_ext      = np.zeros(self.Dim)  # External force
        if self.Dim == 2:
            self.Dim2 = 3
        elif self.Dim == 3:
            self.Dim2 = 6
        else:
            assert False, "Wrong dimension"

        self.stress     = np.zeros(self.Dim2)    # stress
        self.sigma1     = 0.0                    # principal stress
        #self.strain_e   = np.zeros(self.Dim2)    # elastic strain at node
        #self.strain_p   = np.zeros(self.Dim2)    # plastic strain at node
        self.sigmaVM         = np.zeros(self.Dim2)    # VonMises stress



class EdgeAttribute:
    def __init__(self, AdjNode, AdjFacet):
        self.AdjNode  = AdjNode
        self.AdjFacet = AdjFacet

class FacetAttribute:
    def __init__(self, AdjNode, AdjElem):
        self.AdjElem = AdjElem
        self.AdjNode = AdjNode

class ElementAttrib:
    def __init__(self, Connectivity):
        self.Connectivity = Connectivity
        self.NNode        = len(Connectivity)

    #def __del__(self):
        #print("Element ",self.Id," is deleted")

    def SetId(self, Id):
        self.Id = Id

    def InitAttrib(self):
        self.Stiffness   = np.zeros((self.NNode,self.NNode))
        self.B_matrix    = []
        self.Area        = 0.0

    def SetMatProp(self, MatProp):
        self.MatProp = MatProp
        self.E       = MatProp[0]
        self.nu      = MatProp[1]
        E  = self.E
        nu = self.nu
        self.lamda = E*nu / ((1.0+nu)*(1.0-2.0*nu))
        self.K     = E / (3*(1.0-2.0*nu))
        self.mu    = E / (2.0*(1.0+nu))

    def SetElementType(self, Dim=2):
        self.G_Edof       = Dim * self.NNode
        self.Dim          = Dim

        if len(self.Connectivity) == 4 and Dim == 2:
            self.ElmentType = 'q4'
            self.GP = [[-1./np.sqrt(3), -1./np.sqrt(3), 1.],
                        [ 1./np.sqrt(3), -1./np.sqrt(3), 1.],
                        [ 1./np.sqrt(3),  1./np.sqrt(3), 1.],
                        [-1./np.sqrt(3),  1./np.sqrt(3), 1.]]

        elif len(self.Connectivity) == 3:
            self.ElmentType = 't3'
            self.GP = [[1./3., 1./3., 1./2.]]

        elif len(self.Connectivity) == 8:
            self.ElmentType = 'hex8'
            self.GP = [[-1./np.sqrt(3.), -1./np.sqrt(3.),  1./np.sqrt(3), 1.],
                        [ 1./np.sqrt(3.), -1./np.sqrt(3.),  1./np.sqrt(3), 1.],
                        [ 1./np.sqrt(3.),  1./np.sqrt(3.),  1./np.sqrt(3), 1.],
                        [-1./np.sqrt(3.),  1./np.sqrt(3.),  1./np.sqrt(3), 1.],
                        [-1./np.sqrt(3.), -1./np.sqrt(3.), -1./np.sqrt(3), 1.],
                        [ 1./np.sqrt(3.), -1./np.sqrt(3.), -1./np.sqrt(3), 1.],
                        [ 1./np.sqrt(3.),  1./np.sqrt(3.), -1./np.sqrt(3), 1.],
                        [-1./np.sqrt(3.),  1./np.sqrt(3.), -1./np.sqrt(3), 1.]]

        elif len(self.Connectivity) == 6:
            self.ElmentType = 't6'
            self.GP = [[1./6., 1./6., 1./6.],
                        [2./3., 1./6., 1./6.],
                        [1./6., 2./3., 1./6.]]
            assert False,"T6 is not ready yet"

        elif len(self.Connectivity) == 4 and Dim == 3:
            self.ElmentType = 'tet4'
            self.GP = [[1./4., 1./4., 1./4., 1./6.]]

        else:
            print("Current avaiable elemet types are following:")
            print("\t T3 with 1 Gauss quadrature point")
            print("\t Q4 with 4 Gauss quadrature points")
            print("\t Tet4 with 1 Gauss quadrature point")
            print("\t Hex8 with 8 Gauss quadrature points")
            print("\t T6 is under construction")
            assert False, "Check the element type"

class ModelAttrib:
    # pqspace
    # eigenspace
    # tensorspace

    def __init__(self, ModelName):
        self.ModelName = ModelName
        print("\t\tModel "+self.ModelName+" is generated")

        self.NElem: int = 0
        self.NNode: int = 0
        self.Dim: int   = 0
        self.HSP: float = 0.
        self.ReturnMappingMode: str ='eigenspace'
        self.twoD: str = ''
        self.Node    = []
        self.Element = []
        self.GlobalNRStep    = 0
        self.Facet = None
        self.Edge  = None

    #def __del__(self):
        #print("****Mesh is deleted****")

    def GetNodeAtId(self, Id):
        for node in self.Node:
            if node.Id == Id:
                return node

        tmp = "There is not such node with that id:" + str(Id)
        assert False, tmp

    def GetElementAtId(self, Id):
        for elem in self.Element:
            if elem.Id == Id:
                return elem

        tmp = "There is not such elem with that id:" + str(Id)
        assert False, tmp

    def NodeIsValid(self, Node):
        if Node == self.NodeInvalid:
            return False
        else:
            return True

    def ElementIsValid(self, Elem):
        if Elem == self.ElementInvalid:
            return False
        elif Elem == ElementVoid:
            return False
        else:
            return True

    def GetNodeDof(self, Connectivity):
        Ndof = []
        flag = -1
        for ElemNode in Connectivity:
            flag = -1
            for ind, Node in enumerate(self.Node):
                if Node == ElemNode:
                    Ndof.append(ind)
                    flag = 1
                    break
            if flag == -1:
                assert False, "Invalid Node is being handled"
        return Ndof

    def ElementalNodeDof(self, Connectivity):
        ENdof = []
        flag = -1
        for ElemNode in Connectivity:
            flag = -1
            for ind, Node in enumerate(self.Node):
                if Node == ElemNode:
                    for ii in range(self.Dim):
                        ENdof.append(ind*self.Dim+ii)
                    flag = 1
                    break
            if flag == -1:
                assert False, "Invalid Node is being handled"
        return ENdof

    def RemNode(self, SingleNode):
        if self.Facet != None:
            del self.Facet
            del self.Edge
        print("Remove node",SingleNode.Id)

        RemNodeId = -1
        for ii, node in enumerate(self.Node):
            if node == SingleNode:
                RemNodeId = ii
                break

        RemElemId = []
        for ii, elem in enumerate(self.Element):
            flag = 0
            for jj, node in enumerate(elem.Connectivity):
                if node == SingleNode:
                    flag = 1
                    RemElemId.append([ii,jj])
                    break
            if flag == 1:
                for node in elem.Connectivity:
                    node.NElem -= 1

            if SingleNode.NElem == 0:
                break
        del self.Node[RemNodeId]
        for ii in reversed(RemElemId):
            print("Element",self.Element[ii[0]].Id, "is deleted")
            del self.Element[ii[0]]
        del SingleNode
        return 1


    def RemElem(self, SingleElem):
        if self.Facet != None:
            del self.Facet
            del self.Edge
        print("Remove element",SingleElem.Id)

        RemElemId = -1
        for ii, element in enumerate(self.Element):
            if element == SingleElem:
                RemElemId = ii
                break
        RemNode = []
        for node in self.Element[RemElemId].Connectivity:
            node.NElem -= 1
            if node.NElem == 0:
                RemNode.append(node)
        RemNodeId = []

        for ii, node in enumerate(self.Node):
            if node in RemNode:
                RemNodeId.append(ii)

        del SingleElem
        for ii in reversed(RemNodeId):
            del self.Node[ii]
        del self.Element[RemElemId]



    class FacetInvalid:
        Id      = -1

    class NodeInvalid:
        Id      = -1

    class ElementInvalid:
        ##IsBoundary = True # need to think about this
        Id      = -1
    #class ElementVoid: #IsBoundary
        ##IsBoundary = True # need to think about this
        #Id      = -1

    #class CheckConsistency:
        #do something


class ProgramAttrib:
    #def __init__(self, ProgramName):
        #self.ProgramName = ProgramName
        #print("\tModel "+self.ProgramName+" is generated\n")

    def LogGen(self, logname):
        self.file_log = open(logname,'w')

    def LogWrite(self, command):
        self.file_log.write(command)

    def LogClose(self):
        self.file_log.close

    def show(self):
        print("============= Program description =============")
        print("Model")
        print("\ttitle   : ",self.title)
        print("\tresult  : ",self.result,"\n")
        print("Solid")
        print("\tE    = ",self.E,"\tElastic Modulus")
        print("\tv    = ",self.nu,"\t\tPoisson ratio")
        print("\tlamda= ",self.lamda,"\t\t1st lame parameter")
        print("\tmu   = ",self.mu,"\t\t2st lame parameter \n")
        print("============= Program description =============")

    def PrintCommand(self, command, Nindent):
        command.strip('\n')
        for i in range(0,Nindent):
            command = '\t' + command
        print(command)
        self.LogWrite(command+'\n')
        return
