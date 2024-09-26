
import numpy as np
import sys


class NodeInvalid:
    self.IsValid = False
    self.Id      = -1

class ElementInvalid:
    self.IsValid = False
    self.Id      = -1

class NodeAttrib:
    def __init__(self, x):
        self.IsValid = True
        self.Coord  = x
        self.Dim    = len(x)

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
        self.strain_e   = np.zeros(self.Dim2)    # elastic strain at node
        self.strain_p   = np.zeros(self.Dim2)    # plastic strain at node
        self.sigmaVM         = np.zeros(self.Dim2)    # VonMises stress

class ElementAttrib:
    def __init__(self, Connectivity):
        self.IsValid = True
        self.Connectivity = Connectivity
        self.NNode        = len(Connectivity)

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

    def GetNodeAtId(self, Id):
        for Node in self.Node:
            if Node.Id == Id:
                return Node

        tmp = "There is not such node with that id:" + str(Id)
        assert False, tmp




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
