# ---------------------------------------------------------------- 
# Written by: CK, HS in promechanics.org at Columbia University
# ----------------------------------------------------------------       
import numpy as np
import sys
import time as ct
import os.path
from dataclasses import dataclass

# Things to upgrade
#   PlotSetup   : what variable will be plotted how often
#   OutputSetup : what variable will be written how often



def ReadInput(input_name, Fem):
    #input_name = "./input/terzaghi/NonpolarElasticity.inp"
    print('#readinput')
    #Fem = Model()
    Node = NodeAttribute()
    Element = ElementAttribute()

    if not os.path.exists(input_name):
        ErrorMessage("Check the name of *.inp")
        #assert os.path.exists(input_name), "Check the name of *.inp"

    file1 = open(input_name,'r')
    line = file1.readline().strip()

    if (line == "*Title"):
        line = file1.readline().strip()
        Fem.title = line
        Fem.LogGen('./log/' +Fem.title + '.dat')
        Fem.LogWrite('Input name: ' + input_name + '\n\n')
        Fem.LogWrite('#Readinput')
        PrintCommand('*Title',1, Fem)
        PrintCommand(Fem.title,2, Fem)
    else:
        print("First content in *.inp must be *Title")

    line = file1.readline()
    line = file1.readline().strip()
    while line:
        if (line == "*ResultDirectory"):
            PrintCommand(line,1, Fem)
            line = file1.readline().strip()
            Fem.result = line
            PrintCommand(Fem.result, 2, Fem)
            if not os.path.exists(Fem.result):
                ErrorMessage("Check the result directory - *ResultDirectory in *.inp")
        
        elif (line == "*Mesh"):
            PrintCommand(line,1, Fem)
            line = file1.readline().strip()
            NodeName = line + ".NODE"
            if not os.path.exists(NodeName):
                ErrorMessage("Check the mesh name *Mesh in *.inp")
            PrintCommand(NodeName,2, Fem)
            fileNode = open(NodeName,'r')
            linetmp = fileNode.readline().strip()
            Node.NNode = int(linetmp)
            for ind in range(Node.NNode):    
                linetmp = fileNode.readline().strip()
                tmp = linetmp.replace(',', '\t')
                tmp = tmp.split('\t')
                Node.Id.append(int(tmp[0]))
                Node.Coord.append([float(tmp[1]), float(tmp[2])])  
            Fem.Dimension = len(Node.Coord[0])
            
            ElemName = line + ".ELEM"
            if not os.path.exists(ElemName):
                ErrorMessage("Check the mesh name *Mesh in *.inp")
            PrintCommand(ElemName,2, Fem)
            fileNode = open(ElemName,'r')
            linetmp = fileNode.readline().strip()
            Element.NElem = int(linetmp)
            for ind in range(Element.NElem):    
                linetmp = fileNode.readline().strip()
                tmp = linetmp.replace(',', '\t')
                tmp = tmp.split('\t')
                Element.Id.append(int(tmp.pop(0)))
                tmp = list(map(int, tmp))
                Element.Connectivity.append(tmp)
                
            #print(Element.Id)
            #print(Element.Connectivity)

        elif (line == "*ConstitutiveLaw"):
            PrintCommand(line,1, Fem)
            ReadConstitutiveLaw(file1, Fem)

        elif (line == "*LoadingStep"):
            PrintCommand(line,1, Fem)
            line = file1.readline().strip()
            Fem.totalstep = int(line)
            PrintCommand(str(Fem.totalstep),2, Fem)
            #exit(1)
            
        elif (line == "*Plane"):
            PrintCommand(line,1, Fem)
            line = file1.readline().strip()
            line = line.replace(" ", "")
            line = line.replace(',', '\t')
            tmp = line.split('\t')
            Fem.width = float(tmp[1])
            ReadPlane(tmp[0], Fem)
            
        elif (line == "*TimeIntegration"):
            PrintCommand(line,1, Fem)
            line = file1.readline().strip()
            Fem.timeintegration = line
            line = file1.readline().strip()
            line = line.replace(" ", "")
            line = line.replace(',', '\t')
            tmp = line.split('\t')
            Fem.totalstep  = int(tmp[0])
            Fem.totaltime  = float(tmp[1])
            #Fem.dt         = Constant(Fem.totaltime/Fem.totalstep)     # question what is the
            Fem.dt         = (Fem.totaltime/float(Fem.totalstep))     # question what is the Constant function??
            PrintCommand("Total step : " + str(Fem.totalstep ),2, Fem)
            PrintCommand("Total time : " + str(Fem.totaltime),2, Fem)
            PrintCommand("dt         : " + str(Fem.dt),2, Fem)
            
                
        else:
            print(line)
            print("Check command *XXX in *.inp!")
            exit(1)

        line = file1.readline().strip()
        count = 0
        while True:
            if(line == ""):
                line = file1.readline().strip()
                count +=1
            else:
                break
            if count == 5:
                break


    Node.BC_E    = np.empty([Node.NNode*Fem.Dimension])
    Node.BC_E[:] = np.NaN
    Node.u     = np.zeros([Node.NNode*Fem.Dimension])
    Node.u1    = np.zeros([Node.NNode*Fem.Dimension])
    Node.du    = np.zeros([Node.NNode*Fem.Dimension])
    Node.u_e   = np.zeros([Node.NNode*Fem.Dimension])
    Node.u_p   = np.zeros([Node.NNode*Fem.Dimension])
    Node.BC_N  = np.zeros([Node.NNode*Fem.Dimension])
    Node.F_int = np.zeros([Node.NNode*Fem.Dimension])
    Node.F_ext = np.zeros([Node.NNode*Fem.Dimension])
    Node.stress  = np.zeros([Node.NNode*Fem.Dimension,3])
    Node.sigma1  = np.zeros([Node.NNode*Fem.Dimension])
    Node.sigmaVM = np.zeros([Node.NNode*Fem.Dimension])

    Fem.show()

    file1.close
    #Fem.LogClose
    return Node, Element

class NodeAttribute:
    Coord = []
    u: np.array
    u1: np.array
    du: np.array
    u_e: np.array
    u_p: np.array
    BC_E: np.array
    BC_N: np.array
    F_int: np.array # internal force
    F_ext: np.array # external force
    stress: np.array # stress
    sigma1: np.array # principal stress
    # https://www.continuummechanics.org/principalstressesandstrains.html
    sigmaVM:np.array # Von Mises stress
    # https://en.wikipedia.org/wiki/Von_Mises_yield_criterion
    NNode: int
    Id = []
    
class ElementAttribute:
    Connectivity = []
    Stiffness = []
    B_matrix  = []
    NElem : int
    Id = []         # element id
    Area = []
    GPstrain_e = [] # elastic strain
    GPstrain_p = [] # plastic strain
    GPstrain = []
    GPstress = []
    P = []
    Jacc = []
    
class Model:
    title: str = ''
    result: str = ''
    analysis: str = ''
    twoD: str = ''
    dt: float = 0.
    totalstep: int = 1
    step: int = 0
    totaltime: float = 0.
    Dimension: int = 2
    EBC = []
    GP = [[-1./np.sqrt(3), -1./np.sqrt(3), 1.], 
          [ 1./np.sqrt(3), -1./np.sqrt(3), 1.],
          [ 1./np.sqrt(3),  1./np.sqrt(3), 1.],
          [-1./np.sqrt(3),  1./np.sqrt(3), 1.]]
    width: float = 0.
    
    E: float             # Elastic modulus 
    nu: float            # poisson ratio
    lamda: float         # 1st lame parameter
    mu: float            # 2nd lame parameter
    '''
    EBC[0] : Type of boundary condition
        U1 : Displacement x direction
        u2 : Displacement y direction
    EBC[1] : value of the condition (if it is fixed condition, simply 0.0)
    EBC[2:end] : condition ex) left, x[0]>0, x[1]>2.0
    '''
    NBC = []
    '''
    NBC[0] : Type of boundary condition
        Trac1  Trac2
        Force1 Force2
    NBC[1] : value of the condition
    NBC[2:end] : condition ex) left, x[0]>0, x[1]>2.0
    '''
    
    def LogGen(self, logname):
        self.file_log = open(logname,'w')
    
    def LogWrite(self, command):
        self.file_log.write(command)
        
    def LogClose(self):
        self.file_log.close
        
    def show(self):
        print("============= Model description =============")
        print("Model")
        print("\ttitle   : ",self.title)
        print("\tresult  : ",self.result,"\n")
        print("Solid")
        print("\tE    = ",self.E,"\tElastic Modulus")
        print("\tv    = ",self.nu,"\t\tPoisson ratio")
        print("\tlamda= ",self.lamda,"\t\t1st lame parameter")
        print("\tmu   = ",self.mu,"\t\t2st lame parameter \n")
        print("============= Model description =============")

def PrintCommand(command, Nindent, Fem):
    command.strip('\n')
    for i in range(0,Nindent):
        command = '\t' + command
    print(command)
    Fem.LogWrite(command+'\n')
    return 

def ErrorMessage(command):
    command ="\t Error : " + command
    print(command)
    assert 0,command
    exit(1)
    return

def ReadConstitutiveLaw(file1, Fem):
    line = file1.readline().strip()
    PrintCommand(line, 2, Fem)
    while line:
        if(line.lower() == 'solid'):
            line = file1.readline().strip()
            PrintCommand(line, 2, Fem)
            if(line.lower() == 'linearelasticity'):
                line = file1.readline().strip()
                line = line.replace(',', '\t')
                tmp = line.split('\t')
                Fem.E  = float(tmp[0])
                Fem.nu = float(tmp[1])
                PrintCommand("Elastic modulus: " + str(Fem.E),3, Fem)
                PrintCommand("Poisson ratio  : " + str(Fem.nu),3, Fem)
                line = file1.readline().strip()
            elif(line.lower() == 'elastoplasticity'):
                line = file1.readline().strip()
                line = line.replace(',', '\t')
                tmp = line.split('\t')
                Fem.E  = float(tmp[0])
                Fem.nu = float(tmp[1])
                Fem.MatProp  = []
                for ii in tmp:
                    Fem.MatProp.append(float(ii))
                for ii, param in enumerate(Fem.MatProp):
                    PrintCommand(str(ii+1)+"nd parameter : " + str(param),3, Fem)
                line = file1.readline().strip()
                PrintCommand("(1st and 2nd parameter are assumed as E and v)",3, Fem)
            else:
                print("Check Solid at *ConstitutiveLaw in *.inp") 
                print("\tLinearElasticity is available")
                print("\t\t input: E, v (Elastic modulus, Poisson ratio)")
                exit(1)
        else:
            print(line)
            print("Check *ConstitutiveLaw in *.inp") 
            print("\tAnalysisType Solid-> Solid")
            exit(1)
    return

def ReadPlane(line, Fem):
    E = Fem.E
    nu = Fem.nu
    if(line.lower() == 'planestrain'):
        Fem.twoD  = 'planestrain'
        Fem.lamda = E*nu / ((1.0+nu)*(1.0-2.0*nu))
        Fem.K     = E / (3*(1.0-2.0*nu))
        Fem.mu    = E / (2.0*(1.0+nu))
    elif(line.lower() == 'planestress'):
        Fem.twoD  = 'planestress'
        Fem.lamda = E*nu / ((1.0+nu)*(1.0-1.0*nu))
        Fem.K     = E / (2*(1.0-1.0*nu))
        Fem.mu    = E / (2.0*(1.0+nu))
    else:
        print("Check Plane in *inp")
        print("\tPlaneStrain, PlaneStress are avaiable")
        exit(1)
    PrintCommand(line,2, Fem)
