# ---------------------------------------------------------------- 
# Written by: CK, HS in promechanics.org at Columbia University
# ----------------------------------------------------------------       
import numpy as np
import sys
import time as ct
import os.path
from dataclasses import dataclass
from dolfin import *

# Things to upgrade
#   PlotSetup   : what variable will be plotted how often
#   OutputSetup : what variable will be written how often

# Things to check
# PhaseField
# ConstitutiveLaw - Thermo


def ReadInput(input_name):
    #input_name = "./input/terzaghi/NonpolarElasticity.inp"
    print('#readinput')
    Fem = Model()

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
    if (line == "*AnalysisType"):
        PrintCommand(line,1, Fem)
        line = file1.readline().strip()
        ReadAnalysisType(line, Fem)
    else:
        print("Second content in *.inp must be *AnalysisType")
        

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
            if(line == 'Manual'):
                print("Manual function for mesh is not ready yet")  #################### edit
                exit(1)
            PrintCommand(line, 2, Fem)
            if not os.path.exists(line):
                ErrorMessage("Check the .xml file - *Mesh in *.inp")
            mesh = Mesh(line)
            
        elif (line == "*Solid"):
            PrintCommand(line,1, Fem)
            line = file1.readline().strip()
            line = line.replace(',', '\t')
            tmp = line.split('\t')
            Fem.Solid.E = float(tmp[0])
            Fem.Solid.nu = float(tmp[1])
            E  = Fem.Solid.E
            nu = Fem.Solid.nu
            Fem.Solid.lamda = E*nu / ((1.0+nu)*(1.0-2.0*nu))
            Fem.Solid.mu    = E / (2.0*(1.0+nu))
            
            PrintCommand("Elastic modulus: " + str(Fem.Solid.E),2, Fem)
            PrintCommand("Poisson ratio  : " + str(Fem.Solid.nu),2, Fem)
            
        elif (line == "*ConstitutiveLaw"):
            PrintCommand(line,1, Fem)
            ReadConstitutiveLaw(file1, Fem)
            
        elif (line == "*Solver"):
            PrintCommand(line,1, Fem)
            line = file1.readline().strip()
            if(line != 'Monolithic' and line != 'Staggered'):
                print("Check *Solver in *.inp") 
                print("\tAvailable solvers are Monolithic and Staggered") 
                exit(1)
            Fem.solver = line
            PrintCommand(line, 2, Fem)
            
        elif (line == "*Hydro"):
            PrintCommand(line,1, Fem)
            line = file1.readline().strip()
            line = line.replace(" ", "")
            line = line.replace(',', '\t')
            tmp = line.split('\t')
            Fem.Hydro.k_f   = float(tmp[0])
            Fem.Hydro.mu_f  = float(tmp[1])
            Fem.Hydro.poro  = float(tmp[2])
            Fem.Hydro.K_b_f = float(tmp[3])
            PrintCommand("Intrinsic permeability     : " + str(Fem.Hydro.k_f),2, Fem)
            PrintCommand("Dynamic viscosity of water : " + str(Fem.Hydro.mu_f),2, Fem)
            PrintCommand("Porosity                   : " + str(Fem.Hydro.poro),2, Fem)
            PrintCommand("Bulk modulus of fluid      : " + str(Fem.Hydro.K_b_f),2, Fem)
            
        elif (line == "*PhaseField"):
            PrintCommand(line,1, Fem)
            Fem.analysisFlag[3] = 1              # you should check here
            line = file1.readline().strip()
            line = line.replace(" ", "")
            line = line.replace(',', '\t')
            tmp = line.split('\t')
            Fem.PhaseField.Gc = float(tmp[0])
            Fem.PhaseField.ls = float(tmp[1])
            PrintCommand("Fracture Energy        : " + str(Fem.PhaseField.Gc),2, Fem)
            PrintCommand("Length scale parameter : " + str(Fem.PhaseField.ls),2, Fem)
            
        elif (line == "*EssentialBC"):
            PrintCommand(line,1, Fem)
            ReadEssentialBC(file1, Fem)
            
        elif (line == "*NaturalBC"):
            PrintCommand(line,1, Fem)
            ReadNaturalBC(file1, Fem)


        elif (line == "*Plane"):
            PrintCommand(line,1, Fem)
            line = file1.readline().strip()
            ReadPlane(line, Fem)
            
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

    #Fem.show()
    #Fem.Solid.show(Fem.Solid)
    #Fem.Hydro.show(Fem.Hydro)
    #Fem.Hydro.show(Fem.Hydro)

    file1.close
    #Fem.LogClose
    return mesh, Fem
# HSS comment : universial Notation 

class Model:
    title: str = ''
    result: str = ''
    analysis: str = ''
    analysisFlag = [1, 0, 0, 0] # Solid, Hydro, Thermo, PhaseField
    solver: str = ''
    dt: float
    totalstep: int
    totaltime: float 
    EBC = []
    '''
    EBC[0] : Type of boundary condition
        U1 : Displacement x direction
        u2 : Displacement y direction
        T  : Temperature
        P  : Pore pressure
        D  : Damage (initial notch))
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
        print("Model")
        print("\ttitle   : ",self.title)
        print("\tresult  : ",self.result)
        print("\tanalysis: ",self.analysis)
        print("\tsolver  : ",self.solver)
        print("\ttotalstep =  ",self.totalstep)
        print("\ttotaltime =  ",self.totaltime)
        print("\tdt        =  ",self.dt,"\n")
        
    class PhaseField:
        Gc: float           # Fracture energy
        lc: float           # Length scale parameter
        def show(self):
            print("\tGc = ",self.Gc,"\tFracture energy")
            print("\tls = ",self.ls,"\tlength scale parameter\n")
    
    class Solid:
        E: float             # Elastic modulus 
        nu: float            # poisson ratio
        lamda: float         # 1st lame parameter
        mu: float            # 2nd lame parameter
        
        def show(self):
            print("Solid")
            print("\tE    = ",self.E,"\tElastic Modulus")
            print("\tv    = ",self.nu,"\t\tPoisson ratio")
            print("\tlamda= ",self.lamda,"\t\t1st lame parameter")
            print("\tmu   = ",self.mu,"\t\t2st lame parameter \n")
            
        def epsilon(self, u):
            return sym(grad(u))

        def epsilon_v(self, u):
            strain = epsilon(u)
            return tr(strain) 
            # effective stress
            
        def sigma_e(self, u):
            print(self.E)
            print(self.nu)
            print(self.mu)
            print(self.lamda)
            return self.lamda*tr(epsilon(u))*Identity(len(u)) + 2.0* self.mu*epsilon(u)

            # total stress
        def sigma(self, u, p):
            stress_e = sigma_e(u)
            stress   = stress_e - p*Identity(2)

            return stress            
        
    class Hydro:
        k_f: float          # Intrinsic permeability
        mu_f: float         # Dynamic viscosity of fluid
        K_b_f: float        # Bulk modulus of fluid 
        poro: float        # Porosity
        tmp: float
        
        def show(self):
            print("Hydro")
            print("\tk_f   = ",self.k_f,"\tIntrinsic permeability")
            print("\tmu_f  = ",self.mu_f,"\tDynamic viscosity of fluid")
            print("\tporo  = ",self.poro,"\tPorosity")
            print("\tK_b_f = ",self.K_b_f,"\tBulk modulus of fluid\n")
            

        def w_Darcy(self, p):
            w_vel = (-self.k_f/self.mu_f)*grad(p)

            return w_vel
            
            
        
    class Thermo:  # edit
        E: float            # Elastic modulus 
        v: float            # poisson ratio
        
        def show(self):
            print("123123")


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
                Fem.Solid.E  = float(tmp[0])
                Fem.Solid.nu = float(tmp[1])
                E  = Fem.Solid.E
                nu = Fem.Solid.nu
                Fem.Solid.lamda = E*nu / ((1.0+nu)*(1.0-2.0*nu))
                Fem.Solid.mu    = E / (2.0*(1.0+nu))
                
                
                PrintCommand("Elastic modulus: " + str(Fem.Solid.E),3, Fem)
                PrintCommand("Poisson ratio  : " + str(Fem.Solid.nu),3, Fem)
                line = file1.readline().strip()
            else:
                print("Check Solid at *ConstitutiveLaw in *.inp") 
                print("\tLinearElasticity is available")
                print("\t\t input: E, v (Elastic modulus, Poisson ratio)")
                exit(1)
        elif(line.lower() == 'hydro'):
            line = file1.readline().strip()
            PrintCommand(line, 2, Fem)
            if(line.lower() == 'darcylaw'):
                line = file1.readline().strip()
                line = line.replace(',', '\t')
                tmp = line.split('\t')
                Fem.Hydro.k_f   = float(tmp[0])
                Fem.Hydro.mu_f  = float(tmp[1])
                Fem.Hydro.poro  = float(tmp[2])
                Fem.Hydro.K_b_f = float(tmp[3])
                PrintCommand("Intrinsic permeability     : " + str(Fem.Hydro.k_f),3, Fem)
                PrintCommand("Dynamic viscosity of fluid : " + str(Fem.Hydro.mu_f),3, Fem)
                PrintCommand("Porosity                   : " + str(Fem.Hydro.poro),3, Fem)
                PrintCommand("Bulk modulus of fluid      : " + str(Fem.Hydro.K_b_f),3, Fem)
                line = file1.readline().strip()
            else:
                print("Check Hydro at *ConstitutiveLaw in *.inp") 
                print("\tDarycylaw is available")
                print("\t\tinput:")
                print("\t\tIntrinsic permeability")
                print("\t\tDynamic viscosity of fluid")
                print("\t\tPorosity")
                print("\t\tBulk modulus of fluid")
                exit(1)
        elif(line.lower() == 'thermo'):
            line = file1.readline().strip()
            PrintCommand(line, 3, Fem)
            if(line.lower() == 'fourierlaw'):
                print("It should be implemented more edit") # edit!!!
                exit(1);
                #line = file1.readline().strip()
                #line = line.replace(',', '\t')
                #tmp = line.split('\t')
                #Fem.Hydro.k_f   = float(tmp[0])
                #Fem.Hydro.mu_f  = float(tmp[1])
                #Fem.Hydro.poro  = float(tmp[2])
                #PrintCommand("Intrinsic permeability     : " + str(Fem.Thermo.k_f),2, Fem)
                #PrintCommand("Dynamic viscosity of fluid : " + str(Fem.Thermo.mu_f),2, Fem)
                #PrintCommand("Porosity                   : " + str(Fem.Thermo.poro),2, Fem)
                #PrintCommand("Bulk modulus of fluid      : " + str(Fem.Thermo.K_b_f),2, Fem)
                #line = file1.readline().strip()
            else:
                print("Check Thermo at *ConstitutiveLaw in *.inp") 
                print("\tFourierLaw is available")
                print("\t\tinput:")
                print("\t\tIntrinsic permeability")
                exit(1)
        else:
            print("Check *ConstitutiveLaw in *.inp") 
            print("\tAnalysisType Solid-> Solid")
            print("\tAnalysisType UP   -> Solid, Hydro")
            print("\tAnalysisType THM  -> Solid, Hydro, Thermo")
            print("\tAnalysisType TM   -> Solid, Thermo")
            exit(1)
    return

def ReadEssentialBC(file1, Fem):
    line = file1.readline().strip()
    PrintCommand("*type\t*value\t*position", 2, Fem)
    
    while line:
        #EBCType = line
        line = line.replace(',  ', '\t')
        line = line.replace(', ', '\t')
        line = line.replace(',', '\t')
        tmp = line.split('\t')
        tmp1 = []
        for i in tmp:
            tmp1.append(i.lower())
        Fem.EBC.append(tmp1)
        
        if len(tmp) != 3:
            PrintCommand("Check the number of inputs in EBC",2,Fem)
            ErrorMessage("It should be 3 (type, value, position)",2,Fem)
            exit(1)
        
        PrintCommand(tmp[0] + "\t" + tmp[1] + "\t" + tmp[2],2, Fem)
        if (Fem.EBC[-1][0].lower() not in ['u1', 'u2', 't', 'p', 'd'] ): 
            print("Check *EssentialBC in *.inp") 
            ErrorMessage("the condition is one of \
                \n\tU1 : Displacement x direction \
                \n\tU2 : Displacement y direction \
                \n\tT  : Temperature\
                \n\tP  : Pore pressure \
                \n\tD  : Damage (initial notch)")
        line = file1.readline().strip()
    return

def ReadNaturalBC(file1, Fem):
    line = file1.readline().strip()
    PrintCommand("*type\t*value\t*position",2,Fem)
    #
    while line:
        line = line.replace(',  ', '\t')
        line = line.replace(', ', '\t')
        line = line.replace(',', '\t')
        tmp = line.split('\t')
        tmp1 = []
        for i in tmp:
            tmp1.append(i.lower())
        Fem.NBC.append(tmp1)
        
        if len(tmp) < 3:
            PrintCommand("Check the number of inputs in EBC",2,Fem)
            PrintCommand("It should be 3 (type, value, position, (degree))",2,Fem)
            exit(1)
        
        PrintCommand(tmp[0] + "\t" + tmp[1] + "\t" + tmp[2],2, Fem)
        
        
        NBCType = line
        if( Fem.NBC[-1][0].lower() not in ['trac1','trac2','f1','f2'] ):
            print("Check *NaturalBC in *.inp") 
            ErrorMessage("Trac1, Trac2 : Traction \
                \n\tForce     : Force")
        line = file1.readline().strip()
    return


def ReadAnalysisType(line, Fem):
    if(line.lower() == 'solid'):
        Fem.analysisFlag = [1, 0, 0, 0]
    elif (line.lower() == 'thm'):
        Fem.analysisFlag = [1, 1, 1, 0]
    elif (line.lower() == 'up'):
        Fem.analysisFlag = [1, 1, 0, 0]
    elif (line.lower() == 'tm'):
        Fem.analysisFlag = [1, 0, 1, 0]
    else:
        print("Current vailable analysis is following: ")
        print("\tSolid: Solid mechanics")        ########### edit!
        print("\tUP   : Porous solid analysis")
        print("\tTHM  : Thermo-hydro-solid ")
        print("\tTM   : edit ")
        exit(1)
    PrintCommand(line,2, Fem)
    Fem.analysis = line
    
    
'''    
def ReadPlane(line, Fem):
    E = Fem.Solid.E
    nu = Fem.Solid.nu
    if(line.lower() == 'planestrain'):
        lamda = E*nu / ((1.0+nu)*(1.0-2.0*nu))
        mu    = E / (2.0*(1.0+nu))
    elif(line.lower() == 'planestress'):
        alpha = (1.-2.*nu)/(1.-nu)
        lamda = E*nu / ((1.0+nu)*(1.0-2.0*nu))* alpha
        mu    = E / (2.0*(1.0+nu))
    else:
        print("Check Plane in *inp")
        print("\tPlaneStrain, PlaneStress are avaiable")
        exit(1)
    PrintCommand(line,2, Fem)
'''

