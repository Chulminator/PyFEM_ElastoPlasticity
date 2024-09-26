# ---------------------------------------------------------------- 
# 
# Written by: CK, HS in promechanics.org
# ----------------------------------------------------------------       


import numpy as np
import sys
import time as ct
import os.path
from dataclasses import dataclass
from dolfin import *


#def readinput(input_name):


class Fem:
    title: str = ''
    result: str = ''
    analysis: str = ''
    solver: str = ''
    thick: float
    dt: float
    totalstep: int
    totaltime: float 
    company_feed: bool = False
    
    def LogGen(self, logname):
        self.file_log = open(logname,'w')
    
    def LogWrite(self, command):
        self.file_log.write(command)
        
    def LogClose(self):
        self.file_log.close
        
    def show(self):
        print("title   : ",self.title)
        print("result  : ",self.result)
        print("analysis: ",self.analysis)
        print("solver  : ",self.solver)
        print("thickness =  ",self.thick)
        print("totalstep =  ",self.totalstep)
        print("totaltime =  ",self.totaltime)
        print("dt        =  ",self.dt)
        
    class PhaseField:
        Gc: float
        ls: float
        def show(self):
            print("Gc = ",self.Gc,"\tFracture energy")
            print("ls = ",self.ls,"\tlength scale parameter")
    
    class Solid:
        E: float            # Elastic modulus 
        v: float            # poisson ratio
        def show(self):
            print("E = ",self.E,"\tElastic Modulus")
            print("v = ",self.v,"\t\tPoisson ratio")
        
    class Hydro:
        k_w: float          # Intrinsic permeability
        mu_w: float         # Dynamic viscosity of water
        poro: float        # Porosity
        k_b_w: float        # Bulk modulus of water 
        def show(self):
            print("k_w   = ",self.k_w,"\tIntrinsic permeability")
            print("mu_w  = ",self.mu_w,"\tDynamic viscosity of fluid")
            print("poro  = ",self.poro,"\tPorosity")
            print("k_b_w = ",self.k_b_w,"\tBulk modulus of fluid")
        
    class Termo:
        E: float            # Elastic modulus 
        thick: float        # Thickness 
        v: float            # poisson ratio
        
        def show(self):
            print("123123")


def PrintCommand(command, Nindent, fem):
    command.strip('\n')
    for i in range(0,Nindent):
        command = '\t' + command
    print(command)
    fem.LogWrite(command+'\n')
    return 

def ErrorMassage(command):
    command ="\t Error : " + command
    print(command)
    return

def ReadConstitutiveLaw(line, fem):
#*ConstitutiveLaw
#LinearElastic
    if(line != 'Monolithic' and line != 'Staggered'):
        print("Available solvers are Monolithic and Staggered")  #################### edit
        exit(1)
    fem.solver = line
    PrintCommand(line, 2, fem)
    return

def ReadAnaylsisType(line, fem):
    if(line == 'Solid'):  ########### edit!
        print(1)
    elif (line == 'THM'):
        print(2)
    elif (line == 'UP'):
        print(3)
    else:
        print("Current vailable analysis is following: ")
        print("\tSolid: fem.Solid mechanics edit!")        ########### edit!
        print("\tUP: fem.Hydro-XXX couple formulation more description edit!")
        print("\tTHM: Thermo-hydro-XXX formulation more description edit!")
        exit(1)
    PrintCommand(line,2, fem)
    fem.analysis = line
    
    
def ReadPlane(line, fem):
    E = fem.Solid.E
    nu = fem.Solid.v
    if(line == 'PlaneStrain'):
        lamda = E*nu / ((1.0+nu)*(1.0-2.0*nu))
        mu    = E / (2.0*(1.0+nu))
    elif(line == 'PlaneStress'):
        alpha = (1.-2.*nu)/(1.-nu)
        lamda = E*nu / ((1.0+nu)*(1.0-2.0*nu))* alpha
        mu    = E / (2.0*(1.0+nu))
    else:
        print("Check Plane in *inp")
        print("\tPlaneStrain, PlaneStress are avaiable")
        exit(1)
    PrintCommand(line,2, fem)
    



input_name = "./input/terzaghi/NonpolarElasticity.inp"
print('#readinput')
fem = Fem()


if not os.path.exists(input_name):
    ErrorMassage("Check the name of *.inp")
    exit(1)
    #assert os.path.exists(input_name), "Check the name of *.inp"

file1 = open(input_name,'r')
line = file1.readline().strip()

if (line == "*Title"):
    line = file1.readline().strip()
    fem.title = line
    fem.LogGen('./log/' +fem.title + '.dat')
    fem.LogWrite('Input name: ' + input_name + '\n\n')
    fem.LogWrite('#Readinput')
    PrintCommand('*Title',1, fem)
    PrintCommand(fem.title,2, fem)
else:
    print("First content in *.inp must be *Title")
#fem.LogGen('./log/' +fem.title + '.dat')

line = file1.readline()
line = file1.readline().strip()
if (line == "*AnaylsisType"):
    PrintCommand(line,1, fem)
    line = file1.readline().strip()
    ReadAnaylsisType(line, fem)
else:
    print("Second content in *.inp must be *AnaylsisType")
    

line = file1.readline()
line = file1.readline().strip()
while line:
    if (line == "*ResultDirectory"):
        PrintCommand(line,1, fem)
        line = file1.readline().strip()
        fem.result = line
        PrintCommand(fem.result, 2, fem)
        if not os.path.exists(fem.result):
            ErrorMassage("Check the result directory - *ResultDirectory in *.inp")
            exit(1)
    
    elif (line == "*Mesh"):
        PrintCommand(line,1, fem)
        line = file1.readline().strip()
        if(line == 'Manual'):
            print("Manual function for mesh is not ready yet")  #################### edit
            exit(1)
        PrintCommand(line, 2, fem)
        if not os.path.exists(line):
            ErrorMassage("Check the .xml file - *Mesh in *.inp")
            exit(1)
        mesh = Mesh(line)
        
    elif (line == "*Solid"):
        PrintCommand(line,1, fem)
        line = file1.readline().strip()
        line = line.replace(',', '\t')
        tmp = line.split('\t')
        fem.Solid.E = float(tmp[0])
        fem.Solid.v = float(tmp[1])
        fem.thick = float(tmp[2])
        PrintCommand("Elastic modulus: " + str(fem.Solid.E),2, fem)
        PrintCommand("Poisson ratio  : " + str(fem.Solid.v),2, fem)
        PrintCommand("Thickness      : " + str(fem.thick),2, fem)
        
    elif (line == "*ConstitutiveLaw"):
        PrintCommand(line,1, fem)
        line = file1.readline().strip()
        ReadConstitutiveLaw(line, fem)
        
    elif (line == "*Solver"):
        PrintCommand(line,1, fem)
        line = file1.readline().strip()
        if(line != 'Monolithic' and line != 'Staggered'):
            print("Available solvers are Monolithic and Staggered")  #################### edit
            exit(1)
        fem.solver = line
        PrintCommand(line, 2, fem)
        
    elif (line == "*Hydro"):
        PrintCommand(line,1, fem)
        line = file1.readline().strip()
        line = line.replace(" ", "")
        line = line.replace(',', '\t')
        tmp = line.split('\t')
        fem.Hydro.k_w   = float(tmp[0])
        fem.Hydro.mu_w  = float(tmp[1])
        fem.Hydro.poro  = float(tmp[2])
        fem.Hydro.K_b_w = float(tmp[3])
        PrintCommand("Intrinsic permeability     : " + str(fem.Hydro.k_w),2, fem)
        PrintCommand("Dynamic viscosity of water : " + str(fem.Hydro.mu_w),2, fem)
        PrintCommand("Porosity                   : " + str(fem.Hydro.poro),2, fem)
        PrintCommand("Bulk modulus of fluid      : " + str(fem.Hydro.K_b_w),2, fem)
        
    elif (line == "*Plane"):
        PrintCommand(line,1, fem)
        line = file1.readline().strip()
        ReadPlane(line, fem)
        
    elif (line == "*TimeIntegration"):
        PrintCommand(line,1, fem)
        line = file1.readline().strip()
        fem.timeintegration = line
        line = file1.readline().strip()
        line = line.replace(" ", "")
        line = line.replace(',', '\t')
        tmp = line.split('\t')
        fem.totalstep  = int(tmp[0])
        fem.totaltime  = float(tmp[1])
        #fem.dt         = Constant(fem.totaltime/fem.totalstep)     # question what is the
        fem.dt         = (fem.totaltime/float(fem.totalstep))     # question what is the Constant function??
        PrintCommand("Total step : " + str(fem.totalstep ),2, fem)
        PrintCommand("Total time : " + str(fem.totaltime),2, fem)
        PrintCommand("dt         : " + str(fem.dt),2, fem)
        
            
    else:
        print("Check *XXX in *.inp!")
        exit(1)

    line = file1.readline().strip()
    if(line == ""):
        line = file1.readline().strip()

fem.show()
fem.Solid.show(fem.Solid)



file1.close
#fem.LogClose
