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
    thick: float
    dt: float
    TotalStep: int
    TotalTime: float    
    company_feed: bool = False
    
    def LogGen(self, logname):
        self.file_log = open(logname,'w')
    
    def LogWrite(self, command):
        self.file_log.write(command)
        
    def LogClose(self):
        self.file_log.close
    
class Elasticity:
    E: float            # Elastic modulus 
    thick: float        # Thickness 
    v: float            # poisson ratio
    
class Hydro:
    k_w: float          # Intri 
    mu_w: float         # Dynamic viscosity of water
    poro: float        # Porosity
    k_b_w: float        # Bulk modulus of water 
    
class Termo:
    E: float            # Elastic modulus 
    thick: float        # Thickness 
    v: float            # poisson ratio


def PrintCommand(command, Nindent):
    command.strip('\n')
    for i in range(0,Nindent):
        command = '\t' + command
    print(command)
    return 

def ErrorMassage(command):
    command ="\t Error : " + command
    print(command)
    return
    



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
    PrintCommand(line,1)
    line = file1.readline().strip()
    fem.title = line
    PrintCommand(fem.title,2)
else:
    print("First content in *.inp must be *Title")

print(type(file1))
exit(1)


fem.LogGen('./log/' +fem.title + '.dat',w)

line = file1.readline()
line = file1.readline().strip()
if (line == "*AnaylsisType"):
    PrintCommand(line,1)
    line = file1.readline().strip()
    if(line == 'Solid'):  ########### edit!
        solid = Elasticity()
    elif (line == 'THM'):
        solid = Elasticity()
        thermo = Thermo()
        hydro = Hydro()
    elif (line == 'UP'):
        solid = Elasticity()
        hydro = Hydro()
    else:
        print("Current vailable analysis is following: ")
        print("\tSolid: solid mechanics edit!")        ########### edit!
        print("\tUP: Hydro-XXX couple formulation more description edit!")
        print("\tTHM: Thermo-hydro-XXX formulation more description edit!")
        exit(1)
    PrintCommand(line,2)
    Fem.analysis = line
else:
    print("Second content in *.inp must be *AnaylsisType")
    

line = file1.readline()
line = file1.readline().strip()
while line:
    if (line == "*ResultDirectory"):
        PrintCommand(line,1)
        line = file1.readline().strip()
        fem.result = line
        PrintCommand(fem.result, 2)
        if not os.path.exists(fem.result):
            ErrorMassage("Check the result directory - *ResultDirectory in *.inp")
            exit(1)
    
    elif (line == "*Mesh"):
        PrintCommand(line,1)
        line = file1.readline().strip()
        PrintCommand(line, 2)
        if not os.path.exists(line):
            ErrorMassage("Check the .xml file - *Mesh in *.inp")
            exit(1)
        mesh = Mesh(line)
        
    elif (line == "*Elastic"):
        PrintCommand(line,1)
        line = file1.readline().strip()
        line = line.replace(',', '\t')
        tmp = line.split('\t')
        solid.E = float(tmp[0])
        solid.v = float(tmp[1])
        fem.thick = float(tmp[2])
        PrintCommand("Elastic modulus: " + str(solid.E),2)
        PrintCommand("Poisson ratio  : " + str(solid.v),2)
        PrintCommand("Thickness      : " + str(fem.thick),2)
        
    elif (line == "*Hydro"):
        PrintCommand(line,1)
        line = file1.readline().strip()
        line = line.replace(" ", "")
        line = line.replace(',', '\t')
        tmp = line.split('\t')
        Hydro.k_w   = float(tmp[0])
        Hydro.mu_w  = float(tmp[1])
        Hydro.poro  = float(tmp[2])
        Hydro.K_b_w = float(tmp[3])
        PrintCommand("Intrinsic permeability     : " + str(Hydro.k_w),2)
        PrintCommand("Dynamic viscosity of water : " + str(Hydro.mu_w),2)
        PrintCommand("Porosity                   : " + str(Hydro.poro),2)
        PrintCommand("Bulk modulus of fluid      : " + str(Hydro.K_b_w),2)
        
    elif (line == "*Plane"):
        PrintCommand(line,1)
        line = file1.readline().strip()
        E = solid.E
        nu = solid.v
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
        PrintCommand(line,2)
        
    elif (line == "*TimeIntegration"):
        PrintCommand(line,1)
        line = file1.readline().strip()
        fem.timeintegration = line
        line = file1.readline().strip()
        line = line.replace(" ", "")
        line = line.replace(',', '\t')
        tmp = line.split('\t')
        fem.TotalStep  = int(tmp[0])
        fem.TotalTime  = float(tmp[1])
        #fem.dt         = Constant(fem.TotalTime/fem.TotalStep)     # question what is the
        fem.dt         = (fem.TotalTime/float(fem.TotalStep))     # question what is the Constant function??
        PrintCommand("Total step : " + str(fem.TotalStep ),2)
        PrintCommand("Total time : " + str(fem.TotalTime),2)
        PrintCommand("dt         : " + str(fem.dt),2)
        
            
    else:
        print("Check *XXX in *.inp!")
        exit(1)

    line = file1.readline().strip()
    if(line == ""):
        line = file1.readline().strip()
file1.close


#while 


#print(Lines)
#Lines = file1.readline()
#print(Lines)
    
#Line = file1.readline()
#if(Lines == Null):
    #print("zzz")
    

    
    
    
    



  

