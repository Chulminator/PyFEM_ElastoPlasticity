# ---------------------------------------------------------------- 
# Written by: CK in promechanics.org at Columbia University
# ----------------------------------------------------------------       
import numpy as np
import sys
import time as ct
import os.path
from dataclasses import dataclass
from  ElastoPlasticity import *
from MeshHandler import *
from DataStructure import *


# Things to upgrade
#   PlotSetup   : what variable will be plotted how often
#   OutputSetup : what variable will be written how often



def ReadInput(input_name, Program):
    ## The reason Program and model are different is
    ## in case for the multiscale analysis

    print('#readinput')
    #ModelNode    = []
    #ModelElement = []

    if not os.path.exists(input_name):
        ErrorMessage("Check the name of *.inp")

    file1 = open(input_name,'r')
    line = file1.readline().strip()


    if (line == "*Title"):
        line = file1.readline().strip()
        Program.title = line
        Program.LogGen('./log/' +Program.title + '.dat')
        Program.LogWrite('Input name: ' + input_name + '\n\n')
        Program.LogWrite('#Readinput')
        Program.PrintCommand('*Title',1)
        Program.PrintCommand(Program.title,2)
    else:
        print("First content in *.inp must be *Title")

    Model   = ModelAttrib(Program.title)


    line = file1.readline()
    line = file1.readline().strip()

    if (line == "*Mesh"):
        ######################## Node ########################
        Program.PrintCommand(line,1)
        line = file1.readline().strip()
        NodeName = line + ".NODE"
        if not os.path.exists(NodeName):
            print(NodeName)
            ErrorMessage("Check the mesh name *Mesh in *.inp")
        Program.PrintCommand(NodeName,2)
        fileNode = open(NodeName,'r')
        linetmp = fileNode.readline().strip()

        Model.NNode = int(linetmp) # Total Number of node
        for ind in range(Model.NNode):
            linetmp = fileNode.readline().strip()
            tmp = linetmp.replace(',', '\t')
            tmp = tmp.split('\t')

            Id    = int(tmp.pop(0))
            Coord = list(map(float, tmp))

            Node = NodeAttrib(Coord)
            Node.SetId(Id)
            Node.InitAttrib()
            Model.Node.append(Node)
        Model.Dim = Node.Dim
        tmp = "Total number of node: "+str(Model.NNode)
        Program.PrintCommand(tmp, 3)

        ########################  Node   ########################
        ######################## Element ########################
        ElemName = line + ".ELEM"
        if not os.path.exists(ElemName):
            assert False, "Check the mesh name *Mesh in *.inp"

        Program.PrintCommand(ElemName,2)
        fileNode = open(ElemName,'r')
        linetmp = fileNode.readline().strip()
        Model.NElem = int(linetmp)

        for ind in range(Model.NElem):
            linetmp = fileNode.readline().strip()
            tmp = linetmp.replace(',', '\t')
            tmp = tmp.split('\t')
            Id           = int(tmp.pop(0))
            Connectivity = list(map(int, tmp))
            ElemNode = []
            for ii in Connectivity:
                node = Model.GetNodeAtId(ii)
                ElemNode.append( node )
                node.NElem += 1

            Element  = ElementAttrib(ElemNode)
            Element.SetId(Id)
            Element.InitAttrib()
            Element.SetElementType(Model.Dim)
            Model.Element.append(Element)
        tmp = "Total number of element: "+str(Model.NElem)
        Program.PrintCommand(tmp, 3)

        ######################## Element ########################
    else:
        ErrorMessage("Mesh information should be at the second")


    line = file1.readline()
    line = file1.readline().strip()

    while line:
        if (line == "*ResultDirectory"):
            Program.PrintCommand(line,1)
            line = file1.readline().strip()
            Program.result = line
            Program.PrintCommand(Program.result, 2)
            if not os.path.exists(Program.result):
                ErrorMessage("Check the result directory - *ResultDirectory in *.inp")


        elif (line == "*ConstitutiveLaw"):
            Program.PrintCommand(line,1)
            ReadConstitutiveLaw(file1, Program, Model)

            for Element in Model.Element:
                Element.SetMatProp(Model.MatProp)


        elif (line == "*LoadingStep"):
            Program.PrintCommand(line,1)
            line = file1.readline().strip()
            Model.totalstep = int(line)
            Program.PrintCommand(str(Model.totalstep),2)

        elif (line == "*BCFile"):
            Program.PrintCommand(line,1)
            Model.BCFile = True
            line = file1.readline().strip()
            Program.PrintCommand(line,2)
            Model.BCFileDirectory = line

        elif (line == "*Plane"):
            Program.PrintCommand(line,1)
            line = file1.readline().strip()
            line = line.replace(" ", "")
            line = line.replace(',', '\t')
            tmp = line.split('\t')
            Model.width = float(tmp[1])
            ReadPlane(tmp[0], Model)
            Program.PrintCommand(tmp[0].lower(),2)
            for Element in Model.Element:
                ReadPlane(tmp[0], Element)

        elif (line == "*ReturnMappingMode"):
            Program.PrintCommand(line,1)
            line = file1.readline().strip()
            Program.PrintCommand(line,2)
            Model.ReturnMappingMode = line.lower()
            if not (Model.ReturnMappingMode == 'eigenspace' or Model.ReturnMappingMode == 'pqspace' or Model.ReturnMappingMode == 'tensorspace'):
                print("\n\nAvaiable return mapping algorithms are folowing:")
                print("\tPQSpace")
                print("\tEigenSpace")
                print("\tTensorSpace")
                assert False, "\t-Check *ReturnMappingMode in .inp file"


        elif (line == "*PlasticModel"):
            Program.PrintCommand(line,1)
            line = file1.readline().strip()
            Program.PrintCommand(line,2)
            Model.PlasticModel = line
            # add code that this path exists or not


        elif (line == "*InitialPressure"):
            Program.PrintCommand(line,1)
            line = file1.readline().strip()
            Model.InitialStress = 1
            Model.HSP = float(line)
            Program.PrintCommand(str(Model.HSP),2)

            
        elif (line == "*TimeIntegration"): # for dynamic analysis
            Program.PrintCommand(line,1)
            line = file1.readline().strip()
            Model.timeintegration = line
            line = file1.readline().strip()
            line = line.replace(" ", "")
            line = line.replace(',', '\t')
            tmp = line.split('\t')
            Model.totalstep  = int(tmp[0])
            Model.totaltime  = float(tmp[1])
            #Model.dt         = Constant(Model.totaltime/Model.totalstep)     # question what is the
            Model.dt         = (Model.totaltime/float(Model.totalstep))     # question what is the Constant function??
            Program.PrintCommand("Total step : " + str(Model.totalstep ),2)
            Program.PrintCommand("Total time : " + str(Model.totaltime),2)
            Program.PrintCommand("dt         : " + str(Model.dt),2)
                
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

    for Element in Model.Element:
        ConstitutiveTmp = ConstitutiveLaw(Model, Element)
        Element.ConstitutiveModel = ConstitutiveTmp

    ModelFacet    = GenerateFacet(Model)
    Model.NFacet  = len(ModelFacet)
    Model.Facet   = ModelFacet

    ModelEdge     = GenerateEdge(Model)
    Model.NEdge   = len(ModelEdge)
    Model.Edge    = ModelEdge


# Example
## How to Remove Node
#node = Model.GetNodeAtId(1)
#Model.RemNode( Model.GetNodeAtId(2) )
#input("*"*20)

## How to Remove Element
#Model.RemElem( Model.Element[1] )
#Model.RemElem( Model.GetElementAtId(1) )
#input("*"*20)


    #for ii in Model.Node:
        #print(ii.Id, ii.NElem)

    #print("="*30)
    #1. remove node
        #s

    #2. remove element


    #for ii in Model.Element[0].Connectivity:
        ##print( ii == node )
        #print( ii.Id )
    #input("")


    #for ii in Model.Node:
        #del ii

    #node = Model.GetNodeAtId(25)
    #node  = Model.Node.pop(1)
    #print(node.Id)
    #node = Model.Node[24]
    #node = Model.Node[0]
    #Model.Node[0] = Model.NodeInvalid
    #print(node.Id)
    #for Element in Model.Element:
        #for ii in Element.Connectivity:
            #if ii == node:
                #del ii
    #del node

    #for ii in Model.Element[0].Connectivity:
        ##print( ii == node )
        #print( ii.Id )
    ##print( Model.Element[0].Connectivity )
    ##input("")
    #print("="*30)

    #Model.RemNode(node)
    #del node
    #print(node.Id)
    #delattr(Model,'Node')
    #del Model.Node
    #del node
    #print(node.Id)
    #input("1213213")

    #for ii in Model.Node:
        #print(ii.Id)
        ##print(ii == node)
    #print("="*30)
    #input("")


    #Model.RemNode(

    #exit(1)
    ##############################################################################
    # Edge should be also considered "LATER"
    ##############################################################################



    ##############################################################################
    ##here I need to use callback whenever the Node Element are changed But not now
    #Node  = ModelNode.pop(0)
    #del Node
    #print(Node.Id)
    ##############################################################################
    #for ii in ModelFacet.AdjacNode:
        #print(ii[0].Id, ii[1].Id)
    #for ii in ModelFacet.AdjacElem:
        #print(ii[0].Id, ii[1].Id)
    #exit(1)

    file1.close
    return Model


def ErrorMessage(command):
    command ="\t Error : " + command
    print(command)
    assert False, command
    exit(1)
    return


def ReadConstitutiveLaw(file1, Program, Model):
    line = file1.readline().strip()
    Program.PrintCommand(line, 2)
    while line:
        if(line.lower() == 'solid'):
            line = file1.readline().strip()
            Program.PrintCommand(line, 2)
            if(line.lower() == 'linearelasticity'):
                line = file1.readline().strip()
                line = line.replace(',', '\t')
                tmp = line.split('\t')
                Model.E  = float(tmp[0])
                Model.nu = float(tmp[1])
                Program.PrintCommand("Elastic modulus: " + str(Model.E),3)
                Program.PrintCommand("Poisson ratio  : " + str(Model.nu),3)
                line = file1.readline().strip()
            elif(line.lower() == 'elastoplasticity'):
                line = file1.readline().strip()
                line = line.replace(',', '\t')
                tmp = line.split('\t')
                Model.E  = float(tmp[0])
                Model.nu = float(tmp[1])
                Model.MatProp  = []
                for ii in tmp:
                    Model.MatProp.append(float(ii))
                for ii, param in enumerate(Model.MatProp):
                    Program.PrintCommand(str(ii+1)+"nd parameter : " + str(param),3)
                line = file1.readline().strip()
                Program.PrintCommand("(1st and 2nd parameter are assumed as E and v)",3)

                E  = Model.E
                nu = Model.nu
                Model.lamda = E*nu / ((1.0+nu)*(1.0-2.0*nu))
                Model.K     = E / (3*(1.0-2.0*nu))
                Model.mu    = E / (2.0*(1.0+nu))
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


def ReadPlane(line, Model):
    E = Model.E
    nu = Model.nu
    if(line.lower() == 'planestrain'):
        Model.twoD  = 'planestrain'
        Model.lamda = E*nu / ((1.0+nu)*(1.0-2.0*nu))
        Model.K     = E / (3*(1.0-2.0*nu))
        Model.mu    = E / (2.0*(1.0+nu))
    elif(line.lower() == 'planestress'):
        Model.twoD  = 'planestress'
        Model.lamda = E*nu / ((1.0+nu)*(1.0-1.0*nu))
        Model.K     = E / (2*(1.0-1.0*nu))
        Model.mu    = E / (2.0*(1.0+nu))
    else:
        print("Check Plane in *inp")
        print("\tPlaneStrain, PlaneStress are avaiable")
        exit(1)

