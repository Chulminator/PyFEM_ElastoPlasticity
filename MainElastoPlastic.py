# ---------------------------------------------------------------- 
# Written by: Chulmin Kweon  (Chulmin.Kweon@columbia.edu)
# ----------------------------------------------------------------       

import numpy as np
import sys
sys.path.append("./include")
import time as ct

from ReadAnalysis import *
import WriteResult
from  ElastoPlasticity import *

import FemBulk as Bulk
import importlib
import pandas as pd
from MeshHandler import *
from DataStructure import *


max_itr = 20 # Global max iteration
tic = ct.time()
assert len(sys.argv) == 2, "check the input value\n execution example:\n\t python3 MainElastoPlastic.py ./input/PlateWithHole/PlateWithHole_113_93_plastic.inp"


# Model Generation
input_name = sys.argv[1]
Program = ProgramAttrib()


# Read input
Model = ReadInput(input_name, Program)
#exit(1)


print ("FEM analysis is starting in 2 seconds")
ct.sleep(2)
WriteResult.WritePVD(Program, Model)


##########################################################################################
# Fem Nonlinear Analysis
ElementSetUp(Program, Model)


if Model.Dim == 2:
    WriteAttributes = np.zeros((Model.totalstep,10))
elif Model.Dim == 3:
    WriteAttributes = np.zeros((Model.totalstep,17))


if Model.BCFile == True:
    tmp = importlib.import_module(Model.BCFileDirectory)
    tmp.ApplyBC(Model)
else:
    assert False, "Other way to apply boundary is not implemented yet"

#tic1 = ct.time()
for step in range(Model.totalstep):
    #if step == 5:
        #print(ct.time() - tic1)
        #exit(1)
    Model.step = step+1
    Bulk.BC_SetUpAtStep(Model)
    Program.PrintCommand("step: "+str(Model.step),0)


    for ii in range(max_itr):
        Model.GlobalNRStep = ii
        ConstructStiffness(Model)
        Global_K = Bulk.Assembleage(Model)
        loss = Solve(Model, Global_K)
        Program.PrintCommand("\t"+str(loss),0)
        #if loss < 1e-9:
        if loss < 1e-10:
            break


    Bulk.CalculateStrainAtNode(Model)
    Bulk.CalculateStressAtNode(Model)
    Bulk.CalculatePlasticMultiplierAtNode(Model)
    Bulk.CalculateVonMisesStressAtNode(Model)
    Bulk.CalculateNodalDisp(Model)


    if Model.Dim == 2:
        WriteAttributes[step, 0] = np.mean(Model.NodalsigmaVM)
        WriteAttributes[step, 1] = np.mean(Model.Nodalstress[:,0])
        WriteAttributes[step, 2] = np.mean(Model.Nodalstress[:,1])
        WriteAttributes[step, 3] = np.mean(Model.Nodalstress[:,2])
        WriteAttributes[step, 4] = np.mean(Model.Nodalstrain_e[:,0]+ Model.Nodalstrain_p[:,0])
        WriteAttributes[step, 5] = np.mean(Model.Nodalstrain_e[:,1]+ Model.Nodalstrain_p[:,1])
        WriteAttributes[step, 6] = np.mean(Model.Nodalstrain_e[:,2]+ Model.Nodalstrain_p[:,2])
        WriteAttributes[step, 7] = np.mean(Model.NodalPlasticMultiplier)
        WriteAttributes[step, 8] = np.mean(Model.Nodalu[0::2])
        WriteAttributes[step, 9] = np.mean(Model.Nodalu[1::2])
    elif Model.Dim == 3:
        WriteAttributes[step, 0]  = np.mean(Model.NodalsigmaVM)
        WriteAttributes[step, 1]  = np.mean(Model.Nodalstress[:,0])
        WriteAttributes[step, 2]  = np.mean(Model.Nodalstress[:,1])
        WriteAttributes[step, 3]  = np.mean(Model.Nodalstress[:,2])
        WriteAttributes[step, 4]  = np.mean(Model.Nodalstress[:,3])
        WriteAttributes[step, 5]  = np.mean(Model.Nodalstress[:,4])
        WriteAttributes[step, 6]  = np.mean(Model.Nodalstress[:,5])
        WriteAttributes[step, 7]  = np.mean(Model.Nodalstrain_e[:,0]+ Model.Nodalstrain_p[:,0])
        WriteAttributes[step, 8]  = np.mean(Model.Nodalstrain_e[:,1]+ Model.Nodalstrain_p[:,1])
        WriteAttributes[step, 9]  = np.mean(Model.Nodalstrain_e[:,2]+ Model.Nodalstrain_p[:,2])
        WriteAttributes[step, 10] = np.mean(Model.Nodalstrain_e[:,3]+ Model.Nodalstrain_p[:,3])
        WriteAttributes[step, 11] = np.mean(Model.Nodalstrain_e[:,4]+ Model.Nodalstrain_p[:,4])
        WriteAttributes[step, 12] = np.mean(Model.Nodalstrain_e[:,5]+ Model.Nodalstrain_p[:,5])
        WriteAttributes[step, 13] = np.mean(Model.NodalPlasticMultiplier)
        WriteAttributes[step, 14]  = np.mean(Model.Nodalu[0::3])
        WriteAttributes[step, 15]  = np.mean(Model.Nodalu[1::3])
        WriteAttributes[step, 16]  = np.mean(Model.Nodalu[2::3])


    WriteResult.WriteVTU(Program, Model)

    #attribute = np.zeros((Node.NNode, 4))
    #attribute[:,0] = Node.sigmaVM
    #attribute[:,1] = Node.stress[:,0]
    #attribute[:,2] = Node.stress[:,1]
    #attribute[:,3] = Node.stress[:,2]
    #WriteResult.WriteCustomizedAttribute(Fem, Node, Element, attribute)

    ##########################################################################################
    ##### Sandia Inverse Problem Input
    ##### Post processing
    if Model.Dim == 3:
        DispX = Model.Nodalu[0::3]
        DispY = Model.Nodalu[1::3]
        DispZ = Model.Nodalu[2::3]
        Attribute = np.zeros((Model.NNode,3))
        Attribute[:,0] = DispX
        Attribute[:,1] = DispY
        Attribute[:,2] = DispZ
        np.save(Program.result+'/'+Program.title+'_disp_'+str(Model.step)+'.npy', Attribute)
    elif Model.Dim == 2:
        DispX = Model.Nodalu[0::2]
        DispY = Model.Nodalu[1::2]
        Attribute = np.zeros((Model.NNode,5))
        Attribute[:,0] = DispX
        Attribute[:,1] = DispY
        Attribute[:,2] = Model.Nodalstrain_e[:,0]+ Model.Nodalstrain_p[:,0]
        Attribute[:,3] = Model.Nodalstrain_e[:,1]+ Model.Nodalstrain_p[:,1]
        Attribute[:,4] = Model.Nodalstrain_e[:,2]+ Model.Nodalstrain_p[:,2]
        np.save(Program.result+'/'+Program.title+'_info_'+str(Model.step)+'.npy', Attribute)
        np.savetxt(Program.result+'/'+Program.title+'_info_'+str(Model.step)+'.txt', Attribute)
    ##########################################################################################


if Model.Dim == 2:
    output = pd.DataFrame({'VonMises': WriteAttributes[:, 0], \
                           's11'     : WriteAttributes[:, 1], \
                           's22'     : WriteAttributes[:, 2], \
                           's12'     : WriteAttributes[:, 3], \
                           'e11'     : WriteAttributes[:, 4], \
                           'e22'     : WriteAttributes[:, 5], \
                           'e12'     : WriteAttributes[:, 6], \
                           'lamda'   : WriteAttributes[:, 7], \
                           'ux'      : WriteAttributes[:, 8], \
                           'uy'      : WriteAttributes[:, 9]})
if Model.Dim == 3:
    output = pd.DataFrame({'VonMises': WriteAttributes[:, 0], \
                           's11'     : WriteAttributes[:, 1], \
                           's22'     : WriteAttributes[:, 2], \
                           's33'     : WriteAttributes[:, 3], \
                           's12'     : WriteAttributes[:, 4], \
                           's23'     : WriteAttributes[:, 5], \
                           's31'     : WriteAttributes[:, 6], \
                           'e11'     : WriteAttributes[:, 7], \
                           'e22'     : WriteAttributes[:, 8], \
                           'e33'     : WriteAttributes[:, 9], \
                           'e12'     : WriteAttributes[:, 10], \
                           'e23'     : WriteAttributes[:, 11], \
                           'e31'     : WriteAttributes[:, 12], \
                           'lamda'   : WriteAttributes[:, 13], \
                           'ux'      : WriteAttributes[:, 14], \
                           'uy'      : WriteAttributes[:, 15], \
                           'uz'      : WriteAttributes[:, 16]})

savename = Program.result+'/'+Program.title+'_plotfile.csv'
output.to_csv(savename, index = False)


# End of the program
toc = ct.time() - tic
Program.PrintCommand(' ',0)
Program.PrintCommand("Elapsed CPU time: "+str(toc)+"[sec]",0)
Program.PrintCommand("file log is save in ./log/" +Program.title + ".dat",0)
Program.LogClose()

