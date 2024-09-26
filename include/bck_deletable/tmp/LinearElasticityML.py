# ---------------------------------------------------------------- 
# Written by: CK in promechanics.org at Columbia University
# ----------------------------------------------------------------       
import numpy as np
import sys
import time as ct
import os.path
import math
import torch
import torch.nn as nn
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt 
from sklearn.preprocessing import MinMaxScaler

def StrainInvarants(strain):
    #print(strain)
    devT = torch.zeros(3,3)
    ev = strain[0] + strain[1]
    
    devT[0,0] = strain[0] - ev/3.
    devT[0,1] = strain[2]/2
    devT[1,0] = strain[2]/2
    devT[1,1] = strain[1] - ev/3.
    devT[2,2] = - ev/3.
    
    det = torch.sqrt(torch.trace(devT.T @ devT))
    n = devT / det
    es = det * np.sqrt(2./3.)
    
    return ev, es, n

def InvarantsToStress(p,q,n):
    stress_voigt = torch.zeros(3)
    tmp2 = np.sqrt(2./3.)*q*n
    stress = torch.eye(3)*p + tmp2
    stress_voigt[0] = stress[0,0]
    stress_voigt[1] = stress[1,1]
    stress_voigt[2] = stress[0,1]
    return stress_voigt.double()

def ElementalNodeDof(Node, Edof): 
    tmp1 = np.array(Edof)*2 - 2
    tmp2 = np.array(Edof)*2 - 1
    
    c = np.empty((tmp1.size + tmp1.size), dtype =tmp1.dtype)
    c[0::2] = tmp1
    c[1::2] = tmp2
    return c

def AssembleageForce(Fem, Node, Element):
    Global_F = torch.zeros(Node.NNode * Fem.Dimension)
    #Global_F.requires_grad = False
    #Global_F.requires_grad = True
    for EID, Edof in enumerate(Element.Connectivity):
        ENdof = ElementalNodeDof(Node, Edof)
        #print(Element.P[EID].squeeze())
        #print(Element.P[EID])
        #print(Element.P[EID].shape)
        #print(Global_F[ENdof])
        #print(Global_F[ENdof].shape)
        #input(11212)
        Global_F[ENdof] += Element.P[EID]

        
    return Global_F
    

def NormalizeInput(ev,es, model):
    normstrain = torch.zeros(2)
    scale = model.scales_inp.astype(np.float64)
    scale = torch.tensor(scale)
    limit = model.limits_inp.astype(np.float64)
    limit = torch.tensor(limit)
    normstrain[0] = ev * scale[0] + limit[0]
    normstrain[1] = es * scale[1] + limit[1]
    return normstrain.double()

def UnNormalizeInput(grad_pred, model):
    scale = model.scales_grad.astype(np.float64)
    scale = torch.tensor(scale)
    limit = model.limits_grad.astype(np.float64)
    limit = torch.tensor(limit)
    p1 = (grad_pred[0] -limit[0])/scale[0]
    q1 = (grad_pred[1] -limit[1])/scale[1]
    return p1, q1
    
class NNModel2(nn.Module):
    def __init__(self, Ndof):
        super(NNModel2,self).__init__()
        self.weight1 = nn.Parameter(torch.randn(Ndof)/10.)
        #self.weight1 = nn.Parameter(torch.ones(Ndof)/10)
        
    def forward(self, Fem, Node, Element, model):  
        IndexBCN = np.where(np.isnan(Node.BC_E))[0]
        Node.u_torch1 = torch.zeros(Node.NNode*2).double()
        Node.u_torch1[IndexBCN] = self.weight1 * Node.u_torch[IndexBCN]
        
        G_Edof = Fem.Dimension* len(Element.Connectivity[0])
        Element.P = []
        for ind1, Edof in enumerate(Element.Connectivity):
            ENdof = ElementalNodeDof(Node,Edof)
            E_P = torch.zeros(G_Edof)
            B = Element.B_matrix[ind1]
            JaccElem = Element.Jacc[ind1]
            for ind, gp in enumerate(Fem.GP):
                B_GP = B[ind]
                Jacc = JaccElem[ind]
                strain = torch.from_numpy(B_GP) @ Node.u_torch1[ENdof]
                ev,es,n = StrainInvarants(strain)
                #print("ev\n",ev)
                #print("es\n",es)
                #print("n\n",n)
                normstrain = NormalizeInput(ev,es, model)
                model.double()
                #print(normstrain)
                grad_pred = model.forward(normstrain)
                p, q = UnNormalizeInput(grad_pred, model)
                #print("\n\np\n",p)
                #print("\n\nq\n",q)
                stress = InvarantsToStress(p, q, n)
                #print("\n\nstress\n",stress)
                GP_P = torch.tensor(B_GP.T) @ stress * torch.tensor(Fem.width * Jacc * gp[-1])
                E_P += GP_P.T
                #input(111)
            #print(E_P)
            #input(11123123)
            Element.P.append(E_P)
        Global_F = AssembleageForce(Fem, Node, Element)
        return Global_F
    

class NNModel3(nn.Module):
    def __init__(self, Ndof):
        super(NNModel3,self).__init__()
        self.weight1 = nn.Parameter(torch.randn(Ndof))
        
    def forward(self, Fem, Node, Element):
        IndexBCN = np.where(np.isnan(Node.BC_E))[0]
        Node.u_torch1 = torch.zeros(Node.NNode*2).double()
        Node.u_torch1[IndexBCN] = self.weight1 * Node.u_torch[IndexBCN]
        
        G_Edof = Fem.Dimension* len(Element.Connectivity[0])
        Dim    = Fem.Dimension
        Element.P = []
        D = ConstitutiveLaw(Fem)
        Global_F     = torch.zeros(Node.NNode * Fem.Dimension)
        
        for ind1, Edof in enumerate(Element.Connectivity):
            ENdof = ElementalNodeDof(Node,Edof)
            GP_F = torch.tensor(Element.Stiffness[ind1]) @ Node.u_torch1[ENdof]
            Global_F[ENdof] += GP_F
        return Global_F

def ShapeFunction(s,t):
    N_matrix = np.array([(s-1.)*(t-1.)/4.,\
                        (-s-1.)*(t-1.)/4.,\
                        (-s-1.)*(-t-1.)/4.,\
                        (s-1.)*(-t-1.)/4.])

    dN   = np.array([[t-1., 1.-t, 1.+t, -1.-t],
                    [s-1., -1.-s, 1.+s, 1.-s]])
    dN *= 0.25
    return N_matrix, dN 

def ConstitutiveLaw(Fem):
    tmp1 = np.zeros([3,3])
    tmp1[0,0] = 1.
    tmp1[0,1] = 1.
    tmp1[1,0] = 1.
    tmp1[1,1] = 1.
    
    tmp2 = np.zeros([3,3])
    tmp2[0,0] = 2.
    tmp2[1,1] = 2.
    tmp2[2,2] = 1.
    return tmp1 *Fem.lamda + tmp2 * Fem.mu
