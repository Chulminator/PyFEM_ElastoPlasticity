# ---------------------------------------------------------------- 
# Written by: CK, HS in promechanics.org at Columbia University
# ----------------------------------------------------------------       
import numpy as np
import sys
import time as ct
import os.path
from dataclasses import dataclass
from dolfin import *




def DefineSpace(mesh, fem)
#----Element type in FEniCS
#    https://fenicsproject.org/olddocs/dolfin/1.4.0/python/programmers-reference/functions/functionspace/FunctionSpace.html

    u_elem  = VectorElement('CG', mesh.ufl_cell(), 2) 
    
    if ( fem.analysisFlag[1] == 1 ):
        p_elem   = FiniteElement('CG', mesh.ufl_cell(), 1) # pore pressure
        V_up     = FunctionSpace(mesh, MixedElement([u_elem, p_elem]))
        V_u, V_p = V_up.split()
    
        W = FunctionSpace(mesh, 'CG', 1) # porosity

        
    if ( fem.analysisFlag[2] == 1 ):
        Ts_elem = FiniteElement('CG', mesh.ufl_cell(), 1) # solid temperature
        V_Ts = FunctionSpace(mesh, Ts_elem)
    
    if ( fem.analysisFlag[3] == 1 ):
        d_elem  = FiniteElement('CG', mesh.ufl_cell(), 1) # phase field
        V_d = FunctionSpace(mesh, d_elem)
        
        
def DefineBDC(fem)
    
    
        
    return BDC

# strain
def epsilon(u):

  return sym(grad(u))

def epsilon_v(u):

  strain = epsilon(u)

  return tr(strain) 
 

# effective stress
def sigma_e(u):

  return lamda*tr(epsilon(u))*Identity(len(u)) + 2.0*mu*epsilon(u)

# total stress
def sigma(u, p):

  stress_e = sigma_e(u)
  stress   = stress_e - p*Identity(2)

  return stress

# Darcy's velocity
def w_Darcy(p):

  w_vel = (-k_w/mu_w)*grad(p)

  return w_vel
