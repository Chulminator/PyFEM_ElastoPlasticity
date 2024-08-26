# ---------------------------------------------------------------- 
# Written by: CK, HS in promechanics.org at Columbia University
# ----------------------------------------------------------------       
import numpy as np
import sys
import time as ct
import os.path
from dataclasses import dataclass
from dolfin import *
import Read




def SetSpace(mesh, fem):
#----Element type in FEniCS
#    https://fenicsproject.org/olddocs/dolfin/1.4.0/python/programmers-reference/functions/functionspace/FunctionSpace.html
    
    '''
    fem.analysisFlag[0] # Solid
    fem.analysisFlag[1] # Hydro 
    fem.analysisFlag[2] # Temperature
    fem.analysisFlag[3] # Phase Field
    
    # Variable :
        Elem
            u_elem     : 
            p_elem     : 
            theta_elem : 
        Space
    
    '''
    
    u_elem  = VectorElement('CG', mesh.ufl_cell(), 2) 
    if( fem.solver.lower() == 'staggered' or fem.solver == ''):
        u_elem = VectorElement('Lagrange', mesh.ufl_cell(), 2)
        print("not ready")
        assert False
        exit(1)
        
        if ( fem.analysisFlag[1] == 1 ): 
            p_elem       = FiniteElement('CG', mesh.ufl_cell(), 1) # Pore pressure
        if ( fem.analysisFlag[2] == 1 ): 
            theta_elem   = FiniteElement('CG', mesh.ufl_cell(), 1) # Temperature
        if ( fem.analysisFlag[3] == 1 ): 
            W_damage     = FunctionSpace(mesh, 'CG', 1) # Phase field - Damage
            W_history    = FunctionSpace(mesh, 'DG', 0) # History 
            

    elif( fem.solver.lower() == 'monolithic'):
        if ( fem.analysisFlag == [1, 1, 0, 0]):   # Solid Hydro couple
            u_elem     = VectorElement('CG', mesh.ufl_cell(), 2)
            p_elem     = FiniteElement('CG', mesh.ufl_cell(), 1) # Pore pressure
            
            V_up       = FunctionSpace(mesh, MixedElement([u_elem, p_elem]))
            return V_up
            #V_u, V_p   = V_up.split()
            
            #W_porosity = FunctionSpace(mesh, 'CG', 1) #
        elif ( fem.analysisFlag == [1, 1, 1, 0]): # Thm incomplete
            u_elem     = VectorElement('CG', mesh.ufl_cell(), 2)
            p_elem     = FiniteElement('CG', mesh.ufl_cell(), 1) # Pore pressure
            theta_elem = FiniteElement('CG', mesh.ufl_cell(), 1) # Temperature
            # name: thetaS_elem. or t_elem or tS_elem, T_elem
            
            V_thm       = FunctionSpace(mesh, MixedElement([u_elem, p_elem, theta_elem]))
            V_u, V_p   = V_up.split()
            
            print("not ready")
            assert False
            exit(1)            
                
        elif ( fem.analysisFlag == [1, 0, 1, 0]): # Tm incomplete
            print("not ready")
            assert False
            exit(1)            
        elif ( fem.analysisFlag == [1, 1, 0, 1]): # Solid Hydro couple + Phase-field
            print("not ready")
            assert False
            exit(1)
        elif ( fem.analysisFlag == [1, 0, 1, 1]): # Tm  + Phase-field
            print("not ready")
            assert False
            exit(1)
        elif ( fem.analysisFlag == [1, 1, 1, 1]): # Thm + Phase-field
            print("not ready")
            assert False
            exit(1)
    else:
        print("not ready")
        assert False
        exit(1)
    return 
    

# or do I have to classify this with the combination
    '''
        1. Only solid - Linear elasticity small deformation
        2. U - P - u, p 
        3. THM - Theta, p, u
        4. TM  - Theta, u
        Total 8 including phase field on and off
    '''
def SetNaturalBC(V, mesh, NBC_input):
    NBC = []
    print("on the construction")
    assert False
    for NBC_itr in fem.NBC:
        Submesh = GetSubMesh(mesh, NBC_itr[2:end])
        f = Expression(NBC_itr[1],t=0) # I should see this is gonna work if NBC_itr[1] == 0
        if NBC_itr[0].upper() == TRAC1 : 
            NBC.append(DirichletBC(V_u.sub(0), f, Submesh))
        elif NBC_itr[0].upper() == TRAC2 : 
            NBC.append(DirichletBC(V_u.sub(1), f, Submesh))
        elif NBC_itr[0].upper() == F1 : 
            NBC.append(DirichletBC(V_p, f, Submesh))
        elif NBC_itr[0].upper() == F2 :
            NBC.append(DirichletBC(V_theta, f, Submesh))
        else:
            print("other option for NaturalBC is on the construction")
            #print("Check EssentialBD in *.inp")
            print("\tAvaiable option:\n")
            print("SetNaturalBC in DefineProblemSet")
            ErrorMassage("U1\tU2\tP\tT\tD") 
            
    return NBC
    
def SetEssentialBC(V, mesh, EBC_input):
    

    EBC = []
    for EBC_itr in EBC_input:
        Submesh = GetSubMesh(mesh, EBC_itr[2:end])
        
        if EBC_itr[0].upper() == U1 : 
            EBC.append(DirichletBC(V_u.sub(0), Constant(float(EBC_itr[1])), Submesh))
            
        elif EBC_itr[0].upper() == U2 : 
            EBC.append(DirichletBC(V_u.sub(1), Constant(float(EBC_itr[1])), Submesh))
            
        elif EBC_itr[0].upper() == P : 
            EBC.append(DirichletBC(V_p, Constant(float(EBC_itr[1])), Submesh))
            
        elif EBC_itr[0].upper() == T :
            EBC.append(DirichletBC(V_theta, Constant(float(EBC_itr[1])), Submesh))
            
        elif EBC_itr[0].upper() == D : 
            EBC.append(DirichletBC(V_d, Constant(float(EBC_itr[1])), Submesh))
            
        else:
            print("Check EssentialBD in *.inp")
            print("\tAvaiable option:\n")
            ErrorMassage("U1\tU2\tP\tT\tD") 
            
    return EBC
                
            
def GetSubMesh(mesh, option):
    dim = mesh.geometric_dimension()
    mesh_coord = mesh.coordinates()
    if(option[0].lower() == "top" or option.lower() == "up"):
        BDCoord = max(mesh_coord[:,1])
        Submesh     = CompiledSubDomain("near(x[1], BDCoord) && on_boundary", BDCoord = BDCoord)
        
    elif(option[0].lower() == "down" or option.lower() == "bottom"):
        BDCoord = min(mesh_coord[:,1])
        Submesh     = CompiledSubDomain("near(x[1], BDCoord) && on_boundary", BDCoord = BDCoord)
    elif(option[0].lower() == "left"):
        BDCoord = min(mesh_coord[:,0])
        Submesh     = CompiledSubDomain("near(x[0], BDCoord) && on_boundary", BDCoord = BDCoord)
    elif(option[0].lower() == "right"):
        BDCoord = max(mesh_coord[:,0])
        Submesh     = CompiledSubDomain("near(x[0], BDCoord) && on_boundary", BDCoord = BDCoord)
    else:
        print("GetSubMesh in DefineProblemSet.py")
        print("On the construction")
        assert false
        #Submesh     = CompiledSubDomain("near(x[1], BDCoord) && on_boundary", BDCoord = BDCoord)
        
    return Submesh

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
