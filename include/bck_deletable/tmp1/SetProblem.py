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




def SetSpaceAndEBC(mesh, fem):
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
    # EBC
    '''
    
    EBC =[]
    
    #u_elem  = VectorElement('CG', mesh.ufl_cell(), 2) 
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
            V_u, V_p   = V_up.split()
            for EBC_itr in fem.EBC:
                if(EBC_itr[0].lower() == 'u1'):
                    EBC.append(SetEssentialBC(V_u.sub(0), mesh, EBC_itr))
                elif(EBC_itr[0].lower() == 'u2'):
                    EBC.append(SetEssentialBC(V_u.sub(1), mesh, EBC_itr))
                elif(EBC_itr[0].lower() == 'p'):
                    EBC.append(SetEssentialBC(V_p, mesh, EBC_itr))
                else:
                    print(EBC_itr)
                    print("EBC of UP alaysis: U1, U2, P in EssentialBD in *.inp")
                    assert False
                
                
            return V_up, EBC
            
            
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
    '''
        1. Only solid - Linear elasticity small deformation
        2. U - P - u, p 
        3. THM - Theta, p, u
        4. TM  - Theta, u
        Total 8 including phase field on and off
    '''
def SetNaturalBC(V, mesh, NBC_input):
    # MeshFunction
    # https://fenicsproject.org/olddocs/dolfin/1.6.0/python/programmers-reference/cpp/mesh/MeshFunction.html
    #https://scicomp.stackexchange.com/questions/32647/how-to-use-meshfunction-in-fenics-dolfin
    
    dim = mesh.geometric_dimension()
    NBC = [0, 0, 0]
    NBC_u = []
    NBC_th = []
    NBC_p = []
    for NBC_itr in NBC_input:
        
        
        if NBC_itr[0].upper() == "TRAC1" : 
            degree_input=2 # edit here!!
                
            f = Expression(("trac0", "0.0"), trac0 = float(NBC_itr[1]), t=0, degree=degree_input)
            #boundaries = MeshFunction("size_t", mesh, dim-1)
            boundaries = MeshFunction("size_t", mesh, 1)
            boundaries.set_all(0) # Initialize the function to zero
            Submesh = GetSubMesh(mesh, NBC_itr[2:])
            Submesh.mark(boundaries, 1)
            ds = Measure("ds")(subdomain_data=boundaries)
            NBC_u.append([ds(0), f])
            
        elif NBC_itr[0].upper() == "TRAC2" : 
            #print(NBC_itr)
            degree_input=2
            #edit here
                
            f = Expression(("0.0", "trac0"), trac0 = float(NBC_itr[1]), t=0, degree=degree_input)
            #boundaries = MeshFunction("size_t", mesh, dim-1)
            boundaries = MeshFunction("size_t", mesh, 1)
            boundaries.set_all(0) # Initialize the function to zero
            Submesh = GetSubMesh(mesh, NBC_itr[2:])
            Submesh.mark(boundaries, 1)
            ds = Measure("ds")(subdomain_data=boundaries)
            NBC_u.append([ds(1), f])
            #NBC_u.append([ds(1), f])
            '''
            print(ds(1),"\n")
            print(f,"\n")
            print(NBC_u,"\n")
            print(NBC_u[0],"\n")
            tmp = NBC_u[0]
            print(tmp[0],"\n")
            print(tmp[1],"\n")
            exit(1)
            '''
            
        elif NBC_itr[0].upper() == "F1" : 
            print("SetNaturalBC F1")
            print("Under construction")
            exit(1)
        elif NBC_itr[0].upper() == "F2" :
            print("SetNaturalBC F2")
            print("Under construction")
            exit(1)
        else:
            print("other option for NaturalBC is under construction")
            #print("Check EssentialBD in *.inp")
            print("\tAvaiable option:\n")
            print("SetNaturalBC in DefineProblemSet")
            ErrorMassage("U1\tU2\tP\tT\tD") 
    NBC[0] = NBC_u
    NBC[1] = NBC_th
    NBC[2] = NBC_p
    return NBC
    
def SetEssentialBC(V, mesh, EBC_input):
    Submesh = GetSubMesh(mesh, EBC_input[2:])
    # edit
    # expression -> time dependent EBC
    
    EBC = DirichletBC(V, Constant(float(EBC_input[1])), Submesh)
                
    return EBC
                
            
def GetSubMesh(mesh, option):
    #print(option)
    dim = mesh.geometric_dimension()
    mesh_coord = mesh.coordinates()
    if(option[0].lower() == "top" or option[0].lower() == "up"):
        BDCoord = max(mesh_coord[:,1])
        Submesh     = CompiledSubDomain("near(x[1], BDCoord) && on_boundary", BDCoord = BDCoord)
        
    elif(option[0].lower() == "down" or option[0].lower() == "bottom"):
        BDCoord = min(mesh_coord[:,1])
        Submesh     = CompiledSubDomain("near(x[1], BDCoord) && on_boundary", BDCoord = BDCoord)
    elif(option[0].lower() == "left"):
        BDCoord = min(mesh_coord[:,0])
        Submesh     = CompiledSubDomain("near(x[0], BDCoord) && on_boundary", BDCoord = BDCoord)
    elif(option[0].lower() == "right"):
        BDCoord = max(mesh_coord[:,0])
        Submesh     = CompiledSubDomain("near(x[0], BDCoord) && on_boundary", BDCoord = BDCoord)
    else:
        #some other option - specitying x and y
        print("GetSubMesh in DefineProblemSet.py")
        print("under construction")
        assert false
        #Submesh     = CompiledSubDomain("near(x[1], BDCoord) && on_boundary", BDCoord = BDCoord)
        
    return Submesh


def SetVariationalForm_Solid(V, NBC, Solid, fem):
    
    del_x     = TrialFunction(V)
    eta, zeta = TestFunctions(V)

    x_n1 = Function(V)
    u_n1, p_n1 = split(x_n1)

    x_n = Function(V)
    u_n, p_n = split(x_n)


    U, P = V.split()
    return
    
def SetVariationalForm_Hydro(V, NBC, Hydro, fem):
    return
    
    
def SetVariationalForm_Thermo(V, NBC, Thermo, fem):
    return
    
    
def SetVariationalForm_PhaseField(V, NBC, PhaseField, fem):    
    return
    
    
def SetVariationalForm(V, NBC, EBC, Solid, Hydro, Thermo, PhaseField, fem):
    
    del_x     = TrialFunction(V)
    eta, zeta = TestFunctions(V)

    x_n1 = Function(V)
    u_n1, p_n1 = split(x_n1)

    x_n = Function(V)
    u_n, p_n = split(x_n)


    U, P = V.split()
    
    ## Solid
    Res_Solid = inner(grad(eta), Solid.sigma_e(u_n1)) * dx - inner(p_n1, div(eta)) * dx
    
    for NBC_Solid in NBC[0]:
        #print(NBC_Solid)
        #exit(1)
        Res_Solid -= inner (eta,NBC_Solid[1]) * NBC_Solid[0]
        
    #Res = - inner(p_n1, div(eta)) * dx - inner(eta, tmp[0][1]) * tmp[0][0]    
        
    ## Hydro
    
    Res_Hydro = inner(zeta, ((Hydro.poro*(1.+Solid.epsilon_v(u_n1)))/Hydro.K_b_f)*(p_n1 - p_n)/fem.dt) * dx \
        + inner(zeta, div((u_n1 - u_n)/fem.dt)) * dx - inner(grad(zeta), Hydro.w_Darcy(p_n1)) * dx
    
    for NBC_Hydro in NBC[1]:
        Res_Hydro -= inner (zeta,NBC_Hydro[0][1]) * NBC_Hydro[0][0]
        

    ## Thermo
    '''
    Res_Thermo = 
    for NBC_Thermo in NBC[2]:
        Res_Hydro -= inner (zeta,NBC_Thermo[0][1]) * NBC_Thermo[0][0]
    '''
        
    ## PhaseField
        
    Res = Res_Solid + Res_Hydro

    Jac = derivative(Res, x_n1, del_x) # jacobian

    Prob   = NonlinearVariationalProblem(Res, x_n1, EBC, Jac)
    solver = NonlinearVariationalSolver(Prob)

    prm = solver.parameters['newton_solver']
    prm['error_on_nonconvergence'] = False
    prm['relative_tolerance'] = 1e-7
    prm['absolute_tolerance'] = 1e-7
    prm['maximum_iterations'] = 10
    return solver, x_n, x_n1
