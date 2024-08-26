# ---------------------------------------------------------------- 
# Written by: Hyoung Suk Suh (h.suh@columbia.edu)     
# sub by: Chulmin Kweon      (Chulmin.Kweon@columbia.edu)
# ----------------------------------------------------------------       



from dolfin import *
import numpy as np
import sys
sys.path.append("./include")
import time as ct


if __name__ == '__main__':
    
    
    tic = ct.time()
# ---------------------------------------------------------------- 
# Input parameters
# ----------------------------------------------------------------
    import Read
    fem = Read.Model()
    Solid      = Read.Solid()
    Hydro      = Read.Hydro()
    Thermo     = Read.Thermo()
    PhaseField = Read.PhaseField()
    mesh = Read.ReadInput("./input/terzaghi/NonpolarElasticity.inp",fem, Solid, Hydro, Thermo, PhaseField)
    
    fem.show()
    Solid.show()
    Hydro.show()
# ----------------------------------------------------------------
# Define function spaces & Boundary condition
# ----------------------------------------------------------------
# displacement displacement order question
    import SetProblem
    V, EBC = SetProblem.SetSpaceAndEBC(mesh, fem)
    NBC = SetProblem.SetNaturalBC(V, mesh, fem.NBC)
# ---------------------------------------------------------------- 
# Define variational form
# ----------------------------------------------------------------

    del_x     = TrialFunction(V)
    eta, zeta = TestFunctions(V)

    x_n1 = Function(V)
    u_n1, p_n1 = split(x_n1)

    x_n = Function(V)
    u_n, p_n = split(x_n)


    U, P = V.split()

    
    
    tmp = NBC[0]
    #print(tmp[0][0],"\n")
    #print(tmp[0][1],"\n")
    #exit(1)
    #Res = inner(grad(eta), Solid.sigma_e(u_n1)) * dx - inner(p_n1, div(eta)) * dx - inner(eta, f) * ds(1) \
    Res = inner(grad(eta), Solid.sigma_e(u_n1)) * dx - inner(p_n1, div(eta)) * dx - inner(eta, tmp[0][1]) * tmp[0][0] \
        + inner(zeta, ((Hydro.poro*(1.+Solid.epsilon_v(u_n1)))/Hydro.K_b_f)*(p_n1 - p_n)/fem.dt) * dx \
        + inner(zeta, div((u_n1 - u_n)/fem.dt)) * dx - inner(grad(zeta), Hydro.w_Darcy(p_n1)) * dx

    Jac = derivative(Res, x_n1, del_x) # jacobian

    Prob   = NonlinearVariationalProblem(Res, x_n1, EBC, Jac)
    solver = NonlinearVariationalSolver(Prob)

    prm = solver.parameters['newton_solver']
    prm['error_on_nonconvergence'] = False
    prm['relative_tolerance'] = 1e-7
    prm['absolute_tolerance'] = 1e-7
    prm['maximum_iterations'] = 10


# ---------------------------------------------------------------- 
# Solve system & output results
# ----------------------------------------------------------------
    # time stepping
    #time = np.linspace(0., fem.totaltime, fem.totalstep+1)
    time = np.linspace(0, fem.totaltime, fem.totalstep+1)

    # output filew
    xdmf_file = XDMFFile("terzaghi_org.xdmf")
    xdmf_file.parameters["flush_output"] = True
    xdmf_file.parameters["functions_share_mesh"] = True
    xdmf_file.parameters["rewrite_function_mesh"] = False
    
    W = FunctionSpace(mesh, 'CG', 1) # porosity
    for (i, dt) in enumerate(np.diff(time)):

        t = time[i+1]
        print('-----------------------------------------')
        print('>> t =', t, '[sec]')

        # solve system
        solver.solve()

        # update
        x_n.assign(x_n1)

        # update internal variable
        u_n1, p_n1 = x_n1.split()

        poro = project((1.+ Solid.epsilon_v(u_n1))*Hydro.poro, W) # update porosity

        # record output
        u_n1.rename("Displacement", "label")
        p_n1.rename("Pore water pressure", "label")
        poro.rename("Porosity", "label")

        xdmf_file.write(u_n1, t)
        xdmf_file.write(p_n1, t)
        xdmf_file.write(poro, t)



    
    fem.LogClose
    toc = ct.time() - tic
    print('\nElapsed CPU time: ', toc, '[sec]')
