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
    f = Expression((0, fem.NBC[0][1]), t=0, degree=0)
    boundaries = MeshFunction("size_t", mesh, 2)
    boundaries = MeshFunction("size_t", mesh, 1)
    boundaries.set_all(0) # Initialize the function to zero
    Submesh = SetProblem.GetSubMesh(mesh, fem.NBC[0][2:])
    Submesh.mark(boundaries, 1)
    ds = Measure("ds")(subdomain_data=boundaries)
    
    U, P = V.split()
    mesh_coord = mesh.coordinates()
    mesh_xmin  = min(mesh_coord[:,0])
    mesh_xmax  = max(mesh_coord[:,0])
    mesh_ymin  = min(mesh_coord[:,1])
    mesh_ymax  = max(mesh_coord[:,1])
    
    top    = CompiledSubDomain("near(x[1], mesh_ymax) && on_boundary", mesh_ymax = mesh_ymax)
    bottom = CompiledSubDomain("near(x[1], mesh_ymin) && on_boundary", mesh_ymin = mesh_ymin)
    left   = CompiledSubDomain("near(x[0], mesh_xmin) && on_boundary", mesh_xmin = mesh_xmin)
    right  = CompiledSubDomain("near(x[0], mesh_xmax) && on_boundary", mesh_xmax = mesh_xmax)

    # Constrained boundaries
    BC_u_left   = DirichletBC(U.sub(0), Constant(0.0), left)
    BC_u_right  = DirichletBC(U.sub(0), Constant(0.0), right)
    BC_u_bottom = DirichletBC(U.sub(1), Constant(0.0), bottom)

    BC_p_top = DirichletBC(P, Constant(0.0), top)

    BC = [BC_u_left, BC_u_right, BC_u_bottom, BC_p_top]
    #print(fem.NBC[0][1])
    #print(type(fem.NBC[0][1]))
    #exit(1)
    
    #dt     = Constant(fem.totaltime/fem.totalstep)
    t_i    = 0.0    # [sec]
    t_f    = 1000.0 # [sec]
    Nsteps = 100
    dt     = Constant((t_f-t_i)/Nsteps)

    del_x     = TrialFunction(V)
    eta, zeta = TestFunctions(V)

    x_n1 = Function(V)
    u_n1, p_n1 = split(x_n1)

    x_n = Function(V)
    u_n, p_n = split(x_n)


    U, P = V.split()

    W = FunctionSpace(mesh, 'CG', 1) # porosity
    
    tmp = NBC[0]
    #print(tmp[0][0],"\n")
    #print(tmp[0][1],"\n")
    #exit(1)
    #Res = inner(grad(eta), Solid.sigma_e(u_n1)) * dx - inner(p_n1, div(eta)) * dx - inner(eta, tmp[0][1]) * tmp[0][0] \
    Res = inner(grad(eta), Solid.sigma_e(u_n1)) * dx - inner(p_n1, div(eta)) * dx - inner(eta, f) * ds(1) \
        + inner(zeta, ((Hydro.poro*(1.+Solid.epsilon_v(u_n1)))/Hydro.K_b_f)*(p_n1 - p_n)/dt) * dx \
        + inner(zeta, div((u_n1 - u_n)/dt)) * dx - inner(grad(zeta), Hydro.w_Darcy(p_n1)) * dx

    Jac = derivative(Res, x_n1, del_x) # jacobian

    Prob   = NonlinearVariationalProblem(Res, x_n1, BC, Jac)
    solver = NonlinearVariationalSolver(Prob)



# ---------------------------------------------------------------- 
# Solve system & output results
# ----------------------------------------------------------------
    # time stepping
    #time = np.linspace(0., fem.totaltime, fem.totalstep+1)
    time = np.linspace(t_i, t_f, Nsteps+1)

    # output filew
    xdmf_file = XDMFFile("terzaghi_org.xdmf")
    xdmf_file.parameters["flush_output"] = True
    xdmf_file.parameters["functions_share_mesh"] = True
    xdmf_file.parameters["rewrite_function_mesh"] = False

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
