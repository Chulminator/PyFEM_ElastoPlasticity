# ---------------------------------------------------------------- 
# Nonpolar poroelasticity: implementation in FEniCS    
# Written by: Hyoung Suk Suh (h.suh@columbia.edu)     
# ----------------------------------------------------------------       


from ./inlcude import *
from dolfin import *
import numpy as np
import sys
import time as ct


tic = ct.time()


######### option 


######### ReadInput 



# ---------------------------------------------------------------- 
# Input parameters
# ----------------------------------------------------------------
# Mesh and result file names
file_name = 'terzaghi'   # Input directory name

# Applied load
trac0 = -1.0e6 # Applied traction [Pa]

# Material parameters (micropolar poroelasticity)
E     = 1.0e8   # Young's modulus [Pa]
nu    = 0.25    # Poisson's ratio
k_w   = 1.0e-12 # Intrinsic permeability [m2]
mu_w  = 1.0e-3  # Dynamic viscosity of water [Pa.sec]
poro0 = 0.382   # Porosity
K_b_w = 2.2e9   # Bulk modulus of water [Pa]

# Time-stepping parameters
t_i    = 0.0    # [sec]
t_f    = 1000.0 # [sec]
Nsteps = 100
dt     = Constant((t_f-t_i)/Nsteps)





# ---------------------------------------------------------------- 
# Define mesh
# ----------------------------------------------------------------
mesh = Mesh('./'+file_name+'/'+file_name+'.xml')

dim = mesh.geometric_dimension()
mesh_coord = mesh.coordinates()
mesh_xmin  = min(mesh_coord[:,0])
mesh_xmax  = max(mesh_coord[:,0])
mesh_ymin  = min(mesh_coord[:,1])
mesh_ymax  = max(mesh_coord[:,1])

mesh_xcenter = (mesh_xmin + mesh_xmax)/2.





# ---------------------------------------------------------------- 
# Define function spaces
# ----------------------------------------------------------------
u_elem     = VectorElement('CG', mesh.ufl_cell(), 2) # displacement
p_elem     = FiniteElement('CG', mesh.ufl_cell(), 1) # pore water pressure

V    = FunctionSpace(mesh, MixedElement([u_elem, p_elem]))
U, P = V.split()

W = FunctionSpace(mesh, 'CG', 1) # porosity





# ---------------------------------------------------------------- 
# Define boundary conditions
# ----------------------------------------------------------------
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

# Define the loading as an expression depending on t
trac = Expression(("0", "trac0"), trac0=trac0, degree=0)

# Mark boundaries
boundaries = MeshFunction("size_t", mesh, dim-1)
boundaries.set_all(0)

top.mark(boundaries, 1)

ds = Measure("ds")(subdomain_data=boundaries)
n  = FacetNormal(mesh)





# ---------------------------------------------------------------- 
# Define variables & functions
# ----------------------------------------------------------------
lamda = E*nu / ((1.0+nu)*(1.0-2.0*nu))
mu    = E / (2.0*(1.0+nu))

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





# ---------------------------------------------------------------- 
# Define variational form
# ----------------------------------------------------------------
del_x     = TrialFunction(V)
eta, zeta = TestFunctions(V)

x_n1 = Function(V)
u_n1, p_n1 = split(x_n1)

x_n = Function(V)
u_n, p_n = split(x_n)

Res = inner(grad(eta), sigma_e(u_n1)) * dx - inner(p_n1, div(eta)) * dx - dot(eta, trac) * ds(1) \
    + inner(zeta, ((poro0*(1.+epsilon_v(u_n1)))/K_b_w)*(p_n1 - p_n)/dt) * dx + inner(zeta, div((u_n1 - u_n)/dt)) * dx - inner(grad(zeta), w_Darcy(p_n1)) * dx

Jac = derivative(Res, x_n1, del_x) # jacobian

Prob   = NonlinearVariationalProblem(Res, x_n1, BC, Jac)
solver = NonlinearVariationalSolver(Prob)

# Set nonlinear solver parameters
prm = solver.parameters['newton_solver']
prm['error_on_nonconvergence'] = False
prm['relative_tolerance'] = 1e-7
prm['absolute_tolerance'] = 1e-7
prm['maximum_iterations'] = 10





# ---------------------------------------------------------------- 
# Solve system & output results
# ----------------------------------------------------------------
# time stepping
time = np.linspace(t_i, t_f, Nsteps+1)

# output file
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

  poro = project((1.+ epsilon_v(u_n1))*poro0, W) # update porosity

  # record output
  u_n1.rename("Displacement", "label")
  p_n1.rename("Pore water pressure", "label")
  poro.rename("Porosity", "label")

  xdmf_file.write(u_n1, t)
  xdmf_file.write(p_n1, t)
  xdmf_file.write(poro, t)





toc = ct.time() - tic

print(' ')
print('Elapsed CPU time: ', toc, '[sec]')


