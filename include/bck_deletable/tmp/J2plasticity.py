# ==============================================================
# Step 4: run stress point simulation
# written by: Hyoung Suk Suh & Chulmin Kweon (Columbia Univ.)
# ==============================================================

# Import necessary packages and functions
import os
import matplotlib.pyplot as plt
import autograd.numpy as np

from autograd import elementwise_grad as egrad

from util.tensor_operations import *
from util.coordinate_transforms import *

# INPUT --------------------------------------------------------
# Material properties
# >> elastic properties (for both NN and benchmark)
E  = 200.0e3     # Young's modulus [MPa]
nu = 0.3         # Poisson's ratio

# >> yielding properties (for benchmark)
sigma_y = 250.0  # Initial yield stress [MPa]
N       = 2.0    # vonMises parameter
#N       = 1.0    # Hosford parameter
#N       = 1.0    # Hosford parameter
n_hard  = 0.2    # Hardening parameter

# Loading steps & increment
Nstep   = 720    # loading steps
eps_inc = 1.0e-4 # strain increment

# Newton-Raphson parameters
tol     = 1e-9
maxiter = 10

# Define loading direction (strain controlled)
stress_increment = np.array([[eps_inc,  0.0,          0.0],
                             [0.0,     -0.5*eps_inc,  0.0],
                             [0.0,      0.0,         -0.5*eps_inc]])

# File name
file_name = "Tresca"
file_name = "vonMises"
# --------------------------------------------------------------




# Initialize variables
# >> NN-prediction
sigma = np.zeros((3,3)) # stress

eps_e = np.zeros((3,3)) # elastic strain
eps_p = np.zeros((3,3)) # plastic strain
eps   = eps_e + eps_p   # strain

lamda = 0 # plastic multiplier

# >> benchmark
sigma_bm = np.zeros((3,3)) # stress

eps_e_bm = np.zeros((3,3)) # elastic strain
eps_p_bm = np.zeros((3,3)) # plastic strain
eps_bm   = eps_e + eps_p   # strain

lamda_bm = 0 # plastic multiplier




# Define identity tensors
I   = np.eye(3)
II  = identity_4(3)
IxI = tensor_oMult(I,I)




# Define elasticity model
K  = E / (3*(1-2*nu))
mu = E / (2*(1+nu))

a = K + (4/3)*mu
b = K - (2/3)*mu

Ce_principal = np.array([[a, b, b],
                         [b, a, b],
                         [b, b, a]])

dsig1depse1 = a
dsig1depse2 = b
dsig1depse3 = b
dsig2depse1 = b
dsig2depse2 = a
dsig2depse3 = b
dsig3depse1 = b
dsig3depse2 = b
dsig3depse3 = a



# Benchmark - stress point simulation ================================================
# Define yield function
def f_benchmark(sigma1, sigma2, sigma3, lamda):

  phi = (1/2)*(np.abs(sigma1-sigma2)**N + \
               np.abs(sigma2-sigma3)**N + \
               np.abs(sigma3-sigma1)**N)

  kappa = sigma_y*(1 + E*lamda/sigma_y)**n_hard

  return phi**(1/N) - kappa

# Define first order derivatives
get_dfdsig1 = egrad(f_benchmark, 0)
get_dfdsig2 = egrad(f_benchmark, 1)
get_dfdsig3 = egrad(f_benchmark, 2)

def df_benchmark(sigma1, sigma2, sigma3, lamda):

  x, y, z= sigma1, sigma2, sigma3
  denom = np.sqrt(2)*((x-y)**2+(y-z)**2+(z-x)**2)**(1/2)
  dfdsig1 = (2*x-y-z)/ denom
  dfdsig2 = (2*y-z-x)/ denom
  dfdsig3 = (2*z-x-y)/ denom

  #dfdsig1 = get_dfdsig1(sigma1, sigma2, sigma3, lamda)
  #dfdsig2 = get_dfdsig2(sigma1, sigma2, sigma3, lamda)
  #dfdsig3 = get_dfdsig3(sigma1, sigma2, sigma3, lamda)

  return dfdsig1, dfdsig2, dfdsig3

# Define second order derivatives
def df2_benchmark(sigma1, sigma2, sigma3, lamda):

  #dist = 1e-3
  #dfdsig1, dfdsig2, dfdsig3 = df_benchmark(sigma1, sigma2, sigma3, lamda)
  #dfdsig1_s1dist, dfdsig2_s1dist, dfdsig3_s1dist = df_benchmark(sigma1+dist, sigma2, sigma3, lamda)
  #dfdsig1_s2dist, dfdsig2_s2dist, dfdsig3_s2dist = df_benchmark(sigma1, sigma2+dist, sigma3, lamda)
  #dfdsig1_s3dist, dfdsig2_s3dist, dfdsig3_s3dist = df_benchmark(sigma1, sigma2, sigma3+dist, lamda)

  #d2fdsig1dsig1 = (dfdsig1_s1dist - dfdsig1) / dist
  #d2fdsig2dsig2 = (dfdsig2_s2dist - dfdsig2) / dist
  #d2fdsig3dsig3 = (dfdsig3_s3dist - dfdsig3) / dist

  #d2fdsig1dsig2 = (dfdsig1_s2dist - dfdsig1) / dist
  #d2fdsig2dsig3 = (dfdsig2_s3dist - dfdsig2) / dist
  #d2fdsig3dsig1 = (dfdsig3_s1dist - dfdsig3) / dist

  x, y, z= sigma1, sigma2, sigma3
  denom = np.sqrt(2)*((x-y)**2+(y-z)**2+(z-x)**2)**(3/2)
  d2fdsig1dsig1 = 3*(y - z)**2 / denom
  d2fdsig2dsig2 = 3*(z - x)**2 / denom
  d2fdsig3dsig3 = 3*(x - y)**2 / denom
  d2fdsig1dsig2 = 3*(z-x)*(y-z) / denom
  d2fdsig2dsig3 = 3*(x-y)*(z-x) / denom
  d2fdsig3dsig1 = 3*(y-x)*(z-y) / denom
  return d2fdsig1dsig1, d2fdsig2dsig2, d2fdsig3dsig3, d2fdsig1dsig2, d2fdsig2dsig3, d2fdsig3dsig1

# Perform material point simulation
print(":: Stress-point simulation (benchmark) ::")

sigma11_bm = np.zeros(Nstep+1)
eps11_bm   = np.zeros(Nstep+1)
deps = stress_increment

for i in range(Nstep):

  print("Loading step [",i+1,"] ---------------------------------------")

  if i == 50:
    deps = -deps
  elif i == 150:
    deps = -deps
  elif i == 300:
    deps = -deps
  elif i == 500:
    deps = -deps

  # [1] Compute trial strain
  eps_e_tr_bm = eps_e_bm + deps

  eps_e_tr_principal_mag_bm, eps_e_tr_principal_vec_bm = np.linalg.eig(eps_e_tr_bm)

  eps_e_tr1_bm = eps_e_tr_principal_mag_bm[0]
  eps_e_tr2_bm = eps_e_tr_principal_mag_bm[1]
  eps_e_tr3_bm = eps_e_tr_principal_mag_bm[2]

  n1_bm = eps_e_tr_principal_vec_bm[:,0]
  n2_bm = eps_e_tr_principal_vec_bm[:,1]
  n3_bm = eps_e_tr_principal_vec_bm[:,2]

  # [2] Compute trial stress
  sigma_tr_principal_mag_bm = np.inner(Ce_principal, eps_e_tr_principal_mag_bm)

  sigma_tr1_bm = sigma_tr_principal_mag_bm[0]
  sigma_tr2_bm = sigma_tr_principal_mag_bm[1]
  sigma_tr3_bm = sigma_tr_principal_mag_bm[2]

  sigma_tr_bm = sigma_tr1_bm*np.tensordot(n1_bm,n1_bm,axes=0) \
              + sigma_tr2_bm*np.tensordot(n2_bm,n2_bm,axes=0) \
              + sigma_tr3_bm*np.tensordot(n3_bm,n3_bm,axes=0)

  # [3] Check yielding
  f_bm = f_benchmark(sigma_tr1_bm, sigma_tr2_bm, sigma_tr3_bm, lamda_bm)

  # [3.1] If f <= 0, elastic.
  if f_bm <= 0:
    print(">> Elastic!")

    # Update stress & strain
    sigma_bm = sigma_tr_bm

    eps_e_bm = eps_e_tr_bm
    eps_bm   = eps_e_bm + eps_p_bm

  # [3.2] If f > 0, plastic.
  else:
    print(">> Plastic!")

    # Initialize variables
    eps_e_principal_mag_bm, eps_e_principal_vec_bm = np.linalg.eig(eps_e_bm)

    eps_e1_bm = eps_e_principal_mag_bm[0]
    eps_e2_bm = eps_e_principal_mag_bm[1]
    eps_e3_bm = eps_e_principal_mag_bm[2]
    dlamda_bm  = 0

    x_bm = np.zeros(4) # target variables
    x_bm[0] = eps_e1_bm
    x_bm[1] = eps_e2_bm
    x_bm[2] = eps_e3_bm
    x_bm[3] = dlamda_bm

    # Newton-Raphson iteration (return mapping)
    for ii in range(maxiter):

      # Initialize residual and jacobian
      res_bm = np.zeros(4)
      jac_bm = np.zeros((4,4))

      # Current strain
      eps_e1_current_bm = x_bm[0]
      eps_e2_current_bm = x_bm[1]
      eps_e3_current_bm = x_bm[2]

      # Current stress
      sigma1_current_bm = a*eps_e1_current_bm + b*eps_e2_current_bm + b*eps_e3_current_bm
      sigma2_current_bm = b*eps_e1_current_bm + a*eps_e2_current_bm + b*eps_e3_current_bm
      sigma3_current_bm = b*eps_e1_current_bm + b*eps_e2_current_bm + a*eps_e3_current_bm

      # Current lamda
      lamda_current_bm = lamda_bm + x_bm[3]

      # Update derivatives
      # >> First order derivatives
      dfdsig1_bm, dfdsig2_bm, dfdsig3_bm \
        = df_benchmark(sigma1_current_bm, sigma2_current_bm, sigma3_current_bm, lamda_current_bm)

      # >> Second order derivatives
      d2fdsig1dsig1_bm, d2fdsig2dsig2_bm, d2fdsig3dsig3_bm, d2fdsig1dsig2_bm, d2fdsig2dsig3_bm, d2fdsig3dsig1_bm \
        = df2_benchmark(sigma1_current_bm, sigma2_current_bm, sigma3_current_bm, lamda_current_bm)

      # Update residual
      res_bm[0] = x_bm[0] - eps_e_tr1_bm + x_bm[3]*dfdsig1_bm
      res_bm[1] = x_bm[1] - eps_e_tr2_bm + x_bm[3]*dfdsig2_bm
      res_bm[2] = x_bm[2] - eps_e_tr3_bm + x_bm[3]*dfdsig3_bm
      res_bm[3] = f_benchmark(sigma1_current_bm, sigma2_current_bm, sigma3_current_bm, lamda_current_bm)

      # Update Jacobian ***
      jac_bm[0,0] = 1 + x_bm[3]*(d2fdsig1dsig1_bm*dsig1depse1 + d2fdsig1dsig2_bm*dsig2depse1 + d2fdsig3dsig1_bm*dsig3depse1)
      jac_bm[0,1] =     x_bm[3]*(d2fdsig1dsig1_bm*dsig1depse2 + d2fdsig1dsig2_bm*dsig2depse2 + d2fdsig3dsig1_bm*dsig3depse2)
      jac_bm[0,2] =     x_bm[3]*(d2fdsig1dsig1_bm*dsig1depse3 + d2fdsig1dsig2_bm*dsig2depse3 + d2fdsig3dsig1_bm*dsig3depse3)
      jac_bm[0,3] = dfdsig1_bm

      jac_bm[1,0] =     x_bm[3]*(d2fdsig1dsig2_bm*dsig1depse1 + d2fdsig2dsig2_bm*dsig2depse1 + d2fdsig2dsig3_bm*dsig3depse1)
      jac_bm[1,1] = 1 + x_bm[3]*(d2fdsig1dsig2_bm*dsig1depse2 + d2fdsig2dsig2_bm*dsig2depse2 + d2fdsig2dsig3_bm*dsig3depse2)
      jac_bm[1,2] =     x_bm[3]*(d2fdsig1dsig2_bm*dsig1depse3 + d2fdsig2dsig2_bm*dsig2depse3 + d2fdsig2dsig3_bm*dsig3depse3)
      jac_bm[1,3] = dfdsig2_bm

      jac_bm[2,0] =     x_bm[3]*(d2fdsig3dsig1_bm*dsig1depse1 + d2fdsig2dsig3_bm*dsig2depse1 + d2fdsig3dsig3_bm*dsig3depse1)
      jac_bm[2,1] =     x_bm[3]*(d2fdsig3dsig1_bm*dsig1depse2 + d2fdsig2dsig3_bm*dsig2depse2 + d2fdsig3dsig3_bm*dsig3depse2)
      jac_bm[2,2] = 1 + x_bm[3]*(d2fdsig3dsig1_bm*dsig1depse3 + d2fdsig2dsig3_bm*dsig2depse3 + d2fdsig3dsig3_bm*dsig3depse3)
      jac_bm[2,3] = dfdsig3_bm

      jac_bm[3,0] = dfdsig1_bm*dsig1depse1 + dfdsig2_bm*dsig2depse1 + dfdsig3_bm*dsig3depse1
      jac_bm[3,1] = dfdsig1_bm*dsig1depse2 + dfdsig2_bm*dsig2depse2 + dfdsig3_bm*dsig3depse2
      jac_bm[3,2] = dfdsig1_bm*dsig1depse3 + dfdsig2_bm*dsig2depse3 + dfdsig3_bm*dsig3depse3
      jac_bm[3,3] = 0

      # Solve system of equations
      dx_bm = np.linalg.solve(jac_bm, -res_bm) # increment of target variables

      # Update x
      x_bm = x_bm + dx_bm

      # Compute error
      err_bm = np.linalg.norm(dx_bm)

      print(" Newton iter.",ii, ": err =", err_bm)

      if err_bm < tol:
        break

    # Update strain
    eps_bm   = eps_bm + deps
    eps_e_bm = x_bm[0]*np.tensordot(n1_bm,n1_bm,axes=0) + x_bm[1]*np.tensordot(n2_bm,n2_bm,axes=0) + x_bm[2]*np.tensordot(n3_bm,n3_bm,axes=0)
    eps_p_bm = eps_bm - eps_e_bm
    lamda_bm = lamda_bm + x_bm[3]

    # Update stress
    sigma1_bm = a*x_bm[0] + b*x_bm[1] + b*x_bm[2]
    sigma2_bm = b*x_bm[0] + a*x_bm[1] + b*x_bm[2]
    sigma3_bm = b*x_bm[0] + b*x_bm[1] + a*x_bm[2]
    sigma_bm  = sigma1_bm*np.tensordot(n1_bm,n1_bm,axes=0) + sigma2_bm*np.tensordot(n2_bm,n2_bm,axes=0) + sigma3_bm*np.tensordot(n3_bm,n3_bm,axes=0)


  # [4] Record stress and strain
  sigma11_bm[i+1] = sigma_bm[0,0]
  eps11_bm[i+1]   = eps_bm[0,0]
# ====================================================================================




# Plot stress-strain curve
plt.figure(0,figsize=(7,7))
plt.plot(eps11_bm, sigma11_bm, 'b-', linewidth=2.5, label="Benchmark")
#plt.plot(eps11, sigma11, 'r--', linewidth=1.0, label="NN prediction")
plt.xlabel(r'$\epsilon_{11}$', fontsize=15)
plt.ylabel(r'$\sigma_{11}$ [MPa]', fontsize=15)
plt.axhline(0, color = 'k',alpha = 0.5)
plt.axvline(0, color = 'k',alpha = 0.5)
plt.xlim(-0.015, 0.015)
plt.ylim(-600, 600)
plt.show()
