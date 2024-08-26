# ==============================================================
# Run stress point simulation
# written by: Hyoung Suk Suh & Chulmin Kweon (Columbia Univ.)
# ==============================================================

# Import necessary packages and functions
import os
import matplotlib.pyplot as plt
import numpy as np
from J2Plasticity import *
import sys
sys.path.append("../.")
from util.tensor_operations import *
from util.coordinate_transforms import *
from Program import *

# INPUT --------------------------------------------------------
# Material properties
# >> elastic properties (for both NN and benchmark)
E  = 200.0e3     # Young's modulus [MPa]
nu = 0.3         # Poisson's ratio

# Loading steps & increment
Nstep   = 100    # loading steps
eps_inc = 0.03 /100 # strain increment

# Newton-Raphson parameters
tol     = 1e-9
maxiter = 10

# Define loading direction (strain controlled)
stress_increment = np.array([[eps_inc,  0.0,  0.0],
                             [0.0,      0.0,  0.0],
                             [0.0,      0.0,  0.0]])
Fem         = Model()
Fem.MatProp = [E, nu, 250.0, 0.2]
Fem.E       = E
Fem.nu      = nu
Fem.K       = E / (3*(1.0-2.0*nu))
Fem.mu      = E / (2.0*(1.0+nu))

Plasticity  = MyPlasticity(Fem)
# --------------------------------------------------------------
# Initialize variables
# >> NN-prediction
sigma = np.zeros((3,3)) # stress

eps_e = np.zeros((3,3)) # elastic strain
eps_p = np.zeros((3,3)) # plastic strain
eps   = eps_e + eps_p   # strain

lamda = 0 # plastic multiplier

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




# Define yield function & plastic flow direction

# Perform material point simulation
sigma11 = np.zeros(Nstep+1)
eps11   = np.zeros(Nstep+1)

deps = stress_increment

for i in range(Nstep):

  print("Loading step [",i+1,"] ---------------------------------------")

  # [1] Compute trial strain
  eps_e_tr = eps_e + deps

  eps_e_tr_principal_mag, eps_e_tr_principal_vec = np.linalg.eig(eps_e_tr)

  eps_e_tr1 = eps_e_tr_principal_mag[0]
  eps_e_tr2 = eps_e_tr_principal_mag[1]
  eps_e_tr3 = eps_e_tr_principal_mag[2]

  n1 = eps_e_tr_principal_vec[:,0]
  n2 = eps_e_tr_principal_vec[:,1]
  n3 = eps_e_tr_principal_vec[:,2]

  # [2] Compute trial stress
  sigma_tr_principal_mag = np.inner(Ce_principal, eps_e_tr_principal_mag)

  sigma_tr1 = sigma_tr_principal_mag[0]
  sigma_tr2 = sigma_tr_principal_mag[1]
  sigma_tr3 = sigma_tr_principal_mag[2]

  sigma_tr = sigma_tr1*np.tensordot(n1,n1,axes=0) + sigma_tr2*np.tensordot(n2,n2,axes=0) + sigma_tr3*np.tensordot(n3,n3,axes=0)

  # [3] Check yielding
  f = Plasticity.f(sigma_tr1, sigma_tr2, sigma_tr3, lamda)

  # [3.1] If f <= 0, elastic.
  if f <= 0:
    print(">> Elastic!")

    # Update stress & strain
    sigma = sigma_tr

    eps_e = eps_e_tr
    eps   = eps_e + eps_p

  # [3.2] If f > 0, plastic.
  else:
    print(">> Plastic!")

    # Initialize variables
    eps_e_principal_mag, eps_e_principal_vec = np.linalg.eig(eps_e)

    eps_e1 = eps_e_principal_mag[0]
    eps_e2 = eps_e_principal_mag[1]
    eps_e3 = eps_e_principal_mag[2]
    dlamda  = 0

    x = np.zeros(4) # target variables
    x[0] = eps_e1
    x[1] = eps_e2
    x[2] = eps_e3
    x[3] = dlamda

    # Newton-Raphson iteration (return mapping)
    for ii in range(maxiter):

      # Initialize residual and jacobian
      res = np.zeros(4)
      jac = np.zeros((4,4))

      # Current strain
      eps_e1_current = x[0]
      eps_e2_current = x[1]
      eps_e3_current = x[2]

      # Current stress
      sigma1_current = a*eps_e1_current + b*eps_e2_current + b*eps_e3_current
      sigma2_current = b*eps_e1_current + a*eps_e2_current + b*eps_e3_current
      sigma3_current = b*eps_e1_current + b*eps_e2_current + a*eps_e3_current

      sigma1_current = sigma1_current
      sigma2_current = sigma2_current
      sigma3_current = sigma3_current

      # Current lamda
      lamda_current = lamda + x[3]

      # Update derivatives
      # >> First order derivatives
      dfdsig1, dfdsig2, dfdsig3 = Plasticity.df(sigma1_current, sigma2_current, sigma3_current)

      # >> Second order derivatives
      d2fdsig1dsig1, d2fdsig2dsig2, d2fdsig3dsig3, d2fdsig1dsig2, d2fdsig2dsig3, d2fdsig3dsig1 \
        = Plasticity.df2(sigma1_current, sigma2_current, sigma3_current)

      # Update residual
      res[0] = x[0] - eps_e_tr1 + x[3]*dfdsig1
      res[1] = x[1] - eps_e_tr2 + x[3]*dfdsig2
      res[2] = x[2] - eps_e_tr3 + x[3]*dfdsig3
      res[3] = Plasticity.f(sigma1_current, sigma2_current, sigma3_current, lamda_current)

      # Update Jacobian ***
      jac[0,0] = 1 + x[3]*(d2fdsig1dsig1*dsig1depse1 + d2fdsig1dsig2*dsig2depse1 + d2fdsig3dsig1*dsig3depse1)
      jac[0,1] =     x[3]*(d2fdsig1dsig1*dsig1depse2 + d2fdsig1dsig2*dsig2depse2 + d2fdsig3dsig1*dsig3depse2)
      jac[0,2] =     x[3]*(d2fdsig1dsig1*dsig1depse3 + d2fdsig1dsig2*dsig2depse3 + d2fdsig3dsig1*dsig3depse3)
      jac[0,3] = dfdsig1

      jac[1,0] =     x[3]*(d2fdsig1dsig2*dsig1depse1 + d2fdsig2dsig2*dsig2depse1 + d2fdsig2dsig3*dsig3depse1)
      jac[1,1] = 1 + x[3]*(d2fdsig1dsig2*dsig1depse2 + d2fdsig2dsig2*dsig2depse2 + d2fdsig2dsig3*dsig3depse2)
      jac[1,2] =     x[3]*(d2fdsig1dsig2*dsig1depse3 + d2fdsig2dsig2*dsig2depse3 + d2fdsig2dsig3*dsig3depse3)
      jac[1,3] = dfdsig2

      jac[2,0] =     x[3]*(d2fdsig3dsig1*dsig1depse1 + d2fdsig2dsig3*dsig2depse1 + d2fdsig3dsig3*dsig3depse1)
      jac[2,1] =     x[3]*(d2fdsig3dsig1*dsig1depse2 + d2fdsig2dsig3*dsig2depse2 + d2fdsig3dsig3*dsig3depse2)
      jac[2,2] = 1 + x[3]*(d2fdsig3dsig1*dsig1depse3 + d2fdsig2dsig3*dsig2depse3 + d2fdsig3dsig3*dsig3depse3)
      jac[2,3] = dfdsig3

      jac[3,0] = dfdsig1*dsig1depse1 + dfdsig2*dsig2depse1 + dfdsig3*dsig3depse1
      jac[3,1] = dfdsig1*dsig1depse2 + dfdsig2*dsig2depse2 + dfdsig3*dsig3depse2
      jac[3,2] = dfdsig1*dsig1depse3 + dfdsig2*dsig2depse3 + dfdsig3*dsig3depse3
      jac[3,3] = 0

      # Solve system of equations
      dx = np.linalg.solve(jac, -res) # increment of target variables

      # Update x
      x = x + dx

      # Compute error
      err = np.linalg.norm(dx)

      print(" Newton iter.",ii, ": err =", err)

      if err < tol:
        break

    # Update strain
    eps   = eps + deps
    eps_e = x[0]*np.tensordot(n1,n1,axes=0) + x[1]*np.tensordot(n2,n2,axes=0) + x[2]*np.tensordot(n3,n3,axes=0)
    eps_p = eps - eps_e
    lamda = lamda + x[3]

    # Update stress
    sigma1 = a*x[0] + b*x[1] + b*x[2]
    sigma2 = b*x[0] + a*x[1] + b*x[2]
    sigma3 = b*x[0] + b*x[1] + a*x[2]
    sigma  = sigma1*np.tensordot(n1,n1,axes=0) + sigma2*np.tensordot(n2,n2,axes=0) + sigma3*np.tensordot(n3,n3,axes=0)
################### delete ###################
    aAB = np.zeros((3,3))
    mab = np.zeros((3,3,3,3))
    dsigdepse = np.zeros((3,3))
    dsigdepse[0,0] = dsig1depse1
    dsigdepse[0,1] = dsig1depse2
    dsigdepse[0,2] = dsig1depse3
    dsigdepse[1,0] = dsig2depse1
    dsigdepse[1,1] = dsig2depse2
    dsigdepse[1,2] = dsig2depse3
    dsigdepse[2,0] = dsig3depse1
    dsigdepse[2,1] = dsig3depse2
    dsigdepse[2,2] = dsig3depse3

    aAB = np.linalg.solve(jac[0:3,0:3],dsigdepse)

    for ii in range(3):
      for jj in range(3):
        na = eps_e_principal_vec[:,ii].reshape(-1,1)
        nb = eps_e_principal_vec[:,jj].reshape(-1,1)
        mab[ii,jj,:,:] = na @ nb.T
################### delete ###################
  # [4] Record stress and strain
  sigma11[i+1] = sigma[0,0]
  eps11[i+1]   = eps[0,0]

# Plot stress-strain curve
plt.figure(0,figsize=(7,7))
plt.plot(eps11, sigma11, 'r-', linewidth=1.0, label="NN prediction")
plt.xlabel(r'$\epsilon_{11}$', fontsize=15)
plt.ylabel(r'$\sigma_{11}$ [MPa]', fontsize=15)
plt.axhline(0, color = 'k',alpha = 0.5)
plt.axvline(0, color = 'k',alpha = 0.5)
plt.xlim(0.0, 0.03)
#plt.ylim(-600, 600)
plt.show()
