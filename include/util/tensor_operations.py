#import autograd.numpy as np
import numpy as np

def FthTensor2Voigt(Ce):

  ddsdde = np.zeros((6,6))
  for ii in range(3):
    for jj in range(3):
      ddsdde[ii,jj] = Ce[ii,ii,jj,jj]

  ddsdde[0,3] = Ce[0,0,0,1]
  ddsdde[0,4] = Ce[0,0,1,2]
  ddsdde[0,5] = Ce[0,0,2,0]

  ddsdde[1,3] = Ce[1,1,0,1]
  ddsdde[1,4] = Ce[1,1,1,2]
  ddsdde[1,5] = Ce[1,1,2,0]

  ddsdde[2,3] = Ce[2,2,0,1]
  ddsdde[2,4] = Ce[2,2,1,2]
  ddsdde[2,5] = Ce[2,2,2,0]

  ddsdde[3,0] = Ce[0,1,0,0]
  ddsdde[3,1] = Ce[0,1,1,1]
  ddsdde[3,2] = Ce[0,1,2,2]

  ddsdde[4,0] = Ce[1,2,0,0]
  ddsdde[4,1] = Ce[1,2,1,1]
  ddsdde[4,2] = Ce[1,2,2,2]

  ddsdde[5,0] = Ce[2,0,0,0]
  ddsdde[5,1] = Ce[2,0,1,1]
  ddsdde[5,2] = Ce[2,0,2,2]

  ddsdde[3,3] = Ce[0,1,0,1]
  ddsdde[3,4] = Ce[0,1,1,2]
  ddsdde[3,5] = Ce[0,1,2,0]

  ddsdde[4,3] = Ce[1,2,0,1]
  ddsdde[4,4] = Ce[1,2,1,2]
  ddsdde[4,5] = Ce[1,2,2,0]

  ddsdde[5,3] = Ce[2,0,0,1]
  ddsdde[5,4] = Ce[2,0,1,2]
  ddsdde[5,5] = Ce[2,0,2,0]
  return ddsdde

# tensor outer product: C_ijkl = A_ij B_kl
def tensor_oMult(A, B):
  # A, B: 2nd tensors
  # return: 4th tensor
  assert(A.shape == B.shape)
  nDim = A.shape[0]
  res = np.zeros((nDim,nDim,nDim,nDim))
  for i in range(nDim):
    for j in range(nDim):
      for k in range(nDim):
        for l in range(nDim):
          res[i,j,k,l] = A[i,j] * B[k,l]
  return res

# tensor oPlus operation: C_ijkl = A_jl B_ik
def tensor_oPlus(A, B):
  # A, B: 2nd tensors
  # return: 4th tensor
  assert(A.shape == B.shape)
  nDim = A.shape[0]
  res = np.zeros((nDim,nDim,nDim,nDim))
  for i in range(nDim):
    for j in range(nDim):
      for k in range(nDim):
        for l in range(nDim):
          res[i,j,k,l] = A[j,l] * B[i,k]
  return res

# tensor oMinus operation: C_ijkl = A_il B_jk
def tensor_oMinus(A, B):
  # A, B: 2nd tensors
  # return: 4th tensor
  assert(A.shape == B.shape)
  nDim = A.shape[0]
  res = np.zeros((nDim,nDim,nDim,nDim))
  for i in range(nDim):
    for j in range(nDim):
      for k in range(nDim):
        for l in range(nDim):
          res[i,j,k,l] = A[i,l] * B[j,k]
  return res

# compute the 4th order identity tensor II
# such that for any symmetric 2nd order tensor A_ij, II_ijkl A_kl = A_ij
def identity_4(nDim):
  I = np.eye(nDim)
  res = np.zeros((nDim,nDim,nDim,nDim))
  for i in range(nDim):
    for j in range(nDim):
      for k in range(nDim):
        for l in range(nDim):
          res[i,j,k,l] = (I[i,l] * I[j,k] + I[i,k] * I[j,l]) / 2.
  return res
