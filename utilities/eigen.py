'''
Copyright (c) 2021 Angus Graham
All rights reserved.

This source code is licensed under the MIT License found in the
license.md file in the root directory of this source tree.
'''

import numpy as np

'''
qr_givens function
requires Hessenberg matrix as input
uses Givens rotations to decompose matrix into orthogonal matrix times upper
triangular matrix
[ x x x x ]
[ - x x x ]
[ - - x x ]
[ - - - x ]

H = Q.R
'''

def qr_givens(H):
  m = H.shape[0]
  n = H.shape[1]
  R = np.copy(H)
  Q = np.asarray([[ 1.0 if i == j else 0.0 for j in range(n) ] \
    for i in range(m) ])
  for i in range(n-1):
    u = R[i:i+2, i]
    cos = -u[0] / np.linalg.norm(u)  # to opposite axis, theta + pi for accuracy
    sin = -u[1] / np.linalg.norm(u)
    G = np.asarray([[cos, sin], [-sin, cos]])
    R[i:i+2, i:] = np.matmul(G, R[i:i+2, i:])
    Q[i:i+2, :i+2] = np.matmul(G, Q[i:i+2, :i+2])
  R[n-1, n-1] = -R[n-1, n-1]  # given opposite axis, change sign for last
  Q[n-1, :] = -Q[n-1, :]
  return Q.T, R

'''
wilkinson_shift function
finds an eigenvalue of 2x2 lower right block
[a c]
[c d]
'''

def wilkinson_shift(H):
  m = H.shape[0]
  a = H[m-2, m-2]
  c = H[m-1, m-2]
  d = H[m-1, m-1]
  diff = (a - d) / 2
  sign = 1 if diff >= 0 else -1
  w = d - sign * (c ** 2) / (abs(diff) + (diff ** 2 + c ** 2) ** 0.5)
  return w

'''
decomposition function
uses power iteration with Wilkinson shift
to decompose Hessenberg matrix into its eigenvalues & eigenvectors
H = Q.R.QT
'''

def decomposition(H):
  m = H.shape[0]
  n = H.shape[1]
  I = np.asarray([[ 1.0 if i == j else 0.0 for j in range(m) ] \
    for i in range(m) ])
  Q_product = np.copy(I)
  H_topleft = np.copy(H)
  H_full = np.zeros([n,n])
  turn = 0
  for k in range(n,1,-1):  # 512 to 2
    flag = True
    while flag:
      turn += 1
      w = wilkinson_shift(H_topleft)  # returns scalar
      M = H_topleft - w * I[:k, :k]
      Q, R = qr_givens(M)
      Q_full = np.copy(I)
      Q_full[:k, :k] = Q
      Q_product = np.matmul(Q_product, Q_full)
      H_topleft = np.matmul(R, Q) + w * I[:k, :k]
      if abs(H_topleft[k-1, k-2]) < 1e-4:
        flag = False
    if k % 50 == 0:
      print("iter:", turn)
    H_full[:k, :k] = H_topleft
    H_topleft = H_topleft[:k-1, :k-1]
  eigvals = np.asarray([ H_full[i,i] for i in range(n) ])
  return Q_product, eigvals
