'''
Copyright (c) 2021 Angus Graham
All rights reserved.

This source code is licensed under the MIT License found in the
license.md file in the root directory of this source tree.
'''

import numpy as np

'''
hessenberg function
converts matrix to Hessenberg form, with zeros below the first subdiagonal
[ x x x x ]
[ x x x x ]
[ - x x x ]
[ - - x x ]

advantages
- Hessenberg matrix has same eigenvalues (not eigenvectors) as the original matrix
- zeros make further calculations easier

X = Q.H.QT for symmetric matrix
'''

def hessenberg(X):
  m = X.shape[0]
  n = X.shape[1]
  assert m == n, "error, householder reflection here only works for nxn matrix"
  I = np.asarray([[ 1.0 if i == j else 0.0 for j in range(n) ] \
    for i in range(m) ])
  Q = np.copy(I)
  R = np.copy(X)
  for i in range(n-2):
    v = np.copy(R[i+1:, i])  # column vector below diagonal
    sign = 1 if v[0] >= 0 else -1
    v[0] += sign * np.linalg.norm(v)  # as basis vector only affects first component
    H = np.copy(I)
    H[i+1:, i+1:] -= 2 * np.outer(v,v.T) / np.dot(v.T,v)
    R = np.matmul(np.matmul(H, R), H.T)
    Q = np.matmul(H, Q)
  return Q.T, R
