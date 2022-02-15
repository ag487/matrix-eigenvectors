import numpy as np


def hessenberg(X):
    """
    converts matrix to Hessenberg form, with zeros below the first subdiagonal
    [ x x x x ]
    [ x x x x ]
    [ - x x x ]
    [ - - x x ]

    advantages:
    Hessenberg matrix has same eigenvalues (not eigenvectors) as the original matrix
    zeros make further calculations easier

    the function actually transposes the Hessenberg matrix & performs a second reflection
    resulting in a tridiagonal matrix
    [ x x - - ]
    [ x x x - ]
    [ - x x x ]
    [ - - x x ]

    X = Q1.H.Q1T for symmetric matrix X
    """
    m = X.shape[0]
    n = X.shape[1]
    assert m == n, "error, Householder reflection here only works for nxn matrix"
    identity = np.asarray([[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)])
    Q1 = np.copy(identity)
    H = np.copy(X)
    for j in range(n-2):
        vector = np.copy(H[j+1:, j])  # column vector below diagonal
        sign = 1 if vector[0] >= 0 else -1
        vector[0] += sign * pow(np.dot(vector.T, vector), 0.5)  # as basis vector is the norm
        reflect = -1 * np.copy(identity)
        reflect[j+1:, j+1:] += 2 * np.outer(vector, vector.T) / np.dot(vector.T, vector)
        H = np.matmul(np.matmul(H, reflect).T, reflect)
        Q1 = np.matmul(Q1, reflect)
    return Q1, H
