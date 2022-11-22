import numpy as np


def hessenberg(X):
    """ Converts matrix to Hessenberg form, with zeros below the first subdiagonal,
    using Householder reflections.
    [ x x x x ]
    [ x x x x ]
    [ - x x x ]
    [ - - x x ]

    Given the input is symmetrical, actually returns a tridiagonal matrix.
    [ x x - - ]
    [ x x x - ]
    [ - x x x ]
    [ - - x x ]

    X = Q1T.H.Q1

    The advantage is that a Hessenberg matrix has the same eigenvalues as the original matrix, 
    and zeros make further calculations easier.
    """

    dims = X.shape
    assert dims[0] == dims[1], "error, Householder reflection here only works for nxn matrix"
    
    identity = np.asarray([[1.0 if i == j else 0.0 for j in range(dims[1])] for i in range(dims[1])])
    Q1 = np.copy(identity)
    H = np.copy(X).astype("float64")

    for j in range(dims[1] - 2):   # in columns
        vector = np.copy(H[j+1:, j])   # column vector below diagonal
        sign = 1 if vector[0] >= 0 else -1
        vector[0] += sign * pow(np.dot(vector.T, vector), 0.5)   # represents the reflection hyperplane

        reflect = 2 * np.outer(vector, vector.T) / np.dot(vector.T, vector) - identity[j+1:, j+1:]
        # this is the reflection matrix

        H[j+1:, :] = np.matmul(reflect, H[j+1:, :])   # faster than using full nxn equivalent
        H[:, j+1:] = np.matmul(H[:, j+1:], reflect) 
        Q1[j+1:, :] = np.matmul(reflect, Q1[j+1:, :])

    print("converted to Hessenberg form..")
    return Q1.T, H