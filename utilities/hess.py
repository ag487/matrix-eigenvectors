import numpy as np

def hessenberg(X):
    """ Converts matrix to Hessenberg form, with zeros below the first subdiagonal,
    using Householder reflections. Given the input is symmetric, actually returns 
    a tridiagonal matrix. The advantage is that a Hessenberg matrix has the same 
    eigenvalues as the original matrix, and zeros make further calculations easier.
    
    [ x x - - ]
    [ x x x - ]
    [ - x x x ]
    [ - - x x ]   

    X = Q1.H.Q1T 
    """

    dims = X.shape
    assert dims[0] == dims[1], "error, Householder reflection here only works for nxn matrix"
    
    identity = np.eye(dims[1], dtype="float64")
    Q1 = np.copy(identity)
    H = np.copy(X).astype("float64")

    for j in range(dims[1] - 2):
        vector = H[j+1:, j]  # column vector below diagonal

        sign = 1 if vector[0] >= 0 else -1
        hyper = np.copy(vector)
        hyper[0] += sign * pow(np.dot(vector.T, vector), 0.5)  # represents the reflection hyperplane

        # reflection matrix is less than nxn
        reflect = 2 * np.outer(hyper, hyper.T) / np.dot(hyper.T, hyper) - identity[j+1:, j+1:]
        
        H[j+1:, :] = np.matmul(reflect, H[j+1:, :])  # faster than using full nxn multiplication
        H[:, j+1:] = np.matmul(H[:, j+1:], reflect.T)
        Q1[:, j+1:] = np.matmul(Q1[:, j+1:], reflect.T)

    print("converted to Hessenberg form")
    return Q1, H