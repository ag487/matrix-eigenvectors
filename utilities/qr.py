import numpy as np

def givens_rotation(M):
    """ Converts Hessenberg matrix to upper triangular matrix,
    using Givens rotations.
    
    [ x x x x ]
    [ - x x x ]
    [ - - x x ]
    [ - - - x ]

    M = Q.R
    """

    m = M.shape[0]
    Q = np.eye(m, dtype="float64")
    R = np.copy(M).astype("float64")

    for j in range(m - 1):
        vector = R[j:j+2, j]   # the two diagonal and subdiagonal elements of a column

        norm = pow(np.dot(vector.T, vector), 0.5)
        cos = vector[0] / norm
        sin = vector[1] / norm
        rotate = np.array([[cos, sin], [-sin, cos]], dtype="float64")   # 2x2 clockwise rotation matrix

        R[j:j+2, :] = np.matmul(rotate, R[j:j+2, :])   # faster than using full nxn multiplication
        Q[:, j:j+2] = np.matmul(Q[:, j:j+2], rotate.T)

    return Q, R


def wilkinson_shift(M):
    """ Finds an eigenvalue of 2x2 lower right block
    [a c]
    [c d]
    """

    a = M[0, 0]
    c = M[1, 0]
    d = M[1, 1]
    diff = (a - d) / 2
    sign = 1 if diff >= 0 else -1
    w = d - sign * pow(c, 2) / (abs(diff) + pow(pow(diff, 2) + pow(c, 2), 0.5))
    return w


def decomposition(H):
    """ Decompose Hessenberg matrix into its eigenvalues & eigenvectors
    using repeated Givens rotations, with Wilkinson shift to speed convergence.

    H = Q2.right.Q2T
    """
    
    n = H.shape[0]
    identity = np.eye(n, dtype="float64")
    Q2 = np.copy(identity)
    right = np.copy(H).astype("float64")

    count = 0
    for k in range(n, 1, -1):  # 512 .. 2
        while True:
            count += 1

            w = wilkinson_shift(right[k-2:k, k-2:k])   # returns scalar from lower right 2x2 matrix
            right_shifted = right[:k, :k] - w * identity[:k, :k]

            Q, R = givens_rotation(right_shifted)
            Q2[:, :k] = np.matmul(Q2[:, :k], Q)
            right[:k, :k] = np.matmul(R, Q) + w * identity[:k, :k]

            if abs(right[k-1, k-2]) < 1e-4:   # arbitrary stopping value for c in Wilkinson sub matrix
                break

        if (513 - k) % 100 == 0:
            print(f"qr iter: {count}..")

    eigvals = [right[j, j] for j in range(n)]
    return Q2, eigvals
