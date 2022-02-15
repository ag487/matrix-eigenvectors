import numpy as np


def qr_givens(M):
    """
    requires Hessenberg matrix as input
    uses Givens rotations to decompose matrix into upper triangular matrix
    [ x x x x ]
    [ - x x x ]
    [ - - x x ]
    [ - - - x ]
    M = Q.R
    """
    n = M.shape[0]
    R = np.copy(M)
    QT = np.asarray([[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)])
    for j in range(n-1):
        vector = R[j:j+2, j]
        norm = pow(np.dot(vector.T, vector), 0.5)
        cos = -vector[0] / norm  # to opposite axis, theta + pi for improved accuracy
        sin = -vector[1] / norm
        rotate = np.asarray([[cos, sin], [-sin, cos]])
        R[j:j+2, j:] = np.matmul(rotate, R[j:j+2, j:])
        QT[j:j+2, :j+2] = np.matmul(rotate, QT[j:j+2, :j+2])
    R[n-1, n-1] = -R[n-1, n-1]  # given opposite axis, change sign for last
    QT[n-1, :] = -QT[n-1, :]
    return QT.T, R


def wilkinson_shift(M):
    """
    finds an eigenvalue of 2x2 lower right block
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
    """
    decompose Hessenberg matrix into its eigenvalues & eigenvectors
    using repeated Givens rotations, with Wilkinson shift
    H = Q2.right.Q2T
    """
    n = H.shape[0]
    identity = np.asarray([[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)])
    Q2 = np.copy(identity)
    right = np.copy(H)
    turn = 0
    for k in range(n, 1, -1):  # 512 to 2
        flag = True
        while flag:
            turn += 1
            w = wilkinson_shift(right[k-2:k, k-2:k])  # returns scalar from 2x2 matrix block
            right_shifted = right[:k, :k] - w * identity[:k, :k]
            Q, R = qr_givens(right_shifted)
            Q2[:, :k] = np.matmul(Q2[:, :k], Q)  # same as overlaying Identity with Q, then multiplying full Q2
            right[:k, :k] = np.matmul(R, Q) + w * identity[:k, :k]
            if abs(right[k-1, k-2]) < 1e-4:  # the c in wilkinson shift sub matrix
                flag = False
        if k % 50 == 0:
            print("iter:", turn)
    eigvals = np.asarray([right[j, j] for j in range(n)])
    return Q2, eigvals
