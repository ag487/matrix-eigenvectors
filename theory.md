## Overview
This is an explanation of the matrix calculations used in finding the eigenvalues & eigenvectors of a matrix. Throughout, I have assumed the input is a symmetric matrix (which is consistent with the image covariance matrix used in the case study).

## 1. Theory essentials

### 1.1 Matrix determinant
Whether a matrix is invertible is shown by the determinant. First definition, a matrix is "invertible" if it can be multiplied by another matrix to give the identity matrix ie, it has an inverse. Secondly, what is the determinant? Well, let's try to invert a matrix & see what is possible. So, take a 2x2 matrix [ a c | b d ] (formatting here as two column vectors) multiply by 2x2 matrix [ e g | f h ] = 2x2 matrix [ 1 0 | 0 1 ].

Given matrix multiplication a.e + b.g = 1, c.e + d.g = 0. Then solve for e and g.

- g = -c.e / d
- a.e + b.-c.e / d = 1  subs g in other equation
- e.(a.d - b.c) = d
- e = d / (a.d - b.c)
- g = -c / (a.d - b.c)

Given matrix multiplication a.f + b.h = 0, c.f + d.h = 1. Then solve for f and h.

- h = -a.f / b
- c.f + d.-a.f / b = 1  subs
- f.(c.b - a.d) = b
- f = -b / (a.d - b.c)
- h = a / (a.d - b.c)

Then inverse = 2x2 matrix [ d -c | -b a ] / (a.d - b.c). If the determinant a.d - b.c = 0, then the matrix does not have an inverse.

### 1.2 Matrix eigenvalues
Matrix M multiplied by vector u, gives another vector v. The eigenvector or characteristic vector is defined as a non-zero vector u, that gives vector v that is a scalar multiple of u. The scalar multiple is the eigenvalue for that eigenvector. Can then rearrange:

- M. u = eigval. u
- M. u - eigval. u = 0
- (M - eigval. identity). u = 0

If the matrix in brackets above (M - eigval. identity) were invertible, then it could be multiplied by its inverse to give the identity. This would result in 1.u = 0, which contradicts our definition that u is non-zero. Thus the matrix cannot be invertible, meaning its determinant = 0.

### 1.3 Matrix eigenvectors
As an aside, an eigenvector itself is not as important as the relative proportions of its components. For instance, if a matrix has the effect of doubling vector [ 1 1 ], then it will also double 3 times that vector. So there is more than one eigenvector for each eigenvalue, however their proportions are the same. Usually, the representative eigenvector chosen is that with norm 1.

### 1.4 Characteristic polynomial
Take as an example a 2x2 matrix M, with elements a b c d, and calculate the determinant of the above matrix (M - eigval. identity).

- starting matrix is 2x2 [ a c | b d ]
- subtract eigval from the diagonals a & d
- then determinant of this matrix, which was a.d - b.c becomes:
- (a - eigval).(d - eigval) - b.c = 0
- eigval^2^ - (a + d).eigval + (a.d - b.c) = 0

This is the characteristic polynomial of matrix M. The eigenvalues are the zeros of the polynomial. Helpfully, the product of the eigenvalues equals the constant term a.d - b.c, and the sum of the eigenvalues equals the trace of the matrix a + d. This is due to the inherent symmetry of polynomials, for instance (x - 2).(x - 3) = x^2^ - 5.x + 6. The coefficients are the sum and the product of -2 and -3.

Given the above is a 2x2 matrix, the characteristic polynomial is of degree 2, so the quadratic formula gives the solution. For a standard quadratic a.x^2^ + b.x + c = 0 the formula is (-b +- (b^2^ - 4.a.c) ^1/2^) / 2.a. Substituting the characteristic polynomial terms in the quadratic formula:

- (a + d +- ((a + d) ^2^ - 4.(a.d - b.c)) ^1/2^) / 2
- = (a + d)/2 +- (a^2^ + 2.a.d + d^2^ - 4.a.d + 4.b.c) ^1/2^ / 2
- = (a + d)/2 +- (a^2^ - 2.a.d + d^2^ + 4.b.c) ^1/2^ / 2
- = (a + d)/2 +- ((a - d) ^2^ + 4.b.c) ^1/2^ / 2
- = (a + d)/2 +- (((a - d) /2) ^2^ + b.c) ^1/2^

Thie gives an equation for the two eigenvalues, and is used in the Wilkinson shift section below.

## 2. Utilities folder, Hessenberg matrix function
The goal is to convert the symmetric matrix (X) to Hessenberg form (H), that is, a matrix with zeros below the subdiagonal. The advantage is a Hessenberg matrix has the same eigenvalues (not eigenvectors) as the original matrix, and the zeros make further calculations easier (the original eigenvectors can be recovered with the Q matrix in X = QT.H.Q).

Unfortunately, the matrix cannot be taken directly to diagonal form by Householder reflections. The reason is that while the reflection affects the first column, it does not affect the first row, if the subdiagonal is chosen. However, it does affect the first row if the diagonal is chosen. For the LHS multiplication this is fine, but when applied on the RHS this has the effect of disturbing the reflection applied to the first column. Hence, limit to the subdiagonal.

### 2.1 Householder reflection
For a Householder reflection, the portion of each column below the diagonal is to be reflected onto a multiple of the basis vector. For example, a 4x4 matrix would start with rows 2-4 of column 1. This gives a vector x = [ x21 x31 x41 ] using row column indices. The objective is to reflect this onto a multiple of the basis vector [ 1 0 0 ].

Given the symmetry is reflection, both vectors have the same length. The multiple of the basis vector must then be equal to the length of vector x. For instance, length or norm of vector [ 2 4 4 ] = (2^2^ + 2. 4^2^) ^1/2^ = 6. Then the reflected vector must be 6 times the basis vector = [ 6 0 0 ]. In summary:

- vector x is known
- basis vector times multiple is known = norm x. b
- these two vectors are reflected across a hyperplane
- to find the hyperplane, add the two column vectors = [ 8 4 4 ]

### 2.2 Portions parallel & perpendicular
For any two vectors x & v, vector x can be split into the portion in the same direction as v, as well as the portion not in the same direction (ie, perpendicular). Define x = c.v + d.

The portion in same direction: c.v is the portion, c the multiple, so the objective is to find c. Dotted x with v, gives xT.v = c. vT.v + dT.v. Given that dT.v = 0 given orthogonal, then rearranges to c = xT.v / vT.v.

The portion orthogonal: find d, where d = x - c.v. Then d = x - xT.v. v / vT.v.

### 2.3 Reflection matrix
To reflect each column of the matrix across the hyperplane, keep the portion in the same direction as the hyperplane, but reverse the portion orthogonal. Then c.v - (x - c.v) = 2.c.v - x = 2. xT.v.v / vT.v - x. This simplifies to the reflection matrix:

- = 2.v. xT.v / vT.v - x   as scalar. vector = vector. scalar
- = 2.v. vT.x / vT.v - x   as dot product is commutative
- = (2.v.vT / vT.v - identity). x

## 3. Utilities folder, QR decomposition function

### 3.1 Givens rotation function
With the Hessenberg matrix as input, Givens rotations can decompose the matrix into an orthonormal matrix times an upper triangular matrix.

The objective is to rotate the column vector [ a b ] from the diagonal and subdiagonal components, onto a column vector [ r 0 ]. Given the symmetry is rotation, the vectors have the same norm, so r^2^ = a^2^ + b^2^. This brings to mind trigonometry. The angle theta between the vectors has cos = a/r and sin = -b/r. The matrix that rotates a b onto the basis vector is then:

- 2x2 matrix [ cos sin | -sin cos] times 2x1 vector [ a b ] 
- given a^2^/r + b^2^/r = r 
- and -a.b/r + a.b/r = 0
- = 2x1 vector [ r 0 ]

### 3.2 Wilkinson shift function
While the Givens rotations can be applied directly to the Hessenberg matrix, convergence would be slower. To speed convergence the Wilkinson shift essentially subtracts an eigenvalue from the matrix diagonal. The difference is a reduction from about 1250+ iterations without shift, to about 500+ iterations with shift for the example provided.

The eigenvalue formula from the theory section above is used in the denominator of the Wilkinson shift calculation eig = (a + d)/2 +- (((a - d) /2) ^2^ + b.c) ^1/2^. To check the calculation is consistent, rather than re-deriving the formula. As the matrix is symmetric, b = c so 2x2 [ a c | c d ].

- w = d - c^2^ / ((a - d)/2 + (((a - d) /2) ^2^ + c^2^) ^1/2^)
- w = d - c^2^ / (eig - d)   subs from theory section
- w. eig - w. d = d. eig - d^2^ - c^2^
- w. eig = d.(eig + w) - d^2^ - c^2^
- suppose w is eigenvalue, we know trace = eig1 + eig2 = a + d
- eig1. eig2 = d.a + d^2^ - d^2^ - c^2^
- eig1. eig2 = a.d - c^2^

LHS is the product of eigenvalues, RHS is the determinant of the underlying matrix. From the theory section, these are equal. Check complete. Wilkinson shift gives one of the eigenvalues.

### 3.3 Power iteration
Ignoring subtleties of the decomposition function calculations for a moment:

- givens function takes input M & returns Q R where M = Q.R
- so R = QT.M
- decomposition function recombines M2 as R.Q = (QT.M).Q
- iteration continues eg M4 = Q3T.Q2T.QT.M.Q.Q2.Q3

Repeated iteration essentially takes the matrix in the direction of the eigenvectors, until an arbitrary stopping point when c of the Wilkinson shift matrix is small enough.

ENDS