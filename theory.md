This is an explanation of the matrix calculations used in finding the eigenvalues & eigenvectors of a matrix.

I have assumed the input is a symmetric matrix (which reflects the image covariance matrix used in the case study). This also helps with the initial Hessenberg matrix reduction to tridiagonal form.

### Theory basics

#### theory, matrix determinant
Whether a matrix is invertible or not depends on the determinant. What is the determinant? Well let's try to invert a matrix & find e f g h.
- [ a b ].[ e f ] = [ 1 0 ]
- [ c d ] [ g h ]   [ 0 1 ]

We know from matrix multiplication that a.e + b.g = 1, c.e + d.g = 0. Then solve for e g.
- g = -c.e / d
- a.e + b.-c.e / d = 1  subs g in other equation
- e.(a.d - b.c) = d
- e = d / (a.d - b.c)
- g = -c / (a.d - b.c)

We know from matrix multiplication that a.f + b.h = 0, c.f + d.h = 1. Then solve for f h.
- h = -a.f / b
- c.f + d.-a.f / b = 1  subs
- f.(c.b - a.d) = b
- f = -b / (a.d - b.c)
- h = a / (a.d - b.c)

then inverse is.
- [ d -b ] / (a.d - b.c)
- [-c  a ]

If the determinant a.d - b.c = 0, then the matrix does not have an inverse. This means the matrix is not full rank.

#### theory, matrix eigenvalues
Matrix M multiplied by vector u, gives another vector v. We define a characteristic vector or eigenvector as a non-zero vector u, that gives vector v which is a scalar multiple of u. The scalar multiple is the eigenvalue for that eigenvector.

can then rearrange.
- M. u = eigval. u
- M. u - eigval. u = 0
- (M - eigval. identity). u = 0

If the above matrix in brackets (M - eigval. identity) were invertible, we could multiply both sides of the equation by its inverse. This would result in u = 0, which contradicts our definition that u is non-zero. Therefore the matrix cannot be invertible, meaning its determinant = 0.

#### theory, matrix eigenvectors
As an aside, an eigenvector itself is not as important as the relative proportions of its components. For instance, if a matrix has the effect of doubling vector [1 1].T, then it will also double 3 times that vector. So there is more than one eigenvector for each eigenvalue, however their proportions are the same.

#### theory, characteristic polynomial
Take as an example a 2x2 matrix M, with elements a b c d, and calculate the determinant of the above matrix.
- starting matrix is
- [ a  b ]
- [ c  d ]
- subtract eigval from the diagonals a & d
- then determinant of this matrix is
- (a - eigval).(d - eigval) - b.c = 0
- eigval^2 - (a + d).eigval + (a.d - b.c) = 0

This is the characteristic polynomial of matrix M, as the eigenvalues are the zeros of the function. We can see from the polynomial that, helpfully, the product of the eigenvalues equals the constant term of a.d - b.c, and the sum of the eigenvalues the trace of the matrix a + d. For instance, if we could factor to (eigval - p).(eigval - q) = eigval^2 - (p + q).eigval + p.q. Then p q would equal the eigenvalues where both brackets resolve to zero.

Given the 2x2 matrix, the characteristic polynomial is of degree 2, so we can use the quadratic formula to solve. If we had a standard quadratic a.x^2 + b.x + c = 0 then the formula would be (-b +- (b^2 - 4.a.c) ^1/2) / 2.a

Substituting the terms we have gives:
- (a + d +- ((a + d) ^2 - 4.(a.d - b.c)) ^1/2) / 2
- (a + d)/2 +- (a^2 + 2.a.d + d^2 - 4.a.d + 4.b.c) ^1/2 / 2
- (a + d)/2 +- (a^2 - 2.a.d + d^2 + 4.b.c) ^1/2 / 2
- (a + d)/2 +- ((a - d) ^2 + 4.b.c) ^1/2 / 2
- (a + d)/2 +- (((a - d) /2) ^2 + b.c) ^1/2
- which gives an equation for the two eigenvalues
- this is used in the Wilkinson shift section below

### Householder reflection

#### find hyperplane
We use Householder reflections to convert the symmetric matrix to Hessenberg form.

We want to take the portion of each column below the diagonal, and reflect it onto a multiple of a basis vector. Take for example a 4x4 matrix. We would start with rows 2-4 of column 1. This gives us a vector x = [ x21 x31 x41 ] using row column indices. We want to reflect this onto a multiple of the basis vector [ 1 0 0 ].

Given we are using the symmetry of reflection, we know both vectors have the same length. So we know the multiple is the norm of vector x. For instance, norm of vector [ 2 4 4 ] = (2^2 + 2. 4^2) ^1/2 = 6. In this example, we know the reflected vector is 6 times the basis vector = [ 6 0 0 ].

We can visualise in terms of geometry:
- vector x is known
- basis vector times multiple is known = norm x. b
- these two vectors are reflected across a hyperplane
- to find the hyperplane, just add the two vectors = [ 8 4 4 ]

#### split parallel & perpendicular
If we take any two vectors x & v, we can find the portion of vector x in the same direction as v, as well as the portion not in the same direction (ie orthogonal or perpendicular) so define x = c.v + d.

portion in same direction.
- c.v is the portion, c the multiple
- find c
- dot with v, gives xT.v = c. vT.v + dT.v  
- where dT.v = 0 given orthogonal
- rearrange to  c = xT.v / vT.v

portion orthogonal.
- find d
- x - c.v
- x - xT.v. v / vT.v

#### apply reflection
To reflect each column of the matrix across the hyperplane, keep the portion in the same direction as the hyperplane, but reverse the portion orthogonal.
- c.v - (x - c.v)
- 2.c.v - x
- 2. xT.v.v / vT.v - x

we can simplify.
- 2.v. xT.v / vT.v - x  rearrange as scalar. vector = vector. scalar
- 2.v. vT.x / vT.v - x  dot product is commutative
- (2.v.vT / vT.v - identity). x

### Givens rotations
Once we have the Hessenberg form of the matrix, we can apply Givens rotations to decompose the matrix into an orthonormal matrix times an upper triangular matrix.

The objective is to rotate the column vector [a b] onto a column vector [r 0]. The vectors should have the same norm, so r^2 = a^2 + b^2.

This brings to mind trigonometry. The angle theta from the first to second vector has cos = a/r and sin = -b/r. The matrix that would then rotate a b onto the basis vector is.
- [ cos -sin ].[ a ] = [ a^2/r + b^2/r ] = [ r ]
- [ sin  cos ] [ b ]   [-a.b/r + a.b/r ]   [ 0 ]

### Wilkinson shift
While we could apply the Givens rotations directly to the Hessenberg form of the matrix, convergence would be slow. To speed convergence we can apply a Wilkinson shift, essentially subtracting an eigenvalue from the lower right corner of the matrix.

The eigenvalue formula from the theory section above is used in the denominator of the Wilkinson shift calculation. We will check the calculation is consistent, rather than re-deriving the formula.
- matrix is symmetric so b = c
- [ a  c ]
- [ c  d ]

- eig2 = d - c^2 / ((a - d)/2 + (((a - d) /2) ^2 + c^2) ^1/2)
- eig2 = d - c^2 / (eig1 - d)
- eig2. eig1 - eig2. d = d. eig1 - d^2 - c^2
- eig2. eig1 = d.(eig1 + eig2) - d^2 - c^2
- we know trace = eig1 + eig2 = a + d
- eig2. eig1 = d.a + d^2 - d^2 - c^2
- eig2. eig1 = a.d - c^2

LHS is the product of eigenvalues, RHS is the determinant of the underlying matrix. From the theory section, we know these are equal. So check complete. Wilkinson shift gives us one of the eigenvalues.

### Power iteration
Ignoring subtleties of the decomposition function calculations for a moment.
- qr_givens function takes input M & returns Q R where M = Q.R
- so R = QT.M
- decomposition function recombines M2 as R.Q = (QT.M).Q
- iteration continues eg M4 = Q3T.Q2T.QT.M.Q.Q2.Q3

Repeated iteration essentially takes the matrix in the direction of the eigenvectors. We set an arbitrary stopping point when c of the Wilkinson shift matrix is small enough.

ENDS
