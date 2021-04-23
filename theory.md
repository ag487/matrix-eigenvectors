This is an explanation of the matrix calculations used in finding the
eigenvalues & eigenvectors of a matrix.

# Symmetric matrix input

I have assumed the input is a symmetric matrix (which reflects the covariance
matrix used in the case study). The advantage is the eigenvectors of a
symmetric matrix are orthogonal, well actually orthonormal as length is 1.
Refer spectral theorem for proof.
- for non-symmetric matrix X = eigvecs. eigvals. inverse of eigvecs
- however given orthogonal, then eigvecs.eigvecsT = I so transpose is the inverse
- X = eigvecs. eigvals. transpose of eigvecs

# Householder reflection

We use this reflection to convert the matrix X to Hessenberg form.

- [ x11 x12 x13 x14 ]
- [ x21 x22 x23 x24 ]
- [ x31 x32 x33 x34 ]
- [ x41 x42 x43 x44 ]

We want to take the portion of each column below the diagonal, and reflect it
onto a multiple of a basis vector. So, for the first column we want to reflect
vector [ x21 x31 x41 ] onto a multiple of the basis vector [ 1 0 0 ].

Given we are using the symmetry of reflection, we know the vector and multiple
must have the same length, given basis vector has length 1. So we can calculate
multiple as the norm of the vector. For instance, norm of vector [3 4] is 5.

We can visualise in terms of geometry:
- vector x is known
- basis vector times multiple is known = norm x. b
- these two vectors are reflected across a hyperplane (representing matrix X)

Again, as we are using the symmetry of reflection, we know  x - norm x. b is
orthogonal to the hyperplane. So we can use the vector in this direction to
describe the hyperplane. Define vector v = x - norm x. b

## Theory, two vectors orthogonal portion

If we take two vectors x & v, we can find the portion of vector x in the same
direction as v, as well as the portion not in the same direction (ie orthogonal
or perpendicular).

portion in same direction.
- find x = a.v + b
- dot with v, gives x.v = a.v.v + b.v  
- where b.v = 0 as orthogonal to v
- rearrange as  a = x.v / v.v

portion orthogonal.
- x - a.v
- x - x.v / vT.v . v
- x - v.vT.x / vT.v

If instead we were to subtract double the amount from x, this would give a
reflection of x across the portion of x orthogonal to vector v.

## Theory, Householder formula

the theory is that.
- if v = x + norm x. b is non zero, then H(x) = - norm x. b
- else if v = x - norm x. b is non zero, then H(x) = norm x. b

proof for first case.
- from the above, we know the reflection H(x) = x - 2.v.vT.x / vT.v
- can we simplify this?

take vT.v
- (x + norm x. b)T. (x + norm x. b)
- xT.x + norm x. xT. b + norm x. bT. x + norm x ^2. bT. b
- norm x ^2 + norm x. x1 + norm x. x1 + norm x ^2
- = 2. norm x. (norm x + x1)

take v.vT.x
- (x + norm x. b). (x + norm x. b)T. x
- x.xT.x + norm x. x. bT. x + norm x. b. xT. x + norm x ^2. b. bT. x
- norm x ^2. x + norm x. x1. x + norm x ^3. b + norm x ^2. x1. b
- norm x. (norm x. x + x1. x + norm x ^2. b + norm x. x1. b)
- norm x. ((norm x + x1). x + norm x. (norm x + x1). b)
- = norm x. (norm x + x1). (x + norm x. b)

then H(x) simplified.
- x - 2.v.vT.x / vT.v
- x - 2. norm x. (norm x + x1). (x + norm x. b) / (2. norm x. (norm x + x1))
- x - (x + norm x. b)
- = - norm x. b  proven

# Givens rotations

Once we have the Hessenberg form of the matrix, we can apply Givens rotations
to decompose the matrix into an orthonormal matrix times an upper triangular
matrix.

- [ h11 h12 h13 h14 ]
- [ h21 h22 h23 h24 ]
- [ -   h32 h33 h34 ]
- [ -   -   h43 h44 ]

The objective is to rotate the column vector [a b] onto a column vector [r 0].
To preserve the eigenvalues, the column vectors should have the same norm. So
r^2 = a^2 + b^2.

This brings to mind trigonometry. The angle theta between a,b and the basis vector
has cos = a/r and sin = b/r. The matrix that would then rotate a,b onto the
basis vector, when multiplied on the left is.
- [ cos -sin ]
- [ sin cos ]

The other advantage of using cos and sin is that the matrix is orthonormal
as cos ^2 + sin ^2 = 1.

# Wilkinson shift

We could apply the Givens rotations directly to the Hessenberg form of the matrix.
However, convergence would be slow. To speed convergence we can apply a Wilkinson
shift to the lower right corner of the matrix, basically subtracting an eigenvalue.

- [a b]
- [c d]

## Theory, eigenvectors & eigenvalues

A matrix multiplied by vector, gives another vector  A. u = v
- by defn u is eigenvector if u has length 1 & v is scalar multiple of u
- eigenvalue is the scalar multiple

can then rearrange.
- A. u - eigenvalue. u = 0
- (A - eigenvalue. I). u = 0
- has non-zero solution for u if and only if the determinant is 0
- determinant here is the characterisic polynomial
- in geometry, determinant is the scalar area of the parallelogram (0,0) (a,b)
(a+c,b+d) (c,d) = a.d - b.c

can then calculate a.d - b.c
- (a - eig). (d - eig) - b. c
- a.d - a.eig - d.eig + eig ^ 2 - b.c
- eig ^ 2 - (a + d).eig + (a.d - b.c) = 0

use quadratic formula & substitute (-b +- (b.b - 4.a.c) ^ 0.5) / 2.a
- (a + d +- ((a + d) ^ 2 - 4.(a.d - b.c)) ^ 0.5) / 2
- (a + d)/2 +- (a.a + 2.a.d + d.d - 4.a.d + 4.b.c) ^ 0.5 / 2
- (a + d)/2 +- (a.a - 2.a.d + d.d + 4.b.c) ^ 0.5 / 2
- (a + d)/2 +- ((a - d) ^ 2 + 4.b.c) ^ 0.5 / 2
- (a + d)/2 +- (((a - d) / 2) ^ 2 + b.c) ^ 0.5
- gives the eigenvalues

## Theory, Wilkinson shift

The above formula is used in the denominator of the Wilkinson shift calculation.
- matrix is symmetric so b=c
- eig2 = d - sign.(c ^ 2) / (abs((a - d) /2) + (((a - d) /2) ^ 2 + c ^ 2) ^ 0.5)
- eig2 = d - sign.(c ^ 2) / (eig1 - d)
- eig2.eig1 - eig2.d = d.eig1 - d.d - sign.(c ^ 2)
- eig2.eig1 = d.(eig1 + eig2) - d.d - sign.(c ^ 2)
- we know trace = eig1 + eig2 = a + d
- eig2.eig1 = d.a + d.d - d.d - sign.(c ^ 2)
- eig2.eig1 = a.d - sign.(c ^ 2)
- we know RHS is determinant, so proven correct

# Power iteration

Ignoring subtleties of the calculations for a moment.
- qr_givens function takes input H1 & returns Q1,R1 where R1 = Q1T.H1
- decomposition function recombines H2 as R1.Q1 = Q1T.H1.Q1
- iteration continues eg H4 = Q3T.Q2T.Q1T.H1.Q1.Q2.Q3

Repeated iteration essentially takes the matrix in the direction of the
eigenvectors.

ENDS
