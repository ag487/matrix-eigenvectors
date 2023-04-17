## Overview

This is an explanation of the matrix calculations used in finding the eigenvalues & eigenvectors of a matrix. Throughout, I have assumed the input is a symmetric matrix (which is consistent with the image covariance matrix used in the case study).

> 1 Theory essentials
> 1.1 Matrix inverse and determinant
> 1.2 Matrix eigenvalues
> 1.3 Characteristic polynomial
>
> 2 Hessenberg matrix function
> 2.1 Householder reflection
> 2.2 Portions parallel & perpendicular
> 2.3 Reflection matrix
> 2.4 Reduced matrix calculations
> 
> 3 QR decomposition function
> 3.1 Givens rotation function
> 3.2 Reduced matrix calculations, part 2
> 3.3 Wilkinson shift function
> 3.4 Power iteration



## 1. Theory essentials

### 1.1 Matrix inverse and determinant

A matrix may or may not be invertible. A matrix is invertible if it can be multiplied by another matrix to give the identity matrix. That is, it has an inverse. For instance, let's take a 2x2 matrix and assume it has an inverse. Then multiplying the two matrices together would result in the identity matrix.

[ a b ].[ e f ] = [ a.e + b.g   a.f + b.h ]
[ c d ] [ g h ]   [ c.e + d.g   c.f + d.h ]

Given matrix multiplication a.e + b.g = 1, c.e + d.g = 0. Then solve for e and g.

> g = -c.e / d  from second eq
> a.e + b.-c.e / d = 1  subs g in first eq
> e.(a.d - b.c) = d
> e = d / (a.d - b.c)
> g = -c / (a.d - b.c)  subs e in first line

Given matrix multiplication a.f + b.h = 0, c.f + d.h = 1. Then solve for f and h.

> h = -a.f / b  from first eq
> c.f + d.-a.f / b = 1  subs h in second eq
> f.(c.b - a.d) = b
> f = -b / (a.d - b.c)
> h = a / (a.d - b.c)  subs f in first line

Then the inverse of the 2x2 matrix must be a specific permutation of the original matrix, divided by the scalar value of a.d - b.c.

[  d -b ] / (a.d - b.c)
[ -c  a ]

The value a.d - b.c is known as the determinant. If the determinant is zero, then the inverse of the matrix is undefined, in which case the original matrix is not invertible.

To push this point further, what if we started with a matrix where a.d - b.c = 0. Say 3.6 - 2.9 = 0. Then if we ignore the scalar for the moment, but multiply the matrix by the required permutation. In this example, we see the resulting null matrix cannot possibly be divided by a scalar value that would result in the identity matrix.

[ 3  2 ].[  6 -2 ] = [ 0 0 ]
[ 9  6 ] [ -9  3 ]   [ 0 0 ]


### 1.2 Matrix eigenvalues

A matrix multiplied by a vector results in another vector. Say a 3x3 matrix multiplied by a 3x1 vector results in a 3x1 vector. 

If the resulting vector is a scalar multiple of the initial vector, then the initial vector is known as an eigenvector of the original matrix. Let's define this to exclude zero vectors, as a zero vector takes any matrix to another zero vector, so gives a trivial result.

For example.

[ 2 -2  3 ].[ 1 ] = [ 3 ]
[ 5 -2  0 ] [ 1 ]   [ 3 ]
[ 1  3 -1 ] [ 1 ]   [ 3 ]

So an eigenvector here is the column vector 1, and the eigenvalue is the scalar 3. 

As an aside, an eigenvector itself is not as important as the relative proportions of its components. For instance, as the above matrix has the effect of tripling the column vector, it will also triple double that vector. So there is more than one eigenvector for each eigenvalue, however their proportions are the same. Usually, the representative eigenvector chosen is that with length 1.

Rearranging.

> matrix. vector = scalar. vector  (where the scalar in the above example is 3)
> matrix. vector - scalar. vector = 0
> (matrix - scalar. identity). vector = 0

[ -1 -2  3 ].[ 1 ] = [ 0 ]
[  5 -5  0 ] [ 1 ]   [ 0 ]
[  1  3 -4 ] [ 1 ]   [ 0 ]

If the reduced matrix in brackets above (matrix - scalar. identity) were invertible, then it could be multiplied by its inverse to give the identity. This would result in 1.vector = 0, which contradicts our definition that the vector is non-zero. Thus the reduced matrix cannot be invertible, meaning its determinant = 0.


### 1.3 Characteristic polynomial

Let's return to 2x2 matrices to make the calculations simpler, in trying to find the eigenvalues.

[ a b ]
[ c d ]

So then the reduced matrix of (matrix - eigval. identity) is.

[ a-eigval  b        ]
[ c 	    d-eigval ]

We know the determinant of a 2x2 matrix is a.d - b.c. So the determinant for the above is (a-eigval).(d-eigval) - b.c and we know this is zero. Rearranging then gives eigval^2 - (a + d).eigval + (a.d - b.c) = 0. This is the characteristic polynomial for this matrix. 

The eigenvalues are the zeros of the polynomial, so if we can solve the polynomial we will know the eigenvalues. In the case of a 2x2 matrix, the characteristic polynomial is of degree 2, so the quadratic formula gives the solution.

For a standard quadratic a.x^2 + b.x + c = 0 the formula is (-b +- (b^2 - 4.a.c) ^1/2) / 2.a. Substituting the characteristic polynomial terms in the quadratic formula:

> (a + d +- ((a + d) ^2 - 4.(a.d - b.c)) ^1/2) / 2
> = (a + d)/2 +- (a^2 + 2.a.d + d^2 - 4.a.d + 4.b.c) ^1/2 / 2
> = (a + d)/2 +- (a^2 - 2.a.d + d^2 + 4.b.c) ^1/2 / 2
> = (a + d)/2 +- ((a - d) ^2 + 4.b.c) ^1/2 / 2
> = (a + d)/2 +- (((a - d) /2) ^2 + b.c) ^1/2

This gives an equation for the two eigenvalues, and is used in the Wilkinson shift section. 

Also, from the symmetry of polynomial equations we know the product of the eigenvalues = a.d - b.c. For instance, taking a basic polynomial example below, we can see the constant term (which above = a.d - b.c and below = 6) must be equal to the product of the two values that satisfy the equation (which above are the two eigenvalues).

> (x - 2).(x - 3)
> = x^2 - (2 + 3).x + (2.3)
> = x^2 - 5.x + 6



## 2. Hessenberg matrix function

The first step is to convert the symmetric matrix (X) to Hessenberg form (H) with zeros below the subdiagonal. The advantage is a Hessenberg matrix has the same eigenvalues as the original matrix, and the zeros make further calculations easier (the eigenvectors can be recovered with the Q matrix in X = Q.H.QT).

Note the matrix cannot be taken directly to diagonal form by Householder reflections. The reason is that, if we aim for zeros below the subdiagonal, while the reflection affects the first column, it does not affect the first row. However, if we aim for zeros below the diagonal, it does affect the first row. For the LHS multiplication this is fine, but when applied on the RHS this has the effect of disturbing the reflection applied to the first column. Hence, limit to the subdiagonal.


### 2.1 Householder reflection

For a Householder reflection, the portion of each column below the diagonal is to be reflected onto a multiple of the basis vector. For example, a 4x4 matrix would start with rows 2-4 of column 1. This gives a vector x = [ x21 x31 x41 ] using row column indices. The objective is to reflect this onto a multiple of the basis vector [ 1 0 0 ].

Given the symmetry is reflection, both vectors have the same length. The multiple of the basis vector must then be equal to the length of vector x. For instance, length or norm of vector [ 2 4 4 ] = (2^2 + 2. 4^2) ^1/2 = 6. Then the reflected vector must be 6 times the basis vector = [ 6 0 0 ]. In summary:

> vector x is known
> basis vector times multiple is known = norm x. b
> these two vectors are reflected across a hyperplane
> to find the hyperplane, add the two column vectors = [ 8 4 4 ]


### 2.2 Portions parallel & perpendicular

For any two vectors x & v, vector x can be split into the portion in the same direction as v, as well as the portion not in the same direction (ie, perpendicular). Can define this as x = c.v + d.

For the portion in the same direction, c.v is the portion, c the multiple, so the objective is to find c. Dotting x with v, means the sum of x tranposed and multiplied by v, so performing this operation on the equation gives xT.v = c. vT.v + dT.v. Given d and v are perpendicular, then their dot product is zero. This gives xT.v = c. vT.v  which simplifies to c = xT.v / vT.v.

For the portion perpendicular, find d, where d = x - c.v. Then d = x - xT.v. v / vT.v.


### 2.3 Reflection matrix

To reflect each column of the matrix across the hyperplane, keep the portion in the same direction as the hyperplane, but reverse the portion perpendicular. Then c.v - (x - c.v) = 2.c.v - x = 2. xT.v. v / vT.v - x. This simplifies to the reflection matrix:

> = 2.v. xT.v / vT.v - x   as scalar. vector = vector. scalar
> = 2.v. vT.x / vT.v - x   as dot product is commutative
> = (2.v.vT / vT.v - identity). x


### 2.4 Reduced matrix calculations

We could multiply the entire nxn matrix by the reflection matrix embedded in an nxn identity matrix. However, to speed the calculations we can multiply by the smaller reflection matrix. By way of visual demonstration below, where only x influence the outcome. So in the below example we can reduce the calculations to multiplying 3x3 . 3x4 = 3x4.

[ 1 - - - ].[ b b b b ] = [ b b b b ]
[ - x x x ] [ x x x x ]   [ x x x x ]
[ - x x x ] [ x x x x ]   [ x x x x ]
[ - x x x ] [ x x x x ]   [ x x x x ]

Multiplying on the RHS 4x3 . 3x3 = 4x3.

[ a x x x ].[ 1 - - - ] = [ a x x x ]
[ a x x x ] [ - x x x ]   [ a x x x ]
[ a x x x ] [ - x x x ]   [ a x x x ]
[ a x x x ] [ - x x x ]   [ a x x x ]



## 3. QR decomposition function

### 3.1 Givens rotation function

With the Hessenberg matrix as input, Givens rotations can decompose the matrix into an orthonormal matrix times an upper triangular matrix.

The objective is to rotate the column vector [ a b ] from the diagonal and subdiagonal components, onto a column vector [ r 0 ]. Given the symmetry is rotation, the vectors have the same length, so r^2 = a^2 + b^2. This brings to mind trigonometry. The angle between the vectors has cos = a/r and sin = b/r. The matrix that clockwise rotates a b onto the basis vector is then:

[  cos  sin ].[ a ] = [  a/r  b/r ].[ a ] = [ r ]
[ -sin  cos ] [ b ]   [ -b/r  a/r ] [ b ]   [ 0 ]


### 3.2 Reduced matrix calculations, part 2

As before, we can reduce calculations. 

[ 1 - - - - ].[ b b b b b ] = [ b b b b b ]
[ - x x - - ] [ x x x x x ]   [ x x x x x ]
[ - x x - - ] [ - x x x x ]   [ x x x x x ]
[ - - - 1 - ] [ - - b b b ]   [ - - b b b ]
[ - - - - 1 ] [ - - - b b ]   [ - - - b b ]

[ a x x a a ].[ 1 - - - - ] = [ a x x a a ]
[ a x x a a ] [ - x x - - ]   [ a x x a a ]
[ - x x a a ] [ - x x - - ]   [ - x x a a ]
[ - - x a a ] [ - - - 1 - ]   [ - x x a a ]
[ - - - a a ] [ - - - - 1 ]   [ - - - a a ]


### 3.3 Wilkinson shift function

While the Givens rotations can be applied directly to the Hessenberg matrix, convergence would be slow. To speed convergence the Wilkinson shift essentially subtracts an eigenvalue from the matrix diagonal before applying the Givens rotation. The difference is an improvement from about 1250+ iterations without shift, to about 500+ iterations with shift for the example provided.

Instead of re-deriving the formula, we will just check the formula used for Wilkinson shift is consistent with the eigenvalue formula we derived in the theory section. Starting then with the formula used for Wilkinson shift:

> w = d - c^2 / ((a - d)/2 + (((a - d) /2) ^2 + c^2) ^1/2)
> w = d - c^2 / (eig - d)   subs from theory section  eig = (a + d)/2 +- (((a - d) /2) ^2 + b.c) ^1/2
> w. eig - w. d = d. eig - d^2 - c^2
> w. eig = d.(eig + w) - d^2 - c^2
> if w is actually one of the eigenvalues, then trace = sum eigenvalues = a + d
> eig1. eig2 = d.a + d^2 - d^2 - c^2
> eig1. eig2 = a.d - c^2

LHS is the product of eigenvalues, RHS is the determinant of the underlying matrix. From the theory section, these are equal. Check complete. Wilkinson shift gives one of the eigenvalues.


### 3.4 Power iteration

Ignoring subtleties of the decomposition function calculations for a moment:

> Givens function takes input M & returns Q R where M = Q.R
> so R = QT.M
> decomposition function recombines M2 as R.Q = (QT.M).Q
> iteration continues eg M4 = Q3T.Q2T.QT.M.Q.Q2.Q3

Repeated iteration essentially takes the matrix in the direction of the eigenvectors, until an arbitrary stopping point when c of the Wilkinson shift matrix is small enough.

ENDS