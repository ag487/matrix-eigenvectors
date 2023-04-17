This Python3 implementation is in two parts.

## Find eigenvectors of matrix
Calculates eigenvectors & eigenvalues of a symmetric matrix from scratch.
- uses Householder reflections to convert the matrix to Hessenberg form
- uses repeated Givens rotations to QR decompose the matrix, together with Wilkinson shift to speed convergence
- runtime should be a few seconds

I have included a theory document explaining as simply as possible the logic behind these matrix calculations, which was written for my own interest.

## Compress image using eigenvectors
Case study implements the matrix calculations.
- imports 512x512 pixel photo
- decomposes into matrices of 512 eigenvectors & eigenvalues
- takes n most significant eigenvectors (initially n=128)
- projects image onto significant eigenvectors, resulting in principal component matrix, which represents the image in compressed form
- rebuilds an image from the compressed form
- calculates noise difference between original and rebuilt photo - which results in a peak signal to noise ratio of about 38 dB (within the range of 30-50 dB usually considered reasonable for compression)

Part 1 could be replaced with numpy.linalg.eig, if you are only interested in the image decomposition.

ENDS