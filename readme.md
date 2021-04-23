This Python implementation is in two parts.

# Finding eigenvectors of matrix

The main part calculates the eigenvectors & eigenvalues of a symmetric matrix from
scratch.
- uses Householder reflections to convert a matrix to Hessenberg form
- uses Givens rotations to QR decompose the Hessenberg matrix
- then uses power iteration with Wilkinson shift to find the eigenvectors &
eigenvalues
- these can be checked against the numpy.linalg.eig function for correctness
- runtime should be a few seconds

I have included a document explaining as simply as possible the logic behind
these matrix calculations, which was written for my own interest.

# Image decomposition

The second part is a case study implementing the matrix calculations.
- imports 512x512 pixel photo
- decomposes the photo into matrices of 512 eigenvectors & eigenvalues, using
the above
- takes n most significant eigenvalues (initially n=128)
- projects image onto significant eigenvectors, resulting in principal component
matrix, which represents the image in compressed form
- rebuilds the compressed form & displays the image
- calculates noise difference between original and rebuilt photo

If you just want to use the second part without the matrix calculations,
you could use np.linalg.eig in the decompose_image function, in which case
the utilities folder is not needed.

# Results

Using 128 eigenvectors on a 512x512 image results in a compression ratio of about
2.0, giving a rebuilt image with a peak signal to noise ratio of about 37 dB.

This is within the range of 30-50 dB usually considered reasonable for
compression.

ENDS
