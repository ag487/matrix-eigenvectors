'''
Copyright (c) 2021 Angus Graham
All rights reserved.

This source code is licensed under the MIT License found in the
license.md file in the root directory of this source tree.
'''

import math
import numpy as np
from PIL import Image

from utilities import matrix, eigen


# image chosen is 512x512 so if take first 128 eigenvalues
# results in noise of about 37 dB comparing rebuilt image to original
# which is within range of 30-50 dB usually considered reasonable for compression
n_eigs = 128


def import_image():
  print('image reduction by decomposition of covariance matrix')
  img1 = Image.open("Lenna.png")
  m = img1.size[0]
  n = img1.size[1]
  print('importing image..', img1.format, img1.mode, img1.size)
  img2 = np.asarray(img1.convert(mode="L"))  # change to greyscale
  img_norm = np.empty([m,n])
  img_mean = np.zeros([n])
  img_std = np.zeros([n])
  for j in range(n):
    for i in range(m):
      img_mean[j] += img2[i,j]
    img_mean[j] /= m
    for i in range(m):
      img_std[j] += (img2[i,j] - img_mean[j]) ** 2
    img_std[j] = (img_std[j] / m) ** 0.5
    for i in range(m):
      img_norm[i,j] = (img2[i,j] - img_mean[j]) / img_std[j]
  return img2, img_norm, img_mean, img_std

def decompose_image(img_norm):
  print('creating cov matrix..')
  img_cov = np.matmul(img_norm.T, img_norm) / (img_norm.shape[0] - 1)
  print('finding eigenvectors..')
  Q, H = matrix.hessenberg(img_cov)
  Q2, eig_vals = eigen.decomposition(H)
  eig_vecs = np.matmul(Q,Q2)
  # above lines equivalent to   eig_vals, eig_vecs = np.linalg.eig(img_cov)
  abs_vals = np.asarray([ abs(v) for v in eig_vals ])
  v_indices = np.argsort(abs_vals)
  v_count = abs_vals.shape[0]
  vals_sorted = np.empty([n_eigs])
  vecs_sorted = np.empty([v_count,n_eigs])
  for i in range(n_eigs):
    pos = v_indices[v_count-1-i]
    vals_sorted[i] = eig_vals[pos]
    vecs_sorted[:,i] = eig_vecs[:,pos]
  print("sorted by eigenvalue")
  print("vecs,vals:", vecs_sorted.shape, vals_sorted.shape)
  return vals_sorted, vecs_sorted

def PCA(img_norm, vecs):
  pc = np.matmul(img_norm, vecs)  # project image on reduced eigenvectors
  print("principal components:", pc.shape)
  ratio = img_norm.shape[1] / (pc.shape[1] + vecs.shape[1] + 2)
  print("compression ratio:", ratio)  # add 2 to include mean & std vectors
  img_rebuilt = np.matmul(pc, vecs.T)
  return img_rebuilt

def display(img_rebuilt, img_std, img_mean):
  img3_array = np.asarray([[ img_rebuilt[i,j] * img_std[j] + img_mean[j] \
    for i in range(512) ] for j in range(512) ]).T
  img3_array.astype(int)
  img3 = Image.fromarray(img3_array)
  print('rebuilt image:', img3.format, img3.mode, img3.size)
  img3.show()
  return img3_array

def noise(img_orig, img_after):
  m = img_orig.shape[0]
  n = img_orig.shape[1]
  mse = 0
  for i in range(m):
    for j in range(n):
      mse += (img_orig[i,j] - img_after[i,j]) ** 2
  mse /= m * n
  psnr = 10 * math.log(255 ** 2 / mse, 10)  # max pixel value is 255
  print("peak signal to noise ratio in dB:", psnr)
  print("higher the better, typically lossy compression is 30-50 dB")
  return

img_orig, img_norm, img_mean, img_std = import_image()
vals, vecs = decompose_image(img_norm)
img_rebuilt = PCA(img_norm, vecs)
img_after = display(img_rebuilt, img_std, img_mean)
noise(img_orig, img_after)
