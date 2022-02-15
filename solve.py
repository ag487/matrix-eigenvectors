#!/usr/bin/env python3

import math
import numpy as np
from PIL import Image

from utilities import matrix, eigen

"""
image chosen is 512x512 so if take first 128 eigenvalues
results in noise of about 37 dB comparing rebuilt image to original
changing N_EIGS will change the resolution of the rebuilt image
"""

N_EIGS = 128

def import_image():
    print('image reduction by decomposition of covariance matrix')
    with Image.open("Lenna.png") as im:
        print('importing image..', im.format, im.mode, im.size)
        img_bw = np.asarray(im.convert(mode="L"))  # change to greyscale
        m = im.size[0]
        n = im.size[1]
    img_norm = np.empty([m, n])
    cols_mean = np.zeros([n])
    cols_std = np.zeros([n])
    for j in range(n):
        for i in range(m):
            cols_mean[j] += img_bw[i, j]
        cols_mean[j] /= m
        for i in range(m):
            cols_std[j] += pow(img_bw[i, j] - cols_mean[j], 2)
        cols_std[j] = pow(cols_std[j] / m, 0.5)
    for j in range(n):
        for i in range(m):
            img_norm[i, j] = (img_bw[i, j] - cols_mean[j]) / cols_std[j]
    img_cov = np.matmul(img_norm.T, img_norm) / (m - 1)
    return img_bw, img_norm, img_cov, cols_mean, cols_std


def find_eigs(img_cov):
    m = img_cov.shape[0]
    vals_sorted = np.empty([N_EIGS])
    vecs_sorted = np.empty([m, N_EIGS])
    print('finding eigenvectors..')
    Q1, H = matrix.hessenberg(img_cov)
    Q2, eig_vals = eigen.decomposition(H)
    eig_vecs = np.matmul(Q1, Q2)
    # above lines equivalent to   eig_vals, eig_vecs = np.linalg.eig(img_cov)
    abs_vals = np.asarray([abs(v) for v in eig_vals])
    v_indices = np.argsort(abs_vals)
    for i in range(N_EIGS):
        pos = v_indices[m-1-i]
        vals_sorted[i] = eig_vals[pos]
        vecs_sorted[:, i] = eig_vecs[:, pos]
    print("sorted by eigenvalue")
    print("vecs, vals:", vecs_sorted.shape, vals_sorted.shape)
    return vals_sorted, vecs_sorted


def project_components(img_norm, vecs):
    pc = np.matmul(img_norm, vecs)  # project image on reduced eigenvectors
    print("principal components:", pc.shape)
    ratio = img_norm.shape[1] / (pc.shape[1] + vecs.shape[1] + 2)  # add 2 to include mean & std vectors
    print("compression ratio:", ratio)
    img_rebuilt = np.matmul(pc, vecs.T)
    return img_rebuilt


def display(img_rebuilt, cols_std, cols_mean):
    array = []
    for j in range(512):
        line = []
        for i in range(512):
            line.append(img_rebuilt[i, j] * cols_std[j] + cols_mean[j])
        array.append(line)
    img_after = np.asarray(array).T
    img_after.astype(int)
    img2 = Image.fromarray(img_after)
    print('rebuilt image:', img2.format, img2.mode, img2.size)
    img2.show()
    return img_after


def noise(img_orig, img_after):
    m = img_orig.shape[0]
    n = img_orig.shape[1]
    mse = 0
    for i in range(m):
        for j in range(n):
            mse += pow(img_orig[i, j] - img_after[i, j], 2)
    mse /= m * n
    psnr = 10 * math.log(pow(255, 2) / mse, 10)  # max pixel value is 255
    print("peak signal to noise ratio in dB:", psnr)
    print("higher the better, typically lossy compression is 30-50 dB")
    return


def main():
    img_bw, img_norm, img_cov, cols_mean, cols_std = import_image()
    vals, vecs = find_eigs(img_cov)
    img_rebuilt = project_components(img_norm, vecs)
    img_after = display(img_rebuilt, cols_std, cols_mean)
    noise(img_bw, img_after)
    return

if __name__ == '__main__':
    main()
