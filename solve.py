#!/usr/bin/env python3

import cv2
import math
import numpy as np

import sys

from utilities import hess, qr


N_EIGS = 128

# image chosen is 512x512 so if take first 128 eigenvalues
# results in noise of about 38 dB comparing rebuilt image to original
# changing N_EIGS will change the resolution of the rebuilt image


def compress_image(img):
    dims = img.shape
    print(f"imported image..  pixel height {dims[0]}, width {dims[1]}")
        
    col_mean = np.sum(img, axis=0) / dims[0]
    img_norm = img - col_mean   # subtract mean from each row

    col_stdev = np.sum(img_norm * img_norm, axis=0) / dims[0]   # elementwise multiplication
    col_stdev = np.asarray([pow(s, 0.5) for s in col_stdev.tolist()])
    img_norm = img_norm / col_stdev   # divide each row

    img_cov = np.matmul(img_norm.T, img_norm) / (dims[0] - 1)

    print("finding eigenvectors..")

    # below lines could be replaced with   eig_vals, eig_vecs = np.linalg.eig(img_cov)
    Q1, H = hess.hessenberg(img_cov)
    Q2, eig_vals = qr.decomposition(H)    
    eig_vecs = np.matmul(Q1, Q2)
  
    abs_vals = np.asarray([abs(v) for v in eig_vals])
    v_indices = np.argsort(abs_vals)[::-1][:N_EIGS]   # highest N_EIGS, after reversing argsort which returns low to high
    vecs_sorted = eig_vecs[:, v_indices]
    vals_sorted = np.asarray(eig_vals)[v_indices]

    print(f"sorted by {N_EIGS} highest eigenvalues")
    print(f"vecs, vals shape: {vecs_sorted.shape}, {vals_sorted.shape}")

    components = np.matmul(img_norm, vecs_sorted)   # project image on reduced eigenvectors
    print(f"principal components shape {components.shape}")
    return components, vecs_sorted, col_mean, col_stdev


def rebuild_image(components, vecs, cols_mean, cols_std):
    img_after = np.matmul(components, vecs.T) * cols_std + cols_mean   # multiply and add for each row
    img_after = img_after.astype("uint8")
    return img_after


def calc_noise(img, img_after):
    dims = img.shape
    diff = img - img_after
    mse = np.sum(diff * diff) / (dims[0] * dims[1])   # elementwise multiplication
    psnr = 10 * math.log(pow(255, 2) / mse, 10)   # max pixel value is 255

    print(f"peak signal to noise ratio in dB: {psnr:.01f}")
    print("higher the better, typically lossy compression is 30-50 dB")
    return
    

if __name__ == '__main__':
    img = cv2.imread("Lenna.png", 0)   # reads as greyscale

    components, vecs, cols_mean, cols_std = compress_image(img)

    img_after = rebuild_image(components, vecs, cols_mean, cols_std)

    ratio = img.shape[1] / (components.shape[1] + vecs.shape[1] + 2)   # add 2 to include mean & stdev vectors
    print(f"compression ratio = {ratio:.01f}")
    calc_noise(img, img_after)

    cv2.imshow("rebuilt image", img_after)
    cv2.waitKey(0)   # shows image until window is selected and any key pressed