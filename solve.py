#!/usr/bin/env python3

import cv2 as cv
import math
import numpy as np

from utilities import hess, qr


N_EIGS = 128

# image chosen is 512x512 so if take first 128 eigenvalues
# results in noise of about 36 dB comparing rebuilt image to original
# changing N_EIGS will change the resolution of the rebuilt image


def custom_eig(img_cov):
    Q1, H = hess.hessenberg(img_cov)
    Q2, eig_vals = qr.decomposition(H)    
    eig_vecs = np.matmul(Q1, Q2)
    return eig_vals, eig_vecs


def compress_image(img):
    n_rows, n_cols = img.shape
    print(f"imported image..  pixel height {n_rows}, width {n_cols}")
        
    col_mean = img.sum(axis=0) / n_rows
    img_less_mean = img - col_mean   
    col_var = np.square(img_less_mean).sum(axis=0) / n_rows
    col_stdev = np.sqrt(col_var)
    img_norm = img_less_mean / col_stdev
    img_cov = np.matmul(img_norm.T, img_norm) / (n_rows - 1)
    
    print("finding eigenvectors..")
    #alternatively could replace the custom eigenvalue function with the numpy linalg function
    #eig_vals, eig_vecs = np.linalg.eig(img_cov)
    eig_vals, eig_vecs = custom_eig(img_cov)

    print("found eigenvectors")
    eig_vals_np = np.array(eig_vals)
    abs_vals = np.abs(eig_vals_np)
    v_indices = np.argsort(abs_vals)[::-1][:N_EIGS]   
    #argsort returns low to high, then reverse -> high to low, then take highest N_EIGS

    vecs_sorted = eig_vecs[:, v_indices]
    vals_sorted = eig_vals_np[v_indices]

    print(f"sorted by {N_EIGS} highest eigenvalues")
    print(f"vecs, vals shape: {vecs_sorted.shape}, {vals_sorted.shape}")

    components = np.matmul(img_norm, vecs_sorted)   #project image on reduced eigenvectors
    print("compressed photo")
    print(f"principal components shape {components.shape}")
    return components, vecs_sorted, col_mean, col_stdev


def rebuild_image(components, vecs, cols_mean, cols_std):
    img_after = np.matmul(components, vecs.T) * cols_std + cols_mean
    img_after = img_after.astype("uint8")
    return img_after


def calc_noise(img, img_after):
    n_rows, n_cols = img.shape
    diff = img - img_after
    mse = np.square(diff).sum() / (n_rows * n_cols)
    psnr = 10 * math.log(pow(255, 2) / mse, 10)   #max pixel value is 255

    print(f"comparison between photos shows peak signal to noise ratio of {psnr:.01f} dB")
    print("higher the better, typically lossy compression is 30-50 dB")
    return
    

if __name__ == '__main__':
    img = cv.imread("peppers.png", 0)   #reads USC peppers image as greyscale
    components, vecs, cols_mean, cols_std = compress_image(img)

    img_after = rebuild_image(components, vecs, cols_mean, cols_std)

    ratio = img.shape[1] / (components.shape[1] + vecs.shape[1] + 2)   #add 2 to include mean & stdev vectors
    print(f"compression ratio = {ratio:.01f}")
    calc_noise(img, img_after)

    cv.imshow("rebuilt image", img_after)
    cv.waitKey(0)   #shows image until any key pressed

    