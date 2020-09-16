#============================================================================
# Copyright (c) 2018 Nghia T. Vo. All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#============================================================================
# Author: Nghia T. Vo
# E-mail: nghia.vo@diamond.ac.uk
# Description: Methods for removing ring artifacts in reconstructed images by
# combining the polar transform and stripe removal methods.
# Publication date: 05/05/2020
#============================================================================

"""
Module for ring removal methods which work on reconstruction images
by combining the polar transform and stripe removal methods.
"""

import numpy as np
from scipy.ndimage import map_coordinates
import sarepy.prep.stripe_removal_former as srm


def cirle_mask(width, ratio):
    """
    Create a circle mask.

    Parameters
    -----------
    width : int
        Width of a square array.
    ratio : float
        Ratio between the diameter of the mask and the width of the array.

    Returns
    ------
    ndarray
         Square array.
    """
    mask = np.zeros((width, width), dtype=np.float32)
    center = width // 2
    radius = ratio * center
    y, x = np.ogrid[-center:width - center, -center:width - center]
    mask_check = x * x + y * y <= radius * radius
    mask[mask_check] = 1.0
    return mask


def reg_to_polar(width_reg, height_reg, width_pol, height_pol):
    """
    Generate coordinates of a rectangular grid from polar coordinates.

    Parameters
    -----------
    width_reg : int
        Width of the image in the Cartesian coordinate system.
    height_reg : int
        Height of the image in the Cartesian coordinate system.
    width_pol : int
        Width of the image in the polar coordinate system.
    height_pol : int
        Height of the image in the polar coordinate system.

    Returns
    ------
    x_mat : ndarray
         2D array. Broadcast of the x-coordinates.
    y_mat : ndarray
         2D array. Broadcast of the y-coordinates.

    """
    centerx = (width_reg - 1.0) * 0.5
    centery = (height_reg - 1.0) * 0.5
    r_max = np.floor(max(centerx, centery))
    r_list = np.linspace(0.0, r_max, width_pol)
    theta_list = np.arange(0.0, height_pol, 1.0) * 2 * np.pi / (height_pol - 1)
    r_mat, theta_mat = np.meshgrid(r_list, theta_list)
    x_mat = np.float32(
        np.clip(centerx + r_mat * np.cos(theta_mat), 0, width_reg - 1))
    y_mat = np.float32(
        np.clip(centery + r_mat * np.sin(theta_mat), 0, height_reg - 1))
    return x_mat, y_mat


def polar_to_reg(width_pol, height_pol, width_reg, height_reg):
    """
    Generate polar coordinates from grid coordinates.

    Parameters
    -----------
    width_pol : int
        Width of the image in the polar coordinate system.
    height_pol : int
        Height of the image in the polar coordinate system.
    width_reg : int
        Width of the image in the Cartesian coordinate system.
    height_reg : int
        Height of the image in the Cartesian coordinate system.

    Returns
    ------
    r_mat : ndarray
         2D array. Broadcast of the r-coordinates.
    theta_mat : ndarray
         2D array. Broadcast of the theta-coordinates.

    """
    centerx = (width_reg - 1.0) * 0.5
    centery = (height_reg - 1.0) * 0.5
    r_max = np.floor(max(centerx, centery))
    x_list = (1.0 * np.flipud(np.arange(width_reg)) -
              centerx) * width_pol / r_max
    y_list = (1.0 * np.flipud(np.arange(height_reg)) -
              centery) * width_pol / r_max
    x_mat, y_mat = np.meshgrid(x_list, y_list)
    r_mat = np.float32(np.clip(np.sqrt(x_mat**2 + y_mat**2), 0, width_pol - 1))
    theta_mat = np.float32(np.clip((np.pi + np.arctan2(y_mat, x_mat))
                                   * (height_pol - 1) / (2 * np.pi), 0, height_pol - 1))
    return r_mat, theta_mat


def mapping(mat, xmat, ymat):
    """
    Mapping an 2D array at coordinates (xmat, ymat) to rectangular
    grid coordinates.

    Parameters
    ----------
    mat : array_like
        2D array.
    xmat : array_like
        2D array of the x-coordinates.
    ymat : array_like
        2D array of the y-coordinates.

    Returns
    -------
    ndarray
        2D array.
    """
    indices = np.reshape(ymat, (-1, 1)), np.reshape(xmat, (-1, 1))
    mat = map_coordinates(mat, indices, order=1, mode='reflect')
    return mat.reshape(xmat.shape)


def ring_removal_based_wavelet_fft(mat, level=5, sigma=1.0, order=8, pad=150):
    """
    Remove ring artifacts in the reconstructed image by combining the polar
    transform and the wavelet-fft-based method (Ref. [1]).

    Parameters
    ----------
    mat : array_like
        Square array. Reconstructed image
    level : int
        Wavelet decomposition level.
    sigma : int
        Damping parameter. Larger is stronger.
    order : int
        Order of the the Daubechies wavelets.
    pad : int
        Padding for FFT

    Returns
    -------
    ndarray
        Square array. Ring-removed image.

    References
    ----------
    .. [1] https://doi.org/10.1364/OE.17.008567

    """
    (nrow, ncol) = mat.shape
    if nrow != ncol:
        raise ValueError(
            "Width and height of the reconstructed image are not the same")
    mask = cirle_mask(ncol, 1.0)
    (xmat, ymat) = reg_to_polar(ncol, ncol, ncol, ncol)
    (rmat, thetamat) = polar_to_reg(ncol, ncol, ncol, ncol)
    polar_mat = mapping(mat, xmat, ymat)
    polar_mat = srm.remove_stripe_based_wavelet_fft(
        polar_mat, level, sigma, order, pad)
    mat_rec = mapping(polar_mat, rmat, thetamat)
    return mat_rec * mask


def ring_removal_based_fft(mat, u=30, n=8, v=1, pad=150):
    """
    Remove ring artifacts in the reconstructed image by combining the polar
    transform and the fft-based method (Ref. [1]).

    Parameters
    ----------
    mat : array_like
        Square array. Reconstructed image
    u,n : int
        To define the shape of 1D Butterworth low-pass filter.
    v : int
        Number of rows (* 2) to be applied the filter.
    pad : int
        Padding for FFT

    Returns
    -------
    ndarray
        Square array. Ring-removed image.

    References
    ----------
    .. [1] https://doi.org/10.1063/1.1149043
    """
    (nrow, ncol) = mat.shape
    if nrow != ncol:
        raise ValueError(
            "Width and height of the reconstructed image are not the same")
    mask = cirle_mask(ncol, 1.0)
    (xmat, ymat) = reg_to_polar(ncol, ncol, ncol, ncol)
    (rmat, thetamat) = polar_to_reg(ncol, ncol, ncol, ncol)
    polar_mat = mapping(mat, xmat, ymat)
    polar_mat = srm.remove_stripe_based_fft(polar_mat, u, n, v, pad)
    mat_rec = mapping(polar_mat, rmat, thetamat)
    return mat_rec * mask
