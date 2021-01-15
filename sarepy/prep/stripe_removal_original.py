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
# Description: Original implementation of stripe artifact removal methods,
# Nghia T. Vo, Robert C. Atwood, and Michael Drakopoulos, "Superior
# techniques for eliminating ring artifacts in X-ray micro-tomography," Optics
# Express 26, 28396-28412 (2018).
# https://doi.org/10.1364/OE.26.028396
# Publication date: 18th October 2018
#============================================================================

"""
Module for stripe removal methods proposed in:
https://doi.org/10.1364/OE.26.028396
"""

import numpy as np
from scipy.signal.windows import gaussian
from scipy.signal import savgol_filter
from scipy.ndimage import median_filter
from scipy.ndimage import binary_dilation
from scipy.ndimage import uniform_filter1d
from scipy import interpolate
# import scipy.fftpack as fft
import pyfftw.interfaces.scipy_fftpack as fft


def remove_stripe_based_sorting(sinogram, size, dim=1):
    """
    Remove stripe artifacts in a sinogram using the sorting technique,
    algorithm 3 in Ref. [1]. Angular direction is along the axis 0.

    Parameters
    ----------
    sinogram : array_like
        2D array. Sinogram image.
    size : int
        Window size of the median filter.
    dim : {1, 2}, optional
        Dimension of the window.

    Returns
    -------
    ndarray
        2D array. Stripe-removed sinogram.

    References
    ----------
    .. [1] https://doi.org/10.1364/OE.26.028396
    """
    sinogram = np.transpose(sinogram)
    (nrow, ncol) = sinogram.shape
    list_index = np.arange(0.0, ncol, 1.0)
    mat_index = np.tile(list_index, (nrow, 1))
    mat_comb = np.asarray(np.dstack((mat_index, sinogram)))
    mat_sort = np.asarray(
        [row[row[:, 1].argsort()] for row in mat_comb])
    if dim == 2:
        mat_sort[:, :, 1] = median_filter(mat_sort[:, :, 1], (size, size))
    else:
        mat_sort[:, :, 1] = median_filter(mat_sort[:, :, 1], (size, 1))
    mat_sort_back = np.asarray(
        [row[row[:, 0].argsort()] for row in mat_sort])
    return np.transpose(mat_sort_back[:, :, 1])


def remove_stripe_based_filtering(sinogram, sigma, size, dim=1):
    """
    Remove stripe artifacts in a sinogram using the filtering technique,
    algorithm 2 in Ref. [1]. Angular direction is along the axis 0.

    Parameters
    ----------
    sinogram : array_like
        2D array. Sinogram image
    sigma : int
        Sigma of the Gaussian window used to separate the low-pass and
        high-pass components of the intensity profile of each column.
    size : int
        Window size of the median filter.
    dim : {1, 2}, optional
        Dimension of the window.

    Returns
    -------
    array_like
         2D array. Stripe-removed sinogram.

    References
    ----------
    .. [1] https://doi.org/10.1364/OE.26.028396
    """
    pad = min(150, int(0.1 * sinogram.shape[0]))
    sinogram = np.transpose(sinogram)
    sino_pad = np.pad(sinogram, ((0, 0), (pad, pad)), mode='reflect')
    (_, ncol) = sino_pad.shape
    window = gaussian(ncol, std=sigma)
    list_sign = np.power(-1.0, np.arange(ncol))
    sino_smooth = np.copy(sinogram)
    for i, sino_1d in enumerate(sino_pad):
        sino_smooth[i] = np.real(
            fft.ifft(fft.fft(sino_1d * list_sign) * window) * list_sign)[pad:ncol - pad]
    sino_sharp = sinogram - sino_smooth
    if dim == 2:
        sino_smooth_cor = median_filter(sino_smooth, (size, size))
    else:
        sino_smooth_cor = median_filter(sino_smooth, (size, 1))
    return np.transpose(sino_smooth_cor + sino_sharp)


def make_2d_gaussian_window(height, width, sigmax, sigmay):
    """
    Create a 2D Gaussian window.

    Parameters
    ----------
    height : int
        Height of the image.
    width : int
        Width of the image.
    sigmax : int
        Sigma in the x-direction.
    sigmay : int
        Sigma in the y-direction.

    Returns
    -------
    ndarray
        2D array.
    """
    xcenter = (width - 1.0) / 2.0
    ycenter = (height - 1.0) / 2.0
    y, x = np.ogrid[-ycenter:height - ycenter, -xcenter:width - xcenter]
    window = np.exp(-(x ** 2 / (2 * sigmax ** 2) + y ** 2 / (2 * sigmay ** 2)))
    return window


def apply_gaussian_filter(mat, sigmax, sigmay, pad):
    """
    Filtering an image using a 2D Gaussian window.

    Parameters
    ----------
    mat : array_like
        2D array.
    sigmax : int
        Sigma in the x-direction.
    sigmay : int
        Sigma in the y-direction.
    pad : int
        Padding for the Fourier transform.

    Returns
    -------
    ndarray
        2D array. Filtered image.
    """
    mat_pad = np.pad(mat, ((0, 0), (pad, pad)), mode='edge')
    mat_pad = np.pad(mat_pad, ((pad, pad), (0, 0)), mode='mean')
    (nrow, ncol) = mat_pad.shape
    window = make_2d_gaussian_window(nrow, ncol, sigmax, sigmay)
    listx = np.arange(0, ncol)
    listy = np.arange(0, nrow)
    x, y = np.meshgrid(listx, listy)
    mat_sign = np.power(-1.0, x + y)
    mat_filt = np.real(
        fft.ifft2(fft.fft2(mat_pad * mat_sign) * window) * mat_sign)
    return mat_filt[pad:nrow - pad, pad:ncol - pad]


def remove_stripe_based_fitting(sinogram, order, sigmax, sigmay):
    """
    Remove stripes using the fitting technique, algorithm 1 in Ref. [1].
    Angular direction is along the axis 0.

    Parameters
    ----------
    sinogram : array_like
        2D array. Sinogram image
    order : int
        Polynomial fit order.
    sigmax : int
        Sigma of the Gaussian window in the x-direction.
    sigmay : int
        Sigma of the Gaussian window in the y-direction.

    Returns
    -------
    ndarray
        2D array. Stripe-removed sinogram.

    References
    ----------
    .. [1] https://doi.org/10.1364/OE.26.028396
    """
    (nrow, _) = sinogram.shape
    pad = min(150, int(0.1 * sinogram.shape[0]))
    nrow1 = nrow
    if nrow1 % 2 == 0:
        nrow1 = nrow1 - 1
    sino_fit = savgol_filter(sinogram, nrow1, order, axis=0, mode='mirror')
    sino_fit_smooth = apply_gaussian_filter(sino_fit, sigmax, sigmay, pad)
    num1 = np.mean(sino_fit)
    num2 = np.mean(sino_fit_smooth)
    sino_fit_smooth = num1 * sino_fit_smooth / num2
    return (sinogram / sino_fit) * sino_fit_smooth


def detect_stripe(list_data, snr):
    """
    Locate stripe positions using Algorithm 4 in Ref. [1].

    Parameters
    ----------
    list_data : array_like
        1D array. Normalized data.
    snr : float
        Ratio used to segment stripes from background noise.

    Returns
    -------
    ndarray
        1D binary mask.

    References
    ----------
    .. [1] https://doi.org/10.1364/OE.26.028396
    """
    npoint = len(list_data)
    list_sort = np.sort(list_data)
    listx = np.arange(0, npoint, 1.0)
    ndrop = np.int16(0.25 * npoint)
    (slope, intercept) = np.polyfit(listx[ndrop:-ndrop - 1], list_sort[ndrop:-ndrop - 1], 1)
    y_end = intercept + slope * listx[-1]
    noise_level = np.abs(y_end - intercept)
    noise_level = np.clip(noise_level, 1e-6, None)
    val1 = np.abs(list_sort[-1] - y_end) / noise_level
    val2 = np.abs(intercept - list_sort[0]) / noise_level
    list_mask = np.zeros(npoint, dtype=np.float32)
    if val1 >= snr:
        upper_thresh = y_end + noise_level * snr * 0.5
        list_mask[list_data > upper_thresh] = 1.0
    if val2 >= snr:
        lower_thresh = intercept - noise_level * snr * 0.5
        list_mask[list_data <= lower_thresh] = 1.0
    return list_mask


def remove_large_stripe(sinogram, snr, size, drop_ratio=0.1, norm=True):
    """
    Remove large stripes, algorithm 5 in Ref. [1], by: locating stripes,
    normalizing to remove full stripes, and using the sorting technique
    (Ref. [1]) to remove partial stripes. Angular direction is along the
    axis 0.

    Parameters
    ----------
    sinogram : array_like
        2D array. Sinogram image
    snr : float
        Ratio used to segment stripes from background noise.
    size : int
        Window size of the median filter.
    drop_ratio : float, optional
        Ratio of pixels to be dropped, which is used to reduce the false
        detection of stripes.
    norm : bool, optional
        Apply normalization if True.

    Returns
    -------
    ndarray
        2D array. Stripe-removed sinogram.

    References
    ----------
    .. [1] https://doi.org/10.1364/OE.26.028396
    """
    sinogram = np.copy(sinogram)  # Make it mutable
    drop_ratio = np.clip(drop_ratio, 0.0, 0.8)
    (nrow, ncol) = sinogram.shape
    ndrop = int(0.5 * drop_ratio * nrow)
    sino_sort = np.sort(sinogram, axis=0)
    sino_smooth = median_filter(sino_sort, (1, size))
    list1 = np.mean(sino_sort[ndrop:nrow - ndrop], axis=0)
    list2 = np.mean(sino_smooth[ndrop:nrow - ndrop], axis=0)
    list_fact = np.divide(list1, list2,
                         out=np.ones_like(list1), where=list2 != 0)
    list_mask = detect_stripe(list_fact, snr)
    list_mask = np.float32(binary_dilation(list_mask, iterations=1))
    mat_fact = np.tile(list_fact, (nrow, 1))
    if norm is True:
        sinogram = sinogram / mat_fact  # Normalization
    sino_tran = np.transpose(sinogram)
    list_index = np.arange(0.0, nrow, 1.0)
    mat_index = np.tile(list_index, (ncol, 1))
    mat_comb = np.asarray(np.dstack((mat_index, sino_tran)))
    mat_sort = np.asarray(
        [row[row[:, 1].argsort()] for row in mat_comb])
    mat_sort[:, :, 1] = np.transpose(sino_smooth)
    mat_sort_back = np.asarray(
        [row[row[:, 0].argsort()] for row in mat_sort])
    sino_cor = np.transpose(mat_sort_back[:, :, 1])
    listx_miss = np.where(list_mask > 0.0)[0]
    sinogram[:, listx_miss] = sino_cor[:, listx_miss]
    return sinogram


def remove_unresponsive_and_fluctuating_stripe(sinogram, snr, size, residual=False):
    """
    Remove unresponsive or fluctuating stripes, algorithm 6 in Ref. [1], by:
    locating stripes, correcting using interpolation. Angular direction
    is along the axis 0.

    Parameters
    ----------
    sinogram : array_like
        2D array. Sinogram image.
    snr : float
        Ratio used to segment stripes from background noise.
    size : int
        Window size of the median filter.
    residual : bool, optional
        Removing residual stripes if True.

    Returns
    -------
    ndarray
        2D array. Stripe-removed sinogram.

    References
    ----------
    .. [1] https://doi.org/10.1364/OE.26.028396
    """
    sinogram = np.copy(sinogram)  # Make it mutable
    (nrow, _) = sinogram.shape
    sino_smooth = np.apply_along_axis(uniform_filter1d, 0, sinogram, 10)
    list_diff = np.sum(np.abs(sinogram - sino_smooth), axis=0)
    list_diff_bck = median_filter(list_diff, size)
    list_fact = np.divide(list_diff, list_diff_bck,
                         out=np.ones_like(list_diff), where=list_diff_bck != 0)
    list_mask = detect_stripe(list_fact, snr)
    list_mask = np.float32(binary_dilation(list_mask, iterations=1))
    list_mask[0:2] = 0.0
    list_mask[-2:] = 0.0
    listx = np.where(list_mask < 1.0)[0]
    listy = np.arange(nrow)
    mat = sinogram[:, listx]
    finter = interpolate.interp2d(listx, listy, mat, kind='linear')
    listx_miss = np.where(list_mask > 0.0)[0]
    if len(listx_miss) > 0:
        sinogram[:, listx_miss] = finter(listx_miss, listy)
    if residual is True:
        sinogram = remove_large_stripe(sinogram, snr, size)
    return sinogram


def remove_all_stripe(sinogram, snr, la_size, sm_size, drop_ratio=0.1, norm=True, dim=1):
    """
    Remove all types of stripe artifacts by combining algorithm 6, 5, 4,
    and 3 in Ref. [1]. Angular direction is along the axis 0.

    Parameters
    ----------
    sinogram : array_like
        2D array. Sinogram image.
    snr : float
        Ratio used to segment stripes from background noise.
    la_size : int
        Window size of the median filter to remove large stripes.
    sm_size : int
        Window size of the median filter to remove small-to-medium stripes.
    drop_ratio : float, optional
        Ratio of pixels to be dropped, which is used to to reduce
        the possibility of the false detection of stripes.
    norm : bool, optional
        Apply normalization if True.
    dim : {1, 2}, optional
        Dimension of the window.

    Returns
    -------
    ndarray
        2D array. Stripe-removed sinogram.

    References
    ----------
    .. [1] https://doi.org/10.1364/OE.26.028396
    """
    sinogram = remove_unresponsive_and_fluctuating_stripe(
        sinogram, snr, la_size)
    sinogram = remove_large_stripe(sinogram, snr, la_size, drop_ratio, norm)
    sinogram = remove_stripe_based_sorting(sinogram, sm_size, dim)
    return sinogram
