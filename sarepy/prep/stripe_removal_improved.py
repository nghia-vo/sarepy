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
# Description: Improvements of stripe artifact removal methods:
# [1] N. T. Vo, R. C. Atwood, and M. Drakopoulos, "Superior techniques
#     for eliminating ring artifacts in X-ray micro-tomography," Optics
#     Express 26, 28396-28412 (2018). https://doi.org/10.1364/OE.26.028396.
# [2] N. T. Vo, R. C. Atwood, and M. Drakopoulos,"Preprocessing techniques
#     for removing artifacts in synchrotron-based tomographic images,"
#     Proc. SPIE 11113, Developments in X-Ray Tomography XII.
#     https://doi.org/10.1117/12.2530324.
# Publication date: 09th October 2019
#============================================================================

"""
Module for stripe removal methods proposed in:
https://doi.org/10.1117/12.2530324
"""

import numpy as np
from scipy import interpolate
from scipy.signal.windows import gaussian
from scipy.ndimage import median_filter
from scipy.ndimage import binary_dilation
# import scipy.fftpack as fft
import pyfftw.interfaces.scipy_fftpack as fft

from sarepy.prep.stripe_removal_original import remove_stripe_based_sorting
from sarepy.prep.stripe_removal_original import remove_stripe_based_fitting
from sarepy.prep.stripe_removal_original import apply_gaussian_filter
from sarepy.prep.stripe_removal_original import detect_stripe


def remove_stripe_based_filtering_sorting(sinogram, sigma, size, dim=1):
    """
    Removing stripes using the filtering and sorting technique, combination of
    algorithm 2 and algorithm 3 in Ref.[1]. Angular direction is along the axis 0.

    Parameters
    ----------
    sinogram : array_like
        2D array. Sinogram image.
    sigma : int
        Sigma of the Gaussian window used to separate the low-pass and
        high-pass components of the intensity profile of each column.
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
    sino_smooth_cor = np.transpose(
        remove_stripe_based_sorting(np.transpose(sino_smooth), size, dim))
    return np.transpose(sino_smooth_cor + sino_sharp)


def remove_stripe_based_sorting_fitting(sinogram, order, sigmax, sigmay):
    """
    Remove stripes using the sorting and fitting technique, combination of
    algorithm 2 and 1 in Ref. [1]. Angular direction is along the axis 0.

    Parameters
    ----------
    sinogram : array_like
        2D array. Sinogram image.
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
    sinogram = np.transpose(sinogram)
    (nrow, ncol) = sinogram.shape
    list_index = np.arange(0.0, ncol, 1.0)
    mat_index = np.tile(list_index, (nrow, 1))
    mat_comb = np.asarray(np.dstack((mat_index, sinogram)))
    mat_sort = np.asarray(
        [row[row[:, 1].argsort()] for row in mat_comb])
    sino_sort = mat_sort[:, :, 1]
    sino_cor = np.transpose(
        remove_stripe_based_fitting(np.transpose(sino_sort), order, sigmax, sigmay))
    mat_sort[:, :, 1] = sino_cor
    mat_sort_back = np.asarray(
        [row[row[:, 0].argsort()] for row in mat_sort])
    return np.transpose(mat_sort_back[:, :, 1])


def remove_stripe_based_2d_filtering_sorting(sinogram, sigma, size, dim=1):
    """
    Remove stripes using a 2D low-pass filter and the sorting-based technique,
    algorithm in section 3.3.4 in Ref. [1]. Angular direction is along the axis 0.

    Parameters
    ---------
    sinogram : array_like
        2D array. Sinogram image.
    sigma : int
        Sigma of the Gaussian window.
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
    .. [1] https://doi.org/10.1117/12.2530324
    """
    pad = min(150, int(0.1 * sinogram.shape[0]))
    sino_smooth = apply_gaussian_filter(sinogram, sigma, sigma, pad)
    sino_sharp = sinogram - sino_smooth
    sino_cor = remove_stripe_based_sorting(sino_sharp, size, dim)
    return sino_smooth + sino_cor


def remove_stripe_based_interpolation(sinogram, snr, size, drop_ratio=0.1, norm=True):
    """
    Combination of algorithm 4, 5, and 6 in Ref. [1].
    Remove stripes using a detection technique and an interpolation method.
    Angular direction is along the axis 0.

    Parameters
    ----------
    sinogram : array_like
        2D array. Sinogram image
    snr : float
        Ratio used to segment between useful information and noise.
    size : int
        Window size of the median filter used to detect stripes.
    drop_ratio : float, optional
        Ratio of pixels to be dropped, which is used to to reduce
        the possibility of the false detection of stripes.
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
    drop_ratio = np.clip(drop_ratio, 0.0, 0.8)
    sinogram = np.copy(sinogram)
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
        sinogram = sinogram / mat_fact
    list_mask[0:2] = 0.0
    list_mask[-2:] = 0.0
    listx = np.where(list_mask < 1.0)[0]
    listy = np.arange(nrow)
    matz = sinogram[:, listx]
    finter = interpolate.interp2d(listx, listy, matz, kind='linear')
    listx_miss = np.where(list_mask > 0.0)[0]
    if len(listx_miss) > 0:
        sinogram[:, listx_miss] = finter(listx_miss, listy)
    return sinogram
