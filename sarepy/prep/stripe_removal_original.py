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
#from scipy.fftpack import fft, ifft, fft2, ifft2
import pyfftw.interfaces.scipy_fftpack as fft


def remove_stripe_based_sorting(sinogram, size, dim=1):
    """
    Algorithm 3 in the paper. Remove stripes using the sorting technique.
    Work particularly well for removing partial stripes.
    Angular direction is along the axis 0.

    Parameters
    ----------
    sinogram : float
        2D array
    size : int
        Window size of the median filter.
    dim : {1, 2}, optional
        Dimension of the window.

    Returns
    -------
    float
        2D array. Stripe-removed sinogram.
    """
    sinogram = np.transpose(sinogram)
    (nrow, ncol) = sinogram.shape
    listindex = np.arange(0.0, ncol, 1.0)
    matindex = np.tile(listindex, (nrow, 1))
    matcomb = np.asarray(np.dstack((matindex, sinogram)))
    matsort = np.asarray(
        [row[row[:, 1].argsort()] for row in matcomb])
    if dim == 2:
        matsort[:, :, 1] = median_filter(matsort[:, :, 1], (size, size))
    else:
        matsort[:, :, 1] = median_filter(matsort[:, :, 1], (size, 1))
    matsortback = np.asarray(
        [row[row[:, 0].argsort()] for row in matsort])
    return np.transpose(matsortback[:, :, 1])


def remove_stripe_based_filtering(sinogram, sigma, size, dim=1):
    """
    Algorithm 2 in the paper. Remove stripes using the filtering technique.
    Angular direction is along the axis 0.

    Parameters
    ----------
    sinogram : float
        2D array.
    sigma : int
        Sigma of the Gaussian window used to separate the low-pass and
        high-pass components of the intensity profile of each column.
    size : int
        Window size of the median filter.
    dim : {1, 2}, optional
        Dimension of the window.

    Returns
    -------
    float
         2D array. Stripe-removed sinogram.
    """
    pad = 150  # To reduce artifacts caused by FFT
    sinogram = np.transpose(sinogram)
    sinopad = np.pad(sinogram, ((0, 0), (pad, pad)), mode='reflect')
    (_, ncol) = sinopad.shape
    window = gaussian(ncol, std=sigma)
    listsign = np.power(-1.0, np.arange(ncol))
    sinosmooth = np.copy(sinogram)
    for i, sinolist in enumerate(sinopad):
        # sinosmooth[i] = np.real(
        #     ifft(fft(sinolist*listsign)*window)*listsign)[pad:ncol-pad]
        sinosmooth[i] = np.real(
            fft.ifft(fft.fft(sinolist * listsign) * window) * listsign)[pad:ncol - pad]
    sinosharp = sinogram - sinosmooth
    if dim == 2:
        sinosmooth_cor = median_filter(sinosmooth, (size, size))
    else:
        sinosmooth_cor = median_filter(sinosmooth, (size, 1))
    return np.transpose(sinosmooth_cor + sinosharp)


def _2d_window_ellipse(height, width, sigmax, sigmay):
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
    window : float
        2D array.
    """
    centerx = (width - 1.0) / 2.0
    centery = (height - 1.0) / 2.0
    y, x = np.ogrid[-centery:height - centery, -centerx:width - centerx]
    numx = 2.0 * sigmax * sigmax
    numy = 2.0 * sigmay * sigmay
    window = np.exp(-(x * x / numx + y * y / numy))
    return window


def _2d_filter(mat, sigmax, sigmay, pad):
    """
    Filtering an image using a 2D Gaussian window.

    Parameters
    ----------
    mat : float
        2D array.
    sigmax : int
        Sigma in the x-direction.
    sigmay : int
        Sigma in the y-direction.
    pad : int
        Padding for FFT

    Returns
    -------
    mat_filt : float
        2D array. Filtered image.
    """
    matpad = np.pad(mat, ((0, 0), (pad, pad)), mode='edge')
    matpad = np.pad(matpad, ((pad, pad), (0, 0)), mode='mean')
    (nrow, ncol) = matpad.shape
    win2d = _2d_window_ellipse(nrow, ncol, sigmax, sigmay)
    listx = np.arange(0, ncol)
    listy = np.arange(0, nrow)
    x, y = np.meshgrid(listx, listy)
    matsign = np.power(-1.0, x + y)
    # mat_filt = np.real(ifft2(fft2(matpad * matsign) * win2d) * matsign)
    mat_filt = np.real(
        fft.ifft2(fft.fft2(matpad * matsign) * win2d) * matsign)
    return mat_filt[pad:nrow - pad, pad:ncol - pad]


def remove_stripe_based_fitting(sinogram, order, sigmax, sigmay):
    """
    Algorithm 1 in the paper. Remove stripes using the fitting technique.
    Angular direction is along the axis 0.

    Parameters
    ----------
    sinogram : float
        2D array.
    order : int
        Polynomial fit order.
    sigmax : int
        Sigma of the Gaussian window in the x-direction.
    sigmay : int
        Sigma of the Gaussian window in the y-direction.

    Returns
    -------
    float
        2D array. Stripe-removed sinogram.
    """
    (nrow, _) = sinogram.shape
    nrow1 = nrow
    if nrow1 % 2 == 0:
        nrow1 = nrow1 - 1
    sinofit = savgol_filter(sinogram, nrow1, order, axis=0, mode='mirror')
    sinofit_smooth = _2d_filter(sinofit, sigmax, sigmay, 50)
    num1 = np.mean(sinofit)
    num2 = np.mean(sinofit_smooth)
    sinofit_smooth = num1 * sinofit_smooth / num2
    return sinogram / sinofit * sinofit_smooth


def detect_stripe(listdata, snr):
    """
    Algorithm 4 in the paper. Used to locate stripe positions.

    Parameters
    ----------
    listdata : float
        1D normalized array.
    snr : float
        Ratio used to segment between useful information and noise.

    Returns
    -------
    listmask : float
        1D binary mask.
    """
    numdata = len(listdata)
    listsorted = np.sort(listdata)
    xlist = np.arange(0, numdata, 1.0)
    ndrop = np.int16(0.25 * numdata)
    (slope, intercept) = np.polyfit(
        xlist[ndrop:-ndrop - 1], listsorted[ndrop:-ndrop - 1], 1)
    yend = intercept + slope * xlist[-1]
    noiselevel = np.abs(yend - intercept)
    val1 = np.abs(listsorted[-1] - yend) / noiselevel
    val2 = np.abs(intercept - listsorted[0]) / noiselevel
    listmask = np.zeros(numdata, dtype=np.float32)
    if val1 >= snr:
        upper_thresh = yend + noiselevel * snr * 0.5
        listmask[listdata > upper_thresh] = 1.0
    if val2 >= snr:
        lower_thresh = intercept - noiselevel * snr * 0.5
        listmask[listdata <= lower_thresh] = 1.0
    return listmask


def remove_large_stripe(sinogram, snr, size):
    """
    Algorithm 5 in the paper. Remove large stripes by: locating stripes,
    normalizing to remove full stripes, using the sorting technique to
    remove partial stripes.
    Angular direction is along the axis 0.

    Parameters
    ----------
    sinogram : float
        2D array.
    snr : float
        Ratio used to segment between useful information and noise.
    size : int
        Window size of the median filter.

    Returns
    -------
    float
        2D array. Stripe-removed sinogram.
    """
    badpixel_ratio = 0.05
    (nrow, ncol) = sinogram.shape
    ndrop = np.int16(badpixel_ratio * nrow)
    sinosort = np.sort(sinogram, axis=0)
    sinosmooth = median_filter(sinosort, (1, size))
    list1 = np.mean(sinosort[ndrop:nrow - ndrop], axis=0)
    list2 = np.mean(sinosmooth[ndrop:nrow - ndrop], axis=0)
    listfact = list1 / list2
    listmask = detect_stripe(listfact, snr)
    listmask = np.float32(binary_dilation(listmask, iterations=1))
    matfact = np.tile(listfact, (nrow, 1))
    sinogram = sinogram / matfact
    sino_tran = np.transpose(sinogram)
    listindex = np.arange(0.0, nrow, 1.0)
    matindex = np.tile(listindex, (ncol, 1))
    matcomb = np.asarray(np.dstack((matindex, sino_tran)))
    matsort = np.asarray(
        [row[row[:, 1].argsort()] for row in matcomb])
    matsort[:, :, 1] = np.transpose(sinosmooth)
    matsortback = np.asarray(
        [row[row[:, 0].argsort()] for row in matsort])
    sino_cor = np.transpose(matsortback[:, :, 1])
    listxmiss = np.where(listmask > 0.0)[0]
    sinogram[:, listxmiss] = sino_cor[:, listxmiss]
    return sinogram


def remove_unresponsive_and_fluctuating_stripe(sinogram, snr, size):
    """
    Algorithm 6 in the paper. Remove unresponsive or fluctuating stripes by:
    locating stripes, correcting by interpolation.
    Angular direction is along the axis 0.

    Parameters
    ----------
    sinogram : float
        2D array.
    snr : float
        Ratio used to segment between useful information and noise
    size : int
        Window size of the median filter.

    Returns
    -------
    float
        2D array. Stripe-removed sinogram.
    """
    sinogram = np.copy(sinogram)  # Make it mutable
    (nrow, _) = sinogram.shape
    sinosmooth = np.apply_along_axis(uniform_filter1d, 0, sinogram, 10)
    listdiff = np.sum(np.abs(sinogram - sinosmooth), axis=0)
    nmean = np.mean(listdiff)
    listdiff_bck = median_filter(listdiff, size)
    listdiff_bck[listdiff_bck == 0.0] = nmean
    listfact = listdiff / listdiff_bck
    listmask = detect_stripe(listfact, snr)
    listmask = np.float32(binary_dilation(listmask, iterations=1))
    listmask[0:2] = 0.0
    listmask[-2:] = 0.0
    listx = np.where(listmask < 1.0)[0]
    listy = np.arange(nrow)
    matz = sinogram[:, listx]
    finter = interpolate.interp2d(listx, listy, matz, kind='linear')
    listxmiss = np.where(listmask > 0.0)[0]
    if len(listxmiss) > 0:
        matzmiss = finter(listxmiss, listy)
        sinogram[:, listxmiss] = matzmiss
    # Use algorithm 5 to remove residual stripes
    # sinogram = remove_large_stripe(sinogram, snr, size)
    return sinogram


def remove_all_stripe(sinogram, snr, la_size, sm_size, dim=1):
    """
    Remove all types of stripe artifacts by combining algorithm 6, 5, and 3.
    Angular direction is along the axis 0.

    Parameters
    ----------
    sinogram : float
        2D array.
    snr : float
        Ratio used to segment between useful information and noise.
    la_size : int
        Window size of the median filter to remove large stripes.
    sm_size : int
        Window size of the median filter to remove small-to-medium stripes.
    dim : {1, 2}, optional
        Dimension of the window.

    Returns
    -------
    float
        2D array. Stripe-removed sinogram.
    """
    sinogram = remove_unresponsive_and_fluctuating_stripe(
        sinogram, snr, la_size)
    sinogram = remove_large_stripe(sinogram, snr, la_size)
    sinogram = remove_stripe_based_sorting(sinogram, sm_size, dim)
    return sinogram

# ----------------------------------------------------------------------------
# Example of use for a sinogram with w x h ~ 2k x 2k pixel
# sinogram = remove_stripe_based_sorting(sinogram, 31)
# sinogram = remove_stripe_based_filtering(sinogram, 3, 31)
# sinogram = remove_stripe_based_fitting(sinogram, 1, 5, 20)
# sinogram = remove_unresponsive_and_fluctuating_stripe(sinogram1, 3, 81)
# sinogram = remove_large_stripe(sinogram1, 3, 81)
# sinogram = remove_all_stripe(sinogram, 3, 81, 31)
