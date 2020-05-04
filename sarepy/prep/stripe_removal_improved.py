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
from scipy.signal.windows import gaussian
#from scipy.fftpack import fft, ifft, fft2, ifft2
import pyfftw.interfaces.scipy_fftpack as fft

from sarepy.prep.stripe_removal_original import remove_stripe_based_sorting
from sarepy.prep.stripe_removal_original import remove_stripe_based_fitting
from sarepy.prep.stripe_removal_original import _2d_filter


def remove_stripe_based_filtering_sorting(sinogram, sigma, size, dim=1):
    """
    Combination of algorithm 2 and algorithm 3 in [1].
    Removing stripes using the filtering and sorting technique.
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
        # sinosmooth[i] = np.real(ifft(fft(
        #     sinolist * listsign) * window) * listsign)[pad:ncol-pad]
        sinosmooth[i] = np.real(
            fft.ifft(fft.fft(sinolist * listsign) * window) * listsign)[pad:ncol - pad]
    sinosharp = sinogram - sinosmooth
    sinosmooth_cor = np.transpose(
        remove_stripe_based_sorting(np.transpose(sinosmooth), size, dim))
    return np.transpose(sinosmooth_cor + sinosharp)


def remove_stripe_based_sorting_fitting(sinogram, order, sigmax, sigmay):
    """
    Combination of algorithm 3 and 1 in [1].
    Remove stripes using the sorting and fitting technique.
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
    sinogram = np.transpose(sinogram)
    (nrow, ncol) = sinogram.shape
    listindex = np.arange(0.0, ncol, 1.0)
    matindex = np.tile(listindex, (nrow, 1))
    matcomb = np.asarray(np.dstack((matindex, sinogram)))
    matsort = np.asarray(
        [row[row[:, 1].argsort()] for row in matcomb])
    sinosort = matsort[:, :, 1]
    sino_cor = np.transpose(
        remove_stripe_based_fitting(np.transpose(sinosort), order, sigmax, sigmay))
    matsort[:, :, 1] = sino_cor
    matsortback = np.asarray(
        [row[row[:, 0].argsort()] for row in matsort])
    return np.transpose(matsortback[:, :, 1])


def remove_stripe_based_2d_filtering_sorting(sinogram, sigma, size, dim=1):
    """
    Algorithm used in section 3.3.4 in [2].
    Remove stripes using a 2D low-pass filter and the sorting-based technique.
    Angular direction is along the axis 0.

    Parameters
    ---------
    sinogram : float
        2D array.
    sigma : int
        Sigma of the Gaussian window.
    size : int
        Window size of the median filter.
    dim : {1, 2}, optional
        Dimension of the window.

    Returns
    -------
    float
         2D array. Stripe-removed sinogram.
    """
    pad = 150
    sinosmooth = _2d_filter(sinogram, sigma, sigma, pad)
    sinosharp = sinogram - sinosmooth
    sino_cor = remove_stripe_based_sorting(sinosharp, size, dim)
    return sinosmooth + sino_cor
