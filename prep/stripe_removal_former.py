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
# Description: Implementations of former stripe removal methods
# ->  Normalization-based method
# [1] M. Rivers, "Tutorial introduction to X-ray computed microtomography data
#     processing,"
#     http://www.mcs.anl.gov/research/projects/X-ray-cmt/rivers/tutorial.html
# ->  Regularization-based method
# [2] S. Titarenko, P. J. Withers, and A. Yagola, "An analytical formula for
#     ring artefact suppression in X-ray tomography,"
#     Appl. Math. Lett. 23(12), 1489-1495 (2010).
# ->  FFT-based method
# [3] C. Raven, "Numerical removal of ring artifacts in microtomography,"
#     Rev. Sci. Instrum. 69(8), 2978-2980 (1998).
# ->  Wavelet-FFT-based method
# [4] B. Munch, P. Trtik, F. Marone, and M. Stampanoni, "Stripe and ring
#     artifact removal with combined wavelet-Fourier filtering,"
#     Opt. Express 17(10), 8567-8591 (2009).
# Publication date:
#============================================================================

import numpy as np
from scipy.ndimage import gaussian_filter
import pywt
import pyfftw.interfaces.scipy_fftpack as fft
# from scipy.fftpack import fft, ifft, fft2, ifft2


def remove_stripe_based_normalization(sinogram, sigma, num_chunk):
    """
    Remove stripes using the method in [1].
    Angular direction is along the axis 0.
    ---------
    Parameters: - sinogram: 2D array.
                - sigma: sigma of the Gaussian window.
                - num_chunk: number of chunks of rows.
    ---------
    Return:     - stripe-removed sinogram.
    """
    (nrow, _) = sinogram.shape
    listindex = np.array_split(np.arange(nrow), num_chunk)
    for pos in listindex:
        bindex = pos[0]
        eindex = pos[-1] + 1
        listmean = np.mean(sinogram[bindex:eindex], axis=0)
        list_filtered = gaussian_filter(listmean, sigma)
        listcoe = list_filtered - listmean
        matcoe = np.tile(listcoe, (eindex - bindex, 1))
        sinogram[bindex:eindex, :] = sinogram[bindex:eindex, :] + matcoe
    return sinogram


def calculate_reg_mat(width, alpha):
    """
    Calculate coefficients used for the regularization-based method.
    Eq. (7) in [2].
    ---------
    Parameters: - width: width of a square array.
                - alpha: regularization parameter.
    ---------
    Return:     - square array.
    """
    tau = 2.0 * np.arcsinh(np.sqrt(alpha) * 0.5)
    ilist = np.arange(0, width)
    jlist = np.arange(0, width)
    matjj, matii = np.meshgrid(jlist, ilist)
    mat1 = np.abs(matii - matjj)
    mat2 = matii + matjj
    mat1a = np.cosh((width - 1 - mat1) * tau)
    mat2a = np.cosh((width - mat2) * tau)
    matcoe = - \
        (np.tanh(0.5 * tau) / (alpha * np.sinh(width * tau))) * (mat1a + mat2a)
    return matcoe


def remove_stripe_based_regularization(sinogram, alpha, num_chunk):
    """
    Remove stripes using the method in [2].
    Angular direction is along the axis 0.
    ---------
    Parameters: - sinogram: 2D array.
                - alpha: regularization parameter.
                - num_chunk: number of chunks of rows.
    ---------
    Return:     - stripe-removed sinogram.
    """
    sinogram = -np.log(sinogram)
    (nrow, ncol) = sinogram.shape
    sijmat = calculate_reg_mat(ncol, alpha)
    listindex = np.array_split(np.arange(nrow), num_chunk)
    listgrad = np.zeros(ncol, dtype=np.float32)
    matgrad = np.zeros((ncol, ncol), dtype=np.float32)
    for pos in listindex:
        bindex = pos[0]
        eindex = pos[-1] + 1
        listmean = np.mean(sinogram[bindex:eindex], axis=0)
        listgrad[1:-1] = - np.diff(listmean, 2)
        listgrad[0] = listmean[0] - listmean[1]
        listgrad[-1] = listmean[-1] - listmean[-2]
        matgrad[:] = listgrad
        listcoe = np.sum(matgrad * sijmat, axis=1)
        matcoe = np.tile(listcoe, (eindex - bindex, 1))
        sinogram[bindex:eindex, :] = sinogram[bindex:eindex, :] + matcoe
    return np.exp(-sinogram)


def create_2d_window(width, height, u0, v0, n):
    """
    Create a 2d window used for the fft-based method
    ---------
    Parameters: - height, width: shape of the window.
                - u0, n: to define the shape of 1D Butterworth
                        low-pass filter.
                - v0: number of rows (* 2) to be applied the filter.
    ---------
    Return:     - 2D window.
    """
    centerc = np.ceil(width / 2.0) - 1.0
    centerr = np.int16(np.ceil(height / 2.0) - 1)
    listx = np.arange(width) - centerc
    window1d = 1.0 / (1.0 + np.power(listx / u0, 2 * n))
    row1 = centerr - np.int16(v0)
    row2 = centerr + np.int16(v0) + 1
    window2d = np.ones((height, width), dtype=np.float32)
    window2d[row1:row2] = window1d
    return window2d


def remove_stripe_based_fft(sinogram, u0, n, v0, pad):
    """
    Remove stripes using the method in [3].
    Angular direction is along the axis 0.
    ---------
    Parameters: - sinogram: 2D array.
                - u0, n: to define the shape of 1D Butterworth
                        low-pass filter.
                - v0: number of rows (* 2) to be applied the filter.                
                - pad: padding for FFT
    ---------
    Return:     - stripe-removed sinogram.
    """
    if pad > 0:
        sinogram = np.pad(sinogram, ((pad, pad), (0, 0)), mode='mean')
        sinogram = np.pad(sinogram, ((0, 0), (pad, pad)), mode='edge')
    (nrow, ncol) = sinogram.shape
    window2d = create_2d_window(ncol, nrow, u0, v0, n)
    sinogram = fft.ifft2(
        np.fft.ifftshift(np.fft.fftshift(fft.fft2(sinogram)) * window2d))
    return np.real(sinogram[pad:nrow - pad, pad:ncol - pad])


def remove_stripe_based_wavelet_fft(sinogram, level, sigma, order, pad):
    """
    Remove stripes using the method in [4].
    Angular direction is along the axis 0.
    Code adapted from tomopy source code https://github.com/tomopy/tomopy
    with a small improvement of using different ways of padding to
    reduce the side effect of the FFT.
    Angular direction is along the axis 0.
    ---------
    Parameters: - sinogram: 2D array.
                - level: wavelet decomposition level.
                - sigma: Damping parameter. Larger is stronger.
                - order: order of the the Daubechies wavelets.
                - pad: padding for FFT
    ---------
    Return:     - stripe-removed sinogram.
    """
    (nrow, ncol) = sinogram.shape
    if pad > 0:
        sinogram = np.pad(sinogram, ((pad, pad), (0, 0)), mode='mean')
        sinogram = np.pad(sinogram, ((0, 0), (pad, pad)), mode='edge')
    # Wavelet decomposition
    cH = []
    cV = []
    cD = []
    waveletname = "db" + str(order)
    for j in range(level):
        sinogram, (cHt, cVt, cDt) = pywt.dwt2(sinogram, waveletname)
        cH.append(cHt)
        cV.append(cVt)
        cD.append(cDt)
    for j in range(level):
        fcV = np.fft.fftshift(np.fft.fft2(cV[j]))
        nrow1, ncol1 = fcV.shape
        y_hat = (np.arange(-nrow1, nrow1, 2, dtype=np.float32) + 1) / 2
        damp = 1 - np.exp(-np.power(y_hat, 2) / (2 * np.power(sigma, 2)))
        fcV = np.multiply(fcV, np.transpose(np.tile(damp, (ncol1, 1))))
        cV[j] = np.real(np.fft.ifft2(np.fft.ifftshift(fcV)))
    # Wavelet reconstruction.
    for j in range(level)[::-1]:
        sinogram = sinogram[0:cH[j].shape[0], 0:cH[j].shape[1]]
        sinogram = pywt.idwt2((sinogram, (cH[j], cV[j], cD[j])), waveletname)
    if pad > 0:
        sinogram = sinogram[pad:-pad, pad:-pad]
    return sinogram[:nrow, :ncol]
