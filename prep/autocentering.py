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
# Description: Python implementation of a method for calculating center of
# rotation in tomography, Nghia T. Vo, Michael Drakopoulos, Robert C. Atwood,
# and Christina Reinhard, "Reliable method for calculating the center of
# rotation in parallel-beam tomography," Opt. Express 22, 19078-19086 (2014).
# Publication date:
#============================================================================

"""
Module for calculating center of rotation in parallel-beam tomography.

"""

import numpy as np
import pyfftw
import pyfftw.interfaces.scipy_fftpack as fft
import scipy.ndimage as ndi


def make_mask(height, width, radius):
    """
    Make a binary mask to select coefficients outside the double-wedge region.
    Eq.(3) in https://doi.org/10.1364/OE.22.019078
    
    Parameters
    ----------
    height : int
        Image height.
    width : int
        Image width.
    radius : int 
        Radius of an object, in pixel unit.
    
    Returns
    -------
    float
        2D binary mask.
    """
    du = 1.0 / width
    dv = (height - 1.0) / (height * 2.0 * np.pi)
    rdrop = min(20, np.int16(np.ceil(0.05 * height)))
    cen_hei = np.int16(np.ceil(height / 2.0)) - 1
    cen_wid = np.int16(np.ceil(width / 2.0)) - 1
    mask = np.zeros((height, width), dtype=np.float32)
    for i in range(height):
        num1 = np.int16(np.ceil(((i - cen_hei) * dv / radius) / du))
        (p1, p2) = np.int16(np.clip(
            np.sort((-num1 + cen_wid, num1 + cen_wid)), 0, width - 1))
        mask[i, p1:p2 + 1] = 1.0
    mask[cen_hei - rdrop:cen_hei + rdrop + 1, :] = 0.0
    mask[:, cen_wid - 1:cen_wid + 2] = 0.0
    return mask


def coarse_search_based_integer_shift(
        sino_0_180, start_cor=None, stop_cor=None, ratio=0.5):
    """
    Find center of rotation (CoR) using integer shifts of a 180-360 sinogram.
    The 180-360 sinogram is made by flipping horizontally the 0-180 sinogram.
    Angular direction is along the axis 0.
    Note:
    1-pixel shift of the 180-360 sinogram is equivalent to 0.5-pixel shift \
    of CoR. Auto-search is limited to the range of [width/4; width - width/4].
    
    Parameters
    ----------
    sino_0_180 : float 
        Sinogram in the angle range of [0;180].
    start_cor : int 
        Starting point for searching CoR.
    stop_cor : int
        Ending point for searching CoR.
    ratio : float 
        Ratio between a sample and the width of the sinogram. Default value 
        of 0.5 works in most cases (even when the sample is larger than the 
        field_of_view).
    
    Returns
    -------
    float
        Center of rotation with haft-pixel accuracy.
    """
    # Denoising. Should not be used for a downsampled sinogram
    # sino_0_180 = ndi.gaussian_filter(sino_0_180, (3,1), mode='reflect')
    (nrow, ncol) = sino_0_180.shape
    if start_cor is None:
        start_cor = ncol // 4
    if stop_cor is None:
        stop_cor = ncol - ncol // 4 - 1
    start_cor, stop_cor = np.sort((start_cor, stop_cor))
    start_cor = np.int16(np.clip(start_cor, 0, ncol - 1))
    stop_cor = np.int16(np.clip(stop_cor, 0, ncol - 1))
    cen_fliplr = (ncol - 1.0) / 2.0
    sino_180_360 = np.fliplr(sino_0_180)
    comp_sino = np.flipud(sino_0_180)  # Used to avoid local minima
    list_cor = np.arange(start_cor, stop_cor + 0.5, 0.5)
    list_metric = np.zeros(len(list_cor), dtype=np.float32)
    mask = make_mask(2 * nrow, ncol, 0.5 * ratio * ncol)
    for i, cor in enumerate(list_cor):
        shift = np.int16(2.0 * (cor - cen_fliplr))
        sino_shift = np.roll(sino_180_360, shift, axis=1)
        if shift >= 0:
            sino_shift[:, :shift] = comp_sino[:, :shift]
        else:
            sino_shift[:, shift:] = comp_sino[:, shift:]
        mat1 = np.vstack((sino_0_180, sino_shift))
        list_metric[i] = np.mean(
            np.abs(np.fft.fftshift(fft.fft2(mat1))) * mask)
    min_pos = np.argmin(list_metric)
    cor = list_cor[min_pos]
    return cor


def fine_search_based_subpixel_shift(
    sino_0_180, start_cor, search_radius=4, search_step=0.25,
        ratio=0.5):
    """
    Find center of rotation (CoR) with sub-pixel accuracy by shifting a\
    180-360 sinogram around the coarse CoR. The 180-360 sinogram is made\
    by flipping horizontally a 0-180 sinogram.
    Angular direction is along the axis 0.
    
    Parameters
    ----------
    sino_0_180 : float 
        Sinogram in the angle range of [0;180].
    start_cor : float 
        Starting point for searching CoR.
    search_radius : float
        Searching range = (start_cor +/- search_radius)
    search_step : float 
        Searching step.
    ratio : float
        Ratio between the sample and the width of the sinogram. Default value 
        of 0.5 works in most cases (even when a sample is larger than the 
        field_of_view).
    
    Returns
    -------
    float
        Center of rotation.
    """
    # Denoising
    sino_0_180 = ndi.gaussian_filter(sino_0_180, (2, 2), mode='reflect')
    (nrow, ncol) = sino_0_180.shape
    sino_180_360 = np.fliplr(sino_0_180)
    search_radius = np.clip(np.abs(search_radius), 1, ncol // 10 - 1)
    search_step = np.clip(np.abs(search_step), 0.1, 1.1)
    start_cor = np.clip(start_cor, search_radius, ncol - search_radius - 1)
    cen_fliplr = (ncol - 1.0) / 2.0
    list_cor = start_cor + np.arange(
        -search_radius, search_radius + search_step, search_step)
    comp_sino = np.flipud(sino_0_180)  # Used to avoid local minima
    list_metric = np.zeros(len(list_cor), dtype=np.float32)
    mask = make_mask(2 * nrow, ncol, 0.5 * ratio * ncol)
    for i, cor in enumerate(list_cor):
        shift = 2.0 * (cor - cen_fliplr)
        sino_shift = ndi.interpolation.shift(
            sino_180_360, (0, shift), order=3, prefilter=True)
        if shift >= 0:
            shift_int = np.int16(np.ceil(shift))
            sino_shift[:, :shift_int] = comp_sino[:, :shift_int]
        else:
            shift_int = np.int16(np.floor(shift))
            sino_shift[:, shift_int:] = comp_sino[:, shift_int:]
        mat1 = np.vstack((sino_0_180, sino_shift))
        list_metric[i] = np.mean(
            np.abs(np.fft.fftshift(fft.fft2(mat1))) * mask)
    min_pos = np.argmin(list_metric)
    cor = list_cor[min_pos]
    return cor


def _downsample(image, dsp_fact0, dsp_fact1):
    """
    Downsample an image by averaging.
    
    Parameters
    ----------
    image : float
        2D array.
    dsp_fact0 : int
        Downsampling factor along axis 0.
    dsp_fact1 : int
        Downsampling factor along axis 1.
    
    Returns
    -------
    float
        2D array. Downsampled image.
    """
    (height, width) = image.shape
    dsp_fact0 = np.clip(np.int16(dsp_fact0), 1, height // 2)
    dsp_fact1 = np.clip(np.int16(dsp_fact1), 1, width // 2)
    height_dsp = height // dsp_fact0
    width_dsp = width // dsp_fact1
    image_dsp = image[0:dsp_fact0 * height_dsp, 0:dsp_fact1 * width_dsp]
    image_dsp = image_dsp.reshape(
        height_dsp, dsp_fact0, width_dsp, dsp_fact1).mean(-1).mean(1)
    return image_dsp


def find_center_vo(
        sinogram, start_cor=None, stop_cor=None, step_cor=0.25, fine_srange=3,
        ratio=0.5, dsp=True):
    """
    Find center of rotation.
    Note: Downsampling is a trade-off solution to reduce computational cost.
    There's a risk of being stuck in a local minimum by doing that.
    
    Parameters 
    ----------
    sinogram : float
        Sinogram in the angle range of [0;180].
    start_cor : float 
        Starting point for searching CoR.
    stop_cor : float
        Ending point for searching CoR.
    step_cor : float
        Sub-pixel accuracy of searched CoR.
    fine_srange : float
        Searching range for finding CoR with sub-pixel accuracy.
    ratio : float
        Ratio between the sample and the width of the sinogram. Default value 
        of 0.5 works in most cases (even when a sample is larger than the 
        field_of_view).
    dsp : bool 
        Enable/disable downsampling.
    
    Returns
    -------
    float
        Center of rotation.
    """
    (nrow, ncol) = sinogram.shape
    dsp_row = 1
    dsp_col = 1
    if dsp is True:
        if ncol > 2000:
            dsp_col = 4
        if nrow > 2000:
            dsp_row = 2
        # Denoising before downsampling
        sinogram1 = ndi.gaussian_filter(sinogram, (3, 1), mode='reflect')
        sino_dsp = _downsample(sinogram1, dsp_row, dsp_col)
        fine_srange = max(fine_srange, dsp_col)
        if start_cor != None:
            start_cor = np.int16(np.floor(1.0 * start_cor / dsp_col))
        if stop_cor != None:
            stop_cor = np.int16(np.ceil(1.0 * stop_cor / dsp_col))
        off_set = 0.5 * dsp_col
        raw_cor = coarse_search_based_integer_shift(
            sino_dsp, start_cor, stop_cor, ratio)
        fine_cor = fine_search_based_subpixel_shift(
            sinogram, raw_cor * dsp_col + off_set, fine_srange, step_cor, ratio)
    else:
        raw_cor = coarse_search_based_integer_shift(
            sinogram, start_cor, stop_cor, ratio)
        fine_cor = fine_search_based_subpixel_shift(
            sinogram, raw_cor, fine_srange, step_cor, ratio)
    return fine_cor

# #---------------------------------------------------------------------------
# #Example of use
# import multiprocessing
# num_cpu = multiprocessing.cpu_count()
# pyfftw.config.NUM_THREADS = num_cpu - 1
# pyfftw.interfaces.cache.enable()
# pyfftw.interfaces.cache.set_keepalive_time(10.0)
# #---------------------------------------------------------------------------
# #Example 1: Use default parameters
# cor = find_center_vo(sinogram)
# #----------------------------------------------------------------------------
# Example 2: Searching around the horizontal center of a sinogram.
# search_range = 50
# (nrow, ncol) = sinogram.shape
# cor = find_center_vo(sinogram, ncol//2-search_range, ncol//2+search_range)
