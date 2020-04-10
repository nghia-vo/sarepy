#============================================================================
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
# Description: Python codes for reconstructing tomographic data.
#============================================================================

"""
Module for reconstruction algorithms.
"""

import numpy as np
import astra
from scipy.ndimage import interpolation
# import tomopy


def cirle_mask(width, ratio):
    """
    Create a circle mask.
    
    Parameters
    ----------
    width : int
        Width of a square array.
    ratio : float 
        Ratio between the diameter of the mask and the width of the array.
    
    Returns
    -------
    float
        Square array.
    """
    mask = np.zeros((width, width), dtype=np.float32)
    center = width // 2
    radius = ratio * center
    y, x = np.ogrid[-center:width - center, -center:width - center]
    mask_check = x * x + y * y <= radius * radius
    mask[mask_check] = 1.0
    return mask


def recon_astra(sinogram, center, angles=None, ratio=1.0, method="FBP_CUDA", num_iter=1, filter="hann", pad=0):
    """
    Wrapper of reconstruction methods implemented in the astra toolbox package.
    https://www.astra-toolbox.com/docs/algs/index.html
    
    Parameters
    ----------
    sinogram : float 
        2D tomographic data.
    center : float
        Center of rotation.
    angles : float 
        1D array. Tomographic angles in radian.
    ratio : float
        To apply a circle mask to the reconstructed image.
    method : str 
        Reconstruction algorithms. for CPU: 'FBP', 'SIRT', 'SART', 'ART', 
        'CGLS'. for GPU: 'FBP_CUDA', 'SIRT_CUDA', 'SART_CUDA', 'CGLS_CUDA'.
    num_iter : int 
        Number of iterations if using iteration methods.
    filter : str 
        Apply filter if using FBP method. Options: 'hamming', 'hann', 
        'lanczos', 'kaiser', 'parzen',...
    pad : int 
        Padding to reduce the side effect of FFT.
    
    Returns
    -------
    float
        Square array.        
    """
    if pad > 0:
        sinogram = np.pad(sinogram, ((0, 0), (pad, pad)), mode='edge')
        center = center + pad
    (nrow, ncol) = sinogram.shape
    if angles is None:
        angles = np.linspace(0.0, 180.0, nrow) * np.pi / 180.0
    proj_geom = astra.create_proj_geom('parallel', 1, ncol, angles)
    vol_geom = astra.create_vol_geom(ncol, ncol)
    cen_col = (ncol - 1.0) / 2.0
    shift = cen_col - center
    sinogram = interpolation.shift(sinogram, (0, shift), mode='nearest')
    sino_id = astra.data2d.create('-sino', proj_geom, sinogram)
    rec_id = astra.data2d.create('-vol', vol_geom)
    if "CUDA" not in method:
        proj_id = astra.create_projector('line', proj_geom, vol_geom)
    cfg = astra.astra_dict(method)
    cfg['ProjectionDataId'] = sino_id
    cfg['ReconstructionDataId'] = rec_id
    if "CUDA" not in method:
        cfg['ProjectorId'] = proj_id
    if method == "FBP_CUDA" or method =="FBP":
        cfg["FilterType"] = filter
    alg_id = astra.algorithm.create(cfg)
    astra.algorithm.run(alg_id, num_iter)
    rec = astra.data2d.get(rec_id)
    astra.algorithm.delete(alg_id)
    astra.data2d.delete(sino_id)
    astra.data2d.delete(rec_id)
    if pad > 0:
        rec = rec[pad:-pad, pad:-pad]
    if not (ratio is None):
        rec = rec * cirle_mask(rec.shape[0], ratio)
    return rec


def recon_gridrec(sinogram, center, angles=None, ratio=1.0):
    """
    Wrapper of the gridrec method implemented in the tomopy package.
    https://tomopy.readthedocs.io/en/latest/api/tomopy.recon.algorithm.html
     
    Parameters
    ----------
    sinogram : float 
        2D tomographic data.
    center : float
        Center of rotation.
    angles : float
        1D array. Tomographic angles in radian.
    ratio : float
        To apply a circle mask to the reconstructed image.
     
    Returns
    -------
    float
        Square array.
    """
    (nrow, ncol) = sinogram.shape
    if angles is None:
        angles = np.linspace(0.0, 180.0, nrow) * np.pi / 180.0
    sinogram = np.expand_dims(sinogram, 1)
    recont = tomopy.recon(sinogram, angles, center=center, algorithm='gridrec')
    recont = tomopy.circ_mask(recont, axis=0, ratio=ratio)
    return recont[0]
