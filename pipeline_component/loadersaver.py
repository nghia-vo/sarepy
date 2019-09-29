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
#============================================================================

"""
Module for I/O tasks:
- Load data from an image file (tif, png, jpeg).
- Save a 2D array to an image.
"""

import os
import numpy as np
from PIL import Image
import errno

def load_image(file_path):
    """
    Load data from an image.
    ---------
    Parameters - file_path: Path to the file.
    ---------
    Return:    - 2D array.
    """
    if ("\\" in file_path):
        raise ValueError(
            "Please use a file path following the Unix convention")
    mat = None
    try:
        mat = np.asarray(Image.open(file_path), dtype=np.float32)
    except IOError:
        print(("No such file or directory: {}").format(file_path))
        raise
    if len(mat.shape) > 2:
        axis_m = np.argmin(mat.shape)
        mat = np.mean(mat, axis=axis_m)
    return mat

def _create_folder(file_path):
    """
    Create folder if not exists.
    Parameter: file_path.
    """
    file_base = os.path.dirname(file_path)
    if not os.path.exists(file_base):
        try:
            os.makedirs(file_base)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

def _create_file_name(file_path):
    """
    Create file name to avoid overwriting.
    Parameter: File path.
    Return:    Updated file path.
    """
    file_base, file_ext = os.path.splitext(file_path)
    if os.path.isfile(file_path):
        nfile = 1
        check = True
        while check:
            name_add = '0000' + str(nfile)
            file_path = file_base + "_" + name_add[-4:] + file_ext
            if os.path.isfile(file_path):
                nfile = nfile + 1
            else:
                check = False
    return file_path

def save_image(file_path, mat, overwrite=True):
    """
    Save 2D data to an image.
    ---------
    Parameters: - file_path: Path to the file.
                - mat: 2D array.
                - overwrite: Overwrite the existing file if True.
    ---------
    Return:     - Updated file path.
    """
    if ("\\" in file_path):
        raise ValueError(
            "Please use a file path following the Unix convention")
    file_base, file_ext = os.path.splitext(file_path)
    if not ((file_ext == ".tif") or (file_ext == ".tiff")):
        mat = np.uint8(255*(mat-np.min(mat))/(np.max(mat)-np.min(mat)))
    _create_folder(file_path)
    if not overwrite:
        file_path = _create_file_name(file_path)
    image = Image.fromarray(mat)
    try:
        image.save(file_path)
    except IOError:
        print(("Couldn't write to file {}").format(file_path))
        raise
    return file_path
