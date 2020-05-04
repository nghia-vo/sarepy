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
# Description: Python codes for loading and saving data.
#============================================================================

"""
Module for I/O tasks:
- Load data from an image file (tif, png, jpeg) or a hdf file.
- Save a 2D array as a tif image or 2D, 3D array to a hdf file.

"""

import os
import glob
import h5py
import numpy as np
from PIL import Image
import errno


def load_image(file_path):
    """
    Load data from an image.

    Parameters
    ----------
    file_path : str
        Path to the file.

    Returns
    -------
    float
        2D array.
    """
    if "\\" in file_path:
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


def load_hdf_file(file_path, key_path):
    """
    Create a hdf object

    Parameters
    ----------
    file_path : str
        Path to the file
    key_path : str
        Key path to the dataset

    Returns
    -------
    object
        hdf object.
    """
    try:
        ifile = h5py.File(file_path, 'r')
    except IOError:
        print(("Couldn't open file: {}").format(file_path))
        raise
    check = key_path in ifile
    if not check:
        print("Couldn't open object with the key path: {}".format(key_path))
        raise ValueError("!!! Wrong key !!!")
    return ifile[key_path]


def _create_folder(file_path):
    """
    Create folder if not exists.

    Parameters
    ----------
    file_path : str
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

    Parameters
    ----------
    file_path : str

    Returns
    -------
    str
        Updated file path.
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

    Parameters
    ----------
    file_path : str 
        Path to the file.
    mat : int or float 
        2D array.
    overwrite : bool 
        Overwrite the existing file if True.

    Returns
    -------
    str
        Updated file path.
    """
    if "\\" in file_path:
        raise ValueError(
            "Please use a file path following the Unix convention")
    _, file_ext = os.path.splitext(file_path)
    if not (file_ext == ".tif") or (file_ext == ".tiff"):
        mat = np.uint8(255 * (mat - np.min(mat)) / (np.max(mat) - np.min(mat)))
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


def open_hdf_stream(file_path, data_shape, key_path='entry', data_type='float32', overwrite=True):
    """
    Write data to a hdf5 file.

    Parameters
    ----------
    file_path : str
        Path to the file.
    data_shape : tuple of int
        Shape of the data.
    key_path : str
        Key path to the dataset.
    data_type: str
        Type of data.
    overwrite : bool
        Overwrite the existing file if True.

    Returns
    -------
    object
        hdf object.
    """
    file_base, file_ext = os.path.splitext(file_path)
    if not (file_ext == '.hdf' or file_ext == '.h5'):
        file_ext = '.hdf'
    file_path = file_base + file_ext
    _create_folder(file_path)
    if not overwrite:
        file_path = _create_file_name(file_path)
    ofile = None
    try:
        ofile = h5py.File(file_path, 'w')
    except IOError:
        print(("Couldn't write file: {}").format(file_path))
        raise
    data_out = ofile.create_dataset(key_path, data_shape, dtype=data_type)
    return data_out


def create_folder_name(folder_path, name_prefix, zero_prefix=4):
    """
    Create folder name. E.g: Rec_0001, Rec_0002...

    Parameters
    ----------
    folder_path : str
        Path to the parent folder.
    name_prefix : str
        Name prefix
    zero_prefix : int

    Returns
    -------
    str
        Name of the folder.
    """
    scan_name_prefix = name_prefix + "_"
    num_folder_exist = len(
        glob.glob(folder_path + "/" + scan_name_prefix + "*"))
    num_folder_new = num_folder_exist + 1
    name_tmp = "00000" + str(num_folder_new)
    scan_name = scan_name_prefix + name_tmp[-zero_prefix:]
    while os.path.isdir(folder_path + "/" + scan_name):
        num_folder_new = num_folder_new + 1
        name_tmp = "00000" + str(num_folder_new)
        scan_name = scan_name_prefix + name_tmp[-zero_prefix:]
    return scan_name


def find_file(path):
    """
    Search file

    Parameters
    ----------
    path : str
        Path and/or pattern to find files.

    Returns
    -------
    str or list of str
        List of files.
    """
    file_path = glob.glob(path)
    if len(file_path) == 0:
        raise ValueError("No such file {}".format(file_path))
    if len(file_path) == 1:
        file_path = file_path[0]
    return file_path
