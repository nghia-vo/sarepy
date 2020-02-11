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
import h5py
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from collections import OrderedDict
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


def _get_key(name, obj):
    """
    Find the key path to the dataset automatically.
    Use with Group.visititems() method to walk through hdf5 tree.
    """
    wanted_key = 'data'
    if isinstance(obj, h5py.Group):
        for key, val in list(obj.items()):
            if key == wanted_key:
                if isinstance(obj[key], h5py.Dataset):
                    key_path = obj.name + "/" + key
                    return key_path


def load_hdf_file(file_path, key_path=None, index=None, axis=0):
    """
    Load data from a hdf5 file.
    
    Parameters
    ----------
    file_path : str
        Path to the file
    key_path : str
        Key path to the dataset
    index : int 
        Values for slicing data. Can be integer, tuple or list, 
        e.g index=(start,stop,step) or index=(slice1, slice2, slice3,slice4).
    axis : int
        Slice direction
    
    Returns
    -------
    float
        2D array or 3D array.
    """
    ifile = None
    idata = None
    mat = None
    if ("\\" in file_path):
        raise ValueError(
            "Please use a file path following the Unix convention")
    try:
        ifile = h5py.File(file_path, 'r')
    except IOError:
        print(("Couldn't open file: {}").format(file_path))
        raise
    if key_path == None:
        key_path = ifile.visititems(_get_key)  # Find the key automatically
        if key_path == None:
            raise ValueError("Please provide the key path to the dataset!")
    check = key_path in ifile
    if not check:
        print(("Couldn't open object with the key path: {}").format(key_path))
        raise ValueError("!!! Wrong key !!!")
    idata = ifile[key_path]
    shape = idata.shape
    if (len(shape) < 2 or len(shape) > 3):
        raise ValueError("Require 2D or 3D dataset!")
    if len(shape) == 2:
        mat = np.asarray(idata)
    if len(shape) == 3:
        axis = np.clip(axis, 0, 3)
        (depth, height, width) = idata.shape
        if (index == None):
            mat = np.float32(idata[:, :, :])
        else:
            if type(index) == int:
                try:
                    if axis == 0:
                        mat = np.float32(idata[index, :, :])
                    elif axis == 1:
                        mat = np.float32(idata[:, index, :])
                    else:
                        mat = np.float32(idata[:, :, index])
                except ValueError:
                    print("Index out of range")
            if (type(index) == tuple) or (type(index) == list):
                if len(index) == 3:
                    starti = index[0]
                    stopi = index[1]
                    stepi = index[2]
                    list_index = list(range(starti, stopi, stepi))
                elif len(index) == 2:
                    starti = index[0]
                    stopi = index[1]
                    list_index = list(range(starti, stopi))
                else:
                    list_index = list(index)
                try:
                    if axis == 0:
                        mat = np.float32(idata[list_index, :, :])
                    elif axis == 1:
                        mat = np.float32(idata[:, list_index, :])
                    else:
                        mat = np.float32(idata[:, :, list_index])
                except ValueError:
                    print("Index out of range")
            if mat.shape[axis] == 1:
                mat = np.swapaxes(mat, axis, 0)[0]
            if mat.shape[axis] == 0:
                raise ValueError("Empty indices!")
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
    if ("\\" in file_path):
        raise ValueError(
            "Please use a file path following the Unix convention")
    file_base, file_ext = os.path.splitext(file_path)
    if not ((file_ext == ".tif") or (file_ext == ".tiff")):
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


def save_hdf_file(file_path, idata, key_path='entry', overwrite=True):
    """
    Write data to a hdf5 file.
    
    Parameters
    ----------
    file_path : str
        Path to the file.
    idata : data-type 
        Data to be saved.
    key_path : str 
        Key path to the dataset.
    overwrite : bool
        Overwrite the existing file if True.
    
    Returns
    -------
    str
        Updated file path.
    """
    if ("\\" in file_path):
        raise ValueError(
            "Please use a file path following the Unix convention")
    file_base, file_ext = os.path.splitext(file_path)
    if not ((file_ext == '.hdf') or (file_ext == '.h5')):
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
    grp = ofile.create_group(key_path)
    grp.create_dataset("data", data=idata)
    ofile.close()
    return file_path
