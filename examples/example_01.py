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
# Description: Example to show how to load a tomographic data,
# remove stripe artifacts, and perform reconstruction.
#============================================================================

import sys
sys.path.insert(0, "C:/sarepy-master/")

import os
import numpy as np

import sarepy.losa.loadersaver as losa
import sarepy.prep.autocentering as cen
import sarepy.reco.reconstruction as rec

import sarepy.prep.stripe_removal_original as srm1
import sarepy.prep.stripe_removal_improved as srm2
import sarepy.prep.stripe_removal_former as srm3


file_path = "C:/sinogram/sinogram.tif"
out_path = "C:/reconstruction/"

file_path = file_path.replace("\\", "/")
_, file_name = os.path.split(file_path)

# Load a 180-degree sinogram.
# If use a 360-degree sinogram, make sure to update the angles
# and use haft of the sinogram for finding the center of rotation.

sinogram = losa.load_image(file_path)

# Calculate center of rotation
center = 0.0
(nrow, ncol) = sinogram.shape
if center == 0.0:
    search_range = 50
    print("Calculate center of rotation.....")
    center = cen.find_center_vo(
        sinogram, ncol // 2 - search_range, ncol // 2 + search_range)
    # Uncomment the following lines if use a 360-degree sinogram
    # center = cen.find_center_vo(sinogram[0:nrow //2 + 1],
    #                             ncol // 2 - search_range,
    #                             ncol // 2 + search_range)
    print("Center of rotation ---> {}".format(center))
ratio = (min(center, abs(ncol - center))) / (ncol * 0.5)

# Reconstruction without ring artifact removal
rec_image = rec.recon_astra(
    -np.log(sinogram), center, ratio=ratio, pad=50)
# rec_image = rec.recon_gridrec(-np.log(sinogram), center, ratio=ratio)

# If use a 360-degree sinogram
# list_angles = np.linspace(0.0, 360.0, nrow) * np.pi / 180.0
# rec_image = rec.recon_astra(
#     -np.log(sinogram), center, angles=list_angles, ratio=ratio, pad=50)
# rec_image = rec.recon_gridrec(
#     -np.log(sinogram), center, angles=list_angles, ratio=ratio)

losa.save_image(out_path + "/rec_before_correction_" + file_name, rec_image)

# Apply ring artifact methods.
#==========================================================================
# Parameters in the following examples are chosen for sinograms under /data
# Making sure that you change them corresponding to your data, start from
# small parameters.
#==========================================================================
# ->>> Uncomment one of following methods

# Using original methods

sinogram = srm1.remove_all_stripe(sinogram, 3.0, 81, 31)
# sinogram = srm1.remove_unresponsive_and_fluctuating_stripe(sinogram, 3.0, 81)
# sinogram = srm1.remove_stripe_based_sorting(sinogram, 31)
# sinogram = srm1.remove_stripe_based_fitting(sinogram, 2, 10, 60)
# sinogram = srm1.remove_stripe_based_filtering(sinogram, 3, 31)
# sinogram = srm1.remove_large_stripe(sinogram, 3.0, 81)


# Using improved methods

# sinogram = srm2.remove_stripe_based_2d_filtering_sorting(sinogram, 10, 31)
# sinogram = srm2.remove_stripe_based_filtering_sorting(sinogram, 3, 31)
# sinogram = srm2.remove_stripe_based_sorting_fitting(sinogram, 1, 10, 60)
# sinogram = srm2.remove_stripe_based_interpolation(sinogram, 3.0, 31)

# Using former methods

# sinogram = srm3.remove_stripe_based_fft(sinogram, 20, 8, 1, 100)
# sinogram = srm3.remove_stripe_based_normalization(sinogram, 15, 1)
# sinogram = srm3.remove_stripe_based_regularization(sinogram, 0.005, 1)
# sinogram = srm3.remove_stripe_based_wavelet_fft(sinogram, 5, 1.0, 9, 100)

# Reconstruction with ring artifact removal
rec_image = rec.recon_astra(
    -np.log(sinogram), center, ratio=ratio, pad=50)
# rec_image = rec.recon_gridrec(-np.log(sinogram), center, ratio=ratio)

# If use a 360-degree sinogram
# list_angles = np.linspace(0.0, 360.0, nrow) * np.pi / 180.0
# rec_image = rec.recon_astra(
#     -np.log(sinogram), center, angles=list_angles, ratio=ratio, pad=50)
# rec_image = rec.recon_gridrec(
#     -np.log(sinogram), center, angles=list_angles, ratio=ratio)

losa.save_image(out_path + "/rec_after_correction_" + file_name, rec_image)
