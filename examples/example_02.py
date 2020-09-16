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
# Description: Example to show how to remove ring artifacts from a neutron
# image.
#============================================================================

import sys
sys.path.insert(0, "C:/sarepy-master/")

import os
import numpy as np

import sarepy.losa.loadersaver as losa
import sarepy.prep.autocentering as cen
import sarepy.reco.reconstruction as rec

import sarepy.prep.stripe_removal_original as srm1


file_path = "C:/sarepy-master/sarepy/data/sinogram_360_neutron_image.tif"
out_path = "C:/home/reconstruction/"

file_path = file_path.replace("\\", "/")
_, file_name = os.path.split(file_path)

# Load image
# This is a 360-degree sinogram, 16-bit, neutron image.
sinogram = losa.load_image(file_path)
# Apply normalization due to the image is 16-bit.
nmean = np.mean(sinogram[:, 0:30])
sinogram = sinogram / nmean
# Replace values <=0.0
nmean = np.mean(sinogram)
sinogram[sinogram <= 0.0] = nmean

# Calculate center of rotation
center = 0.0
(nrow, ncol) = sinogram.shape
if center == 0.0:
    search_range = 50
    print("Calculate center of rotation.....")
#     center = cen.find_center_vo(
#         sinogram, ncol // 2 - search_range, ncol // 2 + search_range)
    # Uncomment the following lines if use a 360-degree sinogram
    center = cen.find_center_vo(sinogram[0:nrow // 2 + 1],
                                ncol // 2 - search_range,
                                ncol // 2 + search_range)
    print("Center of rotation ---> {}".format(center))
ratio = (min(center, abs(ncol - center))) / (ncol * 0.5)

# Reconstruction without ring artifact removal
list_angles = np.linspace(0.0, 360.0, nrow) * np.pi / 180.0
rec_image = rec.recon_gridrec(
    -np.log(sinogram), center, angles=list_angles, ratio=ratio)
losa.save_image(out_path + "/rec_gridrec_before_" + file_name, rec_image)
# There're streak artifacts introduced by the gridrec method
# which don't shown if using the FBP method.
rec_image = rec.recon_astra(
    -np.log(sinogram), center, angles=list_angles, ratio=ratio)
losa.save_image(out_path + "/rec_fbp_before_" + file_name, rec_image)

# Apply ring artifact methods.
# Parameters used are much smaller than these used in example_01.py.
sinogram = srm1.remove_unresponsive_and_fluctuating_stripe(sinogram, 3.0, 31)
sinogram = srm1.remove_stripe_based_sorting(sinogram, 5)


# Reconstruction with ring artifact removal
rec_image = rec.recon_gridrec(
    -np.log(sinogram), center, angles=list_angles, ratio=ratio)
losa.save_image(out_path + "/rec_after_" + file_name, rec_image)
# Streak artifacts are gone no matter what reconstruction method is used
