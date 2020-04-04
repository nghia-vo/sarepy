# Sarepy
(S)tripe (A)rtifacts (RE)moval in (PY)thon
## Numerical techniques for eliminating ring artifacts in X-ray micro-tomography

**Sarepy** is the Python implementations of methods used for removing ring artifacts in tomography.
 These methods work in sinogram space where artifacts appear as straight lines or stripe artifacts.
 The codes are developed from the original implementations of methods published in Optics Epxress,
 Nghia T. Vo, Robert C. Atwood, and Michael Drakopoulos, *"Superior techniques for eliminating ring artifacts in X-ray micro-tomography,"*
26, 28396-28412 (2018). https://doi.org/10.1364/OE.26.028396

#### Authors:

Nghia Vo - *Diamond Light Source*

#### Contributors:

Daniel S. Hussey - *NIST: National Institute of Standards and Technology* 

How to use
==========
Clone or download the codes to your local machine, then insert the following two lines to your python codes:  

import sys  
sys.path.insert(0, '<path-to-sarepy-pck>')  

make sure that the python libs in the requirements.txt are installed before use.

Features
========
- Methods for cleaning different types of stripe artifacts: full stripes, partial stripes, unresponsive stripes, fluctuating stripes, large stripes, and blurry stripes.
- Various approaches based on the equalization-based methods, i.e equalizing the "response curves" of adjacent pixels, and their combinations.
- A robust stripe detection method.
- Implementations of former methods: a regularization-based method, a normalization-based method, a fft-based method, and a wavelet-fft-based method. 
- Matlab translation of algorithms 3,4,5,6.
- Implementations of a basic pipeline of tomography reconstruction: data loading, automated determination of center of rotation, ring artifact removal, tomographic reconstruction, and data saving.

Update Notes
------------
- 09/10/2019: Add methods for loading and saving data (from vounwarp), methods for calculating center of rotation, wrappers of reconstruction methods (from tomopy and astra). This allows users to try or test ring removal methods easily.
- 11/02/2020: Allow to use 2D kernel in the median filter of the sorting-based correction methods. This is done based on feedbacks from neutron imaging users. Note that this increases the computational cost. Add some sinograms for testing.
