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

Features
========
- Methods for cleaning different types of stripe artifacts: full stripes, partial stripes, unresponsive stripes, fluctuating stripes, large stripes, and blurry stripes.
- Various approaches based on the equalization-based methods, i.e equalizing the "response curves" of adjacent pixels, and their combinations.
- A robust stripe detection method.
- Implementations of former methods: a regularization-based method, a normalization-based method, a fft-based method, and a wavelet-fft-based method. 
- Matlab translation of algorithms 3,4,5,6.
- Implementations of a basic pipeline of tomography reconstruction: data loading, automated determination of center of rotation, ring artifact removal, tomographic reconstruction, and data saving.

