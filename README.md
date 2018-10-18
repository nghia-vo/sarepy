# Sarepy
(S)tripe (A)rtifacts (RE)moval in (PY)thon
## Numerical techniques for eliminating ring artifacts in X-ray micro-tomography


**Sarepy** is the Python implementations of methods used for removing ring artifacts in tomography.
 These methods work in sinogram space where artifacts appear as straight lines or stripe artifacts.
 The codes are the original implementations of methods published in Optics Epxress,
 Nghia T. Vo, Robert C. Atwood, and Michael Drakopoulos, *"Superior techniques for eliminating ring artifacts in X-ray micro-tomography,"*
26, 28396-28412 (2018). https://doi.org/10.1364/OE.26.028396.

Features
========
- Methods for cleaning different types of stripe artifacts: full stripes, partial stripes, unresponsive stripes, and fluctuating stripes.
- Different implementations of the equalization-based methods, i.e equalizing the response curves of adjacent pixels.
- A stripe detection method based on the combination of processing techniques: normalizing, sorting, fitting, and thresholding.

