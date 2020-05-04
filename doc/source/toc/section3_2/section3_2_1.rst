Fft-based methods for removing ring artifacts
=============================================

.. |br| raw:: html

   <br />

Fft-based methods rely on the assumption that ring artifacts are corresponding
to high-frequency components in the Fourier domain. As a result, they can
be removed by damping these components. The following methods are similar to
the preprocessing fft-based methods. However, they work on reconstructed images
instead of sinograms where the polar transform is used to convert ring artifacts to
stripe artifacts. |br| |br|
**Method 1**

    .. autofunction:: sarepy.post.ring_removal_post.ring_removal_based_fft

    How it works
     Input image is transformed into polar coordinates. The Fourier
     transform is applied to the result. A 1D low-pass window is
     multiplied with the row corresponding to the zero frequency in the vertical
     direction and its neighbors to remove the stripes. The resulting image is
     transformed into Cartesian coordinates.

    How to use
     -- The *u* parameter controls the strength of the cleaning capability.
     Smaller is stronger. Recommended starting value: 30. |br|
     -- The *n* parameter defines the shape of the low-pass filter. It's an
     insensitive parameter. Recommended value : 8. |br|
     -- The *v* parameter allows to select how many rows around the
     zero-frequency row to be multiplied with the low-pass window. Larger *v*
     increases void-center artifacts. Recommended value: 1. |br|
     -- The *pad* parameter is needed to reduce the side effects of the Fourier
     transform.

**Method 2**

    .. autofunction:: sarepy.post.ring_removal_post.ring_removal_based_wavelet_fft

    How it works
     It's very similar to the fft-based method. The improvement is that the polar
     transformed image is decomposed using the wavelet transform then the low-pass filter is applied
     to each of the decomposed image.

    How to use
     -- The *level* parameter controls the decomposition level. Higher "level"
     means stronger cleaning capability. It is because applying the low-pass
     filter at a deeper level (corresponding to a smaller-size image) results in
     stronger impact to the recombined image. |br|
     -- The *sigma* parameter also controls the strength of the cleaning
     capability. Larger is stronger, but also increases void-center artifacts. |br|
     -- The *order* parameter is insensitive. Recommended value: 8. |br|
     -- The *pad* parameter is needed to reduce the side effects of the Fourier
     transform.
