Prior methods for removing stripe artifacts
===========================================

.. |br| raw:: html

   <br />

.. _norm_method:

Normalization-based methods
---------------------------

Normalization-based methods rely on the assumption that stripe artifacts are
mainly the full stripe type, i.e the intensity offsets of these stripes are
constant across the angular direction. Different methods use different ways of
calculating these offsets. Sinogram correction is done by adding or multiplying
these offsets to the original sinogram. |br| |br|
**Method 1**

    .. autofunction:: sarepy.prep.stripe_removal_former.remove_stripe_based_normalization

    How it works
      Offsets are calculated by: averaging intensities along the angular
      direction; applying a smoothing filter to the result; subtracting these two.

    How to use
      -- The *sigma* parameter controls the strength of the cleaning capability.
      Larger is stronger. |br|
      -- The *num_chunk* parameter is used to divide a sinogram into multipe chunks of rows.
      This improvement helps to tackle partial stripes, however, it gives rise
      to void-center artifacts. |br|

    How to improve
      -- Can combine with the :ref:`sorting-based method <sorting>` to avoid
      void-center artifacts.

**Method 2**

    .. autofunction:: sarepy.prep.stripe_removal_former.remove_stripe_based_regularization

    How it works
      Firstly, intensities along the angular direction are averaged. Then, the
      1st discrete difference of the result is calculated. The offsets are
      determined by minimizing the differentiated result using a regularizer [2].

    How to use
      -- The *alpha* parameter controls the strength of the cleaning capability.
      Smaller is stronger. Recommended values: 0.01 -> 0.0001 |br|
      -- The *num_chunk* parameter is used to divide the sinogram into multipe chunks of rows.
      This improvement helps to tackle partial stripes, however, it gives rise
      to void-center artifacts. |br|

    How to improve
      -- Can combine with the :ref:`sorting-based method <sorting>` to avoid
      the void-center artifacts.

Fft-based methods
-----------------

Fft-based methods rely on the assumption that stripe artifacts are corresponding
to high-frequency components in the Fourier domain. As a result, they can
be removed by damping these components.  |br| |br|
**Method 1**

  .. autofunction:: sarepy.prep.stripe_removal_former.remove_stripe_based_fft

  How it works
    The Fourier transform is applied to a sinogram. A 1D low-pass window is
    multiplied with the row corresponding to the zero frequency in the vertical
    direction and its neighbors to remove stripes.

  How to use
    -- The *u* parameter controls the strength of the cleaning capability.
    Smaller is stronger. Recommended starting value: 30. |br|
    -- The *n* parameter defines the shape of the low-pass filter. It's an
    insensitive parameter. Recommended value : 8. |br|
    -- The *v* parameter allows to select how many rows around the
    zero-frequency row to be multiplied with the low-pass window. Larger *v*
    increases void-center artifacts. Recommended value: 1. |br|
    -- The *pad* parameter is needed to reduce side effects of the Fourier
    transform.

  How to improve
    -- Can combine with the :ref:`sorting-based method <sorting>` to avoid
    void-center artifacts.

**Method 2**

  .. autofunction:: sarepy.prep.stripe_removal_former.remove_stripe_based_wavelet_fft

  How it works
    It's very similar to the fft-based method. The improvement is that a sinogram is
    decomposed using the wavelet transform then a low-pass filter is applied
    to each of the decomposed image.

  How to use
    -- The *level* parameter controls the decomposition level. Higher "level"
    means stronger cleaning capability. It is because applying the low-pass
    filter at a deeper level (corresponding to a smaller-size image) resulting in
    stronger impact to the recombined image. |br|
    -- The *sigma* parameter also controls the strength of the cleaning
    capability. Larger is stronger, but also increases void-center artifacts. |br|
    -- The *order* parameter is insensitive. Recommended value: 8. |br|
    -- The *pad* parameter is needed to reduce side effects of the Fourier
    transform.

  How to improve
      -- Can combine with the :ref:`sorting-based method <sorting>` to avoid
      void-center artifacts.

**REFERENCES**

1. M. Rivers, "Tutorial introduction to X-ray computed microtomography data processing,"
   http://www.mcs.anl.gov/research/projects/X-ray-cmt/rivers/tutorial.html.
2. S. Titarenko, P. J. Withers, and A. Yagola, "An analytical formula for ring
   artefact suppression in X-ray tomography," Appl. Math. Lett. 23(12), 1489-1495 (2010).
3. C. Raven, "Numerical removal of ring artifacts in microtomography,"
   Rev. Sci. Instrum. 69(8), 2978-2980 (1998).
4. B. Munch, P. Trtik, F. Marone, and M. Stampanoni, "Stripe and ring artifact
   removal with combined wavelet-Fourier filtering," Opt. Express 17(10), 8567-8591 (2009).
