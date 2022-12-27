Tomographic reconstruction pipeline for testing ring removal methods
====================================================================

.. |br| raw:: html

   <br />

To assist users in choosing ring removal methods and turning parameters,
basic data processing modules are added to Sarepy.

**Input-output module**

    This module allows to load/save data from/to an image file (tif, jpg, png...)
    or an hdf file. For example:

    .. code-block:: py

        import sarepy.losa.loadersaver as losa

        # Load an image
        sinogram = losa.load_image("C:/data/sinogram.tif")
        # Save to an image
        losa.save_image("C:/data/image.tif", sinogram)

**Center-of-rotation determination**

    It is used to `calculate the center-of-rotation <https://doi.org/10.1364/OE.22.019078>`_
    from a 180-degree sinogram. If you use a 360-degree sinogram, remember to use
    only a haft of it. For example:

    .. code-block:: py

        import sarepy.prep.autocentering as cen

        (height, width) = sinogram.shape
        # Find the center of rotation by searching around the middle of the image
        # with the radius of 50 pixels
        center = cen.find_center_vo(sinogram, width//2 - 50, width//2 + 50)

**Preprocessing methods for removing stripe artifacts**

    This is the main module of Sarepy with loads of `stripe removal methods <https://doi.org/10.1364/OE.26.028396>`_. Feel free to break them.
    Making sure that there're no negative values in the input.

    .. code-block:: py

        import numpy as np
        import sarepy.prep.stripe_removal_original as srm1
        import sarepy.prep.stripe_removal_improved as srm2
        import sarepy.prep.stripe_removal_former as srm3

        nmean = np.mean(sinogram)
        sinogram[sinogram <= 0.0] = nmean
        # Using original methods
        sinogram1 = srm1.remove_all_stripe(sinogram, 3.0, 31, 11)
        sinogram2 = srm1.remove_unresponsive_and_fluctuating_stripe(sinogram, 3.0, 31)
        sinogram3 = srm1.remove_large_stripe(sinogram, 3.0, 31)
        sinogram4 = srm1.remove_stripe_based_sorting(sinogram, 11)
        sinogram5 = srm1.remove_stripe_based_fitting(sinogram, 2, 5, 60)
        sinogram6 = srm1.remove_stripe_based_filtering(sinogram, 3, 11)

        # Using improved methods
        sinogram7 = srm2.remove_stripe_based_2d_filtering_sorting(sinogram, 11, 11)
        sinogram8 = srm2.remove_stripe_based_filtering_sorting(sinogram, 3, 11)
        sinogram9 = srm2.remove_stripe_based_sorting_fitting(sinogram, 1, 10, 60)

        # Using former methods
        sinogram10 = srm3.remove_stripe_based_fft(sinogram, 20, 8, 1, 100)
        sinogram11 = srm3.remove_stripe_based_normalization(sinogram, 15, 1)
        sinogram12 = srm3.remove_stripe_based_regularization(sinogram, 0.005, 1)
        sinogram13 = srm3.remove_stripe_based_wavelet_fft(sinogram, 5, 1.0, 8, 100)

**Tomographic reconstruction**

    This module is a wrapper of reconstruction methods distributed with the `Tomopy <https://tomopy.readthedocs.io/>`_
    and `Astra toolbox <https://www.astra-toolbox.com/>`_ package.
    If you use a 360_degree sinogram, remember to pass the angles.

    .. code-block:: py

        import sarepy.reco.reconstruction as rec

        ratio = (min(width - center, center))/(0.5 * width)
        # If GPU is available.
        rec_image1 = rec.recon_astra(sinogram, center, angles=None, ratio=ratio)
        # If only CPU is available.
        rec_image2 = rec.recon_gridrec(sinogram, center, angles=None, ratio=ratio)
        # Save results to images
        losa.save_image("C:/data/rec_image1.tif", rec_image1)
        losa.save_image("C:/data/rec_image2.tif", rec_image2)

**Postprocessing methods for removing ring artifacts**

    Users can  have a go with postprocessing methods if need to.

    .. code-block:: py

        import sarepy.post.ring_removal_post as rrp

        rec_image1 = rrp.ring_removal_based_wavelet_fft(rec_image1)
        rec_image2 = rrp.ring_removal_based_fft(rec_image2)
        # Save results to images
        losa.save_image("C:/data/rec_image1.tif", rec_image1)
        losa.save_image("C:/data/rec_image2.tif", rec_image2)

If users would like to develop new approaches for removing ring artifacts, it is highly
recommend that these methods should be tested on challenging `sinograms <https://github.com/nghia-vo/sarepy/tree/master/data>`__
or a full size of tomographic data having all types of ring artifacts `publicly available here <https://doi.org/10.5281/zenodo.1443568>`__.
