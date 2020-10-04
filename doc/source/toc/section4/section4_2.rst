.. _interpolation_method:

Deriving an interpolation-based ring removal method
===================================================

.. |br| raw:: html

   <br />

The methods presented in section 3 mainly use smoothing filters to remove
stripe artifacts. This certainly affects neighboring pixels resulting in
side-effect artifacts. Here, the author derives an method where stripe artifacts
are removed by interpolation with neighbors. The method is a combination of
algorithm 4, 5, and 6 in Ref. [1].

**Python source code**

  .. autofunction:: sarepy.prep.stripe_removal_improved.remove_stripe_based_interpolation

**How it works**
  It works like the method of removing :ref:`large stripes <remove_large_stripe>`.
  However, in step 3 the selective correction is replaced by interpolation, the same as
  used in the method of :ref:`removing unresponsive and fluctuating stripes <remove_dead_stripe>`.

**Demonstration**
  - Remove small-to-medium stripes:

  .. code-block:: py

      import numpy as np
      import sarepy.losa.loadersaver as losa
      import sarepy.prep.autocentering as cen
      import sarepy.reco.reconstruction as rec
      from sarepy.prep.stripe_removal_improved import remove_stripe_based_interpolation

      sinogram = losa.load_image("../data/sinogram_normal.tif")
      (nrow, ncol) = sinogram.shape
      center = cen.find_center_vo(sinogram, ncol // 2 - 30, ncol // 2 + 30, step_cor=0.5)
      ratio = (min(center, abs(ncol - center))) / (ncol * 0.5)
      rec_image = rec.recon_gridrec(-np.log(sinogram), center, ratio=ratio)
      losa.save_image("../data/output/rec1_before.tif", rec_image)
      sinogram = remove_stripe_based_interpolation(sinogram, 3.0, 31)
      rec_image = rec.recon_gridrec(-np.log(sinogram), center, ratio=ratio)
      losa.save_image("../data/output/rec1_after.tif", rec_image)

  .. figure:: section4_2_figs/fig1.jpg
    :figwidth: 100 %

    Figure 1. Reconstructed images before (a) and after (b) removing stripe artifacts.

  - Remove large stripes:

  .. code-block:: py

      sinogram = losa.load_image("../data/sinogram_large_stripe.tif")
      (nrow, ncol) = sinogram.shape
      center = cen.find_center_vo(sinogram, ncol // 2 - 30, ncol // 2 + 30, step_cor=0.5)
      ratio = (min(center, abs(ncol - center))) / (ncol * 0.5)
      rec_image = rec.recon_gridrec(-np.log(sinogram), center, ratio=ratio)
      losa.save_image("../data/output/rec2_before.tif", rec_image)
      sinogram = remove_stripe_based_interpolation(sinogram, 3.0, 51)
      rec_image = rec.recon_gridrec(-np.log(sinogram), center, ratio=ratio)
      losa.save_image("../data/output/rec2_after.tif", rec_image)

  .. figure:: section4_2_figs/fig2.jpg
    :figwidth: 100 %

    Figure 2. Reconstructed images before (a) and after (b) removing stripe artifacts.

    - Remove dead stripes:

    .. code-block:: py

      sinogram = losa.load_image("../data/sinogram_dead_stripe.tif")
      (nrow, ncol) = sinogram.shape
      center = cen.find_center_vo(sinogram, ncol // 2 - 30, ncol // 2 + 30, step_cor=0.5)
      ratio = (min(center, abs(ncol - center))) / (ncol * 0.5)
      rec_image = rec.recon_gridrec(-np.log(sinogram), center, ratio=ratio)
      losa.save_image("../data/output/rec3_before.tif", rec_image)
      sinogram = remove_stripe_based_interpolation(sinogram, 3.0, 51)
      rec_image = rec.recon_gridrec(-np.log(sinogram), center, ratio=ratio)
      losa.save_image("../data/output/rec3_after.tif", rec_image)

    .. figure:: section4_2_figs/fig3.jpg
      :figwidth: 100 %

      Figure 3. Reconstructed images before (a) and after (b) removing stripe artifacts.
