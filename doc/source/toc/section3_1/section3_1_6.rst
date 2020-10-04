Combinations of methods
=======================


.. |br| raw:: html

   <br />

In the previous sections and examples, only a few sinograms with certain types
of stripes were used for demonstration. In practice, one needs to process all
sinograms of a 3D sample where different types of stripes can appear at
different locations or all together in one place. It is impossible for a single
method to remove all of them. Each method has its own advantages and
disadvantages and only works well with the certain types of stripes. As a result,
we have to combine them to remove all types of stripes and reduce extra
artifacts the methods may cause. An efficient combination is that unresponsive
and fluctuating stripes are removed first, then large stripes are removed, and
smaller stripes of the full and partial type are removed last. Removal of large
stripes needs to be performed after removing unresponsive and fluctuating
stripes as explained in :ref:`section 3.1.5 <remove_dead_stripe>`. Methods of
removing smaller stripes need to be employed last because they may enlarge
large stripes. |br|
The following implementation of combining the methods is the simplest and most
efficient because it uses the least number of parameters. More importantly the
same set of parameters can be used for whole dataset.

.. autofunction:: sarepy.prep.stripe_removal_original.remove_all_stripe

How to use
  -- The *snr* parameter controls the sensitivity of the stripe detection
  method. Smaller is more sensitive. Recommended values: 1.1 -> 3.0. |br|
  -- The *la_size* parameter controls the strength of the median filter
  used for removing large stripes, unresponsive stripes, and fluctuating stripes.
  *la_size* can be determined in a straightforward way by the size and the brightness of
  detector defects. Larger is stronger but more computationally expensive. |br|
  -- The *sm_size* parameter controls the strength of the median filter used
  for removing small-to-medium partial and full stripes. |br| |br|

**Another combination of algorithms** which removes stripe artifacts based on interpolation
instead of using a smoothing filter is shown in section 4.2. |br|

  :ref:`4.2 Deriving an interpolation-based ring removal method <interpolation_method>`
