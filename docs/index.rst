f3 Photometry for Kepler’s Full Frame Images
============================================

f3 is a package for doing photometry with Kepler's Full Frame Images, a set of calibration data 
obtained approximately monthly during the primary Kepler mission. While the main mission targeted
about 200,000 stars at 30-minute cadence, there are a total of 4.5 million objects in the Kepler
field that fall in the telescope's field of view but are only accessible through the FFIs.

Our original work with the FFIs can be found at arXiv:1608.01316, while this code is presented
in arXiv:1705.07928.

.. toctree::
   :maxdepth: 7

   install
   photometry
   tutorial



License & Attribution
---------------------

Copyright 2017, Benjamin Montet and contributors.

The source code is made available under the terms of the MIT license.

If you make use of this code, please cite the following paper:

.. code-block:: tex


    @ARTICLE{f3,
       author = {{Montet}, B.~T. and {Tovar}, G. and {Foreman-Mackey}, D.},
        title = "{Long Term Photometric Variability in Kepler Full Frame Images: Magnetic Cycles of Sun-Like Stars}",
      journal = {ArXiv e-prints},
    archivePrefix = "arXiv",
       eprint = {1705.07928},
     primaryClass = "astro-ph.SR",
     keywords = {Astrophysics - Solar and Stellar Astrophysics},
         year = 2017,
        month = may,
       adsurl = {http://adsabs.harvard.edu/abs/2017arXiv170507928M},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

