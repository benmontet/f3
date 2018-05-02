.. image:: https://zenodo.org/badge/90557386.svg
   :target: https://zenodo.org/badge/latestdoi/90557386
.. image:: http://img.shields.io/badge/arXiv-1705.07928-orange.svg?style=flat
        :target: http://arxiv.org/abs/1705.07928

f3
===

**Full Frame Fotometry package for interacting with and measuring photometry
for targets observed in the Kepler Full Frame Images**

You can read the documentation at: `f3.readthedocs.io <http://f3.readthedocs.io>`_.

Installation
------------

Start by installing the dependencies:

1. numpy
2. matplotlib
3. scipy
4. astropy
5. kplr (https://github.com/dfm/kplr) (The version pip will get you is insufficent, grab the one on github)
6. mahotas (http://mahotas.readthedocs.io/en/latest/#) (pip will install this for you along with f3)

Running “pip install f3” or cloning this repository and running setup.py should then
enable you to run the code.


Get the data
------------

To download the Kepler full frame images, run::

    wget -r -nH -nd -np -R index.html -e robots=off https://archive.stsci.edu/pub/kepler/ffi/

Note: this is 44gb worth of data!

f3 will look for this in the “ffidata” subdirectory inside your working directory, but
this can be changed in an argument when you initialize your object.



Use the code
----------------

There is a jupyter notebook in this directory which provides an example use case for this package.


Cite the code
----------------

Please cite `Montet, Tovar, and Foreman-Mackey, 2017, ApJ, 851, 116’ if you find this code
useful in your research. You can find a preprint of the paper at https://arxiv.org/abs/1705.07928.
The bibtex entry for this paper, thanks to NASA ADS, is

.. code-block:: tex

    @ARTICLE{f3,
   author = {{Montet}, B.~T. and {Tovar}, G. and {Foreman-Mackey}, D.},
    title = "{Long-term Photometric Variability in Kepler Full-frame Images: Magnetic Cycles of Sun{\ndash}like Stars}",
  journal = {\apj},
archivePrefix = "arXiv",
   eprint = {1705.07928},
 primaryClass = "astro-ph.SR",
 keywords = {methods: data analysis, stars: activity, stars: solar-type, techniques: photometric },
     year = 2017,
    month = dec,
   volume = 851,
      eid = {116},
    pages = {116},
      doi = {10.3847/1538-4357/aa9e00},
   adsurl = {http://adsabs.harvard.edu/abs/2017ApJ...851..116M},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}


