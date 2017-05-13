f3
===

**Full Frame Fotometry package for interacting with and measuring photometry
for targets observed in the Kepler Full Frame Images**

Installation
------------

Start by installing the dependencies:

1. numpy
2. matplotlib
3. scipy
3. pyfits (I know)
4. kplr (https://github.com/dfm/kplr)
5. mahotas 

Running “pip install f3” or cloning this repository and running setup.py should then
enable you to run the code.


Get the data
------------

To download the Kepler full frame images, run::

    wget -r -nH -nd -np -R index.html -e robots=off https://archive.stsci.edu/pub/kepler/ffi/

f3 will look for this in the “ffidata” subdirectory inside your working directory, but
this can be changed in an argument when you initialize your object.



Use the code
----------------

There is a jupyter notebook in this directory which provides an example use case for this package.

Cite the code
_________________

Please cite `Montet, Tovar, and Foreman-Mackey (to be submitted) if you find this code
useful in your research. If it is after May 18, this code is available on the arXiv and we’ll put a link here when that is true.