Installation
============

f3 is written in Python. It can be installed using `pip <https://pip.pypa.io>`_ 
with the command:

.. code-block:: bash
 
     pip install f3

.. _source:

From Source
-----------

The source code can be downloaded `from GitHub 
<https://github.com/benmontet/f3>`_ by running


.. code-block:: bash

    git clone https://github.com/benmontet/f3.git

and then can be built by running

.. code-block:: bash

    python setup.py install

in the root directory of the source tree.

.. _python-deps:

Dependencies
++++++++++++

In addition to a Python installation, there are a few other dependencies that
are required to run f3:

1. `NumPy <http://www.numpy.org/>`_
2. `Matplotlib <https://matplotlib.org/>`_
3. `SciPy <https://www.scipy.org/>`_
4. `pyfits <http://www.stsci.edu/institute/software_hardware/pyfits>`_
5. `kplr <https://github.com/dfm/kplr>`_
6. `mahotas <http://mahotas.readthedocs.io/en/latest/>`_

mahotas will be automagically installed along with f3 by running setup.py (or through pip). As of June 5, 2017, the kplr version available on pip is not sufficient for running 
f3, you need a version dated June 11, 2016 or newer: the version on github linked above
will be sufficient. 

You will also need the Full Frame Images themselves, which are 44gb and not a part of
the installation package. To get them, run

.. code-block:: bash
 
     wget -r -nH -nd -np -R index.html -e robots=off https://archive.stsci.edu/pub/kepler/ffi/

f3 will by default look for these in a subdirectory from your working directory called “ffidata” but they can be placed anywhere.
