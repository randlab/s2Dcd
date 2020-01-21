Installation
===============================================================================

All the modules included in the "package" are written in Python
(version 3.X), which is platform independent. Therefore, the "package"
should work on Linux, Mac OS X and Windows. Some tests were performed
on all the three platforms, but at the moment the MPS engines
implemented in the `s2Dcd` are available only for Linux and Windows.

Requirements
--------------------------------------------------------------------------------


* `Python <http://python.org>`_ (tested version 3.6.9 on Linux). Older versions should also work.
* If you are on Linux, you will also need the python development packages.
* The Python numerical library `numpy <http://numpy.scipy.org/>`_ and
  `pandas <http://pandas.pydata.org/>`_ and the Python package installer ``pip``.
* The Multiple-point statistics simulation engines `DeeSse` and/or `Impala`.

.. sidebar:: Suggested software (optional)

   Some software that you might find useful:

   * `ParaView <http://www.paraview.org/>`_, useful to visualize the
     `VTK` output of `Impala`.
   * `SGeMS <http://sgems.sourceforge.net/>`_, useful to visualize the
     `GSLIB` output of `DeeSse`.



Python, numpy and pandas installation
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The "easy" way
*********************

If you do not want to install all the Python packages separately,
there are some bundled distributions that allow to install all the
required libraries in "one click". One of this distribution is
`Anaconda <https://www.continuum.io/why-anaconda>`_, which is
available for MS Windows, Linux and Mac OS X (note however that this
will require 200/300 Mbytes on your hard drive). Alternatively, you
can use another "bundled distribution" or install each python package
separately as explained in the following sections.

.. warning:: When you download Anaconda, select the Python *3.X*
             version and **NOT** the one for Python *2.X*.

.. note:: Personally, when working on Linux, I prefer to install Python packages
	  using the OS package manager.
	     
MS Windows
************

If you are on MS Windows, I strongly suggest you to use something like
the aforementioned bundled distribution `Anaconda`.


Linux
*************

Use you package manager, or ``sudo apt-get install`` to install
``python3``, ``python3-dev``, ``python3-numpy`` and
``python3-pandas``.  If `pandas` is not available for your
distribution, please check on the last installation instructions `here
<http://pandas.pydata.org/>`_.


Mac OS X
**************
Python and numpy are available for Mac OS X too. However, at the
moment the binaries of the MPS simulation engines with an interface to
the `s2Dcd` are not available for Mac OS X. You can still use some
utility scripts to work with `GSLIB` files.


The multiple-point simulation engine `DeeSse` 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Actually, the python module `s2Dcd` is only a wrapper for a MPS
simulation engine. For the moment, the implementation allows to use the
MPS simulation engine `DeeSse`.

The executable of your MPS simulation engine (in the case of the
`DeeSse` they are called ``deesse`` and ``deesseOMP``; add ``.exe`` if
you are working on MS Windows), should be located in your ``bin``
directory, in the working directory or in a directory listed in the
``PATH`` environment variable.

.. sidebar:: Adapt the names of the MPS binaries

   If the provided binaries have a different name...
   
   The name of the binary of the simulation code can be changed
   using the variable `s2Dcd.mps_exes`, which is a Python dictionary.

To add a directory (i.e. ``C:\Users\alex\my_dir`` for MS Windows, or ``/home/alex/my_dir``)
to your path variable do:

MS Windows
******************
1. `Start > Control Panel > System and Security > Advances System Settings`
2. Click on the button `Environment variables` at the bottom right.
3. Under `System Variables`, select the line containing the variable `Path`.
4. Edit the variable: add to the end of the variable value  ``;C:\Users\alex\my_dir\``
   (note the comma **;**)

Linux
****************
Add the following line to your `.bash_aliases` file (eventually, create it and make
sure that your `.bashrc` loads it)::
    
    export PATH=${PATH}:/home/alex/my_dir

.. note::
    Before testing the `s2Dcd` module, check if the MPS simulation
    engines work correctly.


Installation of the module `s2Dcd` and other python tools 
--------------------------------------------------------------------------------

Download the sources from Github. From the command line (use the Anaconda prompt if you installed it), move into the directory
``s2Dcd`` (if needed, unpack it)  and type the command::

  pip install -e .

To check if the installation was successful, you can run the python shell and try import the 
main module::

>>> import s2Dcd.s2Dcd




Some details about the main modules:
-------------------------------------

`s2Dcd.py`
    the main module containing for the functions for the `s2Dcd` simulations.
`deesse.py`
    a simple interface to the parameters required by `DeeSse`.
`gslibnumpy.py`
    to convert from numpy and GSLIB and vice versa.
`utili.py`
    some simple utilities...

More details about these in section :ref:`appendices`.

