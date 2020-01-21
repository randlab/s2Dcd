Code purpose
================================================================================

.. .. figure:: ../logo/tutto.gif
   :width: 200
   :align: left   
   :alt: an illustration of the s2Dcd approach


``s2Dcd`` :

    ``s``   sequential

    ``2D``  bi-dimensional (multiple-point statistics simulation with)

    ``cd`` conditioning data

The code `s2Dcd` allows to apply the sequential 2D (multiple-point
simulation) with conditioning data approach described in the paper by
[Comunian2012]_. For a list of publications where the `s2Dcd` was
used, please see the section `publications <publications.html>`_.

.. note:: At the moment the ``s2Dcd`` is it only a **wrapper**
   library: it requires an external MPS simulation engine to work.

.. warning:: Since 2017 only the direct sampling (DS) version of the
             ``s2Dcd`` is supported. Therefore, in the following
             documentation all references to the ``Impala`` MPS
             simulation engine should be considered outdated.
          

Not only ``s2Dcd``...
-------------------------------------------------------------------------------
The package can be used to do other things in addition to the
application of the `s2Dcd` approach. For example, it contains some
modules which are used as interface to the multiple-point simulation
(MPS) codes `DeeSse` and `Impala`. Therefore, if for example you need
to run a number of simulation tests with different parameters, you can
use `Python` and the interface to the two codes to create customized
benchmarks with the flexibility provided by `Python`. See the
documentation provided in the appendices and the examples for more
details.

Also, some simple functions that allows to read and write from `numpy`
to `GSLIB` and `VTK` are provided.

If you have any questions, you want to contribute, suggest new
features, you've found a bug...  you can write me an e-mail:
`alessandro DOT comunian AT gmail DOT com`


MPS simulation engines
----------------------------------

As already mentioned, `s2Dcd` is only a wrapper library and requires
an MPS simulation engine to work.  This version contains and interface
to the DeeSse MPS simulation engine (see `this link
<http://www.randlab.org/research/deesse/>`_ for mode info.)  Another
simulation engine (`Impala`) was supported until 2017, but the current
development of the interface to this engine is discontinued.

Nevertheless, users are encouraged to develop interfaces to their
favorite MPS simulation engine (like for example `SGeMS
<http://sgems.sourceforge.net/>`_ or `MPSLib
<https://github.com/ergosimulation/mpslib>`_).


What's new in this version
----------------------------
* Removed the dependencies from the `VTK` libraries.
* Removed the dependencies from some additional modules derived by
  `Fortran90` subroutines.
* Solved some bugs.
* Added some new examples.  


If you use ``s2Dcd``
----------------------------------
Please:

* Let us know!
* Cite the paper [Comunian2012]_

.. [Comunian2012] Comunian, A.; Renard, P. and Straubhaar, J.
                  *3D multiple-point statistics simulation using 2D
		  training images* Computers & Geosciences, 2012, **40**, 49-65.
		  `doi:10.1016/j.cageo.2011.07.009 <http://dx.doi.org/10.1016/j.cageo.2011.07.009>`_

Acknowledgments
----------------------------

Many thanks to Philippe Renard and Julien Straubhaar for their support
and suggestions, to Andrea Borghi for the fruitful discussions and for
finding some bugs, and to my brother Thomas for the design of the
s2Dcd logo.


