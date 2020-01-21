Notes
================================

Auxiliary variables
-------------------------------------

The information provided in this section are useful to set up the
template file and the *Python* script useful to run the ``s2Dcd``
using auxiliary variables.

Implementation in DeeSse
*********************************

.. note:: This option might work properly, but it is in a early stage
          development.

The use of auxiliary variables with the `DeeSse`
MPS engine could be implemented in a number of ways. Here the
following strategy is selected:

1) The 2D slices (primary and auxiliary variable) are provided in the
   *same* `GSLIB` file using the variable name provided in the
   ``template.in``
   file. In other words, only one TI is defined in the section ``TRAINING
   IMAGE``.
2) The intermediate output files are provided as a unique file
   too. Therefore, in the section ``OUTPUT SETTINGS FOR SIMULATION``
   in the ``template.in`` file the keyword
   ``OUTPUT_SIM_ALL_IN_ONE_FILE`` should be used.


    
    

Implementation in Impala
*******************************

.. warning:: This MPS engine is not supported any more
	     since 2017. This section is left just as a reference for
	     further developments, potentially oriented towards other
	     MPS simulation engines.

There are two main modes that can be implemented to use auxiliary variables:

Full 3D mode
    In this mode a full 3D auxliary variable map must be
    provided. Then the `s2Dcd` automatically slices it along the
    current simulated section to extract a 2D auxiliary variable
    map. This set up is useful when it is relatively easy to find a
    full 3D map for the auxiliary variable.

2D mode
    This should be the default mode, when you use the 2D maps of
    auxiliary variables attached to the 2D TIs provided for the
    simulation.

The `Full 3D mode` can be enabled by setting the
``geostats.s2Dcd.AUX_VAR_FULL3D`` parameter to ``True``. Otherwise,
the default value for the `2D mode` is used.

For the ``Impala`` MPS simulation engine the implementation of usage
of the auxiliary variable is quite univoque and follows a quite
standard workflow.  Simply, the user have to provide (in addition to
the 2D TIs of the main variable):

* a 2D auxiliary variable map (`VTK` files) with the same size of
  the provided 2D TIs, one for each TI.
* a 3D auxiliary variable map with the same size of the 3D simulation
  grid.

.. note::
    1) The format of the variable in the *VTK* files contaning the
       auxiliary variables should be ``float``.
    2) During the simulation process many *VTK* files containing
       sections of the 3D auxiliary variable map will be created. This
       is somehow redundant and in the future should be removed. For
       the moment you have to manually remove manually these files.
    3) Very often the information contained in the 3D auxiliary
       variable sections is not very different from the 2D auxiliary
       variable corresponding to the 2D TI. Some computation time
       could be saved in these cases, but for the moment a complete 3D
       auxliary variable map is required to provide more flexibility.


Licence Issues
---------------------------

The recent versions of the DS codes are running with a licence
manager. Therefore, it is possible that running the code will give an
error code. At the moment the quick and dirty solution is to wait for
some time and retry to run again the simulation of the same
section. You can tune this in the module ``deesse.py``, changing the
values of the variables ``NB_LIC_WAIT`` (number of attempts to contact
the lincese server) and ``LIC_WAIT_TIME`` (pause from one attempt to the
other).

.. warning:: Note that the computing time can be heavily affected by
             this problem!




