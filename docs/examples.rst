Examples
======================================================================
 
Some usage examples of the ``s2Dcd`` module, including a `commented`
example.

.. warning:: Some of the examples contained in the provided directory
             are very work in progress. If the directory related to an
             example contains a file ``README.txt`` where it is
             clearly stated that the example is **work in progress**,
             do not trust too much the corresponding example.

How to run the examples
------------------------------

MS Windows
*****************

.. note:: Here we suppose that you are using the default GUI interface
          ('IDLE') provided with Python. However, if you installed
          Python using `Anaconda`, there is a nicer user friendly GUI
          named `spyderlib` that you can use instead.

If you have installed `Python`, then you should have also the Python
`IDLE` (a GUI interface). From the interface you can open one of the
files in the ``examples`` directory, and from the menu of the `IDLE`
select `run`. Another window should open and show you the advancement
of the simulation.

Alternatively, you can use the `Command Prompt`, like this::

  C:\Users\toto\examples\01_Strebelle> python.exe  s2Dcd_run-ex01.py



.. note:: If you are using `spyderlib` to run
          your scripts (and you don't have any Python installation
          other than the one that comes with the `Anaconda`), it
          can be useful, once you set the ``PYTHONPATH`` environment
          variable from the `spyderlib` menu (inside the `Python Path
          Manager` menu), to push the button `synchronize...` (and
          accept with `Yes`). This should allow you to use `s2Dcd` also
          with an external python shell. Then, to run the `s2Dcd`, go
          to `Run > Configuration per file...` and select `Execute in an
          external system terminal`, and tick `Interact with the python
          console after execution`. Then you should able to run the
          `s2Dcd` with these option without seeing plenty of
          command-prompt windows popping up...
	  Note that you will probably need to restart the Python kernel.
  

Linux
**************
In the directory of the each example, you can simply::

    username@machine$ ./s2Dcd_run-ex01.py

You might need to provide the execution rights to your user
(for example ``chmod +x s2Dcd_run-ex01.py`` ).

If you want, you can redirect the output to some file and submit the process
in background::

   username@machine$ ./s2Dcd_run-ex01.py > somefile.out &


Commented example
------------------------------------------------------------------------
This is a commented version of the file ``s2Dcd_run-ex02.py`` that you
can find in the ``..\examples\02_Strebelle-conditional``
directory.  It runs a simulation using two training images and some
conditioning data, both in the form of data points and in the form of
a conditioning slice.  The example presents the use of the `DeeSse`
MPS simulation engine, but the principles to run a simulation with a
different MPS simulation engine (i.e. `Impala`) are similar.

.. note:: 
    In python, all the text enclosed by ``'''``, ``"""`` and the lines
    starting with ``#`` are comments.

In the first row we define the python interpreter... this is a
standard command that can be included in each python script::

#!/usr/bin/env python3

Import some *standard* python modules::

    import os
    import time
    import sys
    import numpy as np
    import random
    import copy

Import the modules which are part of the `s2Dcd` software, that is::

    import s2Dcd.s2Dcd as s2Dcd

a module that contains all the interface for the `DeeSse` MPS
simulation engine::

    import s2Dcd.deesse as mpds_interface

.. note::
   If you want to use a different MPS engine `Impala`, you need to create your own module for this.

a simple module containing some utilities::
   
    import s2Dcd.utili as utili

a module to read and write output in the GSLIB format::

    import s2Dcd.gslibnumpy as gslibnumpy

If you can use a multi-thread version of `DeeSse`, you can define here
the number of threads::

    s2Dcd.nb_threads = 4

if you don't specify a value for this variable, the default value of 1
is used.

Then, the following line is used to print out some information about
the run and record the start time::

    time_start = utili.print_start() 

A random seed can be set with the following command, in order to allow
to get repeatable results::

    seed = numpy.random.RandomState(456833)

Then you have to load the default definition for the parameters for
the `DeeSse` simulation. In our case, we set all the parameters in a
template file with in mind the resulting 3D simulation. Therefore, the
grid size and the definition of the search template will be provided
as we would do a 3D simulation.  Here all the parameters are contained
in a standard ``.in`` file for `DeeSse`.  We therefore create a the
template of parameters (``par_template``) from the file
``template.in``::

    template_in = "template.in"
    par_template = mpds_interface.Param(file_name=template_in)

.. note:: You could also call ``Param`` without the argument
    ``file_name``, and you will get the default parameters defined in
    the ``mpds_interface``.

Now we define the parameters which are used for the simulation along
the `x`, `y` and `z` axes. We always read the default parameters from
the file defined by the variable ``template_in``, but we could
potentially select a different set of parameters for each direction
(and each step of the sequence). Here we keep the parameters defined
in the default ``template.in`` file. However, we customize them
defining a different training image for each direction.  The value
``None`` must be specified if we don't have a TI for the plane normal
to that simulation direction::

    par_Xnorm = mpds_interface.Param(file_name=template_in)
    par_Xnorm.tis[0].file_name = None
    par_Ynorm = mpds_interface.Param(file_name=template_in)
    par_Ynorm.tis[0].file_name = "ti_250x1x250.gslib"
    par_Znorm = mpds_interface.Param(file_name=template_in)
    par_Znorm.tis[0].file_name = "ti_250x250x1.gslib"

Then you have to specify the maximum number of simulation step.
You can usually put this value to a number bigger that the expected
number of slices required to complete your simulation domain. The
simulation will stop when the simulation domain will be filled.

However, in some cases when you would like to run only the fist step
of the simulation, to check is the simulation and the conditioning are
OK, you can set this value to a smaller integer, for example ``6``.
In this case we select a big value::

    step_max = 3000

Now we define the simulation grid, which is extracted from the size of
the grid contained in the file ``template.in`` and therefore in the
variable ``par_template``::

    simODS = par_template.grid

the grid is printed to the standard output (for double check),
together with other information about the simulation::

    s2Dcd.print_sim_info(simODS, par_Xnorm, par_Ynorm, par_Znorm)

Then, it is required to define a numpy array which will contain the
results of each step of the simulation. The default values of the
array are initialized using the variable ``s2Dcd.no_data``. The
simulation will go on until all the `no_data` values will be simulated
(or until the maximum simulation step ``step_max``)::

    hard_data = s2Dcd.no_data * numpy.ones(
        (simODS.nx, simODS.ny, simODS.nz), 'float')

The default value for the variable ``s2Dcd.no_data`` is -1.

.. note::
    You can redefine the value of the variable ``s2Dcd.no_data`` to
    suit different conversions for the no data values.
    In this example we consider the default value -1, therefore a
    definition is not required. 
    If you want to change the no data value you can do something like::
        
        s2Dcd.no_data = -999999    

The conditioning data, if available, can be loaded in the GSLIB point
data format using the following command::

    s2Dcd.add_gslib_pointdata(["data_points.gslib"], hard_data, simODS)

.. note:: If you have many files containing conditioning data in the
    GSLIB point data format, then you can add them at the same time
    like this::

            s2Dcd.add_gslib_pointdata(["file1.gslib",
	        "file2.gslib","file3.gslib"], hard_data, simODS)

With this command the data points contained in the file
``data_points.gslib`` will be added to the numpy array that contains
all the results. They will also be extracted at each simulation step
if they are located in the simulated section.  You can also add
conditioning data in a structured grid format, for example if you have
access to an outcrop or some data saved as a GSLIB "image".  With the
following command, we load the file ``ti_250x1x250.gslib`` and we
associate its values to the section with `index` 4 along the `y`
coordinate (the exact location depends in the definition of your
grid...)::

    hd_section = gslibnumpy.gslib2numpy("ti_250x1x250.gslib")
    hard_data[:,4,:] = hd_section[0:simODS.nx,0,0:simODS.nz]

Now we can create the simulation sequence. Basically, we alternatively
simulate along the directions defined by the parameters defined by
``par_Xnorm``, ``par_Ynorm``, and ``par_Znorm`` only when a training
image file is present. In this case we didn't provided a value for the
TI file for the sections perpendicular to the plane `yz` (that is the
parameters contained in ``par_Xnorm``). Moreover, for a given
direction we always use the same parameters along the duration of the
simulation sequence. Therefore, to initialize the simulation sequence
we only need::

    seq = s2Dcd.create_seq(simODS, par_Xnorm, par_Ynorm, par_Znorm)

Then comes the true simulation step, which uses the information
collected in the previous ones::

    s2Dcd.sim_run(seq, step_max, hard_data, simODS, par_template, seed)

And finally, compute the running time (approximated!) and print `STOP`::

    utili.print_stop(time_start)

Examples with `DeeSse`
--------------------------------------------------------------------------

.. include:: ../examples/01_Strebelle/s2Dcd_run-ex01.py
   :start-line: 2
   :end-line: 19

.. include:: ../examples/02_Strebelle-conditional/s2Dcd_run-ex02.py
   :start-line: 2
   :end-line: 11

.. include:: ../examples/03_Strebelle-many-realizations/s2Dcd_run-ex03.py
   :start-line: 2
   :end-line: 12

.. include:: ../examples/04_Strebelle-many-param/s2Dcd_run-ex04.py
   :start-line: 2
   :end-line: 17





