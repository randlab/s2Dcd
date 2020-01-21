#!/usr/bin/env python3
'''
Example 1: Strebelle's TI along two directions
*************************************************************************

:file name: `s2Dcd_run-ex01.py`

:directory: `/examples/deesse/01_Strebelle`.

In this example, the `s2Dcd` method is applied to the celebrated
Strebelle's 2D training image to obtain a 3D simulation.
Here we use the same training images along the directions normal to the axis
*z* and *y*.
   
No conditioning data are considered.

.. note::
    Another possibility could be, for example, to use the two TI normal to
    axis *x* and *y*. In this case the hypothesis about the 3D geometry are
    of course completely different.

'''

# Import "standard" modules
import os
import numpy as np

# Import non-standard modules
import s2Dcd.s2Dcd as s2Dcd
import s2Dcd.deesse as mpds_interface
import s2Dcd.utili  as utili

# If you can use a parallel version of the MPS simulation core,
# define here the number of threads.
s2Dcd.nb_threads = 2

# Print header and start counting time
time_start = utili.print_start() 

# The random seed, useful to get repeatable results
seed = np.random.RandomState(456833)

# Load the default simulation parameters from a template file
template_in = "template.in"
par_template = mpds_interface.Param(file_name=template_in)

# Create with the same files some templates for the diffent directions.
# If there is no TI along a given direction, use the keyword "None".
par_Xnorm = mpds_interface.Param(file_name=template_in)
par_Xnorm.tis[0].file_name = None
par_Ynorm = mpds_interface.Param(file_name=template_in)
par_Ynorm.tis[0].file_name = os.path.join("..","data",
                                          "strebelle", "ti_250x1x250.gslib")
par_Znorm = mpds_interface.Param(file_name=template_in)
par_Znorm.tis[0].file_name = os.path.join("..","data", "strebelle", 
                                          "ti_250x250x1.gslib")

# The max number of simulation steps to be performed. You can use it if you
# want to stop the simulation before the 3D domain is completed.
step_max = 100000

# Instantiate an object that contains the dimensions of the simulation grid
simODS = par_template.grid

s2Dcd.print_sim_info(simODS, par_Xnorm, par_Ynorm, par_Znorm)

# Array to store hard data and the result of the simulations. By
# default, the cells that don't contain a simulated (or a default
# conditioning data) are set to the value of the variable
# `s2Dcd.no_data`.
hard_data = s2Dcd.no_data*np.ones(
    (simODS.nx, simODS.ny, simODS.nz), 'float')

# Create the simulation sequence
seq = s2Dcd.create_seq(simODS, par_Xnorm, par_Ynorm, par_Znorm)
    
#
# Simulation
#
s2Dcd.sim_run(seq, step_max, hard_data, simODS, par_template, seed)

# Stop counting time
utili.print_stop(time_start)
