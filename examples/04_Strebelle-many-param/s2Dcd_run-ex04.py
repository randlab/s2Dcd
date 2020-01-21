#!/usr/bin/env python3
'''
Example 4: Changing the simulation parameters
*********************************************************************

:file name: `s2Dcd_run-ex01.py`

:directory: `/examples/04_Strebelle-many-param`.

The modules related to the `s2Dcd` can be used also for other 
purposes, like for example to play with the parameters to be used
for a MPS simulation in an automated way.
You can for example use the modules to run different simulations 
changing automatically the size of the data template.
Here we demonstrate this using the Strebelle's training image.
 
No conditioning data are considered.

'''

# Import "standard" modules
import os
import numpy as np

# Import non-standard modules
import s2Dcd.s2Dcd as s2Dcd
import s2Dcd.deesse as mpds_interface
import s2Dcd.utili as utili 


# If you can use a parallel version of the MPS simulation core,
# define here the number of threads.
# 
s2Dcd.nb_threads = 8

# Print header and start counting time
time_start = utili.print_start() 

# The random seed, useful to get repeatable results
seed = np.random.RandomState(456833)

# Load the default simulation parameters from a template file
template_in = "template.in"
par_template = mpds_interface.Param(file_name=template_in)

# Set up the training images to be used for all the simulations.
# (For other tests, you would want to change it too...)
par_template.tis[0].file_name =  os.path.join("..","..","data","strebelle", 
                                              "ti_250x250x1.gslib")


# An example of how to access the elements of the parameters class.
# Here it is useful to check the order and the ID of the data templates
# defined in the "template.in" file.
for i, dt in enumerate(par_template.dts):
    print("        Data template ID: %d, position %d" %( dt.id, i))


# # Instantiate an object that contains the dimensions of the simulation grid
# simODS = par_template.grid

# s2Dcd.print_sim_info(simODS, par_Xnorm, par_Ynorm, par_Znorm)

# # Array to store hard data and the result of the simulations. By default, the 
# # cells that don't contain a simulated (or a default conditioning data) are
# # set to the value of the variable `s2Dcd.no_data`. 
# hard_data = s2Dcd.no_data * np.ones((simODS.nx, simODS.ny, simODS.nz), 'float')

# # Create the simulation sequence
# seq = s2Dcd.create_seq(simODS, par_Xnorm, par_Ynorm, par_Znorm)
    
# #
# # Simulation
# #
# s2Dcd.sim_run(seq, step_max, hard_data, simODS, par_template, seed)

# Stop counting time
utili.print_stop(time_start)


        

    
    
 
