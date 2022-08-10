'''
:license:
    This file is part of s2Dcd.

    s2Dcd is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
 
    s2Dcd is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with s2Dcd.  If not, see <https://www.gnu.org/licenses/>.

:Purpose:
    A module to apply the `s2Dcd` multiple-point simulation approach.
    At the moment it is implemented with the `DeeSse` MPS simulation engine,
    but could be adapted easily to other MPS engines.
 
:File name:
    ``s2Dcd.py``

:Version:
    1.0 , 2022-04-14 :
        * Adapted to the new version of the DeeSse using the geone interface.
    0.9.9 , 2018-01-23 :
        * Solved one bug that prevented simulation perpendicular
          to axis *x*.

    0.9.8 , 2017-12-15 :
        * Solved one but related to casting some index.
    0.9.6 , 2014-08-12 :
        * Try to include the auxiliary variables treatment.
        * Changed name of the module related to the DS.
    0.9.5 , 2013-11-15 :
        * Converted to Python3 with `2to3`.
        * Corrected a bug in the function `matrioska_interval`.
        * Adapted to the new version of gslibnumpy. 
    0.9.4 , 2012-09-05 :
        * Last version before the movement to Python3
    0.9.3 , 2012-09-04 :
        * Tested on a simple case study with one thread and no
          auxiliary variables.
    0.9.2, 2012-05-02 :
        * Defined some variables to describe the range of integer values
          that the seed should take.
        * Corrected a bug related to bad parenthesis (thank you Andrea!)
    0.9.1,  2012-05-02 :
        * Added an option to better define when a problem is for categorical
          variables or for continuous ones.
        * Solved a bug (forgotten to print the `in` file in the case
          of multiple threads).
    0.8,  2012-04-02 :
        * Modified some structures in order to include them in a
          module containing some tools to deal with structured
          grids... in order to prepare the interface to the MPDS code.
    0.7,  2012-02-24 :
        * Using the variable :py:data:`no_data` the default value for
          the not yet simulated nodes can be re-defined with more
          flexibility.

    0.6,  2012-02-23 :
        * Added a functionality allowing to define a random seed.

        * Added a function to create a simple *standard* simulation
          sequence.

    0.5,  2012-02-16 :
        * added the class :py:class:`SeqStep` to  simplify notation.
    0.4,  2012-02-16 :
        * cleaned up some procedures and documentation.
    0.3, :
        * implemented the treatment of auxiliary variables.
    0.2, :
        * some improvements...
    0.1, :
        * first version.

:Authors:
    Alessandro Comunian.


Usage
-----

See the examples in the corresponding directory and the documentation
of the main functions.


Limitations
-----------
auxiliary variables:
   For the moment they are implemented only for the ``implala`` MPS
   engine.

number of variables:
    `MPDS` can be used with only one variable.

TODO
----

* Include the usage of auxiliary variables when using all the MPS
  engines.
* Implement all the features that can be described in the parameter
  files.
* Deal with the connection problem that can raise when the licence
  cannot be verified or the connection is slow.
* Improve the output and input format, adding for example the SGeMS
  binary input/output format (look for the format in mGstat, for
  example)
* For some simulation engines (like the one based on the snesim/impala engines), 
  that can store information coming from threes/lists in and external file, it can be 
  useful to allow the re-use of these.
* Add some functionality to delete all the output files created,
  something like *DEBUG* or *VERBOSE* mode.


'''


# Import "standard" modules
import glob
import random
import time
import numpy
import os
import copy
import sys
import numpy as np

# Import geone
import geone as gn

# Import "custom" modules
# import s2Dcd.deesse as ds_interface REPLACED BY THE GEONE "deesseinterface"
# import s2Dcd.gslibnumpy as gslibnumpy
import s2Dcd.utili as utili
# import s2Dcd.grid as grid
import s2Dcd.ext as ext


# Name of the executables of the MPS simulation code.

# Note that only the file name (without extension) is required, as the
# extension will be added (Windows only) in a 2nd automated step.
mps_exes = { 
    'XXX': "XXX",
    'XXX_MPI': "None",
    'XXX_OMP': "XXX_OMP",
    'MPDS': "deesse",
    'MPDSOMP': "deesseOMP"}

# The default root name of the result
result_file_root = "result"

# The default number of processing threads to be used by the MPS
# simulation
nb_threads = 1

# Use this variable to select the default parallel engine (if available).
# "MPI" or "OMP" values are allowed. When the number of threads is
# detected and > 1, then the corresponding "thread flavour" is
# implemented.
type_threads = "OMP"

# Default value for the seed, to be used when a `RandomState` is not
# defined.  However, as a genereal rule it is better to always define
# a value for `RandomState`.
seed_default = 321654687235

# The integer value that represents a not yet simulated facie.
no_data = -1

# When the function randint is used to generate a seed, these two
# values define the range of the generated seeds.
min_rand_int = 10e+5
max_rand_int = 10e+8

# A reduction factor user to "transform" the size of the 3D template
# (available only for some simulation engine) from "3D" to "quasi-2D"
HALF_SIZE_RED = 1e-3

# Root name of the current *.in file
IN_FILE_ROOT_SIM="curr_sim"
# Root name of the current *.in file when simulating lists
IN_FILE_ROOT_LIST="curr_list"
# Root name of the current hard data conditioning files
HD_FILE_ROOT="curr_hd"
# Root name of the current mask data
MASK_FILE_ROOT="curr_mask"
# A flag to activate the 3D mode for the auxiliary variable
#AUX_VAR_FULL3D=False
AUX_VAR_FULL3D=True

def all_segms0(segms):
    '''
    Check if a list of segments contains at least one segment with
    length > 0.
    
    Parameters:
        segms:
            A list of tuples representing 1D segments.

    Returns:
        True if all the segments contained into the list have length = 0,
        False if at least one segment has a length > 0.

    Example
    
    >>> all_segms0([(1,4), (2,7), (20,20)])
    False
    >>> all_segms0([(1,1), (7,7), (20,20)])
    True

    '''
    for segm in segms:
        if (segm[1]-segm[0]) > 0:
            return False
    return True


def split_segms( segms, mid_list=[]):
    '''
    A recursive function to print out a "matrioska-like" sequence of
    integers.

    Parameters
        segms : list of tuples of type (a,b)
            A list containing a number of tuples representing a segment.
        mid_list: list of integers, optional
            The list of integers that represents the required output.
            By default this should be called without this argument.

    Returns
        A list containing all the integers in the mixed "matrioska" order.

    Example
        >>> matrioska_interval([(11,26)])
        [11, 26, 18, 14, 22, 12, 16, 20, 24, 13, 15, 17, 19, 21, 23, 25]
    
    .. note:: 
        It is supposed that the coordinates of the segments are all
        integers.

    '''

    if all_segms0(segms):
        # All the segments have a lenght 0. This is the final step.
        return mid_list
    else:
        # At least one segment has a lenght > 0.
        if len(segms) == 1:
            # This is the very first call, where only one segment is
            # provided.  The default values are the extremes of the
            # provided segment.
            mid_list = [segms[0][0], segms[0][1]]

        new_segms=[]
        for segm in segms:
            # Compute the median point
            mid = (segm[1]+segm[0])//2

            # Add it to the final output
            if mid  not in mid_list:
                mid_list.append(mid)
            
            # If the segments are "well formed", then add them to the
            # list of segments.
            if mid >= segm[0]:
                new_segms.append((segm[0],mid))
            if segm[1] >= mid+1:
                new_segms.append((mid+1,segm[1]))
        # Split again the list of segments
        return split_segms(new_segms, mid_list)


def matrioska_interval(points_nb):
    '''
    Computes a "spreaded" list of integers.

    This function is useful to provide a simulation sequence along one
    axis that allows, in priciple, to obtain as much as intersections
    as possible along the other directions.

    Parameters:
        points_nb: integer
            The number of points contained in the list.

    Returns:
        A list of integers containing a "matrioska sequence".
        If there is an error, returns -1.

    Example:
        A typical output of the function is:

           >>> matrioska_interval(8)
           [0, 7, 3, 1, 5, 2, 4, 6]

    .. note::
        * The output sequence always start from 0.
        * The algorithm is probably not efficient, but for small
          number, that is for sizes of simulations (<1000), it can be
          OK.  Moreover, it should be called few times (max 3).
    
    '''

    if points_nb <= 0:
        print('    ERROR in "s2Dcd.matrioska_interval", '
              'points_nb can only be > 0.')
        return -1
    else:
        N = points_nb - 1 # C like numbering

    return split_segms([(0,N)])


class SeqStep(object):
    """
    A class that contains all the information required to create a
    parameters input file for the used MPS engines, that is a *step* of
    the simulation sequence of the approach `s2Dcd`.
    
    NOTE: FOR THE MOMENT ONLY THE DEEESSE IS IMPLEMENTED
    """

    def __init__(self, direct, level, param=None, nthreads=1):
        '''
        Initialization for the class :py:class:`SeqStep`.  All the
        parameters definitions are considered for a given simulation
        step.

        Parameters:
            direct: char, 'x', 'y' or 'z'
                The direction normal to the surface that is simulated.
            level: int
                The corresponding level of the axis defined by
                :py:attr:`direct` where we are simulating. It is an
                integer index in the 3D matrix of the simulation
                domain.
            param: object that contains all the parameters required to 
                run a simulation step of the sequence depending on the 
                given simulation engine. For the moment only the 
                DeeSse is implemented.
                The "None" value is allowed only for testing purposes.
                
        Returns:
            If successful, it creates the object. Otherwise, raises some
            error/warning.

        '''
      
        if direct in ['x','y','z']:
            self.direct = direct
        else:
            raise AttributeError("The value of the attribute *direct* should be in "
                                 "('x', 'y', 'z'). Wrong entered value is: "
                                 "{0}".format(direct))
        self.level = level
        self.param = param
        
        self.ti = param.TI
        self.nthreads = nthreads

        # SEE THE OLDER VERSIONS FOR SUGGESTIONS. THIS COULD BE THE RIGHT
        # PLACE TO LOOK FOR FILES CONTAINING AUXILIARY VARIABLES OR
        # WHERE TO SET UP SOME PROPERTIES (FOR EXAMPLE THE DIMENSIONS OF THE
        # SEARCH NEIGHBORHOOD =0 ALONG THE RIGHT DIRECTION...)

        
    
    def __str__(self):
        """
        Print the content of a sequence step
        """
        out = ("    *** Sequence step info ***\n" +
                "- level {0.direct} = {0.level:5d}".format(self))
        return out


    def simul(self, step, step_max, hard_data, simODS, rcp_lists, verbose=0):
        '''
        Run a MPS simulation for a simulation step of a sequence.

        Parameters:
            step: int
                The current simulation step.
            step_max: int
                The maximum simulation step. Note that the simulation
                domain can be often filled before this step. In this
                case, the simulation procedure is stopped.  Defining
                this parameter smaller than the total number of
                simulations expected to fill the domain can be useful
                for doing some preliminar test.
            hard_data: numpy array
                A numpy array containig all the simulation. The value
                of the global variable :py:data:`no_data` is used in
                locations not yet simulated.
            simODS: object of type :py:class:`Grid`
                Contains all the info concerning the simulation
                domain.
            rcp_lists: flag
                NOTE: THIS IS A LEGACY PARAMETER WHICH WAS USEFUL WHEN
                THE IMPALA MPS ENGINE IS USED.
                Useful if one wants to re-compute the list at each 
                simulation step. This parameter has sense only for the
                `Impala` family of MPS simulation engines.

        Returns:
            * Create some files required for the simulation.
            * Update the content of the :py:class:`geone.img.Img`
              including new simulated nodes.
             
            If the selected section is full of hard data, simply
            returns 1 without creating files and without running the
            MPS simulation.  If there is an error, returns -1.
        
        '''

        
        # This is only for pretty printing purposes...
        intlen = len(str(step_max))
        
        # Check which MPS simulation engine is in use.
        # NOTE THAT AT THE MOMENT THE ONLY MPS SIMULATION ENGINE IMPLEMENTED
        # IS THE DS.
        if isinstance(simODS, gn.deesseinterface.DeesseInput):
            print('    Seq.step {0:{2}d} of {1:{2}d} (max)'.format(
                step+1, step_max, intlen), end=" ")            
        else:
            print("    WARNING: unknown MPS simulation engine!")
            return -1
        
        print("- slice simulated: {0.direct:s} = {0.level:4d}"
              "".format(self))  
    
       
        #   
        # Set up the slice which contains the hard data
        #
        # WARNING: 
        #     1) FOR THE MOMENT THERE IS THE LIMITATION OF USING ONLY ONE TI.
        #        HEREINAFTER THEN THE VALUE "0" SHOULD BE MANAGED ACCORDINGLY
        #
        #     2) WHEN ACCESSING THE TI WITH "VAL", THEN BE...    
        curr_hd = copy.deepcopy(hard_data)
        curr_par = copy.deepcopy(self.param)
        if self.direct == "x":
            # Set the hard data to be used for the simulation
            curr_hd.nx = 1
            curr_hd.val = np.reshape(hard_data.val[0,:,:,self.level],
                                 (1,1,simODS.ny,simODS.nz))
            # Set the parameters for the simulation
            # NOTE: DOUBLE CHECK HERE IF "level" REFERS TO THE TRUE COORDINATE
            #       OR TO THE INDEX WITHIN THE MATRIX. DOUBLE CHECK ALSO FOR
            #       THE OTHER DIRECTIONS.
            curr_par.ox = self.level
            
        elif self.direct == "y":
            curr_hd.ny = 1
            curr_hd.val = np.reshape(hard_data.val[0, :, self.level, :], 
                                 (1,simODS.nx,1,simODS.nz), order="C")
            # Set the parameters for the simulation
            curr_par.y = self.level
        elif self.direct == "z":       
            curr_hd.nz = 1
            curr_hd.val = np.reshape(hard_data.val[0, self.level, :, :],
                                 (1,simODS.nx,simODS.ny,1))
            # Set the parameters for the simulation
            curr_par.z = self.level
            
        # gn.imgplot3d.drawImage3D_surface(curr_hd, excluded_value=np.nan, text="s2Dcd", show_edges=True)            
        # gn.imgplot3d.drawImage3D_slice(curr_hd, slice_normal_y=0,
        #                                show_bounds=True,   # add bounds (axis with graduation)
        #                                text='TI')
        # gn.imgplot3d.drawImage3D_surface(curr_hd, filtering_value=[1], 
                                         # text='TI')
        
        if isinstance(curr_hd, gn.img.Img):
            curr_par.dataImage = np.array([curr_hd])
        else:
            print("ERROR, wrong input type for dataImage!")
        # gn.img.writeImageVtk(curr_hd, "s2Dcd_step{0:06d}.vtk".format(step+1), data_type='int')

        if not isinstance(curr_par, gn.deesseinterface.DeesseInput):
        # LATER ADD A MESSAGE IN THE LOG FILE...
            print("ERROR, wrong input type for DeesseInput parameters!")
        
        # Check if in the slice there are hard data.
        nb_nan = np.sum(np.isnan(curr_hd.val))
        if(nb_nan == 0):
            # No NaN, all the nodes have a valid value, no need to simulate.
            return 1
        # print("ACOMUNIAN")
        # print(curr_par)
        # print(curr_par.dataImage)
        deesse_out = gn.deesseinterface.deesseRun(curr_par, 
                                                  nthreads=self.nthreads,
                                                  verbose=verbose)        
        gn.img.writeImageVtk(curr_hd, "s2Dcd_step{0:06d}.vtk".format(step+1),
                             data_type='int', missing_value=-9999999)
        
        # print(deesse_out)
        
        
        
#         error_code = os.system(runSim)
#         if error_code == ds_interface.LIC_ERR_CODE:
#             for errcd in range(ds_interface.NB_LIC_WAIT):
#                 print("    WARNING: Licence problems with the 'deesse'")
#                 print("        Waiting {0} sec and retry...".
#                       format(ds_interface.LIC_WAIT_TIME))
#                 time.sleep(ds_interface.LIC_WAIT_TIME)
#                 error_code = os.system(runSim)
#                 if error_code == 0:
#                     print("        OK, licence problem solved!")
#                     break
            
            

#         if(error_code != 0):
#             error_message = "    ERROR: problems in %(mps_exe)s, "\
#                 "error_code=%(error_code)d" % locals()
#             sys.exit(error_message)
        
#         # Read and add the results to the storage hard data matrix
#         if mps_core in ["Impala","ImpalaMPI","ImpalaOMP"]:
#             new_hd_file = (os.path.splitext(curr_in_file)[0] + "_00000" +
#                            ext.VTK)
#             new_hd = vn.vtk2numpy(new_hd_file, verbose=False)
#         elif mps_core in ["MPDS","MPDSOMP"]:
#             new_hd_file = os.path.splitext(curr_in_file)[0] + ext.GSLIB
#             new_hd =  gslibnumpy.gslib2numpy_onevar(new_hd_file,
#                                                     verbose=False)
        
#         # Add the new hard data to the "storage" array
#         add_hd(self, new_hd, hard_data, simODS)

            # print(deesse_out) # ACOMUNIAN
        if self.direct == "x":
            # print(deesse_out["sim"][0].val.shape)
            hard_data.val[0,:,:,self.level] = deesse_out["sim"][0].val[0,:,:,0]
        elif self.direct == "y":            
            # print(deesse_out["sim"][0].val.shape)
            
            hard_data.val[0,:,self.level,:] = deesse_out["sim"][0].val[0,:,0,:]
        elif self.direct == "z": 
            # print(deesse_out["sim"][0].val.shape)
            hard_data.val[0,self.level,:,:] = deesse_out["sim"][0].val[0,0,:,:]


#     def create_list(self, simODS, par_template):
#         '''
#         Run Impala to create the lists required along a given
#         direction.

#         Parameters:
#             simODS: object of type :py:class:`Grid`
#                 Contains all the info concerning the simulation
#                 domain.
#             par_template: class :py:class:`impala_interface.Param`.
#                 The parameters for creating the *.in* file for the
#                 simulations.

#         Returns:
#            * Create for each level of multigrid a list in binary format.
#         '''

#         # Check which MPS simulation core is in use
#         if nb_threads > 1:
#             if type_threads == "MPI":
#                 mps_core = 'ImpalaMPI'
#             else:
#                 mps_core = 'ImpalaOMP'
#         elif nb_threads == 1:
#             mps_core = 'Impala'
#         else:
#             print("    Error, wrong number of computation threads (%d)" % 
#                   ( nb_threads ))
#             return -1

#         mps_exe = mps_exes[mps_core]
    
#         in_file_root = IN_FILE_ROOT_LIST # Root name for the files *.in
#         in_file_ext = ext.IMPALA_IN
        
#         print("    Creation of a list along the coordinate %s"
#               %(self.direct))
#         curr_in_file = "%s-%s%s" % (in_file_root, self.direct, in_file_ext)
    
#         # Add the suffix "exe" to the executable if we are running Windows
#         if os.name == 'nt':
#             mps_exe = mps_exe + ext.EXE
        
#         # Create the "*.in" file used for the current simulation step
#         if mps_core == 'Impala':
#             create_in_file4Impala(par_template, self, None, curr_in_file)
#             # ADAPT LATER FOR THE CASE OF AUXILIARY VARIABLES
#             # adaptAuxVarFile(self, simODS)
#             # Create the execution string
#             runSim = ("%s %s %s_ 0 0 0 0 0 1 list_%s   %s.log > %s.out" %  
#                       (mps_exe, curr_in_file, curr_in_file[:-3], 
#                        self.direct, 
#                        curr_in_file[:-3], curr_in_file[:-3]))
#             #print(runSim)
#         elif mps_core == 'ImpalaOMP':
#             create_in_file4Impala(par_template, self, None, curr_in_file)
#             # ADAPT LATER FOR THE CASE OF AUXILIARY VARIABLES
#             # adaptAuxVarFile(self, simODS)
#             # Create the execution string
#             runSim = ("%s %d %s %s_ 0 0 0 0 0 1 list_%s  "
#                       " %s.log > %s.out" %  
#                       (mps_exe, nb_threads, curr_in_file,
#                        curr_in_file[:-3], self.direct , 
#                        curr_in_file[:-3], curr_in_file[:-3]))
#             #print(runSim)
#         elif mps_core == 'ImpalaMPI':
#             # Create the "*.in" file for the current simulation step
#             create_in_file4Impala(in_par, self, curr_hd_file, curr_in_file)
#             # ADAPT LATER FOR THE CASE OF AUXILIARY VARIABLES
#             # adaptAuxVarFile(self, simODS)
#             # Create the execution string
#             runSim = ("mpirun -np %d %s %s %s_ 0 0 0 0 0 1 "
#                       "list_%s  %s.log > %s.out" %  
#                       (nb_threads, mps_exe, curr_in_file,
#                        curr_in_file[:-3], self.direct , 
#                        curr_in_file[:-3], curr_in_file[:-3]))
#         else:
#             print("    ERROR: Wrong MPS engine selected")
                    
#         error_code = os.system(runSim)

#         if(error_code != 0):
#             error_message = ("    ERROR: problems in %(mps_exe)s,"
#                              " error_code=%(error_code)d" % 
#                              locals())
#             sys.exit(error_message)
        


def add_hd(seq_step, new_hd, hard_data, simODS):
    '''
    Add a numpy array containing the data simulated at a given
    sequence step to the array that stores the simulated nodes.

    Parameters:
        seq_step: instance of the class :py:class:`SeqStep`
            Information about the current simulation step.
        new_hd: numpy array (2D)
            The new simulated section.
        hard_data: numpy array (3D)
            The contained of all the simulated data.
        simODS:  object of type :py:class:`Grid`
            Information about the simulation grid.

    Returns:
        Fit the content provided with `new_hd` into the right position
        into the array `hard_data`.

    .. note::
        * A the moment the fact that a gslib file can contain multiple 
          variables is handled with a quick and dirty trick.
    '''

    # WARNING: THIS IS A QUICK AND DIRTY TRICK! BUT SOMEHOW IT MIGHT
    # BE USEFUL TO KEEP THIS OPTIONS FOR FUTURE DEVELOPMENTS TAKING
    # INTO ACCOUNT FOR MULTIPLE VARIABLES...

    if seq_step.direct == "x":
        hard_data[seq_step.level, :, :] = (
            numpy.reshape(new_hd[0], (simODS.ny, simODS.nz)))
    elif seq_step.direct == "y":
        hard_data[:, seq_step.level, :] = (
            numpy.reshape(new_hd[0], (simODS.nx, simODS.nz)))
    elif seq_step.direct == "z":       
        hard_data[:, :, seq_step.level] = (
            numpy.reshape(new_hd[0], (simODS.nx, simODS.ny)))
        
    
     

def sim_run(seq_steps, step_max, hard_data, simODS, nthreads=1,
            rcp_lists=False):
    '''
    Run a `s2Dcd` simulation.

    Parameters:
        seq_steps: list object of type :py:class:`SeqStep`
            A list containing all the information of each simulation step.
        step_max: int
            The maximum simulation step. Note that the simulation
            domain can be often filled before this step.
            In that case, a check is made and the simulation should end in
            advance.
        hard_data: object of type :py:class:`geone.img.Img`
            An empty TI, containing `np.nan` at all the locations where a 
            hard data is not available.
            The final result will be contained here.
        simODS: object of type :py:class:`geone.deesseinterface.DeesseInput`
            Contains all the info concerning the simulation domain.
        in_par: string or a object of type `ds_interface.Param`
            The name of the file containing the template of the input
            parameters file for `Impala` or the object containing all
            the parameters required for a `MPDS` simulation.
        seed: instance of the class :py:class:`RandomState`
            To keep track of the random seed and create simulations
            that can be reproduced using the same seed.
        rcp_lists: flag (default False)
            This flag is useful if one wants to re-compute the list
            at each s2Dcd simulation step. When True it is a waste of time
            if the TI along the same direction is always the same. It is
            only useful when different training images are used for
            the simuations along the same direction.
            

    Returns:
        A VTK or a GSLIB file containing the results of the simulation.
        If the max number of interation is reached, a result file
        containing "no data" values is printed out.
        
    '''

    # Define the name of the output file
    # if isinstance(simODS, geone.deesseinterface.DeesseInput):
    #      res_file = res_file_root + ext.GSLIB
 
    print("\n    *** Simulation starts *** ")
    stp_max = min(len(seq_steps), step_max)
    for i, seq_step in enumerate(seq_steps[0:stp_max]):

        simul_out = seq_step.simul( i, stp_max, hard_data, simODS, rcp_lists)
        # gn.img.writeImageVtk(hard_data, "simul_step{0:06d}.vtk".format(i+1), data_type='int')

        # FOR THE MOMENT IT CONSIDERS THAT THE 3D DOMAIN IS FULL, NO
        # MASK OPTION IS HANDLED

        if simul_out is None:
            continue
        else:
            if simul_out == 1:
                print('    - no need for simulation, section full of'
                      ' hard data.')
                continue
            else:
                print('- ERROR in the simulation. simul returned',
                      simul_out)
                break
            
    # Check if the hard_data matrix is all informed
    # (<=> the simulation is completed)
    
    # Compute the number of np.nan in the results
    sum_nan = np.sum(np.isnan(hard_data.val))
    if ( sum_nan > 0):
        # Could be useful for debug purposes
        print('\n    WARNING: Some nodes were not simulated.')
        print("        Number of NaN nodes: {0}".format(sum_nan))
    else:
        print('\n    All the domain completed!\n')

    print("    *** Saving the result ***")
#     print_result(hard_data, simODS, res_file)
# #            break
 

     

# def create_in_file4Impala(in_par, seq_si, file_cond, file_name_sim,
#                           seed=None):
#     '''
#     Create a parameters input file for *Impala*. A big part of the
#     information is read directly from the template input file for
#     *Impala*.

#     Parameters:
#         in_par: object of type `impala_interface.Param`.
#             All the information contained in the template file.
#         seq_si: object of the class :py:class:`SeqStep`.
#             All the informations required to create a simulation for the
#             current simulation step. See the class :py:class:`SeqStep`
#             for details.
#         file_cond: string or Nonenp.reshape(data_x1["facies"].astype("int").values, (nz, ny, nx))
#             Name of the conditioning file, if None the simulation is
#             considered as non conditional.
#         file_name_sim: string
#             Name of the input paramter file for *Impala*.
#         seed: instance of the class :py:class:`RandomState`, optional
#             To keep track of the random seed and create simulations
#             that can be reproduced using the same seed.  A definition
#             of the seed is not required when the MPS core is called
#             only for the generation of the lists. Therefore, in this
#             case the value of seed can be None, and the value of the
#             module variable `seed_default` is used.
    
#     Returns:
#         A `*.in` file containing the parameters for running `Impala`.

#     .. note::
#         In order to increase the variability of the simulations, a new
#         random seed is generated for each simulation file that is
#         created if a :py:class:`RandomState` instance is provided with
#         the paramter `seed`. Otherwise, the value of the global
#         variable `seed_default` is used.

#     '''
#     curr_par = copy.deepcopy(in_par)
#     curr_par.tis[0].file_name = seq_si.param.tis[0].file_name
#     if seq_si.param.tis[0].nb_aux_var > 0:
#         curr_par.tis[0].aux_var_file_name[0] \
#             = seq_si.param.tis[0].aux_var_file_name[0]

#     #
#     # Set the parameters for this specific direction
#     #
#     if seq_si.direct == 'x':
#         curr_par.grid.ox = curr_par.grid.ox + seq_si.level*curr_par.grid.dx
#         curr_par.grid.nx = 1
#         # Reduce the size of the templates along the direction that is
#         # not simulated
#         for dt in curr_par.dts:
#             dt.half_size[0] = str(HALF_SIZE_RED*float(dt.half_size[0]))

#         if curr_par.nb_aux_var > 0:
#             # Create the auxiliary variable map for the current
#             # simulation step
#             if AUX_VAR_FULL3D:
#                 curr_par.aux_var_files[0] = vn.vtk_slice(
#                     seq_si.param.aux_var_files[0], 'x', seq_si.level) 
#             else:
#                 vn.vtk_trasl(seq_si.param.aux_var_files[0], (curr_par.grid.ox,curr_par.grid.oy,curr_par.grid.oz))

#                 curr_par.aux_var_files[0] = copy.deepcopy(seq_si.param.aux_var_files[0])
#         if curr_par.mask_map_flag == "MASK_MAP_ON":
#             # Create a slice in the main 3D MASK file to be used in
#             # the 2D simulation
#             curr_par.mask_map_file_name = vn.vtk_slice(
#                 seq_si.param.mask_map_file_name, 'x', seq_si.level)

#     elif seq_si.direct == 'y':
#         curr_par.grid.oy = curr_par.grid.oy + seq_si.level*curr_par.grid.dy
#         curr_par.grid.ny = 1
#         for dt in curr_par.dts:
#             dt.half_size[1] = str(HALF_SIZE_RED*float(dt.half_size[1]))
#         if curr_par.nb_aux_var > 0:
#             if AUX_VAR_FULL3D:
#                 curr_par.aux_var_files[0] = vn.vtk_slice(
#                     seq_si.param.aux_var_files[0], 'y', seq_si.level)
#             else:
#                 vn.vtk_trasl(seq_si.param.aux_var_files[0], (curr_par.grid.ox,curr_par.grid.oy,curr_par.grid.oz))
#                 curr_par.aux_var_files[0] = copy.deepcopy(seq_si.param.aux_var_files[0])

#         if curr_par.mask_map_flag == "MASK_MAP_ON":
#             # Create a slice in the main 3D MASK file to be used in
#             # the 2D simulation
#             curr_par.mask_map_file_name = vn.vtk_slice(
#                 seq_si.param.mask_map_file_name, 'y', seq_si.level)

#     elif seq_si.direct == 'z':
#         curr_par.grid.oz = curr_par.grid.oz + seq_si.level*curr_par.grid.dz
#         curr_par.grid.nz = 1
#         for dt in curr_par.dts:
#             dt.half_size[2] = str(HALF_SIZE_RED*float(dt.half_size[2]))
#         if curr_par.nb_aux_var > 0:
#             if AUX_VAR_FULL3D:
#                 curr_par.aux_var_files[0] = vn.vtk_slice(
#                     seq_si.param.aux_var_files[0], 'z', seq_si.level)
#             else:
#                 vn.vtk_trasl(seq_si.param.aux_var_files[0], (curr_par.grid.ox,curr_par.grid.oy,curr_par.grid.oz))
#                 curr_par.aux_var_files[0] = copy.deepcopy(seq_si.param.aux_var_files[0])
#         if curr_par.mask_map_flag == "MASK_MAP_ON":
#             # Create a slice in the main 3D MASK file to be used in
#             # the 2D simulation
#             curr_par.mask_map_file_name = vn.vtk_slice(
#                 seq_si.param.mask_map_file_name, 'z', seq_si.level)

#     else:
#         print('    Error in `create_in_file4Impala`')
#         return -1

#     # Define a value for the seed
#     if seed:
#         curr_par.seed = seed.randint(min_rand_int,max_rand_int)
#     else:
#         curr_par.seed = seed_default

#     if file_cond:
#         curr_par.conditional_flag = "CONDITIONAL_ON"
#         # CONNECTIVITY NOT IMPLEMENTED FOR THE MOMENT
#         curr_par.connectivity_flag  = "CONNECTIVITY_OFF"
#         curr_par.conditional_file_name = file_cond
#       #  curr_par.data_pointset_nb = 1 
#       #  curr_par.data_pointset_files = []
#       #  curr_par.data_pointset_files.append( file_cond )

#     curr_par.out_file = os.path.splitext(file_name_sim)[0] + ext.VTK
#     curr_par.out_report_file = os.path.splitext(file_name_sim)[0] + ext.TXT
#     curr_par.print_file(file_name_sim)


# def create_in_file4MPDS(in_par, seq_si, file_cond, file_name_sim,
#                         seed=None):
#     '''
#     Create a parameters input file for `MPDS` adapted to the current
#     simulation step and the current conditioning data file.


#     Parameters:
#         in_par: object of type `ds_interface.Param`.
#             All the info about the in template file.
#         seq_si: object of the class :py:class:`SeqStep`.
#             (sequence stepA info) All the informations required to
#             create a simulation for the current simulation step. See
#             the class :py:class:`SeqStep` for details.
#         file_cond: string or None
#             Name of the conditioning file, if None the simulation is
#             considered as non conditional.
#         file_name_sim: string
#             Name of the input paramter file for `MPDS`.
#         seed: instance of the class :py:class:`RandomState`, optional
#             To keep track of the random seed and create simulations
#             that can be reproduced using the same seed.  A definition
#             of the seed is not required when the MPS core is called
#             only for the generation of the lists. Therefore, in this
#             case the value of seed can be None, and the value of the
#             module variable `seed_default` is used.

#     Returns:
#         A `*.in` file for the current simulation.

#     .. note::
#         In order to increase the variability of the simulations, a new
#         random seed is generated for each simulation file that is
#         created if a :py:class:`RandomState` instance is provided with
#         the paramter `seed`. Otherwise, the value of the global
#         variable `seed_default` is used.

#     '''

#     curr_par = copy.deepcopy(in_par)

# #    print(curr_par)
#     curr_par.tis[0].file_name = seq_si.param.tis[0].file_name

#     # Set to 0.0 the search radius along the right direction
#     if seq_si.direct == 'x':
#         for var in curr_par.vars:
#             var.search_r_x = 0.0
# #        curr_par.vars[0].search_r_x = 0.0
#         curr_par.grid.ox = seq_si.level
#         curr_par.grid.nx = 1
#     elif seq_si.direct == 'y':
#         for var in curr_par.vars:
#             var.search_r_y = 0.0
#  #       curr_par.vars[0].search_r_y = 0.0
#         curr_par.grid.oy = seq_si.level
#         curr_par.grid.ny = 1
#     elif seq_si.direct == 'z':
#         for var in curr_par.vars:
#             var.search_r_z = 0.0
# #        curr_par.vars[0].search_r_z = 0.0
#         curr_par.grid.oz = seq_si.level
#         curr_par.grid.nz = 1
#     else:
#         print('    Error in `create_in_file4MPDS`')

#     curr_par.seed = seed.randint(min_rand_int, max_rand_int)

#     if file_cond:
#         curr_par.data_pointset_nb = 1 
#         curr_par.data_pointset_files = []
#         curr_par.data_pointset_files.append( file_cond )

#     curr_par.out_file = file_name_sim[:-3] + ext.GSLIB
#     curr_par.out_report_file = file_name_sim[:-3] + ext.TXT
#     curr_par.print_file(file_name_sim)


    
def addVtk2HdArchive(hard_data, facies):
    """
    A function to add all the VTK files contained into the current
    directory into an hard data archive file (that is a 3D matrix
    which once filled in will be the "final" simulation).

    Parameters:
        hard_data: string
            Name of the file which contains all the conditioning data.
        facies: string
            A string containing the facies which are considered.

    .. warning::
        This function is obsolete. It can be time consuming if a lot
        of VTK files are present in the working directory.

    """
    # REMOVE IT?

    # Check if there was a previous simulated vtk which can be converted
    # into hard data (and eventually convert it)
    sim_files = glob.glob("*.vtk")
    if(len(sim_files) == 0):
        print('      No VTK files to be added to the archive available')
    else:
        for sim_file in sim_files:
            print('      Adding the file "%s" to the hard data'
                  ' archive...' % sim_file)
            data, vtkReaderOutput = vn.vtk2numpy(sim_file)
            # Compute where "data" should be added to "hard_data"
            axis, coord = GetWhereToAddData(vtkReaderOutput)


def GetWhereToAddData(vtkReaderOutput):
    '''
    Get where the data file should be added into the hard data
    archive.
    
    Parameters:
        vtkReaderOutput:
            The result of a "GetOutput()" from a
            "vtkStructuredPointsReader" object.

    Returns:
        * The axis normal to the plane where the simulation was
          performed in the format character, that is *x*, *y* or *z*.
        * the "coordinate" (in the reference system of the simulation,
          that is related to the "pixel" which have to be simulated in
          the matrix of hard data) of the slice where the new
          simulated data should be added.

    .. warning::
        This function is obsolete.
    '''
    # REMOVE IT?
    
    axis = ""
    nx, ny, nz = vtkReaderOutput.GetDimensions()
    
    if((nx + ny + nz - max(nx, ny, nz) < 3)):
        "   ERROR: two or more dimensions are =1"
        return axis, -1
    
    if(nx == 1):
        axis = "x"

                   
def numpy2hd4Impala(simODS, hd, hd_file, seq=None):
    '''
    Given a numpy array, prints into an output file all the values in
    the hard data file format of *Impala*.

    Parameters:
        simODS: object of type :py:class:`Grid`
            This object contains all the dimension required to define a
            simulation domain.
        hd: numpy array
            A numpy array containing all the data which to be
            converted into hard data. The array contains
            :py:data:`no_data` where there are not conditioning data.
        hd_file: string
            Name of the file where the hard data will be printed.
        seq: object type :py:class:`SeqStep`, optional
            Some information about the current simulation step. 
            If None, then all the *informed* content of the input
            numpy file is saved.
            
    Returns:
        A file containing all the conditioning data required for a given
        simulation step in the *Impala* hard data format.
    
    '''
    
    nx, ny, nz = hd.shape
    
    # Count the number of values which are different from `no_data`,
    # the "no hard data value
    Nhd = numpy.sum(numpy.where(hd != no_data, 1, 0))
    
    f = open (hd_file, 'w')
    
    # Here we suppose that the dimension is in general 3, even if 
    # the simulation are only 2D.
    print("3   // space dimension", file=f)                
    print(Nhd, "   // number of conditioning data", file=f)
    
    if seq:
        if seq.direct == 'x':
            for j in range(ny):
                for k in range(nz):
                    if(hd[0, j, k] != no_data):
                        print(simODS.dx * (seq.level + 0.5) + simODS.ox,
                              "   ", end=' ', file=f)
                        print(simODS.dy * (j + 0.5) + simODS.oy, "   ",
                              end=' ', file=f)
                        print(simODS.dz * (k + 0.5) + simODS.oz, "   ",
                              end=' ', file=f)
                        print(hd[0, j, k], "   -999", file=f)
        elif seq.direct == 'y':
            for i in range(nx):
                for k in range(nz):
                    if(hd[i, 0, k] != no_data):
                        print(simODS.dx * (i + 0.5) + simODS.ox, "   ",
                              end=' ', file=f)
                        print(simODS.dy * (seq.level + 0.5) + simODS.oy,
                              "   ", end=' ', file=f)
                        print(simODS.dz * (k + 0.5) + simODS.oz, "   ",
                              end=' ', file=f)
                        print(hd[i, 0, k], "   -999", file=f)        
        elif seq.direct == 'z':
            for i in range(nx):
                for j in range(ny):
                    if(hd[i, j, 0] != no_data):
                        print(simODS.dx * (i + 0.5) + simODS.ox, "   ",
                              end=' ', file=f)
                        print(simODS.dy * (j + 0.5) + simODS.oy, "   ",
                              end=' ', file=f)
                        print(simODS.dz * (seq.level + 0.5) + simODS.oz,
                              "   ", end=' ', file=f)
                        print(hd[i, j, 0], "   -999", file=f)
    else:
        # A sequece was not provided, therefore print all data which
        # were simulated
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if(hd[i, j, k] != no_data):
                        print(simODS.dx * (i + 0.5) + simODS.ox, "   ",
                              end=' ', file=f)
                        print(simODS.dy * (j + 0.5) + simODS.oy, "   ",
                              end=' ', file=f)
                        print(simODS.dz * (k + 0.5) + simODS.oz, "   ",
                              end=' ', file=f)
                        print(hd[i, j, k], "   -999", file=f)
         
                
    f.close()
                 

def numpy2hd4MPDS(simODS, hd, seq):
    '''
    Given a numpy array, prints into to a `xyz_data` file to be used
    by the function `gslibnumpy.numpy2gslib_points` to create a file
    with the conditioning points for the simulation for the MPS code
    `MPDS`.

    Parameters:
        simODS: object of type :py:class:`Grid`
            This object contains all the dimension required to define a
            simulation domain.
        hd: numpy array
            A numpy array containing all the data which to be
            converted into hard data. The array contains
            :py:data:`no_data` where there are not conditioning data.
            seq: object type :py:class:`SeqStep` Some information
            about the current simulation step.
            
    Returns: 
        A "xyz_data" file (a tuple of numpy arrays) containing the
        information about the coordinates and the data.
    
    .. note:
        By default the output points are in 3D.
    '''
    
    nx, ny, nz = hd.shape
    
    # Count the number of values which are different from `no_data`,
    # the "no hard data value".
    Nhd = numpy.sum(numpy.where(hd != no_data, 1, 0))

    x = numpy.zeros(Nhd)
    y = numpy.zeros(Nhd)
    z = numpy.zeros(Nhd)
    data = numpy.zeros(Nhd, dtype=hd.dtype)
    iii = 0

    #
    # CHECK IF THE CONVENTION FOR THE COORDINATES IS OK!
    #
    if seq.direct == 'x':
        for j in range(ny):
            for k in range(nz):
                if(hd[0, j, k] != no_data):
                    x[iii] = simODS.dx * (seq.level + 0.5) + simODS.ox
                    y[iii] = simODS.dy * (j + 0.5) + simODS.oy
                    z[iii] = simODS.dz * (k + 0.5) + simODS.oz
                    data[iii] = hd[0, j, k]
                    iii = iii + 1
    elif seq.direct == 'y':
        for i in range(nx):
            for k in range(nz):
                if(hd[i, 0, k] != no_data):
                    x[iii] = simODS.dx * (i + 0.5) + simODS.ox
                    y[iii] = simODS.dy * (seq.level + 0.5) + simODS.oy
                    z[iii] = simODS.dz * (k + 0.5) + simODS.oz
                    data[iii] = hd[i, 0, k]
                    iii = iii + 1
    elif seq.direct == 'z':
        for i in range(nx):
            for j in range(ny):
                if(hd[i, j, 0] != no_data):
                    x[iii] = simODS.dx * (i + 0.5) + simODS.ox
                    y[iii] = simODS.dy * (j + 0.5) + simODS.oy
                    z[iii] = simODS.dz * (seq.level + 0.5) + simODS.oz
                    data[iii] = hd[i, j, 0]
                    iii = iii + 1

    return (x, y, z, data)
                
                
def print_result(hd, simODS, file_name):
    '''
    Print into an ouput file the simulation output.

    Parameters:
        hd: numpy array
            The numpy (3D) array containing all the simulated nodes.
        simOSD: object of type :py:class:`Grid`
            Contains the dimensions of the simulation domain.
        file_name: string
            Name of the output file. The extension can be 
            "vtk" or "gslib".

    Returns:
        Create the output file.
        If some error occurs, returns -1.
    '''

    # NOTE: It could be improved using the interface VTK-Python
    # provided by ParaView... but then there is the problem to be
    # completely dependent on these modules...

    file_name_root, file_ext = os.path.splitext(file_name) 

    #    print("debug file_name:", file_name )
    #    print("debug FILE EXT:", file_ext )
    if file_ext == '.vtk':
        vn.numpy2vtk(hd, file_name, simODS)
    elif file_ext == '.gslib':
        gslibnumpy.numpy2gslib((hd,), file_name)
    else:
        print('    Error in print_result, unknonw file type: %s'
              % (file_ext))
        return -1
    
    
def adaptAuxVarFile(seq, simODS):
    """
    Change the dimension in the VTK file containing the auxiliary
    variable map in order to adapt it for the current 2D simulation.

    Parameters:
        seq: object of type :py:class:`SeqStep`
            Contains the information related to the current simulation
            step.
        simODS: object of type :py:class:`Grid`
            Contains the dimension required to define a simulation domain.

    Returns:
        A file containing a auxiliary variable whith dimensions
        suitable for the current simulation step.
    """
    
    normDir = seq.direct # The direction normal to the current
                         # simulation plane
    currSlice = seq.level # The "level" of the current slice
    auxFile = seq.auxvarti_file # Here the name of the file containing
                                # the aux variable
    
    if not auxFile :
        return 1
    
    nx = simODS.nx
    ny = simODS.ny
    nz = simODS.ny
    
    dx = simODS.dx
    dy = simODS.dy
    dz = simODS.dz
    
    ox = simODS.ox
    oy = simODS.oy
    oz = simODS.oz
    
    data = vn.vtk2numpy(auxFile)[0]
    print(data.shape)
    
    f = open(auxFile, 'w')
    
    print("Aux file modified = ", auxFile)
    
    print("# vtk DataFile Version 3.0", file=f)
    print('Auxiliary variable "depth" Z', file=f)
    print("ASCII", file=f)
    print("DATASET STRUCTURED_POINTS", file=f)
    if normDir == 'x':
        print("DIMENSIONS 1 %d %d" % (ny, nz), file=f)
        print("ORIGIN %f %f %f" % (ox + (currSlice ) * dx, oy, oz), file=f)
        print("SPACING %f %f %f" % (dx, dy, dz), file=f)
        print("POINT_DATA %d" % (ny * nz), file=f)
        print("SCALARS facies float", file=f)
        print("LOOKUP_TABLE default", file=f)
                
        for k in range(nz):
            for j in range(ny):
                print(data[0 , j, k], file=f)        
    elif normDir == 'y':
        print("DIMENSIONS %d 1 %d" % (nx, nz), file=f)
        print("ORIGIN %f %f %f" % (ox, oy + (currSlice ) * dy, oz), file=f)
        print("SPACING %f %f %f" % (dx, dy, dz), file=f)
        print("POINT_DATA %d" % (nx * nz), file=f)
        print("SCALARS facies float", file=f)
        print("LOOKUP_TABLE default", file=f)
                
        for k in range(nz):
            for i in range(nx):
                print(data[i, 0 , k], file=f)        
    elif normDir == 'z':
        
        print("DIMENSIONS %d %d 1" % (nx, ny), file=f)
        print("ORIGIN %f %f %f" % (ox, oy, oz + (currSlice ) * dz), file=f)
        print("SPACING %f %f %f" % (dx, dy, dz), file=f)
        print("POINT_DATA %d" % (nx * ny), file=f)
        print("SCALARS facies float", file=f)
        print("LOOKUP_TABLE default", file=f)
        
        for j in range(ny):
            for i in range(nx):
                print(data[i, j, 0], file=f)
    
    f.close()


def print_sim_info(simODS, par_Xnorm, par_Ynorm, par_Znorm, nthreads):
    '''
    Print some information about the parameters of the simulation.

    Parameters:
        simODS: object of type :py:class:`geone.img.Img`
            Contains info about the simulation grid.
        par_Xnorm: object of type :py:class:`Param` or :py:class:`Param`.
            Information about the `Impala` parameters normal to the
            direction *x*.  (idem for *y* and *z*) and for `MPDS`.
    
            
    .. note:
        * Only one training image per direction.
        * Only one aux var per direction.
    ''' 

    print("\n    Dimension of the simulation domain:"
          " {0.nx:d} x {0.ny:d} x {0.nz:d}" \
          .format(simODS))

    for dire, par_dir in zip(('x','y','z'),(par_Xnorm,
                                            par_Ynorm, par_Znorm)):
        # For each direction, print the TI used
        print("    Plane normal to axis *", dire, "*:", sep='')

        # Check the MPS simulation engine in use
        # print(type(par_dir))
        if par_dir is not None:
            if isinstance(par_dir, gn.deesseinterface.DeesseInput):
                print('        MPS simulation engine: "DeeSse" with {0} threads'.format(nthreads))
        # 
        # If you want to use a different simulation engine, you can select
        # it with something like this:
        #
        # elif type(par_dir) == XXX_interface.Param:
        #     print('        MPS simulation engine: "XXX"')
        #
            else:
                print('        Warning: unknown MPS simulation engine')
        
            # NOTE: HERE ONLY ONE TI IS CONSIDERED.
            print('        TI file:               "', 
                  par_dir.TI[0].name, '"', sep='')
        else:
            # There is no training image along the selected direction
            print('        No TI defined ')
        # if (par_dir.tis[0].file_name is not None and
        #     type(par_dir) == impala_interface.Param):
        #     if par_dir.nb_aux_var > 0:
        #         print('        aux.var file:          "', 
        #               par_dir.aux_var_files[0], '"', sep='')
    
 
def create_lists( simODS, par_template, par_Xnorm, par_Ynorm, par_Znorm):
    """
    Create the lists which will be used by *Impala* for the simulation.

    This preliminary step is required and allows use *Impala* to
    compute the MPS list only once, before starting the simulation,
    and therefore save computing resources.

    Parameters:
        simODS: object of type :py:class:`Grid`
            Contains all the info concerning the simulation domain.
        par_template: impala_interface.Param
            All the information contained in the `template.in` file.
        par_Xnorm: object of type :py:class:`impala_interface.Param`
            Information about the *Impala* parameters normal to the
            direction *x*.
        par_Ynorm: object of type :py:class:`impala_interface.Param`
            Information about the *Impala* parameters normal to the
            direction *y*.
        par_Znorm: object of type :py:class:`impala_interface.Param`
            Information about the *Impala* parameters normal to the
            direction *z*.

    Returns: 
        A list (in binary format) for each multi-grid level and for
        each simulation direction is created in the current directory
        running *Impala*.
        
    """
    print('\n    *** Create the lists ***')

    seq = []

    # Check if there are training images along the given directions
    ti_x = check_ti_file(par_Xnorm)
    ti_y = check_ti_file(par_Ynorm)
    ti_z = check_ti_file(par_Znorm)
    if ti_x :
        seq.append(SeqStep( 'x', 0, par_Xnorm))
    if ti_y :
        seq.append(SeqStep( 'y', 0, par_Ynorm))
    if ti_z :
        seq.append(SeqStep( 'z', 0, par_Znorm))


    for seq_step in seq:
        simul_out = seq_step.create_list( simODS, par_template)




def create_seq(simODS, par_Xnorm, par_Ynorm, par_Znorm, nthreads=1, pseudo3D=0):
    '''
    Creates a simple simulation sequence.

    Parameters:
        simODS: grid definition
            Information about the simulation grid.
        par_Xnorm, par_Ynorm, par_Znorm: :py:class:`Param`
            Information related to the simulation directions.
            See the class documentation for details.
        pseudo3D: integer (default=0)
            Set or not the "pseudo3D" simulation option when this
            value is >0.
            WARNING: FOR THE MOMENT THIS IS NOT IMPLEMENTED HERE.
    Returns:
       * A sequence of :py:class:`SeqStep` objects, with all the
         information needed to run each simulation step.

    .. note::

        This function provide the **basic** definition for a
        simulation sequence when the simulation domain is quite
        simple, like for example when it has a "box" shape (the sizes
        of the simulation grid along the directions *x*, *y* and *z*
        are comparable).  If the dimensions along the different axes
        of your simulation domain are not comparable, you can for
        example explicitly add some customized simulation steps to
        improve the quality of the results. There is one example of
        this in the examples directory.

    '''
    
    # Create "matrioska" simulation intervals only along directions
    # normal to a provided TI.
    ti_x = check_ti_file(par_Xnorm, "yz")
    ti_y = check_ti_file(par_Ynorm, "xz")
    ti_z = check_ti_file(par_Znorm, "xy")
    if ti_x :
        sim_path_Xnorm = matrioska_interval(simODS.nx)
    if ti_y :
        sim_path_Ynorm = matrioska_interval(simODS.ny)
    if ti_z :
        sim_path_Znorm = matrioska_interval(simODS.nz)

    # Create the simulation sequence (a list of objects `SeqStep`).
    # For information about the `SeqStep` object see the documentation.
    seq = []
    for i in range(max(simODS.nx, simODS.ny, simODS.nz)):
        if ( ti_x and i < simODS.nx ):
            try:
                seq.append(SeqStep( 'x', sim_path_Xnorm[i], par_Xnorm, nthreads))
            except IndexError:
                print("\n"
                      "    ERROR in 'create_seq'.\n"
                      "        Check the size of the simulation grid.\n"
                      "        The domain should be 3D.\n"
                      "        Therefore, size along *x* should be >1, "
                      "but '{0}' was provided.\n".format(simODS.nx))
                     
        if ( ti_y and i < simODS.ny ):
            try:
                seq.append(SeqStep( 'y', sim_path_Ynorm[i], par_Ynorm, nthreads))
            except IndexError:
                print("\n"
                      "    ERROR in 'create_seq'.\n"
                      "        Check the size of the simulation grid.\n"
                      "        The domain should be 3D.\n"
                      "        Therefore, size along *y* should be >1, "
                      "but '{0}' was provided.\n".format(simODS.ny))
        if ( ti_z and i < simODS.nz ):
            try:
                seq.append( SeqStep( 'z', sim_path_Znorm[i], par_Znorm, nthreads))
            except IndexError:
                print("\n"
                      "    ERROR in 'create_seq'.\n"
                      "        Check the size of the simulation grid.\n"
                      "        The domain should be 3D.\n"
                      "        Therefore, size along *z* should be >1, "
                      "but '{0}' was provided.\n".format(simODS.nz))
            
    return seq


def check_ti_file(par, plane=None):
    '''
    Check if the parameter info for the simulation along a given 
    direction have a TI assigned.

    Parameters:
        par: geone.deesseinterface.DeesseInput 
            Contains the info about the simulation parameters, including the
            TI.

    Returns:
        True if a TI file name is defined, false if None.

    .. note:
        Check only the 1st TI.
    '''
    if par is not None:
        # NOTE: WORKS ONLY  FOR 1 TI
        if par.TI[0].name:
            return True
        else:
            return False
    else:
        # LATER PUT THIS PROPERLY INTO A LOG FILE.
        print("    WARNING: TI is missing along plane *{0}*.".format(plane))
 

# def add_gslib_pointdata(data_files, hard_data, simODS):
#     '''
#     Add the hard data contained in a number of GSLIB point data files.

#     Parameters:
#         files: list of strings
#             A number of files containing the data. Note that this must
#             be in a list format even if only one file is
#             considered. For example, you should always use a syntax
#             like `["file1.gslib"]` even if you provide only one file.
#         hard_data: 3D numpy array
#             Where all the simulated points are stored.
#         simODS: object of type :py:class:`Grid`
#             Contains all the info concerning the simulation domain.

#     Returns:
#         Update the content of the array `hard_data`.

#     ''' 

#     #
#     # NOTE: I SHOULD CHECK IF THE COORDINATES OF THE POINTS
#     #       CORRESPONDS TO THE TRUE COORDINATES INSIDE THE GRID SYSTEM
#     #       OF IMPALA!
#     #

#     for data_file in data_files:
#         data_all, header = gslibnumpy.gslib_points2numpy(data_file)
#         x = data_all['x']
#         y = data_all['y']
#         z = data_all['z']
#         data = data_all['data']

#         nb_data_added = 0
#         for i in range(len(x)):
#            ix = int(round( ( x[i] - simODS.ox) / simODS.dx))
#            if (ix < 0 or ix > simODS.nx-1):
#                print('    Warning: detected conditioning data '
#                      'outside the ', end=' ')
#                print('simulation domain. Datum ignored.')
#                print('        function: "add_gslib_pointdata"')
#                print('        data file: "', data_file,'"')
#                simODS.print_intervals('x')
#                continue
#            iy = int(round( ( y[i] - simODS.oy) / simODS.dy))
#            if (iy < 0 or iy > simODS.ny-1):
#                print('    Warning: detected conditioning data outside'
#                      ' the ', end=' ')
#                print('simulation domain. Datum ignored.')
#                print('        function: "add_gslib_pointdata"')
#                print('        data file: "', data_file,'"')
#                simODS.print_intervals('y')
#                continue
#            iz = int(round( ( z[i] - simODS.oz) / simODS.dz))
#            if (iz < 0 or iz > simODS.nz-1):
#                print('    Warning: detected conditioning data'
#                      ' outside the ', end=' ')
#                print('simulation domain. Datum ignored.')
#                print('        function: "add_gslib_pointdata"')
#                print('        data file: "', data_file,'"')
#                simODS.print_intervals('z')
#                continue

#            #
#            # NOTE: IF I COULD ACCESS THE ALLOWED FACIES OF IMPALA,
#            #       HERE I COULD CHECK IF THE DATA IS INSIDE THE RANGE
#            #       OF AVAILABLE FACIES.
#            #
#            hard_data[ix,iy,iz] = data[i]
#            nb_data_added = nb_data_added + 1

#         print('    Number of hard data added:', nb_data_added)
           
           



def file_name_sec(file_name, axe):
    """
    Return the name of a file which will be used to store a section
    of the input parameter `file_name`.

    Parameters:
        file_name: string
            The name of the input file which should be sectioned
        axe: string in ('x','y','z')
            The axis perpendicular to the section.
    
    .. note::
        It is supposed that the sectioned file lies in the current 
        directory. Therefore, the path in the name of the input
        file name will be dropped.
    """
    import os
    import ext
    add_ext = {'x':'_yz', 'y':'_xz', 'z':'_xy'}
    return os.path.basename(file_name).rsplit(".",1)+add_ext[axe]+ext.VTK
    
               

               
#
# Here all the testing functions...
#
if __name__ == '__main__':

    import numpy as np    

    #
    # test the function "matrioska_interval"
    #
    print('\n    *** testing the function "matrioska_interval" ***\n')

    points_nb = -5
    print('    points_nb:', points_nb)
    print('    matrioska_interval:', matrioska_interval(points_nb),'\n')

    points_nb = 0
    print('    points_nb:', points_nb)
    print('    matrioska_interval:', matrioska_interval(points_nb),'\n')

    points_nb = 1
    print('    points_nb:', points_nb)
    print('    matrioska_interval:', matrioska_interval(points_nb),'\n')

    points_nb = 7
    print('    points_nb:', points_nb)
    print('    matrioska_interval:', matrioska_interval(points_nb),'\n')

    points_nb = 8
    print('    points_nb:', points_nb)
    print('    matrioska_interval:', matrioska_interval(points_nb),'\n')
    
    
               
    print('    get_mid((1,2)):', get_mid((1,2)) )
    print('    get_mid((7,45)):', get_mid((7,45)) )
    print('    get_mid((1,20)):', get_mid((1,20)) )
    print('    get_mid((6,11)):', get_mid((6,11)) )

    print("input:",  [(10,12),(15, 35)])
    print("out:", split_segms( [(10,12),(15, 35)]))
    #print( fluff([3], [(10,11),(15, 16)]))

