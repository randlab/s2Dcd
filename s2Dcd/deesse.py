"""
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

    A collection of classes and functions to interact with the
    Multiple Point Direct Sampling code (DeeSse).

:File name:
    ``deesse.py``

:Version:
    0.7 , 2020-01-20 :
        * Some clean up before the upload to Github.
    0.6 , 2017-12-15 :
        * Adapted to the 2017 vesion of the code (includes "Pyramids" option)
        * Partially takes the comments for the .in file from an
          external .json file.
    0.5 , 2015-10-15 :
        * Extending to deal with multiple variables simulation.
    0.4 , 2013-06-17 :
        * Some small modifications to make the setting of the default values
          more coherent.
    0.3 , 2013-03-13 :
        * Moved to Python3 using `2to3`.
    0.2 , 2012-07-30 :
        * Moved the function `skip_ccomments` to the module `utili.py`.
        * Last version before the movement to Python3.
    0.1 , 2012-03-30 :
        First version.

:Authors:
    Alessandro Comunian

.. note::
   - A `deesse_file_in_text.json` is required and contains all comment text
     that put in the deesse `.in` file. To create easily a `.json` file that
     respects the intendation and the new lines of the original `.in` file,
     you can use something like::

         gawk '$1=$0' ORS='\\n' file.in

     to replace all the new lines by \\n, at the same time keeping the
     original indentation of the file.

.. warning::

    - The coverage and the support of the simulation options provided
      by DS simulation engine is only partial.


.. todo::
    - Implement the usage of rotations.

"""
import s2Dcd.grid as grid
import s2Dcd.utili as utili
import json
import os


# #########################################################################
# Some default values for some parameters
# #########################################################################
#
# IN THE FUTURE ALL THESE DEFAULT PARAMETERS SHOULD GO INTO AN
# *.ini or *.json file
#
# Size of the search radis along *x* [number of nodes]
SEARCH_R_X = 120.0
# Size of the search radis along *y* [number of nodes]
SEARCH_R_Y = 120.0
# Size of the search radis along *z* [number of nodes]
SEARCH_R_Z = 0.0
# Anisotropy ratio of the search template along *x* [-]
SEARCH_ANIS_X = 1.0
# Anisotropy ratio of the search template along *y* [-]
SEARCH_ANIS_Y = 1.0
# Anisotropy ratio of the search template along *z* [-]
SEARCH_ANIS_Z = 1.0
# Azimuth orientation of the search template [degrees]
SEARCH_ANG_AZ = 0.0
# Dip orientation of the search template [degrees]
SEARCH_ANG_DP = 0.0
# Plunge orientation of the search template [degrees]
SEARCH_ANG_PL = 0.0
# Search "power", for computing weigts according to distance
SEARCH_POW = 0.0
# Maximum number of neighborhood nodes
MAX_NB_NEIGH = 20
# Maximal desity of neighborhing nodes in the search neighborhood
MAX_DENS_NEIGH = 1.0
# Relative distance flag (can be 0 or 1)
REL_DIST_FLAG = 0
# Distance type
DIST_TYPE = 0
# Weight factor for conditioning data
WEIGHT = 1.0
# Distance threshold
DIST_THRE = 0.01
# The default value for the % of the TI to scan
MAX_SCAN = 0.3
# Default variable name
NAME = "facies"
# Default name for the TI file
FILE_NAME = "ti.gslib"
# Default option for the output
OUT_SET = 'OUTPUT_SIM_ALL_IN_ONE_FILE'
# Default name for the output file
OUT_FILE = 'test_simul.gslib'
# Default setting for the print of the report
OUT_REPORT = True
# Default format to print an output variable
FMT = "%10.5E"

#
# These are some settings to deal with problems related to the
# speed of connection with the licence server (I suppose).
#
# Number of attepts to made when there is an error licence
NB_LIC_WAIT = 20
# Duration of the pause when waiting for the licence
LIC_WAIT_TIME = 0.2
# This is the error code of the current version of the deesse in case
# of lincese problems
LIC_ERR_CODE = 40704


# Get the directory containing the module
moddir = os.path.dirname(globals()["__file__"])

# Read the file containing the text to be put in the
# comments of the in file
json_in_text = json.load(
    open(os.path.join(moddir, "deesse_file_in_text.json")))


class VarInfo(object):
    '''
    A class to support the class ParamMPDS. It tries to regroup all
    the information about a variable which have to be simulated.
    '''
    def __init__(self, name=NAME, out_flag=True, fmt=FMT,
                 search_r_x=SEARCH_R_X,
                 search_r_y=SEARCH_R_Y,
                 search_r_z=SEARCH_R_Z,
                 search_anis_x=SEARCH_ANIS_X,
                 search_anis_y=SEARCH_ANIS_Y,
                 search_anis_z=SEARCH_ANIS_Z,
                 search_ang_az=SEARCH_ANG_AZ,
                 search_ang_dp=SEARCH_ANG_DP,
                 search_ang_pl=SEARCH_ANG_PL,
                 search_pow=SEARCH_POW,
                 max_nb_neigh=MAX_NB_NEIGH,
                 max_dens_neigh=MAX_DENS_NEIGH,
                 rel_dist_flag=REL_DIST_FLAG,
                 dist_type=DIST_TYPE,
                 weight=WEIGHT,
                 dist_thre=DIST_THRE,
                 prob_constr=0,
                 block_data=0
                 ):
        '''
        Parameters:
            name: string, optional
                Name of the variable
            prob_constr: int, optional
                This is used for the probability constrains.
                For the moment this is not implemented in the code, and
                therefore only a default value of 0 is allowed
            block_data: int, optional
                This is required for the usage of block data. For the moment
                only the default value of 0 is allowed.
            out_flag: bool, optional
                A flag to determine if the variable
            ...

        .. note:
            See the MPDS documentation for more details.
        '''

        self.name = name
        self.out_flag = out_flag
        self.fmt = fmt
        self.search_r_x = search_r_x
        self.search_r_y = search_r_y
        self.search_r_z = search_r_z
        self.search_anis_x = search_anis_x
        self.search_anis_y = search_anis_y
        self.search_anis_z = search_anis_z
        self.search_ang_az = search_ang_az
        self.search_ang_dp = search_ang_dp
        self.search_ang_pl = search_ang_pl
        self.search_pow = search_pow
        self.max_nb_neigh = max_nb_neigh
        self.max_dens_neigh = max_dens_neigh
        self.rel_dist_flag = rel_dist_flag
        self.dist_type = dist_type
        self.weight = weight
        self.dist_thre = dist_thre
        self.prob_constr = prob_constr
        self.block_data = block_data
    pass


class TiInfo(object):
    '''
    A class to support the class ParamMPDS. Contains all the information
    related to one training image.

    Main attributes:
        file_name: string
            Name of the file containing the training image
        max_scan: float
            Fraction of the TI to scan

    '''

    def __init__(self, file_name=FILE_NAME, max_scan=MAX_SCAN):
        """
        Parameters
            file_name: string
                Name of the file containing the training image
            max_scan: float
                Fraction of the TI to scan
        """
        self.file_name = file_name
        self.max_scan = max_scan

    pass

    def __str__(self):
        """
        Print a TiInfo object
        """
        out = ("*** TiInfo object ***\n"
               "    file_name: {0.file_name}\n"
               "    max_scan:  {0.max_scan}\n".format(
                self))
        return out


class Param(object):
    '''
    A class containing all the parameters contained in the parameters file.
    '''

    def __init__(self, file_name=None,
                 gridin=grid.Grid(),
                 var_nb=1,
                 vars=[VarInfo()],
                 out_set = OUT_SET,
                 out_file = OUT_FILE,
                 out_report = OUT_REPORT,
                 out_report_file = 'test_report.txt',
                 ti_nb = 1,
#                 tis = [TiInfo(file_name=FILE_NAME, max_scan=MAX_SCAN)], 
#                 tis = [], 
                 data_img_nb = 0, 
                 data_img_files = [], 
                 data_pointset_nb = 0, 
                 data_pointset_files = [], 
                 mask = False, 
                 homot_usage = 0, 
                 rot_usage = 0, 
                 cond_data_tol = 0.05, 
                 norm_type = 'NORMALIZING_LINEAR', 
                 sim_type = 'SIM_ONE_BY_ONE', 
                 sim_path = 'PATH_RANDOM', 
                 tol = 0.0, 
                 post_proc_path_max = 1, 
                 post_proc_par = 'POST_PROCESSING_PARAMETERS_DEFAULT',
                 pyramids = 0,
                 seed = 444, 
                 seed_inc = 1, 
                 real_nb = 1 
                 ):
        '''
        The initializer for the class `ParamMPDS`.
        
        The class ParamMPDS contains all the parameters required to
        run a MPDS simulation.
        
        Parameters:
            file_name: string, optional
                The name of the template parameter input file. If it
                is not provided, then the default values of all the
                parameters are used.
            pyramids: int, optional
                Parameter required to turn on the simulation using a pyramyd
                approach. For the moment only the default value 0 is allowed.
        
        .. note::
            If a value for `file_name` is provided, then all the other
            parameters are ignored.

        .. todo::
            - Add some global variables where the default values are
              defined.
        '''

        if file_name:
            self.file_name = file_name
            try:
                par_file = open(self.file_name, 'r')
            except IOError:
                print('   Error reading file "%s"'%(file_name))
                
            # Extract a list containing all the parameters
            par_list = utili.skip_ccomments(par_file)
            
            # Initialize the pointer to the current line of the list
            cur = 0
            
            # Create a grid object
            self.grid = grid.Grid()
            
            # Read the simulation grid
            self.grid.nx = int(par_list[cur][0])
            self.grid.ny = int(par_list[cur][1])
            self.grid.nz = int(par_list[cur][2])
            cur = cur + 1 
            
            # Read the spacing
            self.grid.dx = float(par_list[cur][0])
            self.grid.dy = float(par_list[cur][1])
            self.grid.dz = float(par_list[cur][2])
            cur = cur + 1 
            
            # Read the origin
            self.grid.ox = float(par_list[cur][0])
            self.grid.oy = float(par_list[cur][1])
            self.grid.oz = float(par_list[cur][2])
            cur = cur + 1
            
            # Read the number of simulation variables
            self.var_nb = int(par_list[cur][0])
            cur = cur + 1
            
            # Create a list of variable(s) objects
            self.vars = []
            for i in range(self.var_nb):
                # VERIFY HERE WHAT IT HAPPENS WHEN THE FORMAT IS NOT
                # DEFINED, FOR EXAMPLE IF A VARIABLE IS NOT TO BE
                # PRINTED.
                self.vars.append(
                    VarInfo(name=par_list[cur+i][0],
                            out_flag=bool(int(par_list[cur+i][1])),
                            fmt=par_list[cur+i][2]
                            ))
            cur = cur + self.var_nb
            
            # Read the output setting for the simulation
            self.out_set = par_list[cur][0]
            cur = cur + 1

            # Read the additional information about the output settings
            if self.out_set == 'OUTPUT_SIM_NO_FILE':
                # No output file required
                self.out_file = None
            elif self.out_set == 'OUTPUT_SIM_ALL_IN_ONE_FILE':
                # All the simulations in one file
                self.out_file = par_list[cur][0]
                cur = cur + 1
            elif self.out_set == 'OUTPUT_SIM_ONE_FILE_PER_VARIABLE':
                # One output file per flagged variable
                for var in self.vars:
                    if var.out_flag:
                        var.out_file = par_list[cur][0]
                        cur = cur + 1
            elif self.out_set == 'OUTPUT_SIM_ONE_FILE_PER_REALIZATION':
                # One file for each realization
                self.out_prefix = par_list[cur][0]
                cur = cur + 1
            else:
                print('    Error, wrong output setting value'
                      ' (={0.out_set}),'.format(self))
                return -1

            # Read the flags for the output report
            self.out_report = bool(int(par_list[cur][0]))
            cur = cur + 1
            if self.out_report:
                # Read the file containing the report
                self.out_report_file = par_list[cur][0]
                cur = cur + 1

            # Read the information about the training image
            self.ti_nb = int(par_list[cur][0])
            cur = cur + 1

            # Read the training image(s) files
            self.tis = []
            for i in range(self.ti_nb):
                self.tis.append(TiInfo())
                self.tis[i].file_name = par_list[cur+i][0]
            cur = cur + self.ti_nb

            # Read the data image file(s) related for the Simulation grid
            self.data_img_nb = int(par_list[cur][0])
            cur = cur + 1

            # Read the data image(s) files (An empty list in created
            # in any case to ease the print procedure...)
            self.data_img_files = []
            for i in range(self.data_img_nb):
                self.data_img_files.append(par_list[cur+i][0])
            cur = cur + self.data_img_nb
            
            # Read the data point set file(s) related for the
            # Simulation grid
            self.data_pointset_nb = int(par_list[cur][0])
            cur = cur + 1

            # Read the data point set files (An empty list in created
            # in any case to ease the print procedure...)
            self.data_pointset_files = []
            for i in range(self.data_pointset_nb):
                self.data_pointset_files.append(par_list[cur+i][0])
            cur = cur + self.data_pointset_nb

            # Read the flag concerning the mask
            self.mask = bool(int(par_list[cur][0]))
            cur = cur + 1

            # Read the file containing the mask
            if self.mask:
                self.mask_file = par_list[cur][0]
                cur = cur + 1
         
            # Read all the info concerning homothety
            self.homot_usage = int(par_list[cur][0])
            cur = cur + 1

            # Read the optional parameters about homothety
            if self.homot_usage == 0:
                # No homothety
                pass
            elif self.homot_usage == 1:
                # Homothety without tolerance
                self.homot_img_flag_x = bool(int(par_list[cur][0]))
                if self.homot_img_flag_x:
                    # Read the name of the image
                    self.homot_img_file_x = par_list[cur][1]
                else:
                    # Read the value
                    self.homot_val_x = float(par_list[cur][1])
                    cur = cur + 1
                self.homot_img_flag_y = bool(int(par_list[cur][0]))
                if self.homot_img_flag_y:
                    # Read the name of the image
                    self.homot_img_file_y = par_list[cur][1]
                else:
                    # Read the value
                    self.homot_val_y = float(par_list[cur][1])
                cur = cur + 1
                self.homot_img_flag_z = bool(int(par_list[cur][0]))
                if self.homot_img_flag_z:
                    # Read the name of the image
                    self.homot_img_file_z = par_list[cur][1]
                else:
                    # Read the value
                    self.homot_val_z = float(par_list[cur][1])            
                cur = cur + 1
            elif self.homot_usage == 2:
                # Homothety with tolerance
                self.homot_img_flag_x = bool(int(par_list[cur][0]))
                if self.homot_img_flag_x:
                    # Read the name of the image
                    self.homot_img_file_x = par_list[cur][1]
                else:
                    # Read the value
                    self.homot_lowval_x = float(par_list[cur][1])
                    self.homot_uppval_x = float(par_list[cur][2])
                    cur = cur + 1
                self.homot_img_flag_y = bool(int(par_list[cur][0]))
                if self.homot_img_flag_y:
                    # Read the name of the image
                    self.homot_img_file_y = par_list[cur][1]
                else:
                    # Read the value
                    self.homot_lowval_y = float(par_list[cur][1])
                    self.homot_uppval_y = float(par_list[cur][2])
                cur = cur + 1
                self.homot_img_flag_z = bool(int(par_list[cur][0]))
                if self.homot_img_flag_z:
                    # Read the name of the image
                    self.homot_img_file_z = par_list[cur][1]
                else:
                    # Read the value
                    self.homot_lowval_z = float(par_list[cur][1])          
                    self.homot_uppval_z = float(par_list[cur][2])           
                cur = cur + 1
            else:
                print('    Error, wrong value for the homothety')
                return -1

            # Read info about rotation
            self.rot_usage = int(par_list[cur][0])
            cur = cur + 1

            #
            # TO BE IMPLEMENTED
            #
            if self.rot_usage != 0:
                print('    Error, ROTATION options not implemented yet...')
                return -1

            # Info about consistency of conditioning data
            self.cond_data_tol = float(par_list[cur][0])
            cur = cur + 1

            # Info about the normalization type
            self.norm_type = par_list[cur][0]
            cur = cur + 1

            # Info about seach neighborhood parameters
            for var in self.vars:
                var.search_r_x = float(par_list[cur][0])
                var.search_r_y = float(par_list[cur][1])
                var.search_r_z = float(par_list[cur][2])
                cur = cur + 1
                var.search_anis_x = float(par_list[cur][0])
                var.search_anis_y = float(par_list[cur][1])
                var.search_anis_z = float(par_list[cur][2])
                cur = cur + 1
                var.search_ang_az = float(par_list[cur][0])
                var.search_ang_dp = float(par_list[cur][1])
                var.search_ang_pl = float(par_list[cur][2])
                cur = cur + 1
                var.search_pow = float(par_list[cur][0])
                cur = cur + 1

            # Read the max number of neighboring nodes
            for i, var in enumerate(self.vars):
                var.max_nb_neigh = int(par_list[cur][i])
            cur = cur + 1
                
            # Read the max density of neighboring nodes
            for i, var in enumerate(self.vars):
                var.max_dens_neigh = float(par_list[cur][i])
            cur = cur + 1
   
            # Read the relative distance flag
            for i, var in enumerate(self.vars):
                var.rel_dist_flag = bool(int(par_list[cur][i]))
            cur = cur + 1

            # Read the distance type
            for i, var in enumerate(self.vars):
                var.dist_type = int(par_list[cur][i])
            cur = cur + 1

            # Read the weight factor for conditioning data
            for i, var in enumerate(self.vars):
                var.weight = float(par_list[cur][i])
            cur = cur + 1
                                   
            # Read the simulation type
            self.sim_type = par_list[cur][0]
            cur = cur + 1

            # Read the path type
            self.sim_path = par_list[cur][0]
            cur = cur + 1

            # Read the distance threshold
            for i, var in enumerate(self.vars):
                var.dist_thre = float(par_list[cur][i])
            cur = cur + 1

            # Read the probability constraints
            for var in self.vars:
                var.prob_constr = int(par_list[cur][0])
  #              print("probl constraint", var.prob_constr)
                if var.prob_constr == 1:
                    # WARNING: THIS IS A QUICK AND DIRTY TRICK.
                    # FUTURE IMPLEMENTATIONS SHOULD HANDLE THIS
                    # IN A MORE FLEXIBLE WAY!
                    cur = cur + 1 
                    var.prob_ln1 = " ".join(par_list[cur])
                    cur = cur + 1 
                    var.prob_ln2 = " ".join(par_list[cur])
                    cur = cur + 1 
                    var.prob_ln3 = " ".join(par_list[cur])
                    cur = cur + 1 
                    var.prob_ln4 = " ".join(par_list[cur])
                    cur = cur + 1 
                    var.prob_ln5 = " ".join(par_list[cur])
 #                   print(var.prob_ln1)
 #                   print(var.prob_ln2)
 #                   print(var.prob_ln3)
 #                   print(var.prob_ln4)
 #                   print(var.prob_ln5)
#                    cur = cur + 1
                elif var.prob_constr == 2:
                    print("    WARNING: Probability constraints =2"
                          " not implemented!")
                    print("        Only the default value '0' or '1'"
                          "are allowed for the moment.")
                    print("        Double check if your "
                          "parameter files are")
                    print("        correctly created.")
                cur = cur + 1

            # Read block data information
            for var in self.vars:
                var.block_data = int(par_list[cur][0])
#                print(var.block_d\ata)
                if var.block_data != 0:
                    print("    WARNING: Block data option"
                          " not implemented!")
                    print("        Only the default value '0' "
                          "is allowed for the moment.")
                    print("        (read value=",var.block_data,")")
                    print("        Double check if your "
                          "parameter files are")
                    print("        correctly created.")
                cur = cur + 1
            

            # Read the maximal scan fraction
            for i, ti in enumerate(self.tis):
                ti.max_scan = float(par_list[cur][i])
            cur = cur + 1

            # Read the tolerance
            self.tol = float(par_list[cur][0])
            cur = cur + 1

            # Read the post-processing options
            self.post_proc_path_max = int(par_list[cur][0])
            cur = cur + 1

            # Read the post processing key word
            if self.post_proc_path_max > 0:
                self.post_proc_par = par_list[cur][0]
                cur = cur + 1
                if (self.post_proc_par == 
                    "POST_PROCESSING_PARAMETERS_MANUAL"):
                    print("    ERROR! The option 'POST_PROCESSING_MANUAL'"
                          " is") 
                    print("        not implemented here!")

            # Read the pyramids option
            self.pyramids = int(par_list[cur][0])
            if (self.pyramids != 0):
                print("    ERROR: Pyramid option={0}".format(self.pyramids))
                print("        Only the default value 0 is")
                print("        implemented for the moment.")
            cur = cur + 1

                    
            # Read the seed
            self.seed = int(par_list[cur][0])
            cur = cur + 1

            # Read the seed increment
            self.seed_inc = int(par_list[cur][0])
            cur = cur + 1

            # Read the number of realizations
            self.real_nb = int(par_list[cur][0])
            cur = cur + 1

            # Read the end of file keywork to check if everithing is fine
            if par_list[cur][0] != "END":
                print('    Error in the end of file keyword!')
        else:
            # A parameter file template was not provided, use the
            # default values.

            # Default values for the simulation grid
            self.grid = gridin
      
            # Default number of variables
            self.var_nb = var_nb
            
            # Default value for the variable names and the output flag
            self.vars = vars

            # Default output settings
            self.out_set = out_set
            self.out_file = out_file

            # Default values for the output report
            self.out_report = out_report
            self.out_report_file = out_report_file

            # Default number for the training images
            self.ti_nb = ti_nb
            self.tis = []
            self.tis.append(TiInfo(file_name=FILE_NAME,
                                   max_scan=MAX_SCAN))

            # Default value for the number of data images
            self.data_img_nb = data_img_nb
            self.data_img_files = data_img_files
            self.data_pointset_nb = data_pointset_nb
            self.data_pointset_files = data_pointset_files
            self.mask = mask
            self.homot_usage = homot_usage
            self.rot_usage = rot_usage
            self.cond_data_tol = cond_data_tol
            self.norm_type = norm_type
            self.sim_type = sim_type
            self.sim_path = sim_path
            self.tol = tol
            self.post_proc_path_max = post_proc_path_max
            self.post_proc_par = post_proc_par
            self.pyramids = pyramids
            self.seed = seed
            self.seed_inc = seed_inc
            self.real_nb = real_nb

        
    def __str__(self):
        '''
        A print method for the class ParamMPDS

        .. note::
            * the names of the variables in the section "SIMULATION
              VARIABLES" are not printed according to their original
              order, but ordered alphabetically.
        '''

        out = ('{0}{1}{2}{3}{4}{5}{6}{7}{8}{9}{10}{11}'
               '{12}{13}{14}{15}{16}{17}{18}{19}{20}{21}'
               '{22}{23}{24}{25}{26}{27}{28}{29}{30}'
            ).format(self.print_sim_grid(),
                     self.print_sim_var(),
                     self.print_out_set(),
                     self.print_out_report(),
                     self.print_ti(),
                     self.print_difile(),
                     self.print_dpset(),
                     self.print_maskimage(),
                     self.print_homo(),
                     self.print_rot(),
                     self.print_consi(), 
                     self.print_norm(),
                     self.print_neigh_intro(),
                     self.print_neighs(),
                     self.print_max_nod(),
                     self.print_max_dens(),
                     self.print_reldist(),
                     self.print_disttype(),
                     self.print_weight(),
                     self.print_simpath(),
                     self.print_distthr(),
                     self.print_proconst_intro(),
                     self.print_proconst(),
                     self.print_blockdata_intro(), 
                     self.print_blockdata(),
                     self.print_maxscan(),
                     self.print_tol(),
                     self.print_post(),
                     self.print_pyramids(),
                     self.print_seed(),
                     self.print_realnb()
            )
        return out
    
    def print_sim_grid(self):
        '''
        Print the simulation grid info
        '''
        out = ('/* SIMULATION GRID (SG) */\n'
               '{0.nx:d} {0.ny:d} {0.nz:d} // size in each direction\n'
               '{0.dx} {0.dy} {0.dz} // spacing in each direction\n'
               '{0.ox} {0.oy} {0.oz} // origin\n\n').format(
            self.grid)
        
        return out


    def print_sim_var(self):
        '''
        Print the simulation variables info
        '''

        out = ('/* SIMULATION VARIABLES */\n'
               '/* Number of simulation variable(s), '
               'and for each variable:\n'
               '   variable name and output flag (0 / 1). */\n'
               '{0:d}\n'
               '{1:s}\n\n').format(
            self.var_nb,
            "\n".join([var.name + "    " + str(int(var.out_flag)) + "   " +
                       var.fmt 
                       for var in self.vars]))
        
        return out

    def print_out_set(self):
        '''
        Print the output settings
        '''
        files_string = {}
        files_string["OUTPUT_SIM_NO_FILE"] = ""
        files_string["OUTPUT_SIM_ALL_IN_ONE_FILE"] = (
            self.out_file 
            if self.out_set == "OUTPUT_SIM_ALL_IN_ONE_FILE"
            else ""
            )
        files_string["OUTPUT_SIM_ONE_FILE_PER_VARIABLE"] = (
            "\n".join([var.out_file for var in self.vars])   
            if self.out_set == "OUTPUT_SIM_ONE_FILE_PER_VARIABLE"
            else ""
            )
        files_string["OUTPUT_SIM_ONE_FILE_PER_REALIZATION"] = (
            self.out_prefix
            if self.out_set == "OUTPUT_SIM_ONE_FILE_PER_REALIZATION"
            else ""
            )
        
        out = (
            '/* OUTPUT SETTINGS FOR SIMULATION */\n'
            '/* Key word and required name(s) or prefix, '
            'for output of the realizations:\n'
            '      - OUTPUT_SIM_NO_FILE:\n'
            '           no file in output,\n'
            '      - OUTPUT_SIM_ALL_IN_ONE_FILE:\n'
            '           one file in output,\n'
            '           requires one file name\n'
            '      - OUTPUT_SIM_ONE_FILE_PER_VARIABLE:\n'
            '           one file per variable in output '
            '(flagged as 1 above),\n'
            '           requires as many file name(s) as variable(s)'
            ' flagged as 1 above\n'
            '      - OUTPUT_SIM_ONE_FILE_PER_REALIZATION:\n'
            '           one file per realization,\n'
            '           requires one prefix (for file name) */\n'
            '{0}\n{1}\n\n'
            ).format(
            self.out_set,
            files_string[self.out_set]
            ) 
        
        return out

    def print_out_report(self):
        '''
        Print the output report
        '''

        out = ('/* OUTPUT REPORT */\n'
               '/* Flag (0 / 1), and if 1, output report file. */\n'
               '{0:d}\n'
               '{1:s}\n\n'
            ).format(
            self.out_report,
            self.out_report_file,
            )

        return out



    def print_ti(self):
        '''
        Print the TI info
        '''
        out = (
            '/* TRAINING IMAGE */\n'
            '/* Number of training image(s) (nTI >= 1),'
            ' followed by nTI file(s)\n'
            '   (a file can be replaced by the string "_DATA_" '
            'which means that the\n'
            '   simulation grid itself is taken as training image),'
            'and\n'
            '   if nTI > 1, one pdf image file (for training images,'
            ' nTI variables). */\n'
            '{0:d}\n'
            '{1:s}\n\n'
            ).format(
            self.ti_nb,
            "\n".join([ti.file_name for ti in self.tis])
            )

        return out

    def print_difile(self):
        """
        Print the content of the data image file
        """

        out = ('/* DATA IMAGE FILE FOR SG */\n'
               '/* Number of image file(s) (n >= 0), followed by '
               'n file(s). */\n'
               '{0:d}\n'
               '{1}\n\n'
               ).format(
            self.data_img_nb,
            "\n".join(self.data_img_files) if self.data_img_files else ""
            )
        return out

    def print_dpset(self):
        """
        Print info related to the data point set
        """
        
        out = ('/* DATA POINT SET FILE FOR SG */\n'
               '/* Number of point set file(s) (n >= 0), followed by '
               'n file(s). */\n'
               '{0:d}\n'
               '{1}\n\n').format(
            self.data_pointset_nb,
            "\n".join(self.data_pointset_files) if
            self.data_pointset_files else "",
            )
        return out

    def print_maskimage(self):
        """
        Print infor related to the mask image
        """
        out = ('/* MASK IMAGE */\n'
               '/* Flag (0: mask not used / 1: mask used) and if 1,'
               ' mask image file\n'
               '   (this image contains one variable on the simulation '
               'grid: flag (0 / 1)\n'
               '   for each node of the simulation grid that indicates '
               'if the variable(s)\n'
               '   will be simulated at the corresponding node (flag 1) '
               'or not (flag 0). */\n'
               '{0:d}\n'
               '{1}\n\n').format(
            int(self.mask),
            self.mask_file if self.mask else "",
            )
        return out

    def print_homo(self):
        """
        Print info about the homothety
        """
        out = ('/* HOMOTHETY */\n'
               '/* 1. Homothety usage, integer (homothetyUsage):\n'
               '        - 0: no homothety\n'
               '        - 1: homothety without tolerance\n'
               '        - 2: homothety with tolerance\n'
               '   2a. If homothetyUsage == 1,\n'
               '          then for homothety ratio in each direction,\n'
               '          first for x, then for y, and then for z-axis '
               'direction:\n'
               '             - Flag (0 / 1) indicating if given in an '
               'image file,\n'
               '               followed by\n'
               '                  - one value (real) if flag is 0\n'
               '                  - name of the image file (one variable) '
               'if flag is 1\n'
               '   2b. If homothetyUsage == 2,\n'
               '          then for homothety ratio in each direction,\n'
               '          first for x, then for y, and then for z-axis '
               'direction:\n'
               '             - Flag (0 / 1) indicating if given in an '
               'image file,\n'
               '               followed by\n'
               '                  - two values (lower and upper bounds) '
               '(real) if flag is 0\n'
               '                  - name of the image file (two variables) '
               'if flag is 1\n'
               '*/\n'
               '{0:d}\n\n').format( self.homot_usage)
        return out
        

    def print_rot(self):
        """
        Print info about the rotation
        """
        out = ('/* ROTATION */\n'
               '/* 1. Rotation usage, integer (rotationUsage):\n'
               '        - 0: no rotation\n'
               '        - 1: rotation without tolerance\n'
               '        - 2: rotation with tolerance\n'
               '   2a. If rotationUsage == 1,\n'
               '          then for each angle,\n'
               '          first for azimuth, then for dip, and then for'
               ' plunge:\n'
               '             - Flag (0 / 1) indicating if given in an'
               ' image '
               'file,\n'
               '               followed by\n'
               '                  - one value (real) if flag is 0\n'
               '                  - name of the image file '
               '(one variable) if '
               'flag is 1\n'
               '   2b. If rotationUsage == 2,\n'
               '          then for each angle,\n'
               '          first for azimuth, then for dip, and then for '
               'plunge:\n'
               '             - Flag (0 / 1) indicating if given in an '
               'image file,\n'
               '               followed by\n'
               '                  - two values (lower and upper bounds) '
               '(real) if flag is 0\n'
               '                  - name of the image file (two variables) '
               'if flag is 1\n'
               '*/\n'
               '{0:d}\n\n').format(self.rot_usage)
        return out

    def print_consi(self):
        """
        Print info about the consistency
        """
        out = ('/* CONSISTENCY OF CONDITIONING DATA (TOLERANCE RELATIVELY '
               'TO THE RANGE OF TRAINING VALUES) */\n'
               '/* Maximal accepted expansion in both extremities of the '
               'range of values in training images\n'
               '   for covering the conditioning data values; e.g. if this '
               'number is set to 0.05,\n'
               '   the conditioning data values can be beyond the range of '
               'the values in the training images\n'
               '   (in both extremities) of at most 5 %; this separately '
               'applies to all variables for which\n'
               '   the distance is absolute (not relative, see "relative '
               'distance flag" below) and the\n'
               '   distance type is not 0 (see "distance type" below). */\n'
               '{0}\n\n').format(self.cond_data_tol)
        return out

    def print_norm(self):
        """
        Print info about the normalization
        """

        out = ('/* NORMALIZATION TYPE (FOR VARIABLES FOR WHICH DISTANCE '
               'TYPE IS NOT 0 AND DISTANCE IS ABSOLUTE) */\n'
               '/* Available types:\n'
               '      - NORMALIZING_LINEAR\n'
               '      - NORMALIZING_UNIFORM\n'
               '      - NORMALIZING_NORMAL */\n'
               '{0}\n\n').format( self.norm_type)
        return out

    def print_neigh_intro(self):
        """
        Print some introductory text to the 
        parameters related to the search neighborhood
        """
        out = ('/* SEARCH NEIGHBORHOOD PARAMETERS */\n'
               '/* A search neighborhood is a 3D ellipsoid, defined by:\n'
               '      - search radii (in number of nodes), for each '
               'direction\n'
               '      - anisotropy ratios, for each direction, i.e. '
               'numbers of nodes corresponding\n'
               '        to a distance of one, in each direction; for '
               'example (1.0, 1.0, 2.0) means\n'
               '        that the distance to the central node is the '
               'Euclidean distance where\n'
               '        the unit (distance=1) corresponds to 1, 1 and 2 '
               'nodes for the 1st, 2nd and\n'
               '        3rd direction respectively.\n'
               '      - angles (azimuth, dip and plunge) defining the '
               'rotation that sends the coordinates\n'
               "        system xyz onto the coordinates system x'y'z' in "
               "which the search radii\n"
               '        and the anisotropy ratios are given\n'
               '      - power at which the distance is elevated for '
               'computing the weight of each\n'
               '        node in the search neighborhood\n'
               '   Note that\n'
               '      - the search neighborhood is delimited by the search '
               'radii and the angles\n'
               '      - the anisotropy ratios are used only for computing '
               'the distance to the central\n'
               '        node, from each node in the search neighborhood\n'
               '      - the nodes inside the search neighborhood are '
               'sorted according to their\n'
               '        distance to the central node, from the closest '
               'one to the furthest one */\n')
        return out

    def print_neighs(self):
        """
        Print info related to the search neighborhood
        """
        out = ""

        for i, var in enumerate(self.vars):
            out = out + ( 
                '/* SEARCH NEIGHBORHOOD PARAMETERS FOR VARIABLE #{1} */\n'
                '{0.search_r_x}    '
                '{0.search_r_y}    '
                '{0.search_r_z}    // search radius in each direction\n'
                '{0.search_anis_x}    '
                '{0.search_anis_y}    '
                '{0.search_anis_z}    // anisotropy ratio in each '
                'direction\n'
                '{0.search_ang_az}    '
                '{0.search_ang_dp}    '
                '{0.search_ang_pl}    // angles (azimuth, dip, '
                'plunge in '
                'degrees) for rotation\n'
                '{0.search_pow}     // power for computing weight '
                'according '
                'to distance\n\n'
                ).format(
                var,i
                )
        return out

    def print_max_nod(self):
        """
        Print the max number of neighborihg nodes
        """
        out = (
            '/* MAXIMAL NUMBER OF NEIGHBORING NODES FOR EACH VARIABLE '
            '(as many number(s) as number of variable(s)) */\n'
            '{0}\n\n'
            ).format( 
            " ".join([str(var.max_nb_neigh) for var in self.vars])
            )
        return out

    def print_max_dens(self):
        """
        Print info about the max density of nodes
        """
        out = (
            '/* MAXIMAL DENSITY OF NEIGHBORING NODES IN SEARCH\n'
            '   NEIGHBORHOOD FOR EACH VARIABLE\n'
            '   (as many number(s) as number of variable(s)) */\n'
            '{0}\n\n'
            ).format(
            " ".join([str(var.max_dens_neigh) for var in self.vars])
            )
        return out

    def print_reldist(self):
        """
        Print info relative to the relative distance
        """
        out = (
            '/* RELATIVE DISTANCE FLAG FOR EACH VARIABLE\n'
            '   (as many flag(s) (0 / 1) as number of variable(s) */\n'
            '{0}\n\n'
            ).format(
            " ".join([str(int(var.rel_dist_flag)) for var in self.vars])
            )
        return out

    def print_disttype(self):
        """
        Print info related to the distance type
        """
        out = (
            '/* DISTANCE TYPE FOR EACH VARIABLE (as many number(s) as \n'
            '   number of variable(s)) */\n'
            '/* Available distance (between data events):\n'
            '      - 0: non-matching nodes (typically for categorical '
            'variable)\n'
            '      - 1: L-1 distance\n'
            '      - 2: L-2 distance\n'
            '      - 3: L-p distance, requires the real positive '
            'parameter p\n'
            '      - 4: L-infinity distance */\n'
            '{0}\n\n'
            ).format(
            " ".join([str(var.dist_type) for var in self.vars])
            )
        return out

    def print_weight(self):
        """
        Print info about the weighting factors
        """
        out = (
          '/* WEIGHT FACTOR FOR CONDITIONING DATA, FOR EACH VARIABLE '
            '(as many number(s) as number of variable(s)) */\n'
            '/* For the computation of distance between data events, '
            'if a value is a conditioning\n'
            '   data, its corresponding contribution is multiplied by '
            'the factor given here. */\n'
            '{0}\n\n'
          ).format(
            " ".join([str(var.weight) for var in self.vars])
            )
        return out

    def print_simpath(self):
        """
        Print info about the simulation and the path
        """
        out = (
            '/* SIMULATION AND PATH PARAMETERS */\n'
            '/* Key word for simulation type:\n'
            '      - SIM_ONE_BY_ONE:\n'
            '           successive simulation of one variable at one '
            'node in the simulation grid (4D path)\n'
            '      - SIM_VARIABLE_VECTOR:\n'
            '           successive simulation of all variable(s) at one '
            'node in the simulation grid (3D path) */\n'
            '{0}\n\n'
            '/* Key word for path type:\n'
            '      - PATH_RANDOM:\n'
            '           random path, for simulation type:\n'
            '              - SIM_ONE_BY_ONE     : path visiting all '
            'nodes and variables in a random order\n'
            '              - SIM_VARIABLE_VECTOR: path visiting all '
            'nodes in a random order\n'
            '      - PATH_UNILATERAL:\n'
            '           unilateral path, for simulation type:\n'
            '              - SIM_ONE_BY_ONE: requires a vector of size 4.\n'
            '                   Example: u = (0, -2, 1, 0) means that '
            'the path will visit all nodes:\n'
            '                   randomly in xv-sections, with increasing '
            'z-coordinate, and then decreasing y-coordinate.\n'
            '              - SIM_VARIABLE_VECTOR: requires a vector of '
            'size 3.\n'
            '                   Example: u = (-1, 0, 2) means that the '
            'path will visit all nodes:\n'
            '                   randomly in y-sections, with decreasing '
            'x-coordinate, and then increasing z-coordinate.\n'
            '           This vector must be given after the key word '
            'PATH_UNILATERAL. */\n'
            '{1}\n\n'
            ).format(
            self.sim_type,
            self.sim_path
            )
        return out

    def print_distthr(self):
        """
        Print info about the distance threshold
        """
        out = (
            '/* DISTANCE THRESHOLD FOR EACH VARIABLE (as many number(s) '
            'as number of variable(s)) */\n'
            '{0}\n\n').format(
            " ".join([str(var.dist_thre) for var in self.vars])
            )
        return out

    def print_proconst_intro(self):
        """
        Print some introductory text to the 
        parameters related to the probability constraints

        .. todo::
            This function should be completed with the full text
            
        .. warning::
            Probability constraints are still not implemented in 
            this version of the DeeSse interface.
        """
        out = ('/* PROBABILITY CONSTRAINTS */\n'
               '/* FOR EACH VARIABLE:\n'
               '   ...\n'
               '*/\n')
        return out

    def print_proconst(self):
        """
        Print info related to the probability constraints
        """
        out = ""

        for i, var in enumerate(self.vars):
            if var.prob_constr == 0:
                out = out + ( 
                    '/* PROBABILITY CONSTRAINTS FOR VARIABLE #{1} */\n'
                    '{0}\n\n'
                ).format(
                    var.prob_constr,i
                )
            elif var.prob_constr == 1:
                out = out + ( 
                    '/* PROBABILITY CONSTRAINTS FOR VARIABLE #{1} */\n'
                    '{0}\n{1}\n{2}\n{3}\n{4}\n{5}\n\n'
                ).format(
                    var.prob_constr,
                    var.prob_ln1,
                    var.prob_ln2,
                    var.prob_ln3,
                    var.prob_ln4,
                    var.prob_ln5,
                    i
                )
            else:
                print("    WARNING: Probability constraint set to 0")
                print("    Parameters changes not implemented in the current version.")
                out = out + ( 
                    '/* PROBABILITY CONSTRAINTS FOR VARIABLE #{1} */\n'
                    '{0}\n\n'
                ).format(
                    0,i
                )

        return out


    def print_blockdata_intro(self):
        """
        Print some introductory text to the 
        parameters related to the block data option

        .. todo::
            This function should be completed with the full text
            
        .. warning::
            Probability constraints are still not implemented in 
            this version of the DeeSse interface.
        """
        out = ('/* BLOCK DATA */\n'
               '/* FOR EACH VARIABLE:\n'
               '   1. Block data usage, integer (blockDataUsage):\n'
               '         - 0: no block data\n'
               '         - 1: block data: block mean value\n\n'
               '   2. If blockDataUsage == 1, then\n'
               '         - block data file name\n'
               '*/\n')
        return out


    def print_blockdata(self):
        """
        Print info related to the BLOCK_DATA option
        """
        out = ""

        for i, var in enumerate(self.vars):
            out = out + ( 
                '/* BLOCK DATA FOR VARIABLE #{1} */\n'
                '{0}\n\n'
                ).format(
                var.block_data,i
                )
        return out

    def print_maxscan(self):
        """
        Print info about the percentage of image to scan
        """
        out = (
            '/* MAXIMAL SCAN FRACTION FOR EACH TI (as many number(s) as '
            'number of training image(s)) */\n'
            '{0}\n\n').format(
            self.tis[0].max_scan
            )
        return out
    
    def print_tol(self):
        """
        Print info about the tolerance
        """
        out = (
            '/* TOLERANCE */\n'
            '/* Tolerance t on the threshold value for flagging nodes:\n'
            '   let d(i) be the distance between the data event in the '
            'simulation grid and in the training\n'
            '   image for the i-th variable and s(i) be the distance '
            'threshold for the i-th variable, and let\n'
            '   e(i) = max(0, (d(i)-s(i))/s(i)) be the relative error '
            'for the i-th variable, i.e. the relative\n'
            '   part of the distance d(i) beyond the threshold s(i); '
            'during the scan of the training image, a node\n'
            '   that minimizes e = max (e(i)) is retained (the scan is '
            'stopped if e = 0); if e is greater than the\n'
            '   tolerance t (given here), then the current simulated '
            'node and all non-conditioning nodes of the\n'
            '   data events (one per variable) in the simulation grid '
            'are flagged for resimulation (post-processing). */\n'
            '{0}\n\n'
            ).format(
            self.tol)
        return out

    def print_post(self):
        """
        Print info related to the post-processing procedure
        """
        out = (
            '/* POST-PROCESSING */\n'
            '/* 1. Maximal number of path(s) (npostProcessingPathMax)\n'
            '   2. If npostProcessingPathMax > 0:\n'
            '      key word for post-processing parameters (i. e. '
            'number of neighboring nodes, distance threshold,\n'
            '      maximal scan fraction, and tolerance):\n'
            '         - POST_PROCESSING_PARAMETERS_DEFAULT: for '
            'default parameters\n'
            '         - POST_PROCESSING_PARAMETERS_SAME   : for same '
            'parameters as given above\n'
            '         - POST_PROCESSING_PARAMETERS_MANUAL : for manual '
            'settings\n'
            '   3. If npostProcessingPathMax > 0 and '
            'POST_PROCESSING_PARAMETERS_MANUAL:\n'
            '         MAXIMAL NUMBER OF NEIGHBORING NODES FOR EACH '
            'VARIABLE (as many number(s) as number of variable(s))\n'
            '         MAXIMAL DENSITY OF NEIGHBORING NODES IN SEARCH '
            'NEIGHBORHOOD FOR EACH VARIABLE (as many number(s) as '
            'number of variable(s))\n'
            '         DISTANCE THRESHOLD FOR EACH VARIABLE (as many '
            'number(s) as number of variable(s))\n'
            '         MAXIMAL SCAN FRACTION FOR EACH TI (as many '
            'number(s) as number of training image(s))\n'
            '         TOLERANCE\n'
            '*/\n'
            '{0:d}\n'
            '{1}\n\n'
            ).format(
            self.post_proc_path_max,
            self.post_proc_par if self.post_proc_path_max > 0 else ""
            )
        return out

    def print_pyramids(self):
        """
        Print info related to the pyramids procedure
        """
        out = ('{0}{1}'
               '{2}\n\n'
            ).format(
                json_in_text["pyramids"]["title"],
                json_in_text["pyramids"]["text"],
                self.pyramids
            )
        return out

    def print_seed(self):
        """
        Print info about the seed
        """
        out = (
            '/* SEED NUMBER AND SEED INCREMENT */\n'
            '{0:d}\n'
            '{1:d}\n\n'
            ).format(
            self.seed,
            self.seed_inc
            )
        return out

    def print_realnb(self):
        """
        Print the number of realizations
        """
        out = (
            '/* NUMBER OF REALIZATIONS */\n'
            '{0:d}\n\nEND\n'
            ).format(self.real_nb)
        return out

    def print_file(self, file_name):
        '''
        Print the content of the parameters object into a file.

        Parameters:
            file_name: string
                Name of the file where to print.

        '''
        file_obj = open(file_name, 'w')
        print(self, file=file_obj)
        file_obj.close()

        
    pass
