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
    A module with some utilities to convert from numpy to gslib and
    vice-versa, and to support some GSLIB based software.

    For more details about the GSLIB software, see [Deutsch1988]_ .

:This file:
    ``gslibnumpy.py``


:Version:
    0.9 , 2020-01-20 :
        * Adapted to the inclusion in the s2Dcd package on Github.
    0.8 , 2015-12-02 :
        * Added `gslib_slice`.
        * General clean-up
        * Improved numpy2gslib to allow printing many variables
        * Improved the pandas usage in numpy2gslib
    0.7 , 2014-11-13 :
        * Introduced dependencies from "pandas" to load txt faster
    0.6 , 2013-08-26 :
        * Improved how errors are handled in `gslib2numpy` and the
          extension available.
    0.5 , 2013-07-07 :
        * Improved the function `gslib_points2numpy`.
    0.4 , 2013-03-06 :
        * Converted to Python3 using `2to3`.
        * Improved the flexibility of `gslib2numpy`.
        * Remover the dependence from the module `gslibnumpy_f`.
    0.3 , 2012-11-19 :
        * Last version before the conversion to Python3.
        * Module `gslib2numpy` improved to be able to read also files
          containing multiple variables.
        * merged the function `dat2numpy` into `gslib2numpy`.
    0.2 , 2012-11-18 :
        * Some improvement using the suggestions of Spyder.
        * Deleted the option for the "dtype" in function
          "gslib_points2numpy".
        * included a new function to load GSLIB files coming from Isatis
          (dat2numpy, for files ``.dat``)
    0.1 , 2012-04-27 :
        * Implemented in FORTRAN90 the subroutines `numpy2gslib`.
    0.0 , 2012-03-07:
        * First version.
        
:Authors:
    Alessandro Comunian


References
--------------

.. [Deutsch1988] Deutsch, C. V. & Journel, A.G. (1988) GSLIB:
       Geostatistical Software Library and User's Guide 2nd
       edition. Oxford University Press


Functions and classes
------------------------------
'''

__version__ = 0.8
__author__ = "Alessandro Comunian"

import os
import numpy as np
import pandas as pd

import s2Dcd.utili as utili
import s2Dcd.ext as ext
import s2Dcd.grid as grid


# Default values of some variables
DX = 1.0
DY = 1.0
DZ = 1.0
OX = 0.0
OY = 0.0
OZ = 0.0

# Default values for variables related to variograms
# - Number of structures
STR_NB = 1
# - Nugget effect
NUGGET = 0.0
# Value of the category (or threshold) associated with the variogram
CAT_THR = 0

# Default value for the *c* parameter ("sill")
C_DEF = 1.0
# Default value for the range
A_DEF = 100.0
# Default value for the max lag to plot a variogram
H_MAX = 200.0
# Default value for the variogram structure (1=spherical)
STYPE = 1
# Default value for the sill (c parameter)
C_PAR = 1.0
# Default value for the angles
ANGLES = [0.0, 0.0, 0.0]
# Default value for the ranges
RANGES = [A_DEF, A_DEF, A_DEF]
# The default exponent for the Power variogram model
OMEGA = 1.5
VARIO_TYPE = []
VARIO_TYPE.append("sph.")
VARIO_TYPE.append("exp.")
VARIO_TYPE.append("gau.")
VARIO_TYPE.append("pow.")
VARIO_TYPE.append("hol.")


# Default name for the variables saved in a GSLIB file
VAR_NAMES=("data",)

def numpy2gslib(data, file_name, var_names=VAR_NAMES, grd=None,
                float_format=None, verbose=False):
    '''
    Convert a tuple of numpy arrays into a GSLIB ASCII file.

    Parameters:
        data : tuple of numpy arrays
            The numpy arrays to be saved as GSLIB file.
            
        file_name : string (optional)
            The name of the GSLIB file where to save the data.
        varname : tuple of string, optional
            The names of the variable to be stored in the GSLIB file.
            (one for each numpy array in data).
        grd : grid object, optional
            If available, get the dimension of the grid from a grid object.
        float_format: string (optional, default None)
            This is useful to pass to pandas `to_csv` the float format required.
        verbose : boolean flag, optional
            Print or not some info useful for debugging...

    .. note::
        * Works only with 1D, 2D or 3D numpy arrays.
        * When more than one variable is provided, the order provided
          in "var_names" is respected when printing the columns.
        * `data` should be **always** provided as a tuple, even when
          it contains one array only, otherwise it may result in a corrupted
          file.
    '''

    if not isinstance(data, (tuple, list)):
        print("    WARNING IN numpy2gslib:")
        print("    - Input data must be a tuple or a list of numpy array!")
        print("      (tuple/list should contain at least one element)")

    # Check if we are dealing with 1D, 2D or 3D data
    #
    # HERE FOR THE MOMENT CHECK ONLY THE DIMENSIONS OF THE 1ST NUMPY
    # ARRAY IN THE TUPLE... IN THE FUTURE I SHOULD CHECK THE
    # DIMENSIONS OF ALL THE COMPONENTS OF THE TUPLE
    #
    data_dim = len(data[0].shape)
    #
    # HERE IN THE FUTURE I SHOULD CHECK IF THE LEN OF DATA IS THE SAME
    # AS THE NUMBER OF VARIABLES NAMES.
    #
    data_len = len(var_names)

    # Check the dimensionality of the data, 3D, 2D or 1D
    if data_dim == 3:
        nx, ny, nz = data[0].shape
    elif data_dim == 2:
        nx, ny = data[0].shape
        nz = 1
    elif data_dim == 1:
        nx = data[0].shape
        ny = 1
        nz = 1
    else:
        print("    WARNING: Wrong dimensions of the input data")
        print("        Should be a value in (1,2,3), but *{0}* "
              "found".format(data_dim))

    with open(file_name, 'w') as file_obj:
        # Write the 1st line, which can be a comment, or contain some
        # useful information about the grid.
        if grd is None:
            # Create a grid using the size of data and the default
            # parameters
            grd = grid.Grid(nx=nx, ny=ny, nz=nz,
                            dx=DX, dy=DY, dz=DZ,
                            ox=OX, oy=OY, oz=OZ)
#        print(grd)
        file_obj.write(grd.print4gslib())

        # Print the number of variables
        file_obj.write("{0:d}\n".format( data_len))
        
        for var_name in var_names:
            # Write the name of the variable
            file_obj.write('{0}\n'.format( var_name))

        # Empy dictionary to use in the creation of the data frame
        df_dict = {}
        for (var_name, datum) in zip(var_names, data):
            df_dict[var_name] = np.reshape(datum, nx*ny*nz, order="F")
            
        # Create the dataframe from the dict
        df = pd.DataFrame(df_dict)
        
        # Print on the file
        df.to_csv(file_obj, index=False, header=False, sep="\t", 
                  columns=var_names, float_format=float_format)

    if verbose:
        print('\n    Output file name: "{0}"'.format(file_name))
        print("    Dimensions of the dataset: {0}".format(data[0].shape))
        print("    Number of variables: {0:d}".format(data_len))


def gslib2numpy(file_name, verbose=True, formats=None, names=None,
                source=None, grd=None):
    '''
    Convert a GSLIB ASCII grid into a numpy array.

    Can convert only a subclass of GSLIB ASCII files.
    See note for more details.

    Parameters:
        file_name: string
            Name of the GSLIB input file.
        verbose: bool, optional
            A flag to print reading information or not.
        formats: list of characters, optional (default=None)
            The definition of the input format, :py:func:`numpy.genfromtxt`
            like.
            If nothing is provided, float values are considered for each
            variable.
        names: tuple, optional (default=None)
            The definition of the names of the variables. If not
            provided, the names of the variables contained in the
            GSLIB file are used.

        source: string in ['impala','isatis', 'fluvsim'],
            optional (default=None)
            The header changes a little if the file is to be used with
            "impala", or it is an output from "isatis".  The functions
            tries to extrapolate from the file extension if the file
            comes from "impala" (.gslib) or from "Isatis" (.dat).  If
            this argument is provided, then the format of the header
            is forced according to this (see the note for more
            details).  The output files from FLUVSIM (.out) have a
            header that contains only the dimensions of the grid, and
            therefore they fall as a simple case of "implala" files.
            
    Returns:
        out_dict:
            A numpy dictionary of arrays with the data contained in
            the GSLIB file.
        in_grd:
            A "grid.Grid" object, containing the grid definition (if
            provided).

    .. note::
        Limited to a subsclass of GSLIB ASCII files.
        The files can have an "impala" (of FLUVSIM) format, like::

            nx ny nz [ dx dy dz [ ox oy oz] ]
            nb_var
            var_name[0]
            var_name[1]
            ...
            var_name[nb_var-1]
            var0_value[0] var1_value[0] ...
            var0_value[1] var1_value[1] ...
            ...           ...           ...

        Or, in case the software is "Isatis", then the header will be
        like::
        
            Description line
            nb_var nx ny nz ox oy oz dx dy dz
            ...
   
    .. seealso:: :py:func:`numpy.genfromtxt`
        
    '''

    try:
        in_file = open(file_name, 'r', encoding='utf-8')
        if verbose :
            print('\n    Reading file: "%(file_name)s"' % locals())
    except IOError:
        print('    Error reading file "%s"'%(file_name))
        print('    Check if the file exists...')


    # Try to get the type of header from the file extension
    if source is None:
        in_ext = os.path.splitext(file_name)[1]
        impala_ext = [ext.GSLIB, ext.OUT]
        isatis_ext = [ext.DAT]
        allowed_ext = impala_ext + isatis_ext
        if in_ext in impala_ext:
            source = "impala"
        elif in_ext in isatis_ext:
            source = "isatis"
        else:
            raise ValueError('Invalid extension "{0},"'\
                             'not in {1}'.format(in_ext, allowed_ext))
            
    
    if grd:
        # If a grid definition in provided as input, use that
        in_grd = grd
    else:
        in_grd = grid.Grid()

    # Check if the input format is "Isatis" or "impala", or check the 
    # definition option
    if source == 'impala':
        # Read the first line, containing the size of the data
        line = in_file.readline()
        line_split = line.split()

        if grd is None:
            in_grd.nx = int(  line_split[0])
            in_grd.ny = int(  line_split[1])
            in_grd.nz = int(  line_split[2])
            try:
                in_grd.dx = float(line_split[3])         
                in_grd.dy = float(line_split[4])         
                in_grd.dz = float(line_split[5])
            except IndexError:
                in_grd.dx = DX
                in_grd.dy = DY
                in_grd.dz = DZ
            try:
                in_grd.ox = float(line_split[6])
                in_grd.oy = float(line_split[7])
                in_grd.oz = float(line_split[8])
            except IndexError:
                in_grd.ox = OX
                in_grd.oy = OY
                in_grd.oz = OZ

        # Read the 2nd line, containing the number of variables
        var_nb = int(in_file.readline())

    elif source == 'isatis':
        # Read the header line
        header = in_file.readline()
    
        # Read the line containing the grid description
        grd_descr = in_file.readline()    
        line_split = grd_descr.split()
        
        var_nb = int(line_split[0])
        if grd is None:
            in_grd.nx = int(line_split[1])
            in_grd.ny = int(line_split[2])
            in_grd.nz = int(line_split[3])
            try:
                in_grd.dx = float(line_split[3])         
                in_grd.dy = float(line_split[4])         
                in_grd.dz = float(line_split[5])
            except IndexError:
                in_grd.dx = DX
                in_grd.dy = DY
                in_grd.dz = DZ
            try:
                in_grd.ox = float(line_split[6])
                in_grd.oy = float(line_split[7])
                in_grd.oz = float(line_split[8])
            except IndexError:
                in_grd.ox = OX
                in_grd.oy = OY
                in_grd.oz = OZ


    else:
        print("    ERROR in function `gslib2numpy`")

    
    # Read the 3rd line, containing the name of the variables
    if names is None:
        var_names = [in_file.readline().split()[0] for i in range(var_nb)]
    else:
        var_names = names
        
    if formats is None:
        formats = ','.join(['f' for i in range(var_nb)])
    

    if verbose :
        print('    Grid dimensions: %d x %d x %d' % \
            (in_grd.nx, in_grd.ny, in_grd.nz))
        print('    Number of variables: %(var_nb)d' % locals())
        print('    Variable names: %(var_names)s' % locals())
        print('    Formats: %(formats)s' % locals())
    
    in_file.close()


    # Read the data...
    # NOTE: SOME IMPROVEMENTS MIGHT INCLUDE THE USE OF THE OTHER OPTIONS
    # LIKE usecols and delimiter...

    data = np.genfromtxt(file_name, skip_header=2+var_nb, 
                         dtype=formats, names=var_names)

    out = {}
    for var_name in var_names:
        out[var_name] = np.reshape( data[var_name],
                                    (in_grd.nx, in_grd.ny, in_grd.nz),
                                    order='F')
        
    return out, in_grd

def gslib2numpy_onevar(file_name, verbose=True, formats=None,
                       names=None, source=None):
    '''
    Convert a GSLIB ASCII grid into a numpy array.

    This is a simplified version of ``gslib2numpy`` that can be used
    when you only have one variable and you don't care about the grid
    information contained in the file.

    Can convert only a subclass of GSLIB ASCII files.
    See note for more details.

    Parameters:
        file_name: string
            Name of the GSLIB input file.
        verbose: bool, optional
            A flag to print reading information or not.
        formats: list of characters, optional (default=None)
            The definition of the input format, :py:func:`numpy.genfromtxt`
            like.
            If nothing is provided, float values are considered for each
            variable.
        names: tuple, optional (default=None)
            The definition of the names of the variables. If not provided,
            the names of the variables contained in the GSLIB file are used.
        source: string in ['impala','isatis'], optional (default=None)
            The header changes a little if the file is to be used with 
            "impala", or it is an output from "isatis".
            The functions tries to extrapolate from the file extension if
            the file comes from "impala" (*.gslib) or from "Isatis" (*.dat).
            If this argument is provided, then the format of the header is
            forced according to this (see the note for more details).
            
    Returns:
        out_dict:
            A numpy dictionary of arrays with the data contained in
            the GSLIB file.  NON NON .... ONLY A VARIABLE! UPDATE THIS
            DOCUMENTATION IF IT WORKS.
        in_grd:
            A "grid.Grid" object, containing the grid definition (if
            provided).

    .. note::
        Limited to a subsclass of GSLIB ASCII files.
        The files can have an "impala" format, like::

            nx ny nz [ dx dy dz [ ox oy oz] ]
            nb_var
            var_name[0]
            var_name[1]
            ...
            var_name[nb_var-1]
            var0_value[0] var1_value[0] ...
            var0_value[1] var1_value[1] ...
            ...           ...           ...

        Or, in case the software is "Isatis", then the header will be like::
        
            Description line
            nb_var nx ny nz ox oy oz dx dy dz
            ...
   
    .. seealso:: :py:func:`numpy.genfromtxt`
        

    .. TODO::

        I should carefully revise this function and simplify using
        more directly and only the pandas' function.
    '''

    try:
        in_file = open(file_name, 'r', encoding='utf-8')
        if verbose :
            print('\n    Reading file: "%(file_name)s"' % locals())
    except IOError:
        print('    Error reading file "%s"'%(file_name))
        print('    Check if the file exists...')
        return -1

    in_grd = grid.Grid()

    # Try to get the type of header from the file extension
    if source is None:
        in_ext = os.path.splitext(file_name)[1]
        if in_ext == ".gslib":
            source = "impala"
        elif in_ext == ".dat":
            source = "isatis"
        else:
            print("    ERROR in function `gslib2numpy`,")
            print("    unknown file extension (%s),"%(in_ext))
            print('    valid extensions are ".gslib" or ".dat".')
            return -2
            
    

    # Check if the input format is "Isatis" or "impala", or check the 
    # definition option
    if source == 'impala':
        # Read the first line, containing the size of the data
        line = in_file.readline()
        line_split = line.split()

        in_grd.nx = int(  line_split[0])
        in_grd.ny = int(  line_split[1])
        in_grd.nz = int(  line_split[2])
        try:
            in_grd.dx = float(line_split[3])         
            in_grd.dy = float(line_split[4])         
            in_grd.dz = float(line_split[5])
        except IndexError:
            in_grd.dx = DX
            in_grd.dy = DY
            in_grd.dz = DZ
        try:
            in_grd.ox = float(line_split[6])
            in_grd.oy = float(line_split[7])
            in_grd.oz = float(line_split[8])
        except IndexError:
            in_grd.ox = OX
            in_grd.oy = OY
            in_grd.oz = OZ

        # Read the 2nd line, containing the number of variables
        var_nb = int(in_file.readline())

    elif source == 'isatis':
        # Read the header line
        header = in_file.readline()
    
        # Read the line containing the grid description
        grd_descr = in_file.readline()    
        line_split = grd_descr.split()
        
        var_nb = int(line_split[0])
        in_grd.nx = int(line_split[1])
        in_grd.ny = int(line_split[2])
        in_grd.nz = int(line_split[3])
        in_grd.ox = float(line_split[4])
        in_grd.oy = float(line_split[5])
        in_grd.oz = float(line_split[6])
        in_grd.dx = float(line_split[7])         
        in_grd.dy = float(line_split[8])         
        in_grd.dz = float(line_split[9])

    else:
        print("    ERROR in function `gslib2numpy`")
        return -1
    
    # Read the 3nd line, containing the name of the variables
    if names is None:
        var_names = [in_file.readline().split()[0] for i in range(var_nb)]
    else:
        var_names = names
        
    if formats is None:
        formats = ','.join(['f' for i in range(var_nb)])
    

    if verbose :
        print('    Grid dimensions: %d x %d x %d' % \
            (in_grd.nx, in_grd.ny, in_grd.nz))
        print('    Number of variables: %(var_nb)d' % locals())
        print('    Variable names: %(var_names)s' % locals())
        print('    Formats: %(formats)s' % locals())
    
    in_file.close()


    # Read the data...
    # NOTE: SOME IMPROVEMENTS MIGHT INCLUDE THE USE OF THE OTHER OPTIONS
    # LIKE usecols and delimiter...

    data = pd.read_csv(file_name, skiprows=2).values


    out = np.reshape( data, (in_grd.nx, in_grd.ny, in_grd.nz), order='F')
        
    return out, in_grd





    

def numpy2dat(file_name, data, grd=None, varname=None):
    '''
    Convert a numpy array into a GSLIB ASCII file in the *Isatis* format.

    Parameters:
        file_name : string
            The name of the GSLIB file where to save the data.
        
        data : numpy array
            The numpy array to be saved as GSLIB file.
        grd : grid.Grid object, optional
            The definition of the grid for the dataset.
        varname : string, optional
            The name of the variable to be stored in the GSLIB file.

    .. note::
        * Works only with 1D, 2D or 3D numpy arrays.
        * Only one variable per file.
        * In some cases the input array must be flatten.
    '''
    # THIS FUNCTION COULD BE MERGED WITH numpy2gslib...
    
    err_code = 0
    def_varname = 'data'    
    file_obj = file(file_name, 'w')
    file_obj.write('Created with numpy2dat\n')
    if grd:
        file_obj.write('1 %d %d %d %f %f %f %f %f %f 1\n'% \
            (grd.nx, grd.ny, grd.nz, grd.ox, grd.oy, grd.oz, \
            grd.dx, grd.dy, grd.dz))
    else:
        file_obj.write('1 %d %d %d 0.0 0.0 0.0 1.0 1.0 1.0 1\n'
                       % (data.shape()))

    if not varname:
        varname = def_varname
    file_obj.write('%s\n'%(varname))
    
    if data.dtype in ['float']:
        np.savetxt(file_obj, data, fmt='%12.4e')
    elif data.dtype in ['int']:
        np.savetxt(file_obj, data, fmt='%d')
    else:
        print("    ERROR: Wrong input dtype (only 'int' and "
              "'float' are OK)")
        print("           Value provided: ", data.dtype)
        file_obj.close()
        err_code = -1
        
       
    file_obj.close()
    return err_code
     

        
        
def gslib_points2numpy(file_name, verbose=True, formats=None):
    '''
    Convert to a GSLIB ASCII point data file into a numpy array.

    Parameters:
        file_name : string
            Name of the GSLIB input file.
        verbose : bool, optional
            A flag to print out information or not.
        formats: list of characters, optional (default=None)
            The definition of the input format,
            :py:func:`numpy.genfromtxt` like.  If nothing is provided,
            then `genfromtxt` tries to understand the format
            automatically.
            
    Returns:
        An `ndarray` containing a key for each variable contained in
        the GSLIB file and a string containing the header of the file. 


    Example:
        Read the content of the file `test.gslib`::

        >>> data = gslib_points2numpy("test.gslib")
        >>> x = data['x']

        
    '''


    try:
        in_file = open(file_name, 'r')
        if verbose :
            print('\n    Reading file: "%(file_name)s"' % locals())
    except IOError:
        print('    ERROR reading file "%s"'%(file_name))
        print('    function call aborted.')
        return 


    # Read the first line, containing a header
    header = in_file.readline()

    # Read the 2nd line, containing the number of variables
    var_nb = int(in_file.readline())

    # Read the names of the variables.
    var_names = [in_file.readline().split()[0] for i in range(var_nb)]
    
    in_file.close()

    if verbose:
        print('    Number of variables: %d'%(var_nb))
        print('    ', var_names)

    data = np.genfromtxt(file_name, skip_header=2+var_nb, \
                                dtype=formats, names=var_names, unpack=True)

    return data, header

    
def numpy2gslib_points(xyz_data, file_name, varname="data"):
    '''
    Convert a numpy array into a GSLIB points ASCII file.

    Parameters:
        xyz_data : tuple of numpy arrays
            The tuple is made of four numpy arrays of the same lenght:
            the coordinates x, y, z and the point data *data*.
            Depending on the size of the tuple, the points are considered
            1D, 2D or 3D.
        file_name : string
            The name of the GSLIB file where to save the data.
        varname : string, optional
            The name of the variable.

    Returns:
        A gslib point file containing the coordinates of the points
        and the data values.  If the size of the input array is <1,
        then the value -1 is returned and a warning message is printed
        out.

    .. note::
        * The default names for the variables are *x*, *y*, *z* and *data*.
        * The data type is automatically detected from the data type of the 
          data.
    '''

    data_size = xyz_data[0].size
    data_dim = len(xyz_data)

    if data_size < 1:
        print('    Warning: the numpy array to be converted into'
              ' gslib points is empty!')
        return -1

    header = ("%i\n") % (data_size)

    data = np.column_stack(xyz_data)

    if data_dim == 2 :
        # This is a 1D array
        header = header + "2\nx\n%s\n"%(varname)
        fmt = '%.{0}e {1}'.format(utili.SIG_FIG, utili.dtype_fmt(xyz_data[1]))
    elif data_dim == 3 :
        # This is a 2D array
        header = header + "3\nx\ny\n%s\n"%(varname)
        fmt = '%.{0}e %.{0}e {1}'.format(utili.SIG_FIG,
                                         utili.dtype_fmt(xyz_data[2]))     
    elif data_dim == 4 :
        # This is a 3D array
        header = header + "4\nx\ny\nz\n%s\n"%(varname)
        fmt = '%.{0}e %.{0}e %.{0}e {1}'.format(utili.SIG_FIG,
                                         utili.dtype_fmt(xyz_data[3]))     
    else :
        print('    Error in "numpy2gslib_points"')
        print('    Must be between 2 and 4 '
              '(but is is %(data_dim)d.' % locals())

    with open(file_name, 'wb') as fid:
        fid.write(header.encode('utf-8'))
        np.savetxt(fid, data, fmt=fmt)


def spherical(h, c=C_DEF, a=A_DEF):
    '''
    Definition of the model of a Spherical variogram.
    
    Parameters:
        h: float
            The values of lag :math:`h` where to plot the variogram
        c: float
            The sill value :math:`c`
        a: float
            The actual range :math:`a`

    Returns:
        The values of the variogram :math:`\gamma(h)`
    '''

    return(np.where(h <= a, c*(1.5*h/a-0.5*(h/a)**3), c))

def exponential(h, c=C_DEF, a=A_DEF):
    '''
    Definition of the model of a Exponential variogram.
    
    Parameters:
        h: float
            The values of lag :math:`h` where to plot the variogram
        c: float
            The sill value :math:`c`
        a: float
            The actual range :math:`a`

    Returns:
        The values of the variogram :math:`\gamma(h)`
    '''

    return(c*(1.0-np.exp(-3*h/a)))


def gaussian(h, c=C_DEF, a=A_DEF):
    '''
    Definition of the model of a Gaussian variogram.
    
    Parameters:
        h: float
            The values of lag :math:`h` where to plot the variogram
        c: float
            The sill value :math:`c`
        a: float
            The actual range :math:`a`

    Returns:
        The values of the variogram :math:`\gamma(h)`
    '''

    return(c*(1.0-np.exp(-(3*h)**2/a**2)))


def power(h, c=C_DEF, omega=OMEGA):
    '''
    Definition of the model of a power variogram.
    
    Parameters:
        h: float
            The values of lag :math:`h` where to plot the variogram
        c: float
            The sill value :math:`c`
        omega: float
            The exponent of the model :math:`\omega`

    Returns:
        The values of the variogram :math:`\gamma(h)`
    '''

    return(c*h**omega)


class VarioStruct(object):
    '''
    A class to contain info about one structure of a variogram
    '''

    def __init__(self, stype=STYPE, c_par=C_DEF, angles=ANGLES,
                 ranges=RANGES):
        self.stype = stype
        self.c_par = c_par
        self.angles = angles
        self.ranges = ranges
        pass

    def __str__(self):
        '''
        Print a string with the information contained in the class.
                
        '''
        out =  "        - type: {0.stype}\n".format(self)
        out += "        - c par.: {0.c_par}\n".format(self)
        out += "        - angles: " + ' '.join(
            str(elem) for elem in self.angles)+"\n"
        out += "        - ranges: "+' '.join(
            str(elem) for elem in self.ranges)+"\n"

        return out

    def print4par(self):
        '''
        Print a string with the information contained in the class
        useful in a parameter file.
                
        '''
        out =  "{0.stype} ".format(self)
        out += "{0.c_par} ".format(self)
        out += (' '.join(str(elem) for elem in self.angles) +
                "\\    it,cc,ang1,ang2,ang3\n")
        out += (' '.join(str(elem) for elem in self.ranges) +
                "\\    a_hmax,a_hmin,a_vert\n")

        return out

    def print4csv(self):
        '''
        Print a string with the information contained in the class
        useful in a CSV file.
                
        '''
        out =  "{0.stype} ".format(self)
        out += "{0.c_par} ".format(self)
        out += ' '.join(str(elem) for elem in self.angles)+' '
        out += ' '.join(str(elem) for elem in self.ranges)

        return out
        

class Variogram(object):
    '''
    A class to contain all the parameters related to a GSLIB variogram
    model.
    '''

    def __init__(self, str_nb=STR_NB, nugget=NUGGET, cat_thr=CAT_THR):
        '''
        Initializer class for Variogram
        
        Parameters:
            str_nb: int
                Default number of structures
        '''
        # The number of structures
        self.str_nb = str_nb
        self.nugget = nugget
        self.cat_thr = cat_thr

        self.structs = []
        for struct in range(self.str_nb):
            self.structs.append(VarioStruct())
            
        

    def __str__(self):
        '''
        Print the info about a variogram
        '''
        
        out =  "    Variogram information\n"
        out += "    ---------------------\n"
        out += "    - category/threshold: {0.cat_thr}\n".format(self)
        out += "    - nb. of structures : {0.str_nb}\n".format(self)
        out += "    - nugget            : {0.nugget}\n".format(self)
        for i_st, struct in enumerate(self.structs):
            out += "    - Structure {0}\n".format(i_st+1)
            out += str(struct)
        
        return out

    def print4par(self):
        '''
        Print the info about a variogram for a parameter file
        '''
        out = ("{0.str_nb} {0.nugget} \{0.cat_thr}   nst,"
               " nugget effect\n".format(self))
        for i_st, struct in enumerate(self.structs):
             out += struct.print4par()
        
        return out

    def print4csv(self):
        '''
        Print the info about a variogram for a CSV file
        '''
        out = "{0.str_nb} {0.nugget} ".format(self)
        for i_st, struct in enumerate(self.structs):
            out += struct.print4csv()+' '
        
        return out





    def plot(self, file_name, h_max=H_MAX):
        '''
        Plot the variogram model using matplotlib
        
        Parameters:
            h_max: float
                Maximum value of the lag
        
        '''
        import pylab as pl

        # Use latex for the labels...
        # pl.rc('text', usetex=True)
        # pl.rc('font', family='serif')
        h = np.arange(h_max)
    
        gamma = np.zeros(h.shape)

        for struct in self.structs:
            if struct.stype == 1:
                # Spherical model
                if any(t > 0 for t in struct.angles):
                    print("    WARNING, angles > 0 not implemented")
                # PER IL MOMENTO USO SOLO IL PRIMO RANGE...
                gamma += spherical(h, struct.c_par, struct.ranges[0])
            elif struct.stype == 2:
                # Exponential model
                gamma += exponential(h, struct.c_par, struct.ranges[0])
            elif struct.stype == 3:
                # Gaussian model
                gamma += gaussian(h, struct.c_par, struct.ranges[0])
            else: 
                print("    ERROR: Structure type not implemented")

        pl.plot(h, gamma)
        pl.ylim((0.0, 1.05*max(gamma)))
        pl.xlabel('h')
        pl.ylabel('gamma(h)')
        pl.savefig(file_name)


    def plot_3dir(self, file_name, h_max=None, suptitle=None):
        '''
        Plot the variogram model using matplotlib along all the three
        direction for all the available ranges.
        
        Parameters:
            h_max: float
                Maximum value of the lag
            suptitle: string
                A string for the super-title.

        .. note::
            Use by default the range of the 1st structure to define the 
            upper limit of the :math:`x` axis.
        '''
        import pylab as pl

        # A multiplier to the max range, better to define it as a
        # parameter...
        h_mult = 2
        n_points = 100

        # Use latex for the labels...
        #        pl.rc('text', usetex=True)
        #        pl.rc('font', family='serif')
        if h_max is None:
            
            h_max = max(self.structs[0].ranges)
        h = np.arange(0.0, h_mult*h_max, h_mult*h_max/n_points)
        
        gamma = np.zeros((h.shape[0], 3))

        for struct in self.structs:
            if struct.stype == 1:
                # Spherical model
                if any(t > 0 for t in struct.angles):
                    print("    WARNING, angles > 0 not implemented")
                # PER IL MOMENTO USO SOLO IL PRIMO RANGE...
                for i in range(3):
                    gamma[:,i] += spherical(h, struct.c_par, 
                                            struct.ranges[i])
            elif struct.stype == 2:
                # Exponential model
                for i in range(3):
                    gamma[:,i] += exponential(h, struct.c_par,
                                              struct.ranges[i])
            elif struct.stype == 3:
                # Gaussian model
                for i in range(3):
                    gamma[:,i] += gaussian(h, struct.c_par,
                                           struct.ranges[i])
            else: 
                print("    ERROR: Structure type not implemented")

        figu, (ax1, ax2, ax3) = pl.subplots(1, 3, sharey=True,
                                            figsize=(9,3))

        if suptitle:
            pl.suptitle(suptitle)

        ax1.plot(h, gamma[:,0])
        ax1.set_title('horizontal, max.')

        ax2.plot(h, gamma[:,1])
        ax2.set_title('horizontal, min.')
        
        ax3.plot(h, gamma[:,2])
        ax3.set_title('vertical')

        ax1.set_ylim((0.0, 1.05*np.max(gamma)))
        ax1.set_xlabel('h')
        ax1.set_ylabel('gamma(h)')
        pl.tight_layout()
        pl.savefig(file_name)

        

        

        

        


def align_parfile(par_file):
    '''
    Read a string as a parfile, detect where the character "\" is and
    align everything.

    Parameters:
        par_file: string
            A string containing the parameter file to be indented.

    Returns:
        A string containing the parameter file indented.
    '''

    par_split = par_file.split("\n")

    # Create a list containing the positions of the \\
    bar_pos = []
    for line in par_split:
        bar_pos.append(line.find("\\") )
    
    max_pos = max(bar_pos)
    
    out = ""
    for line in par_split:
        # Add spaces wher the position of \ is smalles than the max position
        bar_pos_curr = line.find("\\")
        shift = max_pos - bar_pos_curr
        if bar_pos_curr > 0:
            out += line[:bar_pos_curr]+" "*shift+line[bar_pos_curr:]+"\n"
        else:
            out += line+"\n"

    return out



def gslib_slice(file_in, axis, level=None):
    """
    Cut a 2D slice within a 3D `GSLIB` file perpendicular to a given
    *axis* at a given *level*. The `GSLIB` file can contain more than
    one variable.

    Parameters:
        file_in : string,
            The name of the input file.
        axis : char, in ['x', 'y', 'z']
            The axis coordinate which is used to cut the slice.
        level : int, optional (default=None)
            The index (**not** the true coordinate) along the
            coordinate where `axis` is defined where to cut the slice.
            If *None* then an index in the middle of the input
            domain is selected.

    Returns:
        A file with the same name as `file_in` and a suffix
        `-<axis><level>.gslib` is created, and 0 if successful.
        If the return value is < 0 then there was some problem...

    .. warning::
        If in the target directory there is a file with the same name as the
        file created as return value this is overwritten.

    """

    data, grid = gslib2numpy(file_in)

    file_out = (file_in[:-len(ext.GSLIB)] + "-" + axis + str(level) +
                ext.GSLIB)

    if axis == 'x':
        if level is None:
            level = int(0.5*grid.nx)
        grid.ox = grid.ox + grid.dx*level
        grid.nx=1
        numpy2gslib(tuple(data[key][level, :, :] for key in data.keys()),
                    file_out, tuple(key for key in data.keys()), grid)
    elif axis == 'y':
        if level is None:
            level = int(0.5*grid.ny)
        grid.oy = grid.oy + grid.dy*level
        grid.ny=1
        numpy2gslib(tuple(data[key][:,level, :] for key in data.keys()),
                    file_out, tuple(key for key in data.keys()), grid)
    elif axis == 'z':
        if level is None:
            level = int(0.5*grid.nz)
        grid.oz = grid.oz + grid.dz*level
        grid.nz=1
        numpy2gslib(tuple(data[key][:,:,level] for key in data.keys()),
                    file_out, tuple(key for key in data.keys()), grid)
    else:
        print("    ERROR in `gslib_slice`:")
        print("    Value provided for `axis` (%s)" % (axis), end=' ')
        print(" not in ['x', 'y', 'z'].")

    return file_out



#
# Here all the testing functions...
#
if __name__ == '__main__':


    import numpy as np
    
    # 
    # test the function `numpy2gslib`
    #
    print("\n"
          "    **************************************\n"
          "    ** Testing the function numpy2gslib **\n"
          "    ** with one variable (test1)        **\n"
          "    **************************************\n"
          )
    data_int = np.random.randint(10, size=(4,5,3))
    numpy2gslib((data_int,), "./test/gslibnumpy_test1.gslib", verbose=True)

    print("\n"
          "    ***************************************\n"
          "    ** Testing the function numpy2gslib  **\n"
          "    ** with two variables (test2)        **\n"
          "    ***************************************\n"
          )
    data_float = np.random.uniform( size=(4,5,3))
    numpy2gslib((data_int,data_float), 
                "./test/gslibnumpy_test2.gslib",
                ("data_int", "data_float"),
                verbose=True)

    print("\n"
          "    ***************************************\n"
          "    ** Testing the function gslib_slice  **\n"
          "    **     (test3)                       **\n"
          "    ***************************************\n"
          )
    gslib_slice("./test/gslibnumpy_test2.gslib", 'x', level=2)


#     data, header = gslib_points2numpy("./tests/data_ref/gslib_points2numpy-2D.gslib")

#     print(header)
#     print(data)

#     vario_test = Variogram(str_nb = 2)
#     vario_test.structs[1].stype = 1
# #    print(vario_test)
#     print(vario_test.print4csv())

#     vario_test.plot("vario_test.pdf")
#     vario_test.plot_3dir("vario_test_3dir.pdf")

    

    
    

    
    

    
    

    
