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
    A module containing some classes and function useful to work with
    structured grids.

:This file:
   `grid.py`

:Version:
    0.3 , 2013-08-11 :
        * Adapted to the case of grids made of cells or grid made of
          points.
    0.2 , 2013-03-06 :
        * Converted to Python3.X with `2to3`
    0.1 , 2012-04-02 :
        * First version.
 
:Authors:
    Alessandro Comunian

.. note::
   The grids for the moment are always considered as 3D.

'''

import numpy as np

class Grid:
    '''
    A simple class that contains the *Origin*, *Delta* and *Size* of a
    simulation.  It can also be used as container for some information
    contained in a VTK structured grid file.
  
    .. note::
        * By default, the size of the grid is considered in term of points.

    .. seealso::
        :py:mod:`vtknumpy`

    '''
    
    def __init__(self, \
                 ox=0.0, oy=0.0, oz=0.0, \
                 dx=1.0, dy=1.0, dz=1.0, \
                 nx=200, ny=200, nz=1, gtype='points'):
        """
        Define a structured grid, with some default values. 

        Parameters:

            ox, oy, oz: float, optional
               
                Coordinates of the origin.

            dx, dy, dz: float, optional

                Distance between points of the structured grid or side
                of the cells.

            nx, ny, nz: int, optional

                Number of points. The number of cells is nx-1, ny-1
                and nz-1.
            
            gtype: string in ('points', 'cells')

                Define the grid type, for grids made of points or
                grids made of cells.

        """
        self.ox = ox
        self.oy = oy
        self.oz = oz
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.nx = nx        
        self.ny = ny        
        self.nz = nz
        self.gtype = gtype

        if self.gtype == 'points':
            self.points = self.nx*self.ny*self.nz
        elif self.gtype == 'cells':
            self.cells = ( self.nx-1 if self.nx> 1 else 1)*(
                self.ny-1 if self.ny> 1 else 1
                )*( self.nz-1 if self.nz> 1 else 1)
                           
        # Compute the size of the grid
        self._lx = None #self.dx * self.nx - self.ox
        self._ly = None #self.dy * self.ny - self.oy 
        self._lz = None #self.dz * self.nz - self.oz

      
    
    def __str__(self):
        """
        Print some of the informations provided in the class.
        """
        out = ("Grid info:\n" +
            "    Origin: ( {0.ox:f}, {0.oy:f}, {0.oz:f})\n".format(self) + 
            "    Delta:  ( {0.dx:f}, {0.dy:f}, {0.dz:f})\n".format(self) +
            "    Size:   ( {0.lx:f}, {0.ly:f}, {0.lz:f})\n".format(self) +
            "    N:      ( {0.nx:d}, {0.ny:d}, {0.nz:d})\n".format(self) +
            "    type:     {0.gtype}\n".format(self) +
            "    points:   {0.points}\n".format(self) +
            "    cells:    {0.cells}\n".format(self) )
        return out


    def print_intervals(self, axis='xyz'):
        '''
        Print the intervals that constitute the simulation domain in a
        format like::

            [ ox, ox+nx*dx] [ oy, oy+ny*dy] [ oz, oz+nz*dz]

        where *ox* is the origin, *nx* is the number of points and *dx*
        is the delta between points (*idem* for *y* and *z*).

        Parameters:
            axis: string containing ['x','y','z'], optional
                If the default value "xyz" is used, then all the intevals
                are printed.
        
        '''
        
        ox = self.ox; nx = self.nx; dx = self.dx
        oy = self.oy; ny = self.ny; dy = self.dy
        oz = self.oz; nz = self.nz; dz = self.dz


        if 'x' in axis:
            print("        interval along *x* [ ", ox, ", ", nx*dx + ox, "] ", end=' ')
        if 'y' in axis:
            print("        interval along *y* [ ", oy, ", ", ny*dy + oy, "] ", end=' ')
        if 'z' in axis:
            print("        interval along *z* [ ", oz, ", ", nz*dz + oz, "] ", end=' ') 

        print('')


    def spacing(self):
        """
        To print out the spacing of the grid as a tuple

        Parameters:
            self : an instance of the Grid class

        Returns:
            A tuple  containing the spacing defined in the grid.
        """
        return (self.dx, self.dy, self.dz)


    def origin(self):
        """
        To print out the origin of the grid as a tuple

        Parameters:
            self : an instance of the Grid class

        Returns:
            A tuple containing the origin defined in the grid.
        """
        return (self.ox, self.oy, self.oz)


    def shape(self):
        """
        To print out the shape of the grid as a tuple

        Parameters:

            self:

                an instance of the Grid class

            gtype: string in ('points', 'cells')
            
                String to decide to print the shape in terms of points
                or in terms of cells.

        Returns:
            A tuple containing the origin defined in the grid.
        """
        if self.gtype == 'points':
            out = (self.nx, self.ny, self.nz)
        elif self.gtype == 'cells':
            # By default, nx, ny and nz contains the size of the grid
            # in terms of points. When one dimension is 1, then along
            # that direction number of cells == number of points.
            out = (self.nx-1 if self.nx > 1 else 1,
                   self.ny-1 if self.ny > 1 else 1,
                   self.nz-1 if self.nz > 1 else 1)
        return out


    def compute_max(self):
        """
        Compute the max values for *x*, *y* and *z* of the grid.
        """
        self.x_max = self.ox + self.dx*self.nx
        self.y_max = self.oy + self.dy*self.ny
        self.z_max = self.oz + self.dz*self.nz


    def get_lx(self):
        """
        Provide as output a tuple containing the size of a
        grid. Useful for the implementation of ``property``.
        """
        return self.dx * self.nx - self.ox

    def get_ly(self):
        """
        Provide as output a tuple containing the size of a
        grid. Useful for the implementation of ``property``.
        """
        return self.dy * self.ny - self.oy


    def get_lz(self):
        """
        Provide as output a tuple containing the size of a
        grid. Useful for the implementation of ``property``.
        """
        return self.dz * self.nz - self.oz


    def set_lx(self, val=None):
        self._lx = self.dx * self.nx - self.ox

    def set_ly(self):
        self._ly = self.dy * self.ny - self.oy

    def set_lz(self):
        self._lz = self.dz * self.nz - self.oz

    def del_lx(self):
        del self._lx

    def del_ly(self):
        del self._ly

    def del_lz(self):
        del self._lz

    lx = property(get_lx, set_lx, del_lx, "'size' along *x* of the grid.")
    ly = property(get_ly, set_ly, del_ly, "'size' along *y* of the grid.")
    lz = property(get_lz, set_lz, del_lz, "'size' along *z* of the grid.")


    def get_points(self):
        """
        Update the number of points for a points grid
        """
        return self.nx*self.ny*self.nz

    def set_points(self, val=None):
        self._points = self.nx*self.ny*self.nz

    def del_points(self):
        del self._points

    points = property(get_points, set_points, del_points,
                      "Number of points")
    
    def get_size(self):
        """
        Compute the total number of points/cells in the grid.
        """
        return self.nx*self.ny*self.nz

    def get_cells(self):
        """
        Update the number of cells for a cells grid
        """
        return \
            (self.nx-1 if self.nx>1 else 1)* \
            (self.ny-1 if self.ny>1 else 1)* \
            (self.nz-1 if self.nz>1 else 1)

    

    def set_cells(self, val=None):
        self._cells = \
            (self.nx-1 if self.nx>1 else 1)* \
            (self.ny-1 if self.ny>1 else 1)* \
            (self.nz-1 if self.nz>1 else 1)

    def del_cells(self):
        del self._cells

    cells = property(get_cells, set_cells, del_cells, "Number of cells")
    

    def print4gslib(self):
        """
        A function to print out the grid information in the header of
        a GSLIB file.
        """
        out = ('{0.nx:d} {0.ny:d} {0.nz:d} '
               '{0.dx:f} {0.dy:f} {0.dz:f} '
               '{0.ox:f} {0.oy:f} {0.oz:f} '
               '\n'.format(self))
        return out
        
    def get_center(self):
        """
        Returns the center of the grid
        """
        x_c = 0.5*(self.ox + self.nx*self.dx)
        y_c = 0.5*(self.oy + self.ny*self.dy)
        z_c = 0.5*(self.oz + self.nz*self.dz)
        return x_c, y_c, z_c
    
def sg_info(data):
    """
    Collect some information about a structured grid
    contained in a numpy array
    
    
    """
    if data.dtype == np.int:
        unq = np.unique(data)
        tot = data.size
        tot_ge0 = np.sum(data>=0)
        props = np.bincount(data[data>=0])/tot_ge0

        
    return(unq, props, tot_ge0)
    
if __name__ == '__main__':

    # Create a grid with the default parameters
    grid1 = Grid()

    print(grid1)

    grid1.print_intervals('xz')

    print("    Spacing:", grid1.spacing())
    print("    Origin: ", grid1.origin())
    print("    Shape:  ", grid1.shape())
    print("    Max:  ", grid1.compute_max())
    

    print("Before update")
    print(grid1)
    grid1.dx = 1.1
    print("After update")
    print(grid1)
    
    
