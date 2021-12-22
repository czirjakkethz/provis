import pyvista as pv
import numpy as np
import trimesh

from provis.src.processing.surface_handler import SurfaceHandler
from provis.utils.name_checker import check_name

class Surface:
    """
    The Surface class is used to visualize the surface information of the given molecule.
    
    All member functions of the class are a "two in one" function:
    - plot the surface created by the MSMS binary. 
    - create a good approximation of the surface natively with o3d and trimesh and plot it. MSMS does not have to be installed for this option. (It is slightly slower and less precise)
    
    Choose between the two by setting the msms Boolean variable (True corresponding to the MSMS binary option is default)
    
    """
    def __init__(self, name, dens=None):
        """
        Initialize Surface class with given filename. Creates internal data structures; a DataHandler to extract basic surface information and stores it in self._atmsurf (this is a list of Spheres for each atom roughly equating the Van-der-Waals radius).
        
        :param name: dens - sampling density used in msms binary. Also needed to load the face and vert files, as their (file)names include the density
        :param type: float
        :param name: name - name of file to be loaded
        :param type: str
        """
        self._path, self._out_path = check_name(name)
        if dens:
            self._density = dens
        self._sh = SurfaceHandler(name, dens=dens)

    def plot_surface(self, msms=True, outname=0):
        """
        Plot the surface of protein.
        
        :param name: outname - save image of plot to specified filename. Will appear in data/img/ directory. default: data/img/{self._out_path}_surface
        :param type: string
        
        :return: void - plot
        """
        if msms:
            mesh = self._sh.return_mesh_and_color(simple=True)
        
        else:
            mesh = self._sh.native_mesh()
        
        #plot
        pl = pv.Plotter(lighting=None)
        pl.background_color = 'grey'
        pl.enable_3_lights()

        style = 'surface'
        pl.add_mesh(mesh, color="white", smooth_shading=True, style=style, show_edges=False)        
        # save a screenshot
        if not outname:
            new_name = self._out_path.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = 'data/img/' + new_name + '_surface.png'
        pl.show(screenshot=outname)

    def plot_hydrophob(self, outname="hydrophob"):
        """
        Plot the hydrophobic features of a protein.

        :param name: outname - Name of output file. Defaults to "hydrophob".
        :param type: str, optional
        """
        
        mesh, cas = self._sh.return_mesh_and_color(feature="hydrophob")

        # plot surface with feature visualization
        pl = pv.Plotter()
        pl.add_mesh(mesh, scalars=cas, cmap='RdBu')
        pl.background_color = 'grey'
        pl.camera_position = 'xy'
        pl.show(screenshot=outname + '.jpeg')

    def plot_shape(self, outname="shape"):
        """
        Plot the shape features of a protein.

        :param name: outname - Name of output file. Defaults to "shape".
        :param type: str, optional
        """
        
        mesh, cas = self._sh.return_mesh_and_color(feature="shape")

        # plot surface with feature visualization
        pl = pv.Plotter()
        pl.add_mesh(mesh, scalars=cas, cmap='RdBu')
        pl.background_color = 'grey'
        pl.camera_position = 'xy'
        pl.show(screenshot=outname + '.jpeg')

    def plot_charge(self, outname="charge"):
        """
        Plot the charge features of a protein.

        :param name: outname - Name of output file. Defaults to "charge".
        :param type: str, optional
        """
        
        mesh, cas = self._sh.return_mesh_and_color(feature="charge")

        # plot surface with feature visualization
        pl = pv.Plotter()
        pl.add_mesh(mesh, scalars=cas, cmap='RdBu')
        pl.background_color = 'grey'
        pl.camera_position = 'xy'
        pl.show(screenshot=outname + '.jpeg')