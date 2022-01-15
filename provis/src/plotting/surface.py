import pyvista as pv
import numpy as np
import trimesh

from provis.src.processing.surface_handler import SurfaceHandler
from provis.src.processing.name_checker import NameChecker

class Surface:
    """
    The Surface class is used to visualize the surface information of the given molecule.
    
    All member functions of the class are a "two in one" function:
    - plot the surface created by the MSMS binary. 
    - create a good approximation of the surface natively with o3d and trimesh and plot it. MSMS does not have to be installed for this option. (It is slightly slower and less precise)
    
    Choose between the two by setting the msms Boolean variable (True corresponding to the MSMS binary option is default)
    
    """
    def __init__(self, nc, sh=None, density=None, msms=False, notebook=False, dh=None):
        """
        Initialize Surface class with given filename. Creates internal data structures; a DataHandler to extract basic surface information and stores it in self._atmsurf (this is a list of Spheres for each atom roughly equating the Van-der-Waals radius).
        
        :param name: nc - Instance of a NameChecker class. Used to pass the pdb file name and paths.
        :param type: NameChecker
        :param name: msms - If True plot msms binary version of surface. If False plot the native (non-binary) surface. Default: False.
        :param type: bool, optional
        :param name: density - sampling density used in msms binary. Also needed to load the face and vert files, as their (file)names include the density
        :param type: float, optional
        :param name: notebook - Needs to be set to true to work in a notebook environment. Defualts to False.
        :param type: bool, optional 
        :param name: dh - Instance of DataHandler. To be passed to self._sh (SurfaceHandler), not needed otherwise. Default: None.
        """
        
        self._path, self._out_path, self._base_path = nc.return_all()
        if density:
            self._density = density
        if not sh:
            self._sh = SurfaceHandler(nc, density=density, dh=dh)
        else:
            self._sh = sh
        self._msms = msms
        self._notebook = notebook
        self._shading = not self._notebook


    def plot_surface(self, outname=None, feature=None, patch=None):
        """
        Plot the surface of protein.
        
        :param name: outname - save image of plot to specified filename. Will appear in data/img/ directory. Default: data/img/{self._out_path}_surface.
        :param type: string, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot.
        """
        mesh, cas = self._sh.return_mesh_and_color(self._msms, feature, patch)
        
        #plot
        pl = pv.Plotter(notebook=self._notebook)
        pl.background_color = 'grey'
        pl.enable_3_lights()

        pl.add_mesh(mesh, scalars=cas, cmap='RdBu', smooth_shading=self._shading, show_edges=False)
        # pl.add_mesh(mesh, color="white", smooth_shading=True, style=style, show_edges=False)        
        
        # save a screenshot
        if not outname:
            new_name = self._out_path.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_surface.png'
        pl.show(screenshot=outname, title="Surface")

    def plot_hydrophob(self, outname=None):
        """
        Plot the hydrophobic features of a protein.

        :param name: outname - Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_hydrophob.png.
        :param type: string, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot.
        """
        mesh, cas = self._sh.return_mesh_and_color(self._msms, feature="hydrophob")

        # plot surface with feature visualization
        pl = pv.Plotter(notebook=self._notebook)
        pl.add_mesh(mesh, scalars=cas, cmap='RdBu', smooth_shading=self._shading, show_edges=False)
        pl.background_color = 'grey'
        pl.camera_position = 'xy'
        
        # save a screenshot
        if not outname:
            new_name = self._out_path.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_hydrophob.png'
        pl.show(screenshot=outname, title="Hydrophob")

    def plot_shape(self, outname=None):
        """
        Plot the shape features of a protein.

        :param name: outname - Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_shape.png.
        :param type: string, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot.
        """
        
        mesh, cas = self._sh.return_mesh_and_color(self._msms, feature="shape")

        # plot surface with feature visualization
        pl = pv.Plotter(notebook=self._notebook)
        pl.add_mesh(mesh, scalars=cas, cmap='RdBu', smooth_shading=self._shading, show_edges=False)
        pl.background_color = 'grey'
        pl.camera_position = 'xy'
        
        # save a screenshot
        if not outname:
            new_name = self._out_path.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_shape.png'
        pl.show(screenshot=outname, title="Shape")

    def plot_charge(self, outname=None):
        """
        Plot the charge features of a protein.

        :param name: outname - Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_charge.png.
        :param type: string, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot.
        """
        
        mesh, cas = self._sh.return_mesh_and_color(self._msms, feature="charge")

        # plot surface with feature visualization
        pl = pv.Plotter(notebook=self._notebook)
        pl.add_mesh(mesh, scalars=cas, cmap='RdBu', smooth_shading=self._shading, show_edges=False)
        pl.background_color = 'grey'
        pl.camera_position = 'xy'
        
        # save a screenshot
        if not outname:
            new_name = self._out_path.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_charge.png'
        pl.show(screenshot=outname, title="Charge")