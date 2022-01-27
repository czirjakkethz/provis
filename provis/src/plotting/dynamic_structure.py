import pyvista
import numpy as np
import pyvista as pv
import time

from provis.src.processing.file_converter import FileConverter
from provis.src.plotting.structure import Structure
from provis.src.plotting.surface import Surface
from provis.src.processing.name_checker import NameChecker
from provis.src.processing.data_handler import DataHandler
from provis.src.processing.surface_handler import SurfaceHandler

class DynamicStructure:
    """
    The Dynamic Structure class, similarly to the Protein class, encapsulates every other class and creates a user friendly way to plot your desired dynamic structure of a protein molecules.
    
    While the class is built similarly to the Protein class it does not use the Protein class itself. This is due to the fact that the Protein class is a rigid class made for a single molecule and
    """
  
    def __init__(self, pdb_name, base_path=None, density=3.0, plot_solvent=False, msms=False, notebook=False):
        """
        Initialize the class with the name of the pdb file and you are ready for plotting!
        
        Set the optional variables to get the most out of provis.

        :param name: pdb_name - Name of the pdb file with or without the extension. If the file is not stored in the data/pdb directory the full path has to be specified.
        :param type: str
        :param name: base_path - Path to the root directory of the required directory structure (see: https://pro-vis.readthedocs.io/en/latest/setup.html#). Only needs to be specified if it is not the current working directory. Default: None.
        :param type: str, optional
        :param name: density - Default: 3.0.
        :param type: float, optional
        :param name: plot_solvent - If True solvent atoms will also be plotted. Default: False.
        :param type: bool, optional
        :param name: msms - Set to True if you want to compute the surface information using the msms binary. If False surface will be computed natively. Default: False.
        :param type: bool, optional
        :param name: notebook - Set to True when using running in a Jupyter Notebook environment. Default: False.
        :param type: bool, optional
        """
        
            
#        Protein(pdb_name, base_path=None, density=3.0, plot_solvent=False, msms=False, notebook=False)
        import warnings
        warnings.filterwarnings("ignore")
        self._solvent = plot_solvent
        if notebook:
            pyvista.set_jupyter_backend('pythreejs')
        self._notebook = notebook
        self._name_checker = NameChecker(pdb_name, base_path)
        self._path, self._out_path, self._base_path, self._mesh_path = self._name_checker.return_all()
        self.file_converter = FileConverter(self._name_checker, density=density)
        self._data_handler =  DataHandler(self._name_checker, fc=self.file_converter)
        self._surface_handler = SurfaceHandler( self._name_checker, fc=self.file_converter, dh=self._data_handler, density=density)
                
        self.structure = Structure(self._name_checker, dh=self._data_handler, plot_solvent=plot_solvent, notebook=notebook)
        self.surface = Surface(self._name_checker, sh=self._surface_handler, msms=msms, notebook=notebook)
        print("Initialized DynamicStructure class")
        
    def plot_atoms(self, box=False, res=None, outname=None, camera=None, title="Atoms"):
        """
        Plot the dynamic atom cloud.

        :param name: box, optional - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_stick_point.png.
        :param type: string, optional
        :param name: camera - Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to 'xy'. Default: None.
        :param type: pyvista.Camera, optional
        :param name: title - Title of the plotting window. Default: "Atoms".
        :param type: str, optional
        """
        atom_data = self._data_handler.get_atoms(show_solvent=self._solvent, model_id=0)
            
        # Create a plotter object and initialize first mesh
        plotter = pv.Plotter(notebook=self._notebook, off_screen=False)
        plotter.smooth_shading = True
        plotter.background_color = 'grey'
        _atoms, _, _atom_names = self._data_handler.get_atom_mesh(atom_data, vw=0) # second arg: 1 = showvw spheres instead of "normal" radius
        bigmesh = pv.PolyData()
        colors = []
        for j, mesh in enumerate(_atoms):
            bigmesh += mesh
            colors += [_atom_names[j]] * len(mesh.points)
        plotter.add_mesh(bigmesh, scalars=colors, lighting=True, smooth_shading=True)

        # Open the movie file
        plotter.open_movie("animation.mp4", 1)


        # Update mesh
        i = 0
        plotter.enable_eye_dome_lighting()
        plotter.render()
        plotter.write_frame()
        i = 0
        for model in self._data_handler._structure:
            atom_data = self._data_handler.get_atoms(show_solvent=self._solvent, model_id=i) # second arg: 1 = show_solvent
            
            
            plotter.clear()
            plotter.smooth_shading = True
            _atoms, _, _atom_names = self._data_handler.get_atom_mesh(atom_data, vw=0) # second arg: 1 = showvw spheres instead of "normal" radius
            bigmesh = pv.PolyData()
            colors = []
            for j, mesh in enumerate(_atoms):
                bigmesh += mesh
                colors += [_atom_names[j]] * len(mesh.points)
            plotter.add_mesh(bigmesh, scalars=colors, lighting=True, smooth_shading=True)
            
            # must update normals when smooth shading is enabled
            plotter.mesh.compute_normals(cell_normals=False, inplace=True)
            plotter.render()
            plotter.write_frame()
            time.sleep(0.5)
            i+=1

    def plot_surface(self, feature=None, patch=False, title="Surface", box=None, res=None, outname=None, camera=None):
        """
        Plot the dynamic atom cloud.


        :param name: feature - Pass which feature (coloring) you want to plot. Options: hydrophob, shape, charge. Default: None (uniform coloring).
        :param type: str, optional
        :param name: patch - If True then coloring will be read in from "root directory"/data/tmp/{pdb_id}.pth file. Default: False.
        :param type: bool, optional
        :param name: box, optional - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_stick_point.png.
        :param type: string, optional
        :param name: camera - Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to 'xy'. Default: None.
        :param type: pyvista.Camera, optional
        :param name: title - Title of the plotting window. Default: "Surface".
        :param type: str, optional
        """
  
        # Create and structured surface
        mesh, cas = self._surface_handler.return_mesh_and_color(False, feature=feature, patch=patch, model_id=0)
        if not cas:
            np.zeros(len(mesh.points))
            
        # Create a plotter object and initialize first mesh
        plotter = pv.Plotter(notebook=self._notebook, off_screen=False)
        plotter.add_mesh(mesh, scalars=cas, cmap='RdBu', lighting=True, smooth_shading=True, show_edges=False)

        # Open the movie file
        plotter.open_movie("animation.mp4", 1)


        # Update mesh
        i = 0
        plotter.enable_eye_dome_lighting()
        plotter.render()
        plotter.write_frame()
        for model in self._data_handler._structure:
            mesh, cas = self._surface_handler.return_mesh_and_color(False, feature=feature, patch=patch, model_id=i)
            if not cas:
                np.zeros(len(mesh.points))
            
            i += 1
            
            plotter.clear()
            plotter.smooth_shading = True
            plotter.add_mesh(mesh, scalars=cas, cmap='RdBu', lighting=True, smooth_shading=True, show_edges=False)

            # must update normals when smooth shading is enabled
            plotter.mesh.compute_normals(cell_normals=False, inplace=True)
            plotter.render()
            plotter.write_frame()
            time.sleep(0.5)
            

        # Closes and finalizes movie
        plotter.close()
        
