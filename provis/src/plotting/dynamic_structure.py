import pyvista
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
        if notebook:
            pyvista.set_jupyter_backend('panel')
        self._name_checker = NameChecker(pdb_name, base_path)
        self.file_converter = FileConverter(self._name_checker, density=density)
        self._data_handler =  DataHandler(self._name_checker, fc=self.file_converter)
        self._surface_handler = SurfaceHandler( self._name_checker, fc=self.file_converter, dh=self._data_handler, density=density)
                
        self.structure = Structure(self._name_checker, dh=self._data_handler, plot_solvent=plot_solvent, notebook=notebook)
        self.surface = Surface(self._name_checker, sh=self._surface_handler, msms=msms, notebook=notebook)
        
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
        i = 0
        for model in self._data_handler._structure:
            self.structure.plot(atoms=1, box=box, vw=0, bonds=0, residues=0, res=res, outname=outname, title=title, camera=camera, model_id=i, dynamic=True)
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
        i = 0
        for model in self._data_handler._structure:
            self.surface.plot(feature=feature, title=title, patch=patch, box=box, res=res, outname=outname, camera=camera, model_id=i)
            i+=1
