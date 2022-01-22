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
  
    def __init__(self):
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
        
            
        Protein(pdb_name, base_path=None, density=3.0, plot_solvent=False, msms=False, notebook=False):
        if notebook:
            pyvista.set_jupyter_backend('panel')
        self._name_checker = NameChecker(pdb_name, base_path)
        self.file_converter = FileConverter(self._name_checker, density=density)
        self._data_handler =  DataHandler(self._name_checker, fc=self.file_converter)
        self._surface_handler = SurfaceHandler( self._name_checker, fc=self.file_converter, dh=self._data_handler, density=density)
        for model in self._data_handler._structure:
            self.structure = Structure(self._name_checker, dh=self._data_handler, plot_solvent=plot_solvent, notebook=notebook)
            self.surface = Surface(self._name_checker, sh=self._surface_handler, msms=msms, notebook=notebook)

