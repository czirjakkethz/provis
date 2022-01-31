import imp
import pyvista
from provis.src.processing.file_converter import FileConverter
from provis.src.plotting.structure import Structure
from provis.src.plotting.surface import Surface
from provis.src.processing.name_checker import NameChecker
from provis.src.processing.data_handler import DataHandler
from provis.src.processing.surface_handler import SurfaceHandler

class Protein:
    """
    The protein class encapsulates every other class in the provis.src.processing package. 
    
    This class contains all the necessary information need for plotting.
    """
  
    def __init__(self, pdb_name, base_path=None, density=3.0, model_id=0):
        """
        Initialize the class with the name of the pdb file.
        
        Check if the molecule of the given model_id exists for this given .pdb file. If not set model_id to 0.

        Parameters:
            pdb_name: str
                Name of the pdb file with or without the extension. If the file is not stored in the data/pdb directory the full path has to be specified.
            base_path: str, optional
                Path to the root directory of the required directory structure (see: https://pro-vis.readthedocs.io/en/latest/setup.html#). Only needs to be specified if it is not the current working directory. Default: None.
            density: float, optional
                Default: 3.0.
            model_id: int, optional
                Scpecify the model id of the desired molecule from a trajectory file. If the passed model id is not in the trajectory then model_id will be automatically set to 0. Count starts at 0. Leave default value for static molecules. Default: 0.
        """
        self._model_id = model_id
        self._name_checker = NameChecker(pdb_name, base_path)
        self.file_converter = FileConverter(self._name_checker, density=density)
        self._data_handler =  DataHandler(self._name_checker, fc=self.file_converter)
        self._surface_handler = SurfaceHandler( self._name_checker, fc=self.file_converter, dh=self._data_handler, density=density)

        i = 0
        for model in self._data_handler._structure:
            i += 1
        self._num_models = i
        
        # If model id too large reset to 0
        if self._model_id >= self._num_models:
            self._model_id = 0