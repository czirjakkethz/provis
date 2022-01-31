from email.mime import base
from genericpath import exists
import os
import provis

class NameChecker:
    """
    This class provides uniform names and path locations to all the other classes of provis.
    NameChecker has internal variables and a method to return these variables. 
    
    Class member variables: 
    self._pdb_name - Full path to the pdb file without the .pdb extension. Usually PROVIS_PATH/data/pdb/{pdb_id}.
    self._out_path - Full path to the temporary files. The names of all temporary files are derived from this variable. Usually PROVIS_PATH/data/tmp/{pdb_id}.
    self._base_path - Full path of the provis directory or any directory that has the following directory structure within: {path}/data/data, {path}/data/img, {path}/data/tmp, {path}/data/pdb. 
    """
    
    def __init__(self, name, base_path: str=None):
        """
        Checks if the provided name corresponds to an existing pdb file.
        The "name" variable is a path to the pdb file with or without the pdb extension. One can also simply pass the pdb ID (name without extension) if the file is saved under: "root directory"/data/pdb.
        Always sets "out_path" to default output location: "root directory"/data/tmp.
        
        This class is passed to all other classes and ensures that all the plotting classes plot the same pdb file.

        The code in words:
        Find the base directory (the root of the special directory structure specified in: https://pro-vis.readthedocs.io/en/latest/setup.html#setting-up-provis)
        Deduce the path to the .pdb file, the data/tmp and data/meshes directories.
        Return everything.
        
        Parameters:
            name: str
                Name and path to location of the pdb file. If simply pdb id (file name without .pdb extension) is provided the default directory - provis/data/pdb - will be searched for the file. If no file found the input is returned.
            base_path: str, optional
                Easiest if you leave this alone, clone the provis git repo and run everything in the provis root directory. Path to the base directory with the following file structure: {base_path}/data/data, {base_path}/data/img, {base_path}/data/tmp, {base_path}/data/pdb AND {base_path}/binaries that contains the binaries required to run provis. Defaults to the "provis" packages path.
        """
        # disable warning from external binaries. Provis works regardless and clutters output.
        import warnings
        warnings.filterwarnings("ignore")

        if not base_path:
            path = os.getcwd() #path.dirname(provis.__file__)
            base_path = path
            if base_path.split("/")[-2] == "provis":
                len_f = len(base_path.split("/")[-1])
                base_path = base_path[0:-len_f]
            
        if base_path[-1] != "/":
            base_path += "/"
        
        name = name.split(".")[0]
        name_pdb = name + ".pdb"

        unkown_name = name.split("/")[-1]
        id = unkown_name.split(".")[0]

        
        if not (exists(name_pdb)):
            # name = name.split(".")[0]
            default = base_path + "data/pdb/" + id
            default_file = default + ".pdb"
            if exists(default_file):
                name = default
            else:
                print("The provided pdb file could not be found. Please make sure you have it downloaded and provided the correct file path.")


        out_path = base_path + "data/tmp/" 
        out_path += id
        mesh_path = base_path + "data/meshes/" + id
        self._pdb_name = name 
        self._out_path = out_path
        self._base_path = base_path
        self._mesh_path = mesh_path
    
    def return_all(self):
        """
        Return all class variables. This function is used by all other provis classes to retrieve the paths to the files needed.
        
        Returns: 
            str
                path to the input pdb id (file without the extension)
            str
                path to the output pdb id - data/tmp/{pdb_id}
            str
                path to the root of the provis directory. Following file structure HAS TO exist within: {path}/data/data, {path}/data/img, {path}/data/tmp, {path}/data/pdb.
        """
        
        return self._pdb_name, self._out_path, self._base_path, self._mesh_path
