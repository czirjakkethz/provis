from genericpath import exists
import os
import provis

class NameChecker:
    """
    The point of this class is to provide uniform naming and path locations to all the other classes of provis.
    NameChecker has class variables and a class method to return these variables. 
    
    Variables: 
    pdb_name - Full path to the pdb file without the .pdb extension. Usually PROVIS_PATH/data/pdb/{pdb_id}.
    out_path - Full path to the temporary files. The names of all temporary files are derived from this variable. Usually PROVIS_PATH/data/tmp/{pdb_id}.
    base_path - Full path of the provis directory or any directory that has the following directory structure within: {path}/data/data, {path}/data/img, {path}/data/tmp, {path}/data/pdb. 
    """
    pdb_name = ""
    out_path = ""
    base_path = ""
    
    def __init__(self, name, base_path: str=None):
        """
        Checks if the provided name corresponds to an existing pdb file.
        Works if "name" points to a pdb file in the working directory or a pdb file in the default pdb file location: provis/data/pdb.
        Always sets "out_path" to default output location: provis/data/tmp.

        :param name: name - Name and path to location of the pdb file. If simply pdb id (file name without .pdb extension) is provided the default directory - provis/data/pdb - will be searched for the file. If no file found the input is returned.
        :param type: str
        :param name: base_path - Easiest if you leave this alone, clone the provis git repo and run everything in the provis root directory. Path to the base directory with the following file structure: {base_path}/data/data, {base_path}/data/img, {base_path}/data/tmp, {base_path}/data/pdb AND {base_path}/binaries that contains the binaries required to run provis. Defaults to the "provis" packages path.
        :param type: str, optional
        """
        if not base_path:
            path = os.path.dirname(provis.__file__)
            # path points to provis/provis, we only want provis/
            base_path = path[0:-6]
        
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
        NameChecker.pdb_name = name 
        NameChecker.out_path = out_path
        NameChecker.base_path = base_path
    
    @classmethod
    def return_all(cls):
        """
        Return all class variables.
        
        :return: str - path to the input pdb id (file without the extension)
        :return: str - path to the output pdb id - data/tmp/{pdb_id}
        :return: str - path to the root of the provis directory. Following file structure HAS TO exist within: {path}/data/data, {path}/data/img, {path}/data/tmp, {path}/data/pdb.
        """
        
        return NameChecker.pdb_name, NameChecker.out_path, NameChecker.base_path
