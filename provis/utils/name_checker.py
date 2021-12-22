from genericpath import exists
import os

def check_name(name):
    """
    Checks if the provided name corresponds to an existing pdb file.
    Works if "name" points to a pdb file in the working directory or a pdb file in the default pdb file location: provis/data/pdb.
    Always sets "out_path" to default output location: provis/data/tmp.

    :param name: name - Name and path to location of the pdb file. If simply pdb id (file name without .pdb extension) is provided the default directory - provis/data/pdb - will be searched for the file. If no file found the input is returned.
    :param type: str
    
    :return: str - path to the input pdb id (file without the extension)
    :return: str - path to the output pdb id - data/tmp/{pdb_id}
    """
    name_pdb = name + ".pdb"
    default_loc_name = "data/pdb/" + name
    default_loc_name_pdb = default_loc_name + ".pdb"

    unkown_name = name.split("/")[-1]
    id = unkown_name.split(".")[0]

    if not (exists(name) or exists(name_pdb)):
        default = "data/pdb/" + id
        default_file = default + ".pdb"
        if exists(default_file):
            name = default
        else:
            print("The provided pdb file could not be found. Please make sure you have it downloaded and provided the correct file path.")

    out_path = "data/tmp/" 
    # full_path = os.path.abspath(out_path)
    # os.mkdir(full_path.join(Path("/" + id)))
    out_path += id #+ "/" + id 
    return name, out_path