from genericpath import exists
import os

def check_name(name):
    unkown_name = name.split("/")[-1]
    id = unkown_name.split(".")[0]
    pdb_path = name.split(".")[0]
    pdb_name = pdb_path + ".pdb"
    if exists(pdb_name):
        name = pdb_path
    else:
        name = "FILE NOT FOUND"
    out_path = "data/tmp/" 
    # full_path = os.path.abspath(out_path)
    # os.mkdir(full_path.join(Path("/" + id)))
    out_path += id #+ "/" + id 
    return name, out_path