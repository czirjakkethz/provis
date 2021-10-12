import subprocess

def pdb_to_xyzr_script(name):
    """
    Run the pdb_to_xyzr script for given filename
    Works only on linux
    (could easily be extended with other parameters)
    :param name: name - Name of file
    :param type: str
    :return: void - xyzr file
    """
    # print("./pdb_to_xyzr %4s.pdb > %4s.xyzr" % (name, name)
    rc = subprocess.call("./pdb_to_xyzr %4s.pdb > %4s.xyzr" % (name, name), shell=True)

