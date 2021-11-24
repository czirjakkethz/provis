import subprocess
import platform

def pdb_to_xyzr_script(name, solv, bash=0):
    """
    Run the pdb_to_xyzr script for given filename
    Works on both linux and windows
    (could easily be extended with other parameters)
    :param name: name - Name of file
    :param type: str
    :param name: solv - If True also convert solvent atoms, False only structural atoms
    :param type: bool
    :param name: bash - If you want to run original script, set to 1
    :param type: bool
    :return: void - xyzr file
    """
    # print("python pdb_to_xyz.py %1s %4s.pdb > %4s.xyzr" % (str(int(solv)), name, name))
    if not bash:
        ac = subprocess.call("python scripts/pdb_to_xyzr.py %1s %4s.pdb > %4s.xyzr" % (str(int(solv)), name, name), shell=True)
    else:
        ac = subprocess.call("./scripts/pdb_to_xyzr %4s.pdb > %4s.xyzr" % (name, name), shell=True)


# def main():
#     pdb_to_xyzr_script("2fd7", 1)

# if __name__ == "__main__":
#     main()
