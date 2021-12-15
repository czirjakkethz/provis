import subprocess
import platform
import os.path

def msms_script(name, dens):
    """
    Run the msms script for given filename and density.
    Can work on all operating systems if you have msms.exe binary installed.
    (could easily be extended with other parameters for binary, currently only dens)
    
    :param name: name - Name of file
    :param type: str
    :param name: dens - Density of triangulation
    :param type: float
    
    :return: void - face and vert files
    """
    # print("./msms.exe -if %4s.xyzr -of %4s_out_%s -density %s" % (name, name, str(int(dens)), str(dens)))
    if os.path.isfile(".scripts/msms.exe"):
        if platform.system() == "Windows":
            ac = subprocess.call("scripts/msms.exe -if %4s.xyzr -of %4s_out_%s -density %s" % (name, name, str(int(dens)), str(dens)), shell=True)
        else:
            ac = subprocess.call("./scripts/msms.exe -if %4s.xyzr -of %4s_out_%s -density %s" % (name, name, str(int(dens)), str(dens)), shell=True)


# def main():
#     msms_script("2fd7", 3.0)

# if __name__ == "__main__":
#     main()
