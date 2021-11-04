import subprocess
import platform

def msms_script(name, dens):
    """
    Run the msms script for given filename and density
    Works on both linux and windows
    (could easily be extended with other parameters)
    :param name: name - Name of file
    :param type: str
    :param name: dens - Density of triangulation
    :param type: float
    :return: void - face and vert files
    """
    # print("./msms.exe -if %4s.xyzr -of %4s_out_%s -density %s" % (name, name, str(int(dens)), str(dens)))
    if platform.system() == "Windows":
        ac = subprocess.call("msms.exe -if %4s.xyzr -of data/%4s_out_%s -density %s" % (name, name, str(int(dens)), str(dens)), shell=True)
    else:
        ac = subprocess.call("./msms.exe -if %4s.xyzr -of data/%4s_out_%s -density %s" % (name, name, str(int(dens)), str(dens)), shell=True)


# def main():
#     msms_script("2fd7", 3.0)

# if __name__ == "__main__":
#     main()
