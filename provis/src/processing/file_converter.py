import subprocess
import os.path
from provis.src.processing.name_checker import NameChecker
from provis.utils.surface_utils import output_pdb_as_xyzrn

class FileConverter():
    """
    Class to create and destroy necessary files required in other parts of code.
    It has a bunch of member functions to call the binaries and scripts to convert the files.
    
    Also has a cleanup function that removes everything in the "root_directory"/data/tmp (and data/img if specified) directory. Best practice is to call this function at the end of your main file.
    However if you want to plot the same protein many times, then it is benificial to keep the temporary (data/tmp) files as if they exist provis will not recompute them.
    """
    
    def __init__(self, nc, density=3.0, convert_all=False):
        """
        Can be constructed empty. 
        If "convert_all" set to True conversions instant. Creates xyzr and mol2 files in every case and pqr, face and vert files if appropriate binary exists.
        
        :param name: nc - Instance of a NameChecker class. Used to pass the pdb file name and paths.
        :param type: str, optional
        :param name: density - Density of triangles. Default: 3.0.
        :param type: float, optional
        :param name: base_path - Path to "working directory" according to the rules of NameChecker().
        :param type: str, optional
        :param name: convert_all - Set to True if you want to convert all necessairy files on initialization. Default: False.
        :param type: bool, optional
        """

        self._path, self._out_path, self._base_path = nc.return_all()
        if convert_all:
            self.pdb_to_xyzrn(self._path, self._out_path)
            self.msms(self._path, density)
            self.pdb_to_mol2(self._path, self._out_path)
            self.pdb_to_pqr(self._path, self._out_path)
            
    
    def pdb_to_xyzrn(self, path, output):
        """
        Converts .pdb to .xyzrn file.
        
        :param name: path - Name of input (pdb) file (without extension)
        :param type: str
        :param name: output - Name of output (xyzrn) file (without extension)
        :param type: str
        
        :return: void - xyzrn file
        """

        # create filepath+filename.xyzrn for output file
        xyzrn = output + ".xyzrn"
        name = path + ".pdb"

        if not os.path.exists(xyzrn):
            print("Converting xyzrn file...")
            output_pdb_as_xyzrn(name, xyzrn)
        else:
            print("xyzrn files already exist. No conversion needed")
        
    def msms(self, path, dens):
        """
        Run the msms binary for given filename.
        
        It takes the {path}.xyzrn file as input and output is written to {path}_out_{dens}
        Binary path is read in from environment variable: MSMS_BIN. If environment variable does not exist binary will be looked up in provis/binaries/msms.
        
        :param name: path - Path to the/Name of the .xyzrn file to be converted.
        :param type: str
        :param name: dens - Density of triangulation
        :param type: float
        
        :return: void - face and vert files
        """
        # print("./msms.exe -if %4s.xyzr -of %4s_out_%s -density %s" % (name, name, str(int(dens)), str(dens)))

        # Now run MSMS on xyzrn file       
        try:
            MSMS_BIN = os.environ['MSMS_BIN']
        except:
            pass
        
        MSMS_BIN = ""
        if not os.path.exists(MSMS_BIN):
            MSMS_BIN = self._base_path + 'binaries/msms'

        file_base = f"{path}_out_{int(dens * 10)}"

        FNULL = open(os.devnull, 'w')
        args = [MSMS_BIN, "-density", f"{dens}", "-hdensity", "3.0", "-probe",\
                        "1.5", "-if", f"{path}.xyzrn", "-of", file_base, "-af", file_base]
        
        output_exists = os.path.exists(file_base + ".vert") and os.path.exists(file_base + ".face") + os.path.exists(file_base + ".area")
        if os.path.isfile(MSMS_BIN):# and not output_exists:
            if not output_exists:
                print("Converting face and vert files...")
                ac = subprocess.call(args)
            else:
                print("face and vert files already exist. No conversion needed")
        else:
            print("MSMS Binary not found under: ", MSMS_BIN)
            
    def pdb_to_mol2(self, path, outpath):
        """
        Run openbabel, to convert pdb to mol2
        
        :param name: path - Name of input (pdb) file (without extension)
        :param type: str        
        :param name: outpath - Name of desired output file (without extension). It will add .mol2 to the given path.
        :param type: str
        """
        
        if not os.path.exists(outpath + ".mol2"):
            print("Converting mol2 file...")
            ac = subprocess.call("obabel %s.pdb -O %s.mol2" % (path,  outpath), shell=True)
        else:
            print("mol2 file already exist. No conversion needed")

    def pdb_to_pqr(self, path, outpath, forcefield="swanson"):
        """
        Run pdb2pqr, to convert pdb to pqr
        
        Binary path is read in from environment variable: PDB2PQR_BIN. If environment variable does not exist binary will be looked up in binaries/pdb2pqr/pdb2pqr.
        
        :param name: path - Name of input (pdb) file (without extension)
        :param type: str        
        :param name: outpath - Name of desired output file (without extension). It will add .pqr to the given path.
        :param type: str
        :param name: forcefield - Force field used for charge computation, by binary. Default: swanson. Options: amber, charmm, parse, tyl06, peoepb and swanson
        :param type: str, optional
        """
        
        try:
            PDB2PQR_BIN  = os.environ['PDB2PQR_BIN']
        except:
            pass
            
        PDB2PQR_BIN = ""
        if not os.path.exists(PDB2PQR_BIN):
            PDB2PQR_BIN = self._base_path + "binaries/pdb2pqr/pdb2pqr"
            
        #possibilities: amber, charmm, parse, tyl06, peoepb and swanson
        args = [PDB2PQR_BIN , f"--ff={forcefield}", f"{path}.pdb", f"{outpath}.pqr"]

        output_exists = os.path.exists(f"{outpath}.pqr")        
        
        if os.path.isfile(PDB2PQR_BIN):# and not output_exists:
            if not output_exists:
                print("Converting pqr file...")
                ac = subprocess.call(args) 
            else:
                print("pqr file already exist. No conversion needed")
        else:
            print("PDB2PQR Binary not found under: ", PDB2PQR_BIN)
       

    def cleanup(self, delete_img: bool=False):
        """
        Deletes all files from data/tmp (and data/img) directories.

        :param name: delete_img - If True data/img directory also cleared, else just the data/tmp. Default: False.
        :param type: bool, optional
        """

        import os, shutil
        tmp = self._base_path + "data/tmp/"
        folders = [tmp]
        if delete_img:
            img = self._base_path + "data/img/"
            folders.append(img)
        for folder in folders:
            for filename in os.listdir(folder):
                file_path = os.path.join(folder, filename)
                try:
                    if os.path.isfile(file_path) or os.path.islink(file_path):
                        os.unlink(file_path)
                    elif os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                except Exception as e:
                    print('Failed to delete %s. Reason: %s' % (file_path, e))

