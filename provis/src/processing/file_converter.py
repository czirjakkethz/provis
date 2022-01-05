import subprocess
import os.path
from provis.utils.name_checker import NameChecker
from provis.utils.surface_utils import output_pdb_as_xyzrn

class FileConverter():
    """
    Class to create and destroy necessary files required in other parts of code.
    It has a bunch of classmethod functions to call the binaries and scripts to convert the files.
    
    Best practice is to initialize a FileConverter in a given file before calling any function of provis and call the cleanup method after the last plotting function is called. This will keep your directories decluttered. 
    However if you want to plot the same protein many times, then it is benificial to keep the temporary (data/tmp) files as if they exist provis will not recompute them.
    """
    
    _base_path = ""
    def __init__(self, name=None, dens=None, solv=0, bash=0, base_path=None):
        """
        Can be constructed empty. If called with arguments conversions instant. Creates xyzr and mol2 files in every case and face and vert files if msms binary exists.

        :param name: name - Name of file to be loaded. Passed to NameChecker() to get usable paths. Default: None.
        :param type: str, optional
        :param type: float, optional
        :param name: solv - set to True if you want to plot solvent atoms as well. Default: False.
        :param type: bool, optional
        :param name: bash - set to True if you want to run binary version of pdb_to_xyzr (comes with msms). Default: False.
        :param type: bool, optional
        :param name: base_path - Path to "working directory" according to the rules of NameChecker().
        :param name: dens - Density of triangles. Default: None.
        """

        self._solv = solv
        self._bash = bash
        if name:
            NameChecker(name, base_path)
            self._path, self._out_path, FileConverter._base_path = NameChecker.return_all()
            self.pdb_to_mol2(self._path, self._out_path)
            self.pdb_to_pqr(self._path, self._out_path)
            self.pdb_to_xyzrn(self._path, self._out_path)
        if dens:
            self._dens = dens
        if name and dens:
            self.msms(self._out_path, self._dens)
        pass
    
    @classmethod
    def pdb_to_xyzrn(cls, path, output):
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
        
    @classmethod
    def msms(cls, path, dens):
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
        
        if not os.path.exists(MSMS_BIN):
            MSMS_BIN = FileConverter._base_path + 'binaries/msms'

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
            
    @classmethod
    def pdb_to_mol2(cls, path, outpath):
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

    @classmethod
    def pdb_to_pqr(cls, path, outpath, forcefield="swanson"):
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
            
        if not os.path.exists(PDB2PQR_BIN):
            PDB2PQR_BIN = FileConverter._base_path + "binaries/pdb2pqr/pdb2pqr"
            
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
       

    @classmethod
    def cleanup(cls, delete_img: bool=False):
        """
        Deletes all files from data/tmp (and data/img) directories.

        :param name: delete_img - If True data/img directory also cleared, else just the data/tmp. Default: False.
        :param type: bool, optional
        """

        import os, shutil
        tmp = FileConverter._base_path + "data/tmp/"
        folders = [tmp]
        if delete_img:
            img = FileConverter._base_path + "data/img/"
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

