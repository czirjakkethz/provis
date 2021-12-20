import subprocess
import os.path
from provis.utils.name_checker import check_name
from provis.utils.surface_utils import output_pdb_as_xyzrn

class FileConverter():
    """
    Class to create and destroy necessary files required in other parts of code.
    It has a bunch of staticmethod functions to call the binaries and scripts to convert the files.
    """
    def __init__(self, name=None, dens=None, solv=0, bash=0):
        """
        Can be constructed empty. If called with arguments conversions instant. Creates xyzr and mol2 files in every case and face and vert files if msms binary exists.

        :param name: name - name of file. With or without .pdb extension. Assumes it to be in one of the following directories: /provis or /provis/data/pdb
        :param type: str
        :param name: dens - density of triangles
        :param type: float
        :param name: solv - set to True if you want to plot solvent atoms as well. Default: false
        :param type: bool
        :param name: bash - set to True if you want to run binary version of pdb_to_xyzr (comes with msms). Default: false
        :param type: bool
        """

        self._solv = solv
        self._bash = bash
        if name:
            self._path, self._out_path = check_name(name)
            self.pdb_to_mol2(self._path, self._out_path)
            self.pdb_to_pqr(self._path, self._out_path)
            self.pdb_to_xyzrn(self._path, self._out_path)
        if dens:
            self._dens = dens
        if name and dens:
            self.msms(self._out_path, self._dens)
        pass
    
    @staticmethod
    def pdb_to_xyzrn(path, output):
        """
        Converts .pdb to .xyzrn file.
        
        :param name: path - Name of input (pdb) file (without extension)
        :param type: str
        :param name: output - Name of output (xyzrn) file (without extension)
        
        :return: void - xyzrn file
        """
        print("Converting xyzrn file...")

        # create filepath+filename.xyzrn for output file
        xyzrn = output + ".xyzrn"
        name = path + ".pdb"
        
        output_pdb_as_xyzrn(name, xyzrn)
        
    @staticmethod
    def msms(path, dens):
        """
        Run the msms binary for given filename.
        
        It takes the {path}.xyzrn file as input and output is written to {path}_out_{dens}
        
        :param name: path - Path to the/Name of the .xyzrn file to be converted.
        :param type: str
        :param name: dens - Density of triangulation
        :param type: float
        
        :return: void - face and vert files
        """
        print("Converting face and vert files...")
        # print("./msms.exe -if %4s.xyzr -of %4s_out_%s -density %s" % (name, name, str(int(dens)), str(dens)))

        # Now run MSMS on xyzrn file
        MSMS_BIN = os.environ['MSMS_BIN']

        file_base = f"{path}_out_{int(dens * 10)}"

        FNULL = open(os.devnull, 'w')
        args = [MSMS_BIN, "-density", f"{dens}", "-hdensity", "3.0", "-probe",\
                        "1.5", "-if", f"{path}.xyzrn", "-of", file_base, "-af", file_base]
        if os.path.isfile(MSMS_BIN):
            ac = subprocess.call(args) 
        else:
            print("MSMS Binary not found under: ", MSMS_BIN)
            
    @staticmethod
    def pdb_to_mol2(path, outpath):
        """
        Run openbabel, to convert pdb to mol2
        
        :param name: path - Name of input (pdb) file (without extension)
        :param type: str        
        :param name: outpath - Name of desired output file (without extension). It will add _out.pqr to the given path.
        :param type: str
        """
        print("Converting mol2 file...")
        ac = subprocess.call("obabel %4s.pdb -O %4s.mol2" % (path,  outpath), shell=True)

    @staticmethod
    def pdb_to_pqr(path, outpath):
        """
        Run pdb2pqr, to convert pdb to pqr
        
        :param name: path - Name of input (pdb) file (without extension)
        :param type: str        
        :param name: outpath - Name of desired output file (without extension). It will add _out.pqr to the given path.
        :param type: str
        """
        print("Converting pqr file...")

        PDB2PQR_BIN  = os.environ['PDB2PQR_BIN']

        args = [PDB2PQR_BIN , "--clean", f"{path}.pdb", f"{outpath}_out.pqr"]

        if os.path.isfile(PDB2PQR_BIN ):
            ac = subprocess.call(args)
        else:
            print("MSMS Binary not found under: ", PDB2PQR_BIN )

    @staticmethod
    def cleanup(delete_img: bool=False):
        """
        Deletes all files from data/tmp (and data/img) directories.

        :param name: delete_img - If True data/img directory also cleared, else just the data/tmp
        :param type: bool
        """

        import os, shutil
        folders = ["data/tmp/"]
        if delete_img:
            folders.append("data/img/")
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
        # files_to_remove = [
        #     xyzrn_file, file_base + ".vert", file_base + ".area",
        #     file_base + ".face"
        # ]

        # for file in files_to_remove:
        #     if os.path.exists(file):
        #         os.remove(file)
