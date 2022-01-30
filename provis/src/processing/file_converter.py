from genericpath import exists
import subprocess
import os.path
from provis.src.processing.name_checker import NameChecker
from provis.utils.surface_utils import output_pdb_as_xyzrn

class FileConverter():
    """
    Class to create and destroy necessary files required in other parts of code.
    The member functions call the binaries and scripts to convert the files.
    
    Also has a cleanup function that removes everything in the "root_directory"/data/tmp (and data/img if specified) directory. Best practice is to call this function at the end of your main file.
    However if you want to plot the same protein many times, then it is benificial to keep the temporary (data/tmp) files as if they exist provis will not recompute them.
    """
    
    def __init__(self, nc, density=3.0, convert_all=False):
        """
        If "convert_all" set to True conversions instant. Creates xyzr and mol2 files in every case and pqr, face and vert files if appropriate binary exists.
        
        Parameters:
            nc: str, optional
                Instance of a NameChecker class. Used to pass the pdb file name and paths.
            density: float, optional
                Density of triangles. Default: 3.0.
            base_path: str, optional
                Path to "working directory" according to the rules of NameChecker().
            convert_all: bool, optional
                Set to True if you want to convert all necessairy files on initialization. Default: False.
        """

        self._path, self._out_path, self._base_path, mesh = nc.return_all()
        if convert_all:
            self.pdb_to_xyzrn(self._path, self._out_path)
            self.msms(self._path, density)
            self.pdb_to_mol2(self._path, self._out_path)
            self.pdb_to_pqr(self._path, self._out_path)
            self.decompose_traj(self._path)
            
    
    def pdb_to_xyzrn(self, path, output):
        """
        Converts .pdb to .xyzrn file.
        
        Parameters:
            path: str
                Name of input (pdb) file (without extension)
            output: str
                Name of output (xyzrn) file (without extension)
        
        Returns: 
            void
                xyzrn file
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
        
        Parameters:
            path: str
                Path to the/Name of the .xyzrn file to be converted.
            dens: float
                Density of triangulation
        
        Returns: 
            void
                face and vert files
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
        print(args)
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
        
        Parameters:
            path: str        
                Name of input (pdb) file (without extension)
            outpath: str
                Name of desired output file (without extension). It will add .mol2 to the given path.
               
        Returns: 
            void
                mol2 file
        """
        
        if not os.path.exists(outpath + ".mol2"):
            print("Converting mol2 file...")
            ac = subprocess.call("obabel %s.pdb -O %s.mol2" % (path,  outpath), shell=True)
        else:
            print("mol2 file already exist. No conversion needed")

    def pdb_to_pqr(self, path, outpath, forcefield="swanson"):
        """
        Run the pdb2pqr binary for given filename.
        
        It takes the {path}.pdb file as input and output is written to {outpath}.pqr.
        Binary path is read in from environment variable: MSMS_BIN. If environment variable does not exist binary will be looked up in provis/binaries/msms.

        
        Binary path is read in from environment variable: PDB2PQR_BIN. If environment variable does not exist binary will be looked up in binaries/pdb2pqr/pdb2pqr.
        
        Parameters:
            path: str        
                Name of input (pdb) file (without extension)
            outpath: str
                Name of desired output file (without extension). It will add .pqr to the given path.
            forcefield: str, optional
                Force field used for charge computation, by binary. Default: swanson. Options: amber, charmm, parse, tyl06, peoepb and swanson
              
        Returns: 
            void
                pqr file
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
       
    def decompose_traj(self, path):
        """
        As the MSMS binary is unable to work with trajectory .pdb files the large .pdb file containing all models needs to be decomposed to one file per model.
        This function completes this exact task.

        Only executes decomposition if the first file ("{pdb_file_name}_0.pdb") does not exists to avoid unnecessairy recomputation.
        
        Parameters:
            path: str   
                Name of input (pdb) file (without extension)

        Returns: 
            int
                Number of models in the trajectory.
        """
        print("Decomposing trajectory pdb file")
        # Writing to a file
        traj_name = path + ".pdb"
        traj = open(traj_name, 'r')
        model_id = 0
        
        model_0 = path + "_0.pdb"
        if os.path.exists(model_0):
            from Bio.PDB import PDBParser
            parser = PDBParser()
            file_name = path + ".pdb"
            structure = parser.get_structure(path, file_name)
            i = 0
            for model in structure:
                i += 1

            return i
        
        while True:
            new_name = path + "_" + str(model_id) + ".pdb"

            new_file = open(new_name, 'w')
        
            # Get next line from file
            line = traj.readline()
        
            # if line is empty
            # end of file is reached
            if not line:
                new_file.close()
                os.unlink(new_name)
                break
            
            while line != 'TER\n':
                while line == "ENDMDL\n" or line[0:6] == "REMARK":
                    line = traj.readline()
                new_file.writelines(line)
                line = traj.readline()
                if not line:
                    break
                
            new_file.close()
            model_id += 1
            
        traj.close()
        
        return model_id
        
    def cleanup(self, delete_img: bool=False, delete_meshes: bool=False):
        """
        Deletes all files related to the current .pdb id from data/tmp (and if specified, the data/img and data/meshes) directories.

        CAUTION: provis does not recompute existing files. So if you have a molecule that you want to plot multiple times then do not delete the temporary files.
        WARNING: as provis does not recompute existing files it might occur that an old version of the file is stored in the temporary directories and this might cause provis to fail. If this is the case simply delete all temporary files (as well as meshes).

        Parameters:
            delete_img: bool, optional
                 If True all files of the form {pdb_id}_{*} will be deleted from the data/meshes directory. (pdb_id is the name of the .pdb file without the .pdb extension and {*} represents that "anything"). Default: False.
            delete_meshes: bool, optional
                If True all files of the form {pdb_id}_{*} will be deleted from the data/meshes directory. (pdb_id is the name of the .pdb file without the .pdb extension and {*} represents that "anything"). Default: False
        """

        import os, shutil

        pdb_id = self._path.split("/")[-1].split("_")[0]
        tmp = self._base_path + "data/tmp/"
        folders = [tmp]
        if delete_img:
            img = self._base_path + "data/img/"
            folders.append(img)
        for folder in folders:
            for filename in os.listdir(folder):
                if filename.split(".")[0].split("_")[0] == pdb_id:
                    file_path = os.path.join(folder, filename)
                    try:
                        if os.path.isfile(file_path) or os.path.islink(file_path):
                            os.unlink(file_path)
                        elif os.path.isdir(file_path):
                            shutil.rmtree(file_path)
                    except Exception as e:
                        print('Failed to delete %s. Reason: %s' % (file_path, e))

