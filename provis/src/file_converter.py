import subprocess

class FileConverter():
    """
    Class to create necessary files required in other parts of code.
    It has a bunch of staticmethod functions to call the binaries and scripts to convert the files.
    """
    def __init__(self, name=None, dens=None, solv=0, bash=0):
        """
        Can be constructed empty. If called with arguments conversions instant. Creates xyzr and mol2 files in every case and face and vert files if msms binary exists.

        :param name: name - name of file
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
            self._name = name
            # self.pdb_to_mol2(self._name)
            # self.pdb_to_pqr(self._name)
            # self.pdb_to_xyzr(self._name, self._solv, self._bash)
        if dens:
            self._dens = dens
        # if name and dens:
            # self.msms(self._name, self._dens)
        pass
    
    @staticmethod
    def pdb_to_xyzr(name, solv, bash):
        """
        Run the pdb_to_xyzr script for given filename.
        Works on all operating systems, as it is a python script.
        (could easily be extended with other parameters)
        
        :param name: name - Name of file
        :param type: str
        :param name: solv - If True also convert solvent atoms, False only structural atoms
        :param type: bool
        :param name: bash - If you want to run original script, set to 1
        :param type: bool
        
        :return: void - xyzr file
        """
        print("Converting xyzr file...")
        
        
        # This is calls the pdb_to_xyzr.py script, that is based on the original msms bash script
        # # print("python pdb_to_xyz.py %1s %4s.pdb > %4s.xyzr" % (str(int(solv)), name, name))
        # if not bash:
        #     ac = subprocess.call("python scripts/pdb_to_xyzr.py %1s %4s.pdb > %4s.xyzr" % (str(int(solv)), name, name), shell=True)
        # else:
        #     ac = subprocess.call("./binaries/pdb_to_xyzr %4s.pdb > %4s.xyzr" % (name, name), shell=True)

    @staticmethod
    def msms(name, dens):
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
        import platform
        import os.path
        print("Converting face and vert files...")
        # print("./msms.exe -if %4s.xyzr -of %4s_out_%s -density %s" % (name, name, str(int(dens)), str(dens)))

        # Now run MSMS on xyzrn file
        MSMS_BIN = os.environ['MSMS_BIN']
        from subprocess import Popen, PIPE

        FNULL = open(os.devnull, 'w')
        args = [MSMS_BIN, "-density", f"{dens}", "-hdensity", "3.0", "-probe",\
                        "1.5", "-if", f"{name}.xyzr", "-of", f"{name}_out_{int(dens)}", "-af", f"{name}_out_{int(dens)}"]
        #print msms_bin+" "+`args`
        # p2 = Popen(args, stdout=PIPE, stderr=PIPE)
        # stdout, stderr = p2.communicate()
        if os.path.isfile(".binaries/msms.exe") or os.path.isfile("binaries/msms"):
            ac = subprocess.call(args) # TODO clean this up
            # if platform.system() == "Windows":
            #     ac = subprocess.call("binaries/msms.exe -if %4s.xyzr -of %4s_out_%s -density %s" % (name, name, str(int(dens)), str(dens)), shell=True)
            # else:
            #     ac = subprocess.call("./binaries/msms -if %4s.xyzr -of %4s_out_%s -density %s" % (name, name, str(int(dens)), str(dens)), shell=True)

    @staticmethod
    def pdb_to_mol2(name):
        """
        Run openbabel, to convert pdb to mol2
        
        :param name: name - name of file (without extension)
        :param type: str
        """
        print("Converting mol2 file...")
        ac = subprocess.call("obabel %4s.pdb -O %4s.mol2" % (name, name), shell=True)

    @staticmethod
    def pdb_to_pqr(name):
        """
        Run pdb2pqr, to convert pdb to pqr
        
        :param name: name - name of file (without extension)
        :param type: str
        """
        print("Converting pqr file...")
        ac = subprocess.call("binaries/pdb2pqr --clean %4s.pdb %4s_out.pqr" % (name, name), shell=True)

