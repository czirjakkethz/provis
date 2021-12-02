import subprocess

from provis.utils.msms_caller import msms_script
from provis.utils.pdb_to_xyzr_caller import pdb_to_xyzr_script

class FileConverter():
    """
    Class to create necessary files required in other parts of code.
    """
    def __init__(self, name=None, dens=None, solv=0, bash=None):
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
    
        if name:
            self._name = name
        if dens:
            self._dens = dens
        self._solv = solv
        if bash:
            self._bash = bash
        if name and bash:
            self.pdb_to_xyzr(self._name, self._solv, self._bash)
        if name and dens:
            self.msms(self._name, self._dens)
        if name:
            self.pdb_to_mol2(self._name)
        pass
    
    @staticmethod
    def pdb_to_xyzr(name, solvent, bash=0):
        """
        Run the pdb_to_xyzr script for given filename
        
        :param name: name - Name of file
        :param type: str
        :param name: bash - set to True if you want to run binary version of pdb_to_xyzr (comes with msms). Default: false
        :param type: bool
        
        :return: void - xyzr file
        """
        pdb_to_xyzr_script(name, solvent, bash)

    @staticmethod
    def msms(name, dens):
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
        msms_script(name, dens)

    @staticmethod
    def pdb_to_mol2(name):
        """
        Run openbabel, to convert pdb to mol2
        
        :param name: name - name of file (without extension)
        :param type: str
        """
        ac = subprocess.call("obabel %4s.pdb -O %4s.mol2" % (name, name), shell=True)


