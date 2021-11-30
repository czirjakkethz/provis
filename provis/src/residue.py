
class Residue():
    """
    Residue class is mainly used for plotting a bounding box around the given residue.
    """
    def __init__(self, id=None, chain=0, padding=0):
        """
        Initialize class. Can be empty, but better initialized with residue id and chain id.
        
        :param name: id - Residue id. Count starting at 0.
        :param type: int
        :param name: chain - Chain id. Count starting at 0. Specify which chain the residue is on. Default 0 (in case of single chain)
        :param type: int
        :param name: padding - Optional padding to bounding box.
        :param type: int
        """
        if id == None:
            self._id_list = []
            self._chain_list = []
        else:
            self._id_list = [id]
            self._chain_list = [chain]
        self._pad = padding
        
    def add_residue(self, id, chain=0):
        """
        Add a new residue to the internal list
        
        :param name: id - id of residue
        :param type: int
        :param name: chain - chain id of residue
        :param type: int
        
        :return: list - list of the current residues
        """
        self._id_list.append(id)
        self._chain_list.append(chain)
    
    def get_res_info(self):
        """
        Returns all internal information of class
        
        :return: list - list of the current residues
        :return: int - padding for bounding box
        """
        return self._id_list, self._chain_list, self._pad
    
