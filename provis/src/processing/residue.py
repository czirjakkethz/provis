
class Residue():
    """
    Residue class is mainly used for plotting a bounding box around the given residue.
    """
    def __init__(self, id=None, chain=0, padding=0):
        """
        Initialize class. Can be empty, but better initialized with residue id and chain id.
        
        :param name: id - Residue id. Count starting at 0. Default: None.
        :param type: int, optional
        :param name: chain - Chain id. Count starting at 0. Specify which chain the residue is on. Default 0 (in case of single chain). Default: 0.
        :param type: int, optional
        :param name: padding - Optional padding to bounding box. Default: 0.
        :param type: int, optional
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
        Add a new residue to the internal list of residues.
        
        :param name: id - Residue id. Count starting at 0. Default: None.
        :param type: int, optional
        :param name: chain - Chain id. Count starting at 0. Specify which chain the residue is on. Default 0 (in case of single chain). Default: 0.
        :param type: int, optional
        """
        self._id_list.append(id)
        self._chain_list.append(chain)
    
    def get_res_info(self):
        """
        Returns all internal information of class
        
        :returns: list - list of the current residues
        :returns: list - list of chain ID's corresponding to residues
        :returns: int - padding for bounding box
        """
        return self._id_list, self._chain_list, self._pad
    
    def remove_residue(self, id, chain=0):
        """
        Remove speciefied residue from internal list.

        :param name: id - Residue id. Count starting at 0. Default: None.
        :param type: int, optional
        :param name: chain - Chain id. Count starting at 0. Specify which chain the residue is on. Default 0 (in case of single chain). Default: 0.
        :param type: int, optional
        """
        
        
        curr_idx = self._id_list.index(id)

        len_ = len(self._id_list) 
        while(self._chain_list[curr_idx] != chain and curr_idx < len_):
            # if first instance of residue id does not correspond to the specified chain then search the rest of the list
            curr_idx = curr_idx + 1 + self._id_list[curr_idx + 1:].index(id)
        
        if curr_idx < len_:
            del self._id_list[curr_idx]
            del self._chain_list[curr_idx]
        else: 
            print(f"No such residue found: id {id}, chain: {chain}")
    
