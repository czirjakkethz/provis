
class Residue():
    """
    Residue class is used for plotting a bounding box around the specified residue.
    """
    def __init__(self, id=None, chain=0, padding=0):
        """
        Initialize class. Can be empty, but better initialized with residue id and chain id.

        Chain specifies which chain the residue is on and id is the id of the residue on the given chain.
        
        Parameters:
            id: int, optional
                Residue id. Count starting at 0. Default: None.
            chain: int, optional
                Chain id. Count starting at 0. Specify which chain the residue is on. Default 0 (in case of single chain). Default: 0.
            padding: int, optional
                Optional padding to bounding box. Default: 0.
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
        
        Parameters:
            id: int, optional
                Residue id. Count starting at 0. Default: None.
            chain: int, optional
                Chain id. Count starting at 0. Specify which chain the residue is on. Default 0 (in case of single chain). Default: 0.
        """
        self._id_list.append(id)
        self._chain_list.append(chain)
    
    def get_res_info(self):
        """
        Returns all internal information of class
        
        Returns: 
            list:
                list of the current residues
            list:
                list of chain ID's corresponding to residues
            int:
                padding for bounding box
        """
        return self._id_list, self._chain_list, self._pad
    
    def remove_residue(self, id, chain=0):
        """
        Remove speciefied residue from internal list.

        Parameters:
            id: int, optional
                Residue id. Count starting at 0. Default: None.
            chain: int, optional
                Chain id. Count starting at 0. Specify which chain the residue is on. Default 0 (in case of single chain). Default: 0.
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
    
