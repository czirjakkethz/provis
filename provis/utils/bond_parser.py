import numpy as np
import pandas as pd

def bond_parser(filename):
    """
    Calculates bond information from mol2 file.
    
    :param name: filename - Name of file to open
    :param type: str
    
    :returns: pandas.dataFrame - Pandas df with bonds
    """
    f = open(filename,'r')
    f_text = f.read()
    f.close()
    bond_start = f_text.find('@<TRIPOS>BOND')
    bond_end = f_text[bond_start:].replace('@<TRIPOS>BOND','').find('@')
    l = int(len(f_text[bond_start:bond_end].replace('@<TRIPOS>BOND\n','').replace('\n',' ').strip().split())/4)
    df_bonds = pd.DataFrame(np.array(f_text[bond_start:bond_end].replace('@<TRIPOS>BOND\n','').replace('\n',' ').strip().split()).reshape((l,4)), #
            columns=['bond_id', 'atom1', 'atom2', 'bond_type'])
    df_bonds.set_index(['bond_id'], inplace=True)
    return df_bonds

