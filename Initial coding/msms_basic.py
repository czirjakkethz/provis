from Bio.PDB import *
from Bio import SeqIO
import numpy as np
import pyvista as pv
import ssbio.utils
from ssbio.protein.structure.utils.structureio import StructureIO

# ssbio.protein.structure.properties.msms.get_msms_df_on_file(pdb_file, outfile=None, outdir=None, outext='_msms.df', force_rerun=False)

# parser = PDBParser() 
# load_pdb = "2fd7.pdb"
# structure = parser.get_structure("2fd7", load_pdb) 

# ssbio.protein.structure.StructProp("2fd7")#, description=None, root_dir=None, pdb_file_type='pdb')

# StructProp(ident, 
# description=None, 
# chains=None,
# mapped_chains=None,
# is_experimental=False,
# structure_path=None,
# file_type=None)
# ssbio.protein.structure - get_residue_depths

from ssbio.databases.pdb import PDBProp
from ssbio.databases.uniprot import UniProtProp

import sys
import logging

# Create logger
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)  # SET YOUR LOGGING LEVEL HERE #

# Other logger stuff for Jupyter notebooks
handler = logging.StreamHandler(sys.stderr)
formatter = logging.Formatter('[%(asctime)s] [%(name)s] %(levelname)s: %(message)s', datefmt="%Y-%m-%d %H:%M")
handler.setFormatter(formatter)
logger.handlers = [handler]

