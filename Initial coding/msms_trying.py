import sys
import logging

# Import the Protein class
from ssbio.core.protein import Protein

# Printing multiple outputs per cell
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"

# Create logger
logger = logging.getLogger()
logger.setLevel(logging.INFO)  # SET YOUR LOGGING LEVEL HERE #

# Other logger stuff for Jupyter notebooks
handler = logging.StreamHandler(sys.stderr)
formatter = logging.Formatter('[%(asctime)s] [%(name)s] %(levelname)s: %(message)s', datefmt="%Y-%m-%d %H:%M")
handler.setFormatter(formatter)
logger.handlers = [handler]

# SET FOLDERS AND DATA HERE
import tempfile
ROOT_DIR = tempfile.gettempdir()

PROTEIN_ID = 'SRR1753782_00918'
PROTEIN_SEQ = 'MSKQQIGVVGMAVMGRNLALNIESRGYTVSVFNRSREKTEEVIAENPGKKLVPYYTVKEFVESLETPRRILLMVKAGAGTDAAIDSLKPYLEKGDIIIDGGNTFFQDTIRRNRELSAEGFNFIGTGVSGGEEGALKGPSIMPGGQKDAYELVAPILTKIAAVAEDGEPCVTYIGADGAGHYVKMVHNGIEYGDMQLIAEAYSLLKGGLNLSNEELANTFTEWNNGELSSYLIDITKDIFTKKDEDGNYLVDVILDEAANKGTGKWTSQSALDLGEPLSLITESVFARYISSLKAQRVAASKVLSGPKAQPAGDKAEFIEKVRRALYLGKIVSYAQGFSQLRAASDEYHWDLNYGEIAKIFRAGCIIRAQFLQKITDAYAENADIANLLLAPYFKKIADEYQQALRDVVAYAVQNGIPVPTFSAAVAYYDSYRAAVLPANLIQAQRDYFGAHTYKRTDKEGIFHTEWLE'

# Create the Protein object
my_protein = Protein(ident='my_prot', root_dir='2fd7', pdb_file_type='pdb')

# Load the protein sequence
# This sets the loaded sequence as the representative one
my_protein.load_manual_sequence(seq=PROTEIN_SEQ, ident='WT', write_fasta_file=True, set_as_representative=True)

# Mapping using BLAST
my_protein.blast_representative_sequence_to_pdb(seq_ident_cutoff=0.9, evalue=0.00001)
my_protein.df_pdb_blast.head()


# from Bio.PDB import *
# from Bio import SeqIO
# import numpy as np
# import pyvista as pv

# import ssbio.utils
# from ssbio.protein.structure.utils.structureio import StructureIO

# ssbio.protein.structure.properties.msms.get_msms_df_on_file("2fd7.pdb", outfile="out", outdir=None, outext='_msms.df', force_rerun=False)
