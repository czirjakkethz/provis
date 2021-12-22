
import os
from Bio.PDB import PDBParser, Selection

parser = PDBParser()

#Kyte-Doolittle scale for hydrophobicity
KD_SCALE = {}
KD_SCALE["ILE"] = 4.5
KD_SCALE["VAL"] = 4.2
KD_SCALE["LEU"] = 3.8
KD_SCALE["PHE"] = 2.8
KD_SCALE["CYS"] = 2.5
KD_SCALE["MET"] = 1.9
KD_SCALE["ALA"] = 1.8
KD_SCALE["GLY"] = -0.4
KD_SCALE["THR"] = -0.7
KD_SCALE["SER"] = -0.8
KD_SCALE["TRP"] = -0.9
KD_SCALE["TYR"] = -1.3
KD_SCALE["PRO"] = -1.6
KD_SCALE["HIS"] = -3.2
KD_SCALE["GLU"] = -3.5
KD_SCALE["GLN"] = -3.5
KD_SCALE["ASP"] = -3.5
KD_SCALE["ASN"] = -3.5
KD_SCALE["LYS"] = -3.9
KD_SCALE["ARG"] = -4.5

def get_residues(pdb_file):
    global parser
    pdb_file = os.path.abspath(pdb_file)
    base_file = pdb_file.split("/")[-1]  # Remove any full path prefixes
    pdb_id = base_file.split(".")[0]

    struct = parser.get_structure(file=pdb_file, id=pdb_id)
    residues = Selection.unfold_entities(struct, "R")
    # Remove heteroatoms
    residues = [res for res in residues if res.get_full_id()[3][0] == " "]
    return residues

