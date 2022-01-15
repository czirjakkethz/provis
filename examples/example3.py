"""
This file was created by Kristof Czirjak as part of his Bachelor's Thesis - provis
"""
from provis.src.processing.protein import Protein
from provis.src.processing.residue import Residue

def main():
    """
    This is an example file to showcase some functionalities of provis.
    
    We showcase the easiest way to run provis. For this you should have this file in the root directory of the special directroy structure specified in the setup section of the documentation.
    The path to the "root directory"/data/tmp will automatically be found.
    This way you can have your pdb files nicely organized in the data/pdb directory (or simply have them in the root directory). 
    Your temporary files will be in the data/tmp directory and the screenshots of the plots in the data/img directory.
    
    First:
    Define variables needed later:
    """
    name = "2fd7"
    density = 3.0
    solvent = 0
    
    prot = Protein(name, density=density, plot_solvent=solvent)
    
    # prot.structure.plot_backbone()
    # prot.structure.plot_atoms()
    # prot.structure.plot_bonds()
    # prot.structure.plot_vw()
    # prot.structure.plot_stick_point()
    # r = Residue(1)
    # r.add_residue(3)
    # r. add_residue(1, 1)
    # r.remove_residue(1, 1)
    # prot.structure.plot(atoms=1, box=1, bonds=1, vw=0, residues=0, res=r, bb=0)
    
    prot.surface.plot_surface()    
    prot.surface.plot_hydrophob()
    prot.surface.plot_shape()
    prot.surface.plot_charge()

    prot._file_converter.cleanup()

if __name__ == "__main__":
    main()

