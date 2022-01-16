"""
This file was created by Kristof Czirjak as part of his Bachelor's Thesis - provis
"""
from provis.src.processing.file_converter import FileConverter
from provis.src.processing.residue import Residue
from provis.src.plotting.structure import Structure
from provis.src.processing.name_checker import NameChecker

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
    
    """
    Second:
    Initialize the FileConverter class. This class creates the necessairy files (xyzrn, pqr, face, vert, etc.) for the plotting. 
    If you initialize the class with variables, like done here, conversion is instant. (You can also call member functions of the class if you only need a specific conversion).
    
    The FileConverter class also initializes the NameChecker class, thus storing the name of the pdb file and output path to be used globally.
    """
    nc = NameChecker(name)
    fc = FileConverter(nc, density, solvent, convert_all=True) 
    
    """
    Third:
    Structure is a class that handles all plotting not related to surfaces. This includes simple stick and point plots, the visualization of the bonds or the backbone of the protein.
    """
    sp = Structure(nc)
    sp.plot_backbone()
    sp.plot_atoms()
    sp.plot_bonds()
    sp.plot_vw()
    sp.plot_stick_point()
    r = Residue(1)
    r.add_residue(3)
    r. add_residue(1, 1)
    r.remove_residue(1, 1)
    sp.plot(atoms=1, box=1, bonds=1, vw=0, residues=0, res=r, bb=0)

    """
    Finally we can clean up all the temporary files from the "root direcotry"/data/tmp (and data/img) folders.
    
    It is worth to mention however, that sometimes it is worth to keep the temporary files if you want to visualize this protein later on.
    If you keep the temporary files they will not have to be computed again and visualization will be quicker.
    """
    fc.cleanup(delete_img=0)

if __name__ == "__main__":
    main()

