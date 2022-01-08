"""
This file was created by Kristof Czirjak as part of his Bachelor's Thesis - provis
"""
from provis.utils.name_checker import NameChecker
from provis.src.processing.residue import Residue
from provis.src.plotting.surface import Surface

def main():
    """
    This is an example file to showcase some functionalities of provis.
    
    We showcase how to run provis when you are outside the special directroy structure specified in the setup section of the documentation.
    The full path to the pdb file has to be passed, as well as the full path to the above mentioned special directory structure.
        
    First:
    Define variables needed later:
    """
    name = "/home/kczi/Documents/provis/data/pdb/selected_prediction"
    base_path='/home/kczi/Documents/provis/'
    density = 3.0
    solvent = 0
    bash = 0
    
    """
    Second:
    Initializing the NameChecker class will store the name of the pdb file and output path.
    This class is usually initialized when by the FileConverter class, but if you have the necessairy files you do not need to call the FileConverter class at all.
    """
    nc = NameChecker(name, base_path)
    
    """
    Third:
    Surface is a class that handles all plotting related to surfaces. This includes basic surface visualization as well as the surface properties.
    """
    s = Surface(dens=density)
    s.plot_surface()    
    s.plot_hydrophob()
    s.plot_shape()
    s.plot_charge()


if __name__ == "__main__":
    main()

