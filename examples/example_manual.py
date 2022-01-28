"""
This file was created by Kristof Czirjak as part of his Bachelor's Thesis - provis
"""
from ctypes import Structure
from provis.src.processing.data_handler import DataHandler
from provis.src.processing.name_checker import NameChecker
from provis.src.processing.residue import Residue
from provis.src.plotting.surface import Surface

def main():
    """
    This is an example file to showcase some functionalities of provis.
    
    We showcase how to run provis when you are outside the special directroy structure specified in the setup section of the documentation.
    The full path to the pdb file has to be passed, as well as the full path to the above mentioned special directory structure.
        
    Define variables needed later:
    """
    name = "/home/kczi/Documents/provis/data/pdb/2fd7"
    base_path='/home/kczi/Documents/provis/'
    density = 3.0
    
    """
    The NameChecker class will find and store the path of the pdb file and output path.
    """
    nc = NameChecker(name, base_path)
    
    """
    A DataHandler class can be initialzied explicitly, by passing the NameChecker class instance.
    However, this is actually not obligatory, as the higher level plotting classes can create their own DataHandler classes.
    It is still best practice to initialize all classes explicitly and pass them as arguments to other classes to avoid unnecessairy duplication in memory and unnecessairy computations.
    """
    dh = DataHandler(nc)
    
    """
    Structure is a class that handles all plotting related to the structure / general layout of the molecule. This includes the atom cloud, the bonds or the backbone of the molecule, amongst others.    
    
    As mentioned above we can initialize the Structure class with a DataHandler.
    """

    struct = Structure(nc, dh=dh)
    struct.plot_backbone()
    struct.plot_atoms()
    struct.plot_bonds()
    struct.plot_vw()
    struct.plot_stick_point()
    struct.plot_residues()
    r = Residue(29)
    r.add_residue(50)
    r. add_residue(1, 1)
    r.remove_residue(1, 1)
    struct.plot(atoms=1, box=1, bonds=1, vw=0, residues=0, res=r, bb=0)    
    
    """
    Surface is a class that handles all plotting related to surfaces. This includes basic surface visualization as well as the surface properties.    
    
    As seen here the Surface class can be initialized with only a NameChecker instance. Ideally, however one would pass a SurfaceHandler class instance on initialization.
    """
    s = Surface(nc, density=density)
    s.plot()    
    s.plot_hydrophob()
    s.plot_shape()
    s.plot_charge()


if __name__ == "__main__":
    main()

