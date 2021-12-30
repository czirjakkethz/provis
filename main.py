"""
This file was created by Kristof Czirjak as part of his Bachelor's Thesis - provis
"""
from provis.src.processing.file_converter import FileConverter
from provis.src.processing.residue import Residue
from provis.src.plotting.stick_point import StickPoint
from provis.src.plotting.surface import Surface

# TODO: RUN Sphinx and push                                    
def main():
    """
    This is an example file to showcase all functionalities of provis.
    
    Provis has two main functionalities. Plotting the atoms and bonds as a stick and point model as well as plotting the surface information of the protein.
    """
    
    # define variables needed later
    name = "2fd7" # "1a3n" # "data/pdb/2fd7" # "data/pdb/1a3n" # "data/pdb/7nkd" #
    density = 3.0
    solvent = 0
    bash = 0
    
    # Create all the necessary files from the .pdb file
    fc = FileConverter(name, density, solvent, bash)
    
    # Plot stick point
    sp = StickPoint()
    sp.plot_atoms()
    sp.plot_bonds()
    sp.plot_vw()
    sp.plot_stick_point()
    r = Residue(1)
    r.add_residue(3)
    r. add_residue(1, 1)
    r.remove_residue(1, 1)
    sp.plot(atoms=1, box=1, bonds=1, vw=0, residues=0, res=r)


    # Plot surface
    s = Surface(dens=density)
    s.plot_hydrophob()
    s.plot_surface(msms=False)    
    s.plot_shape()
    s.plot_charge()

    # Clean up directories
    fc.cleanup(delete_img=0)

if __name__ == "__main__":
    main()

