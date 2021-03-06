"""
This file was created by Kristof Czirjak as part of his Bachelor's Thesis - provis
"""
from provis.src.plotting.plotter import Plotter
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
    name = "traj"
    density = 3.0
    msms = True
    notebook = False
    """
    Second:
    Initialize the Protein class with the desired molecule. 
    
    It is possible to choose a specific model from a trajectory to plot it as a single protein. 
    """
    
    prot = Protein(name, base_path=None, density=density)
    prot2 = Protein(name, base_path=None, density=density, model_id=30)

    """
    Third:
    Plot!
    
    Use the Plotter class to plot. At least one Protein instance has to be passed. If two proteins are passed then they will be plotted side-by-side.
    It is also possible to add more proteins using the Plotter.add_protein(Protein) method.
    """
    plot = Plotter(prot, prot2, msms=msms, notebook=notebook)
    
    plot.plot_backbone()
    plot.plot_atoms()
    plot.plot_bonds()
    plot.plot_vw()
    plot.plot_stick_point()
    plot.plot_residues()
    r = Residue(29)
    r.add_residue(50)
    r. add_residue(1, 1)
    r.remove_residue(1, 1)
    plot.plot_structure(atoms=1, box=1, bonds=1, vw=0, residues=0, res=r, bb=0)

    plot.plot_surface(res=r)
    plot.plot_hydrophob(res=r)
    plot.plot_shape(res=r)
    plot.plot_charge(res=r)

    """
    And finally clean up everything with the "cleanup" function of the prot.file_converter (FileConverter class) member variable.    

    CAUTION: it is not advised to call the cleanup function if you plan to plot the same exact molecule once agian, since it is more efficient to load information from a temporary file than to recalculate it.
    """
    # prot.file_converter.cleanup()

if __name__ == "__main__":
    main()

