"""
This file was created by Kristof Czirjak as part of his Bachelor's Thesis - provis
"""
from provis.src.processing.protein import Protein
from provis.src.plotting.dynamic_plotter import DynamicPlotter
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
    base_path = None
    density = 3.0
    msms = True
    plot_solvent=False
    notebook=False
    """
    Second:
    Initialize the Protein class with the desired molecule. 
    
    It is possible to choose a specific model from a trajectory to plot it as a single protein. 
    """

    prot = Protein(name, base_path=None, density=density)

    """
    Third:
    Plot!
    
    Use the DynamicPlotter class to plot the full trajectory of an protein.
    """
    dp = DynamicPlotter(prot, msms=msms, notebook=notebook, plot_solvent=plot_solvent)

    r = Residue(29)
    r.add_residue(50)
    r. add_residue(1, 1)
    r.remove_residue(1, 1)

    dp.plot_backbone()
    dp.plot_atoms(res=r)
    dp.plot_bonds()
    dp.plot_vw()
    dp.plot_stick_point()
    dp.plot_residues()
    dp.plot_structure(atoms=1, box=1, bonds=1, vw=0, residues=0, res=r, bb=0)

    dp.plot_surface()
    dp.plot_hydrophob()
    dp.plot_shape()
    dp.plot_charge()
    
    """
    And finally clean up everything with the "cleanup" function of the prot.file_converter (FileConverter class) member variable.
    """
    # prot.file_converter.cleanup()

if __name__ == "__main__":
    main()

