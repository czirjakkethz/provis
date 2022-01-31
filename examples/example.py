"""
This file was created by Kristof Czirjak as part of his Bachelor's Thesis - provis
"""
from provis.src.plotting.dynamic_plotter import DynamicPlotter
from provis.src.plotting.plotter import Plotter
from provis.src.processing.protein import Protein

def main():
    """
    This is an example to showcase briefly showcase the functionalities of provis through a few examples.
    
    This example showcases the easiest way to run provis. For this you should have this file in the root directory of the special directroy structure 
    specified in the setup section of the documentation.
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
    Define a Protein class. This class instance is to be passed to the Plotter classes.
    Initialize with pdb name
    """
    prot = Protein(name, base_path=None, density=density, plot_solvent=False, msms=True, notebook=False)

    """
    Initializing the Plotter and DynamicPlotter classes.
    Use the methods of the classes to plot a single protein or a dynamic trajectory of the protein.   
    """
    
    pl = Plotter(prot, msms=msms, notebook=notebook, plot_solvent=plot_solvent)
    pl.plot_charge()
    
    dp = DynamicPlotter(prot, msms=msms, notebook=notebook, plot_solvent=plot_solvent)
    dp.plot_atoms()

    """
    And finally clean up everything with the "cleanup" function of the prot.file_converter  or ds.file_converter (FileConverter class) member variable.
    """
    prot.file_converter.cleanup()

if __name__ == "__main__":
    main()

