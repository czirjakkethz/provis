"""
This file was created by Kristof Czirjak as part of his Bachelor's Thesis - provis
"""
from provis.src.plotting.dynamic_structure import DynamicStructure
from provis.src.processing.residue import Residue
from provis.src.processing.protein import Protein

def main():
    """
    This is an example to showcase briefly showcase the functionalities of provis through a few examples.
    
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
    Define a Protein class. This class is used for static plotting. Initialize it with the pdb name and then plot using its internal members.
    
    Use the methods of prot.structure (Structure class) and prot.surface (Surface class) to plot.
    """
    prot = Protein(name, base_path=None, density=density, plot_solvent=False, msms=True, notebook=False)
    prot.surface.plot_shape()

    """
    Initializing the DynamicStructure class is used for dynamic plotting.
   
    Use the methods of the DynamicStructure class to plot the dynamic trajectory of the protein.   
    """
    
    ds = DynamicStructure(name, base_path=base_path, density=density, plot_solvent=plot_solvent, msms=msms, notebook=notebook)
    ds.plot_atoms()


    """
    And finally clean up everything with the "cleanup" function of the prot.file_converter  or ds.file_converter (FileConverter class) member variable.
    """
    # ds.file_converter.cleanup()

if __name__ == "__main__":
    main()

