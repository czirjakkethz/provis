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
    """
    Second:
    Initializing the Protein. This class instance needs to be passed to the Plotter.
    """
    
    prot = Protein(name, base_path=None, density=density)

    """
    Third:
    Plot!
    
    Use the Plotter to plot. Create the necessairy lists of meshes with the DataHandler object.
    """
    plot = Plotter(prot, msms=msms, notebook=False)
    
    atom_data = prot._data_handler.get_atoms(show_solvent=False, model_id=prot._model_id) 
    atoms, col_a, _ = prot._data_handler.get_atom_mesh(atom_data, vw=False) # second arg = True => show vw spheres instead of "normal" radius
   
        
    plot.manual_plot(atoms=atoms, col_a=col_a)    
  
    """
    And finally clean up everything with the "cleanup" function of the prot.file_converter (FileConverter class) member variable.
    """
    # prot.file_converter.cleanup()

if __name__ == "__main__":
    main()

