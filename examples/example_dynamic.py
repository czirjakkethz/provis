"""
This file was created by Kristof Czirjak as part of his Bachelor's Thesis - provis
"""
from provis.src.plotting.dynamic_structure import DynamicStructure
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
    Initializing the DynamicStructure class will prepare everything for plotting. It creates all the necessairy classes in the background and you are already good to go!
    """
    
    ds = DynamicStructure(name, base_path=base_path, density=density, plot_solvent=plot_solvent, msms=msms, notebook=notebook)
    """
    Third:
    Plot!
    
    Use the member functions of the DynamicStructure to plot the dynamic trajectory of the protein.
    """
    ds.plot_surface()#feature="hydrophob")

    r = Residue(29)
    r.add_residue(50)
    r. add_residue(1, 1)
    r.remove_residue(1, 1)

    ds.plot_backbone()
    ds.plot_atoms(res=r)
    ds.plot_bonds()
    ds.plot_vw()
    ds.plot_stick_point()
    ds.plot_residues()
    ds.plot(atoms=1, box=1, bonds=1, vw=0, residues=0, res=r, bb=0)


    """
    And finally clean up everything with the "cleanup" function of the prot.file_converter (FileConverter class) member variable.
    """
    # ds.file_converter.cleanup()

if __name__ == "__main__":
    main()

