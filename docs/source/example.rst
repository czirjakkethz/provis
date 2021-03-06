Example
===========

This is an example file to showcase the easiest way to run provis in particular how to plot a single protein. For this you should have this file in the root directory of the special directroy structure specified in the setup section of the documentation. Otherwise set the **base_path** variable of the **NameChecker** object.

If this is fullfilled path to the *"root directory"/data/tmp* will automatically be found.
This way you can have your pdb files nicely organized in the data/pdb directory (or simply have them in the root directory). 
Your temporary files will be in the *data/tmp* directory and the screenshots of the plots in the *data/img* directory.
    

Import the necessairy files.

.. code-block:: python

	from provis.src.processing.protein import Protein
	from provis.src.processing.residue import Residue

First:
Define variables needed later:

.. code-block:: python

    name = "2fd7"
    density = 3.0
    plot_solvent = False
    msms = True
    notebook = False


Second:
Create a **Protein** class instance. Initialize it with your *.pdb* file name and other parameters.

If you want to plot multiple proteins (or different models of the same trajectory) this is also possible. Simply create a second **Protein** class instance and pass it to the **Plotter**

.. code-block:: python

    prot = Protein(name, base_path=None, density=density)
    prot2 = Protein(name, base_path=None, density=density, model_id=30)

Initialize the **Plotter** class. This creates all the necessairy classes in the background and you are already good to go!

.. code-block:: python

    plot = Plotter(prot, prot2, msms=msms, notebook=notebook, plot_solvent=plot_solvent)

Third:
Plot!

Use the **Plotter** to plot.

.. code-block:: python

    plotter.plot_backbone()
    plotter.plot_atoms()
    plotter.plot_bonds()
    plotter.plot_vw()
    plotter.plot_stick_point()
    plotter.plot_residues()
    r = Residue(1)
    r.add_residue(3)
    r. add_residue(1, 1)
    r.remove_residue(1, 1)
    plotter.plot_structure(atoms=1, box=1, bonds=1, vw=0, residues=0, res=r, bb=0)
    
    plotter.plot_surface()    
    plotter.plot_hydrophob()
    plotter.plot_shape()
    plotter.plot_charge()


And finally clean up everything with the "cleanup" function of the **Protein.file_converter** (**FileConverter** class) member variable.

.. code-block:: python

    prot.file_converter.cleanup()


The following image shows the hydrophobicity of the 2fd7 protein.

.. image:: images/2fd7_hydrophob.png
  :width: 600
  :alt: The following image shows the hydrophobicity of the 2fd7 protein.

