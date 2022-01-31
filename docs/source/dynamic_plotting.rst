
.. _dynamic_plotting:

Dynamic Structures
====================

Dynamic structures showcase the transformation of a molecule in time. This trajectory of change can also be plotted with provis.

The **DynamicPlotter** class achieves exactly that. Unlike the **Protein** class, the **DynamicPlotter** class has its own member functions used for plotting. This is due to the fact that each model (molecule at a given time) has to be plotted seperately and this is achieved by a loop within the member function.

For large molecules and a long trajectory the surface meshes can take a while to compute. It is advised to first precompute the meshes (they will be stored in *.obj* files in the **"root directory"/data/meshes** directory). Once the *.obj* files exist plotting will be seamless.

The result of the plot will be saved to the **"root directory"/data/img** folder as an *.mp4* file.

How to use the DynamicPlotter class?
++++++++++++++++++++++++++++++++++++++++++

The **DynamicPlotter** class is very similar to the **Protein** class.

Import the necessairy files.

First:
Define variables needed later:

.. code-block:: python

    name = "2fd7"
    density = 3.0


Second:
Create a **Protein** class instance. Initialize it with your *.pdb* file name and other parameters.

.. code-block:: python

    prot = Protein(name, base_path=None, density=density)


Initialize the **DynamicPlotter** class. This creates all the necessairy classes in the background and you are already good to go!

.. code-block:: python

    dp = DynamicPlotter(prot, msms=msms, notebook=notebook, plot_solvent=plot_solvent)

Third:
Plot!

This class has its own plotting methods. Here is a complete list:

.. code-block:: python

    dp.plot_backbone()
    dp.plot_atoms()
    dp.plot_bonds()
    dp.plot_vw()
    dp.plot_stick_point()
    dp.plot_residues()
    dp.plot_structure(atoms=1, box=1, bonds=1, vw=0, residues=0, res=r, bb=0)

    dp.plot_surface()
    dp.plot_hydrophob()
    dp.plot_shape()
    dp.plot_charge()


And finally, cleaning up is also possible with the "cleanup" function of the Protein.file_converter (**FileConverter** class) member variable.

.. code-block:: python

    prot.file_converter.cleanup()

