.. _tutorial:

******************
How to use Provis
******************

This section will explain how to use provis and how to get your desired plots.

In short you have to initialize a **Protein** class instance with the desired *.pdb* file. Then you have to create a **Plotter** or a **DynamicPlotter** class instance and pass the previously created **Protein** to it. Once the **(Dynamic)Plotter** is initialized you can plot using its plotting methods.

Loading a pdb
###############

Provis uses the information in *.pdb* files to plot your desired protein. The most straightforward way is to pass the full path of the file to provis, as the *.pdb* file can be saved anywhere if you pass the full path. For example save the path to the **name** variable:

 .. code-block:: python

   name = "/home/username/provis/data/pdb/2fd7.pdb"

As explained in 
:ref:`setup`, provis requires a specific directory structure. If you have your *.pdb* file stored in the **/data/pdb** directory you do not have to specify a full path to the pdb file, but simply the name of the file:


 .. code-block:: python

   name = "2fd7.pdb" # "2fd7" also works

Initializing the Protein class
#######################################

The Protein class encompasses and combines all classes of the **provis.src.processing** package.

It first figures out the location of the *.pdb* file using the **NameChecker** class. Then it instantiates a **FileConverter** class so the necessairy temporary files can be converted from our *.pdb* file.

Then a **DataHandler** class is created. This class calculates all the necessairy non-surface related meshes, such as the Spheres corresponding to each atom or the backbone of the molecule. Next, a **SurfaceHandler** class handles all the surface related computation.

.. code-block:: python

    prot = Protein(name, base_path=None, density=3.0)

Specify the **name** of the *.pdb* file.

The path to the special direcotry explained in 
:ref:`setup` has to also be provided. This is passed in the **base_path** variable and should point to the root directory of the **/data** and **/binaries** directories.


Plotter class
###############

Use the Plotter class to plot. At least one Protein has to be passed. If two proteins are passed then they will be plotted side-by-side.
It is possible to add more proteins later using the Plotter.add_protein(Protein) method.

.. code-block:: python

    prot2 = Protein(name, model_id=30)
    plot = Plotter(prot, prot2, msms=msms, notebook=notebook)

MSMS
++++++

You will have to also specify if you want to plot the **msms** binary version of the surface or the simpler native mesh. If the **msms** option is chosen you can also specify the **density** of the triangulation (to be passed to the **msms** binary).

Solvent atoms can also be plotted by setting the **plot_solvent** variable to *True*.

And finally if you are working in a Jupyter Notebook like environment then set the **notebook** varaible to *True*.

Plotting
+++++++++

Plotting can be achieved by calling the member functions of the **Structure** and the **Surface** classes. For example for the **prot** class instance defined above the bonds of the molecule can be plotted as follows:

.. code-block:: python

	plot.plot_bonds()
	
All plotting functions have the following input variables:
 - box (bool): If True bounding box will be plotted around molecule.
 - res (Residue): Specified residue will be marked with a red box.
 - outname (str): Save image of the plot to the file passed in this variable (otherwise saved in data/img).
 - camera (pyvista.camera): A pyvista camera object to be make it easier to set a fixed camera position to compare two molecules.
 
 Some of the plotting functions have additional input variables. One example; **plot_bonds()**:
 
 - colorful (bool): If True different bond types will be plotted in different colors.
 
  	Single bonds: white
        Double bonds: blue
        Triple bonds: green
        Amide bonds: red
        Aromatic bonds: purple
        Undefined/Anything else: black


Camera
++++++++

Setting a good camera position is very important. By default the camera portion is set to **[0, max * 3, 0]**, where **max** is the largest deviation of the coordinates from the center of the molecule. This ensures that the whole molecule is visible in the plot window and that the camera will always face the same direction when plotting dynamically.

To set the camera position manually the **DataHandler** class' instance variables named **DataHandler._cam_pos** (the default camera position) and **DataHandler._max_coords** (the maximum deviation, as explained above) might be helpful.
        
The output
+++++++++++

The output will be an interactive **vtk.Window** window. The output will also be saved as an image to *"root directory"/data/img*.

The following image are the bonds of the 1st and 31st model of a given dynamic trajectory.

.. image:: images/traj_30_bonds_msms.png
  :width: 600
  :alt: The following images are the bonds of the 1st and 31st model of a given dynamic trajectory.

The following image are the atoms of the 1st and 31st model of a given dynamic trajectory.

.. image:: images/traj_30_atoms_msms.png
  :width: 600
  :alt: The following images are the atoms of the 1st and 31st model of a given dynamic trajectory.

