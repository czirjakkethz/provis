.. _tutorial:

******************
How to use Provis
******************

This section will explain how to use provis and how to get your desired plots.

Loading a pdb
###############

Provis uses the information in *.pdb* files to plot your desired protein. The most straightforward way is to pass the full path of the file to provis, as the *.pdb* file can be saved anywhere if you pass the full path. For example save the path to the **name** variable:

 .. code-block:: python

   name = "/home/username/provis/data/pdb/2fd7.pdb"

As explained in 
:ref:`setup`, provis requires a specific directory structure. If you have your *.pdb* file stored in the **/data/pdb** directory you do not have to specify a full path to the pdb file, but simply the name of the file:


 .. code-block:: python

   name = "2fd7.pdb" # "2fd7" also works

Setting Protein class input variables
#######################################

.. code-block:: python

	prot = Protein(name, base_path=None, density=density, plot_solvent=solvent, msms=False, notebook=False)

Setting the input variables of the Protein class is the next step to take. You have to pass the **name** variable from before.

The path to the special direcotry explained in 
:ref:`setup` has to also be provided. This is passed in the **base_path** variable and should point to the root directory of the **/data** and **/binaries** directories.

MSMS
++++++

You will have to also specify if you want to plot the **msms** binary version of the surface or the simpler native mesh. If the **msms** option is chosen you can also specify the **density** of the triangulation (to be passed to the **msms** binary).

Solvent atoms can also be plotted by setting the **plot_solvent** variable to *True*.

And finally if you are working in a Jupyter Notebook like environment then set the **notebook** varaible to *True*.

How the Protein class works
#############################

The Protein class encompasses and combines all other classes and files of provis.

It first figures out the location of the *.pdb* file using the **NameChecker** class. Then it instantiates a **FileConverter** class so the necessairy temporary files can be converted from our *.pdb* file.

Then a **DataHandler** class is created. This class calculates all the necessairy non-surface related meshes, such as the Spheres corresponding to each atom or the backbone of the molecule. Next, a **SurfaceHandler** class handles all the surface related computation.

Finally, a **Structure** and a **Surface** class is created. These classes are instantiated with the above mentioned classes as input variables to ensure that no data duplication occurs and that the two will handle the same protein. 

Plotting
############

Plotting can be achieved by calling the member functions of the **Structure** and the **Surface** classes. For example for the **prot** class instance defined above the bonds of the molecule can be plotted as follows:

.. code-block:: python

	prot.structure.plot_bonds()
	
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

Setting a good camera position is very important. By default the camera portion is set to **[0, max * 4, 0]**, where **max** is the largest deviation of the coordinates from the center of the molecule. This ensures that the whole molecule is visible in the plot window and that the camera will always face the same direction when plotting dynamically.

To set the camera position manually the **DataHandler** class' instance variables named **DataHandler._cam_pos** (the default camera position) and **DataHandler._max_coords** (the maximum deviation, as explained above) might be helpful.
        

