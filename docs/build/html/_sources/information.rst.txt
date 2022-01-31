General information about provis
=====================================

This section describes the general code architecture and design decisions of the **provis** library.

MSMS
+++++

As explained in the Getting started section, there are two ways to compute the surface. Using the **msms** binary or natively. 

Therefore, it is very important to specify which version you want to use. This can be done by setting the **msms** input variable in the constructor of the class that handles the surface plotting. These classes are: **Surface**, **Protein** and **DynamicPlotter** (all three classes ultimately rely on the **Surface** class in the background).

Classes
++++++++++

Protein
^^^^^^^^^

The **Protein** class encapsulates the whole processing package of the **provis** library. You can initialize it with a simple pdb filename and the class can then be used to retrieve all the necessairy information for plotting.


NameChecker
^^^^^^^^^^^^^^

The **NameChecker** class is responsible for storing the paths to all the files used by **provis**. All other classes require a **NameChecker** class instance to be passed. Passing this **NameChecker** instance ensures that all other classes will work with the same molecule.

FileConverter
^^^^^^^^^^^^^^

The **FileConverter** class is responsible for creating the required files from the *.pdb* file. This class is hidden in other classes and only called when a file is needed/missing.

IMPORTANT: the **FileConverter** class will not overwrite any existing files. This decision was made to not create the exact same file again and again. But this also means that if you have temporary file (used by the **FileConverter** class) that is corrupted or has bad information it will be used by **provis**. If you have any such problems you can always just delete all temporary files and then they will simply be newly created.

DataHandler
^^^^^^^^^^^^

The **DataHandler** class processes the structural information of the protein and creates the meshes for plotting.

SurfaceHandler
^^^^^^^^^^^^^^^^

The **SurfaceHandler** class processes the structural information of the protein and creates the meshes for plotting.


Plotter
^^^^^^^^^^

The **Plotter** class is responsible for static plotting. It plots the internal data of the passed **Protein** class(es). If multiple **Protein** classes are passed, then they will be plotted in a side-by-side manner. 

DynamicPlotter
^^^^^^^^^^^^^^^^^

Similar to the **Plotter** class, the **DynamicPlotter** class plots instances of the **Protein** class. This class however, focuses on dynamic plotting. It plots the trajectory of a molecule.

See:
:ref:`dynamic_plotting`.


Residue
^^^^^^^^

The **Residue** class is a very basic class. With this class you can specify a given residue on a given chain that you want to mark in the plot. It is passed to the plotting classes and a red box is drawn around the specified residue.



Design decisions
+++++++++++++++++

Since creating meshes requires a lot of computation it was decided to store the results of these computations to files. This way if a protein is analysed with **provis** for the second time, then results of the expensive computation will already be known. This way one can simply load the information from the files and plot right away.

Some classes of **provis** need the same information and need access to the same methods. For this reason, most classes of **provis** can be initialised with already existing instances of the internal variables. This way the exact same computations (for example: initialising the passed class) will be mitigated. Passing initialised classes also insures that the same molecule is being considered throughout the whole program.

To ensure that only necessary computations are run the existence of the object to be computed is always checked. Both files and class instance variables are checked. If the given object already exists (and possible other checks are passed) then this object will not be recomputed.
**WARNING: this decision can result in an error if a new file with a previously existing name is analysed. **Provis** will not know that the molecule is different and errors might occur. To fix this simply delete (or move) all temporary files related to this molecule.** This is a very rare problem, if you use **provis** normally and do not modify existing files then you will not run into this problem.

Code architecture
++++++++++++++++++

In the following will show how the above mentioned classes depend on one another.

One can easily observe that all classes depend on the NameChecker class. And since some classes have other dependencies as well, duplication of the exact same class instance would seem quite likely. This is the reason (as described above) why all classes can be initialised with pre-existing instances of the necessary classes, so instead of duplication they will be shared. (Disclaimer: the classes can also be initialised empty, but then the above mentioned duplication occurs.)

::

	FileConverter
	   └── NameChecker          


::

	DataHandler
	   ├── NameChecker
	   └── FileConverter          

::

	SurfaceHandler
	   ├── NameChecker
	   ├── FileConverter
	   └── DataHandler      


::

	Protein
	   ├── NameChecker
	   ├── FileConverter
	   ├── DataHandler
	   └── SurfaceHandler      


::

	Plotter
	   └── Protein
	         
::

	DynamicPlotter
	   └── Protein
	         

