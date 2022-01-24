Additional information about provis
=====================================

This section will try to explain the general code architecture and design decisions of the **provis** library.

Classes
++++++++++

Protein
^^^^^^^^^

The **Protein** class encapsulates the whole **provis** library and is the easiest way to plot your (static) proteins. You can initialize it with a simple pdb filename in a single line of code and plot in the next.

FileConverter
^^^^^^^^^^^^^^

The **FileConverter** class is responsible for creating the required files from the *.pdb* file. This class is hidden in other classes and only called when a file is needed/missing.

IMPORTANT: the **FileConverter** class will not overwrite any existing files. This decision was made to not create the exact same file again and again. But this also means that if you have temporary file (used by the **FileConverter** class) that is corrupted or has bad information it will be used by **provis**. If you have any such problems you can always just delete all temporary files and then they will simply be newly created.

DataHandler
^^^^^^^^^^^^

The **DataHandler** class processes the structural information of the protein and creates the meshes for plotting. The meshes created by this class are passed to the **Structure** class for plotting.

SurfaceHandler
^^^^^^^^^^^^^^^^

The **SurfaceHandler** class processes the structural information of the protein and creates the meshes for plotting. The meshes created by this class are passed to the **Surface** class for plotting.

Structure
^^^^^^^^^^

The Structure class handles every non-surface related plotting.

The atoms can be plotted as a point cloud, bonds or the backbone of the molecule can be visualized and so can the Van-der-Waals radius of the atoms, as well as the residues of the protein.

Surface
^^^^^^^^

The Surface class handles surface related plotting.

Apart from the standard surface plot, interesting surface features - hydrophobicity, shape index, charges - can also be visualized.

Residue
^^^^^^^^

The **Residue** class is a very basic class. With this class you can specify a given residue on a given chain that you want to mark in the plot. It is passed to the plotting classes and a red box is drawn around the specified residue.

DynamicStructure
^^^^^^^^^^^^^^^^^

See:
:ref:`dynamic_plotting`.
