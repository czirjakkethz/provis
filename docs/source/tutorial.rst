******************
How to use Provis
******************

This section will explain how to use provis and how to get your desired plots.

Loading a pdb
###############

Provis uses the information in .pdb files to plot your desired protein. The most straightforward way is to pass the full path of the file to provis. For example save the path to the 'name' variable:

 .. code-block:: python

   name = "/home/username/provis/data/pdb/2fd7.pdb"

As explained in 
:ref:`setup`, provis requires a specific directory structure. If you have your .pdb file stored in the 'data/pdb/' directory you do not have to specify a full path to the pdb file, but simply the name of the file:


 .. code-block:: python

   name = "2fd7.pdb" # "2fd7" also works



You will need to provide a .pdb file describing the protein you want to work with. To load this file, you will have to tell provis where the file is located.


