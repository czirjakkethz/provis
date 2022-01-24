.. _setup:

Setting up Provis
=========================

Once you have every library installed, binary downloaded - and of course provis installed (either pip or with the setup.py in the root directory) - you will need to create an environment where provis feels at home. When this is set up you are ready to use provis and run it from anywhere.

Binaries
---------

After installing the binaries you can either set environment variables to specify the path or you can manually move the binaries to where provis can find them.

Binaries directory
^^^^^^^^^^^^^^^^^^^^

If you are running Ubuntu (20.04.3 LTS) and installed provis by cloning the github repository then you are all set.

Otherwise you should move your binary files to the 'binaries' directory as explained below ::refMore versitile option.

Setting environment variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Alternatively you can simply set the environment variables to point to the binary files.

You should set MSMS_BIN, PDB2PQR_BIN and APBS_BIN to the location where the appropriate binary files are located.

Example:

.. code-block::

	export MSMS_BIN='/home/username/Downloads/msms'

Directory structure
--------------------
In order to keep your directories uncluttered and make it easy to create a pipline, provis requires a specific directory structure.See below:

::

	provis
	   ├── data
	   │   ├── data/
	   │   ├── pdb/
	   │   ├── img/
	   │   ├── meshes/
	   │   └── tmp/
	   └── binaries          
	       ├── apbs
	       ├── msms
	       └── pdb2pqr
		    └── pdb2pqr

The directories are needed as some parts of file loading is hard coded so if the directory structure is not present errors will occur.

Easy option
^^^^^^^^^^^^
If you clone provis from `github <https://github.com/czirjakkethz/provis>`_ and simply use this (provis) git directory as your base directory.

More versitile option
^^^^^^^^^^^^^^^^^^^^^^^

You will need to create a data/ and (potentially) a binaries/ directory within your root directory.

 -- If you set the environment variables for the binaries (provis will then use these to find the binaries) then the binaries/ directory is not needed. Otherwise the binaries from the Requirements 
 :ref:`req` section will all have to be copied into the binaries/ directory. --

The data/ directory needs the following subdirectiories:
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
::

    data
    ├── data
    ├── pdb
    ├── img   
    ├── meshes
    └── tmp

The pdb directory is the location to store the pdb files to convert, as if a *.pdb* file is stored here then it is enough to pass the pdb id (filename without extension) to provis and you do not have to pass a full path to the *.pdb* file. 

The img directory stores all the screenshots of the outputted plots.

The tmp directory stores all temporary files created by provis, such as the *.face* and *.vert* files of MSMS or the *.mol2* files needed for the bonds.


The binaries/ directory needs the following executables and subdirectiories:
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
::

    binaries
    ├── apbs
    ├── msms       
    └── pdb2pqr     
         └── pdb2pqr

Where **apbs** and **msms** are executables/binaries and the *pdb2pqr* is the directory downloaded from the official website, containing the **pdb2pqr** binary.
