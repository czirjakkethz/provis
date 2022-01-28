
Download Provis
================

Provis can be downloaded by pip or from github.

**Pip:**

.. code-block::

	pip install provis
	
**Github:** https://github.com/czirjakkethz/provis

If provis was downoaded from the git repository then the following command has to be run: (run it from the root directory)

.. code-block::
	
	python3 setup.py develop
	
This will set up provis as a python library on your machine, but also download all the python (pip) dependencies.

.. _setup:

Setting up Provis
=========================

Some file paths are hard coded in provis in order to maintain an uncluttered and organized directory structure.

However, this also means, that for provis to work this specific directory structure has to exist.

Directory structure
--------------------

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



Easy option
^^^^^^^^^^^^
Clone provis from `github <https://github.com/czirjakkethz/provis>`_ and simply use this git directory (provis) as the base directory.

More versitile option
^^^^^^^^^^^^^^^^^^^^^^^

A data/ and (potentially) a binaries/ directory within the root directory will have to be created.

 -- If the environment variables for the binaries are set then the binaries/ directory is not needed (as provis will then use these variables to find the binaries). Otherwise the binaries from the Requirements 
 :ref:`req` section will all have to be copied into the binaries/ directory. --

Subdirectiories of the data/ directory:
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
::

    data
    ├── data
    ├── pdb
    ├── img   
    ├── meshes
    └── tmp

The pdb directory is the location to store the pdb files. If a *.pdb* file is stored here then it is enough to pass the pdb id (filename without extension) to provis. Otherwise the full path to the *.pdb* file needs to be passed. 

The img directory stores all the screenshots of the outputted plots.

The tmp directory stores all temporary files created by provis, such as the *.face* and *.vert* files of MSMS or the *.mol2* files needed for the bonds.

Binaries
---------

After installing the binaries either set the environment variables to specify the path of their location or manually move the binaries to where provis can find them.


Binaries directory
^^^^^^^^^^^^^^^^^^^^

If you are running Ubuntu (20.04.3 LTS) and installed provis by cloning the github repository then you are all set.

Otherwise move the binary files to the 'binaries' directory as explained below.


Subdirectiories and executables of the binaries/ directory:
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
::

    binaries
    ├── apbs
    ├── msms       
    └── pdb2pqr     
         └── pdb2pqr

**apbs** and **msms** are executables and *pdb2pqr* is the directory downloaded from the official website, containing the **pdb2pqr** binary.


Setting environment variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Alternatively the environment variables can be set to point to the binary files.

Set the MSMS_BIN, PDB2PQR_BIN and APBS_BIN variables to the full path to their appropriate binary files.

Example:

.. code-block::

	export MSMS_BIN='/home/username/Downloads/msms'




