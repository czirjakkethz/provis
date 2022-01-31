# provis - A python based PROtein VISualization toolkit

Welcome to the repository of Provis - a python based protein visualization library.

[Documentation](https://pro-vis.readthedocs.io/en/latest/), [installation instructions](https://pro-vis.readthedocs.io/en/latest/req.html).

Installation
===============

The easiest way to install is with ``` pip install provis ``` or by cloning this github repo and running: "python3 setup.py develop" from the base directory.

The library has a few dependencies that have to be installed seperately.

**FOLLOW THE INSTRUCTIONS IN "Requirements for Provis" AND "Setting up Provis"** Otherwise **provis** will not work.

Requirements for Provis
=========================

Provis is based on the idea of not to reinvent the wheel, so it requires quite a few third party packages and binaries.

Binaries
----------------------

### OpenBabel ###

Easiest to install on Linux is by calling:

``` sudo apt install openbabel ```

Alterantively: http://openbabel.org/wiki/Main_Page


#### PDB2PQR and APBS ####

These binaries are **required**, as they are used to extract surface inforamtion.

Download from: https://www.poissonboltzmann.org/


#### MSMS ####

 MSMS is optional. It is used to compute the surface, but a native method for the surface computation also exists in provis.

 Download MSMS form:
 https://ccsb.scripps.edu/mgltools/downloads/

 This tutorial might help:
 http://biskit.pasteur.fr/install/applications/deprecated/msms

Pip
----------------------

 All of the following can be dowloaded using pip and should be automatically installed when installing pip with provis:

* BioPython
* Trimesh
* PyVista
* Biopandas
* Torch
* Pyvtk
* Open3d
* rTree
* Panel
* Imageio-ffmpeg
   
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

