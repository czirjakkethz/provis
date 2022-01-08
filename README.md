# provis - A python based PROtein VISualization toolkit

Welcome to the repository of Provis - a python based protein visualization library.

The easiest way to install is with ``` pip install provis ```, but the library has a few dependencies that have to be installed seperately.
Please visit the  [documentation](https://pro-vis.readthedocs.io/en/latest/) for [installation instructions](https://pro-vis.readthedocs.io/en/latest/req.html) or continue reading below.

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
   
Jupyter Notebook
---------------------

If you want to run a fully interactive provis in the Jupyter Notebook environment you will have to install the following additional packages. Failing to do so will simply produce snapshot images instead of 3D plots.

``` pip install ipyvtklink ```

Setting up Provis
=========================

Once you have every library installed, binary downloaded and of course provis installed (either pip or with the setup.py in the root directory) you will need to create a directory structure where provis feels at home. When this is set up you are ready to use provis and run it from anywhere.

Directory structure
--------------------
You will need a directory with a specific folder structure.
This way your directories will remain uncluttered and it makes it easy to create a pipline around provis.

The directories are needed as some parts of file loading is hard coded in provis and if the directory structure is not present errors will occur.

### Easy option ###

The easiest is if you clone provis from [github](https://github.com/czirjakkethz/provis) and simply use this github directory as your base directory.

### More versitile option ###

You will need a data/ and (potentially) a binaries/ directory as well.

 -- If you set the environment variables for the binaries (provis will then use these to find the binaries) then the binaries/ directory is not needed, but otherwise it is and the binaries from the Requirements section will all have to be copied into there. --

#### The data/ directory needs the following subdirectiories: ####


    data
    ├── data
    ├── pdb
    ├── img        
    └── tmp

The pdb directory is the location to store the pdb files to convert, as if a pdb file is stored here then it is enough to pass the pdb id (filename without extension) to provis and you do not have to pass a full path to the pdb file. 

The img directory stores all the screenshots of the outputted plots.

The tmp directory stores all temporary files created by provis, such as the .face and .vert files of MSMS or the .mol2 files needed for the bonds.


#### The binaries/ directory needs the following executables and subdirectiories: ####


    data
    ├── apbs
    ├── msms       
    └── pdb2pqr

Where apbs and msms are executables/binaries and the pdb2pqr is the directory downloaded from the official website, containing the pdb2pqr binary.
