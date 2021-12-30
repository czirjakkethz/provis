Setting up Provis
=========================

Once you have every library installed, binary downloaded and of course provis installed (either pip or with the setup.py in the root directory) you will need to create a directory structure where provis feels at home. When this is set up you are ready to use provis and run it from anywhere.

Directory structure
--------------------
You will need a directory with a specific folder structure.
This way your directories will remain uncluttered and it makes it easy to create a pipline around provis.

The directories are needed as some parts of file loading is hard coded in provis and if the directory structure is not present errors will occur.

Easy option
^^^^^^^^^^^^
The easiest is if you clone provis from `github <https://github.com/czirjakkethz/provis>`_ and simply use this github directory as your base directory.

More versitile option
^^^^^^^^^^^^^^^^^^^^^^^
You will need a data/ and (potentially) a binaries/ directory as well.

 -- If you set the environment variables for the binaries (provis will then use these to find the binaries) then the binaries/ directory is not needed, but otherwise it is and the binaries from the Requirements section will all have to be copied into there. --

The data/ directory needs the following subdirectiories:
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
::

    data
    ├── data
    ├── pdb
    ├── img        
    └── tmp

The pdb directory is the location to store the pdb files to convert, as if a pdb file is stored here then it is enough to pass the pdb id (filename without extension) to provis and you do not have to pass a full path to the pdb file. 

The img directory stores all the screenshots of the outputted plots.

The tmp directory stores all temporary files created by provis, such as the .face and .vert files of MSMS or the .mol2 files needed for the bonds.


The binaries/ directory needs the following executables and subdirectiories:
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
::

    data
    ├── apbs
    ├── msms       
    └── pdb2pqr

Where apbs and msms are executables/binaries and the pdb2pqr is the directory downloaded from the official website, containing the pdb2pqr binary.
