.. _req:


Requirements for Provis
=========================

Provis is based on the idea of not to reinvent the wheel, so quite a few third party packages and binaries are required to run it.

Binaries
----------------------

Binaries are 3rd party, ready to use programs, that are provided to the user as is. Here is a list of the required external programs and installation instructions.

OpenBabel
^^^^^^^^^^^^^^^

Easiest to install on Linux is by calling:

 .. code-block:: bash

   sudo apt install openbabel

Alterantively: http://openbabel.org/wiki/Main_Page

OpenBabel is needed to create the *mol.2* files. These files store the bond information.

PDB2PQR and APBS
^^^^^^^^^^^^^^^^^^^^^

Download from: https://www.poissonboltzmann.org/

These binaries are required for the surface feature plotting. The binaries create the 

MSMS
^^^^^^^^^^^^^^^^^^^^^^

 MSMS is optional. It is used to compute the surface, but a native method for the surface computation also exists in provis (although it is slower and chemically less accurate).

 Download MSMS form:
 https://ccsb.scripps.edu/mgltools/downloads/

 This tutorial might help:
 http://biskit.pasteur.fr/install/applications/deprecated/msms

Pip
----------------------

If provis was downloaded via pip (and not from the `github <https://github.com/czirjakkethz/provis>`_) then all of the following packages should be installed. 

If not, then run the following command in the root directory of provis:

.. code-block:: bash
	
	python3 setup.py develop
    
This should pip install everything, including the provis package.


Here is a list of provis' pip dependencies:

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
   
    

