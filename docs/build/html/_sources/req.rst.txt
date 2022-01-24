.. _req:


Requirements for Provis
=========================

Provis is based on the idea of not to reinvent the wheel, so it requires quite a few third party packages and binaries.

Binaries
----------------------

OpenBabel
^^^^^^^^^^^^^^^

Easiest to install on Linux is by calling:

 .. code-block:: bash

   sudo apt install openbabel

Alterantively: http://openbabel.org/wiki/Main_Page


PDB2PQR and APBS
^^^^^^^^^^^^^^^^^^^^^

These binaries are **required**, as they are used to extract surface inforamtion.

Download from: https://www.poissonboltzmann.org/


MSMS
^^^^^^^^^^^^^^^^^^^^^^
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
   
Jupyter Notebook
---------------------

If you want to run a fully interactive provis in the Jupyter Notebook environment you will have to install the following additional packages. Failing to do so will simply produce snapshot images instead of 3D plots.

.. code-block:: bash

    pip install ipyvtklink
    

