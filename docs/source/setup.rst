Requirements for Provis
=========================

Pip
----------------------

Provis is based on the idea of not to reinvent the wheel, so it requires quite a few third party packages. All of the following should be downloadable with pip

* BioPython
* PyVista
* Pandas
* Numpy (included in Pandas)
* Biopandas
* Trimesh

Binaries
----------------------

DSSP
^^^^^^^^^^^^^^^^^^^^^^

 DSSP is **required**. It is used to extract surface information. 

 One can easily downloaded using the `conda
 <https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html>`_ installer.

 .. code-block:: bash

   conda install -c salilab dssp

MSMS
^^^^^^^^^^^^^^^^^^^^^^
 MSMS is optional. It is used to compute the surface, but a native method for the surface computation also exists in provis.

 Download MSMS form:
 https://ccsb.scripps.edu/mgltools/downloads/

 This tutorial might help:
 http://biskit.pasteur.fr/install/applications/deprecated/msms

Jupyter Notebook
---------------------

If you want to run a fully interactive provis in the Jupyter Notebook environment you will have to install the following additional packages. Failing to do so will simply produce snapshot images instead of 3D plots.

.. code-block:: bash

    pip install ipyvtklink