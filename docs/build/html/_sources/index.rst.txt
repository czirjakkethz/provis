.. provis documentation master file, created by
   sphinx-quickstart on Thu Nov 25 15:42:07 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to provis's Wiki!
==================================

Provis is a protein visualization library based on python. As the name suggests it is used to visualize proteins from their .pdb data formats (from the 3D coordinates). One can visualize the protein as a stick and point model or as a surface, and many different options can be specified to show the desired properties. 

The library and this documentation were created by Kristof Czirjak as his Bachelor's Thesis at ETH Zurich.

You can find installation instructions under
:ref:`setup`, a tutorial on 
:ref:`tutorial`, as well as the documentation of the `source code <https://github.com/czirjakkethz/provis>`_.

What **provis** is capable of:

.. raw:: html

	<video width="320" height="240" controls>
	  <source src="images/traj_hydrophob_animation.mp4" type="video/mp4">
	Your browser does not support the video tag.
	</video>

The following video shows the hydrophobic trajectory of a given protein.

|

.. toctree::
   :maxdepth: 2
   :caption: Getting started:
   
   req
   setup
   

.. toctree::
   :maxdepth: 3
   :caption: Documentation:

   provis

.. toctree::
   :maxdepth: 3
   :caption: Use and examples:
   
   information
   tutorial
   example
   dynamic_plotting
   
.. toctree::
   :maxdepth: 2
   :caption: Author(s)

   authors

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
