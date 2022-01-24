
.. _dynamic_plotting:

Dynamic Structures
====================

Dynamic structures showcase the transformation of a molecule in time. This trajectory of change can also be plotted with provis.

The **DynamicStructure** class is used for this purpose. Unlike the **Protein** class, the **DynamicStructure** class has its own member functions used for plotting. This is due to the fact that each model (molecule at a given time) has to be plotted seperately and this is achieved by a loop within the member function.

For large molecules and a long trajectory the surface meshes can take a while to compute. It is advised to first precompute the meshes (they will be stored in *.obj* files in the **"root directory"/data/meshes** directory). Once the *.obj* files exist plotting will be seamless.



