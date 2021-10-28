
from Bio.PDB import *
from Bio import SeqIO
import numpy as np
import pyvista

# load atoms into a python list #TODO: probably faster in numpy array
parser = PDBParser() 
structure = parser.get_structure("2fd7", "2fd7.pdb") 
atoms = structure.get_atoms()
atoms_list = list(atoms)


# Create a numpy array with the coordinates of all atoms
num_atoms = len(atoms_list)
#print("Number of atoms in protein is:", num_atoms)
coord_array = np.empty((num_atoms, 3))
radius_array = [None] * num_atoms
color_array = [None] * num_atoms


# Create different spheres for the different atoms
carbon = pyvista.Sphere(radius=6.7, phi_resolution=10, theta_resolution=10)
oxygen = pyvista.Sphere(radius=4.8, phi_resolution=10, theta_resolution=10)
nitrogen = pyvista.Sphere(radius=5.6, phi_resolution=10, theta_resolution=10)
sulfur = pyvista.Sphere(radius=8.8, phi_resolution=10, theta_resolution=10)
calcium = pyvista.Sphere(radius=19.4, phi_resolution=10, theta_resolution=10)

# Create list for geoms
geoms = [None] * num_atoms

i = 0
for atom in atoms_list:
    coord_array[i][0:3] = atom.get_coord() * 10
    element = atom.element
    radius = 0
    if element == "C":
        radius = 0.67
        geoms[i] = carbon
        color_array[i] = '#000000'
    elif element == "O":
        radius = 0.48
        geoms[i] = oxygen
        color_array[i] = '#FFFFFF'
    elif element == "N":
        radius = 0.56
        geoms[i] = nitrogen
        color_array[i] = '#007DEE'
    elif element == "S":
        radius = 0.88
        geoms[i] = sulfur
        color_array[i] = '#FFDC61'
    elif element == "CA":
        radius = 1.94
        geoms[i] = calcium
        color_array[i] = '#38761D'
    radius_array[i] = radius
    i = i+1


#print(radius_array)
"""
# ALTERNATE PLOTTING, spheres are points not full faced (if add_points NOT add_mesh)
# Can not figure out coloring, probably have to add with each, individually
# create pyvista data structure, mesh(?)
pdata = pyvista.PolyData(coord_array)


# create pyvista plotter
pl = pyvista.Plotter()

# create many spheres from the point cloud
#sphere = pyvista.Sphere(radius=0.5, phi_resolution=10, theta_resolution=10)
pc = pdata.glyph(scale=False, geom=geoms)
pl.add_mesh(pc, specular=1, specular_power=150,
                 smooth_shading=True, show_scalar_bar=False)

pl.show()
#pc.plot()#cmap='Reds')

"""
# test_coord = np.array([[10, 0, 0],[100, 0, 0]]) 
# test_geoms = [carbon, oxygen]
# test_color = [0, 1]
# ALTERNATE PLOTTING, shperes are solid, can't figure out color
# create pyvista data structure
pdata = pyvista.PolyData(coord_array)
pdata['point_color'] = color_array

# create many spheres from the point cloud
#sphere = pyvista.Sphere(radius=0.5, phi_resolution=10, theta_resolution=10)
pc = pdata.glyph(scale=False ,geom=geoms, factor=1)

pc.plot(scalars='point_color')#cmap='Reds')

#"""