
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
coord_array = np.empty((num_atoms, 3))
# Create list for colors
color_array = [None] * num_atoms
# Create list for geoms
geoms = [None] * num_atoms

# Create different spheres for the different atoms
# Set resolution variables so only need to change once
phi_res = 10
theta_res = 10
carbon = pyvista.Sphere(radius=0.67, phi_resolution=phi_res, theta_resolution=theta_res)
oxygen = pyvista.Sphere(radius=0.48, phi_resolution=phi_res, theta_resolution=theta_res)
nitrogen = pyvista.Sphere(radius=0.56, phi_resolution=phi_res, theta_resolution=theta_res)
sulfur = pyvista.Sphere(radius=0.88, phi_resolution=phi_res, theta_resolution=theta_res)
calcium = pyvista.Sphere(radius=1.94, phi_resolution=phi_res, theta_resolution=theta_res)

i = 0
for atom in atoms_list:
    coord_array[i][0:3] = atom.get_coord() 
    element = atom.element
    if element == "C":
        geoms[i] = carbon
        color_array[i] = '#333333'
    elif element == "O":
        geoms[i] = oxygen
        color_array[i] = '#EA2128'
    elif element == "N":
        geoms[i] = nitrogen
        color_array[i] = '#007DEE'
    elif element == "S":
        geoms[i] = sulfur
        color_array[i] = '#FFDC61'
    elif element == "CA":
        geoms[i] = calcium
        color_array[i] = '#38761D'
    i = i+1

# create pyvista data structure
pdata = pyvista.PolyData(coord_array)
pdata['point_color'] = color_array

# create the point cloud with many shperes from geoms
pc = pdata.glyph(geom=geoms)

# plot with colors
pc.plot(scalars='point_color')

#"""