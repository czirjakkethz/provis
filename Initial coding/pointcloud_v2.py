
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
# Create list for radius
radius_array = [None] * num_atoms
# Create list for color
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
    radius = 0
    if element == "C":
        radius = 0.67
        geoms[i] = carbon
        color_array[i] = element #'#333333'
    elif element == "O":
        radius = 0.48
        geoms[i] = oxygen
        color_array[i] = element #'#EA2128'
    elif element == "N":
        radius = 0.56
        geoms[i] = nitrogen
        color_array[i] = element #'#007DEE'
    elif element == "S":
        radius = 0.88
        geoms[i] = sulfur
        color_array[i] = element #'#FFDC61'
    elif element == "CA":
        radius = 1.94
        geoms[i] = calcium
        color_array[i] = element #'#38761D'
    radius_array[i] = radius
    i = i+1


# Can not figure out coloring, probably have to add with each, individually
# create pyvista data structure, mesh(?)
pdata = pyvista.PolyData(coord_array)
pdata['point_color'] = color_array


# create pyvista plotter
pl = pyvista.Plotter()

# create the point cloud with geoms
pc = pdata.glyph(geom=geoms)
pl.add_mesh(pc, specular=1, specular_power=150,
                 smooth_shading=True, show_scalar_bar=True, scalars='point_color', clim=[-1,20])

pl.show()

