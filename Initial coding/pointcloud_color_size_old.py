"""
This file was created by Kristof Czirjak as part of his Bachelor's Thesis
"""

from Bio.PDB import *
from Bio import SeqIO
import numpy as np
import pyvista as pv


# Creates a dictionary of atomic coordinates by atom type
def load_atoms(file_name):

    # load atoms into a python list
    parser = PDBParser() 
    load_pdb = file_name + ".pdb"
    structure = parser.get_structure(file_name, load_pdb) 
    atoms = structure.get_atoms()
    atoms_list = list(atoms)

    # initialize dictionary, then fill up iteratively
    atom_data = dict()
    for atom in atoms_list:
        type_name = atom.element
        if type_name not in atom_data:
            atom_data[type_name] = [atom.get_coord()]
        else:
            atom_data[type_name].append(atom.get_coord())

    # return the 3D positions of atoms organized by atom type        
    return atom_data


# preprocessing for plotting by atom type
def pre_plot_atoms(atom_data):

    # Set resolution variables so only need to change once
    phi_res = 10
    theta_res = 10
    # these two dictionaries have to be manually created
    # also, if more/other atoms present in protein it will not work
    atoms_size_dict = {'C': pv.Sphere(radius=0.67, phi_resolution=phi_res, theta_resolution=theta_res),
                    'O': pv.Sphere(radius=0.48, phi_resolution=phi_res, theta_resolution=theta_res),
                    'N': pv.Sphere(radius=0.56, phi_resolution=phi_res, theta_resolution=theta_res),
                    'S': pv.Sphere(radius=0.88, phi_resolution=phi_res, theta_resolution=theta_res),
                    'CA': pv.Sphere(radius=1.94, phi_resolution=phi_res, theta_resolution=theta_res)}
    color_dict = {'C': '#333333',
                    'O': '#EA2128',
                    'N': '#007DEE',
                    'S': '#FFDC61',
                    'CA': '#38761D'}

    # create glyphs (spherical) to represent each atom
    atoms_spheres = []
    colors_spheres = []
    for atoms_type in atom_data:

        # create a mesh with each atoms position
        mesh = pv.PolyData(np.array(atom_data[atoms_type]))

        # place a specific sphere at given position
        glyphs = mesh.glyph(geom=atoms_size_dict[atoms_type])
        atoms_spheres.append(glyphs)

        #color the sphere according to 'CPK' standard
        colors_spheres.append(color_dict[atoms_type])

    # return list of spheres and colors representing each atom
    return atoms_spheres, colors_spheres

# plot point cloud of atoms
def plot_atoms(atoms_spheres, colors_spheres):
    # Use 3 lights because it's a bit brighter
    pl = pv.Plotter(lighting=None)
    pl.background_color = 'grey'
    pl.enable_3_lights()

    # we have to add atoms types one at a type since it seems we can only
    # apply one texture per atoms type.  This may complicate things with complicated or mixed atoms types.
    j = 0
    for mesh in atoms_spheres:
        pl.add_mesh(mesh, color=colors_spheres[j])
        j=j+1

    # this returns  camera position that you can use later for automated plots
    # here, we also save a png screenshot
    cpos = pl.show(screenshot='test.png')

def main():
    
    atom_data = load_atoms("2fd7")
    atoms, colors = pre_plot_atoms(atom_data)
    plot_atoms(atoms, colors)


if __name__ == "__main__":
    main()

# https://github.com/pv/pv-support/issues/374