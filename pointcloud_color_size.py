"""
This file was created by Kristof Czirjak as part of his Bachelor's Thesis
"""

from Bio.PDB import *
from Bio import SeqIO
import numpy as np
import pyvista as pv

# load structure form pdb file
def load_structure(file_name):

    # parse file and return structure
    parser = PDBParser() 
    load_pdb = file_name + ".pdb"
    structure = parser.get_structure(file_name, load_pdb) 
    return structure

# Creates a dictionary of atomic coordinates by atom type

    # load file into a python list
    # initialize dictionary, then fill up iteratively
    atom_data = dict()

    # return the 3D positions of atoms organized by atom type        


# Creates a dictionary of coordinates by residues
def get_residues(structure):

    # load file into a python list
    residues = structure.get_residues()
    res_list = list(residues)

    # initialize dictionary, then fill up iteratively
    res_data = dict()
    for res in res_list:
        type_name = res.get_resname()
        center_coord = [0., 0., 0.]
        c = 0.
        for atom in res:    
            center_coord = center_coord + atom.get_coord()
            c += 1.
        center_coord = center_coord / c
        if type_name not in res_data:
            res_data[type_name] = [center_coord]
        else:
            res_data[type_name].append(center_coord)

    # return the 3D positions of atoms organized by atom type        
    return res_data


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


# preprocessing for plotting by atom type
def pre_plot_residues(res_data):

    # Set resolution variables so only need to change once
    phi_res = 25
    theta_res = 25
    # these two dictionaries have to be manually created
    # could not find volume information for CYS, HIS, LYS, THR, TYR
    # for above mentioned apporximation by "closest" available
    res_size_dict = { "ALA": 2.80, "ARG": 3.77, "ASN": 3.18, "ASP": 1.51, 
        "CYS": 2.87, "GLN": 3.37, "GLU": 1.68, "GLY": 2.51, "HIS": 3.65,
        "ILE": 3.43, "LEU": 3.42, "LYS": 3.77, "MET": 3.44, "PHE": 3.65, 
        "PRO": 3.13, "SER": 2.87, "THR": 2.87, "TRP": 3.54, "TYR": 3.65, 
        "VAL": 3.24, "HOH":0}


    color_dict = {"ALA": '#13B6E2', "ARG": '#23DEFE', "ASN": '#0EFF57', "ASP": '#822DD2', 
        "CYS": '#A22282', "GLN": '#60C7B0', "GLU": '#6913FE', "GLY": '#8E41D0', "HIS": '#7089FC',
        "ILE": '#0CAFF9', "LEU": '#806769', "LYS": '#8B7928', "MET": '#68D139', "PHE": '#8BA695', 
        "PRO": '#9FEBA4', "SER": '#BBD7EB', "THR": '#D1A67A', "TRP": '#F93713', "TYR": '#E5613D', 
        "VAL": '#128033', "HOH": 'w' }

    # create glyphs (spherical) to represent each res
    res_spheres = []
    colors_spheres = []
    for res_type in res_data:

        shere = pv.Sphere(radius=res_size_dict[res_type], phi_resolution=phi_res, theta_resolution=theta_res)

        # create a mesh with each residues position
        mesh = pv.PolyData(np.array(res_data[res_type]))

        # place a specific sphere at given position
        glyphs = mesh.glyph(geom=shere)
        res_spheres.append(glyphs)

        #color the sphere according to 'CPK' standard
        colors_spheres.append(color_dict[res_type])

    # return list of spheres and colors representing each res
    return res_spheres, colors_spheres

def pre_plot_bonds(atom_list):
    ns = NeighborSearch(atom_list)

    bonds = []
    for target in atom_list:
        close_atoms = ns.search(target.coord, _cutoff_dist)
        my_coord = target.get_coord()
        for close_atom in close_atoms:
            # cylinder = pyvista.Cylinder(center=[1, 2, 3], direction=[1, 1, 1], 
            #     radius=1, height=2)
            c_coord = close_atom.get_coord()
            line = pv.Line(my_coord, c_coord, resolution=5)
            bonds.append(line)
    
    return bonds

# plot point cloud of atoms
    # Use 3 lights because it's a bit brighter
    pl = pv.Plotter(lighting=None)
    pl.background_color = 'grey'
    pl.enable_3_lights()

    # adding the spheres (by atom type) one at a time
    j = 0
    for mesh in atoms_spheres:
        pl.add_mesh(mesh, color=colors_spheres[j])
        j=j+1

    # adding the bonds one at a time
    if bonds != 0:
        l = 0
        for b in bonds:
            pl.add_mesh(b, color='w', render_lines_as_tubes=True, line_width=5)
    # adding the spheres (by residue) one at a time
    # only executes if residue information provided
        k = 0
        for mesh in r:
            pl.add_mesh(mesh, color=c[j], opacity=0.2)
            k+=1
    
    # if specified add bounding box
    if box:
        pl.add_bounding_box(color='white', corner_factor=0.5, line_width=1)

    # save a screenshot
    pl.show(screenshot='test.png')

def main():
    
    struct = load_structure("2fd7")
    res_data = get_residues(struct)
    atoms, colors = pre_plot_atoms(atom_data)
    ress, colors_r = pre_plot_residues(res_data)


if __name__ == "__main__":
    main()

# https://github.com/pv/pv-support/issues/374