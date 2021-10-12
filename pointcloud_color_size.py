"""
This file was created by Kristof Czirjak as part of his Bachelor's Thesis
"""

from Bio.PDB import *
from Bio import SeqIO
from Thesis_l.msms_script import msms_script
import numpy as np
import pyvista as pv
import os
import trimesh
import pdb_to_xyzr
import msms_script

class FileConverter():
    def __init__(self, *args):
        """
        Can be called empty. If called with arguments conversion instant
        :param name: args[0] - name of file
        :param type: str
        :param name: args[1] - density of triangles
        :param type: float
        :return: void - converted files
        """
        if len(args) > 0:
            self._name = args[0]
            self._dest = args[1]
            self.pdb_to_xyzr(self._name, self._dest)
            self.msms(self._name)
        pass

    def pdb_to_xyzr(name, dest):
        """
        Run the pdb_to_xyzr script for given filename
        Works only on linux
        (could easily be extended with other parameters)
        :param name: name - Name of file
        :param type: str
        :return: void - xyzr file
        """
        pdb_to_xyzr.pdb_to_xyzr_script(name, dest)

    def msms(name):
        """
        Run the msms script for given filename and density
        Works on both linux and windows
        (could easily be extended with other parameters)
        :param name: name - Name of file
        :param type: str
        :param name: dens - Density of triangulation
        :param type: float
        :return: void - face and vert files
        """
        msms_script.msms_script(name)


class StickPoint:
    def __init__():
        pass

    def load_structure(file_name):
        """
        load structure form pdb file (by parsing file) and return structure
        :param name: file_name - Name of file (without filetype) to be loaded
        :param type: str
        :return: structure - Biopython structure representation of the file
        """

        parser = PDBParser() 
        load_pdb = file_name + ".pdb"
        structure = parser.get_structure(file_name, load_pdb) 
        
        return structure

    def get_atoms(structure, show_solvent=False):
        """
        Creates a dictionary of atomic coordinates from structure
        :param name: structure - Biopython structure object of desired molecule
        :param type: structure
        :param name: show_solvent - If True solvent molecules also added to retrun dictionary
        :param type: bool
        :return: dict - Dictionary of atomic coordinates by atom type
        """
        # load file into a python list
        residues = structure.get_residues()
        residues_list = list(residues)
        
        # initialize dictionary, then fill up iteratively
        atom_data = dict()
        for residue in residues_list:
            if residue.get_resname() == "HOH" and not show_solvent:
                continue
            for atom in residue:
                type_name = atom.element
                if type_name not in atom_data:
                    atom_data[type_name] = [atom.get_coord()]
                else:
                    atom_data[type_name].append(atom.get_coord())

        # return the 3D positions of atoms organized by atom type        
        return atom_data


    def get_residues(structure):
        """
        Creates a dictionary of coordinates by residues from structure object
        :param name: structure - Biopython structure object of desired molecule
        :param type: structure
        :return: dict - Dictionary of atomic coordinates by residue type
        """

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


    def pre_plot_atoms(atom_data):
        """
        Create a list of Shperes representing each atom for plotting
        :param name: atom_data - Dictionary of atoms and their coordinates, by atom type
        :param type: dict
        :return: list - List of pyvista Shperes representing each atom 
        :return: list - List of colors for each atom
        """
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


    def pre_plot_residues(res_data):
        """
        Create a list of Shperes representing each residue for plotting
        :param name: res_data - Dictionary of residues and their coordinates, by residue type
        :param type: dict
        :return: list - List of pyvista Shperes representing each residue 
        :return: list - List of colors for each residue
        """

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

    def pre_plot_bonds(struct):
        """
        Determine bonds from 3D information
        :param name: struct - Biopython structure object of desired molecule
        :param type: structure
        :return: list - List of pyvista lines representing each bond 
        """
        # create a list of non-solvent atoms
        residues = struct.get_residues()
        residues_list = list(residues)
        atom_list = []
        for residue in residues_list:
            if residue.get_resname() == "HOH":
                continue
            for atom in residue:
                atom_list.append(atom)

        ns = NeighborSearch(atom_list)
        _cutoff_dist = 3 # not sure how to choose this number

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

    def plot_atoms(atoms, col_s, box=0, bonds=0, res=0, col_r=0):
        """
        Plot stick and point model
        :param name: atoms - List of pyvista Shperes representing each atom 
        :param type: list
        :param name: col_s - List of colors for each atom
        :param type: list
        :param name: box - If True bounding box also visualized, default 0
        :param type: bool        
        :param name: bonds - List of pyvista Lines representing each bond, default 0 
        :param type: list
        :param name: res - List of pyvista Shperes representing each residue, default 0 
        :param type: list
        :param name: col_r - List of colors for each residue, default 0
        :param type: list
        :return: void - Window with interactive plot
        """
        # Use 3 lights because it's a bit brighter
        pl = pv.Plotter(lighting=None)
        pl.background_color = 'grey'
        pl.enable_3_lights()

        # adding the spheres (by atom type) one at a time
        j = 0
        for mesh in atoms:
            pl.add_mesh(mesh, color=col_s[j])
            j=j+1

        # adding the bonds one at a time
        if bonds != 0:
            l = 0
            for b in bonds:
                pl.add_mesh(b, color='w', render_lines_as_tubes=True, line_width=5)
        # adding the spheres (by residue) one at a time
        # only executes if residue information provided
        if res != 0:
            k = 0
            for mesh in res:
                pl.add_mesh(mesh, color=col_r[k], opacity=0.2)
                k+=1
        
        # if specified add bounding box
        if box:
            pl.add_bounding_box(color='white', corner_factor=0.5, line_width=1)

        # save a screenshot
        pl.show(screenshot='test.png')

class Surface:
    def __init__():
        pass

    def load(self, file_name):
        """
        Load surface information from both face and vert files
        :param name: file_name - Name of input files
        :param type: str
        :return: list - list of face data
        :return: list - list of vert data
        """
        
        face = self.load_fv(file_name, ".face", 0)
        vert = self.load_fv(file_name, ".vert", 1)
        return face, vert


    def load_fv(file_name, end, vorf):
        """
        Load surface information from face or vert file
        :param name: file_name - Name of input file
        :param type: str
        :param name: end - Type of input file
        :param type: str
        :param name: vorf - Specify face or vertex
        :param type: bool
        :return: list - list of data
        """
        
        outfile = open(file_name + end,"r")
        data = outfile.readlines()
        l3 = str.split(data[2])
        numlines = int(l3[0])
        numspheres = int(l3[1])
        density = float(l3[2])
        probe = float(l3[3])
        ret = [[] for x in range(numlines)]
        i = 0
        for line in data[3:]:
            line_split = str.split(line)
            k = 0
            for entry in line_split[:3]:
                if vorf:
                    ret[i].append(float(entry))
                else:
                    ret[i].append(int(entry)-1)
            i+=1

        outfile.close()
        return ret

    def plot_surface(self, filename):
        """
        Plot surface from face and vert files
        :param name: self - The surface object
        :param type: Surface
        :param name: file_name - Name of input file
        :param type: str
        :return: void - plot
        """

        face, vertice = self.load(self, filename)
        vertices = np.array(vertice)
        faces = np.hstack(face)
        pl = pv.Plotter(lighting=None)
        pl.background_color = 'grey'
        pl.enable_3_lights()
        tmesh = trimesh.Trimesh(vertice, faces=face, process=False)
        mesh = pv.wrap(tmesh)
        pl.add_mesh(mesh)
        pl.show()


def main():
    s = Surface
    s.plot_surface(s, "2fd7_fine")

    sp = StickPoint
    struct = sp.load_structure("2fd7")
    atom_data = sp.get_atoms(struct) # second arg: 1 = showsolvent
    bonds = sp.pre_plot_bonds(struct)
    res_data = sp.get_residues(struct)
    atoms, colors = sp.pre_plot_atoms(atom_data)
    ress, colors_r = sp.pre_plot_residues(res_data)
    sp.plot_atoms(atoms, colors, 1, bonds, ress, colors_r)


if __name__ == "__main__":
    main()

# https://github.com/pv/pv-support/issues/374