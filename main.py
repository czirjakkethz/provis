"""
This file was created by Kristof Czirjak as part of his Bachelor's Thesis
"""

from Bio.PDB import *
from Bio import SeqIO
import numpy as np
import pyvista as pv
import os
import trimesh
from scripts.msms_script import msms_script
from scripts.pdb_to_xyzr_script import pdb_to_xyzr_script
import platform
from scripts.atminfo import import_atm_size_info, import_atm_mass_info
from biopandas.mol2 import PandasMol2
import re
import pandas as pd
from os.path import exists
import subprocess


def bond_parser(filename):
    f = open(filename,'r')
    f_text = f.read()
    f.close()
    bond_start = f_text.find('@<TRIPOS>BOND')
    bond_end = f_text[bond_start:].replace('@<TRIPOS>BOND','').find('@')
    l = int(len(f_text[bond_start:bond_end].replace('@<TRIPOS>BOND\n','').replace('\n',' ').strip().split())/4)
    df_bonds = pd.DataFrame(np.array(f_text[bond_start:bond_end].replace('@<TRIPOS>BOND\n','').replace('\n',' ').strip().split()).reshape((l,4)), # 
            columns=['bond_id', 'atom1', 'atom2', 'bond_type'])
    df_bonds.set_index(['bond_id'], inplace=True)
    return df_bonds


class FileConverter():
    def __init__(self, *args):
        """
        Can be constructed empty. If called with arguments conversion instant
        :param name: args[0] - name of file
        :param type: str
        :param name: args[1] - density of triangles
        :param type: float
        :return: void - converted files
        """
        if len(args) > 0:
            self._name = args[0]
            self._dest = args[1]
            self._solv = args[2]
            self._bash = args[3]
            self.pdb_to_xyzr(self._name, self._solv, self._bash)
            self.msms(self._name, self._dest)
            self.pdb_to_mol2(self._name)
        pass

    def pdb_to_xyzr(self, name, solvent, bash=0):
        """
        Run the pdb_to_xyzr script for given filename
        Works only on linux
        (could easily be extended with other parameters)
        :param name: name - Name of file
        :param type: str
        :return: void - xyzr file
        """
        pdb_to_xyzr_script.pdb_to_xyzr_script(name, solvent, bash)

    def msms(self, name, dest):
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
        msms_script.msms_script(name, dest)

    def pdb_to_mol2(self, name):
        ac = subprocess.call("obabel %4s.pdb -O %4s.mol2" % (name, name), shell=True)


class DataHandler:
    def __init__(self, name):
        """
        load structure form pdb file (by parsing file) and return structure
        :param name: file_name - Name of file (without filetype) to be loaded
        :param type: str
        :return: structure - Biopython structure representation of the file
        """
        self._name = name
        parser = PDBParser()
        file_name = name + ".pdb"
        structure = parser.get_structure(name, file_name)
        
        self._structure = structure

    def get_atoms(self, show_solvent=False):
        """
        Creates a dictionary of atomic coordinates from structure
        :param name: show_solvent - If True solvent molecules also added to retrun dictionary
        :param type: bool
        :return: dict - Dictionary of atomic coordinates by atom type
        """
        # load file into a python list
        residues = self._structure.get_residues()
        residues_list = list(residues)
        
        # initialize dictionary, then fill up iteratively
        atom_data = dict()
        for residue in residues_list:
            if residue.get_resname() == "HOH" and not show_solvent:
                continue
            for atom in residue:
                type_name = atom.element
                # If atom not in dictionary, add it as key with coords in list
                if type_name not in atom_data:
                    atom_data[type_name] = [atom.get_coord()]
                # If atom already in dictionary, append its coordinates to list
                else:
                    atom_data[type_name].append(atom.get_coord())

        # return the 3D positions of atoms organized by atom type        
        return atom_data


    def get_residues(self):
        """
        Creates a dictionary of coordinates by residues from structure object
        :return: dict - Dictionary of atomic coordinates by residue type
        """

        # load file into a python list
        residues = self._structure.get_residues()
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

    def get_residue_info(self, res, chain, option):
        """
        Calculates information about specified residue from mol2 file
        If exact residue specified eg. "THR1" only that COM returned
        Else list of COM of all residues of that name returned TODO
        :param name: res - Residue name to be looked at
        :param type: str
        :param name: option - choose what property of residue you want
        :param type: str - com for Centre Of Mass, ch for charge
        :return: list - list of COM coords of given (exact) residue
        """
        
        # load info for given residue
        fname = self._name + ".mol2"
        pmol = PandasMol2().read_mol2(fname)
        my = pmol.df[pmol.df['subst_id'] == res]
        
        if len(my) == 0:
            print("ERROR: No such residue: ", res)
            return 1
        
        fid = my.iloc[0]['atom_id'] -1
        curr_chain = 0
        actual = pd.DataFrame(columns = list(my))
        for i, row in my.iterrows():
            if row['atom_id'] != fid + 1:
                curr_chain += 1
            if curr_chain == chain:
                actual = actual.append(row)
            fid = row['atom_id']
                
        # if com option calculate Centre of mass
        if option.upper() == "COM":
            COM_sum = [0.0, 0.0, 0.0]
            mass_dict = import_atm_mass_info()
            a = 0
            for i in range(len(actual)):
                m = mass_dict[actual.iloc[i]['atom_type'][0]]
                COM_sum[0] += actual.iloc[i]['x'] * m
                COM_sum[1] += actual.iloc[i]['y'] * m
                COM_sum[2] += actual.iloc[i]['z'] * m
                a += m
            
            COM = [i/a for i in COM_sum]
            return COM
            
        # if charge option calculate residue charge
        elif option.upper() == "CH" or option.upper() == "CHARGE":
            charge = 0.0
            for i in range(len(my)):
                charge += my.iloc[i]['charge']

            return charge
    
    def load_fv(self, file_name, end, vorf):
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
        
    def get_structure(self):
        """
        Return the loaded structure object
        """
        return self._structure


class StickPoint:
    def __init__(self, file_name):
        self._name = file_name
        self._dh = DataHandler(file_name)
        # these two dictionaries have to be manually created
        # could not find volume information for CYS, HIS, LYS, THR, TYR
        # for above mentioned apporximation by "closest" available
        self._res_size_dict = { "ALA": 2.80, "ARG": 3.77, "ASN": 3.18, "ASP": 1.51,
            "CYS": 2.87, "GLN": 3.37, "GLU": 1.68, "GLY": 2.51, "HIS": 3.65,
            "ILE": 3.43, "LEU": 3.42, "LYS": 3.77, "MET": 3.44, "PHE": 3.65,
            "PRO": 3.13, "SER": 2.87, "THR": 2.87, "TRP": 3.54, "TYR": 3.65,
            "VAL": 3.24, "HOH":1.375, "HEM": 33.43 } #HEM vol 156536


        self._res_color_dict = {"ALA": '#13B6E2', "ARG": '#23DEFE', "ASN": '#0EFF57', "ASP": '#822DD2',
            "CYS": '#A22282', "GLN": '#60C7B0', "GLU": '#6913FE', "GLY": '#8E41D0', "HIS": '#7089FC',
            "ILE": '#0CAFF9', "LEU": '#806769', "LYS": '#8B7928', "MET": '#68D139', "PHE": '#8BA695',
            "PRO": '#9FEBA4', "SER": '#BBD7EB', "THR": '#D1A67A', "TRP": '#F93713', "TYR": '#E5613D',
            "VAL": '#128033', "HOH": 'w', "HEM": 'r'}
        self._atoms_size_dict, self._atoms_color_dict, self._vw_dict = import_atm_size_info(1)



    def pre_plot_atoms(self, atom_data, vw=0, probe=0):
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

        # create glyphs (spherical) to represent each atom
        atoms_spheres = []
        colors_spheres = []
        rad = probe / 2.0
        for atoms_type in atom_data:

            # create a mesh with each atoms position
            mesh = pv.PolyData(np.array(atom_data[atoms_type]))
            # place a specific sphere at given position
            if not vw:
                sphere = pv.Sphere(radius=self._atoms_size_dict[atoms_type], phi_resolution=phi_res, theta_resolution=theta_res)
            else:
                sphere = pv.Sphere(radius=(self._vw_dict[atoms_type] + rad), phi_resolution=phi_res, theta_resolution=theta_res)
            glyphs = mesh.glyph(geom=sphere)
            atoms_spheres.append(glyphs)

            try:
                # color the sphere according to 'CPK' standard
                colors_spheres.append(self._atoms_color_dict[atoms_type])
            except KeyError:
                # if atom is unrecognized, color it pink
                colors_spheres.append('#DD77FF')

        # return list of spheres and colors representing each atom
        return atoms_spheres, colors_spheres


    def pre_plot_residues(self, res_data):
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

        # create glyphs (spherical) to represent each res
        res_spheres = []
        colors_spheres = []
        for res_type in res_data:

            sphere = pv.Sphere(radius=0.0, phi_resolution=phi_res, theta_resolution=theta_res)
            try:    
                #color the sphere according to 'CPK' standard
                colors_spheres.append(self._res_color_dict[res_type])

                sphere = pv.Sphere(radius=self._res_size_dict[res_type], phi_resolution=phi_res, theta_resolution=theta_res)
            except KeyError:
                colors_spheres.append('#DD77FF')



            # create a mesh with each residues position
            mesh = pv.PolyData(np.array(res_data[res_type]))

            # place a specific sphere at given position
            glyphs = mesh.glyph(geom=sphere)
            res_spheres.append(glyphs)


        # return list of spheres and colors representing each res
        return res_spheres, colors_spheres

    def pre_plot_bonds_old(self, struct):
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

    def pre_plot_bonds(self, name):
        """
        Determine bonds from 3D information
        :param name: fname - file name without extension
        :param type: structure
        :return: list - List of pyvista lines representing each bond 
        """
        fname = name + ".mol2"

        # Check if mol2 file exists. If not convert it from pdb
        file_exists = exists(fname)
        fc = FileConverter()
        if not file_exists:
            fc.pdb_to_mol2(name)
        pmol = PandasMol2().read_mol2(fname)
        bonds_in = bond_parser(fname)

        # Find endpoints of bonds, create Line() and add to list of lines
        bonds = []
        for i in range(len(bonds_in)):    
            my = pmol.df.iloc[int(bonds_in.iloc[i]['atom1'])-1][['x', 'y', 'z']]
            c = pmol.df.iloc[int(bonds_in.iloc[i]['atom2'])-1][['x', 'y', 'z']]
            
            my_coord = my.values.tolist()
            c_coord = c.values.tolist()
            if bonds_in.iloc[i]['bond_type'] == '1':
                line = pv.Line(my_coord, c_coord, resolution=5)
            bonds.append(line)
        
        return bonds

    def plot_stick_ball(self, atoms, col_s, box=0, bonds=0, vw=0, residues=0, col_r=0, res=0):
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
        opacity = 1 - vw*0.4
        style = 'surface'
        if vw:
            style='wireframe'
        for j, mesh in enumerate(atoms):
            pl.add_mesh(mesh, color=col_s[j], opacity=opacity, smooth_shading=True, style=style)

        # adding the bonds one at a time
        if bonds != 0:
            for b in bonds:
                pl.add_mesh(b, color='w', render_lines_as_tubes=True, line_width=5)
        # adding the spheres (by residue) one at a time
        # only executes if residue information provided
        if residues != 0:
            for k, mesh in enumerate(residues):
                pl.add_mesh(mesh, color=col_r[k], opacity=0.2)
        
        # if specified add bounding box
        if box:
            pl.add_bounding_box(color='white', corner_factor=0.5, line_width=1)
        if res:
            res_list, chain_list, pad = res.get_res_info()
            for i, r in enumerate(res_list):
                chain = chain_list[i]
                residues = self._dh.get_structure().get_residues()
                residues_list = list(residues)
                res_name = residues_list[r + 1].get_resname()
                d = (self._res_size_dict[res_name] + pad) * 2
                x, y, z = d,d,d
                pl.add_mesh(pv.Cube(center=self._dh.get_residue_info(r, chain,'com'), x_length=x, y_length=y, z_length=z), style='wireframe', show_edges=1, line_width=5)
        # save a screenshot
        pl.show(screenshot='data/test.png')

class Surface:
    def __init__(self, name):
        self._dh = DataHandler(name)
        

    def load(self, file_name):
        """
        Load surface information from both face and vert files
        :param name: file_name - Name of input files
        :param type: str
        :return: list - list of face data
        :return: list - list of vert data
        """
        
        face = self._dh.load_fv(file_name, ".face", 0)
        vert = self._dh.load_fv(file_name, ".vert", 1)
        return face, vert


    def plot_surface(self, filename):
        """
        Plot surface from face and vert files
        :param name: self - The surface object
        :param type: Surface
        :param name: file_name - Name of input file
        :param type: str
        :return: void - plot
        """
        face, vertice = self.load(filename)
        vertices = np.array(vertice)
        faces = np.hstack(face)
        pl = pv.Plotter(lighting=None)
        pl.background_color = 'grey'
        pl.enable_3_lights()
        tmesh = trimesh.Trimesh(vertice, faces=face, process=False)
        mesh = pv.wrap(tmesh)
        pl.add_mesh(mesh)
        pl.show()

    def new_surface(self, atoms, atmsrf, col):
        
#        by_type = list(atoms.values())
#        points_ = []
#        for i in range(len(by_type)):
#            points_.extend(by_type[i])
#        points = np.array(points_)
#        hull = Delaunay(points)
#        indices = hull.simplices
#        faces = np.hstack(indices)
#        vertices = points[indices]
#        vertices_ = vertices.squeeze()

        #############======================

        pl = pv.Plotter(lighting=None)
        pl.background_color = 'grey'
        pl.enable_3_lights()

        # adding the spheres (by atom type) one at a time
        j = 0
        style = 'surface'
        mesh_ = pv.wrap(atmsrf[0])
        for mesh in atmsrf[1:]:
            mesh_ = mesh_ + (mesh)
            
        vol = mesh_.delaunay_3d(alpha=1.4)
        shell = vol.extract_surface().reconstruct_surface()
        # Compute the normals in-place and use them to warp the globe
#        shell.compute_normals(inplace=True)  # this activates the normals as well
#
#        # Now use those normals to warp the surface
#        warp = shell.warp_by_scalar(factor=0.5e-5)
        pl.add_mesh(shell, color="white", smooth_shading=True, style=style, show_edges=True)#, culling='back')
#        # save a screenshot
##        pl.show(screenshot='test.png')
##
##        mesh = pv.wrap(atomss[0])
##        for r in range(1,len(atomss)):
##            print(str(r) + " is triangulated: " + str(atomss[r].is_all_triangles()))
##            mesh = mesh.boolean_union(atomss[r])
##            print(str(r) + " is done")
##        mesh.plot()
##        pl = pv.Plotter(lighting=None)
##        pl.background_color = 'grey'
##        pl.enable_3_lights()
##        for mesh in res:
##            pl.add_mesh(mesh, color='w')
        pl.show()
        
        ###########===========================
        
#from scipy.spatial import Delaunay, ConvexHull
#import numpy as np
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.tri as mtri
#from scipy.spatial import Delaunay
#import numpy as np
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.tri as mtri

#        # points is a 3D numpy array (n_points, 3) coordinates of a sphere
#        cloud = pv.PolyData(points)
#        volume = cloud.delaunay_3d(alpha=2.7)
#        shell = volume.extract_geometry()
#        shell.plot()
        
        
        
#        X = np.array(points)
#        # My data points are strictly positive. This doesn't work if I don't center about the origin.
#        X -= X.mean(axis=0)
#
#        rad = np.linalg.norm(X, axis=1)
#        zen = np.arccos(X[:,-1] / rad)
#        azi = np.arctan2(X[:,1], X[:,0])
#
#        tris = mtri.Triangulation(zen, azi)
#
#        fig = plt.figure()
#        ax  = fig.add_subplot(111, projection='3d')
#        ax.plot_trisurf(X[:,0], X[:,1], X[:,2], triangles=tris.triangles, cmap=plt.cm.bone)
#        plt.show()


#        # colorful mess, topological surface
#        from Bio.PDB.DSSP import DSSP
#        p = PDBParser()
#        structure = p.get_structure("2fd7", "data/2fd7.pdb")
#        model = structure[0]
#        dssp = DSSP(model, "data/2fd7.pdb")
#        a_key = list(dssp.keys())[2]
#
#        print(dssp[a_key])

#        fig = plt.figure()
#        ax = fig.add_subplot(1, 1, 1, projection='3d')
#        # The triangles in parameter space determine which x, y, z points are
#        # connected by an edge
#        ax.plot_trisurf(points[:,0], points[:,1], points[:,2], cmap=plt.cm.Spectral)
#        # ax.plot_trisurf(points[:,0], points[:,1], points[:,2], triangles=hull.simplices, cmap=plt.cm.Spectral)
#        plt.show()

         ##############JUST COORDS, Kinda works#################
#
#        pl = pv.Plotter(lighting=None)
#        pl.background_color = 'grey'
#        pl.enable_3_lights()
#         ###########According to actual######################
#
#        #mesh = pv.PolyData(new)
#        tmesh = trimesh.Trimesh(points, faces=indices, process=False)
#
##         ##############adjusted for some reason####################
##        new_v = []
##        for quads in vertices:
##            for point in quads:
##                new_v.append(point)
##        new_f = []
##        for el in range(len(indices)):
##            new_f.append([])
##            for xyz in range(len(indices[el])):
##                new_f[el].append(indices[el][xyz] - 1)
##
##        tmesh = trimesh.Trimesh(new_v, faces=new_f, process=False)
## ############################################
#
#        mesh = pv.wrap(tmesh)
#        pl.add_mesh(mesh)
#
#        pl.show(screenshot='surf_2.png')

class Residue():
    def __init__(self, id=None, chain=0, padding=0):
        if id == None:
            self._id_list = []
            self._chain_list = []
        else:
            self._id_list = [id]
            self._chain_list = [chain]
        self._pad = padding
        
    def add_residue(self, id, chain=0):
        """
        :param name: id - id of residue
        :param type: int
        :param name: chain - chain id of residue
        :param type: int
        :return: list - list of the current residues
        """
        self._id_list.append(id)
        self._chain_list.append(chain)
    
    def get_res_info(self):
        """
        :return: list - list of the current residues
        :return: int - padding for bounding box
        """
        return self._id_list, self._chain_list, self._pad
    
                                    
def main():
    name = "data/2fd7" # "data/1a3n" #
    density = 3.0
    solvent = 0
    bash = 0
    show_box = 1

    dh = DataHandler(name)
    ## for given residue (ID) calculate COM or Charge
#   (dh.get_residue_info(3,'com'))
    # TODO: change from name to residue id DONE
    # TODO: add bounding box to residue DONE
    # TODO: use Charlotte surface color code
    
    ## plot stick point
    sp = StickPoint(name)
    atom_data = dh.get_atoms() # second arg: 1 = showsolvent
    ## return list of spheres (meshes) and colors for the spheres
    atoms, colors = sp.pre_plot_atoms(atom_data, vw=0) # second arg: 1 = showvw spheres instead of "normal" radius
    ## return list of lines (meshes)
#    bonds = sp.pre_plot_bonds(name)
    ## return list of spheres (meshes) and colors for the spheres
    res_data = dh.get_residues()
    ress, colors_r = sp.pre_plot_residues(res_data)
    ## Create residue for plotting with bounding box
    r = Residue(6)
    ## plot stick and balls model
    sp.plot_stick_ball(atoms=atoms, col_s=colors, box=0, vw=0, bonds=0, residues=0, col_r=0, res=r)
#
#
#    ## plot surface
#    fc = FileConverter(name, density, solvent, bash)
#    s = Surface(name)
#    out_name = name + "_out_" + str(int(density))
#    s.plot_surface(out_name)
#    atmsurf, col = sp.pre_plot_atoms(atom_data, vw=1, probe=1)
#    s.new_surface(atom_data, atmsurf, col)

if __name__ == "__main__":
    main()

# https://github.com/pv/pv-support/issues/374
