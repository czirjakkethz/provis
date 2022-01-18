from os import fchdir
from Bio.PDB import *
from Bio import SeqIO
import pyvista as pv
import numpy as np
from os.path import exists
import pandas as pd
from biopandas.mol2 import PandasMol2
import trimesh

from provis.utils.atminfo import import_atm_size_info, import_atm_mass_info
from provis.utils.bond_parser import bond_parser
from provis.src.processing.name_checker import NameChecker
from provis.src.processing.file_converter import FileConverter


class DataHandler:
    """
    The 'brain' of provis, when it comes to handling atomic positions. 
    
    This class loads information from a variety of files and creates meshes to be plotted. 
    Upper level classes - eg. Structure - have their own DataHandler objects that do all the work.
    """
    def __init__(self, nc, fc=None):
        """
        Load structure form pdb file (by parsing file) and save Biopython structure representation of the protein. Also loads all size and color helper dictionaries.
        
        :param name: nc - Instance of a NameChecker class. Used to pass the pdb file name and paths.
        :param type: NameChecker
        :param name: fc - Instance of a FileConverter class. Needed if temporary files have not been created before instantiating this class. Default: None.
        :param type: FileConverter, optional
        """
        self._path, self._out_path, self._base_path = nc.return_all()
        if not fc:
            self._fc = FileConverter(nc)
        else:
            self._fc = fc
        
        self._fc.pdb_to_xyzrn(self._path, self._out_path)
        
        parser = PDBParser()
        file_name = self._path + ".pdb"
        self._structure = parser.get_structure(self._path, file_name)
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
        self._atoms_size_dict, self._atoms_color_dict, self._vw_dict = import_atm_size_info(self._base_path, 1)


    def get_atoms(self, show_solvent=False):
        """
        Creates a dictionary of atomic coordinates from structure
        
        :param name: show_solvent - If True solvent molecules also added to retrun dictionary. Default: False.
        :param type: bool, optional
        
        :returns: dict - Dictionary of atomic coordinates by atom type.
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
        
            
    def get_atoms_IDs(self):
        """
        Get atomic coordinates and residue IDs (format from output_pdb_as_xyzrn()) from the xyzrn file.

        :returns: dict - Dictionary of atomic coordinates by atom type.
        :returns: list - List of unique residue IDs (format from output_pdb_as_xyzrn())
        :returns: list - Atomic coordinates (in same order as the residue IDs)
        """
        xyzrnfile = open(self._out_path + ".xyzrn")
        meshdata = (xyzrnfile.read().rstrip()).split("\n")
        xyzrnfile.close()
        
        lenm = len(meshdata)
        res_id = [""] * lenm
        
        atom_data = dict()
        atom_coords = []
        for vi in range(lenm):
            fields = meshdata[vi].split()
            vertices = [None] * 3
            vertices[0] = float(fields[0])
            vertices[1] = float(fields[1])
            vertices[2] = float(fields[2])
            res_id[vi] = fields[5]
            res = res_id[vi].split("_")
            atmtype = res[4][0]
            if atmtype not in atom_data:
                atom_data[atmtype] = [vertices]
            # If atom already in dictionary, append its coordinates to list
            else:
                atom_data[atmtype].append(vertices)
            atom_coords.append(vertices)


        # return the 3D positions of atoms organized by atom type
        return atom_data, res_id, atom_coords

    def get_residues(self):
        """
        Creates a dictionary of coordinates by residues from structure object
        
        :returns: dict - Dictionary of atomic coordinates by residue type.
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
        
        :param name: res - Residue number of specified residue be looked at.
        :param type: str
        :param name: chain - Chain number of corresponding to residue be looked at.
        :param name: option - choose what property of residue you want. com for Centre Of Mass, ch for charge
        :param type: str - options: com, ch
        
        :returns: list - List of COM coords of given (exact) residue.
        """
        
        # load info for given residue
        fname = self._out_path + ".mol2"
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
        
    def get_structure(self):
        """
        Return the loaded structure object
        
        :returns: structure
        """
        return self._structure
        
    
    def get_atom_mesh(self, atom_data, vw=0, probe=0, phi_res=10, theta_res = 10):
        """
        Create a list of Shperes and colors representing each atom for plotting. Can later be added to a mesh for plotting.
        
        :param name: atom_data - Dictionary of atoms and their coordinates, by atom type.
        :param type: dict
        :param name: vw, optional - When set to True Van-der-Waals atomic radii used instead of empirical radii. Default: False.
        :param type: bool, optional
        :param name: probe - size of probe (representing the solvent size) needed for surface calculation. Default: 0.
        :param type: int, optional
        :param name: phi_res - pyvista phi_resolution for Sphere objects representing atoms. Default: 10.
        :param type: int, optional
        :param name: theta_res - pyvista theta_resolution for Sphere objects representing atoms. Default: 10.
        :param type: int, optional
        
        :returns: list - List of pyvista Shperes representing each atom
        :returns: list - List of colors corresponding to each atom
        """
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

    def get_atom_trimesh(self, atom_data, vw=False, probe=0):
        """
        Create a list of shperes and colors representing each atom for plotting in a Trimesh format. Used for feature computation in the surface_handler class.
        
        :param name: atom_data - Dictionary of atoms and their coordinates, by atom type.
        :param type: dict
        :param name: vw, optional - When set to True Van-der-Waals atomic radii used instead of empirical radii. Default: False.
        :param type: bool, optional
        :param name: probe - size of probe (representing the solvent size) needed for surface calculation. Default: 0.
        :param type: int, optional
        
        :returns: list - List of pyvista Shperes representing each atom
        :returns: list - List of colors for each atom
        """
        # these two dictionaries have to be manually created
        # also, if more/other atoms present in protein it will not work

        # create glyphs (spherical) to represent each atom
        atoms_spheres = []
        colors_spheres = []
        rad = probe / 2.0
        

        for atoms_type in atom_data:

            # create a mesh with each atoms position
            # mesh = trimesh.points.PointCloud()#np.array(atom_data[atoms_type]))
            # print(mesh)
            # place a specific sphere at given position

            spheres = []
            for atom in atom_data[atoms_type]:
                if not vw:
                    spheres.append(trimesh.primitives.Sphere(radius=self._atoms_size_dict[atoms_type], center=atom))
                else:
                    spheres.append(trimesh.primitives.Sphere(radius=(self._vw_dict[atoms_type] + rad), center=atom))

            mesh = trimesh.util.concatenate(spheres)
            atoms_spheres.append(mesh)

            try:
                # color the sphere according to 'CPK' standard
                colors_spheres.append(self._atoms_color_dict[atoms_type])
            except KeyError:
                # if atom is unrecognized, color it pink
                colors_spheres.append('#DD77FF')

        # return list of spheres and colors representing each atom
        return atoms_spheres, colors_spheres

    def get_residue_mesh(self, res_data, phi_res = 25, theta_res = 25):
        """
        Create a list of Shperes and colors representing each residue for plotting
        
        :param name: res_data - Dictionary of residues and their coordinates by residue type.
        :param type: dict
        :param name: phi_res - pyvista phi_resolution for Sphere objects representing atoms. Defaul: 25.
        :param type: int, optional
        :param name: theta_res - pyvista theta_resolution for Sphere objects representing atoms. Default:25.
        :param type: int, optional
        
        :returns: list - List of pyvista Shperes representing each residue
        :returns: list - List of colors for each residue
        """

        # create glyphs (spherical) to represent each res
        res_spheres = []
        colors_spheres = []
        for res_type in res_data:

            sphere = pv.Sphere(radius=0.0, phi_resolution=phi_res, theta_resolution=theta_res)
            try:
                sphere = pv.Sphere(radius=self._res_size_dict[res_type], phi_resolution=phi_res, theta_resolution=theta_res)
                
                #color the sphere according to 'CPK' standard
                colors_spheres.append(self._res_color_dict[res_type])
                
            except KeyError:
                # color unkown to light purple
                colors_spheres.append('#DD77FF')

            # create a mesh with each residues position
            mesh = pv.PolyData(np.array(res_data[res_type]))

            # place a specific sphere at given position
            glyphs = mesh.glyph(geom=sphere)
            res_spheres.append(glyphs)


        # return list of spheres and colors representing each res
        return res_spheres, colors_spheres

    def get_bond_mesh(self):
        """
        Determine bonds from 3D information.
        
        The color information is as follows:
            White for all single bonds, 
            Blue for all double bonds,
            Green for all triple bonds,
            Red for all amide bonds,
            Purple for all aromatic bonds,
            Black for everything else.
        
        :returns: list - List of pyvista lines representing each bond.
        :retruns: list - List of colors corresponding to the lines in the above list
        """
        fname = self._out_path + ".mol2"

        # Check if mol2 file exists. If not convert it from pdb
        file_exists = exists(fname)
        if not file_exists:
            self._fc.pdb_to_mol2(self._path, self._out_path)
            
        pmol = PandasMol2().read_mol2(fname)
        bonds_in = bond_parser(fname)

        # Find endpoints of bonds, create Line() and add to list of lines
        bonds = []
        col = []
        for i in range(len(bonds_in)):
            my = pmol.df.iloc[int(bonds_in.iloc[i]['atom1'])-1][['x', 'y', 'z']]
            c = pmol.df.iloc[int(bonds_in.iloc[i]['atom2'])-1][['x', 'y', 'z']]
            
            my_coord = my.values.tolist()
            c_coord = c.values.tolist()
            if bonds_in.iloc[i]['bond_type'] == '1': #TODO: differentiate between different bonds.
                line = pv.Line(my_coord, c_coord, resolution=5)
                col.append('w')
            elif bonds_in.iloc[i]['bond_type'] == '2':
                line = pv.Line(my_coord, c_coord, resolution=1)
                col.append('b')
            elif bonds_in.iloc[i]['bond_type'] == '3':
                line = pv.Line(my_coord, c_coord, resolution=1)
                col.append('g')
            # amide bond in red
            elif bonds_in.iloc[i]['bond_type'] == 'am':
                line = pv.Line(my_coord, c_coord, resolution=1)
                col.append('r')
            # aromatic bond in purple
            elif bonds_in.iloc[i]['bond_type'] == 'ar':
                line = pv.Line(my_coord, c_coord, resolution=1)
                col.append('m')
            else:
                line = pv.Line(my_coord, c_coord, resolution=1)
                col.append('k')
            bonds.append(line)
        
        return bonds, col
    
    def get_backbone_mesh(self):
        """
        Creates and returns a Spline object representing the backbone of the protein.
        
        :returns: numpy.ndarray - List of coordinates representing the centre of mass of each residue. To be used in stick_point to create a Spline.
        """
        
        residues = self._structure.get_residues()
        
        res_list = []
        for res in residues:
            if res.get_resname() != "HOH":
                com = [0, 0, 0]
                count = 0
                for atom in res: 
                    coords = atom.get_coord()
                    com += coords
                    count += 1
                com /= count
                res_list.append(com)
                
        return np.array(res_list)