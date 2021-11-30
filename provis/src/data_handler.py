from Bio.PDB import *
from Bio import SeqIO
import pyvista as pv
import numpy as np
from os.path import exists
import pandas as pd
from biopandas.mol2 import PandasMol2


from provis.utils.atminfo import import_atm_size_info, import_atm_mass_info
from provis.utils.bond_parser import bond_parser
from provis.src.file_converter import FileConverter


class DataHandler:
    """
    Class that does most computation.
    """
    def __init__(self, name):
        """
        Load structure form pdb file (by parsing file) and return structure. Load all size and color helper dictionaries.
        
        :param name: file_name - Name of file (without filetype) to be loaded
        :param type: str
        
        :return: structure - Biopython structure representation of the file
        """
        self._name = name
        parser = PDBParser()
        file_name = name + ".pdb"
        self._structure = parser.get_structure(name, file_name)
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
        
        :param name: res - Residue number of specified residue be looked at
        :param type: str
        :param name: chain - Chain number of corresponding to residue be looked at
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
        
    def get_structure(self):
        """
        Return the loaded structure object
        
        :return: structure
        """
        return self._structure
        
    
    def get_atom_mesh(self, atom_data, vw=0, probe=0, phi_res=10, theta_res = 10):
        """
        Create a list of Shperes and colors representing each atom for plotting. Can later be added to a mesh for plotting.
        
        :param name: atom_data - Dictionary of atoms and their coordinates, by atom type
        :param type: dict
        :param name: vw, optional - When set to True Van-der-Waals atomic radii used instead of empirical radii; default false
        :param type: bool
        :param name: probe - size of probe (representing the solvent size) needed for surface calculation
        :param type: int
        :param name: phi_res - pyvista phi_resolution for Sphere objects representing atoms
        :param type: int
        :param name: theta_res - pyvista theta_resolution for Sphere objects representing atoms
        :param type: int
        
        :return: list - List of pyvista Shperes representing each atom
        :return: list - List of colors for each atom
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


    def get_residue_mesh(self, res_data, phi_res = 25, theta_res = 25):
        """
        Create a list of Shperes and colors representing each residue for plotting
        
        :param name: res_data - Dictionary of residues and their coordinates, by residue type
        :param type: dict
        :param name: phi_res - pyvista phi_resolution for Sphere objects representing atoms
        :param type: int
        :param name: theta_res - pyvista theta_resolution for Sphere objects representing atoms
        :param type: int
        
        :return: list - List of pyvista Shperes representing each residue
        :return: list - List of colors for each residue
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
        Determine bonds from 3D information
        
        :return: list - List of pyvista lines representing each bond
        """
        fname = self._name + ".mol2"

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

