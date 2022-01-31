from os import fchdir
from Bio.PDB import *
import pyvista as pv
import numpy as np
from os.path import exists
import pandas as pd
from biopandas.mol2 import PandasMol2
import trimesh
from biopandas.mol2 import split_multimol2
import copy

from provis.utils.atminfo import *
from provis.src.processing.name_checker import NameChecker
from provis.src.processing.file_converter import FileConverter


class DataHandler:
    """
    The 'brain' of provis, when it comes to handling atomic positions. 
    
    This class loads information from a variety of files and creates meshes to be plotted. 
    Upper level classes - eg. Protein - use DataHandler objects to create the meshes.
    
    The DataHandler class loads atom-positional information from a .pdb file and from this information computes the necessairy molecular structure mesh.
    It also loads pre-defined dictionaries from atminfo.py, that encode the size, coloring and mass of a given atom or residue.
    The member functions range from loading the atoms from the .pdb file and storing them by type, to creating the meshes from this information, as well as calculating bonds or the backbone.
    """
    def __init__(self, nc, fc=None):
        """
        Load structure form pdb file (by parsing file) and save Biopython structure representation of the protein.
        It also loads, above mentioned, pre-defined dictionaries from atminfo.py, that encode the size, coloring and mass of a given atom or residue.
        
        It also calculates the center of the atomic positions (as they appear in the pdb file). This will be used to re-center the atoms.

        Parameters:
            nc: NameChecker
                Instance of a NameChecker class. Used to pass the pdb file name and paths.
            fc: FileConverter, optional
                Instance of a FileConverter class. Needed if temporary files have not been created before instantiating this class. Default: None.
        """
        self._path, self._out_path, self._base_path, mesh_path = nc.return_all()
        if not fc:
            self._fc = FileConverter(nc)
        else:
            self._fc = fc
        
        self._fc.pdb_to_xyzrn(self._path, self._out_path)
        
        parser = PDBParser()
        file_name = self._path + ".pdb"
        self._structure = parser.get_structure(self._path, file_name)

        self._res_size_dict, self._res_color_dict = import_res_size_info()
        self._atoms_size_dict, self._atoms_color_dict, self._vw_dict = import_atm_size_info(True)
        
        #calculate the center and bounds of the atom cloud
        self._cam_pos = [0, 0, 0]
        self._max_coords = [0, 0, 0]
        self._centroid = [0, 0, 0]
        num_atoms = 0
        for model in self._structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() == "HOH":
                        continue
                    for atom in residue:
                        num_atoms += 1
            break

        i = 0
        s = int(num_atoms / 3)
        atoms = []
        for model in self._structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() == "HOH":
                        continue
                    for atom in residue:
                        atoms.append(atom.get_coord())
                        i += 1
            break
        center = np.array(atoms)
        self._centroid = center.mean(axis=0)
        self._max_coords = center.max(axis=0)
        self._cam_pos = [0, max(self._max_coords) * 3, 0]

    def get_structure(self):
        """
        Return the loaded structure object
        
        Returns: 
            structure
        """
        return self._structure

    def get_atoms(self, show_solvent=False, model_id=0):
        """
        Creates a dictionary that stores the 3D coordinates for each atom.
        The dictionary keys are the atom names. For each atom type in the given molecule the coordinates of the atoms of this type are stored in a list within the dictionary.
        
        The code in words:
        The .pdb file (loaded in __init__()) is iterated through. For each atom it is checked if the type of this atom is already in the dictionary. 
        If not then a new list is created with the coordinates of this atom and added to the dictionary with the name of the atom as the key.
        If the name of this atom is already present then the coordinates of the current atom are added to the list of coordinates of this same type of atom.
        The dictionary is returned.
        
        Parameters:
            show_solvent: bool, optional
                If True solvent molecules also added to retrun dictionary. Default: False.
            model_id: int, optional
                The dynamic model ID of the desired molecule. Count starts at 0. Leave default value for static molecules. Default: 0.
        
        Returns: 
            dict
                Dictionary of atomic coordinates by atom type.
            int
                Maximum coordinate in y axis. (Used to create default camera)
        """
        calculate = True
        for i in range(3):
            if self._centroid[i] != 0:
                calculate *= False

        if calculate:
            num_atoms = 0
            for model in self._structure:
                if model.id == model_id:
                    for chain in model:
                        for residue in chain:
                            if residue.get_resname() == "HOH" and not show_solvent:
                                continue
                            for atom in residue:
                                num_atoms += 1

            i = 0
            s = int(num_atoms / 3)
            atoms = []
            for model in self._structure:
                if model.id == model_id:
                    for chain in model:
                        for residue in chain:
                            if residue.get_resname() == "HOH" and not show_solvent:
                                continue
                            for atom in residue:
                                atoms.append(atom.get_coord())
                                i += 1
            center = np.array(atoms)
            self._centroid = center.mean(axis=0)
            self._max_coords = center.max(axis=0)
            self._cam_pos = [0, max(self._max_coords) * 3, 0]
        i = 0
        atom_data = dict()
        for model in self._structure:
            if model.id == model_id:
                for chain in model:
                    for residue in chain:
                        if residue.get_resname() == "HOH" and not show_solvent:
                            continue
                        for atom in residue:
                            type_name = atom.element
                            # If atom not in dictionary, add it as key with coords in list
                            if type_name not in atom_data:
                                coords = copy.deepcopy(atom.get_coord())
                                for j in range(3):
                                    coords[j] -= self._centroid[j]
                                atom_data[type_name] = [coords]
                            # If atom already in dictionary, append its coordinates to list
                            else:
                                coords = copy.deepcopy(atom.get_coord())
                                for j in range(3):
                                    coords[j] -= self._centroid[j]
                                atom_data[type_name].append(coords)
            i += 1

        # return the 3D positions of atoms organized by atom type
        return atom_data
        
            
    def get_atoms_IDs(self, model_id=0):
        """
        Get dictionary of atomic coordinates (same as get_atoms()) and residue IDs (format from output_pdb_as_xyzrn()) from the xyzrn file.
        Also return a list of all the atomic coordinates in a list in the same order as in the .pdb file.
        
        Used by the Surface class.
        
        The code in words:
        The .xyzrn file is loaded and is iterated through. The relevant fields are stored to temporary variables.
        For each atom it is checked if the type of this atom is already in the dictionary. 
        If not then a new list is created with the coordinates of this atom and added to the dictionary with the name of the atom as the key.
        If the name of this atom is already present then the coordinates of the current atom are added to the list of coordinates of this same type of atom.
        Regardless of the atom already being present in the dictionary the residue id and the coordinates are added to the two lists.
        The dictionary and the two lists are returned.
        
        Parameters:
            model_id: int, optional
                The dynamic model ID of the desired molecule. Count starts at 0. Leave default value for static molecules. Default: 0.
        
        Returns: 
            dict
                Dictionary of atomic coordinates by atom type.
            list
                List of unique residue IDs (format from output_pdb_as_xyzrn())
            list
                Atomic coordinates (in same order as the residue IDs)
        """
        xyzrnfile = open(self._out_path + ".xyzrn")
        meshdata = (xyzrnfile.read().rstrip()).split("\n")
        xyzrnfile.close()
       
        res_id = []
        atom_data = dict()
        atom_coords = []
        lenm = len(meshdata)
        for vi in range(lenm):
            fields = meshdata[vi].split()
            if int(fields[5]) == model_id:
                vertices = [None] * 3
                vertices[0] = float(fields[0]) - self._centroid[0]
                vertices[1] = float(fields[1]) - self._centroid[1]
                vertices[2] = float(fields[2]) - self._centroid[2]
                res_id.append(fields[6])
                res = fields[6].split("_")
                atmtype = res[4][0]
                if atmtype not in atom_data:
                    atom_data[atmtype] = [vertices]
                # If atom already in dictionary, append its coordinates to list
                else:
                    atom_data[atmtype].append(vertices)
                atom_coords.append(vertices)

        # return the 3D positions of atoms organized by atom type
        return atom_data, res_id, atom_coords

    def get_residues(self, model_id=0, show_solvent=False):
        """
        Creates a dictionary of coordinates by residues from structure object
        
        The code in words:
        The .pdb file (loaded in __init__()) is iterated through. For each residue it is checked if the type of this residue is already in the dictionary. 
        If not then a new list is created with the coordinates of this residue and added to the dictionary with the type of the residue as the key.
        If the type of this residue is already present then the coordinates of the current residue are added to the list of coordinates of this same type of residue.
        The coordinates of the residue are calculated as the arithmetic center of the coordinates.
        The dictionary is returned.
        
        Parameters:
            model_id: int, optional
                The dynamic model ID of the desired molecule. Count starts at 0. Leave default value for static molecules. Default: 0.
            show_solvent: bool, optional
                If True solvent molecules also added to retrun dictionary. Default: False.
        
        Returns: 
            dict
                Dictionary of atomic coordinates by residue type.
        """

        # load file into a python list
        residues = self._structure.get_residues()
        res_list = list(residues)

        # initialize dictionary, then fill up iteratively
        res_data = dict()
        for model in self._structure:
            if model.id == model_id:
                for chain in model:
                    for res in chain:
                        if res.get_resname() == "HOH" and not show_solvent:
                            continue
                        else:
                            type_name = res.get_resname()
                            center_coord = [0., 0., 0.]
                            c = 0.
                            for atom in res:
                                center_coord = center_coord + atom.get_coord()
                                center_coord -= self._centroid
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
        Calculates information about specified residue from mol2 file.
        Depending on what is specified, either the center of mass (COM) or the charge is computed.
        
        Parameters:
            res: str
                Residue number of specified residue be looked at.
            chain:                choose what property of residue you want. com for Centre Of Mass, ch for charge
                Chain number of corresponding to residue be looked at. str - options: com, ch
        
        Returns: 
            list
                List of COM coords of given (exact) residue.
        """
        # load info for given residue
        fname = self._out_path + ".mol2"#TODO: make so it works with trajectories
        if not exists(fname):
            self._fc.pdb_to_mol2(self._path, self._out_path)
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
            if a == 0:
                print("ERROR: residue empty -> center of mass could not be found")
                return 1
            COM_sum -= self._centroid * a
            COM = [i/a for i in COM_sum]
            return COM
            
        # if charge option calculate residue charge
        elif option.upper() == "CH" or option.upper() == "CHARGE":
            charge = 0.0
            for i in range(len(my)):
                charge += my.iloc[i]['charge']

            return charge

    
    def get_atom_mesh(self, atom_data, vw=0, probe=0, phi_res=10, theta_res = 10):
        """
        Create a list of Shperes and colors representing each atom for plotting. Can later be added to a mesh for plotting.
        
        The code in words:
        Iterates through the atom_data dictionary by atom type (from get_atoms()).
        It creates uniform Spheres (same size and color) in the position specified by the coordinates list for each atom type.
        Also differentiates between Van der Waals and normal radii and handles unkown atoms.
        
        Parameters:
            atom_data: dict
                Dictionary of atoms and their coordinates, by atom type.
            vw: bool, optional
                ptional - When set to True Van-der-Waals atomic radii used instead of empirical radii. Default: False.
            probe: int, optional
                size of probe (representing the solvent size) needed for surface calculation. Default: 0.
            phi_res: int, optional
                pyvista phi_resolution for Sphere objects representing atoms. Default: 10.
            theta_res: int, optional
                pyvista theta_resolution for Sphere objects representing atoms. Default: 10.
        
        Returns: 
            list
                List of pyvista Shperes representing each atom
            list
                List of colors corresponding to each atom
            list
                List of atom ID's for each atom
        """
        # these two dictionaries have to be manually created
        # also, if more/other atoms present in protein it will not work

        # create glyphs (spherical) to represent each atom
        atoms_spheres = []
        colors_list = []
        name_list = []
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
                colors_list.append(self._atoms_color_dict[atoms_type])
                name_list.append(atoms_type)
            except KeyError:
                # if atom is unrecognized, color it pink
                colors_list.append('#DD77FF')
                name_list.append("Unkown")

        # return list of spheres and colors representing each atom
        return atoms_spheres, colors_list, name_list

    def get_atom_trimesh(self, atom_data, vw=False, probe=0):
        """
        Create a list of shperes and colors representing each atom for plotting in a Trimesh format. Used for feature computation in the surface_handler class.
        
        The code in words:
        Iterates through the atom_data dictionary by atom type (from get_atoms()).
        It creates uniform Spheres (same size and color) in the position specified by the coordinates list for each atom type.
        Also differentiates between Van der Waals and normal radii and handles unkown atoms.
        
        Parameters:
            atom_data: dict
                Dictionary of atoms and their coordinates, by atom type.
            vw: bool, optional
                ptional - When set to True Van-der-Waals atomic radii used instead of empirical radii. Default: False.
            probe: int, optional
                size of probe (representing the solvent size) needed for surface calculation. Default: 0.
        
        Returns: 
            list
                List of pyvista Shperes representing each atom
            list
                List of colors for each atom
            list
                List of atom ID's for each atom
        """
        # these two dictionaries have to be manually created
        # also, if more/other atoms present in protein it will not work

        # create glyphs (spherical) to represent each atom
        atoms_spheres = []
        colors_list = []
        name_list = []
        rad = probe / 2.0
        

        for atoms_type in atom_data:

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
                colors_list.append(self._atoms_color_dict[atoms_type])
                name_list.append(atoms_type)
            except KeyError:
                # if atom is unrecognized, color it pink
                colors_list.append('#DD77FF')
                name_list.append("Unkown")

        # return list of spheres and colors representing each atom
        return atoms_spheres, colors_list, name_list

    def get_residue_mesh(self, res_data, phi_res = 25, theta_res = 25):
        """
        Create a list of Shperes and colors representing each residue for plotting. Can later be added to a mesh for plotting.
        
        The code in words:
        Iterates through the res_data dictionary by atom type (from get_residues()).
        It creates uniform Spheres (same size and color) in the position specified by the coordinates list for each residue type.
        Also differentiates between Van der Waals and normal radii and handles unkown residues.
        
        Parameters:
            res_data: dict
                Dictionary of residues and their coordinates by residue type.
            phi_res: int, optional
                pyvista phi_resolution for Sphere objects representing residues. Defaul: 25.
            theta_res: int, optional
                pyvista theta_resolution for Sphere objects representing residues. Default:25.
        
        Returns: 
            list
                List of pyvista Shperes representing each residue
            list
                List of colors for each residue
            list
                List of residue names
        """

        # create glyphs (spherical) to represent each res
        res_spheres = []
        colors_list = []
        res_names = []
        for res_type in res_data:

            sphere = pv.Sphere(radius=0.0, phi_resolution=phi_res, theta_resolution=theta_res)
            try:
                sphere = pv.Sphere(radius=self._res_size_dict[res_type], phi_resolution=phi_res, theta_resolution=theta_res)
                
                #color the sphere according to 'CPK' standard
                colors_list.append(self._res_color_dict[res_type])
                res_names.append(res_type)
                
            except KeyError:
                # color unkown to light purple
                colors_list.append('#DD77FF')
                res_names.append("Unkown")

            # create a mesh with each residues position
            mesh = pv.PolyData(np.array(res_data[res_type]))

            # place a specific sphere at given position
            glyphs = mesh.glyph(geom=sphere)
            res_spheres.append(glyphs)


        # return list of spheres and colors representing each res
        return res_spheres, colors_list, res_names

    def get_bond_mesh(self, model_id=0):
        """
        Determine bonds from 3D information.
        
        The color information is as follows:
            White for all single bonds, 
            Blue for all double bonds,
            Green for all triple bonds,
            Red for all amide bonds,
            Purple for all aromatic bonds,
            Black for everything else.
              
        The code in words:
        Parse mol2 file (also works on multi model file).
        Find where the boundaries of the current molecule are in the file.
        Extract the atomic and bond information by creating DataFrame.
        From the DataFrames get the information corresponding to the current bond: create a pyvista.Line() and store the bond type.
        Return the compiled lists of Lines and bond types (colors).
        
        Parameters:
            model_id: int, optional
                The dynamic model ID of the desired molecule. Count starts at 0. Leave default value for static molecules. Default: 0.
        
        Returns: 
            list
                List of pyvista lines representing each bond.
            list
                List of colors corresponding to the lines in the above list
            list
                List of names of the bonds: single, double, triple, amide, aromatic, unkown
        """
        calculate = True
        for i in range(3):
            if self._centroid[i] != 0:
                calculate *= False

        if calculate:
            num_atoms = 0
            for model in self._structure:
                if model.id == model_id:
                    for chain in model:
                        for residue in chain:
                            if residue.get_resname() == "HOH" :
                                continue
                            for atom in residue:
                                num_atoms += 1

            i = 0
            s = int(num_atoms / 3)
            atoms = []
            for model in self._structure:
                if model.id == model_id:
                    for chain in model:
                        for residue in chain:
                            if residue.get_resname() == "HOH":
                                continue
                            for atom in residue:
                                atoms.append(atom.get_coord())
                                i += 1
            center = np.array(atoms)
            self._centroid = center.mean(axis=0)
            self._max_coords = center.max(axis=0)
            self._cam_pos = [0, max(self._max_coords) * 3, 0]

        # Start of actual bond info calculation
        fname = self._out_path + ".mol2"
        # Check if mol2 file exists. If not convert it from pdb
        file_exists = exists(fname)
        if not file_exists:
            self._fc.pdb_to_mol2(self._path, self._out_path)
            

        k = 0
        bonds = []
        col = []
        names = []
        for mol2 in split_multimol2(fname):
            # only do computation on desired model_id
            keep_molecule = k == model_id
            if keep_molecule:
                # easy to handle the mol2 content as text -> convert it from list
                str_mol2 = ''.join(mol2[1])
                
                # find start and end of the @<TRIPOS>ATOM section
                start_atom = str_mol2.find('@<TRIPOS>ATOM')
                # 13 = len('@<TRIPOS>ATOM')
                end_atom = start_atom + 13 + str_mol2[start_atom:].replace('@<TRIPOS>ATOM','').find('@')
                # figure out number of atoms
                num_atoms = int(len(str_mol2[start_atom:end_atom].replace('@<TRIPOS>ATOM\n','').replace('\n',' ').strip().split())/9)
                # create a dataframe of atoms
                df_atoms = pd.DataFrame(np.array(str_mol2[start_atom:end_atom].replace('@<TRIPOS>ATOM\n','').replace('\n',' ').strip().split()).reshape((num_atoms,9)), columns=['atom_id', 'atom_name', 'x', 'y', 'z', 'atom_type', 'subst_id', 'subst_name', 'charge'])
                
                # find start of the @<TRIPOS>BOND section (end not needed, goes to end of string)
                start_bond = str_mol2.find('@<TRIPOS>BOND')
                # figure out number of bonds
                num_bonds = int(len(str_mol2[start_bond:].replace('@<TRIPOS>BOND\n','').replace('\n',' ').strip().split())/4)
                # create a dataframe of bonds
                df_bond = pd.DataFrame(np.array(str_mol2[start_bond:].replace('@<TRIPOS>BOND\n','').replace('\n',' ').strip().split()).reshape((num_bonds,4)), columns=['bond_id', 'atom1', 'atom2', 'bond_type'])
                
                # find each the 3D coords for each two atoms in a bond
                for i in range(num_bonds):
                    my = df_atoms.iloc[int(df_bond.iloc[i]['atom1'])-1][['x', 'y', 'z']]
                    c = df_atoms.iloc[int(df_bond.iloc[i]['atom2'])-1][['x', 'y', 'z']]
                    
                    my_coord = my.values.tolist()
                    c_coord = c.values.tolist()
                    my_coord = [float(i) for i in my_coord]
                    c_coord = [float(i) for i in c_coord]
                    my_coord -= self._centroid
                    c_coord -= self._centroid
                    if df_bond.iloc[i]['bond_type'] == '1':
                        line = pv.Line(my_coord, c_coord, resolution=5)
                        col.append('w')
                        names.append('single')
                    elif df_bond.iloc[i]['bond_type'] == '2':
                        line = pv.Line(my_coord, c_coord, resolution=1)
                        col.append('b')
                        names.append('double')
                    elif df_bond.iloc[i]['bond_type'] == '3':
                        line = pv.Line(my_coord, c_coord, resolution=1)
                        col.append('g')
                        names.append('triple')
                    # amide bond in red
                    elif df_bond.iloc[i]['bond_type'] == 'am':
                        line = pv.Line(my_coord, c_coord, resolution=1)
                        col.append('r')
                        names.append('amide')
                    # aromatic bond in purple
                    elif df_bond.iloc[i]['bond_type'] == 'ar':
                        line = pv.Line(my_coord, c_coord, resolution=1)
                        col.append('m')
                        names.append('aromatic')
                    else:
                        line = pv.Line(my_coord, c_coord, resolution=1)
                        col.append('k')
                        names.append('unkown')
                    bonds.append(line)
                break
            k += 1
        
        return bonds, col, names
    
    def get_backbone_mesh(self, model_id=0):
        """
        Creates and returns a Spline object representing the backbone of the protein.
          
        The code in words:
        Iterates through the res_data dictionary by atom type (from get_residues()).
        Calculates the center of each residue and returns these points as a numpy array.
        (Later used to create a Spline.)
             
        Parameters:         
            model_id: int, optional
                The dynamic model ID of the desired molecule. Count starts at 0. Leave default value for static molecules. Default: 0.
        
        Returns: 
            pyvista.Spline
                Spline running through coordinates representing the centre of mass of each residue.
        """
        
        i = 0
        res_list = []
        for model in self._structure:
            if model.id == model_id:
                for chain in model:
                    for res in chain:
                        if res.get_resname() != "HOH":
                            com = [0, 0, 0]
                            count = 0
                            for atom in res:
                                coords = atom.get_coord()
                                com += coords
                                com -= self._centroid
                                count += 1
                            com /= count
                            res_list.append(com)
            i += 1
        
        return pv.Spline(np.array(res_list))
