
# imports
from genericpath import sameopenfile
import os

import numpy as np
from numpy.lib.twodim_base import tri
import torch
import pyvista as pv
from pyvtk import PolyData, PointData, CellData, Scalars, Vectors, VtkData
import trimesh
from subprocess import PIPE, Popen
import open3d as o3d
import enum
from os.path import exists

from provis.utils.surface_utils import get_surface, compute_normal, prepare_trimesh, fix_trimesh, find_nearest_atom
from provis.utils.surface_feat import compute_surface_features
from provis.src.processing.file_converter import FileConverter
from provis.src.processing.data_handler import DataHandler

class SurfaceHandler:
    """
    The 'brain' of provis, when it comes to handling surface information. 
    
    This class loads information from a variety of files and creates meshes to be plotted. 
    Upper level classes - eg. StickPoint - have their own AtomHandler objects that do all the work.    
    """
  
    def __init__(self, nc, fc=None, dh=None, density=3.0):
        """
        Initializes SurfaceHandler class memeber variables: 
        _dh: DataHandler - needed for native mesh creation.
        _fc: FileConverter - needed to create necessairy files.   
        _path: str - path to pdb file (without extension).
        _out_path: str - path to output files (without extension).
        _base_path: str - path to root directory.
        _density: float - density of surface mesh.
        _features: tuple(numpy.ndarray) - A collection of arrays representing the feature information for the surface of the protein, saved in a tuple.
        _mesh: Mesh - Surface mesh of protein.
        _col: numpy.ndarray - Coloring map of surface. To be passed to pyvista.PolyData.add_mesh(scalars=).
        _res_id: list - List of residue IDs (format from output_pdb_as_xyzrn())
        _atom_coords: list - List of atomic coordinates. Just the coordinates in the same order as in _res_id.

        :param name: nc - Instance of a NameChecker class. Used to pass the pdb file name and paths.
        :param type: NameChecker
        :param name: fc - Instance of a FileConverter class. Used to create the necessairy files (face, vert, pqr, ...) if they do not exist yet. Default: None.
        :param type: FileConverter
        :param name: dh - Instance of a DataHandler class. Used to retrieve atom-positional information when calculating the surface of the protein natively. Default: None. If None a new DataHandler variable will be initialized with "nc".
        :param type: DataHandler, optional
        :param name: density - Density needed for msms. Default: 3.0.
        :param type: float, optional
        """
        
        if not fc:
            self._fc = FileConverter(nc)
        else:
            self._fc = fc
        if not dh:
            self._dh = DataHandler(nc, fc)
        else:
            self._dh = dh
        self._path, self._out_path, self._base_path = nc.return_all()
        self._density = density
        self._features = None
        self._mesh = None
        self._col = None
        self._res_id = None
        self._atom_coords = None
        
    def get_assignments(self):
        """
        Get assignments (coloring) for the mesh. File has to exist, no way to produce it with provis. 
        Loads "root directory"/data/tmp/{pdb_id}.pth and returns it.

        :returns: PyTourch object - Coloring of surface.
        """
        filename = self._out_path.upper() + '.pth'
        assigment = torch.load(filename)
        
        return assigment

    def get_surface_features(self, mesh, feature,  res_id=None):
        """
        Get the coloring corresponding to a specific feature.

        :param name: mesh - The mesh
        :param type: Trimesh
        :param name: feature - Name of feature we are interested in. Options: hydrophob, shape, charge.
        :param type: str
        :param name: native - Set to True when using native mesh computation. Default: False.
        :param type: bool, optional
        
        :raises: NotImplementedError - If unkown feature specified error is raised

        returns: numpy.ndarray - Array of coloring corresponding to surface.
        """

        # get surface
        pdb_file = self._path + '.pdb'
        
        pqr_file = self._out_path + '.pqr'
        if not exists(pqr_file):
            self._fc.pdb_to_pqr(self._path, self._out_path)
        
        if not res_id:
            surface = get_surface(self._out_path, density=self._density)
        else:
            full_res_id = find_nearest_atom(self._atom_coords, res_id, mesh.vertices)
            
            surface = [mesh.vertices, mesh.faces, mesh.face_normals, full_res_id, mesh.area_faces]
            
        # compute features of surface
        if not self._features:
            self._features = compute_surface_features(surface, pdb_file, self._out_path, mesh, pdb_id = self._path)

        if feature == 'hydrophob':
            return self._features[2]
        elif feature == 'shape':
            return self._features[0]
        elif feature == 'charge':
            return self._features[3]
        else:
            raise NotImplementedError

    def msms_mesh_and_color(self, feature=None, patch=False):
        """
        Return the mesh and coloring ready for plotting.

        :param name: feature - Name of feature, same as in get_surface_features. Options: hydrophob, shape, charge. Defaults to None.
        :param type: str, optional
        :param name: patch - Set coloring of mesh manually, from a file. If set to True get_assignments() will be called. Defaults to False.
        :param type: bool, optional

        :returns: Trimesh - The mesh, always returned
        :returns: numpy.ndarray - Color map corresponding to specified feature. Only returned if a feature is specified as an input. This map can be passed to a pyvista.Plotter.add_mesh() function as the 'scalars' argument to get the encoded coloring.
        """
        
        path = f"{self._out_path}_out_{int(self._density * 10)}"
        face = path + '.face'
        vert = path + '.vert'
        from os.path import exists
        # Check if face file exists. If not convert it from pdb
        file_exists = exists(face) and exists(vert)
        if not file_exists:
            self._fc.msms(self._out_path, self._density)
        
        if not self._mesh:
            surface = get_surface(out_path=self._out_path, density=self._density)
            self._mesh = prepare_trimesh(vertices=surface[0], faces=surface[1], normals=surface[2], 
                            resolution=1.5, apply_fixes=True)
        
        if feature:    
            if patch:
                self._col = self.get_assignments()
            if not patch:
                self._col = self.get_surface_features(self._mesh, feature)
        else:
            self._col = None

    def native_mesh_and_color(self, feature=None):
        """
        Returns a mesh without the need for the MSMS binary.
        Always returns a pyvista.PolyData mesh of the surface and if a feature is specified it also returns the coloring according to that feature.

        :param name: feature - Name of feature, same as in get_surface_features. Options: hydrophob, shape, charge. Defaults to "".
        :param type: str, optional

        :returns: pyvista.PolyData - The mesh, always returned
        :returns: numpy.ndarray - Color map corresponding to specified feature. Only returned if a feature is specified as an input. This map can be passed to a pyvista.Plotter.add_mesh() function as the 'scalars' argument to get the encoded coloring.
        """
        if not self._mesh:
            # create rough surface by combining vw radii of each atom
            atom_data, self._res_id, self._atom_coords = self._dh.get_atoms_IDs()
            self._atmsurf, col = self._dh.get_atom_mesh(atom_data, vw=1, probe=0.1)

            # adding the spheres (by atom type) one at a time
            j = 0
            mesh_ = pv.wrap(self._atmsurf[0])
            for mesh in self._atmsurf[1:]:
                mesh_ = mesh_ + (mesh)

            # blur the spheres and extract edges (pyvista)
            mesh_ = mesh_.delaunay_3d(alpha=1.5).extract_feature_edges()
            # create trimesh mesh
            new = trimesh.Trimesh(mesh_.points)
            cloud = o3d.geometry.PointCloud()
            cloud.points = o3d.utility.Vector3dVector(new.vertices)
            cloud.normals = o3d.utility.Vector3dVector(new.vertex_normals)
            # reconstruct the surface mesh as a trimesh mesh
            radii = [0.1, 0.3, 0.45, 0.6, 0.75, 0.9,1.2, 1.5]
            tri_mesh= o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(cloud, o3d.utility.DoubleVector(radii))#, depth=depth, width=width, scale=scale, linear_fit=linear_fit)
            v = np.asarray(tri_mesh.vertices)
            f = np.array(tri_mesh.triangles)
            f = np.c_[np.full(len(f), 3), f]
            mesh = pv.PolyData(v, f)
            shell =  mesh.clean().reconstruct_surface()#.clean()
            shell.compute_normals(inplace=True)

            # parse list of PolyData (format) faces and convert them into Trimesh (format) faces 
            len_f = len(shell.faces)
            faces_ = []
            i = 0
            while(i < len_f):
                curr = shell.faces[i]
                temp = [None] * curr
                idxx = 0
                j = i + 1
                while(curr):
                    temp[idxx] = shell.faces[j + idxx]
                    idxx += 1
                    curr -= 1
                faces_.append(temp)
                i += idxx + 1
            
            tri_mesh = trimesh.Trimesh(shell.points, faces=faces_)
            self._mesh = tri_mesh

        # only needed if feature is specified (-> color map needs to be calculated)
        if feature:
            self._col = self.get_surface_features(self._mesh, feature, self._res_id)
        else:
            self._col = None

    
    def return_mesh_and_color(self, msms=False, feature=None, patch=False):
        """
        Wrapper function to choose between the msms surface visualization vs the native surface visualization. 
        If you could not download the MSMS binary leave this variable false and you will not run into problems.
        
        :param name: msms - If true, surface generated by msms binary is returned, else the native mesh. Default: False.
        :param type: bool, optional
        :param name: feature - Name of feature, same as in get_surface_features. Options: hydrophob, shape, charge. Defaults to None.
        :param type: str, optional
        :param name: patch - Set coloring of mesh manually. If set to True get_assignments() will be called. Defaults to False.
        :param type: bool, optional
        
        :returns: trimesh.Trimesh - The mesh corresponding to the surface of the protein.
        :returns: numpy.ndarray - Coloring map corresponding to specified feature.
        """
        
        # run appropriate function to compute the mesh and the coloring
        if msms:
            self.msms_mesh_and_color(feature, patch)
        else:
            self.native_mesh_and_color(feature)
          
        return self._mesh, self._col
