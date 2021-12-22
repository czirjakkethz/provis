
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

from provis.utils.surface_utils import get_surface, compute_normal, prepare_trimesh, fix_trimesh
from provis.utils.name_checker import check_name
from provis.src.processing.data_handler import DataHandler
from provis.src.processing.surface_feat import compute_surface_features

class SurfaceHandler:
    """
    The 'brain' of provis, when it comes to handling surface information. 
    
    This class loads information from a variety of files and creates meshes to be plotted. 
    Upper level classes - eg. StickPoint - have their own AtomHandler objects that do all the work.    
    """
  
    def __init__(self, name, dens=None):
        """
        Initializes SurfaceHandler class memeber variables: 
        _dh: DataHandler - needed for native mesh creation
        _path: str - path to pdb file (without extension)
        _out_path: str - path to output files (without extension)
        _density: float - density of surface mesh
        _features: #TODO
        _mesh: Mesh - Surface mesh of protein

        :param name: name - name to be passed to name_checker.check_name()
        :param type: str
        :param name: dens - Density needed for msms. Defaults to None.
        :param type: float, optional
        """
        self._dh = DataHandler(name)
        self._path, self._out_path = check_name(name)
        if dens:
            self._density = dens
        self._features = None
        self._mesh = None
        
    def get_assignments(self):
        """
        Get assignments (coloring) for the mesh. File has to exist, no way to produce it with provis. 
        Loads provis/data/tmp/{pdb_id}.pth and returns it.

        :returns: [type] - [description] #TODO
        """
        filename = self._out_path.upper() + '.pth'
        assigment = torch.load(filename)
        
        print(assigment)
        print(type(assigment))
        return assigment

    def get_surface_features(self, mesh, feature):
        """
        Get the coloring corresponding to a specific feature.

        Args:
            mesh (Trimesh): The mesh
            feature (str): Name of feature we are interested in. Options: hydrophob, shape, charge

        Raises:
            NotImplementedError: If unkown feature specified error is raised

        Returns:
            [type]: [description] # TODO
        """

        # get surface
        pdb_file = self._path + '.pdb'
        surface = get_surface(self._out_path, density=self._density)

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

    def return_mesh_and_color(self, feature="", patch=False, simple=False):
        """
        Return the mesh and coloring ready for plotting.

        :param name: feature - Name of feature, same as in get_surface_features. Options: hydrophob, shape, charge. Defaults to "".
        :param type: str, optional
        :param name: patch - Set coloring of mesh manually. If set to True get_assignments() will be called. Defaults to False.
        :param type: bool, optional
        :param name: simple - If set to True only mesh is returned without coloring. Used for surface plotting. Defaults to False.
        :param type: bool, optional

        :returns: Trimesh - The mesh, always returned
            #TODO: Coloring map for mesh corresponding to the specified feature 
        """
        
        if not self._mesh:
            surface = get_surface(out_path=self._out_path, density=self._density)
            self._mesh = prepare_trimesh(vertices=surface[0], faces=surface[1], normals=surface[2], 
                            resolution=1.5, apply_fixes=True)
        
        if not simple:    
            if patch:
                cas = self.get_assignments()
            if not patch:
                cas = self.get_surface_features(self._mesh, feature)
            
            return self._mesh, cas

        else:
            return self._mesh

    def native_mesh(self):
        """
        Returns a mesh without the need for the MSMS binary.

        :returns: pyvista.PolyData - The mesh
        """

        atom_data = self._dh.get_atoms()
        # change vw and probe to get a finer/better/differently refined mesh
        self._atmsurf, col = self._dh.get_atom_mesh(atom_data, vw=1, probe=0.1)

        # adding the spheres (by atom type) one at a time
        j = 0
        mesh_ = pv.wrap(self._atmsurf[0])
        for mesh in self._atmsurf[1:]:
            mesh_ = mesh_ + (mesh)


        mesh_ = mesh_.delaunay_3d(alpha=1.5).extract_feature_edges()
        # self._atmsurf, col = self._dh.get_atom_trimesh(atom_data, vw=1, probe=0.1)
        new = trimesh.Trimesh(mesh_.points)#trimesh.util.concatenate(mesh_)
        cloud = o3d.geometry.PointCloud()
        cloud.points = o3d.utility.Vector3dVector(new.vertices)
        cloud.normals = o3d.utility.Vector3dVector(new.vertex_normals)
        radii = [0.1, 0.3, 0.45, 0.6, 0.75, 0.9,1.2, 1.5]
        tri_mesh= o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(cloud, o3d.utility.DoubleVector(radii))#, depth=depth, width=width, scale=scale, linear_fit=linear_fit)
        v = np.asarray(tri_mesh.vertices)
        f = np.array(tri_mesh.triangles)
        f = np.c_[np.full(len(f), 3), f]
        mesh = pv.PolyData(v, f)
        shell =  mesh.clean().reconstruct_surface()#.clean()
        
        shell.compute_normals(inplace=True)
        
        # Work in progress; trying to figure out how to plot features with native mesh
        # tri_mesh = trimesh.Trimesh(shell.points, shell.faces)
        # fix_trimesh(tri_mesh)
        
        # self._mesh = prepare_trimesh(tri_mesh.vertices, tri_mesh.faces)

        return shell
