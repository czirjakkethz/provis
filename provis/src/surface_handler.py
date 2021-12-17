
# imports
import os

import numpy as np
import torch
# import pymesh
import pyvista as pv
from pyvtk import PolyData, PointData, CellData, Scalars, Vectors, VtkData
from provis.src.surface_utils import get_surface, compute_normal, prepare_trimesh
from provis.src.surface_feat import compute_surface_features
import trimesh
from subprocess import PIPE, Popen

ROOT_DIR = os.environ['PROT']
MSMS_BIN = os.environ['MSMS_BIN']

APBS_BIN = os.environ['APBS_BIN'] # "binaries/apbs/bin/apbs"
PDB2PQR_BIN = os.environ['PDB2PQR_BIN'] # "binaries/pdb2pqr/pdb2pqr"
MULTIVALUE_BIN = os.environ['MULTIVALUE_BIN'] # "binaries/apbs/share/apbs/tools/bin/multivalue"

class SurfaceHandler:
    """
    The 'brain' of provis, when it comes to handling surface information. 
    
    This class loads information from a variety of files and creates meshes to be plotted. 
    Upper level classes - eg. StickPoint - have their own AtomHandler objects that do all the work.    
    """
  
    def __init__(self, name, dens=None):
        self._name = name
        if dens:
            self._density = dens

    def save_vtk(self, fname, xyz, triangles=None, values=None, vectors=None,
                triangle_values=None):
        """
        Saves a point cloud or triangle mesh as a .vtk file.

        Files can be opened with Paraview or displayed using the PyVista library.
        Args:
            fname (string): filename.
            xyz (Tensor): (N,3) point cloud or vertices.
            triangles (integer Tensor, optional): (T,3) mesh connectivity.
            values (Tensor, optional): (N,D) values, supported by the vertices.
            vectors (Tensor, optional): (N,3) vectors, supported by the vertices.
            triangle_values (Tensor, optional): (T,D) values, supported by triangles.
        """

        # encode the points/vertices as a VTK structure:
        if triangles is None:
            # point cloud
            structure = PolyData(points=xyz)
        else:
            # surface mesh
            structure = PolyData(points=xyz, polygons=triangles)

        data = [structure]
        pointdata, celldata = [], []

        # point values - one channel per column of the `values` array
        if values is not None:
            values = values
            if len(values.shape) == 1:
                values = values[:, None]
            features = values.T
            pointdata += [
                Scalars(f,
                        name=f"features_{i:02d}") for i, f in enumerate(features)
            ]

        # point vectors - one vector per point
        if vectors is not None:
            pointdata += [Vectors(vectors, name="vectors")]

        # store in the VTK object:
        if pointdata != []:
            pointdata = PointData(*pointdata)
            data.append(pointdata)

        # triangle values - one channel per column of the `triangle_values` array
        if triangle_values is not None:
            triangle_values = triangle_values
            if len(triangle_values.shape) == 1:
                triangle_values = triangle_values[:, None]
            features = triangle_values.T
            celldata += [
                Scalars(f,
                        name=f"features_{i:02d}") for i, f in enumerate(features)
            ]

            celldata = CellData(*celldata)
            data.append(celldata)

        #  write to hard drive
        vtk = VtkData(*data)
        vtk.tofile(fname)


    def load_mesh_from_file(self, mesh_file: str):

        mesh = trimesh.load(mesh_file)
        tmp_mesh = trimesh.Trimesh(mesh.vertices, mesh.faces)
        
        normals = compute_normal(tmp_mesh.vertices, tmp_mesh.faces)

        mesh = trimesh.Trimesh(tmp_mesh.vertices, tmp_mesh.faces, vertex_normals=normals)
        return mesh

    def get_assignments(self, pdb):
        # set data path
        # packagedir = holoprot.__path__[0]
        # filename = os.path.join(ROOT_DIR,
        #                         'datasets/assignments/pdbbind',
        #                         pdb.upper() + '.pth')
        # load file
        filename = pdb.upper() + '.pth'
        assigment = torch.load(filename)

        return assigment


    def get_mesh(self, pdb):
        # set data path
        # packagedir = holoprot.__path__[0]
        # filename = os.path.join(ROOT_DIR,
        #                         'datasets/surface_mesh/processed/pdbbind',
        #                         pdb + '.obj')

        filename = pdb + '.obj'
        # load file
        mesh = self.load_mesh_from_file(filename)
        return mesh


    def get_pdb_file(self, pdb):
        # set data path
        # packagedir = holoprot.__path__[0]
        # filename = os.path.join(ROOT_DIR,
        #                         f'datasets/raw/pdbbind/pdb_files/{pdb}',
        #                         pdb + '.pdb')
        filename = pdb + '.pdb'
        return filename


    def get_surface_features(self, pdb, feature):

        # get surface
        pdb_file = self.get_pdb_file(pdb)
        surface = get_surface(pdb_file, msms_bin=MSMS_BIN)

        # get mesh
        mesh = self.get_mesh(pdb)
        
        # compute features of surface
        features = compute_surface_features(surface, pdb_file, mesh, pdb_id = pdb)

        if feature == 'hydrophob':
            return features[2]
        elif feature == 'shape':
            return features[0]
        elif feature == 'charge':
            return features[3]
        else:
            raise NotImplementedError

    def mesh_color(self, feature="", patch=0):

        surface = get_surface(pdb_file=self._name, msms_bin=MSMS_BIN, density=self._density)
        mesh = prepare_trimesh(vertices=surface[0], faces=surface[1], normals=surface[2], 
                        resolution=1.5, apply_fixes=True)
        if patch:
            cas = self.get_assignments(self._name)
        if not patch:
            cas = self.get_surface_features(self._name, feature)

        self.save_vtk('tmp.vtk', mesh.vertices, mesh.faces.astype(int))
        mesh = pv.PolyData('tmp.vtk')
        return mesh, cas

