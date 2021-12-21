
# imports
from genericpath import sameopenfile
import os

import numpy as np
import torch
# import pymesh
import pyvista as pv
from pyvtk import PolyData, PointData, CellData, Scalars, Vectors, VtkData
from provis.utils.surface_utils import get_surface, compute_normal, prepare_trimesh
from provis.src.surface_feat import compute_surface_features
import trimesh
from subprocess import PIPE, Popen
import open3d as o3d

from provis.utils.name_checker import check_name
from provis.src.data_handler import DataHandler

MSMS_BIN = os.environ['MSMS_BIN'] # 'binaries/msms

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
        self._dh = DataHandler(name)
        self._path, self._out_path = check_name(name)
        if dens:
            self._density = dens
        self._features = None
        self._mesh = None
        
    def get_assignments(self):
        # set data path
        # packagedir = holoprot.__path__[0]
        # filename = os.path.join(ROOT_DIR,
        #                         'datasets/assignments/pdbbind',
        #                         pdb.upper() + '.pth')
        # load file
        filename = self._path.upper() + '.pth'
        assigment = torch.load(filename)

        return assigment

    def get_surface_features(self, mesh, feature):

        # get surface
        pdb_file = self._path + '.pdb'
        surface = get_surface(self._out_path, density=self._density)

        # # get mesh from file
        # obj_file = pdb + '.obj'
        # mesh = trimesh.load(obj_file)
        # tmp_mesh = trimesh.Trimesh(mesh.vertices, mesh.faces)
        # normals = compute_normal(tmp_mesh.vertices, tmp_mesh.faces)
        # mesh = trimesh.Trimesh(tmp_mesh.vertices, tmp_mesh.faces, vertex_normals=normals)

        # compute features of surface
        if not self._features:
            self._features = compute_surface_features(surface, pdb_file, self._out_path, mesh, pdb_id = self._path) # TODO save features so as to not recompute all the time

        if feature == 'hydrophob':
            return self._features[2]
        elif feature == 'shape':
            return self._features[0]
        elif feature == 'charge': #TODO charge plots weirdly
            return self._features[3]
        else:
            raise NotImplementedError

    def return_mesh_and_color(self, feature="", patch=0):
        
        if not self._mesh:
            surface = get_surface(out_path=self._out_path, density=self._density)
            self._mesh = prepare_trimesh(vertices=surface[0], faces=surface[1], normals=surface[2], 
                            resolution=1.5, apply_fixes=True)
        
    
        if patch:
            cas = self.get_assignments()
        if not patch:
            cas = self.get_surface_features(self._mesh, feature)

        return self._mesh, cas

    def load_forv(self, file_name, end, vorf):
        """
        Load surface information from face or vert file
        
        :param name: file_name - Name of input file
        :param type: str
        :param name: end - Type of input file
        :param type: str
        :param name: vorf - Vertex or face file. "v" for vertex, "f" for face
        :param type: str
        
        :return: list - list of data
        """
        
        outfile = open(file_name + end, "r")
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
                if vorf == "v":
                    ret[i].append(float(entry))
                elif vorf == "f":
                    ret[i].append(int(entry)-1)
            i+=1

        outfile.close()
        return ret


    def poisson_mesh(self, pc, depth=8, width=0, scale=1.1, linear_fit=False):
        """`pc` is a `pyvista.PolyData` point cloud. The default arguments are abitrary"""
        cloud = o3d.geometry.PointCloud()
        poly_pc = pc.extract_all_edges()
        cloud.points = o3d.utility.Vector3dVector(poly_pc.points)
        cloud.normals = o3d.utility.Vector3dVector(poly_pc["norms"])
        trimesh, _ = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(cloud, depth=depth, width=width, scale=scale, linear_fit=linear_fit)
        v = np.asarray(trimesh.vertices)
        f = np.array(trimesh.triangles)
        f = np.c_[np.full(len(f), 3), f]
        mesh = pv.PolyData(v, f)
        return mesh.clean()

    def native_mesh(self):

        atom_data = self._dh.get_atoms()
        # change vw and probe to get a finer/better/differently refined mesh
        self._atmsurf, col = self._dh.get_atom_mesh(atom_data, vw=1, probe=0.1)

        # adding the spheres (by atom type) one at a time
        j = 0
        mesh_ = pv.wrap(self._atmsurf[0])
        for mesh in self._atmsurf[1:]:
            mesh_ = mesh_ + (mesh)


        # shell = mesh_
        # # shell = trimesh.smoothing.filter_taubin(mesh_)
        # pcd = mesh_.sample_points_poisson_disk(750)
        # o3d.visualization.draw_geometries([pcd])
        # alpha = 0.03
        # print(f"alpha={alpha:.3f}")
        # mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd, alpha)
        # mesh.compute_vertex_normals()
        # o3d.visualization.draw_geometries([mesh], mesh_show_back_face=True)
       
        # FLIPPING WORKDSSSSSSSSSS
        mesh_ = mesh_.delaunay_3d(alpha=1.5).extract_feature_edges()
        # self._atmsurf, col = self._dh.get_atom_trimesh(atom_data, vw=1, probe=0.1)
        new = trimesh.Trimesh(mesh_.points)#trimesh.util.concatenate(mesh_)
        cloud = o3d.geometry.PointCloud()
        poly_pc = new.vertices[new.edges_unique]
        # poly_pc = pc.extract_all_edges()
        cloud.points = o3d.utility.Vector3dVector(new.vertices)
        cloud.normals = o3d.utility.Vector3dVector(new.vertex_normals)
        radii = [0.1, 0.3, 0.45, 0.6, 0.75, 0.9,1.2, 1.5]
        tri_mesh= o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(cloud, o3d.utility.DoubleVector(radii))#, depth=depth, width=width, scale=scale, linear_fit=linear_fit)
        v = np.asarray(tri_mesh.vertices)
        f = np.array(tri_mesh.triangles)
        f = np.c_[np.full(len(f), 3), f]
        mesh = pv.PolyData(v, f)
        shell =  mesh.clean().reconstruct_surface()#.clean()

        # shell = self.poisson_mesh(mesh_) 

        
        # my_tubes = dict()
        # for i in range(len(self._atmsurf)):
        #     my_tubes[i] = self._atmsurf[i]

        # blocks = pv.MultiBlock(my_tubes)
        # merged = blocks.combine()
        # merged # this is now a single unstructured grid containing all geometry
        # shell = merged#mesh_.extract_surface().extract_all_edges()

        # prepare_trimesh()
            
        # create one mesh out of many spheres
        # vol = mesh_.delaunay_3d(alpha=1.5)
        # # extract surface from new mesh
        # # # shell = self.poisson_mesh(vol)
        # shell = vol.extract_feature_edges(12).reconstruct_surface().clean()
        # shell = shell.reconstruct_surface(sample_spacing=1.3)
        # shell = vol.extract_surface().reconstruct_surface(sample_spacing=1.2)
        # shell = shell.extract_all_edges().delaunay_3d(alpha=1.4)

        return shell
