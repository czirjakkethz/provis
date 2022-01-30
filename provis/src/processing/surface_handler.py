
# imports
from genericpath import sameopenfile
import os
from os.path import exists
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
import pickle

from provis.utils.surface_utils import get_surface, compute_normal, prepare_trimesh, fix_trimesh, find_nearest_atom
from provis.utils.surface_feat import compute_surface_features
from provis.src.processing.file_converter import FileConverter
from provis.src.processing.data_handler import DataHandler

class SurfaceHandler:
    """
    The 'brain' of provis, when it comes to handling surfaces. 
    
    This class loads information from a variety of files and creates meshes to be plotted. 
    Upper level classes - eg. Surface - use SurfaceHandler objects to create the meshes.
    
    The SurfaceHandler class loads atom-positional information from a .pdb file and from this information computes the surface mesh.
    The surface can be computed by the MSMS binary or natively. The MSMS version is chemically more accurate and faster, but the MSMS binary has to be downloaded for it to work.
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

        Parameters:
            nc: NameChecker
                Instance of a NameChecker class. Used to pass the pdb file name and paths.
            fc: FileConverter
                Instance of a FileConverter class. Used to create the necessairy files (face, vert, pqr, ...) if they do not exist yet. Default: None.
            dh: DataHandler, optional
                Instance of a DataHandler class. Used to retrieve atom-positional information when calculating the surface of the protein natively. Default: None. If None a new DataHandler variable will be initialized with "nc".
            density: float, optional
                Density needed for msms. Default: 3.0.
        """
        
        if not fc:
            self._fc = FileConverter(nc)
        else:
            self._fc = fc
        if not dh:
            self._dh = DataHandler(nc, fc)
        else:
            self._dh = dh
        self._path, self._out_path, self._base_path, self._mesh_path = nc.return_all()
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

        Returns: 
            PyTourch object
                Coloring of surface.
        """
        filename = self._out_path.upper() + '.pth'
        assigment = torch.load(filename)
        
        return assigment

    def get_surface_features(self, mesh, feature,  res_id=None):
        """
        Get the coloring corresponding to a specific feature.

        The code in words:
        Creates the pqr file if it does not exist.
        If the res_id variable is not empty retrieve the surface structure (list of specific surface-related information) using the get_surface() function.
        Else compile the surface structure from the mesh and the find_nearest_atom() function.
        Next, using the above mentioned surface structure compute all surface feature information.
        Finally return the coloring corresponding to the specified feature.

        Parameters:
            mesh: Trimesh
                The mesh
            feature: str
                Name of feature we are interested in. Options: hydrophob, shape, charge, hbonds.
            res_id: bool, optional
                List of unique residue IDs (format from output_pdb_as_xyzrn()). Default: False.
        
        Returns: 
            numpy.ndarray
                Array of coloring corresponding to surface.

        Raises:
            NotImplementedError
                If unkown feature specified error is raised
        """

        # get surface
        pdb_file_path = self._path 
        
        pqr_file_path = self._out_path
        if self._color_needed and self._dynamic:
            pdb_file_path = self._path + "_" + str(self._model_id) 
            pqr_file_path = self._out_path + "_" + str(self._model_id)
        if not exists(pqr_file_path + '.pqr'):
            self._fc.pdb_to_pqr(pdb_file_path, pqr_file_path)
        pdb_file = pdb_file_path + '.pdb'
        if not res_id:
            surface = get_surface(pqr_file_path, density=self._density)
        else:
            full_res_id = find_nearest_atom(self._atom_coords, res_id, mesh.vertices)
            
            surface = [mesh.vertices, mesh.faces, mesh.face_normals, full_res_id, mesh.area_faces]
            
        # compute features of surface
        if not self._features or self._color_needed:
            self._features = compute_surface_features(surface, pdb_file, pqr_file_path, mesh, pdb_id = pdb_file_path)

        if feature == 'hydrophob':
            return self._features[2]
        elif feature == 'shape':
            return self._features[0]
        elif feature == 'charge':
            return self._features[3]
        elif feature == 'hbonds':
            return self._features[1]
        else:
            raise NotImplementedError

    def msms_mesh_and_color(self, feature=None, patch=False):
        """
        Return the mesh and coloring created by the MSMS binary.

        The code in words:
        If self._mesh_needed is set to True - if the mesh could not be loaded from a file - compute the mesh using the .face and .vert files created by the MSMS binary.
        If self._dynamic is False - if the mesh is needed for a static molecule - and the .face and .vert files do not exist compute these files using the FileConverter class.
        Do the same for self._dynamic == True.
        If the mesh is needed compute the mesh from the .face and .vert files.
        Finally, if a feature is specified check if the color information could be loaded from a file (self._color_needed) and compute it if needed.
        If no feature specified set the variable that stores color information to None. This will result in a white mesh.

        Parameters:
            feature: str, optional
                Name of feature, same as in get_surface_features. Options: hydrophob, shape, charge, hbonds. Defaults to None.
            patch: bool, optional
                Set coloring of mesh manually, from a file. If set to True get_assignments() will be called. Defaults to False.
        """
        print("MSMS mesh calculation for model id: ", self._model_id)

        path = f"{self._out_path}_out_{int(self._density * 10)}"
        face = path + '.face'
        vert = path + '.vert'
        # Check if face file exists. If not convert it from pdb
        file_exists = exists(face) and exists(vert)
        if self._mesh_needed:
            if not self._dynamic:
                if not file_exists:
                    self._fc.msms(self._out_path, self._density)
            
            new_xyzrn_path = self._out_path
            if self._dynamic:
                new_xyzrn_path = self._out_path + "_" + str(self._model_id) 
                
                # convert current mesh
                path = f"{new_xyzrn_path}_out_{int(self._density * 10)}"
                face = path + '.face'
                vert = path + '.vert'
                # Check if face file exists. If not convert it from pdb
                file_exists = exists(face) and exists(vert)
                if not file_exists:
                    print("Creating face and vert files for current model...")
                    new_name = self._path + "_" + str(self._model_id) 
                    
                    self._fc.pdb_to_xyzrn(new_name, new_xyzrn_path)
                    self._fc.msms(new_xyzrn_path, self._density)
                    path = f"{new_xyzrn_path}_out_{int(self._density * 10)}"
                
            print(" - Get surface")
            surface = get_surface(out_path=new_xyzrn_path, density=self._density)
            self._mesh = prepare_trimesh(vertices=surface[0], faces=surface[1], normals=surface[2], 
                            resolution=1.5, apply_fixes=True)
        
            meshname = self._mesh_path + "_msms_" + str(self._model_id) + '.obj'
            self._mesh.export(meshname)
        
        if feature:
            if self._color_needed:
                print(" - Feature: ", feature)
                if patch:
                    self._col = self.get_assignments()
                if not patch:
                    self._col = self.get_surface_features(self._mesh, feature)
                fname = self._mesh_path + "_msms_" + feature + "_" + str(self._model_id)
                with open(fname, 'wb') as f:
                    pickle.dump(self._col, f)
        else:
            self._col = None

    def native_mesh_and_color(self, feature=None):
        """
        Returns a mesh without the need for the MSMS binary. The mesh and coloring is also saved to the following file names:
        Mesh: "root directory"/data/meshes/{pdb_id}_{model_id}.obj
        Color: "root directory"/data/meshes/{pdb_id}_{feature}_{model_id}
        Always returns a pyvista.PolyData mesh of the surface and if a feature is specified it also returns the coloring according to that feature.

        Parameters:
            feature: str, optional
                Name of feature, same as in get_surface_features. Options: hydrophob, shape, charge, hbonds. Defaults to "".
        """
        print("Native mesh calculation for model id: ", self._model_id)
        if self._mesh_needed:

            print(" - Feature: ", feature)
            # create rough surface by combining vw radii of each atom
            atom_data, self._res_id, self._atom_coords = self._dh.get_atoms_IDs(model_id=self._model_id)
            self._atmsurf, col, _ = self._dh.get_atom_mesh(atom_data, vw=1, probe=0.1)

            # adding the spheres (by atom type) one at a time
            j = 0
            mesh_ = pv.wrap(self._atmsurf[0])
            for mesh in self._atmsurf[1:]:
                mesh_ = mesh_ + (mesh)
            
            # blur the spheres and extract edges (pyvista)
            mesh_ = mesh_.extract_surface()#.extract_feature_edges(1)#.delaunay_3d(alpha=a).extract_feature_edges(2)
            # create trimesh mesh
            new = trimesh.Trimesh(mesh_.points)
            cloud = o3d.geometry.PointCloud()
            cloud.points = o3d.utility.Vector3dVector(new.vertices)
            # cloud.normals = o3d.utility.Vector3dVector(new.vertex_normals)
            cloud.estimate_normals()

            # to obtain a consistent normal orientation
            cloud.orient_normals_towards_camera_location(cloud.get_center())

            # or you might want to flip the normals to make them point outward, not mandatory
            cloud.normals = o3d.utility.Vector3dVector( - np.asarray(cloud.normals))
            # reconstruct the surface mesh as a trimesh mesh
            tri_mesh, _ = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(cloud, depth=5)
            print(" - Poisson completed")
            v = np.asarray(tri_mesh.vertices)
            f = np.array(tri_mesh.triangles)

            tri_mesh = trimesh.Trimesh(v, faces=f)

            self._mesh = tri_mesh
            meshname = self._mesh_path + "_" + str(self._model_id) + '.obj'
            tri_mesh.export(meshname)
            print("Mesh calculation done.")
            

        # only needed if feature is specified (-> color map needs to be calculated)
        if feature:
            if self._color_needed:
                self._col = self.get_surface_features(self._mesh, feature, self._res_id)
                
                fname = self._mesh_path + "_" + feature + "_" + str(self._model_id)
                with open(fname, 'wb') as f:
                    pickle.dump(self._col, f)
        else:
            self._col = None

    
    def return_mesh_and_color(self, msms=False, feature=None, patch=False, model_id=0, num_models=0):
        """
        Wrapper function to choose between the msms surface visualization vs the native surface visualization. 
        If you could not download the MSMS binary leave the msms variable as False.
        
        If the mesh and appropriate coloring has already been stored to a file, then this information will be loaded and no computation will be done.
        
        The code in words:
        First, set a few class member variables that are used in the {*}_mesh_and_color() member functions.
        Next check if the mesh and color information has already been computed. If the files already exist load them and return. No computation will be done.
        Otherwise, depending on what the msms input variable is set to, compute either the msms_mesh_and_color() ir the native_mesh_and_color().
        Return the mesh and color information.

        Parameters:
            msms: bool, optional
                If true, surface generated by msms binary is returned, else the native mesh. Default: False.
            feature: str, optional
                Name of feature, same as in get_surface_features. Options: hydrophob, shape, charge, hbonds. Defaults to None.
            patch: bool, optional
                Set coloring of mesh manually. If set to True get_assignments() will be called. Defaults to False.
            model_id: int, optional
                The dynamic model ID of the desired molecule. Count starts at 0. Leave default value for static molecules. Default: 0.
            num_models: int, optional
                Number of models in a trajectory. Only important for dynamic structures. Default: 0.
        
        Returns: 
            trimesh.Trimesh
                The mesh corresponding to the surface of the protein.
            numpy.ndarray
                Coloring map corresponding to specified feature.
        """
        self._color_needed = True
        self._mesh_needed = True
        self._msms = msms
        self._model_id = model_id
        self._dynamic = num_models > 0
        connect = "_"
        if msms:
            connect = "_msms_"
            
        meshname = self._mesh_path + connect + str(model_id) + '.obj'
        mesh_exists = exists(meshname)
        if feature:
            fname = self._mesh_path + connect + feature + "_" + str(model_id)
            col_exists = exists(fname)
            if mesh_exists and col_exists:
                print("Loading mesh from: ", meshname)
                print("Loading features from: ", fname)
                self._mesh = trimesh.load_mesh(meshname)
                self._mesh_needed = False
                with open(fname, 'rb') as f:
                    self._col = pickle.load(f)
                self._color_needed = False
        else:
            if mesh_exists:
                self._mesh = trimesh.load_mesh(meshname)
                print("Loading mesh from: ", meshname)
                self._mesh_needed = False
           
        if self._color_needed or self._mesh_needed:
            # run appropriate function to compute the mesh and the coloring
            if msms:
                self.msms_mesh_and_color(feature, patch)
            else:
                self.native_mesh_and_color(feature)

        return self._mesh, self._col
