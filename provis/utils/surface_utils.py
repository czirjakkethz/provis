"""
https://github.com/bunnech/holoprot/blob/main/holoprot/utils/surface.py is the base for this file. Modifications were made.

Utilities for preparing and computing features on molecular surfaces.
"""
import os
import numpy as np
from numpy.core.numeric import full
from numpy.matlib import repmat
from Bio.PDB import PDBParser, Selection
from subprocess import Popen, PIPE
from typing import Tuple, List, Dict
import trimesh 

from provis.utils import RADII, POLAR_HYDROGENS

eps = 1e-6
Surface = Tuple[np.ndarray, np.ndarray, np.ndarray, List[str], Dict[str, str]]


def read_msms(
        file_root: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str]]:
    """Read surface constituents from output files generated using MSMS.

    :param name: file_root - Root name for loading .face and .vert files (produced by MSMS). Default location is data/tmp/{pdb_id}s
    :param type: str
   
    :returns: numpy.ndarray - vertices
    :returns: numpy.ndarray - faces
    :returns: numpy.ndarray - vertex normals
    :returns: list - List of res_id's from output_pdb_as_xyzrn()
    """
    vertfile = open(file_root + ".vert")
    meshdata = (vertfile.read().rstrip()).split("\n")
    vertfile.close()

    # Read number of vertices.
    count = {}
    header = meshdata[2].split()
    count["vertices"] = int(header[0])
    ## Data Structures
    vertices = np.zeros((count["vertices"], 3))
    normalv = np.zeros((count["vertices"], 3))
    atom_id = [""] * count["vertices"]
    res_id = [""] * count["vertices"]
    for i in range(3, len(meshdata)):
        fields = meshdata[i].split()
        vi = i - 3
        vertices[vi][0] = float(fields[0])
        vertices[vi][1] = float(fields[1])
        vertices[vi][2] = float(fields[2])
        normalv[vi][0] = float(fields[3])
        normalv[vi][1] = float(fields[4])
        normalv[vi][2] = float(fields[5])
        atom_id[vi] = fields[7]
        res_id[vi] = fields[9]
        count["vertices"] -= 1

    # Read faces.
    facefile = open(file_root + ".face")
    meshdata = (facefile.read().rstrip()).split("\n")
    facefile.close()

    # Read number of faces.
    header = meshdata[2].split()
    count["faces"] = int(header[0])
    faces = np.zeros((count["faces"], 3), dtype=int)
    normalf = np.zeros((count["faces"], 3))

    for i in range(3, len(meshdata)):
        fi = i - 3
        fields = meshdata[i].split()
        faces[fi][0] = int(fields[0]) - 1
        faces[fi][1] = int(fields[1]) - 1
        faces[fi][2] = int(fields[2]) - 1
        count["faces"] -= 1

    assert count["vertices"] == 0
    assert count["faces"] == 0

    return vertices, faces, normalv, res_id


def output_pdb_as_xyzrn(pdb_file: str, xyzrn_file: str) -> None:
    """
    Converts a .pdb file to a .xyzrn file.

    :param name: pdb_file - path to PDB File to convert (with extension)
    :param type: str
    :param name: xyzrn_file - path to the xyzrn File (with extension)
    :param type: str
    """
    # get absolute path for xyzrn_file, needed for open
    xyzrn_file = os.path.join(os.getcwd(), xyzrn_file)

    filename = pdb_file.split("/")[-1]
    pdb_id = filename.split(".")[0]

    parser = PDBParser()
    struct = parser.get_structure(id=pdb_id, file=pdb_file)
    outfile = open(xyzrn_file, "w")

    for atom in struct.get_atoms():
        name = atom.get_name()
        residue = atom.get_parent()
        # Ignore hetatms.
        if residue.get_id()[0] != " ":
            continue
        resname = residue.get_resname()
        reskey = residue.get_id()[1]
        chain = residue.get_parent().get_id()
        atomtype = name[0]

        color = "Green"
        coords = None
        if atomtype in RADII and resname in POLAR_HYDROGENS:
            if atomtype == "O":
                color = "Red"
            if atomtype == "N":
                color = "Blue"
            if atomtype == "H":
                if name in POLAR_HYDROGENS[resname]:
                    color = "Blue"  # Polar hydrogens
            coords = "{:.06f} {:.06f} {:.06f}".format(atom.get_coord()[0],
                                                      atom.get_coord()[1],
                                                      atom.get_coord()[2])
            insertion = "x"
            if residue.get_id()[2] != " ":
                insertion = residue.get_id()[2]
            # this is the weird ID; the 9th column of the xyzrn file
            full_id = "{}_{:d}_{}_{}_{}_{}".format(
                chain,
                residue.get_id()[1], insertion, resname, name, color)
        if coords is not None:
            outfile.write(coords + " " + RADII[atomtype] + " 1 " + full_id +
                          "\n")


def get_surface(out_path: str, density: float):
    """
    Wrapper function that reads in the output from the MSMS executable to build the protein surface.

    :param name: out_path - path to output (output path from namechecker) directory. Usually data/tmp
    :param type: str
    :param name: density - Need to pass same density as used by the MSMS binary, as the face and vert files have the density included in their names. The variable is needed for loading these files.
    :param type: bool
        
    :returns: numpy.ndarray - vertices
    :returns: numpy.ndarray - faces
    :returns: numpy.ndarray - vertex normals
    :returns: list - List of res_id's from output_pdb_as_xyzrn()
    :returns: dict - Dictionary: residues as keys, areas as values
    """
    file_base = os.path.abspath(out_path)
    file_base = f"{file_base}_out_{int(density * 10)}"

    # file_base = pdb_path.split(".")[0]
    # file_base = f"{pdb_path}_out_{int(density * 10)}"

    vertices, faces, normals, names = read_msms(file_base)
    areas = {}
    ses_file = open(file_base + ".area")
    next(ses_file)  # ignore header line
    for line in ses_file:
        fields = line.split()
        areas[fields[3]] = fields[1]

    return vertices, faces, normals, names, areas


def compute_normal(vertices: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """
    Compute normals for the vertices and faces

    :param name: vertices - Vertices of the mesh
    :param type: np.ndarray 
    :param name: faces - Faces of the mesh
    :param type: np.ndarray 
    
    :returns: np.ndarray - Normals of the mesh
    """
    vertex = vertices.T
    face = faces.T
    nface = np.size(face, 1)
    nvert = np.size(vertex, 1)
    normal = np.zeros((3, nvert))
    # unit normals to the faces
    normalf = crossp(
        vertex[:, face[1, :]] - vertex[:, face[0, :]],
        vertex[:, face[2, :]] - vertex[:, face[0, :]],
    )
    sum_squares = np.sum(normalf**2, 0)
    d = np.sqrt(sum_squares)
    d[d < eps] = 1
    normalf = normalf / repmat(d, 3, 1)
    # unit normal to the vertex
    normal = np.zeros((3, nvert))
    for i in np.arange(0, nface):
        f = face[:, i]
        for j in np.arange(3):
            normal[:, f[j]] = normal[:, f[j]] + normalf[:, i]

    # normalize
    d = np.sqrt(np.sum(normal**2, 0))
    d[d < eps] = 1
    normal = normal / repmat(d, 3, 1)
    # enforce that the normal are outward
    vertex_means = np.mean(vertex, 0)
    v = vertex - repmat(vertex_means, 3, 1)
    s = np.sum(np.multiply(v, normal), 1)
    if np.sum(s > 0) < np.sum(s < 0):
        # flip
        normal = -normal
        normalf = -normalf
    return normal.T


def crossp(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """
    Creates the cross product of two numpy arrays

    :param name: x - Array 1
    :param type: np.ndarray
    :param name: y - Array 2
    :param type: np.ndarray

    :returns: np.ndarray - (Array 1) x (Array 2)
    """
    z = np.zeros((x.shape))
    z[0, :] = np.multiply(x[1, :], y[2, :]) - np.multiply(x[2, :], y[1, :])
    z[1, :] = np.multiply(x[2, :], y[0, :]) - np.multiply(x[0, :], y[2, :])
    z[2, :] = np.multiply(x[0, :], y[1, :]) - np.multiply(x[1, :], y[0, :])
    return z

def prepare_trimesh(vertices: np.ndarray,
                 faces: np.ndarray,
                 normals: np.ndarray = None,
                 resolution: float = 1.0,
                 apply_fixes: bool = False):
    """
    Prepare the mesh surface given vertices and faces. Optionally, compute
    normals and apply fixes to mesh.

    :param name: vertices - Surface vertices
    :param type: np.ndarray
    :param name: faces - Triangular faces on the mesh
    :param type: np.ndarray
    :param name: normals - Normals for each vertex
    :param type: np.ndarray
    :param name: apply_fixes - Optional application of fixes to mesh. Check fix_mesh for details on fixes. Default: False,
    :param type: bool

    :returns: trimesh.Trimesh - Mesh
    """
    tmp_mesh = trimesh.Trimesh(vertices, faces)
    
    if apply_fixes:
        tmp_mesh.process(validate=1)

    if apply_fixes or normals is None:
        normals = compute_normal(tmp_mesh.vertices, tmp_mesh.faces)

    mesh = trimesh.Trimesh(tmp_mesh.vertices, tmp_mesh.faces, vertex_normals=normals)
    return mesh



def fix_trimesh(mesh, resolution: float = 1.0):
    """
    Applies a predefined set of fixes to the mesh, and converts it to a
    specified resolution. These fixes include removing duplicated vertices wihin
    a certain threshold, removing degenerate triangles, splitting longer edges to
    a given target length, and collapsing shorter edges.

    :param name: mesh - Mesh
    :param type: trimesh.Trimesh
    :param name: resolution - Maximum size of edge in the mesh
    :param tyoe: float

    :returns: trimesh.Trimesh - mesh with all fixes applied
    """
    target_len = resolution
    mesh, _ = trimesh.remove_duplicated_vertices(mesh, 0.001)

    count = 0
    mesh, __ = trimesh.remove_degenerated_triangles(mesh, 100)
    mesh, __ = trimesh.split_long_edges(mesh, target_len)
    num_vertices = mesh.num_vertices
    while True:
        mesh, __ = trimesh.collapse_short_edges(mesh, 1e-6)
        mesh, __ = trimesh.collapse_short_edges(
            mesh, target_len, preserve_feature=True)
        mesh, __ = trimesh.remove_obtuse_triangles(mesh, 150.0, 100)
        if mesh.num_vertices == num_vertices:
            break

        num_vertices = mesh.num_vertices
        #print("#v: {}".format(num_vertices));
        count += 1
        if count > 10: break

    mesh = trimesh.resolve_self_intersection(mesh)
    mesh, __ = trimesh.remove_duplicated_faces(mesh)
    mesh = trimesh.compute_outer_hull(mesh)
    mesh, __ = trimesh.remove_duplicated_faces(mesh)
    mesh, __ = trimesh.remove_obtuse_triangles(mesh, 179.0, 5)
    mesh, __ = trimesh.remove_isolated_vertices(mesh)
    mesh, _ = trimesh.remove_duplicated_vertices(mesh, 0.001)

    return mesh

def find_nearest_atom(coords, res_id, new_verts):
    j = len(res_id)
    el = res_id[j-1]
    full_res_id = []
    
    from scipy import spatial
    #airports = new_verts
    tree = spatial.KDTree(coords)
    
    for i in range(len(new_verts)):
        dist, loc = tree.query(new_verts[i])
        lenr = len(res_id)
        full_res_id.append(res_id[loc])
        
    return full_res_id