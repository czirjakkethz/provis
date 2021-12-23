"""
https://github.com/bunnech/holoprot/blob/main/holoprot/feat/surface.py is the base for this file. Modifications were made.

Functions to compute features for a patch on the protein surface.

Some of these are borrowed from https://github.com/LPDI-EPFL/masif under the https://github.com/LPDI-EPFL/masif/blob/master/LICENSE license
"""
import numpy as np
from Bio.PDB import Selection
from Bio.PDB.Residue import Residue
from subprocess import PIPE, Popen
import os
from scipy.spatial import KDTree
from typing import List, Tuple, Dict
import trimesh

from provis.src.processing import KD_SCALE, get_residues
from provis.utils.charges_utils import (compute_hbond_helper, compute_satisfied_CO_HN,
                    normalize_electrostatics)
from provis.utils.surface_utils import prepare_trimesh


Surface = Tuple[np.ndarray, np.ndarray, np.ndarray, List[str], Dict[str, str]]

def compute_shape_index(mesh) -> np.ndarray:
    """
    Computes shape index for the patches. Shape index characterizes the shape
    around a point on the surface, computed using the local curvature around each
    point. These values are derived using Trimesh's available geometric
    processing functionality.

    :param name: mesh: The mesh is constructed using information about vertices and faces.
    :param type: Trimesh

    :returns: np.ndarray - Shape index for each vertex
    """
    normals = mesh.vertex_normals
    # n1 = mesh.get_attribute("vertex_nx")
    # n2 = mesh.get_attribute("vertex_ny")
    # n3 = mesh.get_attribute("vertex_nz")
    # normals = np.stack([n1, n2, n3], axis=1)
    # mesh.add_attribute("vertex_mean_curvature")
    # H = mesh.get_attribute("vertex_mean_curvature")
    # mesh.add_attribute("vertex_gaussian_curvature")
    # K = mesh.get_attribute("vertex_gaussian_curvature")
    H = trimesh.curvature.discrete_mean_curvature_measure(mesh, mesh.vertices, 1)
    K = trimesh.curvature.discrete_gaussian_curvature_measure(mesh, mesh.vertices, 1)
    elem = np.square(H) - K
    # In some cases this equation is less than zero, likely due to the method
    # that computes the mean and gaussian curvature. set to an epsilon.
    elem[elem < 0] = 1e-8
    k1 = H + np.sqrt(elem)
    k2 = H - np.sqrt(elem)
    # Compute the shape index
    si = (k1 + k2) / (k1 - k2)
    si = np.arctan(si) * (2 / np.pi)
    return si


def compute_hbonds(vertices: np.ndarray, residues: List[Residue],
                   names: List[str]) -> np.ndarray:
    """Compute H-bond (hydrogen-bond) induced charges at every vertex.
    
    :param name: vertices - Vertices of mesh
    :param type: np.ndarray
    :param name: residues - List of residues to compute
    :param type: List[Residue]
    :param name: names - List of custom names created by output_pdb_as_xyzrn()
    :param type: List[str]
    
    :returns: np.ndarray - Array of bonds, by vertex
    """
    
    residue_dict = {}
    for res in residues:
        chain_id = res.get_parent().get_id()
        if chain_id == "":
            chain_id = " "
        residue_dict[(chain_id, res.get_id())] = res

    atoms = Selection.unfold_entities(residues, "A")
    satisfied_CO, satisfied_HN = compute_satisfied_CO_HN(atoms)

    hbond = np.array([0.0] * len(vertices))
    # Go over every vertex
    for ix, name in enumerate(names):
        fields = name.split("_")
        chain_id = fields[0]
        if chain_id == "":
            chain_id = " "
        if fields[2] == "x":
            fields[2] = " "
        res_id = (" ", int(fields[1]), fields[2])
        aa = fields[3]
        atom_name = fields[4]
        # Ignore atom if it is BB and it is already satisfied.
        if atom_name == "H" and res_id in satisfied_HN:
            continue
        if atom_name == "O" and res_id in satisfied_CO:
            continue
        # Compute the charge of the vertex
        hbond[ix] = compute_hbond_helper(
            atom_name, residue_dict[(chain_id, res_id)], vertices[ix])

    return hbond


def assign_props_to_new_mesh(new_vertices,
                             old_vertices: np.ndarray,
                             old_props: np.ndarray,
                             feature_interpolation: bool = True) -> np.ndarray:
    """
    Assign properties to vertices in modified mesh given the initial mesh. The
    assignment is carried using a KDTree data structure to query nearest points.

    :param name: new_vertices - Vertices on the modified mesh
    :param type: np.ndarray 
    :param name: old_vertices - Vertices on the original mesh
    :param type: np.ndarray 
    :param name: old_props - Property values for each vertex on the original mesh
    :param type: np.ndarray 
    :param name: feature_interpolation - If set to True interpolates features to new vertices.
    :param type: bool, (default True)

    :returns: np.ndarray - Property values for vertices on the modified mesh
    """
    dataset = old_vertices
    testset = new_vertices
    new_props = np.zeros(len(new_vertices))
    if feature_interpolation:
        num_inter = 4  # Number of interpolation features
        # Assign k old vertices to each new vertex.
        kdt = KDTree(dataset)
        dists, result = kdt.query(testset, k=num_inter)
        # Square the distances (as in the original pyflann)
        dists = np.square(dists)
        # The size of result is the same as new_vertices
        for vi_new in range(len(result)):
            vi_old = result[vi_new]
            dist_old = dists[vi_new]
            # If one vertex is right on top, ignore the rest.
            if dist_old[0] == 0.0:
                new_props[vi_new] = old_props[vi_old[0]]
                continue

            total_dist = np.sum(1 / dist_old)
            for i in range(num_inter):
                new_props[vi_new] += (
                    old_props[vi_old[i]] * (1 / dist_old[i]) / total_dist)
    else:
        # Assign k old vertices to each new vertex.
        kdt = KDTree(dataset)
        dists, result = kdt.query(testset)
        new_props = old_props[result]
    return new_props


def compute_hydrophobicity(names: List[str]) -> np.ndarray:
    """
    Compute hydrophobicity value for all vertices on the surface. Each surface
    vertex has a mapping to the corresponding residue from the original protein.
    This is used to assign a hydrophobicity value to each vertex using the Kyte-
    Doolittle scale.

    :param name: names - Identifier names for each vertex in the surface
    :param type: List[str]
    
    :returns: np.ndarray - Hydrophobicity values for each surface vertex
    """
    hp = np.zeros(len(names))
    for ix, name in enumerate(names):
        res_name = name.split("_")[3]
        hp[ix] = KD_SCALE[res_name] / 4.5
    return hp


def compute_charges(vertices: np.ndarray, pdb_id: str, path: str) -> np.ndarray:
    """
    Computes electrostatics for the surface vertices. The function first calls
    the PDB2PQR executable to prepare the pdb file for electrostatics. Poisson-
    Boltzmann electrostatics are computed using APSB executable. Multivalue,
    provided within APSB suite is used to assign charges to each vertex. The
    charges are further normalized.

    :param name: vertices - Surface vertex coordinates
    :param type: np.ndarray
    :param name: pdb_id - PDB ID of the protein
    :param type: str
    :param name: path - Path of the pqr file in the form {path}.pqr
    :param type: str
    
    :returns: np.ndarray - Charge values for each vertex
    """

    charge_file = f"{path}.pqr"
    if not os.path.exists(charge_file):
        raise ValueError(f"Charges cannot be computed. Missing file {pdb_id}.pqr. {pdb_id}")

    chargefile = open(charge_file)
    charges = np.array([0.0] * len(vertices))
    for ix, line in enumerate(chargefile.readlines()):
        try:
            charges[ix] = float(line.split()[8]) # "," for some reson pqr does not have commas, but whitespaces
        except:
            pass
    charges = normalize_electrostatics(charges)
    return charges


def compute_surface_features(surface: Surface,
                             pdb_file: str,
                             path: str,
                             mesh = None,
                             fix_mesh: bool = False,
                             return_mesh: bool = False,
                             pdb_id: str = None) -> Tuple[np.ndarray]:
    """Computes all patch features.

    :param name: surface: Tuple of attributes characterizing the surface. These include vertices,
        faces, normals to each vertex, areas, fdue identifiers for vertices.
    :param type: Surface
    :param name: pdb_file - PDB File containing the atomic coordinates
    :param type: str
    :param name: path - Path to files without extensions. Usually data/tmp/{pdb_id}
    :param type: str
    :param name: mesh - The mesh.
    :param type: trimesh.Trimesh
    :param name: fix_mesh - Whether to fix the mesh by collapsing nodes and edges. Default: False.
    :param type: bool, optional
    :param name: return_mesh - Whether to return the mesh. Default: False.
    :param type: bool, optional
    :param name: pdb_id - PDB id of the associated protein. Default: None
    :param type: str, optional

    :returns: np.ndarray - Shape index
    :returns: np.ndarray - Hydrogen bond induced charges
    :returns: np.ndarray - Hydrophobicity of each surface vertex
    :returns: np.ndarray - Electrostatics of each surface vertex
    """
    residues = get_residues(pdb_file)
    vertices, faces, normals, names, areas = surface
    hbonds = compute_hbonds(vertices, residues, names)
    hydrophob = compute_hydrophobicity(names)
    if mesh is None:
        mesh = prepare_trimesh(vertices, faces, normals=normals, apply_fixes=fix_mesh)
    si = compute_shape_index(mesh)
    
    if fix_mesh or (len(mesh.vertices) != len(vertices)):
        hbonds = assign_props_to_new_mesh(
            new_vertices=mesh.vertices,
            old_vertices=vertices,
            old_props=hbonds,
            feature_interpolation=True)
        hydrophob = assign_props_to_new_mesh(
            new_vertices=mesh.vertices,
            old_vertices=vertices,
            old_props=hydrophob,
            feature_interpolation=True)

    charges = compute_charges(mesh.vertices, pdb_id, path)

    if si is None or not np.isfinite(si).all():
        raise ValueError(f"{pdb_id}: Shape index failed nan check.")

    if hbonds is None or not np.isfinite(hbonds).all():
        raise ValueError(f"{pdb_id}: HBonds failed nan check.")

    if hydrophob is None or not np.isfinite(hydrophob).all():
        raise ValueError(f"{pdb_id}: Hydrophob failed nan check.")

    if charges is None or not np.isfinite(charges).all():
        raise ValueError(f"{pdb_id}: Charges failed nan check.")

    node_feats = (si, hbonds, hydrophob, charges)
    
    if return_mesh:
        return node_feats, mesh

    return node_feats
