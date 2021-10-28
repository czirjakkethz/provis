
from Bio.PDB import *
from Bio import SeqIO
import numpy as np


parser = PDBParser() # PERMISSIVE=1

structure = parser.get_structure("2fd7", "2fd7.pdb") #("7o86", "7o86.pdb")

# import nglview as nv

# view = nv.show_biopython(structure)

# view

polypeptide_builder = CaPPBuilder()
counter = 1
for polypeptide in polypeptide_builder.build_peptides(structure):
    seq = polypeptide.get_sequence()
    print(f"Sequence: {counter}, Length: {len(seq)}")
    print(seq)
    counter += 1

from Bio.SeqUtils.ProtParam import ProteinAnalysis
print("\nAnalyzed sequence")
analyzed_seq = ProteinAnalysis(str(seq))
print(analyzed_seq)

print("\nGravy")
print(analyzed_seq.gravy())

print("\nMol weight")
print(analyzed_seq.molecular_weight())

print("\nAmino acid count")
print(analyzed_seq.count_amino_acids())

for atom in structure.get_atoms():
    print(atom)

atoms = structure.get_atoms()
atoms_list = list(atoms)

for atom in atoms:
    print(atom.coord)

    
residue = atom.get_parent()
chain = residue.get_parent()
atom_list = Selection.unfold_entities(chain, "A")
print("Different method")
print(atom_list[0].coord)

num_atoms = len(atoms_list)
coord_array = np.empty((num_atoms, 3))

print('hi')
i = 0
for atom in atoms_list:
    coord_array[i] = atom.coord
    i = i+1
    

#coord_array = np.array([atom.coord for atom in atoms])  

print(coord_array)


import pyvista

#point_cloud = np.random.random((100, 3))
coord_array
pdata = pyvista.PolyData(coord_array)

#pdata['orig_sphere'] = np.arange(num_atoms)

# create many spheres from the point cloud
sphere = pyvista.Sphere(radius=0.2, phi_resolution=10, theta_resolution=10)
pc = pdata.glyph(scale=False, geom=sphere)

pc.plot()#cmap='Reds')

def lines_from_points(points):
    """Given an array of points, make a line set"""
    poly = pyvista.PolyData()
    poly.points = points
    cells = np.full((len(points)-1, 3), 2, dtype=np.int_)
    cells[:, 1] = np.arange(0, len(points)-1, dtype=np.int_)
    cells[:, 2] = np.arange(1, len(points), dtype=np.int_)
    poly.lines = cells
    return poly


line = lines_from_points(coord_array)
line

line["scalars"] = np.arange(line.n_points)
tube = line.tube(radius=0.1)
tube.plot(smooth_shading=True)

"""
https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
https://towardsdatascience.com/visualizing-and-analyzing-proteins-in-python-bd99521ccd 
https://docs.pyvista.org/index.html
https://docs.pyvista.org/examples/00-load/create-spline.html#sphx-glr-examples-00-load-create-spline-py
"""