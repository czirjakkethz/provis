import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD, GRO, XTC
from provis.utils.atminfo import import_atm_size_info

print(mda.__version__)
my_univ = mda.Universe("data/2fd7.mol2")
print(my_univ)

print(my_univ.atoms.bonds)
for bond in my_univ.atoms.bonds:
    for atom in bond:
        print(atom.position)

# print(my_univ.trajectory[0])

size, col ,vw = import_atm_size_info(1)
# print(vw)
vw["O.3"] = vw["O"] 
vw["N.3"] = vw["N"]
vw["S.3"] = vw["S"]
vw["C.ar"] = vw["C"]
vw["N.pl3"] = vw["N"]
vw["C.2"] = vw["C"]
vw["N.am"] = vw["N"]
vw["C.cat"] = vw["C"]
vw["C.3"] = vw["C"]
vw["O.2"] = vw["O"]
vw["O.co2"] = vw["O"]
bonds = my_univ.atoms.guess_bonds(vdwradii=vw)
print(bonds)