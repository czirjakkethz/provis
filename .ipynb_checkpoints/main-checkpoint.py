"""
This file was created by Kristof Czirjak as part of his Bachelor's Thesis - provis
"""
from provis.src.file_converter import FileConverter
from provis.src.stick_point import StickPoint
from provis.src.surface import Surface
from provis.src.residue import Residue
                                    
def main():
    name =  "data/2fd7" # "data/1a3n" # "data/7nkd" #
    density = 3.0
    solvent = 0
    bash = 0
    show_box = 1
    
    fc = FileConverter(name, density, solvent, bash)
    
#    ## plot stick point
#    sp = StickPoint(name)
#    sp.plot_atoms()
#    sp.plot_bonds()
#    sp.plot_vw()
#    sp.plot_stick_point()
#    r = Residue(1)
#    sp.plot(atoms=1, box=1, bonds=1, vw=0, residues=0, res=0)
#

    ## plot surface
    s = Surface(name)
    s.plot_msms_surface(dens=density)
    s.plot_surface()

if __name__ == "__main__":
    main()

# https://github.com/pv/pv-support/issues/374