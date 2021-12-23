import pyvista as pv

from provis.src.processing.data_handler import DataHandler
from provis.utils.name_checker import NameChecker

class StickPoint:
    """
    The StickPoint class is used to visualize the structural information of the given molecule. One can easily plot the atoms, residues, bonds or any combination of these.
    """
    def __init__(self):
        """
        Initializes StickPoint class with given filename.
        Creates internal data structures; a DataHandler class instance and loads all the atomic information required for plotting.
        """
        self._id, self._name, self._base_path = NameChecker.return_all()
        # create "brain" of plotting class
        self._dh = DataHandler()
          
        atom_data = self._dh.get_atoms() # second arg: 1 = showsolvent
        ## return list of spheres (meshes) and colors for the spheres
        self._atoms, self._col_a = self._dh.get_atom_mesh(atom_data, vw=0) # second arg: 1 = showvw spheres instead of "normal" radius
        ## return list of lines (meshes)
        self._bonds = self._dh.get_bond_mesh()
        ## return list of spheres (meshes) and colors for the spheres
        res_data = self._dh.get_residues()
        self._residues, self._col_r = self._dh.get_residue_mesh(res_data)
        
        self._atoms_vw, self._col_vw = self._dh.get_atom_mesh(atom_data, vw=1)
        

    def manual_plot(self, atoms=0, col_a=0, box=0, bonds=0, vw=0, residues=0, col_r=0, res=0, outname=0):
        """
        Plot stick and point model
        
        :param name: atoms - List of pyvista Shperes representing each atom, default: 0.
        :param type: list, optional
        :param name: col_a - List of colors for each atom, default: 0.
        :param type: list, optional
        :param name: box - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: bonds - List of pyvista Lines representing each bond, default: 0.
        :param type: list, optional
        :param name: vw - If True styling for Van-der-Waals plotting set. Vw atomic objects still have to be passed under 'atoms' variable.
        :param type: bool, optional
        :param name: res - List of pyvista Shperes representing each residue, default: 0.
        :param type: list, optional
        :param name: col_r - List of colors for each residue, default: 0.
        :param type: list, optional
        :param name: res - Specified residues will be plotted with a bounding box around them.
        :param type: Residue, optional
        :param name: outname - save image of plot to specified filename. Will appear in data/img directory. default: data/img/{self._name}_stick_point.
        :param type: string, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot
        """
        # Use 3 lights because it's a bit brighter
        pl = pv.Plotter(lighting=None)
        pl.background_color = 'grey'
        pl.enable_3_lights()

        # adding the spheres (by atom type) one at a time
        opacity = 1 - vw*0.4
        style = 'surface'
        if atoms != 0:
            if vw:
                style='wireframe'
            for j, mesh in enumerate(atoms):
                pl.add_mesh(mesh, color=col_a[j], opacity=opacity, smooth_shading=True, style=style)

        # adding the bonds one at a time
        if bonds != 0:
            for b in bonds:
                pl.add_mesh(b, color='w', render_lines_as_tubes=True, line_width=5)
        # adding the spheres (by residue) one at a time
        # only executes if residue information provided
        if residues != 0:
            for k, mesh in enumerate(residues):
                pl.add_mesh(mesh, color=col_r[k], opacity=0.2)
        
        # if specified add bounding box
        if box:
            pl.add_bounding_box(color='white', corner_factor=0.5, line_width=1)
        if res:
            res_list, chain_list, pad = res.get_res_info()
            for i, r in enumerate(res_list):
                chain = chain_list[i]
                residues = self._dh.get_structure().get_residues()
                residues_list = list(residues)
                res_name = residues_list[r + 1].get_resname()
                d = (self._dh._res_size_dict[res_name] + pad) * 2
                x, y, z = d,d,d
                pl.add_mesh(pv.Cube(center=self._dh.get_residue_info(r, chain,'com'), x_length=x, y_length=y, z_length=z), style='wireframe', show_edges=1, line_width=5)
        
        # save a screenshot
        if not outname:
            new_name = self._name.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_stick_point.png'
        pl.show(screenshot=outname)


    def plot(self, atoms=0, box=0, bonds=0, vw=0, residues=0, res=None, outname=0, title=None):
        """
        Plot stick and point model
        
        :param name: atoms - Plot atoms, default: 0.
        :param type: bool, optional
        :param name: box, optional - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: bonds, optional - Plot bond, default: 0.
        :param type: bool, optional
        :param name: vw - Plot Wan-der-Waals radii instead of atomic radii.
        :param type: bool, optional
        :param name: residues, optional - Plot residue, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_stick_point.png.
        :param type: string, optional
        :param name: title - Title of the plot window. Defaults to None.
        :param type: str, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot
        """
        # Use 3 lights because it's a bit brighter
        pl = pv.Plotter(lighting=None)
        pl.background_color = 'grey'
        pl.enable_3_lights()

        # adding the spheres (by atom type) one at a time
        opacity = 1 - vw*0.4
        style = 'surface'
        if atoms:
            if vw:
                style='wireframe'
                for j, mesh in enumerate(self._atoms_vw):
                    pl.add_mesh(mesh, color=self._col_vw[j], opacity=opacity, smooth_shading=True, style=style)
            else:
                for j, mesh in enumerate(self._atoms):
                    pl.add_mesh(mesh, color=self._col_a[j], opacity=opacity, smooth_shading=True, style=style)

        # adding the bonds one at a time
        if bonds:
            for b in self._bonds:
                pl.add_mesh(b, color='w', line_width=5, render_lines_as_tubes=True)
                # TODO: figure out why  render_lines_as_tubes=True, crashes
        # adding the spheres (by residue) one at a time
        # only executes if residue information provided
        if residues:
            for k, mesh in enumerate(self._residues):
                pl.add_mesh(mesh, color=self._col_r[k], opacity=0.2)
        
        # if specified add bounding box
        if box:
            pl.add_bounding_box(color='white', corner_factor=0.5, line_width=1)
        
        if res:
            res_list, chain_list, pad = res.get_res_info()
            for i, r in enumerate(res_list):
                chain = chain_list[i]
                residues = self._dh.get_structure().get_residues()
                residues_list = list(residues)
                res_name = residues_list[r + 1].get_resname()
                d = (self._dh._res_size_dict[res_name] + pad) * 2
                x, y, z = d,d,d
                pl.add_mesh(pv.Cube(center=self._dh.get_residue_info(r, chain,'com'), x_length=x, y_length=y, z_length=z), style='wireframe', show_edges=1, line_width=5)

        # save a screenshot
        if not outname:
            new_name = self._name.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_stick_point.png'
        pl.show(screenshot=outname, title=title)


    def plot_stick_point(self, box=0, r=None, outname=0):
        """
        Plot stick and point model
        
        :param name: box, optional - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_stick_point.png.
        :param type: string, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot.
        """

        self.plot(atoms=1, box=box, vw=0, bonds=1, residues=0, res=r, outname=outname, title="Stick Point")
        
    def plot_atoms(self, box=0, r=0, outname=0):
        """
        Plot stick and point model
        
        :param name: box, optional - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_atoms.png.
        :param type: string, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot.
        """
        # save a screenshot
        if not outname:
            new_name = self._name.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_atoms.png'
        self.plot(atoms=1, box=box, vw=0, bonds=0, residues=0, res=r, outname=outname, title="Atoms")
        
    def plot_vw(self, box=0, r=0, outname=0):
        """
        Plot Van-der-Waals atoms
        
        :param name: box, optional - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_atoms.png.
        :param type: string, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot.
        """
        # save a screenshot
        if not outname:
            new_name = self._name.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_vw.png'
        self.plot(atoms=1, box=box, vw=1, bonds=0, residues=0, res=r, title="Van-der-Waals")
        
    def plot_bonds(self, box=0, r=0, outname=0):
        """
        Plot bonds only
        
        :param name: box, optional - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_atoms.png.
        :param type: string, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot.
        """
        # save a screenshot
        if not outname:
            new_name = self._name.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_bonds.png'
        self.plot(atoms=0, box=box, vw=0, bonds=1, residues=0, res=r, outname=outname, title="Bonds")

