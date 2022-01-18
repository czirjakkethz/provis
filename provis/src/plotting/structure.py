import pyvista as pv

from provis.src.processing.data_handler import DataHandler
from provis.src.processing.name_checker import NameChecker

class Structure:
    """
    The Structure class is used to visualize the structural information of the given molecule. One can easily plot the atoms, residues, bonds or any combination of these structures.
    """
    def __init__(self, nc, dh=None, plot_solvent=False, notebook=False):
        """
        The constructor creates the internal data structures and loads all the atomic information required for plotting.
        
        A NameChecker object is required for initialization, as this is how the program finds the desired pdb file/molecule.
        If nothing else is passed the NameChecker object will be used to initialize the other internal objects of the structure class.
        
        Apart from the required NameChecker object one can also pass a DataHandler for even more control.
        
        :param name: nc - Instance of a NameChecker class. Used to pass the pdb file name and paths.
        :param type: NameChecker
        :param name: dh - Instance of a DataHandler class. Used to retrieve atom-positional information. Default: None. If None a new DataHandler variable will be initialized with "nc".
        :param type: DataHandler, optional
        :param name: plot_solvent - If True solvent molecules will also be plotted. Default: False.
        :param type: bool, optional
        :param name: notebook - Needs to be set to true to work in a notebook environment. Defualts to False.
        :param type: bool, optional 
        """
        self._notebook = notebook
        self._shading = not self._notebook

        self._id, self._name, self._base_path = nc.return_all()
        # create "brain" of plotting class
        if not dh:
            self._dh = DataHandler(nc)
        else:
            self._dh = dh
        atom_data = self._dh.get_atoms(show_solvent=plot_solvent) # second arg: 1 = show_solvent
        ## return list of spheres (meshes) and colors for the spheres
        self._atoms, self._col_a = self._dh.get_atom_mesh(atom_data, vw=0) # second arg: 1 = showvw spheres instead of "normal" radius
        ## return list of lines (meshes)
        self._bonds, self._bond_col = self._dh.get_bond_mesh()
        ## return list of spheres (meshes) and colors for the spheres
        res_data = self._dh.get_residues()
        self._residues, self._col_r = self._dh.get_residue_mesh(res_data)
        
        self._atoms_vw, self._col_vw = self._dh.get_atom_mesh(atom_data, vw=1)
        

    def manual_plot(self, box=0, res=0, outname=0, atoms=0, col_a=0, bonds=0, vw=0, residues=0, col_r=0, bb=0, camera=None):
        """
        Plot stick and point model. In this function one can pass all the desired meshes to be plotted. One can get these meshes from the DataHandler class.
        
        :param name: box - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - List of pyvista Shperes representing each residue, default: 0.
        :param type: list, optional
        :param name: outname - save image of plot to specified filename. Will appear in data/img directory. default: data/img/{self._name}_stick_point.
        :param type: string, optional
        :param name: atoms - List of pyvista Shperes representing each atom, default: 0.
        :param type: list, optional
        :param name: col_a - List of colors for each atom, default: 0.
        :param type: list, optional
        :param name: bonds - List of pyvista Lines representing each bond, default: 0.
        :param type: list, optional
        :param name: vw - If True styling for Van-der-Waals plotting set. Vw atomic objects still have to be passed under 'atoms' variable.
        :param type: bool, optional
        :param name: col_r - List of colors for each residue, default: 0.
        :param type: list, optional
        :param name: res - Specified residues will be plotted with a bounding box around them.
        :param type: Residue, optional
        :param name: bb - List of coordinates describing the back-bone of the protein, default: 0.
        :param type: bool, optional
        :param name: camera - Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to 'xy'. Default: None.
        :param type: pyvista.Camera, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot
        """
        
        pl = pv.Plotter(notebook=self._notebook)
        pl.background_color = 'grey'
        pl.enable_3_lights()

        # adding the spheres (by atom type) one at a time
        opacity = 1 - vw*0.4
        style = 'surface'
        if atoms != 0:
            if vw:
                style='wireframe'
            for j, mesh in enumerate(atoms):
                pl.add_mesh(mesh, color=col_a[j], opacity=opacity, smooth_shading=self._shading, style=style)

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
                pl.add_mesh(pv.Cube(center=self._dh.get_residue_info(r, chain,'com'), x_length=x, y_length=y, z_length=z), style='wireframe', show_edges=1, line_width=5, color='r')
        
        
        if bb:
            pl.add_mesh(pv.Spline(bb), render_lines_as_tubes=True, smooth_shading=self._shading, line_width=10)
        
        # save a screenshot
        if not outname:
            new_name = self._name.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_stick_point.png'
        
        if camera: 
            pl.camera = camera
        else:
            pl.camera.position = 'xy'
        pl.show(screenshot=outname, title='Provis')


    def plot(self, box=0, res=None, outname=0, atoms=0, bonds=0, vw=0, residues=0, bb=0, title=None, camera=None):
        """
        This member function is called by all the others. Using this function you can plot any combination of the results gotten from the specialized member functions. For example you could plot the atoms and the backbone of the protein in the same plot.
        
        All information to be plotted is already computed. This function simply dictates what is to be plotted.
        
        :param name: box, optional - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_stick_point.png.
        :param type: string, optional
        :param name: atoms - Plot atoms, default: 0.
        :param type: bool, optional
        :param name: bonds, optional - Plot bond. If zero or undefined then it does not plot the bonds, if 1 it plots all bonds uniformly, if 2 it plots colorful bonds (see data_handler). Default: 0.
        :param type: int, optional
        :param name: vw - Plot Wan-der-Waals radii instead of atomic radii.
        :param type: bool, optional
        :param name: residues, optional - Plot residue, default: 0.
        :param type: bool, optional
        :param name: bb - If True backbone of protein is plotted. Default: False.
        :param type: bool, optional
        :param name: title - Title of the plot window. Defaults to None.
        :param type: str, optional
        :param name: camera - Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to 'xy'. Default: None.
        :param type: pyvista.Camera, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot
        """
        pl = pv.Plotter(notebook=self._notebook)
        pl.background_color = 'grey'
        pl.enable_3_lights()

        # adding the spheres (by atom type) one at a time
        opacity = 1 - vw*0.4
        style = 'surface'
        if atoms:
            if vw:
                style='wireframe'
                for j, mesh in enumerate(self._atoms_vw):
                    pl.add_mesh(mesh, color=self._col_vw[j], opacity=opacity, smooth_shading=self._shading, style=style)
            else:
                for j, mesh in enumerate(self._atoms):
                    pl.add_mesh(mesh, color=self._col_a[j], opacity=opacity, smooth_shading=self._shading, style=style)

        # adding the bonds one at a time
        if bonds == 1:
            for b in self._bonds:
                pl.add_mesh(b, color="w", line_width=5, render_lines_as_tubes=True)
        if bonds == 2:
            for b, c in zip(self._bonds, self._bond_col):
                pl.add_mesh(b, color=c, line_width=5, render_lines_as_tubes=True)
                
                
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
                pl.add_mesh(pv.Cube(center=self._dh.get_residue_info(r, chain,'com'), x_length=x, y_length=y, z_length=z), style='wireframe', show_edges=1, line_width=5, smooth_shading=self._shading, color='r')

        if bb:
            bb_mesh = self._dh.get_backbone_mesh()
            pl.add_mesh(pv.Spline(bb_mesh), render_lines_as_tubes=True, line_width=10)
        # save a screenshot
        if not outname:
            new_name = self._name.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_stick_point.png'
           
        if camera: 
            pl.camera.position = camera
        pl.show(screenshot=outname, title=title)


    def plot_stick_point(self, box=0, res=None, outname=0, camera=None):
        """
        Plot stick and point model of the protein. Atoms are spheres, bonds are tubes.
        
        :param name: box, optional - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_stick_point.png.
        :param type: string, optional
        :param name: camera - Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to 'xy'. Default: None.
        :param type: pyvista.Camera, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot.
        """

        self.plot(atoms=1, box=box, vw=0, bonds=1, residues=0, res=res, outname=outname, title="Stick Point", camera=camera)
        
    def plot_atoms(self, box=0, res=0, outname=0, camera=None):
        """
        Plot the atoms as spheres.
        
        Consult https://en.wikipedia.org/wiki/CPK_coloring for the coloring. 
        
        :param name: box, optional - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_atoms.png.
        :param type: string, optional
        :param name: camera - Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to 'xy'. Default: None.
        :param type: pyvista.Camera, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot.
        """
        # save a screenshot
        if not outname:
            new_name = self._name.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_atoms.png'
        self.plot(atoms=1, box=box, vw=0, bonds=0, residues=0, res=res, outname=outname, title="Atoms", camera=camera)
        
    def plot_vw(self, box=0, res=0, outname=0, camera=None):
        """
        Plot Van-der-Waals radius of atoms as wireframe spheres.
        
        :param name: box, optional - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_atoms.png.
        :param type: string, optional
        :param name: camera - Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to 'xy'. Default: None.
        :param type: pyvista.Camera, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot.
        """
        # save a screenshot
        if not outname:
            new_name = self._name.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_vw.png'
        self.plot(atoms=1, box=box, vw=1, bonds=0, residues=0, res=res, title="Van-der-Waals", camera=camera)
        
    def plot_bonds(self, box=0, res=0, outname=0, colorful=False, camera=None):
        """
        Plot only the bonds. By default all bonds will be plotted uniformly. If you want to view the difference in bonds you can set the colorful variable to True.
        
        Single bonds: white
        Double bonds: blue
        Triple bonds: green
        Amide bonds: red
        Aromatic bonds: purple
        Undefined/Anything else: black
        
        :param name: box - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_atoms.png.
        :param type: string, optional
        :param name: colorful - If True bonds will be plotted in a colorful manner. If False all bonds are white. Default: False
        :param type: bool, optional
        :param name: camera - Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to 'xy'. Default: None.
        :param type: pyvista.Camera, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot.
        """
        b = bool(colorful) * 1 + 1
        # save a screenshot
        if not outname:
            new_name = self._name.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_bonds.png'
        self.plot(atoms=0, box=box, vw=0, bonds=b, residues=0, res=res, outname=outname, title="Bonds", camera=camera)
        
        
    def plot_backbone(self, box=0, res=0, outname=0, camera=None):
        """
        Plots the backbone (roughly the amide bonds) of the protein.
        
        :param name: box - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_backbone.png.
        :param type: string, optional
        :param name: camera - Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to 'xy'. Default: None.
        :param type: pyvista.Camera, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot.
        """
        # save a screenshot
        if not outname:
            new_name = self._name.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_backbone.png'
        self.plot(atoms=0, box=box, vw=0, bonds=0, residues=0, res=res, bb=1, outname=outname, title="Backbone", camera=camera)

