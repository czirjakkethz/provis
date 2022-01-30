import pyvista as pv
from os.path import exists
from provis.src.processing.data_handler import DataHandler
from provis.src.processing.name_checker import NameChecker

class Structure:
    """
    The Structure class is used to visualize the structural information of the given molecule. One can easily plot the atoms, residues, bonds or any combination of these structures.
    """
    def __init__(self, nc, dh=None, plot_solvent=False, notebook=False, msms=False):
        """
        In this constructor no mesh information is computed, simply preparations are made. Meshes are loaded in the plotting function.
        
        A NameChecker object is required for initialization, as this is how the program finds the desired pdb file, molecule.
        If nothing else is passed the NameChecker object will be used to initialize the other internal objects of the structure class.
        
        Apart from the required NameChecker object one can also pass a DataHandler for even more control.
        
        Parameters:
            nc: NameChecker
                Instance of a NameChecker class. Used to pass the pdb file name and paths.
            dh: DataHandler, optional
                Instance of a DataHandler class. Used to retrieve atom-positional information. Default: None. If None a new DataHandler variable will be initialized with "nc".
            plot_solvent: bool, optional
                If True solvent molecules will also be plotted. Default: False.
            notebook: bool, optional 
                Needs to be set to true to work in a notebook environment. Defualts to False.
            msms: bool, optional
                Set to True when running the msms version. Only used to save image with "msms" at end of filename ({root directory}/data/img/{pdb_id}_{plot type}_msms.png). Default: False.
        """
        self._msms = msms
        self._notebook = notebook
        self._shading = not self._notebook
        self._solvent = plot_solvent
        self._path, self._out_path, self._base_path, self._mesh_path = nc.return_all()
        # create "brain" of plotting class
        if not dh:
            self._dh = DataHandler(nc)
        else:
            self._dh = dh
        self._cam_pos = [0, 0, 0]

    def manual_plot(self, box=0, res=0, outname=0, atoms=0, col_a=0, bonds=0, vw=0, residues=0, col_r=0, bb=0, camera=None):
        """
        Plot stick and point model. In this function one can pass all the desired meshes to be plotted. One can get these meshes from the DataHandler class.
        
        Parameters:
            box: bool, optional
                If True bounding box also visualized, default: 0.
            res: list, optional
                List of pyvista Shperes representing each residue, default: 0.
            outname: string, optional
                save image of plot to specified filename. Will appear in data/img directory. default: data/img/{self._out_path}_stick_point.
            atoms: list, optional
                List of pyvista Shperes representing each atom, default: 0.
            col_a: list, optional
                List of colors for each atom, default: 0.
            bonds: list, optional
                List of pyvista Lines representing each bond, default: 0.
            vw: bool, optional
                If True styling for Van-der-Waals plotting set. Vw atomic objects still have to be passed under 'atoms' variable.
            col_r: list, optional
                List of colors for each residue, default: 0.
            res: Residue, optional
                Specified residues will be plotted with a bounding box around them.
            bb: bool, optional
                List of coordinates describing the back-bone of the protein, default: 0.
            camera: pyvista.Camera, optional
                Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to [0, 4 * "max distance from the center", 0] (see: https://pro-vis.readthedocs.io/en/latest/tutorial.html for more detail). Default: None.
        
        Returns: 
            Pyvista.Plotter window
                Window with interactive plot
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
            pl.add_mesh(bb, render_lines_as_tubes=True, smooth_shading=self._shading, line_width=10)
        
        # save a screenshot
        if not outname:
            new_name = self._out_path.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_stick_point.png'
        
        if camera: 
            pl.camera = camera
        else:
            self._cam_pos = self._sh._dh._cam_pos
            pl.camera.position = self._cam_pos
        pl.show(screenshot=outname, title='Provis')


    def plot(self, box=0, res=None, outname=0, atoms=0, bonds=0, vw=0, residues=0, bb=0, title=None, camera=None, model_id=0, dynamic=False):
        """
        This member function is called by all the others. Using this function you can plot any combination of the results gotten from the specialized member functions. For example you could plot the atoms and the backbone of the protein in the same plot.
        
        All information to be plotted is already computed. This function simply dictates what is to be plotted.
        
        Parameters:
            box: bool, optional
                ptional - If True bounding box also visualized, default: 0.
            res: Residue, optional
                Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
            outname: string, optional
                Save image of plot to specified filename. Will appear in data/img directory. Defaults to {root directory}/data/img/{pdb_id}_{model_id}_stick_point.png. If Structure class was initialized with msms=True then output will have "_msms.png" as the ending.
            atoms: bool, optional
                Plot atoms, default: 0.
            bonds: int, optional
                ptional - Plot bond. If zero or undefined then it does not plot the bonds, if 1 it plots all bonds uniformly, if 2 it plots colorful bonds (see data_handler). Default: 0.
            vw: bool, optional
                Plot Wan-der-Waals radii instead of atomic radii.
            residues: bool, optional
                ptional - Plot residue, default: 0.
            bb: bool, optional
                If True backbone of protein is plotted. Default: False.
            title: str, optional
                Title of the plot window. Defaults to None.
            camera: pyvista.Camera, optional
                Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to [0, 4 * "max distance from the center", 0] (see: https://pro-vis.readthedocs.io/en/latest/tutorial.html for more detail). Default: None.
            model_id: int, optional
                The dynamic model ID of the desired molecule. Count starts at 0. Leave default value for static molecules. Default: 0.
            dynamic: bool, optional
                Set to True if you are plotting a dynamic model. Default: False.
        
        Returns: 
            Pyvista.Plotter window
                Window with interactive plot
        """
        atom_data = self._dh.get_atoms(show_solvent=self._solvent, model_id=model_id) # second arg: 1 = show_solvent
        self._cam_pos = self._dh._cam_pos
        pl = pv.Plotter(notebook=self._notebook)
        print("Structure plotter created for model id: ", model_id)
        pl.background_color = 'grey'
        pl.enable_3_lights()
            
        if camera: 
            pl.camera = camera
            print("Camera added...")
        else:
            pl.camera_position = self._cam_pos
              
        # adding the spheres (by atom type) one at a time
        opacity = 1 - vw*0.4
        style = 'surface'
        if atoms:
            if vw:
                _atoms_vw, _col_vw, _ = self._dh.get_atom_mesh(atom_data, vw=1)
                style='wireframe'
                for j, mesh in enumerate(_atoms_vw):
                    pl.add_mesh(mesh, color=_col_vw[j], opacity=opacity, smooth_shading=self._shading, style=style)
                print("Van-der-Waals atoms added...")
            else:
                ## return list of spheres (meshes) and colors for the spheres
                _atoms, _col_a, _ = self._dh.get_atom_mesh(atom_data, vw=0) # second arg: 1 = showvw spheres instead of "normal" radius
                for j, mesh in enumerate(_atoms):
                    pl.add_mesh(mesh, color=_col_a[j], opacity=opacity, smooth_shading=self._shading, style=style)
                print("Atoms added...")

        if bonds:
            # adding the bonds one at a time
            bond_mesh_name = self._mesh_path + "_" + str(model_id) + "_bonds"
            if bonds == 1:
                bond_mesh_name += ".vtk"
                if exists(bond_mesh_name):
                    mesh = pv.PolyData(bond_mesh_name)
                else:
                    print("Calculating Bond mesh")
                    ## return list of lines (meshes)
                    _bonds, _bond_col, _ = self._dh.get_bond_mesh()
                    mesh = pv.PolyData()
                    for b in _bonds:
                        mesh += b
                    mesh.save(bond_mesh_name)
                pl.add_mesh(mesh, color="w", line_width=5, render_lines_as_tubes=True)
            if bonds == 2:
                _bonds, _bond_col, _ = self._dh.get_bond_mesh()
                for b, c in zip(_bonds, _bond_col):
                    pl.add_mesh(b, color=c, line_width=5, render_lines_as_tubes=True)
                # colors = []
                # mesh = pv.PolyData()
                # bond_col_name = bond_mesh_name + "_col"
                # bond_mesh_name += "_col.stl"         
                # if exists(bond_mesh_name) and exists(bond_col_name):
                #     mesh = pv.PolyData(bond_mesh_name)
                #     with open(bond_col_name, 'rb') as f:
                #         colors = pickle.load(f)
                # else:
                #     # # return list of lines (meshes)
                #     # _bonds, _bond_col, _ = self._dh.get_bond_mesh()
                #     # mesh = pv.PolyData()

                #     # for b, c in zip(_bonds, _bond_col):
                #     #     pl.add_mesh(b, color=c, line_width=5, render_lines_as_tubes=True)
                    
                #     # return list of lines (meshes)
                    
                #     _bonds, _bond_col, _bond_names = self._dh.get_bond_mesh()
                #     print(len(_bond_names))
                #     mesh = pv.PolyData()
                #     colors = []
                #     for j, b in enumerate(_bonds):
                #         mesh += b
                #         colors += [_bond_names[j]] * len(b.points)
                #     import pickle
                #     print(len(mesh.points))
                #     with open(bond_col_name, 'wb') as f:
                #         pickle.dump(colors, f)
                # pl.add_mesh(b, scalars=colors, line_width=5, render_lines_as_tubes=True)
            print("Bonds added...")
                
                
        # adding the spheres (by residue) one at a time
        # only executes if residue information provided
        if residues:
            ## return list of spheres (meshes) and colors for the spheres
            res_data = self._dh.get_residues()
            _residues, _col_r, _ = self._dh.get_residue_mesh(res_data)
            
            for k, mesh in enumerate(_residues):
                pl.add_mesh(mesh, color=_col_r[k], opacity=0.2)
            print("Residues added...")

        if bb:
            bb_mesh = self._dh.get_backbone_mesh()
            pl.add_mesh(bb_mesh, render_lines_as_tubes=True, line_width=10)
            print("Back-bone added...")
             
        if res:
            res_list, chain_list, pad = res.get_res_info()
            for i, r in enumerate(res_list):
                chain = chain_list[i]
                residues = self._dh.get_structure().get_residues()
                residues_list = list(residues)
                res_name = residues_list[r + 1].get_resname()
                d = (self._dh._res_size_dict[res_name] + pad) * 2
                x, y, z = d,d,d
                
                center = self._dh.get_residue_info(r, chain,'com')
                # if residue not found 1 is returned. Otherwise the coordinates
                if center != 1:
                    pl.add_mesh(pv.Cube(center=center, x_length=x, y_length=y, z_length=z), style='wireframe', show_edges=1, line_width=5, smooth_shading=self._shading, color='r')            
            print("Residues marked...")
       
        # if specified add bounding box 
        if box:
            pl.add_bounding_box(color='white', corner_factor=0.5, line_width=1)
            print("Bounding box added...")
            
        # save a screenshot
        if not outname or outname[0] == '_':
            ending = ".png"
            if self._msms:
                ending = "_msms" + ending
                
            new_name = self._out_path.split('/')
            new_name = new_name[-1].split('.')[0]
            
            ident = '_stick_point'
            if outname:
                ident = outname
            outname = self._base_path + 'data/img/' + new_name + "_" + str(model_id) + ident + ending 
          
        print("Showing plot")  
        pl.show(screenshot=outname, title=title)


    def plot_stick_point(self, box=0, res=None, outname=0, camera=None):
        """
        Plot stick and point model of the protein. Atoms are spheres, bonds are tubes.
        
        Parameters:
            box: bool, optional
                ptional - If True bounding box also visualized, default: 0.
            res: Residue, optional
                Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
            outname: string, optional
                Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_stick_point.png.
            camera: pyvista.Camera, optional
                Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to [0, 4 * "max distance from the center", 0] (see: https://pro-vis.readthedocs.io/en/latest/tutorial.html for more detail). Default: None.
        
        Returns: 
            Pyvista.Plotter window
                Window with interactive plot.
        """

        self.plot(atoms=1, box=box, vw=0, bonds=1, residues=0, res=res, outname=outname, title="Stick Point", camera=camera)
        
    def plot_atoms(self, box=0, res=0, outname=None, camera=None):
        """
        Plot the atoms as spheres. Each atom has a radius proportianal to its calculated atomic radius.
        
        Consult https://en.wikipedia.org/wiki/CPK_coloring for the coloring. 
        
        Parameters:
            box: bool, optional
                ptional - If True bounding box also visualized, default: 0.
            res: Residue, optional
                Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
            outname: string, optional
                Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_atoms.png.
            camera: pyvista.Camera, optional
                Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to [0, 4 * "max distance from the center", 0] (see: https://pro-vis.readthedocs.io/en/latest/tutorial.html for more detail). Default: None.
        
        Returns: 
            Pyvista.Plotter window
                Window with interactive plot.
        """
        if not outname:
            outname = '_atoms'
        self.plot(atoms=1, box=box, vw=0, bonds=0, residues=0, res=res, outname=outname, title="Atoms", camera=camera)
        
    def plot_residues(self, box=0, res=0, outname=0, camera=None):
        """
        Plot the residues as Spheres. Each sphere is the approximate size of the radius of the given residue. This plot should only be used to get a general feel for the layout of the protein.
        
        For coloring information please visit: http://acces.ens-lyon.fr/biotic/rastop/help/colour.htm
        
        Parameters:
            box: bool, optional
                ptional - If True bounding box also visualized, default: 0.
            res: Residue, optional
                Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
            outname: string, optional
                Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_residues.png.
            camera: pyvista.Camera, optional
                Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to [0, 4 * "max distance from the center", 0] (see: https://pro-vis.readthedocs.io/en/latest/tutorial.html for more detail). Default: None.
        
        Returns: 
            Pyvista.Plotter window
                Window with interactive plot.
        """
        
        if not outname:
            outname = '_residues'
        self.plot(atoms=0, box=box, vw=0, bonds=0, residues=1, res=res, title="Residues", outname='_residues', camera=camera)
        
    def plot_vw(self, box=0, res=0, outname=0, camera=None):
        """
        Plot Van-der-Waals radius of atoms as spheres. Spheres have a wireframe style to be able to view inner structure as well.
        To plot Van-der-Waals radii as solid spheres use the manual_plot() member function.
        
        Parameters:
            box: bool, optional
                ptional - If True bounding box also visualized, default: 0.
            res: Residue, optional
                Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
            outname: string, optional
                Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_atoms.png.
            camera: pyvista.Camera, optional
                Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to [0, 4 * "max distance from the center", 0] (see: https://pro-vis.readthedocs.io/en/latest/tutorial.html for more detail). Default: None.
        
        Returns: 
            Pyvista.Plotter window
                Window with interactive plot.
        """
        if not outname:
            outname = '_vw'
        self.plot(atoms=1, box=box, vw=1, bonds=0, residues=0, res=res, title="Van-der-Waals", camera=camera)
        
    def plot_bonds(self, box=0, res=0, outname=0, colorful=False, camera=None):
        """
        Plot only the bonds. By default all bonds will be plotted uniformly. 
        
        If the difference in bond types is of interest set the "colorful" variable to True.
        Coloring:
        - Single bonds: white
        - Double bonds: blue
        - Triple bonds: green
        - Amide bonds: red
        - Aromatic bonds: purple
        - Undefined/Anything else: black
        
        Parameters:
            box: bool, optional
                If True bounding box also visualized, default: 0.
            res: Residue, optional
                Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
            outname: string, optional
                Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_atoms.png.
            colorful: bool, optional
                If True bonds will be plotted in a colorful manner. If False all bonds are white. Default: False
            camera: pyvista.Camera, optional
                Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to [0, 4 * "max distance from the center", 0] (see: https://pro-vis.readthedocs.io/en/latest/tutorial.html for more detail). Default: None.
        
        Returns: 
            Pyvista.Plotter window
                Window with interactive plot.
        """
        b = bool(colorful) * 1 + 1
        
        if not outname:
            outname = '_bonds'
        self.plot(atoms=0, box=box, vw=0, bonds=b, residues=0, res=res, outname=outname, title="Bonds", camera=camera)
        
        
    def plot_backbone(self, box=0, res=0, outname=0, camera=None):
        """
        Plots the backbone (roughly the amide bonds) of the protein.
        
        Parameters:
            box: bool, optional
                If True bounding box also visualized, default: 0.
            res: Residue, optional
                Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
            outname: string, optional
                Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_backbone.png.
            camera: pyvista.Camera, optional
                Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to [0, 4 * "max distance from the center", 0] (see: https://pro-vis.readthedocs.io/en/latest/tutorial.html for more detail). Default: None.
        
        Returns: 
            Pyvista.Plotter window
                Window with interactive plot.
        """
        
        if not outname:
            outname = '_backbone'
        self.plot(atoms=0, box=box, vw=0, bonds=0, residues=0, res=res, bb=1, outname=outname, title="Backbone", camera=camera)

