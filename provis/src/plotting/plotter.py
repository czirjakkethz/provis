import enum
from provis.src.processing.protein import Protein
import pyvista as pv
import numpy as np
from os.path import exists

from provis.src.processing.surface_handler import SurfaceHandler
from provis.src.processing.name_checker import NameChecker

class Plotter:
    """
    The Plotter class is used for plotting structural and surface information of one or more proteins.
    
    Add Proteins to the Plotter to plot them next to one another. While it is possible to add more than two proteins to one Plotter class it is discouraged, as the window will get cluttered.
    
    The class can visualize two kinds of surfaces:
     - a chemically accurate surface created by the MSMS binary.
     - a good approximation of the surface computed natively with o3d and trimesh. (MSMS does not have to be installed for this option. It is fast, but less precise.)
    
    Choose between the two by setting the msms Boolean variable. (Default: True, corresponding to the MSMS binary option.)
    
    """
    def __init__(self, prot: Protein, prot2=None, msms=True, notebook=False, plot_solvent=False):
        """
        In this constructor no mesh information is computed, simply preparations are made. Meshes are loaded in the plotting function.
        
        A NameChecker object is required for initialization, as this is how the program finds the desired pdb file, molecule.
        If nothing else is passed the NameChecker object will be used to initialize the other internal objects of the surface class.
        
        Apart from the required NameChecker object one can also pass a SurfaceHandler for even more control.
        
        Parameters:
            prot: Protein
                Instance of a Protein class. All information to be plotted is taken from this class. 
            prot2: Protein
                Second instance of a Protein class.
            msms: bool, optional
                If True plot msms binary version of surface. If False plot the native (non-binary) surface. Default: True. 
            density: float, optional
                sampling density used in msms binary. Also needed to load the face and vert files, as their (file)names include the density 
            notebook: bool, optional 
                Needs to be set to true to work in a notebook environment. Defualts to False. 
            msms: bool, optional
                Set to True when running the msms version. Only used to save image with "msms" at end of filename ({root directory}/data/img/{pdb_id}_{plot type}_msms.png). Default: False. 
            plot_solvent: bool, optional
                If True solvent molecules will also be plotted. Default: False.
        """
        if notebook:
            pv.set_jupyter_backend('panel')
        self._proteins = [prot]
        if prot2:
            self._proteins.append(prot2)
            
        self._path, self._out_path, self._base_path, self._mesh_path = prot._name_checker.return_all()        
        self._msms = msms
        self._notebook = notebook
        self._shading = not self._notebook
        self._solvent = plot_solvent

    def add_protein(self, protein: Protein):
        """
        Add another Protein to the (internal list of the) Plotter object.

        Parameters:
            protein: Protein 
                Instance of Protein class.
        """
        
        self._proteins.append(protein)

    def manual_plot(self, box=False, res=None, outname=None, atoms=None, col_a=None, bonds=None, vw=0, residues=None, col_r=None, bb=None, camera=None):
        """
        Plots list of meshes directly. One can get these meshes from the DataHandler class.
        
        Parameters:
            box: bool, optional
                If True bounding box also visualized. Default: False.
            res: list, optional
                List of pyvista Shperes representing each residue. Default: None.
            outname: string, optional
                save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{self._out_path}_stick_point.
            atoms: list, optional
                List of pyvista Shperes representing each atom. Default: None.
            col_a: list, optional
                List of colors for each atom. Default: None.
            bonds: list, optional
                List of pyvista Lines representing each bond. Default: None.
            vw: bool, optional
                If True styling for Van-der-Waals plotting set. Vw atomic objects still have to be passed under 'atoms' variable.
            col_r: list, optional
                List of colors for each residue. Default: None.
            res: Residue, optional
                Specified residues will be plotted with a bounding box around them.
            bb: pyvista.Spline, optional
                Spline describing the back-bone of the protein. Default: None.
            camera: pyvista.Camera, optional
                Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to [0, 4 * "max distance from the center", 0] (see: https://pro-vis.readthedocs.io/en/latest/tutorial.html for more detail). Default: None.
        
        Returns: 
            Pyvista.Plotter window
                Window with interactive plot
        """
        
        # plot
        plot_size = len(self._proteins)
        pl = pv.Plotter(notebook=self._notebook, shape=(1, plot_size))
        print("Surface plotter created...")
        for p_id, prot in enumerate(self._proteins):
            pl.subplot(0, p_id)
            pl.background_color = 'grey'
            pl.enable_3_lights()

            if camera: 
                pl.camera = camera
            else:
                self._cam_pos = prot._data_handler._cam_pos
                pl.camera.position = self._cam_pos
                
            # adding the spheres (by atom type) one at a time
            opacity = 1 - vw*0.4
            style = 'surface'
            if atoms != 0:
                if vw:
                    style='wireframe'
                for j, mesh in enumerate(atoms):
                    pl.add_mesh(mesh, color=col_a[j], opacity=opacity, smooth_shading=self._shading, style=style)

            # adding the bonds one at a time
            if bonds:
                for b in bonds:
                    pl.add_mesh(b, color='w', render_lines_as_tubes=True, line_width=5)
            # adding the spheres (by residue) one at a time
            # only executes if residue information provided
            if residues:
                for k, mesh in enumerate(residues):
                    pl.add_mesh(mesh, color=col_r[k], opacity=0.2)
            
            
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
            
            # if specified add bounding box
            if box:
                pl.add_bounding_box(color='white', corner_factor=0.5, line_width=1)
        
        # save a screenshot
        if not outname:
            new_name = self._out_path.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_stick_point.png'
        
        print("Showing plot")
        pl.show(screenshot=outname, title='Provis')


    def plot_structure(self, box=0, res=None, outname=0, atoms=0, bonds=0, vw=0, residues=0, bb=0, title="Structure", camera=None, dynamic=False):
        """
        This member function is called by all the others. Using this function you can plot any combination of the results gotten from the specialized member functions. For example you could plot the atoms and the backbone of the protein in the same plot.
        
        All information to be plotted is already computed. This function simply dictates what is to be plotted.
        
        Parameters:
            box: bool, optional
                ptional - If True bounding box also visualized. Default: None.
            res: Residue, optional
                Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
            outname: string, optional
                Save image of plot to specified filename. Will appear in data/img directory. Defaults to {root directory}/data/img/{pdb_id}_{model_id}_stick_point.png. If Structure class was initialized with msms=True then output will have "_msms.png" as the ending.
            atoms: bool, optional
                Plot atoms. Default: None.
            bonds: int, optional
                ptional - Plot bond. If zero or undefined then it does not plot the bonds, if 1 it plots all bonds uniformly, if 2 it plots colorful bonds (see data_handler). Default: None.
            vw: bool, optional
                Plot Wan-der-Waals radii instead of atomic radii.
            residues: bool, optional
                ptional - Plot residue. Default: None.
            bb: bool, optional
                If True backbone of protein is plotted. Default: False.
            title: str, optional
                Title of the plot window. Default: Structure.
            camera: pyvista.Camera, optional
                Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to [0, 4 * "max distance from the center", 0] (see: https://pro-vis.readthedocs.io/en/latest/tutorial.html for more detail). Default: None.
            dynamic: bool, optional
                Set to True if you are plotting a dynamic model. Default: False.
        
        Returns: 
            Pyvista.Plotter window
                Window with interactive plot
        """
        
        # plot
        plot_size = len(self._proteins)
        pl = pv.Plotter(notebook=self._notebook, shape=(1, plot_size))
        print("Surface plotter created...")
        for p_id, prot in enumerate(self._proteins):
            pl.subplot(0, p_id)
            pl.background_color = 'grey'
            pl.enable_3_lights()
             
            atom_data = prot._data_handler.get_atoms(show_solvent=self._solvent, model_id=prot._model_id) # second arg: 1 = show_solvent
            self._cam_pos = prot._data_handler._cam_pos
            print("Structure plotter created for model id: ", prot._model_id)
                
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
                    _atoms_vw, _col_vw, _ = prot._data_handler.get_atom_mesh(atom_data, vw=1)
                    style='wireframe'
                    for j, mesh in enumerate(_atoms_vw):
                        pl.add_mesh(mesh, color=_col_vw[j], opacity=opacity, smooth_shading=self._shading, style=style)
                    print("Van-der-Waals atoms added...")
                else:
                    ## return list of spheres (meshes) and colors for the spheres
                    _atoms, _col_a, _ = prot._data_handler.get_atom_mesh(atom_data, vw=0) # second arg: 1 = showvw spheres instead of "normal" radius
                    for j, mesh in enumerate(_atoms):
                        pl.add_mesh(mesh, color=_col_a[j], opacity=opacity, smooth_shading=self._shading, style=style)
                    print("Atoms added...")

            if bonds:
                # adding the bonds one at a time
                bond_mesh_name = self._mesh_path + "_" + str(prot._model_id) + "_bonds"
                if bonds == 1:
                    bond_mesh_name += ".vtk"
                    if exists(bond_mesh_name):
                        mesh = pv.PolyData(bond_mesh_name)
                    else:
                        print("Calculating Bond mesh")
                        ## return list of lines (meshes)
                        _bonds, _bond_col, _ = prot._data_handler.get_bond_mesh(model_id=prot._model_id)
                        mesh = pv.PolyData()
                        for b in _bonds:
                            mesh += b
                        mesh.save(bond_mesh_name)
                    pl.add_mesh(mesh, color="w", line_width=5, render_lines_as_tubes=True)
                if bonds == 2:
                    _bonds, _bond_col, _ = prot._data_handler.get_bond_mesh(model_id=prot._model_id)
                    for b, c in zip(_bonds, _bond_col):
                        pl.add_mesh(b, color=c, line_width=5, render_lines_as_tubes=True)
                    
                print("Bonds added...")
                    
                    
            # adding the spheres (by residue) one at a time
            # only executes if residue information provided
            if residues:
                ## return list of spheres (meshes) and colors for the spheres
                res_data = prot._data_handler.get_residues(model_id=prot._model_id)
                _residues, _col_r, _ = prot._data_handler.get_residue_mesh(res_data)
                
                for k, mesh in enumerate(_residues):
                    pl.add_mesh(mesh, color=_col_r[k], opacity=0.2)
                print("Residues added...")

            if bb:
                bb_mesh = prot._data_handler.get_backbone_mesh(model_id=prot._model_id)
                pl.add_mesh(bb_mesh, render_lines_as_tubes=True, line_width=10)
                print("Back-bone added...")
                
            if res:
                res_list, chain_list, pad = res.get_res_info()
                for i, r in enumerate(res_list):
                    chain = chain_list[i]
                    residues_ = prot._data_handler.get_structure().get_residues()
                    residues_list = list(residues_)
                    res_name = residues_list[r + 1].get_resname()
                    d = (prot._data_handler._res_size_dict[res_name] + pad) * 2
                    x, y, z = d,d,d
                    
                    center = prot._data_handler.get_residue_info(r, chain,'com')
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
            
            ident = '_structure'
            if outname:
                ident = outname
            outname = self._base_path + 'data/img/' + new_name + "_" + str(prot._model_id) + ident + ending 
            
        print("Showing plot")  
        pl.show(screenshot=outname, title=title)


    def plot_stick_point(self, box=0, res=None, outname=0, camera=None):
        """
        Plot stick and point model of the protein. Atoms are spheres, bonds are tubes.
        
        Parameters:
            box: bool, optional
                ptional - If True bounding box also visualized. Default: None.
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

        self.plot_structure(atoms=1, box=box, vw=0, bonds=1, residues=0, res=res, outname=outname, title="Stick Point", camera=camera)
        
    def plot_atoms(self, box=0, res=0, outname=None, camera=None):
        """
        Plot the atoms as spheres. Each atom has a radius proportianal to its calculated atomic radius.
        
        Consult https://en.wikipedia.org/wiki/CPK_coloring for the coloring. 
        
        Parameters:
            box: bool, optional
                ptional - If True bounding box also visualized. Default: None.
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
        self.plot_structure(atoms=1, box=box, vw=0, bonds=0, residues=0, res=res, outname=outname, title="Atoms", camera=camera)
        
    def plot_residues(self, box=0, res=0, outname=0, camera=None):
        """
        Plot the residues as Spheres. Each sphere is the approximate size of the radius of the given residue. This plot should only be used to get a general feel for the layout of the protein.
        
        For coloring information please visit: http://acces.ens-lyon.fr/biotic/rastop/help/colour.htm
        
        Parameters:
            box: bool, optional
                ptional - If True bounding box also visualized. Default: None.
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
        self.plot_structure(atoms=0, box=box, vw=0, bonds=0, residues=1, res=res, title="Residues", outname='_residues', camera=camera)
        
    def plot_vw(self, box=0, res=0, outname=0, camera=None):
        """
        Plot Van-der-Waals radius of atoms as spheres. Spheres have a wireframe style to be able to view inner structure as well.
        To plot Van-der-Waals radii as solid spheres use the manual_plot_structure() member function.
        
        Parameters:
            box: bool, optional
                ptional - If True bounding box also visualized. Default: None.
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
        self.plot_structure(atoms=1, box=box, vw=1, bonds=0, residues=0, res=res, title="Van-der-Waals", camera=camera)
        
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
                If True bounding box also visualized. Default: None.
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
        self.plot_structure(atoms=0, box=box, vw=0, bonds=b, residues=0, res=res, outname=outname, title="Bonds", camera=camera)
        
        
    def plot_backbone(self, box=0, res=0, outname=0, camera=None):
        """
        Plots the backbone (roughly the amide bonds) of the protein.
        
        Parameters:
            box: bool, optional
                If True bounding box also visualized. Default: None.
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
        self.plot_structure(atoms=0, box=box, vw=0, bonds=0, residues=0, res=res, bb=1, outname=outname, title="Backbone", camera=camera)


    def plot_surface(self, feature=None, title="Surface", patch=False, box=None, res=None, outname=None, camera=None):
        """
        Plot the surface of the protein. If the surface has already been computed and saved to the default file, then the surface will automatically be loaded from there.
        The surface can be computed either using the msms binary or natively. The msms binary is chemically accurate surface, while the native one is only for visualization purposes.
        
        If you run into any sort of error concerning array size mismatching or of the sort delete all the temporary files and the mesh ({root directory}/data/meshes/{pdb_id}_{model_id}.obj).
        This will force everything to be recomputed and the dimension mismatch should disappear.
        
        Parameters:
            feature: str, optional
                Pass which feature (coloring) you want to plot. Options: hydrophob, shape, charge. Default: None (uniform coloring). 
            title: str, optional
                Title of the plot window. Default: Surface. 
            patch: bool, optional
                If True then coloring will be read in from "root directory"/data/tmp/{pdb_id}.pth file. Default: False. 
            box, optional: bool, optional
                If True bounding box also visualized. Default: None. 
            res: Residue, optional
                Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None. 
            outname: string, optional
                Save image of plot to specified filename. Will appear in data/img directory. Defaults to {root directory}/data/img/{pdb_id}_{model_id}_surface.png. If Surface class was initialized with msms=True then output will have "_msms.png" as the ending. 
            camera: pyvista.Camera, optional
                Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to [0, 3 * "max distance from the center", 0] (see: https://pro-vis.readthedocs.io/en/latest/tutorial.html for more detail). Default: None. 
            
        Returns: 
            Pyvista.Plotter window
                Window with interactive plot.
        """
        print("Calculating surface mesh")
        
        # plot
        plot_size = len(self._proteins)
        pl = pv.Plotter(notebook=self._notebook, shape=(1, plot_size))
        for p_id, prot in enumerate(self._proteins):
            print("Surface plotter created for model id:", prot._model_id)
            mesh, cas = prot._surface_handler.return_mesh_and_color(self._msms, feature=feature, patch=patch, model_id=prot._model_id)
            pl.subplot(0, p_id)
            pl.background_color = 'grey'
            pl.enable_3_lights()
            pl.add_mesh(mesh, scalars=cas, cmap='RdBu', smooth_shading=self._shading, show_edges=False)
            print("Mesh added to plotter...")
            # if specified add bounding box
            
            if res:
                res_list, chain_list, pad = res.get_res_info()
                for i, r in enumerate(res_list):
                    chain = chain_list[i]
                    residues_ = prot._data_handler.get_structure().get_residues()
                    residues_list = list(residues_)
                    res_name = residues_list[r + 1].get_resname()
                    d = (prot._data_handler._res_size_dict[res_name] + pad) * 2
                    x, y, z = d,d,d
                    center = prot._data_handler.get_residue_info(r, chain,'com')
                    # if residue not found 1 is returned. Otherwise the coordinates
                    if center != 1:
                        pl.add_mesh(pv.Cube(center=center, x_length=x, y_length=y, z_length=z), style='wireframe', show_edges=1, line_width=5, smooth_shading=self._shading, color='r')
                print("Residues marked...")

            if box:
                pl.add_bounding_box(color='white', corner_factor=0.5, line_width=1)
                print("Bounding box added...")
        
            if camera: 
                pl.camera = camera
                print("Camera added...")
            else:
                self._cam_pos = prot._data_handler._cam_pos
                pl.camera.position = self._cam_pos
            
            # save a screenshot
            if not outname or outname[0] == '_':
                ending = ".png"
                if self._msms:
                    ending = "_msms" + ending
                    
                new_name = self._path.split('/')
                new_name = new_name[-1].split('.')[0]
                
                ident = '_surface'
                if outname:
                    ident = outname
                outname = self._base_path + 'data/img/' + new_name + "_" + str(prot._model_id) + ident + ending
       
        pl.show(screenshot=outname, title=title)

    def plot_hydrophob(self, box=None, res=None, outname=None, camera=None):
        """
        Plot the hydrophobic features of a protein.

        Parameters:
            box, optional: bool, optional
                If True bounding box also visualized. Default: None. 
            res: Residue, optional
                Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None. 
            outname: string, optional
                save image of plot to specified filename. Will appear in data/img/ directory. Default: data/img/{self._out_path}_surface. 
            camera: pyvista.Camera, optional
                Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to [0, 3 * max_distance_from_center, 0]. Default: None. 
        
        Returns:
            Pyvista.Plotter window
                Window with interactive plot.
        """
        if not outname:
            outname = '_hydrophob'
        self.plot_surface(feature="hydrophob", title="Hydrophob", box=box, res=res, outname=outname, camera=camera)

    def plot_shape(self, box=None, res=None, outname=None, camera=None):
        """
        Plot the shape features of a protein.

        Parameters:
            box, optional: bool, optional
                If True bounding box also visualized. Default: None. 
            res: Residue, optional
                Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None. 
            outname: string, optional
                save image of plot to specified filename. Will appear in data/img/ directory. Default: data/img/{self._out_path}_surface. 
            camera: pyvista.Camera, optional
                Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to [0, 3 * "max distance from the center", 0] (see: https://pro-vis.readthedocs.io/en/latest/tutorial.html for more detail). Default: None. 
        
        Returns:
            Pyvista.Plotter window
                Window with interactive plot.
        """
        
        if not outname:
            outname = '_shape'
        self.plot_surface(feature="shape", title="Shape", box=box, res=res, outname=outname, camera=camera)

    def plot_charge(self, box=None, res=None, outname=None, camera=None):
        """
        Plot the charge features of a protein.

        Parameters:
            box, optional: bool, optional
                If True bounding box also visualized. Default: None. 
            res: Residue, optional
                Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None. 
            outname: string, optional
                save image of plot to specified filename. Will appear in data/img/ directory. Default: data/img/{self._out_path}_surface. 
            camera: pyvista.Camera, optional
                Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to [0, 3 * "max distance from the center", 0] (see: https://pro-vis.readthedocs.io/en/latest/tutorial.html for more detail). Default: None. 
        
        Returns:
            Pyvista.Plotter window
                Window with interactive plot.
        """
        if not outname:
            outname = '_charge'
        self.plot_surface(feature="charge", title="Charge", box=box, res=res, outname=outname, camera=camera)

