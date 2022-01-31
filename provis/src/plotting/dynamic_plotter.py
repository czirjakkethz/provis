import pyvista
import numpy as np
import pyvista as pv
import time

from provis.src.processing.protein import Protein

class DynamicPlotter:
    """
    The DynamicPlotter class, similarly to the Protein class, encapsulates every other class and creates a user friendly way to plot your desired dynamic structure of a protein molecules.
    
    While the class is built similarly to the Protein class it does not use the Protein class itself. This is due to the fact that the Protein class is a rigid class made for a single molecule and
    """
        
    def __init__(self, prot: Protein, msms=True, notebook=False, plot_solvent=False):
        """
        In this constructor no mesh information is computed, simply preparations are made. Meshes are loaded in the plotting function.
        
        A NameChecker object is required for initialization, as this is how the program finds the desired pdb file, molecule.
        If nothing else is passed the NameChecker object will be used to initialize the other internal objects of the surface class.
        
        Apart from the required NameChecker object one can also pass a SurfaceHandler for even more control.
        
        Parameters:
            prot: Protein
                Instance of a Protein class. All information to be plotted is taken from this class. 
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
        self._msms = msms
        self._notebook = notebook
        self._shading = not self._notebook
        self._solvent = plot_solvent
        self._path, self._out_path, self._base_path, self._mesh_path = prot._name_checker.return_all()

        if notebook:
            pv.set_jupyter_backend('panel')
        self._protein = prot
             
        if msms:
            self._num_models = prot.file_converter.decompose_traj(self._path)
        else: 
            i = 0
            for model in prot._data_handler._structure:
                i += 1
            self._num_models = i
        self._cam_pos = [0, 0, 0]
        print("Initialized DynamicPlotter class")
       
    def plot_structure(self, box=False, res=None, outname=None, camera=None, title="Structure", atoms=0, bonds=0, vw=0, residues=0, bb=0):
        """
        Plot the dynamic atom cloud.

        Parameters:
            box: bool, optional
                ptional - If True bounding box also visualized, default: 0.
            res: Residue, optional
                Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
            outname: string, optional
                Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_stick_point.png.
            camera: pyvista.Camera, optional
                Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to [0, 4 * "max distance from the center", 0] (see: https://pro-vis.readthedocs.io/en/latest/tutorial.html for more detail). Default: None.
            title: str, optional
                Title of the plotting window. Default: "Atoms".
        
        Returns: 
            Pyvista.Plotter window
                Window with interactive plot
            str - only if self._notebook is True
                The full path to where the .mp4 video file is stored.

        """            
        # Create a plotter object and initialize first mesh
        plotter = pv.Plotter(notebook=self._notebook, off_screen=False)
        plotter.smooth_shading = True
        plotter.background_color = 'grey'
        plotter.add_title(title)
        # save a screenshot
        if not outname or outname[0] == '_':
            ending = "_animation.mp4"
                
            new_name = self._path.split('/')
            new_name = new_name[-1].split('.')[0]
            
            ident = '_structure'
            if outname:
                ident = outname
            outname = self._base_path + 'data/img/' + new_name + ident + ending 
          
        # Open the movie file
        plotter.open_movie(outname, 1)


        # Update mesh
        i = 0
        plotter.enable_eye_dome_lighting()
        plotter.render()
        plotter.write_frame()
        for model in self._protein._data_handler._structure:

            plotter.clear()
            plotter.smooth_shading = True
            plotter.background_color = 'grey'
            plotter.add_title(title)

            atom_data = self._protein._data_handler.get_atoms(show_solvent=self._solvent, model_id=i) # second arg: 1 = show_solvent
            self._cam_pos = self._protein._data_handler._cam_pos
  
            if camera: 
                plotter.camera = camera
                print("Camera added for model id: ", i)
            else:
                plotter.camera.position = self._cam_pos

            if atoms: 
                # adding the spheres (by atom type) one at a time
                opacity = 1 - vw*0.4
                style = 'surface'
                if vw:
                    _atoms_vw, _, _atom_names = self._protein._data_handler.get_atom_mesh(atom_data, vw=1)
                    style='wireframe'
                    bigmesh = pv.PolyData()
                    colors = []
                    for j, mesh in enumerate(_atoms_vw):
                        bigmesh += mesh
                        colors += [_atom_names[j]] * len(mesh.points)
                    plotter.add_mesh(bigmesh, scalars=colors, style=style)
                    print("Van-der-Waals atoms added for model id: ", i)
                else:
                    ## return list of spheres (meshes) and colors for the spheres
                    _atoms, _, _atom_names= self._protein._data_handler.get_atom_mesh(atom_data, vw=0) 
                    
                    bigmesh = pv.PolyData()
                    colors = []
                    for j, mesh in enumerate(_atoms):
                        bigmesh += mesh
                        colors += [_atom_names[j]] * len(mesh.points)
                    plotter.add_mesh(bigmesh, scalars=colors, style=style)
                    print("Atoms added for model id: ", i)
            

            if bonds:
                ## return list of lines (meshes)
                _bonds, _, _bond_names = self._protein._data_handler.get_bond_mesh(model_id=i)
                # adding the bonds one at a time
                if bonds == 1:

                    bond_mesh_name = self._mesh_path + "_" + str(i) + "_bonds"
                    bond_mesh_name += ".vtk"
                    from os.path import exists
                    if exists(bond_mesh_name):
                        bigmesh = pv.PolyData(bond_mesh_name)
                    else:
                        ## return list of lines (meshes)
                        bigmesh = pv.PolyData()
                        for b in _bonds:
                            bigmesh += b
                        bigmesh.save(bond_mesh_name)

                    plotter.add_mesh(bigmesh, color="w", line_width=5)

                if bonds == 2:
                    bigmesh = pv.PolyData()
                    colors = []
                    for j, mesh in enumerate(_bonds):
                        bigmesh += mesh
                        colors += [_bond_names[j]] * len(mesh.points)
                    plotter.add_mesh(bigmesh, scalars=colors)
                print("Bonds added for model id: ", i)

            # adding the spheres (by residue) one at a time
            # only executes if residue information provided
            if residues:
                ## return list of spheres (meshes) and colors for the spheres
                res_data = self._protein._data_handler.get_residues(model_id=i)
                _residues, _col_r, _res_names = self._protein._data_handler.get_residue_mesh(res_data)
                
                bigmesh = pv.PolyData()
                colors = []
                for j, mesh in enumerate(_residues):
                    bigmesh += mesh
                    colors += [_res_names[j]] * len(mesh.points)
                plotter.add_mesh(bigmesh, scalars=colors)
                print("Residues added for model id: ", i)
             
            if bb:
                bb_mesh = self._protein._data_handler.get_backbone_mesh(model_id=i)
                
                plotter.add_mesh(bb_mesh, line_width=10)
                print("Back-bone added for model id: ", i)   
               
            if res:
                res_list, chain_list, pad = res.get_res_info()

                bigmesh = pv.PolyData()
                colors = []
                for j, r in enumerate(res_list):
                    chain = chain_list[j]
                    residues_ = self._protein._data_handler.get_structure().get_residues()
                    residues_list = list(residues_)
                    res_name = residues_list[r + 1].get_resname()
                    d = (self._protein._data_handler._res_size_dict[res_name] + pad) * 2
                    x, y, z = d,d,d
                    res_exists = (self._protein._data_handler.get_residue_info(r, chain,'com') != 1)
                    if res_exists:
                        bigmesh += pv.Cube(center=self._protein._data_handler.get_residue_info(r, chain,'com'), x_length=x, y_length=y, z_length=z)
                
                center = self._protein._data_handler.get_residue_info(r, chain,'com')
                # if residue not found 1 is returned. Otherwise the coordinates
                if center != 1:
                    plotter.add_mesh(pv.Cube(center=center, x_length=x, y_length=y, z_length=z), style='wireframe', show_edges=1, line_width=5, color='r')
                print("Residues marked for model id: ", i)
 
            # if specified add bounding box
            if box:
                plotter.add_bounding_box(color='white', corner_factor=0.5, line_width=1)
                print("Bounding box added for model id: ", i)
              
            # must update normals when smooth shading is enabled
            plotter.mesh.compute_normals(cell_normals=False, inplace=True)
            plotter.render()
            plotter.write_frame()
            time.sleep(0.5)
            i+=1
        plotter.close()
        if self._notebook:
            return outname 

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

        self.plot_structure(atoms=1, box=box, vw=0, bonds=1, residues=0, res=res, outname=outname, title="Stick Point", camera=camera)
        
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
        self.plot_structure(atoms=1, box=box, vw=0, bonds=0, residues=0, res=res, outname=outname, title="Atoms", camera=camera)
        
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
        self.plot_structure(atoms=0, box=box, vw=0, bonds=0, residues=1, res=res, title="Residues", outname='_residues', camera=camera)
        
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
        self.plot_structure(atoms=0, box=box, vw=0, bonds=b, residues=0, res=res, outname=outname, title="Bonds", camera=camera)
        
        
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
        self.plot_structure(atoms=0, box=box, vw=0, bonds=0, residues=0, res=res, bb=True, outname=outname, title="Backbone", camera=camera)


    def plot_surface(self, feature=None, patch=False, title="Surface", box=None, res=None, outname=None, camera=None):
        """
        Plot the dynamic surface of the molecule.

        Parameters:
            feature: str, optional
                Pass which feature (coloring) you want to plot. Options: hydrophob, shape, charge. Default: None (uniform coloring).
            patch: bool, optional
                If True then coloring will be read in from "root directory"/data/tmp/{pdb_id}.pth file. Default: False.
            box: bool, optional
                ptional - If True bounding box also visualized, default: 0.
            res: Residue, optional
                Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
            outname: string, optional
                Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_stick_point.png.
            camera: pyvista.Camera, optional
                Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to [0, 4 * "max distance from the center", 0] (see: https://pro-vis.readthedocs.io/en/latest/tutorial.html for more detail). Default: None.
            title: str, optional
                Title of the plotting window. Default: "Surface".
                
        Returns: 
            Pyvista.Plotter window
                Window with interactive plot
            str - only if self._notebook is True
                The full path to where the .mp4 video file is stored.
        """
  
        # Create and structured surface
        mesh, cas = self._protein._surface_handler.return_mesh_and_color(msms=self._msms, feature=feature, patch=patch, model_id=0, num_models=self._num_models)
        
            
        # Create a plotter object and initialize first mesh
        plotter = pv.Plotter(notebook=self._notebook, off_screen=False)
        plotter.smooth_shading = True
        plotter.background_color = 'grey'
        plotter.add_title(title)

        plotter.add_mesh(mesh, scalars=cas, cmap='RdBu', show_edges=False)
        
        # save a screenshot
        if not outname or outname[0] == '_':
            ending = "_animation.mp4"
                
            new_name = self._path.split('/')
            new_name = new_name[-1].split('.')[0]
            
            ident = '_surface'
            if outname:
                ident = outname
            outname = self._base_path + 'data/img/' + new_name + ident + ending 
        # Open the movie file
        plotter.open_movie(outname, 1)


        # Update mesh
        i = 0
        plotter.enable_eye_dome_lighting()
        plotter.render()
        plotter.write_frame()
        for model in self._protein._data_handler._structure:
            i += 1
            if i == self._num_models:
                break
            
            mesh, cas = self._protein._surface_handler.return_mesh_and_color(msms=self._msms, feature=feature, patch=patch, model_id=i, num_models=self._num_models)            
            self._cam_pos = self._protein._data_handler._cam_pos
            plotter.clear()
            plotter.camera.position = self._cam_pos
            plotter.smooth_shading = True
            plotter.background_color = 'grey'
            plotter.add_title(title)

            plotter.add_mesh(mesh, scalars=cas, cmap='RdBu', show_edges=False)

            # must update normals when smooth shading is enabled
            plotter.mesh.compute_normals(cell_normals=False, inplace=True)
            plotter.render()
            plotter.write_frame()
            time.sleep(0.5)
            
        plotter.render()

        # Closes and finalizes movie
        plotter.close()
        if self._notebook:
            return outname
        

    def plot_hydrophob(self, box=None, res=None, outname=None, camera=None):
        """
        Plot the hydrophobic features of a protein.

        Parameters:
            box, optional: bool, optional
                If True bounding box also visualized, default: 0. 
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
                If True bounding box also visualized, default: 0. 
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
                If True bounding box also visualized, default: 0. 
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