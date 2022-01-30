import pyvista as pv
import numpy as np
import trimesh

from provis.src.processing.surface_handler import SurfaceHandler
from provis.src.processing.name_checker import NameChecker

class Surface:
    """
    The Surface class is used to visualize the surface information of the given molecule.
    
    The class can visualize two kinds of surfaces:
     - a chemically accurate surface created by the MSMS binary.
     - a good approximation of the surface computed natively with o3d and trimesh. (MSMS does not have to be installed for this option. It is fast, but less precise.)
    
    Choose between the two by setting the msms Boolean variable. (Default: True, corresponding to the MSMS binary option.)
    
    """
    def __init__(self, nc, sh=None, msms=True, density=3.0, notebook=False):
        """
        In this constructor no mesh information is computed, simply preparations are made. Meshes are loaded in the plotting function.
        
        A NameChecker object is required for initialization, as this is how the program finds the desired pdb file, molecule.
        If nothing else is passed the NameChecker object will be used to initialize the other internal objects of the surface class.
        
        Apart from the required NameChecker object one can also pass a SurfaceHandler for even more control.
        
        Parameters:
            nc: NameChecker
                Instance of a NameChecker class. Used to pass the pdb file name and paths. 
            sh: SurfaceHandler, optional
                Instance of SurfaceHandler. This class computes the required surface meshes for the plots. The "brain" of the Surface class. Default: None :
                a SurfaceHandler class will be initialized with the NameChecker object. 
            msms: bool, optional
                If True plot msms binary version of surface. If False plot the native (non-binary) surface. Default: True. 
            density: float, optional
                sampling density used in msms binary. Also needed to load the face and vert files, as their (file)names include the density 
            notebook: bool, optional 
                Needs to be set to true to work in a notebook environment. Defualts to False. 
            msms: bool, optional
                Set to True when running the msms version. Only used to save image with "msms" at end of filename ({root directory}/data/img/{pdb_id}_{plot type}_msms.png). Default: False. 
        """
        self._msms = msms
        self._path, self._out_path, self._base_path, self._mesh_path = nc.return_all()
        if not sh:
            self._sh = SurfaceHandler(nc, density=density)
        else:
            self._sh = sh
        self._msms = msms
        self._notebook = notebook
        self._shading = not self._notebook

    def plot(self, feature=None, title="Surface", patch=False, box=None, res=None, outname=None, camera=None, model_id=0):
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
                If True bounding box also visualized, default: 0. 
            res: Residue, optional
                Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None. 
            outname: string, optional
                Save image of plot to specified filename. Will appear in data/img directory. Defaults to {root directory}/data/img/{pdb_id}_{model_id}_surface.png. If Surface class was initialized with msms=True then output will have "_msms.png" as the ending. 
            camera: pyvista.Camera, optional
                Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to [0, 3 * "max distance from the center", 0] (see: https://pro-vis.readthedocs.io/en/latest/tutorial.html for more detail). Default: None. 
            model_id: int, optional
                The dynamic model ID of the desired molecule. Count starts at 0. Leave default value for static molecules. Default: 0. 
        
        Returns: 
            Pyvista.Plotter window
                Window with interactive plot.
        """
        print("Calculating surface mesh")
        # get appropriate mesh and coloring
        mesh, cas = self._sh.return_mesh_and_color(self._msms, feature=feature, patch=patch)
        # plot
        pl = pv.Plotter(notebook=self._notebook)
        print("Surface plotter created...")
        pl.background_color = 'grey'
        pl.enable_3_lights()
        pl.add_mesh(mesh, scalars=cas, cmap='RdBu', smooth_shading=self._shading, show_edges=False)
        print("Mesh added to plotter...")
        # if specified add bounding box
            
        if res:
            res_list, chain_list, pad = res.get_res_info()
            for i, r in enumerate(res_list):
                chain = chain_list[i]
                residues = self._sh._dh.get_structure().get_residues()
                residues_list = list(residues)
                res_name = residues_list[r + 1].get_resname()
                d = (self._sh._dh._res_size_dict[res_name] + pad) * 2
                x, y, z = d,d,d
                center = self._sh._dh.get_residue_info(r, chain,'com')
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
            self._cam_pos = self._sh._dh._cam_pos
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
            outname = self._base_path + 'data/img/' + new_name + "_" + str(model_id) + ident + ending
       
        pl.show(screenshot=outname, title=title)

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
        self.plot(feature="hydrophob", title="Hydrophob", box=box, res=res, outname=outname, camera=camera)

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
        self.plot(feature="shape", title="Shape", box=box, res=res, outname=outname, camera=camera)

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
        self.plot(feature="charge", title="Charge", box=box, res=res, outname=outname, camera=camera)
