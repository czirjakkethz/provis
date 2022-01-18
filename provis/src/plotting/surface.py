import pyvista as pv
import numpy as np
import trimesh

from provis.src.processing.surface_handler import SurfaceHandler
from provis.src.processing.name_checker import NameChecker

class Surface:
    """
    The Surface class is used to visualize the surface information of the given molecule.
    
    The class can visualize two kinds of surfaces:
    - the "correct" surface created by the MSMS binary. 
    - a good approximation of the surface computed natively with o3d and trimesh. (MSMS does not have to be installed for this option. It is slightly slower and less precise)
    
    Choose between the two by setting the msms Boolean variable. (True corresponding to the MSMS binary option is default)
    
    """
    def __init__(self, nc, sh=None, msms=False, density=3.0, notebook=False):
        """
        In this constructor no mesh information is computed, simply preparations are made. Meshes are loaded in the plotting function.
        
        A NameChecker object is required for initialization, as this is how the program finds the desired pdb file/molecule.
        If nothing else is passed the NameChecker object will be used to initialize the other internal objects of the surface class.
        
        Apart from the required NameChecker object one can also pass a SurfaceHandler for even more control.
        
        :param name: nc - Instance of a NameChecker class. Used to pass the pdb file name and paths.
        :param type: NameChecker
        :param name: sh - Instance of SurfaceHandler. This class computes the required surface meshes for the plots. The "brain" of the Surface class. Default: None - a SurfaceHandler class will be initialized with the NameChecker object.
        :param type: SurfaceHandler, optional
        :param name: msms - If True plot msms binary version of surface. If False plot the native (non-binary) surface. Default: False.
        :param type: bool, optional
        :param name: density - sampling density used in msms binary. Also needed to load the face and vert files, as their (file)names include the density
        :param type: float, optional
        :param name: notebook - Needs to be set to true to work in a notebook environment. Defualts to False.
        :param type: bool, optional 
        """
        
        self._path, self._out_path, self._base_path = nc.return_all()
        if not sh:
            self._sh = SurfaceHandler(nc, density=density)
        else:
            self._sh = sh
        self._msms = msms
        self._notebook = notebook
        self._shading = not self._notebook

    def plot(self, feature=None, title="Surface", patch=False, box=None, res=None, outname=None, camera=None):
        """
        Plot the surface of the protein.
        
        :param name: feature - Pass which feature (coloring) you want to plot. Options: hydrophob, shape, charge. Default: None (uniform coloring).
        :param type: str, optional
        :param name: title - Title of the plot window. Default: Surface.
        :param type: str, optional
        :param name: patch - If True then coloring will be read in from "root directory"/data/tmp/{pdb_id}.pth file. Default: False.
        :param type: bool, optional
        :param name: box, optional - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - save image of plot to specified filename. Will appear in data/img/ directory. Default: data/img/{self._out_path}_surface.
        :param type: string, optional
        :param name: camera - Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to 'xy'. Default: None.
        :param type: pyvista.Camera, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot.
        """
        # get appropriate mesh and coloring
        mesh, cas = self._sh.return_mesh_and_color(self._msms, feature=feature, patch=patch)
        
        # plot
        pl = pv.Plotter(notebook=self._notebook)
        pl.background_color = 'grey'
        pl.enable_3_lights()
        pl.add_mesh(mesh, scalars=cas, cmap='RdBu', smooth_shading=self._shading, show_edges=False)
        
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

        
        # save a screenshot
        if not outname:
            new_name = self._out_path.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_surface.png'
            
        if camera: 
            pl.camera = camera
        else:
            pl.camera_position = 'xy'
        pl.show(screenshot=outname, title=title)

    def plot_hydrophob(self, box=None, res=None, outname=None, camera=None):
        """
        Plot the hydrophobic features of a protein.

        :param name: box, optional - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - save image of plot to specified filename. Will appear in data/img/ directory. Default: data/img/{self._out_path}_surface.
        :param type: string, optional
        :param name: camera - Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to 'xy'. Default: None.
        :param type: pyvista.Camera, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot.
        """

        # save a screenshot
        if not outname:
            new_name = self._out_path.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_hydrophob.png'
        self.plot(feature="hydrophob", title="Hydrophob", box=box, res=res, outname=outname, camera=camera)

    def plot_shape(self, box=None, res=None, outname=None, camera=None):
        """
        Plot the shape features of a protein.

        :param name: box, optional - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - save image of plot to specified filename. Will appear in data/img/ directory. Default: data/img/{self._out_path}_surface.
        :param type: string, optional
        :param name: camera - Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to 'xy'. Default: None.
        :param type: pyvista.Camera, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot.
        """
        
        # save a screenshot
        if not outname:
            new_name = self._out_path.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_shape.png'
        self.plot(feature="shape", title="Shape", box=box, res=res, outname=outname, camera=camera)

    def plot_charge(self, box=None, res=None, outname=None, camera=None):
        """
        Plot the charge features of a protein.

        :param name: box, optional - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - save image of plot to specified filename. Will appear in data/img/ directory. Default: data/img/{self._out_path}_surface.
        :param type: string, optional
        :param name: camera - Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to 'xy'. Default: None.
        :param type: pyvista.Camera, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot.
        """
    
        # save a screenshot
        if not outname:
            new_name = self._out_path.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_charge.png'
        self.plot(feature="charge", title="Charge", box=box, res=res, outname=outname, camera=camera)