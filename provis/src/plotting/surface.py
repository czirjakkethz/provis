import pyvista as pv
import numpy as np
import trimesh

from provis.src.processing.surface_handler import SurfaceHandler
from provis.src.processing.name_checker import NameChecker

class Surface:
    """
    The Surface class is used to visualize the surface information of the given molecule.
    
    All member functions of the class are a "two in one" function:
    - plot the surface created by the MSMS binary. 
    - create a good approximation of the surface natively with o3d and trimesh and plot it. MSMS does not have to be installed for this option. (It is slightly slower and less precise)
    
    Choose between the two by setting the msms Boolean variable (True corresponding to the MSMS binary option is default)
    
    """
    def __init__(self, nc, dh=None, sh=None, density=None, msms=False, notebook=False):
        """
        Initialize Surface class with given filename. Creates internal data structures; a DataHandler to extract basic surface information and stores it in self._atmsurf (this is a list of Spheres for each atom roughly equating the Van-der-Waals radius).
        
        :param name: nc - Instance of a NameChecker class. Used to pass the pdb file name and paths.
        :param type: NameChecker
        :param name: msms - If True plot msms binary version of surface. If False plot the native (non-binary) surface. Default: False.
        :param type: bool, optional
        :param name: density - sampling density used in msms binary. Also needed to load the face and vert files, as their (file)names include the density
        :param type: float, optional
        :param name: notebook - Needs to be set to true to work in a notebook environment. Defualts to False.
        :param type: bool, optional 
        :param name: dh - Instance of DataHandler. To be passed to self._sh (SurfaceHandler), not needed otherwise. Default: None.
        """
        
        self._path, self._out_path, self._base_path = nc.return_all()
        if density:
            self._density = density
        if not sh:
            self._sh = SurfaceHandler(nc, density=density, dh=dh)
        else:
            self._sh = sh
        self._msms = msms
        self._notebook = notebook
        self._shading = not self._notebook

    def plot(self, feature=None, title=None, patch=False, box=None, res=None, outname=None):
        """
        Plot the surface of protein.
        
        :param name: feature - Pass which feature (coloring) you want to plot. Options: hydrophob, shape, charge. Default: None (uniform coloring).
        :param type: str, optional
        :param name: title - Title of the plot window. Defaults to None.
        :param type: str, optional
        :param name: patch - If True then coloring will be read in from "root directory"/data/tmp/{pdb_id}.pth file. Default: False.
        :param type: bool, optional
        :param name: box, optional - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - save image of plot to specified filename. Will appear in data/img/ directory. Default: data/img/{self._out_path}_surface.
        :param type: string, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot.
        """
        # get appropriate mesh and coloring
        mesh, cas = self._sh.return_mesh_and_color(self._msms, feature=feature, patch=patch)
        
        # plot
        pl = pv.Plotter(notebook=self._notebook)
        pl.background_color = 'grey'
        pl.enable_3_lights()
        pl.camera_position = 'xy'
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
        pl.show(screenshot=outname, title=title)

    def plot_hydrophob(self, box=None, res=None, outname=None):
        """
        Plot the hydrophobic features of a protein.

        :param name: box, optional - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - save image of plot to specified filename. Will appear in data/img/ directory. Default: data/img/{self._out_path}_surface.
        :param type: string, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot.
        """

        # save a screenshot
        if not outname:
            new_name = self._out_path.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_hydrophob.png'
        self.plot("hydrophob", "Hydrophob", box=box, res=res, outname=outname)

    def plot_shape(self, box=None, res=None, outname=None):
        """
        Plot the shape features of a protein.

        :param name: box, optional - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - save image of plot to specified filename. Will appear in data/img/ directory. Default: data/img/{self._out_path}_surface.
        :param type: string, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot.
        """
        
        # save a screenshot
        if not outname:
            new_name = self._out_path.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_shape.png'
        self.plot("shape", "Shape", box=box, res=res, outname=outname)

    def plot_charge(self, box=None, res=None, outname=None):
        """
        Plot the charge features of a protein.

        :param name: box, optional - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - save image of plot to specified filename. Will appear in data/img/ directory. Default: data/img/{self._out_path}_surface.
        :param type: string, optional
        
        :return: Pyvista.Plotter window - Window with interactive plot.
        """
    
        # save a screenshot
        if not outname:
            new_name = self._out_path.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = self._base_path + 'data/img/' + new_name + '_charge.png'
        self.plot("charge", "Charge", box=box, res=res, outname=outname)