import pyvista as pv
import numpy as np
import trimesh

from provis.src.data_handler import DataHandler
from provis.src.surface_handler import SurfaceHandler
from provis.utils.name_checker import check_name

class Surface:
    """
    The Surface class is used to visualize the surface information of the given molecule.
    """
    def __init__(self, name, dens=None):
        """
        Initialize Surface class with given filename. Creates internal data structures; a DataHandler to extract basic surface information and stores it in self._atmsurf (this is a list of Spheres for each atom roughly equating the Van-der-Waals radius).
        
        :param name: dens - sampling density used in msms binary. Also needed to load the face and vert files, as their (file)names include the density
        :param type: float
        :param name: name - name of file to be loaded
        :param type: str
        """
        self._path, self._out_path = check_name(name)
        if dens:
            self._density = dens
        self._sh = SurfaceHandler(name, dens=dens)

    def plot_msms_surface(self, outname=0):
        """
        Plot surface from face and vert files
        
        :param name: outname - save image of plot to specified filename. Will appear in data/img/ directory. default: data/img/{self._out_path}_surface
        :param type: string
        
        :return: void - plot
        """
        filename = self._out_path + '_out_' + str(int(self._density * 10))
        # get faces and vertices
        face = self._sh.load_forv(filename, ".face", "f")
        vert = self._sh.load_forv(filename, ".vert", "v")
        vertices = np.array(vert)
        faces = np.hstack(face)
        
        #plot
        pl = pv.Plotter(lighting=None)
        pl.background_color = 'grey'
        pl.enable_3_lights()
        faces.resize(int(len(faces)/3), 3)
        mesh = trimesh.Trimesh(vertices, faces=faces, process=False)
        # mesh = pv.wrap(tmesh)
        mesh = pv.wrap(mesh)
        pl.add_mesh(mesh)
        
        # save a screenshot
        if not outname:
            new_name = self._out_path.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = 'data/img/' + new_name + '_msms_surf.png'
        pl.show(screenshot=outname)

    def plot_surface(self, outname=0):
        """
        Plot surface natively, without binaries.
        
        :param name: outname - save image of plot to specified filename. Will appear in data/img/ directory. default: data/img/{self._out_path}_surface
        :param type: string
        
        :returns: plot
        """

        shell = self._sh.native_mesh()

        pl = pv.Plotter(lighting=None)
        pl.background_color = 'grey'
        pl.enable_3_lights()

        style = 'surface'
        pl.add_mesh(shell, color="white", smooth_shading=True, style=style, show_edges=False)#, culling='back')
        # save a screenshot
        if not outname:
            new_name = self._out_path.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = 'data/img/' + new_name + '_surface.png'
        pl.show(screenshot=outname)

    def plot_hydrophob(self, outname="hydrophob", patch=0):
        
        mesh, cas = self._sh.return_mesh_and_color(feature="hydrophob")

        # plot surface with feature visualization
        pl = pv.Plotter()
        pl.add_mesh(mesh, scalars=cas, cmap='RdBu')
        pl.background_color = 'white'
        pl.camera_position = 'xy'
        pl.show(screenshot=outname + '.jpeg')

    def plot_shape(self, outname="shape", patch=0):
        
        mesh, cas = self._sh.return_mesh_and_color(feature="shape")

        # plot surface with feature visualization
        pl = pv.Plotter()
        pl.add_mesh(mesh, scalars=cas, cmap='RdBu')
        pl.background_color = 'white'
        pl.camera_position = 'xy'
        pl.show(screenshot=outname + '.jpeg')

    def plot_charge(self, outname="charge", patch=0):
        
        mesh, cas = self._sh.return_mesh_and_color(feature="charge")

        # plot surface with feature visualization
        pl = pv.Plotter()
        pl.add_mesh(mesh, scalars=cas, cmap='RdBu')
        pl.background_color = 'white'
        pl.camera_position = 'xy'
        pl.show(screenshot=outname + '.jpeg')