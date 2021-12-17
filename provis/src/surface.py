import pyvista as pv
import numpy as np
import trimesh

from provis.src.data_handler import DataHandler
from provis.src.surface_handler import SurfaceHandler

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
        self._name = name
        if dens:
            self._density = dens
        self._dh = DataHandler(name)
        self._sh = SurfaceHandler(name, dens=dens)
        atom_data = self._dh.get_atoms()
        self._atmsurf, col = self._dh.get_atom_mesh(atom_data, vw=1, probe=0.1)


    def plot_msms_surface(self, outname=0):
        """
        Plot surface from face and vert files
        
        :param name: outname - save image of plot to specified filename. Will appear in data/output/ directory. default: data/output/{self._name}_surface
        :param type: string
        
        :return: void - plot
        """
        filename = self._name + '_out_' + str(int(self._density))
        # get faces and vertices
        face = self._dh.load_forv(filename, ".face", "f")
        vert = self._dh.load_forv(filename, ".vert", "v")
        vertices = np.array(vert)
        faces = np.hstack(face)

        #plot
        pl = pv.Plotter(lighting=None)
        pl.background_color = 'grey'
        pl.enable_3_lights()
        tmesh = trimesh.Trimesh(vertices, faces=faces, process=False)
        mesh = pv.wrap(tmesh)
        pl.add_mesh(mesh)
        
        # save a screenshot
        if not outname:
            new_name = self._name.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = 'data/output/' + new_name + '_msms_surf.png'
        pl.show(screenshot=outname)

    def plot_surface(self, outname=0):
        """
        Plot surface natively, without binaries.
        
        :param name: outname - save image of plot to specified filename. Will appear in data/output/ directory. default: data/output/{self._name}_surface
        :param type: string
        
        :returns: plot
        """

        # adding the spheres (by atom type) one at a time
        j = 0
        style = 'surface'
        mesh_ = pv.wrap(self._atmsurf[0])
        for mesh in self._atmsurf[1:]:
            mesh_ = mesh_ + (mesh)
            
        # create one mesh out of many spheres
        vol = mesh_.delaunay_3d(alpha=1.4)
        # extract surface from new mesh
        shell = vol.extract_surface().reconstruct_surface(sample_spacing=1.2)

        pl = pv.Plotter(lighting=None)
        pl.background_color = 'grey'
        pl.enable_3_lights()

        pl.add_mesh(shell, color="white", smooth_shading=True, style=style, show_edges=False)#, culling='back')
        # save a screenshot
        if not outname:
            new_name = self._name.split('/')
            new_name = new_name[-1].split('.')[0]
            outname = 'data/output/' + new_name + '_surface.png'
        pl.show(screenshot=outname)

    def plot_hydrophob(self, outname="hydrophob", patch=0):
        
        mesh, cas = self._sh.mesh_color(feature="hydrophob")

        # plot surface with feature visualization
        pl = pv.Plotter()
        pl.add_mesh(mesh, scalars=cas, cmap='RdBu')
        pl.background_color = 'white'
        pl.camera_position = 'xy'
        pl.show(screenshot=outname + '.jpeg')
