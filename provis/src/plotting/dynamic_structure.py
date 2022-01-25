import pyvista
from provis.src.processing.file_converter import FileConverter
from provis.src.plotting.structure import Structure
from provis.src.plotting.surface import Surface
from provis.src.processing.name_checker import NameChecker
from provis.src.processing.data_handler import DataHandler
from provis.src.processing.surface_handler import SurfaceHandler

class DynamicStructure:
    """
    The Dynamic Structure class, similarly to the Protein class, encapsulates every other class and creates a user friendly way to plot your desired dynamic structure of a protein molecules.
    
    While the class is built similarly to the Protein class it does not use the Protein class itself. This is due to the fact that the Protein class is a rigid class made for a single molecule and
    """
  
    def __init__(self, pdb_name, base_path=None, density=3.0, plot_solvent=False, msms=False, notebook=False):
        """
        Initialize the class with the name of the pdb file and you are ready for plotting!
        
        Set the optional variables to get the most out of provis.

        :param name: pdb_name - Name of the pdb file with or without the extension. If the file is not stored in the data/pdb directory the full path has to be specified.
        :param type: str
        :param name: base_path - Path to the root directory of the required directory structure (see: https://pro-vis.readthedocs.io/en/latest/setup.html#). Only needs to be specified if it is not the current working directory. Default: None.
        :param type: str, optional
        :param name: density - Default: 3.0.
        :param type: float, optional
        :param name: plot_solvent - If True solvent atoms will also be plotted. Default: False.
        :param type: bool, optional
        :param name: msms - Set to True if you want to compute the surface information using the msms binary. If False surface will be computed natively. Default: False.
        :param type: bool, optional
        :param name: notebook - Set to True when using running in a Jupyter Notebook environment. Default: False.
        :param type: bool, optional
        """
        
            
#        Protein(pdb_name, base_path=None, density=3.0, plot_solvent=False, msms=False, notebook=False)
        if notebook:
            pyvista.set_jupyter_backend('pythreejs')
        self._notebook = notebook
        self._name_checker = NameChecker(pdb_name, base_path)
        self._path, self._out_path, self._base_path, self._mesh_path = self._name_checker.return_all()
        self.file_converter = FileConverter(self._name_checker, density=density)
        self._data_handler =  DataHandler(self._name_checker, fc=self.file_converter)
        self._surface_handler = SurfaceHandler( self._name_checker, fc=self.file_converter, dh=self._data_handler, density=density)
                
        self.structure = Structure(self._name_checker, dh=self._data_handler, plot_solvent=plot_solvent, notebook=notebook)
        self.surface = Surface(self._name_checker, sh=self._surface_handler, msms=msms, notebook=notebook)
        
    def plot_atoms(self, box=False, res=None, outname=None, camera=None, title="Atoms"):
        """
        Plot the dynamic atom cloud.

        :param name: box, optional - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_stick_point.png.
        :param type: string, optional
        :param name: camera - Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to 'xy'. Default: None.
        :param type: pyvista.Camera, optional
        :param name: title - Title of the plotting window. Default: "Atoms".
        :param type: str, optional
        """
        i = 0
        for model in self._data_handler._structure:
            self.structure.plot(atoms=1, box=box, vw=0, bonds=0, residues=0, res=res, outname=outname, title=title, camera=camera, model_id=i, dynamic=True)
            i+=1

    def plot_surface(self, feature=None, patch=False, title="Surface", box=None, res=None, outname=None, camera=None):
        """
        Plot the dynamic atom cloud.


        :param name: feature - Pass which feature (coloring) you want to plot. Options: hydrophob, shape, charge. Default: None (uniform coloring).
        :param type: str, optional
        :param name: patch - If True then coloring will be read in from "root directory"/data/tmp/{pdb_id}.pth file. Default: False.
        :param type: bool, optional
        :param name: box, optional - If True bounding box also visualized, default: 0.
        :param type: bool, optional
        :param name: res - Residues passed in 'res' will be plotted with a bounding box around them. Defaults to None.
        :param type: Residue, optional
        :param name: outname - Save image of plot to specified filename. Will appear in data/img directory. Defaults to data/img/{pdb_id}_stick_point.png.
        :param type: string, optional
        :param name: camera - Pass a Pyvista Camera https://docs.pyvista.org/api/core/camera.html to manually set the camera position. If nothing/None is passed then the camera position will be set to 'xy'. Default: None.
        :param type: pyvista.Camera, optional
        :param name: title - Title of the plotting window. Default: "Surface".
        :param type: str, optional
        """
        import numpy as np

        import pyvista as pv


        # x = np.arange(-10, 10, 0.25)
        # y = np.arange(-10, 10, 0.25)
        # x, y = np.meshgrid(x, y)
        # r = np.sqrt(x ** 2 + y ** 2)
        # z = np.sin(r)

        # # Create and structured surface
        # grid = pv.StructuredGrid(x, y, z)

        # # Create a plotter object and set the scalars to the Z height
        # plotter = pv.Plotter(notebook=False, off_screen=False)
        # plotter.add_mesh(grid, scalars=z.ravel(), smooth_shading=True)

        # # Open a gif
        # plotter.open_gif("wave.gif")

        # pts = grid.points.copy()

        # # Update Z and write a frame for each updated position
        # nframe = 15
        # for phase in np.linspace(0, 2 * np.pi, nframe + 1)[:nframe]:
        #     z = np.sin(r + phase)
        #     pts[:, -1] = z.ravel()
        #     print(pts)
        #     pts *= 0.9
        #     plotter.update_coordinates(pts, render=False)
        #     plotter.update_scalars(z.ravel(), render=False)

        #     # must update normals when smooth shading is enabled
        #     plotter.mesh.compute_normals(cell_normals=False, inplace=True)
        #     plotter.render()
        #     plotter.write_frame()

        # # Closes and finalizes movie
        # plotter.close()
        ##################################
        # Create and structured surface
        import time
        mesh, cas = self._surface_handler.return_mesh_and_color(False, feature=feature, patch=patch, model_id=0)
        if not cas:
            np.zeros(len(mesh.vertices))
            
        # Create a plotter object and set the scalars to the Z height
        plotter = pv.Plotter(notebook=self._notebook, off_screen=False)
        plotter.add_mesh(mesh, scalars=cas, cmap='RdBu', lighting=True, smooth_shading=True, show_edges=False)

        # Open a gif
        plotter.open_movie("animation.mp4", 1)


        # Update Z and write a frame for each updated position
        i = 0
        plotter.enable_eye_dome_lighting()
        plotter.render()
        plotter.write_frame()
        for model in self._data_handler._structure:
            mesh, cas = self._surface_handler.return_mesh_and_color(False, feature=feature, patch=patch, model_id=i)
            if not cas:
                np.zeros(len(mesh.vertices))
            
            i += 1
            # print(cas)
            
            plotter.clear()
            plotter.smooth_shading = True
            # plotter.enable_eye_dome_lighting()
            plotter.add_mesh(mesh, scalars=cas, cmap='RdBu', lighting=True, smooth_shading=True, show_edges=False)
            # plotter.update_coordinates(points, mesh=mesh, render=False)
            # plotter.update_scalars(cas, render=False)

            # must update normals when smooth shading is enabled
            plotter.mesh.compute_normals(cell_normals=False, inplace=True)
            plotter.render()
            plotter.write_frame()
            time.sleep(0.5)
            

        # Closes and finalizes movie
        plotter.close()
        ############################
        # import mayavi.mlab

        # mesh, cas = self._surface_handler.return_mesh_and_color(False, feature=feature, patch=patch, model_id=0, dynamic=True)

        # fig = mayavi.mlab.figure(bgcolor=(0, 0, 0), size=(640, 360))
        # print(mesh.vertices)
        # mayavi.mlab.points3d(mesh.vertices,
        #                     cas,          # Values used for Color
        #                     mode="point",
        #                     colormap='spectral', # 'bone', 'copper', 'gnuplot'
        #                     # color=(0, 1, 0),   # Used a fixed (r,g,b) instead
        #                     figure=fig,
        #                     )
        # mayavi.mlab.show()
        ##########################################
        # from threading import Thread
        # import time
        # import numpy as np
        # import pyvista as pv
        # import pyvistaqt as pvqt
        # from pyvista import examples

        # OFF_SCREEN = False
        # mesh_, cas = self._surface_handler.return_mesh_and_color(False, feature=feature, patch=patch, model_id=0, dynamic=True)
        # mesh = pv.wrap(mesh_)
        # p = pyvista.Plotter(off_screen=OFF_SCREEN)
        # func = lambda box: box # Does nothing
        # p.add_mesh(mesh)
        # mesh.plot()
        # #p.add_box_widget(callback=func)
        # p.close()

        # p = pyvista.Plotter(off_screen=OFF_SCREEN)
        # func = lambda box, widget: box # Does nothing
        # p.add_mesh(mesh)
        # p.add_box_widget(callback=func, pass_widget=True)
        # p.close()

        # p = pyvista.Plotter(off_screen=OFF_SCREEN)
        # p.add_mesh_clip_box(mesh)
        # p.close() 
        ###################################################################
        # globe = examples.load_globe()
        # globe.point_data['scalars'] = np.random.rand(globe.n_points)
        # globe.set_active_scalars('scalars')
        
        # plotter = pvqt.BackgroundPlotter()
        # plotter.add_mesh(globe, lighting=False, show_edges=True, texture=True, scalars='scalars')
        # plotter.view_isometric()

        # # shrink globe in the background
        # def shrink():
        #     for i in range(50):
        #         globe.points *= 0.95
        #         # Update scalars
        #         globe.point_data['scalars'] = np.random.rand(globe.n_points)
        #         time.sleep(0.5)

        # thread = Thread(target=shrink)
        # thread.start()
        
        
        
        # i = 0
        # import pyvista as pv
        # self._msms = False
        # self._shading = True
        # pl = pv.Plotter(notebook=False)
        # pl.background_color = 'grey'
        # pl.enable_3_lights()
        # for model in self._data_handler._structure:
        #     #self.surface.plot(feature=feature, title=title, patch=patch, box=box, res=res, outname=outname, camera=camera, model_id=i)
        #     mesh, cas = self._surface_handler.return_mesh_and_color(self._msms, feature=feature, patch=patch, model_id=i, dynamic=True)
        #     # plot

        #     print("Adding mesh...")
        #     pl.add_mesh(mesh, scalars=cas, cmap='RdBu', smooth_shading=self._shading, show_edges=False)
        #     print("Mesh added to plotter")
        #     # if specified add bounding box
        #     if box:
        #         pl.add_bounding_box(color='white', corner_factor=0.5, line_width=1)
            
        #     if res:
        #         res_list, chain_list, pad = res.get_res_info()
        #         for i, r in enumerate(res_list):
        #             chain = chain_list[i]
        #             residues = self._data_handler.get_structure().get_residues()
        #             residues_list = list(residues)
        #             res_name = residues_list[r + 1].get_resname()
        #             d = (self._data_handler._res_size_dict[res_name] + pad) * 2
        #             x, y, z = d,d,d
        #             pl.add_mesh(pv.Cube(center=self._data_handler.get_residue_info(r, chain,'com'), x_length=x, y_length=y, z_length=z), style='wireframe', show_edges=1, line_width=5, smooth_shading=self._shading, color='r')

            
        #     # save a screenshot
        #     if not outname:
        #         new_name = self._out_path.split('/')
        #         new_name = new_name[-1].split('.')[0]
        #         outname = self._base_path + 'data/img/' + new_name + '_surface.png'
                
        #     if camera: 
        #         pl.camera = camera
        #     else:
        #         pl.camera_position = 'xy'
        #     # pl.close()
        #     pl.show(screenshot=outname, title=title)
        #     pl.clear()
        #     i+=1
