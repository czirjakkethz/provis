import pyvista as pv
import numpy as np
import trimesh


from provis.src.data_handler import DataHandler

class Surface:
    def __init__(self, name):
        self._name = name
        self._dh = DataHandler(name)
        atom_data = self._dh.get_atoms()
        self._atmsurf, col = self._dh.get_atom_mesh(atom_data, vw=1, probe=0.1)
        
        
    def load_forv(self, file_name, end, vorf):
        """
        Load surface information from face or vert file
        
        :param name: file_name - Name of input file
        :param type: str
        :param name: end - Type of input file
        :param type: str
        :param name: vorf - Vertex or face file. "v" for vertex, "f" for face
        :param type: str
        
        :return: list - list of data
        """
        
        outfile = open(file_name + end,"r")
        data = outfile.readlines()
        l3 = str.split(data[2])
        numlines = int(l3[0])
        numspheres = int(l3[1])
        density = float(l3[2])
        probe = float(l3[3])
        ret = [[] for x in range(numlines)]
        i = 0
        for line in data[3:]:
            line_split = str.split(line)
            k = 0
            for entry in line_split[:3]:
                if vorf == "v":
                    ret[i].append(float(entry))
                elif vorf == "f":
                    ret[i].append(int(entry)-1)
            i+=1

        outfile.close()
        return ret
        

    def load_fv(self, file_name):
        """
        Load surface information from both face and vert files
        
        :param name: file_name - Name of input files
        :param type: str
        
        :return: list - list of face data
        :return: list - list of vert data
        """
        
        face = self.load_forv(file_name, ".face", "f")
        vert = self.load_forv(file_name, ".vert", "v")
        return face, vert


    def plot_msms_surface(self, filename):
        """
        Plot surface from face and vert files
        
        :param name: file_name - Name of face and vertex input file(s) without extension. Need to have the same name
        :param type: str
        
        :return: void - plot
        """
        face, vertice = self.load_fv(filename)
        vertices = np.array(vertice)
        faces = np.hstack(face)
        pl = pv.Plotter(lighting=None)
        pl.background_color = 'grey'
        pl.enable_3_lights()
        tmesh = trimesh.Trimesh(vertice, faces=face, process=False)
        mesh = pv.wrap(tmesh)
        pl.add_mesh(mesh)
        pl.show()

    def plot_surface(self):
        """
        Plot surface natively, without binaries.
        
        :returns: plot
        """

        pl = pv.Plotter(lighting=None)
        pl.background_color = 'grey'
        pl.enable_3_lights()

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

        pl.add_mesh(shell, color="white", smooth_shading=True, style=style, show_edges=False)#, culling='back')
        # save a screenshot
        pl.show(screenshot='test.png')

