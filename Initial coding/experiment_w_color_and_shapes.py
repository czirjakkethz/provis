import pyvista as pv
import numpy as np

# get dataset for the glyphs: supertoroids in xy plane
# use N random kinds of toroids over a mesh with 27 points
N = 5
values = np.arange(N)  # values for scalars to look up glyphs by


# taken from:
# rng = np.random.default_rng()
# params = rng.uniform(0.5, 2, size=(N, 2))  # (n1, n2) parameters for the toroids
params = np.array([[1.56821334, 0.99649769],
                   [1.08247844, 1.83758874],
                   [1.49598881, 0.83495047],
                   [1.52442129, 0.89600688],
                   [1.92212387, 0.78096621]])

geoms = [pv.ParametricSuperToroid(n1=n1, n2=n2) for n1, n2 in params]

# get dataset where to put glyphs
x,y,z = np.mgrid[:3, :3, :3]
mesh = pv.StructuredGrid(x, y, z)

# add random scalars
# rng_int = rng.integers(0, N, size=x.size)
scaling = [4, 1, 2, 0, 4, 0, 1, 4, 3, 1, 1, 3, 3, 4, 3, 4, 4,
                    3, 3, 2, 2, 1, 1, 1, 2, 0, 3]
rng_int = np.array([4, 1, 2, 0, 4, 0, 1, 4, 3, 1, 1, 3, 3, 4, 3, 4, 4,
                    3, 3, 2, 2, 1, 1, 1, 2, 0, 3])
mesh.point_data['scalars'] = rng_int

# construct the glyphs on top of the mesh; don't scale by scalars now
glyphs = mesh.glyph(geom=geoms, indices=values, scale=scaling,  rng=(0, N-1))

# create plotter and add our glyphs with some nontrivial lighting
plotter = pv.Plotter()
plotter.add_mesh(glyphs, specular=1, specular_power=15,
                 smooth_shading=True, show_scalar_bar=False)

plotter.save_graphic("exp.pdf", title='PyVista Export', raster=True, painter=True)

plotter.show()
