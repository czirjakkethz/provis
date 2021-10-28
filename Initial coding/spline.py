import pyvista

import numpy as np

points = np.random.random((3, 3))
print(points)

spline = pyvista.Spline(points, 10)
print(splin)
print(spline.lines)