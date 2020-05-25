import pyfada
import numpy as np
import time, sys

#-----------------------------------------------------------------#
#def vista(nlevelmax, nlevels, n, plot=True):
def vista(umg, plot=True):
#  print("umg", umg)
  import pyvista
#  nlevels = nlevelmax
#  op = pyfada.Operator(nlevelmax, nlevels, n)
#  print("umg.bounds()",umg.bounds())
#  sys.exit(1)
#  print("umg.dx()",umg.dx())
#  n = umg.n()
#  print("n", n)
#  print("umg.n()",umg.n())
  op = pyfada.Operator(umg.n(), umg.bounds())
  iter = op.testsolve()
  plotter = pyvista.Plotter()
  plotter.set_background([0.9,0.9,0.9])
  print(f"dims = {op.get_dimensions().flat}")
  grid = pyvista.UniformGrid(op.get_dimensions())
  u = op.get_solution()
  print(f"u: {u.mean()} {u.max()}")
  if not plot: return
  grid.point_arrays['u'] = u
  plotter.add_mesh(grid, show_edges=True, scalars='u', opacity=0.3)
  plotter.add_mesh_isovalue(grid, show_edges=True, scalars='u', widget_color=[0.1,0.1,0.1], line_width=3)
  plotter.view_isometric()
  if umg.dim==2: plotter.view_xy()
  plotter.show()


#=================================================================#
if __name__ == '__main__':
  umg = pyfada.UniformMultiGrid(6, 3, [3,3,3])
  vista(umg)
#  vista(6, 3, [3,3,3])
  
