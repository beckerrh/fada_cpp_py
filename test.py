import pyfada
import numpy as np
import time

#-----------------------------------------------------------------#
def vista(n=[3,3,3], nlevelmax=6, plot=True):
  import pyvista
  op = pyfada.Operator(nlevelmax, n)
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
  if len(n)==2: plotter.view_xy()
  plotter.show()


#=================================================================#
if __name__ == '__main__':
  vista(n=[3,3,3], nlevelmax=6)
