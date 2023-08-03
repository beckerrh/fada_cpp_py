import pyfada
import numpy as np
import time, sys
import pyvista

#-----------------------------------------------------------------#
def vista(umg, plot=True):
  solver = pyfada.SolverLaplace(umg, {"stenciltype":"Trapez", "matrix":"matrix", "smoother":"GS"})
  print("solver", solver)
  application="DirichletRhsOne"
  application="Linear"
  iter = solver.testsolve(application=application)
  plotter = pyvista.Plotter()
  plotter.set_background([0.9,0.9,0.9])
  print(f"dims = {umg.get_dimensions().flat}")
  grid = pyvista.UniformGrid(dimensions=umg.get_dimensions().flat)
  u = solver.get_solution()
  print(f"u: {u.mean()} {u.max()}")
  if not plot: return
  grid.point_data['u'] = u
  plotter.add_mesh(grid, show_edges=True, scalars='u', opacity=0.3)
  plotter.add_mesh_isovalue(grid, show_edges=True, scalars='u', widget_color=[0.1,0.1,0.1], line_width=3)
  plotter.view_isometric()
  if umg.dim()==2: plotter.view_xy()
  plotter.show()


#=================================================================#
if __name__ == '__main__':
  umg = pyfada.UniformMultiGrid(6, 3, [3,3,3])
  umg = pyfada.UniformMultiGrid(6, 3, [3,3])
  # print("umg", umg)
  vista(umg)
