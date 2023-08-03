import numpy as np
import h5py
import pyvista
import pyfada
import pathlib

#-----------------------------------------------------------------#
def visu(filename_grid="grid.hdf", filename_solution="solution.hdf", point_data=True):

  f = h5py.File(filename_solution, 'r')
  # print(f"{list(f.keys())=}")
  # dset = f['dataset']
  # print(f"{dset.shape}=")

  # g = h5py.File(filename_grid, 'r')
  # print(f"{list(g.keys())=}")
  # print(f"{g['n'].shape=}=")
  # print(f"{g['dx'].shape=}=")
  # print(f"{g['n'][:]=}=")


  ugrid = pyfada.UniformGrid()
  ugrid.loadhdf5(filename_grid)

  plotter = pyvista.Plotter()
  plotter.set_background([0.9,0.9,0.9])
  # print(f"dims = {ugrid.get_dimensions().flat}")
  grid = pyvista.UniformGrid(dimensions=ugrid.get_dimensions().flat)
  u = np.array(f['dataset'][:])
  print(f"u: {u.mean()} {u.max()}")
  if point_data:
      grid.point_data['u'] = u.T
  else:
      grid.cell_data['u'] = u.T
  plotter.add_mesh(grid, show_edges=False, scalars='u', opacity=0.5)
  plotter.add_mesh_isovalue(grid, show_edges=False, scalars='u', widget_color=[0.1,0.1,0.1], line_width=2)
  plotter.view_isometric()
  if ugrid.dim()==2: plotter.view_xy()
  plotter.show()


#=================================================================#
if __name__ == '__main__':
  visu()
