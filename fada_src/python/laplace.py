import pathlib, sys, os
SCRIPT_DIR = str(pathlib.Path().absolute())
sys.path.insert(0,SCRIPT_DIR)
import pyfada
import numpy as np
import time, sys
import pyvista

#-----------------------------------------------------------------#
def run(plot=True):
    application="Sinus_dir"
    solver = pyfada.SolverLaplace(parameters={"stenciltype":"Trapez", "matrix":"matrix", "smoother":"GS", "application":application})
    print("solver", solver)
    info = solver.testsolve()
    print(f"{info=}")
    plotter = pyvista.Plotter()
    plotter.set_background([0.9,0.9,0.9])
    dims = solver.get_grid().get_dimensions().flat
    if len(dims)==2: dims = [dims[0], dims[1], 1]
    grid = pyvista.UniformGrid(dimensions=dims)
    u = solver.get_solution_nodes()
    print(f"u: {u.mean()} {u.max()}")
    if not plot: return
    grid.point_data['u'] = u
    plotter.add_mesh(grid, show_edges=True, scalars='u', opacity=0.3)
    plotter.add_mesh_isovalue(grid, show_edges=True, scalars='u', widget_color=[0.1,0.1,0.1], line_width=3)
    plotter.view_isometric()
    if dims[2]==1: plotter.view_xy()
    plotter.show()


#=================================================================#
if __name__ == '__main__':
    run()
