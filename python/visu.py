import pathlib, sys, os
SCRIPT_DIR = str(pathlib.Path().absolute())
sys.path.insert(0,SCRIPT_DIR)
import numpy as np
import h5py
import pyvista as pv
import pyfada
import datatk


class visu_backend:
    def def_plotters(self):
        return [self.plot1, self.plot2]
    def init_from_directory(self, datadir):
        self.datadir = datadir
    def plot(self, fig):
        filename_grid = str(self.datadir / "grid.hdf")
        filename_solution = self.datadir / "solution.hdf"
        f = h5py.File(filename_solution, 'r')
        keys = list(f.keys())
        keys.remove("n")
        data = {key:np.array(f[key][:]).T for key in keys}
        ugrid = pyfada.UniformGrid()
        ugrid.loadhdf5(filename_grid)
        
        n, dx = ugrid.n(), ugrid.dx()
        print(f"{n=} {dx=}")
        
        
        # pyvista
        pyvista = False
        if pyvista:
            plotter = pv.Plotter()
            plotter.set_background([0.9,0.9,0.9])
            # print(f"dims = {ugrid.get_dimensions().flat}")
            dims = ugrid.get_dimensions().flatten()
            if len(dims)==2: dims = [dims[0], dims[1], 1]
            grid = pv.UniformGrid(dimensions=dims)
            # u = np.array(f['dataset'][:])
            # print(f"u: {u.mean()} {u.max()}")
            if point_data:
                grid.point_data['u'] = data['u']
            else:
                grid.cell_data['u'] = u.T
            plotter.add_mesh(grid, show_edges=False, scalars='u', opacity=0.5)
            plotter.add_mesh_isovalue(grid, show_edges=False, scalars='u', widget_color=[0.1,0.1,0.1], line_width=2)
            plotter.view_isometric()
            if ugrid.dim()==2: plotter.view_xy()
            plotter.show()


#=================================================================#
if __name__ == '__main__':
    root = tk.Tk()
    dtk = datatk.DataTk(root=root, app=visu_backend())
    dtk.open()
    tk.mainloop()