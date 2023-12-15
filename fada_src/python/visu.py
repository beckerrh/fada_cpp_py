import pathlib, sys, os
SCRIPT_DIR = str(pathlib.Path().absolute())
sys.path.insert(0,SCRIPT_DIR)
import numpy as np
import h5py
import pyfada
import tkinter as tk
import fadalib.datatk
import matplotlib.gridspec
import matplotlib.pyplot as plt


class visu_backend:
    def def_plotters(self):
        return [self.plot]
    def init_from_directory(self, datadir):
        self.datadir = datadir
    def get_data(self):
        filename_grid = str(self.datadir / "grid.hdf")
        filename_solution = self.datadir / "solution.hdf"
        ugrid = pyfada.UniformGrid()
        ugrid.loadhdf5(filename_grid)
        n, dx = ugrid.n().flatten(), ugrid.dx().flatten()
        f = h5py.File(filename_solution, 'r')
        keys = list(f.keys())
        keys.remove("n")
        # print(f"{keys=} {n=}")
        return n, dx, {key:np.array(f[key][:]).reshape(*n) for key in keys}
    def plot(self, fig):
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        n, dx, data = self.get_data()
        x = np.arange(n[0])*dx[0]
        y = np.arange(n[1])*dx[1]
        keys = data.keys()
        nplots = len(keys)
        # gs = matplotlib.gridspec.GridSpec(1, nplots, wspace=0.4, hspace=0.3)
        # axs = [fig.add_subplot(gs[0, j]) for j in range(nplots)]
        axs = fig.subplots(nrows=1, ncols=nplots)
        if nplots==1: axs = [axs]
        # print(f"{nplots=}")
        for i,k in enumerate(keys):
            axs[i].set_title(f"{k}")
            CS = axs[i].contourf(x, y, data[k], levels=12)
            CS2 = axs[i].contour(CS, levels=CS.levels[::2], colors='k')
            divider = make_axes_locatable(axs[i])
            cax = divider.append_axes('right', size='5%', pad=0.4)
            cbar = plt.colorbar(CS, cax=cax, orientation='vertical')
            cbar.add_lines(CS2)

    def plot_pyvista(self, fig):
        
        # pyvista
        pyvista = False
        if pyvista:
            import pyvista as pv
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
    dtk = fadalib.datatk.DataTk(root=root, app=visu_backend())
    dtk.open()
    tk.mainloop()