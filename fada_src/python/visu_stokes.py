import pathlib, sys, os
SCRIPT_DIR = str(pathlib.Path().absolute())
sys.path.insert(0,SCRIPT_DIR)
import numpy as np
import h5py
import pyvista as pv
import pyfada

#-----------------------------------------------------------------#
def visu(datadir=None, filename_grid="grid.hdf", filename_p="solution_p.hdf", filenames_v=["solution_v0.hdf","solution_v1.hdf"]):

    if datadir is not None:
        filename_grid = str(datadir / "grid.hdf")
        filename_solution = datadir / "solution.hdf"
        
        
    ugrid = pyfada.UniformGrid()
    ugrid.loadhdf5(filename=filename_grid)
    dims = ugrid.get_dimensions().flatten()
    print(f"{dims=}")
    dim = len(dims)

    if datadir is not None:
        fp = h5py.File(filename_solution, 'r')
        # print(f"{list(fp.keys())=}")
        n = np.array(fp['n'][:])
        # print(f"{n=}")
        p = np.array(fp['p'][:]).reshape(dims).T.flatten()
        v0 = np.array(fp['v0'][:]).reshape(dims).T.flatten()
        v1 = np.array(fp['v1'][:]).reshape(dims).T.flatten()
        if "v2" in list(fp.keys()):        
            v2 = np.array(fp['v2'][:]).reshape(dims).T.flatten()
        else:
            v2 = np.zeros_like(v1)
    else:
        fp = h5py.File(filename_p, 'r')
        fv = [h5py.File(filename, 'r') for filename in filenames_v]
        p = np.array(fp['dataset'][:]).reshape(dims).T.flatten()    
        v0 = np.array(fv[0]['dataset'][:]).reshape(dims).T.flatten()
        v1 = np.array(fv[1]['dataset'][:]).reshape(dims).T.flatten()
        if len(fv) ==3:
            v2 = np.array(fv[2]['dataset'][:]).reshape(dims).T.flatten()
        else:
            v2 = np.zeros_like(v1)
    # print(f"{list(fp.keys())=}")
    # dset = fp['dataset']
    # print(f"{dset.shape}=")

    # g = h5py.File(filename_grid, 'r')
    # print(f"{list(g.keys())=}")
    # print(f"{g['n'].shape=}=")
    # print(f"{g['dx'].shape=}=")
    # print(f"{g['n'][:]=}=")
    # print(f"{g['dx'][:]=}=")


    print(f"p: {p.min()} {p.mean()} {p.max()} v0: {v0.min()} {v0.mean()} {v0.max()} v1: {v1.min()}  {v1.mean()} {v1.max()} v2: {v2.min()}{v2.mean()} {v2.max()}")

    v = np.vstack( [v0, v1, v2] ).T
    # v = np.hstack( [v0.T, v1.T, v2.T] )
    print(f"{v.shape=}")
    


    plotter = pv.Plotter(shape=(2, 2))
    plotter.set_background([0.9,0.9,0.9])
    
    if len(dims)==2: dims = [dims[0], dims[1], 1]
    grid = pv.UniformGrid(dimensions=dims)
    grid.point_data['p'] = p
    grid.point_data['v0'] = v0
    grid.point_data['v1'] = v1
    scale = 4/v.max()
    grid['vectors'] = scale*v
    grid.set_active_vectors('vectors')
    
    plotter.subplot(0,0)
    plotter.add_text("p", font_size=20, color='black')
    plotter.add_mesh(grid, show_edges=False, scalars='p', opacity=0.5, copy_mesh=True)
    plotter.add_mesh_isovalue(grid, show_edges=False, scalars='p', line_width=3)
    # plotter.view_isometric()
    if dim==2: plotter.view_xy()


    plotter.subplot(0,1)
    plotter.add_text("v", font_size=20, color='black')
    plotter.add_mesh(grid.arrows, lighting=False, color='black')
    # plotter.view_isometric()
    if dim==2: plotter.view_xy()
    
    plotter.subplot(1,0)
    plotter.add_text("v0", font_size=20, color='black')
    plotter.add_mesh(grid, show_edges=False, scalars='v0', opacity=0.5, copy_mesh=True)
    plotter.add_mesh_isovalue(grid, show_edges=False, scalars='v0', line_width=3)
    # plotter.view_isometric()
    if dim==2: plotter.view_xy()
    
    plotter.subplot(1,1)
    plotter.add_text("v1", font_size=20, color='black')
    plotter.add_mesh(grid, show_edges=False, scalars='v1', opacity=0.5, copy_mesh=True)
    plotter.add_mesh_isovalue(grid, show_edges=False, scalars='v1', line_width=3)
    # plotter.view_isometric()
    if dim==2: plotter.view_xy()
    
    plotter.show()


#=================================================================#
if __name__ == '__main__':
    visu(datadir=pathlib.Path("datadir_Sinus"))
