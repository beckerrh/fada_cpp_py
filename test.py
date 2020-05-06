import pyfada
import numpy as np
import time

#-----------------------------------------------------------------#
def vista(n=[3,3,3], nlevelmax=6):
  import pyvista
  op = pyfada.Operator(nlevelmax, n)
  iter = op.testsolve()
  plotter = pyvista.Plotter()
  plotter.set_background("gray")
  print(f"dims = {op.get_dimensions().flat}")
  grid = pyvista.UniformGrid(op.get_dimensions())
  u = op.get_solution()
  print(f"u: {u.mean()} {u.max()}")
  grid.point_arrays['u'] = u
  plotter.add_mesh(grid, show_edges=True, scalars='u', opacity=0.3)
  plotter.add_mesh_isovalue(grid, show_edges=True, scalars='u')
  plotter.view_isometric()
  plotter.show()

#-----------------------------------------------------------------#
def testsmoothers(n, nlevelsmax=12):
  import matplotlib.pyplot as plt
  smoothers = ['jac', 'gs1', 'gs2']
  nlevels = np.arange(3, nlevelsmax)
  times, iters, Ns = {}, {}, []
  for smoother in smoothers:
    times[smoother] = []
    iters[smoother] = []
  for nlevel in nlevels:
    op = pyfada.Operator(nlevel, n)
    Ns.append(op.nall())
    for smoother in smoothers:
      op.smoother = smoother;
#      op.optmem = 0;
      t0 = time.time()
      iter = op.testsolve(print=False)
      t1 = time.time()
      iters[smoother].append(iter)
      times[smoother].append(t1-t0)
  fig = plt.figure()
  ax = fig.add_subplot(211)
  ax.set_xlabel(r'$\log_{10}(N)$')
  ax.set_ylabel(r'it')
  ax.set_title(f"iter")
  for smoother in smoothers:
      ax.plot(np.log10(Ns), iters[smoother], '-x', label=f"{smoother}")
  ax.legend()
  ax = fig.add_subplot(212)
  ax.set_xlabel(r'$\log_{10}(N)$')
  ax.set_ylabel(r't')
  ax.set_title(f"time")
  for smoother in smoothers:
      ax.plot(np.log10(Ns), np.log10(times[smoother]), '-x', label=f"{smoother}")
  ax.legend()
  plt.show()



#=================================================================#
if __name__ == '__main__':
#  testsmoothers(np.array([3,3]), nlevelsmax=8)
#  testsmoothers(np.array([3,3,3]), nlevelsmax=8)
  vista()
