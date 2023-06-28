import pyfada
import numpy as np
import time, os, psutil
#-----------------------------------------------------------------#
def test(n, nlevelmax=12, stenciltype="Trapez", matrixtype="stencil", smoothertype="GS"):
  nlevels = nlevelmax
  umg = pyfada.UniformMultiGrid(nlevelmax, nlevels, n)
  solver = pyfada.SolverLaplace(umg, stenciltype, matrixtype, smoothertype)
  t0 = time.time()
  iter = solver.testsolve(problem="DirichletRhsOne")
  t1 = time.time()
  process = psutil.Process(os.getpid())
  mem = process.memory_info().rss//10**6
  print(f"No. Iterations {iter:3d}  8 (N = {umg.nall():8d} dim = {umg.dim():1d})")
  print(f"Total time: {t1-t0:6.2f} Memory: {mem//10**3:3d}.{mem%10**3:3d} GB")



#=================================================================#
if __name__ == '__main__':
  test(np.array([3,3]), nlevelmax=12)
  test(np.array([3,3,3]), nlevelmax=8)
