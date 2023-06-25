import pyfada
import numpy as np
import time

#-----------------------------------------------------------------#
def testsparameters(nlevelmin, nlevelmax, parameters, paramname, getMgSolver, title):
  import matplotlib.pyplot as plt
  nlevelmaxs = np.arange(nlevelmin, nlevelmax)
  times, iters, Ns = {}, {}, []
  for param in parameters:
    times[param] = []
    iters[param] = []
  for nlevelmax in nlevelmaxs:
    for param in parameters:
      umg, solver = getMgSolver(nlevelmax, param)
      print("nlevelmax", nlevelmax, "nall", umg.nall())
      t0 = time.time()
      iter = solver.testsolve(print=False)
      t1 = time.time()
      iters[param].append(iter)
      times[param].append(t1-t0)
    Ns.append(umg.nall())
#  print("Ns", Ns)
  fig = plt.figure()
  ax = fig.add_subplot(211)
  ax.set_xlabel(r'$\log_{10}(N)$')
  ax.set_ylabel(r'it')
  ax.set_title(f"iter")
  for param in parameters:
      ax.plot(np.log10(Ns), iters[param], '-x', label=f"{param}")
  ax.legend()
  ax = fig.add_subplot(212)
  ax.set_xlabel(r'$\log_{10}(N)$')
  ax.set_ylabel(r't')
  ax.set_title(f"time")
  for param in parameters:
      ax.plot(np.log10(Ns), times[param], '-x', label=f"{param}")
#      ax.plot(np.log10(Ns), np.log10(times[param]), '-x', label=f"{param}")
  ax.legend()
  fig.suptitle(f"{title} {umg.dim()}d param={paramname}")
  plt.show()

#-----------------------------------------------------------------#
def testsmoothers(n, nlevelmax=12, femtype="Q1", matrixtype="Full"):
  def get(nlevelmax, param):
    nlevels = nlevelmax
    umg = pyfada.UniformMultiGrid(nlevelmax, nlevels, n)
    solver = pyfada.SolverLaplace(umg, femtype, matrixtype, param)
    return umg, solver
  smoothers = ['Jac', 'GS1', 'GS2', 'GS']
  pname = "smoother"
  title = f"fem={femtype} {matrixtype}"
  testsparameters(nlevelmin=3, nlevelmax=nlevelmax, parameters=smoothers, paramname=pname, getMgSolver=get, title=title)


#-----------------------------------------------------------------#
def testcoarsesolve(n, nlevelmin, nlevelmax, femtype="Q1", matrixtype="Full", smoothertype="Jac"):
  def get(nlevelmax, param):
    nlevels = param
    if param==-1: nlevels=nlevelmax
    umg = pyfada.UniformMultiGrid(nlevelmax, nlevels, n)
    solver = pyfada.SolverLaplace(umg, femtype, matrixtype, smoothertype)
    return umg, solver
  nlevels = np.arange(1,nlevelmin)
  nlevels = [1, -1]
  pname = "nlevels"
  title = f"fem={femtype} {matrixtype}"
  testsparameters(nlevelmin=nlevelmin, nlevelmax=nlevelmax, parameters=nlevels, paramname=pname, getMgSolver=get, title=title)


#=================================================================#
if __name__ == '__main__':
#  testcoarsesolve(np.array([3,3]), nlevelmin=6, nlevelmax=11)
#  testcoarsesolve(np.array([3,3,3]), nlevelmin=3, nlevelmax=6)
#  testsmoothers(np.array([3,3]), nlevelmax=9)
  testsmoothers(np.array([3,3]), nlevelmax=9, matrixtype="Trapez")
#  testsmoothers(np.array([3,3,3]), nlevelmax=7)
#  testsmoothers(np .array([3,3,3]), nlevelmax=7, matrixtype="Q1Trapez")
