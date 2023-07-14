import pyfada
import numpy as np
import time

#-----------------------------------------------------------------#
def testsparameters(pname, parameters, methods, getMgSolver, title):
  import matplotlib.pyplot as plt
  times, iters, Ns = {}, {}, []
  for method in methods:
    times[method] = []
    iters[method] = []
  for parameter in parameters:
    for method in methods:
      umg, solver = getMgSolver(parameter, method)
      print("parameter", parameter, "n_fine", umg.n_fine())
      t0 = time.time()
      iter = solver.testsolve(print=False)
      t1 = time.time()
      iters[method].append(iter)
      times[method].append(t1-t0)
    Ns.append(umg.n_fine())
  if pname=="N":
      plabel = r'$\log_{10}(N)$'
      pplot = np.log10(Ns)
  else:
      plabel = pname
      pplot = np.array(parameters)
  fig = plt.figure()
  ax = fig.add_subplot(211)
  ax.set_xlabel(plabel)
  ax.set_ylabel(r'it')
  ax.set_title(f"iter")
  for method in methods:
      ax.plot(pplot, iters[method], '-x', label=f"{method}")
  ax.legend()
  ax = fig.add_subplot(212)
  ax.set_xlabel(plabel)
  ax.set_ylabel(r't')
  ax.set_title(f"time")
  for method in methods:
      ax.plot(pplot, np.log10(times[method]), '-x', label=f"{method}")
  ax.legend()
  fig.suptitle(f"{title} {umg.dim()}d")
  plt.show()

#-----------------------------------------------------------------#
def testsmoothers(n=np.array([3,3]), nlevelmax=10, femtype="Q1", matrixtype="Trapez"):
  def get(nlevelmax, method):
    nlevels = nlevelmax
    umg = pyfada.UniformMultiGrid(nlevelmax, nlevels, n)
    solver = pyfada.SolverLaplace(umg, femtype, matrixtype, method)
    return umg, solver
  smoothers = ['Jac', 'GS1', 'GS2', 'GS']
  title = f"smoother comapre fem={femtype} {matrixtype}"
  parameters = np.arange(4, nlevelmax)
  testsparameters(pname="N", parameters=parameters, methods=smoothers, getMgSolver=get, title=title)


#-----------------------------------------------------------------#
def testcoarsesolve(n=np.array([3,3,3]), nlevels=7, stenciltype="Trapez", smoothertype="GS"):
  def get(param, method):
    umg = pyfada.UniformMultiGrid(nlevels, param, n)
    solver = pyfada.SolverLaplace(umg=umg, parameters={"stenciltype":stenciltype, "matrixtype":method, "smoothertype":smoothertype})
    return umg, solver
  parameters = np.arange(nlevels-3, nlevels)
  title = f"coarse solver {stenciltype=}"
  testsparameters(pname="nlevels", parameters=parameters, methods=["matrix"], getMgSolver=get, title=title)


#=================================================================#
if __name__ == '__main__':
    testcoarsesolve()
    # testsmoothers()
