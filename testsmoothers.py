import pyfada
import numpy as np
import time

#-----------------------------------------------------------------#
def testsmoothers(n, nlevelsmax=12, matrixtype="Q1"):
  import matplotlib.pyplot as plt
  smoothers = ['jac', 'gs1', 'gs2', 'gs']
  nlevels = np.arange(3, nlevelsmax)
  times, iters, Ns = {}, {}, []
  for smoother in smoothers:
    times[smoother] = []
    iters[smoother] = []
  for nlevel in nlevels:
    op = pyfada.Operator(nlevel, n, matrixtype)
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
  fig.suptitle(f"{matrixtype} {op.dim()}d")
  plt.show()



#=================================================================#
if __name__ == '__main__':
  testsmoothers(np.array([3,3]), nlevelsmax=9)
  testsmoothers(np.array([3,3]), nlevelsmax=9, matrixtype="Q1Trapez")
  testsmoothers(np.array([3,3,3]), nlevelsmax=7)
  testsmoothers(np.array([3,3,3]), nlevelsmax=7, matrixtype="Q1Trapez")
