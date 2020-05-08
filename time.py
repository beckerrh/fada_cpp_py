import pyfada
import numpy as np
import time, os, psutil
#-----------------------------------------------------------------#
def test(n, nlevels=12):
  op = pyfada.Operator(nlevels, n)
  t0 = time.time()
  iter = op.testsolve()
  t1 = time.time()
  process = psutil.Process(os.getpid())
  mem = process.memory_info().rss//10**6
  print(f"No. Iterations {iter:3d}  8 (N = {op.nall():8d})")
  print(f"Total time: {t1-t0:6.2f} Memory: {mem//10**3:3d}.{mem%10**3:3d} GB")



#=================================================================#
if __name__ == '__main__':
  test(np.array([3,3]), nlevels=12)
  test(np.array([3,3,3]), nlevels=8)
