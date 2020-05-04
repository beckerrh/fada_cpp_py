import pyfada
import numpy as np
import time

print(f"{dir(pyfada)}")

nlevels = 12

n = np.array([3,3])
op = pyfada.Operateur(nlevels, n)

print(f"{dir(op)}")

op.smoother = "jac";

t0 = time.time()
iter = op.testsolve()

print(f"smoother = {op.smoother} iter = {iter} t={time.time()-t0:6.2f}")

u = op.get_solution()
print(f"u = {u.shape}")
