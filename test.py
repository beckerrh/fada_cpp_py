import pyfada
import numpy as np

print(f"{dir(pyfada)}")

nlevels = 5

n = np.array([3,3])
op = pyfada.Operateur(nlevels, n)

print(f"{dir(op)}")

op.smoother = "gs1";

iter = op.testsolve()

print("hallo", op.smoother, iter)
