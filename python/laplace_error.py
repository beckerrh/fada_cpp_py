import pathlib, sys, os
SCRIPT_DIR = str(pathlib.Path().absolute())
sys.path.insert(0,SCRIPT_DIR)
import pyfada
import numpy as np
import time
import matplotlib.pyplot as plt
from scipy.stats import linregress
import h5py
import matplotlib.pyplot as plt
from scipy.stats import linregress



#-----------------------------------------------------------------#
def test():
    nlevels = np.arange(1,8)
    errs, ns = [], []
    method = "Q1_0"
    for nlevel in nlevels:
        solver = pyfada.SolverLaplace(parameters={"method": method, "application": "Sinus_dir", "n0":"3", "nlevels":str(nlevel), "mgtimer":str(False)})
        info = solver.testsolve(print=False)
        errs.append(info.err)
        ns.append(solver.n_gridpoints())
        print(f"{info=}")
        solver.save_for_visu()
    res = linregress(np.log(ns), np.log(errs))
    plt.plot(np.log(ns), np.log(errs), '-x', label="err")
    plt.plot(np.log(ns), res.intercept + res.slope * np.log(ns), 'b--', label=f"r={-res.slope:4.1f}", linewidth=1)
    plt.legend()
    plt.grid()
    plt.show()


#=================================================================#
if __name__ == '__main__':
    test()
