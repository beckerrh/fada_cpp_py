import pathlib, sys, os
SCRIPT_DIR = str(pathlib.Path().absolute())
sys.path.insert(0,SCRIPT_DIR)
import pyfada
import numpy as np
import time
import matplotlib.pyplot as plt
from scipy.stats import linregress
import h5py


#-----------------------------------------------------------------#
def test():
    nlevels = np.arange(1,8)
    errs, ns = [], []
    method = "Q1_0"
    application = "Sinus_dir"
    application = "Sinus_per"
    for nlevel in nlevels:
        solver = pyfada.SolverLaplace(parameters={"method": method, "application": application, "n0":"4", "nlevels":str(nlevel), "mgtimer":str(0)})

#-----------------------------------------------------------------#
def test():
    application = "SinusBis_per"
    errs_v, errs_p, ns = [], [], []

    # n0s = 2**np.arange(2,8)
    # for n0 in n0s:
        # solver = pyfada.SolverStokes(parameters={"application": application, "n0":str(n0), "nlevels":str(1), "mgtimer":str(0)})
    nlevels = np.arange(1,12)
    for nlevel in nlevels:
        solver = pyfada.SolverStokes(parameters={"application": application, "n0":"4", "nlevels":str(nlevel), "mgtimer":str(0)})

        info = solver.solve_stationary(print=False)
        print(f"{info=}")
        errs_v.append(info.err_v)
        errs_p.append(info.err_p)
        ns.append(solver.n_gridpoints())
        solver.save_for_visu()
    # print(f"{ns=}")
    # print(f"{errs_v=}")
    # print(f"{errs_p=}")
    logns, logerrs_v, logerrs_p = np.log10(ns), np.log10(errs_v), np.log10(errs_p)
    res = linregress(logns, logerrs_v)
    plt.plot(logns, logerrs_v, '-x', label="err_v")
    plt.plot(logns, res.intercept + res.slope * logns, 'b--', label=f"r={-res.slope:4.1f}", linewidth=1)
    res = linregress(logns, logerrs_p)
    plt.plot(logns, logerrs_p, '-x', label="err_p")
    plt.plot(logns, res.intercept + res.slope * logns, 'b--', label=f"r={-res.slope:4.1f}", linewidth=1)
    plt.legend()
    plt.grid()
    plt.show()


#=================================================================#
if __name__ == '__main__':
    test()
