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
    n0s = 4*np.arange(1,8)
    for n0 in n0s:
        solver = pyfada.SolverStokes(parameters={"n0":str(n0), "nlevels":str(1), "mgtimer":str(False)})
        info = solver.solve_stationary(print=False)
        print(f"{info=}")
        solver.save_for_visu()


#=================================================================#
if __name__ == '__main__':
    test()
