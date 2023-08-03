import pathlib, sys
SCRIPT_DIR = str(pathlib.Path().absolute())
sys.path.insert(0,SCRIPT_DIR)
import pyfada
import numpy as np
import time
import matplotlib.pyplot as plt
from scipy.stats import linregress

#=================================================================#
if __name__ == '__main__':
    
    
    nlevels = np.arange(2,7)
    factors = np.power(2,np.arange(-2,1, dtype=float))
    
    nl, nf = len(nlevels), len(factors)
    niter = np.empty(shape=(nl,nf), dtype=int)
    err_p = np.empty(shape=(nl,nf), dtype=float)
    err_v = np.empty(shape=(nl,nf), dtype=float)
    niter_mean_v = np.empty(shape=(nl,nf), dtype=float)
    niter_mean_p = np.empty(shape=(nl,nf), dtype=float)
    dts = np.empty(shape=(nl,nf), dtype=float)
    
    for inl, nlevel in enumerate(nlevels):
        print(f"{inl}({nl})")
        for i, factor in enumerate(factors):
            h = 1/(2**nlevel)
            dt = 2.*factor*h**2
            # dt = 10.0*factor*h
            # dt = 1.0*factor
            dts[inl, i] = dt
            params = {'dt':str(dt), 'nlevels':str(nlevel), 'mgtimer':'false'}
            solver = pyfada.SolverStokes(params)
            io = solver.chorin_stationary(print=False)
            if io.niter==-1:
                print(f"no convergence for {nlevel=} {factor=}")
            niter[inl, i] = io.niter
            err_p[inl, i] = io.err_p
            err_v[inl, i] = io.err_v
            niter_mean_v[inl, i] = io.niter_mean_v
            niter_mean_p[inl, i] = io.niter_mean_p
            

    fig, axs = plt.subplots(2, 2, sharex=True)
    axs[0,0].set_title("niter")
    axs[0,1].set_title("niter_mean_v")
    axs[1,0].set_title("err_p")
    axs[1,1].set_title("err_v")
    for i, factor in enumerate(factors):
        axs[0,0].plot(nlevels, niter[:,i], label=f"{factor}")
        axs[0,1].plot(nlevels, niter_mean_v[:,i], label=f"{factor}")
        axs[1,0].plot(nlevels, np.log(err_p[:,i]), label=f"{factor}")
        res = linregress(nlevels, np.log(err_p[:,i]))
        axs[1,0].plot(nlevels, res.intercept + res.slope * nlevels, 'b--', label=f"r={-res.slope:4.1f}", linewidth=1)
        axs[1,1].plot(nlevels, np.log(err_v[:,i]), label=f"{factor}")
        res = linregress(nlevels, np.log(err_v[:,i]))
        axs[1,1].plot(nlevels, res.intercept + res.slope * nlevels, 'b--', label=f"r={-res.slope:4.1}", linewidth=1)
    for i in range(2): 
        for j in range(2): 
            axs[i,j].legend()
            axs[i,j].grid()
    plt.show()
