import pathlib, sys, os
SCRIPT_DIR = str(pathlib.Path().absolute())
sys.path.insert(0,SCRIPT_DIR)
import pyfada
import numpy as np
import time
import matplotlib.pyplot as plt
from scipy.stats import linregress
import h5py



def get_data(dirname):
    dsplit = dirname.split('_')
    niter_out, nt, dt, q, end_time = int(dsplit[-1]), int(dsplit[-2]), float(dsplit[-3]), float(dsplit[-4]), float(dsplit[-5])
    statusfile = pathlib.Path(dirname).joinpath('status')
    with statusfile.open() as f: 
        print(f.readline()) 
        nsamples, nt2, ngrid = f.readline().split()
    nsamples, ngrid = int(nsamples), int(ngrid)
    
    data_v, data_p = np.empty(shape=(nt,2,ngrid)), np.empty(shape=(nt,ngrid))
    err_p, err_v, dts = np.zeros(nt-1), np.zeros(nt-1), np.zeros(nt-1)
    for isample in range(nsamples):
        for iter_out in range(niter_out):
            for it in range(nt):
                filename = dirname + f"/solution_{isample:06d}_{iter_out:06d}_{it:06d}_v0.hdf"
                file = h5py.File(filename, 'r')
                data_v[it,0,:] = np.array(file['dataset'][:])
                filename = dirname + f"/solution_{isample:06d}_{iter_out:06d}_{it:06d}_v1.hdf"
                file = h5py.File(filename, 'r')
                data_v[it,1,:] = np.array(file['dataset'][:])
                filename = dirname + f"/solution_{isample:06d}_{iter_out:06d}_{it:06d}_p.hdf"
                file = h5py.File(filename, 'r')
                data_p[it,:] = np.array(file['dataset'][:])
            for it in range(nt-1):
                data_v[it,:,:] -= data_v[nt-1,:,:]
                data_p[it,:] -= data_p[nt-1,:]
                err_v[it] += np.sum(data_v[it,:,:]**2)
                err_p[it] += np.sum(data_p[it,:]**2)
    for it in range(nt-1):
        err_v[it] = np.sqrt(err_v[it])/nsamples/ngrid/niter_out
        err_p[it] = np.sqrt(err_p[it])/nsamples/ngrid/niter_out
        dts[it] = dt/np.power(2,it);

    return err_v, err_p, dts, ngrid, nsamples 
                       

def plot(dt, err_p, err_v, ngrid, nsamples):        
    fig, axs = plt.subplots(1, 2, sharex=True)
    res = linregress(np.log(dt), np.log(err_p))
    fig.suptitle(f"N={ngrid} M={nsamples}")
    axs[0].set_title("err_p")
    axs[0].plot(np.log(dt), np.log(err_p), '-x')
    axs[0].plot(np.log(dt), res.intercept + res.slope * np.log(dt), 'b--', label=f"r={-res.slope:4.1f}", linewidth=1)
    axs[0].set_xlabel("log(dt)")
    axs[0].set_ylabel("log(err_p)")
    res = linregress(np.log(dt), np.log(err_v))
    axs[1].set_title("err_v")
    axs[1].plot(np.log(dt), np.log(err_v), '-x')
    axs[1].plot(np.log(dt), res.intercept + res.slope * np.log(dt), 'b--', label=f"r={-res.slope:4.1f}", linewidth=1)
    axs[1].set_ylabel("log(err_v)")    
    for j in range(2): 
        axs[j].legend()
        axs[j].grid()
    plt.show()
    

#=================================================================#
if __name__ == '__main__':

    err_v, err_p, dts, ngrid, nsamples  = get_data("datadir_Sinus_1_0.1_0.0625_8_2")
    print(f"{err_v=}")
    print(f"{err_p=}")
    plot(dts, err_p, err_v, ngrid, nsamples)
    # _application = "Sinus";
    # if(_datadir=="")
    # {
    #     _datadir = "datadir_" +_application;
    # }
    # std::stringstream ss;
    # ss << "," << _end_time << "," << _q << "," << _dt << "," << _nt << "," << _niter_out << "," ;
    # _datadir += ss.str();
    # if(not _directoryExists(_datadir))
    # {
    #     std::string command = "mkdir -p " + _datadir;
    #     system( command.c_str() );
    # }
        
    # os.environ["OMP_NUM_THREADS"] = "11"
    #
    # params = {'dt':0.0625, 'nt':8, 'end_time':1, 'niter_out':2, 'n0':13, 'nlevels':4, 'mgtimer':'false', 'nsamples':10, 'q':0.1}
    # params[]
    #
    #
    # # params = {'dt':1/16, 'nt':12, 'end_time':1, 'niter_out':2, 'n0':7, 'nlevels':2, 'mgtimer':'false', 'nsamples':1, 'q':0.01}
    # # params = {'dt':1/16, 'nt':5, 'end_time':1, 'niter_out':3, 'nlevels':6, 'mgtimer':'false', 'nsamples':2, 'q':0.01}
    #
    # params = {'dt':1/16, 'nt':8, 'end_time':1, 'niter_out':3, 'n0':4, 'nlevels':3, 'mgtimer':'false', 'nsamples':3, 'q':0.1}
    #
    # solver = pyfada.SolverStokes({k:str(v) for k,v in params.items()})
    # dt, err_p, err_v = solver.chorin_sde(print=True)
    # print(dt, err_p, err_v)
    # dt, err_p, err_v = dt.flat, err_p.flat, err_v.flat
    #
            

