#!/usr/bin/env python3

import pathlib, shutil, sys
from distutils.spawn import find_executable

neededs = ['cmake']
for needed in neededs:
    if find_executable(needed) is None:
        print("executable '{}' not found".format(needed))
        sys.exit(1)

# datadir =  pathlib.Path.home().joinpath( 'data_dir', datadir_def_name)
runname = pathlib.Path().absolute().name + "_run"
rundir = pathlib.Path().absolute().parent / runname
if rundir.is_dir():
    rd = input(f"{rundir} exists. Remove? [yN]")
    if rd=='y':
        shutil.rmtree(rundir)
    else:
        sys.exit(1)
rundir.mkdir()
shutil.copy(pathlib.Path().absolute()/"fada.py", rundir)
