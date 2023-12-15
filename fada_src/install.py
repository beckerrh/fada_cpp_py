#!/usr/bin/env python3

import pathlib, shutil, sys, os, subprocess
from distutils.spawn import find_executable

neededs = ['cmake']
for needed in neededs:
    if find_executable(needed) is None:
        print("executable '{}' not found".format(needed))
        sys.exit(1)

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
os.chdir(rundir)
out = subprocess.run(["./fada.py", "cmake"], capture_output=True)
if out.returncode:
    print(f"{out.stdout=}")
    print(f"{out.stderr=}")
#subprocess.run(["./fada.py", "compile"], capture_output=True)