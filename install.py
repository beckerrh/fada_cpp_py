import pathlib, shutil
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
    rd = input(f"{rundir} exists remove[yN]")
    if rd=='y':
        shutil.rmtree(rundir)
rundir.mkdir(parents=True, exist_ok=True)
