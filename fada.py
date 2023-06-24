import sys, os, shutil, subprocess

#------------------------------------------------------------------------
def runsubprocess(command, options=None):
    if isinstance(command, list):
        commandstring = ' '.join( map(lambda x: str(x), command))
        commandlist = command
    else:
        commandstring = command
        commandlist = command.split()
    # print("command=%s\ncommandstring=%s\ncommandlist=%s" %(command, commandstring, commandlist))
    if options==None:
        return subprocess.check_call(commandstring, shell=True)
    elif options=="file":
        try:
            stdout_file = open('./simfem.stdout', 'a+')
            stderr_file = open('./simfem.stderr', 'a+')
        except:
            raise IOError('cannot open files')
        proc = subprocess.Popen(args=commandlist, stdout=stdout_file, stderr=stderr_file, shell=False)
        returncode = proc.wait()
        if returncode:
            raise RuntimeError("error in command '%s'" %commandstring)
        return returncode
    else:
        p = subprocess.Popen(commandlist, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout, stderr = p.communicate()
        if p.returncode != 0:
            raise RuntimeError("error in command '%s' (running in %s) errorcode=%d\n\nstdout=\x1b[6;30;42m %s \x1b[0m" %(commandstring, os.getcwd(), p.returncode,stdout.decode('ascii')))
        return p.returncode

#------------------------------------------------------------------------
def compile(sourcedir, installdir, build_type, verbose=False, cleanlib=False):
    builddir = os.path.join(startdirup, 'simfemsrc.compile')
    builddir = os.path.join(builddir, build_type)
    if verbose:
        print ('compile: sourcedir', sourcedir)
        print ('compile: installdir', installdir)
        print ('compile: builddir', builddir)
    if cleanlib:
        shutil.rmtree(builddir, ignore_errors=True)
    try:
        os.makedirs(builddir)
    except:
        pass
    startdir = os.getcwd()
    os.chdir(builddir)
    cmakeoptions = " -DCMAKE_BUILD_TYPE="+build_type + " -DCMAKE_INSTALL_PREFIX="+installdir
    command = "cmake " + sourcedir + cmakeoptions
    returncode = runsubprocess(command)
    command = "make -j4"
    returncode = runsubprocess(command)
    command = "make install"
    returncode = runsubprocess(command)
    os.chdir(startdir)


#------------------------------------------------------------------------
from argparse import ArgumentParser
parser = argparse.ArgumentParser(
    # description='utility',
    usage='''fada <command> [<args>]
            commands are:
            compile     compile library and/or project
            test        run test
            ''')
parser.add_argument('-t', default = self.CMAKE_BUILD_TYPES[0], help='build type', choices=self.CMAKE_BUILD_TYPES)
parser.add_argument('--cleanlib', default = False, action="store_true", help='clean library (build)')
parser.add_argument('--verbose', default = False, action="store_true", help='compile blabbing')
args = vars(parser.parse_args(sysargs))
build_type = args['t']
cleanlib = args['cleanlib']
verbose = args['verbose']
if verbose:
    print ('ArgumentsParser() args', args)
