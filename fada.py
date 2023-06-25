#!/usr/bin/env python3

import sys, os, shutil, pathlib, subprocess
import argparse

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
class FadaRun:
    def __init__(self, sourcedir, installdir):
        self.CMAKE_BUILD_TYPES = ['Release','RelWithDebInfo','Debug','DebugFull','Profile']
        self.sourcedir, self.installdir = sourcedir, installdir
        self.builddir = installdir
        parser = argparse.ArgumentParser(
            # description='utility',
            usage='''fada <command> [<args>]
                    commands are:
                    compile     compile library and/or project
                    test        run test
                    ''')
        parser.add_argument('command', help='subcommand to run')
        parser.add_argument('-t', default = self.CMAKE_BUILD_TYPES[0], help='build type', choices=self.CMAKE_BUILD_TYPES)
        parser.add_argument('--cleanlib', default = False, action="store_true", help='clean library (build)')
        parser.add_argument('--verbose', default = False, action="store_true", help='compile blabbing')
        args = vars(parser.parse_args(sys.argv[1:]))
        if args['verbose']:
            print ('ArgumentsParser() args', args)
        if args['command']=='cmake':
            self.cmake(build_type=args['t'], verbose=args['verbose'], cleanlib=args['cleanlib'])
        elif args['command']=='compile':
            self.compile(build_type=args['t'], verbose=args['verbose'])
        else:
            raise ValueError(f"unknown subcommand {args['command']}")

    def cmake(self, build_type=None, verbose=False, cleanlib=False):
        if build_type is None: build_type=self.CMAKE_BUILD_TYPES[0]
        localbuilddir = self.builddir / build_type
        if verbose:
            print ('compile: sourcedir', self.sourcedir)
            print ('compile: installdir', self.installdir)
            print ('compile: builddir', self.builddir)
            print ('compile: localbuilddir', localbuilddir)
        if cleanlib:
            shutil.rmtree(localbuilddir, ignore_errors=True)
        try:
            os.makedirs(localbuilddir)
        except:
            pass
        startdir = os.getcwd()
        os.chdir(localbuilddir)
        cmakeoptions = " -DCMAKE_BUILD_TYPE="+build_type + " -DCMAKE_INSTALL_PREFIX="+str(self.installdir)
        command = "cmake " + str(self.sourcedir) + cmakeoptions
        returncode = runsubprocess(command)
    def compile(self, build_type=None, verbose=False):
        if build_type is None: build_type=self.CMAKE_BUILD_TYPES[0]
        localbuilddir = self.builddir / build_type
        if not localbuilddir.exists(): self.cmake(build_type=build_type, verbose=verbose)
        if verbose:
            print ('compile: sourcedir', self.sourcedir)
            print ('compile: installdir', self.installdir)
            print ('compile: builddir', self.builddir)
            print ('compile: localbuilddir', localbuilddir)
        startdir = os.getcwd()
        os.chdir(localbuilddir)
        command = "make -j4"
        returncode = runsubprocess(command)
        command = "make install"
        returncode = runsubprocess(command)
        os.chdir(startdir)


#------------------------------------------------------------------------

#ici
installdir = pathlib.Path().absolute()
sourcedir = installdir.parent / 'fada_cpp_py'

fr = FadaRun(sourcedir=sourcedir, installdir=installdir)
