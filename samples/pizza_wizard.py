#!/usr/bin/env python
from __future__ import print_function
import argparse
import os
import shutil
import sys
import subprocess as sp
import numpy as np
import unittest
import NonRotRa2e3.unitTest
import RotE1e3EkPump.unitTest
import RotInteg.unitTest
import TestRestart.unitTest
import TestChebMap.unitTest
import TimeStepChange.unitTest
import MultistepSchemes.unitTest
import DIRKSchemes.unitTest
import GalerkinBases.unitTest
import InhomogeneousHeatFlux.unitTest
import FingConv.unitTest

__version__ = '1.0'


def getParser():
    """
    Get script option
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s '+__version__, 
                        help="Show program's version number and exit.")
    parser.add_argument('--level', action='store', dest='test_level', type=int,
                        default=-1, help='Test level')
    parser.add_argument('--use-debug-flags', action='store_true', 
                        dest='use_debug_flags', 
                        default=False, help='Use compilation debug flags')
    parser.add_argument('--use-mkl', action='store_true', dest='use_mkl', 
                        default=False, 
                        help='Use the MKL for FFTs and Lapack calls')
    parser.add_argument('--use-precond', action='store', dest='use_precond', 
                        type=bool, default=True, 
                        help='Use matrix preconditioning')
    parser.add_argument('--nranks', action='store', dest='nranks', type=int,
                        default=4, help='Specify the number of MPI ranks')
    parser.add_argument('--mpicmd', action='store', dest='mpicmd', type=str,
                        default='mpirun', help='Specify the mpi executable')

    return parser


def wizard():
    print('\n')
    print("pizza's auto-tests suite")
    print('\n')

    print("                         ___              ")  
    print("                      |  ~~--.            ")
    print("                      |%=@%%/             ")
    print("                      |o%%%/              ")
    print("                   __ |%%o/               ")
    print("             _,--~~ | |(_/ ._             ")
    print("          ,/'  m%%%%| |o/ /  `\.          ")
    print("         /' m%%o(_)%| |/ /o%%m `\         ")
    print("       /' %%@=%o%%%o|   /(_)o%%% `\       ")
    print("      /  %o%%%%%=@%%|  /%%o%%@=%%  \      ")
    print("     |  (_)%(_)%%o%%| /%%%=@(_)%%%  |     ")
    print("     | %%o%%%%o%%%(_|/%o%%o%%%%o%%% |     ")
    print("     | %%o%(_)%%%%%o%(_)%%%o%%o%o%% |     ")
    print("     |  (_)%%=@%(_)%o%o%%(_)%o(_)%  |     ")
    print("      \ ~%%o%%%%%o%o%=@%%o%%@%%o%~ /      ")
    print("       \. ~o%%(_)%%%o%(_)%%(_)o~ ,/       ")
    print("         \_ ~o%=@%(_)%o%%(_)%~ _/         ")
    print("           `\_~~o%%%o%%%%%~~_/'           ")
    print("              `--..____,,--'              ")  
    print('\n')


def cmake(args, startdir, execDir):
    """
    Run cmake
    """
    if not os.path.exists(execDir):
        os.mkdir(execDir)
    else:
        shutil.rmtree(execDir)
        os.mkdir(execDir)
    os.chdir(execDir)

    if args.use_debug_flags:
         build_type='-DCMAKE_BUILD_TYPE=Debug'
    else:
         build_type='-DCMAKE_BUILD_TYPE=Release'

    if args.use_precond:
        precond_opt = '-DUSE_PRECOND=yes'
    else:
        precond_opt = '-DUSE_PRECOND=no'

    if args.use_mkl:
        mkl_opt = '-DUSE_LAPACKLIB=MKL'
    else:
        mkl_opt = '-DUSE_LAPACKLIB=LAPACK'

    # Compilation
    cmd = 'cmake %s/.. %s %s %s' % (startdir, build_type,
                                    precond_opt, mkl_opt)
    print('  '+cmd)
    print('\n')
    sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))


def compile():
    """
    Compile the code
    """
    cmd = 'make -j'
    print('  '+cmd)
    print('\n')
    sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
            stderr=open(os.devnull, 'wb'))


def get_env(args):
    """
    Display the environment variables
    """
    if 'FC' in os.environ:
        fortran_comp = os.environ['FC']
    else:
        fortran_comp = 'FC is not defined. Default Fortran compiler will be used!'
    if 'CC' in os.environ:
        c_comp = os.environ['CC']
    else:
        c_comp = 'CC is not defined. Default C compiler will be used!'
        
    print('  FC        : %s' % fortran_comp)
    print('  CC        : %s' % c_comp)
    print('  nranks    : %i' % args.nranks)
    print('  mpi exec  : %s' % args.mpicmd)
    print('  MKL       : %r' % args.use_mkl)
    print('\n')


def get_exec_cmd(args, execDir):
    """
    Determine execution command
    """
    pizzaExec = '%s/pizza.exe' % execDir

    os.environ['I_MPI_PIN_PROCESSOR_LIST'] = 'allcores'

    execCmd = '%s -n %i %s' % (args.mpicmd, args.nranks, pizzaExec)

    return execCmd


def clean_exec_dir(execDir):
    """
    Remove build directory
    """
    shutil.rmtree(execDir)


def getSuite(startdir, cmd, precision, args):
    """
    Construct test suite
    """
    suite = unittest.TestSuite()

    if args.test_level in [-1, 0]:
        # Non-rotating annulus convection
        suite.addTest(NonRotRa2e3.unitTest.NonRotRa2e3('outputFileDiff',
                                                  '%s/NonRotRa2e3' \
                                                  % startdir, 
                                                  execCmd=cmd,
                                                  precision=precision))
        # QG case with Ekman pumping
        suite.addTest(RotE1e3EkPump.unitTest.RotE1e3EkPump('outputFileDiff',
                                                  '%s/RotE1e3EkPump' % startdir, 
                                                  execCmd=cmd,
                                                  precision=precision))
        # Fingering convection
        suite.addTest(FingConv.unitTest.FingConv('outputFileDiff',
                                                 '%s/FingConv' % startdir, 
                                                 execCmd=cmd,
                                                 precision=precision))
        # QG case with Cheb integration method
        suite.addTest(RotInteg.unitTest.RotInteg('outputFileDiff',
                                                 '%s/RotInteg' % startdir, 
                                                 execCmd=cmd,
                                                 precision=precision))
        # Inhomogeneous heat flux at the outer Boundary
        suite.addTest(InhomogeneousHeatFlux.unitTest.InhomogeneousHeatFlux( \
                                              'outputFileDiff',
                                              '%s/InhomogeneousHeatFlux' % startdir,
                                              execCmd=cmd,
                                              precision=precision))
        # Test restart from a checkpoint
        suite.addTest(TestRestart.unitTest.TestRestart('outputFileDiff',
                                                  '%s/TestRestart' % startdir, 
                                                  execCmd=cmd,
                                                  precision=precision))
        # Test restart from a checkpoint and remap to a cheb grid
        suite.addTest(TestChebMap.unitTest.TestChebMap('outputFileDiff',
                                                  '%s/TestChebMap' % startdir, 
                                                  execCmd=cmd,
                                                  precision=precision))
        # Test several multistep schemes
        suite.addTest(MultistepSchemes.unitTest.TestMultistepSchemes('outputFileDiff',
                                                  '%s/MultistepSchemes' % startdir, 
                                                  execCmd=cmd,
                                                  precision=precision))
        # Test several DIRK schemes
        suite.addTest(DIRKSchemes.unitTest.TestDIRKSchemes('outputFileDiff',
                                                  '%s/DIRKSchemes' % startdir, 
                                                  execCmd=cmd,
                                                  precision=precision))
        # Test a timestep change
        suite.addTest(TimeStepChange.unitTest.TimeStepChange('outputFileDiff',
                                                  '%s/TimeStepChange' % startdir, 
                                                  execCmd=cmd,
                                                  precision=precision))
        # Test several Galerkin Bases
        suite.addTest(GalerkinBases.unitTest.TestGalerkinBases('outputFileDiff',
                                                  '%s/GalerkinBases' % startdir, 
                                                  execCmd=cmd,
                                                  precision=precision))

    return suite


if __name__ == '__main__':
    precision = 1e-8 # relative tolerance between expected and actual result
    startdir = os.getcwd()
    execDir = '%s/tmp' % startdir # where pizza will be built

    parser = getParser()
    args = parser.parse_args()

    # Initialisation
    wizard()

    # Display environment variables
    print('0.    Environment     ')
    print('----------------------')
    get_env(args)

    # Run cmake
    print('1. cmake configuration')
    print('----------------------')
    cmake(args, startdir, execDir)

    # Compile the code
    print('2.    compilation     ')
    print('----------------------')
    compile()

    # Determine the execution command
    cmd = get_exec_cmd(args, execDir)

    # Run the auto-test suite
    print('3.   Auto-tests       ')
    print('----------------------')
    suite = getSuite(startdir, cmd, precision, args)
    runner = unittest.TextTestRunner(verbosity=0)
    ret = not runner.run(suite).wasSuccessful()

    # Clean build directory
    clean_exec_dir(execDir)

    sys.exit(ret)
