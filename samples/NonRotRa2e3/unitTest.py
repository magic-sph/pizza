from __future__ import print_function
import unittest
import numpy as np
import glob
import os
import time
import shutil
import subprocess as sp

def cleanDir(dir):
    for f in glob.glob(f'{dir}/*.test'):
        os.remove(f)
    if os.path.exists(f'{dir}/stdout.out'):
        os.remove(f'{dir}/stdout.out')
    for f in glob.glob(f'{dir}/*.pyc'):
        os.remove(f)
    if os.path.exists(f'{dir}/__pycache__'):
        shutil.rmtree(f'{dir}/__pycache__')


def readData(file):
    return np.loadtxt(file)


class NonRotRa2e3(unittest.TestCase):

    def __init__(self, testName, dir, execCmd='mpirun -n 8 ../tmp/pizza.exe', 
                 precision=1e-8):
        super(NonRotRa2e3, self).__init__(testName)
        self.dir = dir
        self.precision = precision
        self.execCmd = execCmd
        self.startDir = os.getcwd()
        self.description = "Annulus non-rotating convection, Ra=2e3"

    def list2reason(self, exc_list):
        if exc_list and exc_list[-1][0] is self:
            return exc_list[-1][1]

    def setUp(self):
        # Cleaning when entering
        print(f'\nDirectory   :           {self.dir}')
        print(f'Description :           {self.description}')
        self.startTime = time.time()
        cleanDir(self.dir)
        os.chdir(self.dir)
        cmd = f'{self.execCmd} {self.dir}/input.nml'
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

    def tearDown(self):
        # Cleaning when leaving
        os.chdir(self.startDir)
        cleanDir(self.dir)

        t = time.time()-self.startTime
        st = time.strftime("%M:%S", time.gmtime(t))
        print(f'Time used   :                            {st}')

        if hasattr(self, '_outcome'): # python 3.4+
            if hasattr(self._outcome, 'errors'):  # python 3.4-3.10
                result = self.defaultTestResult()
                self._feedErrorsToResult(result, self._outcome.errors)
            else:  # python 3.11+
                result = self._outcome.result
        else:  # python 2.7-3.3
            result = getattr(self, '_outcomeForDoCleanups', 
                             self._resultForDoCleanups)

        error = self.list2reason(result.errors)
        failure = self.list2reason(result.failures)
        ok = not error and not failure

        if ok:
            print('Validating results..                     OK')
        else:
            if error:
                print('Validating results..                     ERROR!')
                print('\n')
                print(result.errors[-1][-1])
                print('\n')
            if failure:
                print('Validating results..                     FAIL!')
                print('\n')
                print(result.failures[-1][-1])
                print('\n')

    def outputFileDiff(self):
        # Kinetic energy
        datRef = readData(f'{self.dir}/reference.out')
        datTmp = readData(f'{self.dir}/e_kin.test')
        np.testing.assert_allclose(datRef, datTmp, rtol=self.precision, atol=1e-20)
