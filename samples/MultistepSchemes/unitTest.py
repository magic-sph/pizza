from __future__ import print_function
import unittest
import numpy as np
import glob
import os
import time
import shutil
import subprocess as sp

def cleanDir(dir):
    for f in glob.glob('{}/*.test'.format(dir)):
        os.remove(f)
    if os.path.exists('{}/stdout.out'.format(dir)):
        os.remove('{}/stdout.out'.format(dir))
    for f in glob.glob('{}/*.pyc'.format(dir)):
        os.remove(f)
    if os.path.exists('{}/__pycache__'.format(dir)):
        shutil.rmtree('{}/__pycache__'.format(dir))

def readStack(file):
    f = open(file, 'r')
    out = np.array([])
    for line in f.readlines():
        cut = line.split()
        dat = np.asarray(cut, dtype=np.float64)
        out = np.append(out, dat)
    return out

class TestMultistepSchemes(unittest.TestCase):

    def __init__(self, testName, dir, execCmd='mpirun -n 8 ../tmp/pizza.exe', 
                 precision=1e-8):
        super(TestMultistepSchemes, self).__init__(testName)
        self.dir = dir
        self.precision = precision
        self.execCmd = execCmd
        self.startDir = os.getcwd()
        self.description = "Test several multistep schemes"

    def list2reason(self, exc_list):
        if exc_list and exc_list[-1][0] is self:
            return exc_list[-1][1]

    def setUp(self):
        # Cleaning when entering
        print('\nDirectory   :           {}'.format(self.dir))
        print('Description :           {}'.format(self.description))
        self.startTime = time.time()
        cleanDir(self.dir)
        for f in glob.glob('{}/*.start'.format(self.dir)):
            os.remove(f)
        for f in glob.glob('{}/*.first_continue'.format(self.dir)):
            os.remove(f)
        for f in glob.glob('{}/*.second_continue'.format(self.dir)):
            os.remove(f)
        for f in glob.glob('{}/*.third_continue'.format(self.dir)):
            os.remove(f)
        for f in glob.glob('{}/*.fourth_continue'.format(self.dir)):
            os.remove(f)
        for f in glob.glob('{}/*.final'.format(self.dir)):
            os.remove(f)

        os.chdir(self.dir)
        cmd = '{} {}/input.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))
        cmd = '{} {}/input_first_restart.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))
        cmd = '{} {}/input_second_restart.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))
        cmd = '{} {}/input_third_restart.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))
        cmd = '{} {}/input_fourth_restart.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))
        cmd = '{} {}/input_final.nml'.format(self.execCmd, self.dir)
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))
        cmd = 'cat e_kin_3D.start e_kin_3D.first_continue e_kin_3D.second_continue e_kin_3D.third_continue e_kin_3D.fourth_continue e_kin_3D.final > e_kin_3D.test'
        sp.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))

    def tearDown(self):
        # Cleaning when leaving
        os.chdir(self.startDir)
        cleanDir(self.dir)
        for f in glob.glob('{}/*.start'.format(self.dir)):
            os.remove(f)
        for f in glob.glob('{}/*.first_continue'.format(self.dir)):
            os.remove(f)
        for f in glob.glob('{}/*.second_continue'.format(self.dir)):
            os.remove(f)
        for f in glob.glob('{}/*.third_continue'.format(self.dir)):
            os.remove(f)
        for f in glob.glob('{}/*.fourth_continue'.format(self.dir)):
            os.remove(f)
        for f in glob.glob('{}/*.final'.format(self.dir)):
            os.remove(f)

        t = time.time()-self.startTime
        st = time.strftime("%M:%S", time.gmtime(t))
        print('Time used   :                            {}'.format(st))

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
            if failure:
                print('Validating results..                     FAIL!')
                print('\n')
                print(result.failures[-1][-1])

    def outputFileDiff(self):
        datRef = readStack('{}/reference.out'.format(self.dir))
        datTmp = readStack('{}/e_kin_3D.test'.format(self.dir))
        np.testing.assert_allclose(datRef, datTmp, rtol=self.precision,
                                   atol=1e-20)
