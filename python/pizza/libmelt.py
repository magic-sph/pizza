# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import copy
from .log import PizzaSetup
from .libpizza import scanDir
from .balance import avg_std
from scipy.interpolate import interp1d
import os
import re


class PizzaMelt(PizzaSetup):
    """
    This module is used to load and display the output files that
    contain the time evolution of the melt radius rmelt.TAG

    >>> m = PizzaMelt(all=True)
    >>> m.plot()
    """

    def __init__(self, datadir='.', tag=None, endian='l', iplot=False,
                 all=False):
        """
        :param datadir: working directory
        :type datadir: str
        :param tag: a specific trailing tag, default is None
        :type tag: str
        :param endian: endianness of the file ('B' or 'l')
        :type endian: str
        :param iplot: boolean to toggle the display (False by default)
        :type iplot: bool
        :param all: when set to True, the complete time series is reconstructed
                    by stacking all the corresponding files from the working
                    directory (False by default)
        :type all: bool
        """
        pattern = os.path.join(datadir, 'log.*')
        logFiles = scanDir(pattern)

        if tag is not None:
            pattern = os.path.join(datadir, 'rmelt.{}'.format(tag))
            files = scanDir(pattern)

            #  Either the log.tag directly exists and the setup is easy
            #  to obtain
            if os.path.exists(os.path.join(datadir, 'log.{}'.format(tag))):
                PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.{}'.format(tag))
            #  Or the tag is a bit more complicated and we need to find
            #  the corresponding log file
            else:
                mask = re.compile(r'{}/rmelt\.(.*)'.format(datadir))
                if mask.match(files[-1]):
                    ending = mask.search(files[-1]).groups(0)[0]
                    pattern = os.path.join(datadir, 'log.{}'.format(ending))
                    if os.path.exists(pattern):
                        PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                            nml='log.{}'.format(ending))

            for k, file in enumerate(files):
                print('reading {}'.format(file))
                if k == 0:
                    self.phi, self.time, self.rmelt, self.dt_rmelt = \
                        self._read(file, endian)
                else:
                    phi, time, rmelt, dt_rmelt = self._read(file, endian)
                    self.add(phi, time, rmelt, dt_rmelt)

        # If no tag is specified, the most recent is plotted
        elif not all:
            if len(logFiles) != 0:
                PizzaSetup.__init__(self, quiet=True, nml=logFiles[-1])
                name = 'rmelt.{}'.format(self.tag)
                filename = os.path.join(datadir, name)
                print('reading {}'.format(filename))
                self.phi, self.time, self.rmelt, self.dt_rmelt = \
                    self._read(filename, endian)
            else:
                dat = scanDir('rmelt.*')
                filename = dat[-1][1]
                print('reading {}'.format(filename))
                self.phi, self.time, self.time, self.rmelt = \
                    self._read(filename, endian)

        # If no tag is specified but all=True, all the directory is plotted
        else:
            if len(logFiles) != 0:
                PizzaSetup.__init__(self, quiet=True, nml=logFiles[-1])
            pattern = os.path.join(datadir, 'rmelt.*')
            files = scanDir(pattern)
            for k, file in enumerate(files):
                print('reading {}'.format(file))
                if k == 0:
                    self.phi, self.time, self.rmelt, self.dt_rmelt = \
                        self._read(file, endian)
                else:
                    phi, time, rmelt, dt_rmelt = self._read(file, endian)
                    self.add(phi, time, rmelt,  dt_rmelt)

        if iplot:
            self.plot()

    def __add__(self, new):
        """
        Built-in functions to sum two vphi force balances
        """
        out = copy.deepcopy(new)
        if self.n_phi_max != new.n_phi_max:
            ip = interp1d(self.phi, self.rmelt, axis=1, fill_value='extrapolate')
            self.rmelt = ip(self.rmelt, new.phi)
            ip = interp1d(self.phi, self.dt_rmelt, axis=1, fill_value='extrapolate')
            self.dt_rmelt = ip(self.dt_rmelt, new.phi)
            self.phi = new.phi

        if new.time[0] == self.time[-1]:
            out.time = np.concatenate((self.time, new.time[1:]), axis=0)
            out.rmelt = np.concatenate((self.rmelt, new.rmelt[1:, :]), axis=0)
            out.dt_rmelt = np.concatenate((self.dt_rmelt, new.dt_rmelt[1:, :]),
                                          axis=0)
        else:
            out.time = np.concatenate((self.time, new.time), axis=0)
            out.rmelt = np.concatenate((self.rmelt, new.rmelt), axis=0)
            out.dt_rmelt = np.concatenate((self.dt_rmelt, new.dt_rmelt), axis=0)

        return out

    def add(self, phi, time, rmelt, dt_rmelt):
        if self.n_phi_max != len(phi):
            ip = interp1d(self.phi, self.rmelt, axis=1, fill_value='extrapolate')
            self.rmelt = ip(self.rmelt, phi)
            ip = interp1d(self.phi, self.dt_rmelt, axis=1, fill_value='extrapolate')
            self.dt_rmelt = ip(self.dt_rmelt, phi)
            self.phi = phi

        if time[0] == self.time[-1]:
            self.time = np.concatenate((self.time, time[1:]), axis=0)
            self.rmelt = np.concatenate((self.rmelt, rmelt[1:, :]), axis=0)
            self.dt_rmelt = np.concatenate((self.dt_rmelt, dt_rmelt[1:, :]), axis=0)
        else:
            self.time = np.concatenate((self.time, time), axis=0)
            self.rmelt = np.concatenate((self.rmelt, rmelt), axis=0)
            self.dt_rmelt = np.concatenate((self.dt_rmelt, dt_rmelt), axis=0)

    def _read(self, filename, endian):
        """
        Read one rmelt.TAG file
        """

        with open(filename, 'rb') as file:
            version = np.fromfile(file, dtype=np.int32, count=1)[0]
            self.n_phi_max = np.fromfile(file, dtype=np.int32, count=1)[0]
            phi = np.fromfile(file, dtype=np.float64, count=self.n_phi_max)
            data = np.fromfile(file, dtype=np.float64, count=8)
            self.ra, self.ek, self.pr, self.radratio, self.raxi, \
                self.sc, self.stef, self.tmelt = data
            data = np.fromfile(file, dtype=np.float64)
            nsteps = len(data)//(2*self.n_phi_max+1)  # Modify if two fields

            data = data.reshape((nsteps, 2*self.n_phi_max+1))
            time = data[:, 0]
            rmelt = data[:, 1:self.n_phi_max+1]
            dt_rmelt = data[:, self.n_phi_max+1:2*self.n_phi_max+1]

        return phi, time, rmelt, dt_rmelt

    def plot(self, nstep=1):
        """
        Display when ``iplot = True``

        :param nstep: time interval (by default 1)
        :type nstep: int
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax1 = ax.twinx()
        rmelt_m, rmelt_std = avg_std(self.time[::nstep], self.rmelt[::nstep])
        ax.fill_between(self.phi, rmelt_m-rmelt_std, rmelt_m+rmelt_std,
                        alpha=0.1)
        ax.plot(self.phi, rmelt_m)
        dt_rmelt_m, dt_rmelt_std = avg_std(self.time[::nstep],
                                           self.dt_rmelt[::nstep])
        ax1.fill_between(self.phi, dt_rmelt_m-dt_rmelt_std, dt_rmelt_m+dt_rmelt_std,
                         alpha=0.1, color='C1')
        ax1.plot(self.phi, dt_rmelt_m, color='C1')
        ax.set_xlabel('Longitude')
        ax.set_ylabel('r melt')
        ax1.set_ylabel('dT(r melt)')
        ax.set_xlim(self.phi[0], self.phi[-1])
        fig.tight_layout()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.contourf(self.time[::nstep], self.phi, self.rmelt[::nstep, :].T, 65,
                    cmap=plt.get_cmap('terrain'))

        ax.set_xlabel('Time')
        ax.set_ylabel('Longitude')

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.contourf(self.time[::nstep], self.phi, self.dt_rmelt[::nstep, :].T, 65,
                    cmap=plt.get_cmap('plasma'))

        ax.set_xlabel('Time')
        ax.set_ylabel('Longitude')

        fig.tight_layout()


