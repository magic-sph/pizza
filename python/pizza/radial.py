# -*- coding: utf-8 -*-
import os, re
import matplotlib.pyplot as plt
import numpy as np
from .log import PizzaSetup
from .libpizza import fast_read, scanDir


class PizzaRadial(PizzaSetup):
    """
    This class can be used to read and display the time and 
    horizontally averaged files radial_profiles.TAG

    >>> rad = PizzaRadial() # display the content of radial_profiles.tag
    """

    def __init__(self, datadir='.', iplot=True, tag=None, all=False):
        """
        :param datadir: working directory
        :type datadir: str
        :param iplot: to plot the output, default is True
        :type iplot: bool
        :param tag: a specific tag, default is None
        :type tag: str
        :param all: if all=True, then all the radial profiles from the directory
                    are stacked together and averaged by their respective time
                    span
        :type all: bool
        """

        name = 'radial_profiles'

        if not all:
            if tag is not None:
                pattern = os.path.join(datadir, '%s.%s' % (name, tag))
                files = scanDir(pattern)

                # Either the log.tag directly exists and the setup is easy to obtain
                if os.path.exists(os.path.join(datadir, 'log.%s' % tag)):
                    PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                        nml='log.%s' % tag)
                # Or the tag is a bit more complicated and we need to find 
                # the corresponding log file
                else:
                    mask = re.compile(r'%s\.(.*)' % name)
                    if mask.match(files[-1]):
                        ending = mask.search(files[-1]).groups(0)[0]
                        if logFiles.__contains__('log.%s' % ending):
                            PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                                nml='log.%s' % ending)

                # Sum the files that correspond to the tag
                mask = re.compile(r'%s\.(.*)' % name)
                for k, file in enumerate(files):
                    tag = mask.search(file).groups(0)[0]
                    nml = PizzaSetup(nml='log.%s' % tag, datadir=datadir, quiet=True)
                    filename = file
                    if k == 0:
                        self.tstart = nml.start_time
                        data = fast_read(filename)*(nml.stop_time-nml.start_time)
                    else:
                        if os.path.exists(filename):
                            data += fast_read(filename)*(nml.stop_time-nml.start_time)
                self.tstop = nml.stop_time
                data /= (self.tstop-self.tstart)

            else:
                pattern = os.path.join(datadir, '%s.*'% name)
                files = scanDir(pattern)
                filename = files[-1]
                # Determine the setup
                mask = re.compile(r'%s\.(.*)' % name)
                ending = mask.search(files[-1]).groups(0)[0]
                if os.path.exists('log.%s' % ending):
                    try:
                        PizzaSetup.__init__(self, datadir=datadir, quiet=True, 
                                        nml='log.%s' % ending)
                    except AttributeError:
                        pass

                data = fast_read(filename, skiplines=0)
        else:
            self.nsteps = 0
            pattern = os.path.join(datadir, '%s.*' % name)
            files = scanDir(pattern)

            # Determine the setup
            mask = re.compile(r'%s\.(.*)' % name)
            for k, file in enumerate(files):
                tag = mask.search(file).groups(0)[0]
                nml = PizzaSetup(nml='log.%s' % tag, datadir=datadir, quiet=True)
                filename = file
                if k == 0:
                    self.tstart = nml.start_time
                    data = fast_read(filename)*(nml.stop_time-nml.start_time)
                else:
                    if os.path.exists(filename):
                        data += fast_read(filename)*(nml.stop_time-nml.start_time)
            self.tstop = nml.stop_time
            data /= (self.tstop-self.tstart)
            PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                nml='log.%s' % tag)

        self.radius = data[:, 0]
        self.us2_mean = data[:, 1]
        self.us2_std = data[:, 2]
        self.up2_mean = data[:, 3]
        self.up2_std = data[:, 4]
        self.enst_mean = data[:, 5]
        self.enst_std = data[:, 6]
        self.uphi_mean = data[:, 7]
        self.uphi_std = data[:, 8]
        self.temp_mean = data[:, 9]
        self.temp_std = data[:, 10]
        self.nushell_mean = data[:, 11]
        self.nushell_std = data[:, 12]

        if iplot:
            self.plot()

    def plot(self):
        """
        Display the result when ``iplot=True``
        """

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.fill_between(self.radius, self.temp_mean-self.temp_std, self.temp_mean+\
                        self.temp_std, alpha=0.1)
        ax.plot(self.radius, self.temp_mean)
        ax.set_xlabel('Radius')
        ax.set_ylabel('Temperature')
        ax.set_ylim(self.temp_mean.min(), self.temp_mean.max())
        ax.set_xlim(self.radius[-1], self.radius[0])
        fig.tight_layout()

        if abs(self.uphi_mean).max() > 0.:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.fill_between(self.radius, self.uphi_mean-self.uphi_std, self.uphi_mean+\
                            self.uphi_std, alpha=0.1)
            ax.plot(self.radius, self.uphi_mean)
            ax.set_xlabel('Radius')
            ax.set_ylabel('Uphi')
            ax.set_xlim(self.radius[-1], self.radius[0])
            fig.tight_layout()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.fill_between(self.radius, self.us2_mean-self.us2_std, self.us2_mean+\
                        self.us2_std, alpha=0.1)
        ax.plot(self.radius, self.us2_mean, label='us**2')
        ax.fill_between(self.radius, self.up2_mean-self.up2_std, self.up2_mean+\
                        self.up2_std, alpha=0.1)
        ax.plot(self.radius, self.up2_mean, label='up**2')
        ax.fill_between(self.radius, self.enst_mean-self.enst_std, self.enst_mean+\
                        self.enst_std, alpha=0.1)
        ax.plot(self.radius, self.enst_mean, label='omega**2')
        ax.legend(loc='best', frameon=False)
        ax.set_xlabel('Radius')
        ax.set_xlim(self.radius[-1], self.radius[0])
        ax.set_yscale('log')
        fig.tight_layout()
