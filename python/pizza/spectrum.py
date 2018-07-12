# -*- coding: utf-8 -*-
import os, re
import matplotlib.pyplot as plt
import numpy as np
from .log import PizzaSetup
from .libpizza import scanDir, fast_read
from .npfile import *


class PizzaSpectrum(PizzaSetup):
    """
    This class can be used to read and display the spectra 'spec_#.TAG'
    or 'spec_avg.TAG'

    >>> # display the content of spec_1.TAG
    >>> # where TAG is the most recent file in the current directory
    >>> sp = PizzaSpectrum(ispec=1)
    """

    def __init__(self, datadir='.', iplot=True, ispec=None, ave=True, tag=None):
        """
        :param iplot: display the output plot when set to True (default is True)
        :type iplot: bool
        :param ispec: the number of the spectrum you want to plot
        :type ispec: int
        :param tag: file suffix (tag), if not specified the most recent one in
                    the current directory is chosen
        :type tag: str
        :param ave: plot a time-averaged spectrum when set to True
        :type ave: bool
        :param datadir: current working directory
        :type datadir: str
        """

        if ispec is not None:
            ave = False

        self.ave = ave

        if self.ave:
            self.name = 'spec_avg'
        else:
            self.name = 'spec_'

        if tag is not None:
            if ispec is not None:
                file = '%s%i.%s' % (self.name, ispec, tag)
                filename = os.path.join(datadir, file)
            else:
                pattern = os.path.join(datadir, '%s*%s' % (self.name, tag))
                files = scanDir(pattern)
                if len(files) != 0:
                    filename = files[-1]
                else:
                    print('No such tag... try again')
                    return

            if os.path.exists(os.path.join(datadir, 'log.%s' % tag)):
                PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % tag)
        else:
            if ispec is not None:
                pattern = os.path.join(datadir, '%s%i*' % (self.name, ispec))
                files = scanDir(pattern)
                filename = files[-1]
            else:
                pattern = os.path.join(datadir, '%s*' % self.name)
                files = scanDir(pattern)
                filename = files[-1]
            # Determine the setup
            mask = re.compile(r'.*\.(.*)')
            ending = mask.search(files[-1]).groups(0)[0]
            if os.path.exists(os.path.join(datadir, 'log.%s' % ending)):
                PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % ending)

        if not os.path.exists(filename):
            print('No such file')
            return

        data = fast_read(filename)

        self.index = data[:, 0]
        if not self.ave:
            self.us2 = data[:, 1]
            self.up2 = data[:, 2]
            self.enst = data[:, 3]
        else:
            self.us2_mean = data[:, 1]
            self.us2_std = data[:, 2]
            self.up2_mean = data[:, 3]
            self.up2_std = data[:, 4]
            self.enst_mean = data[:, 5]
            self.enst_std = data[:, 6]

        if iplot:
            self.plot()

    def plot(self):
        """
        Plotting function
        """
        if not self.ave:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            if abs(self.up2[0]) > 0.:
                ax.loglog(self.index[1:]+1, self.us2[1:], label='us**2')
                ax.loglog(self.index+1, self.up2, label='up**2')
                ax.loglog(self.index+1, self.enst, label='omega**2')
            else:
                ax.loglog(self.index[1:]+1, self.us2[1:], label='us**2')
                ax.loglog(self.index[1:]+1, self.up2[1:], label='up**2')
                ax.loglog(self.index[1:]+1, self.enst[1:], label='omega**2')
            ax.set_xlabel('m+1')
            ax.set_xlim(1, self.index[-1])
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()
        else:
            fig = plt.figure()
            ax = fig.add_subplot(111)

            if abs(self.up2_mean[0]) > 0.:
                ax.fill_between(self.index[1:]+1, self.us2_mean[1:]-self.us2_std[1:], \
                                self.us2_mean[1:]+self.us2_std[1:], alpha=0.1)
                ax.plot(self.index+1, self.us2_mean, label='us**2')

                ax.fill_between(self.index+1, self.up2_mean-self.up2_std, \
                                self.up2_mean+self.up2_std, alpha=0.1)
                ax.plot(self.index+1, self.up2_mean, label='up**2')

                ax.fill_between(self.index+1, self.enst_mean-self.enst_std, \
                                self.enst_mean+self.enst_std, alpha=0.1)
                ax.plot(self.index+1, self.enst_mean, label='omega**2')
            else:
                ax.fill_between(self.index[1:]+1, self.us2_mean[1:]-self.us2_std[1:], \
                                self.us2_mean[1:]+self.us2_std[1:], alpha=0.1)
                ax.plot(self.index[1:]+1, self.us2_mean[1:], label='us**2')

                ax.fill_between(self.index[1:]+1, self.up2_mean[1:]-self.up2_std[1:], \
                                self.up2_mean[1:]+self.up2_std[1:], alpha=0.1)
                ax.plot(self.index[1:]+1, self.up2_mean[1:], label='up**2')

                ax.fill_between(self.index[1:]+1, self.enst_mean[1:]-self.enst_std[1:], \
                                self.enst_mean[1:]+self.enst_std[1:], alpha=0.1)
                ax.plot(self.index[1:]+1, self.enst_mean[1:], label='omega**2')

            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.set_xlabel('m+1')
            ax.set_xlim(1, self.index[-1])
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()

