#-*- coding: utf-8 -*-
import os, re
import matplotlib.pyplot as plt
import numpy as np
from .log import PizzaSetup
import glob
from .libpizza import fast_read, scanDir


class PizzaTs(PizzaSetup):
    """
    This python class is used to read and plot the different time series 
    written by the code: 

       * Kinetic energy: :ref:`e_kin.TAG <secEkinFile>`
       * Heat transfer: :ref:`heat.TAG <secHeatFile>`
       * Power budget: :ref:`power.TAG <secpowerFile>`

    Here are a couple of examples of how to use this function.

    >>> # plot the most recent e_kin.TAG file found in the directoy
    >>> PizzaTs(field='e_kin')
    >>>
    >>> # stack **all** the power.TAG file found in the directory
    >>> ts = PizzaTs(field='power', all=True)
    >>> print(ts.time, ts.buoPower) # print time and buoyancy power 
    >>>
    >>> # If you only want to read the file ``heat.N0m2z``
    >>> ts = PizzaTs(field='heat', tag='N0m2z', iplot=False)
    """

    def __init__(self, datadir='.', field='e_kin', iplot=True, all=False, tag=None):
        """
        :param datadir: working directory
        :type datadir: str
        :param field: the file you want to plot
        :type field: str
        :param iplot: when set to True, display the plots (default True)
        :type iplot: bool
        :param all: when set to True, the complete time series is reconstructed by
                    stacking all the corresponding files from the working directory
                    (default False)
        :type all: bool
        :param tag: read the time series that exactly corresponds to the specified tag
        :type tag: str
        """
        self.field = field
        pattern = os.path.join(datadir, 'log.*')
        logFiles = scanDir(pattern)

        if tag is not None:
            pattern = os.path.join(datadir, '%s.%s' % (self.field, tag))
            files = scanDir(pattern)

            # Either the log.tag directly exists and the setup is easy to obtain
            if os.path.exists(os.path.join(datadir, 'log.%s' % tag)):
                PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % tag)
            # Or the tag is a bit more complicated and we need to find 
            # the corresponding log file
            else:
                mask = re.compile(r'%s\.(.*)' % self.field)
                if mask.match(files[-1]):
                    ending = mask.search(files[-1]).groups(0)[0]
                    if logFiles.__contains__('log.%s' % ending):
                        PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                            nml='log.%s' % ending)

            # Concatenate the files that correspond to the tag
            for k,file in enumerate(files):
                filename = file
                datanew = fast_read(filename)
                if k == 0:
                    data = datanew.copy()
                    ncolRef = data.shape[1]
                else:
                    ncol = datanew.shape[1]
                    if ncol == ncolRef:
                        data = np.vstack((data, datanew[1:,:]))
                    else: # If the number of columns has changed
                        data = np.vstack((data, datanew[1:, 0:ncolRef]))

        # If no tag is specified, the most recent is plotted
        elif not all:
            if len(logFiles) != 0:
                PizzaSetup.__init__(self, quiet=True, nml=logFiles[-1])
                name = '%s.%s' % (self.field, self.tag)
                filename = os.path.join(datadir, name)
                data = fast_read(filename)
            else:
                mot = '%s.*' % (self.field)
                dat = [(os.stat(i).st_mtime, i) for i in glob.glob(mot)]
                dat.sort()
                filename = dat[-1][1]
                data = fast_read(filename)

        # If no tag is specified but all=True, all the directory is plotted
        else:
            if len(logFiles) != 0:
                PizzaSetup.__init__(self, quiet=True, nml=logFiles[-1])
            pattern = os.path.join(datadir, '%s.*' % (self.field))
            files = scanDir(pattern)
            for k, file in enumerate(files):
                filename = file
                datanew = fast_read(filename)
                if k == 0:
                    data = datanew.copy()
                    ncolRef = data.shape[1]
                else:
                    if datanew.shape[0] != 0: # In case the file is empty
                        ncol = datanew.shape[1]
                        if ncol == ncolRef:
                            data = np.vstack((data, datanew[1:,:]))
                        else: # If the number of columns has changed
                            data = np.vstack((data, datanew[1:, 0:ncolRef]))

        if self.field == 'e_kin':
            self.time = data[:, 0]
            self.us2 = data[:, 1]
            self.up2 = data[:, 2]
            self.up2_axi = data[:, 3]
            self.ekin = self.us2+self.up2
        elif self.field == 'heat':
            self.time = data[:, 0]
            self.topnuss = data[:, 1]
            self.botnuss = data[:, 2]
            self.toptemp = data[:, 3]
            self.bottemp = data[:, 4]
            self.beta = data[:, 5]
        elif self.field == 'reynolds':
            self.time = data[:, 0]
            self.rey = data[:, 1]
            self.rey_zon = data[:, 2]
            self.rey_fluct = data[:, 3]
        elif self.field in ('power'):
            self.time = data[:, 0]
            self.buoPower = data[:, 1]
            self.viscDiss = data[:, 2]
        if iplot:
            self.plot()

    def plot(self):
        """
        Plotting subroutines. Only called if 'iplot=True'
        """
        if self.field == 'e_kin':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.us2, ls='-', c='#30a2da',
                    label='us**2')
            ax.plot(self.time, self.up2, ls='-', c='#fc4f30',
                    label='up**2')
            ax.plot(self.time, self.up2_axi, ls='--', c='#fc4f30',
                    label='up_axi**2')
            ax.plot(self.time, self.ekin, ls='-', c='#31363B')
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Ekin')
        elif self.field == 'reynolds':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.rey, label='Re')
            ax.plot(self.time, self.rey_fluct, label='Re fluct')
            ax.plot(self.time, self.rey_zon, label='Re zon')
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Reynolds')
        elif self.field == 'heat':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.time, self.topnuss, label='Top Nusselt')
            ax.plot(self.time, self.botnuss, label='Bottom Nusselt')
            ax.legend(loc='lower right', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Nusselt number')
        elif self.field in ('power'):
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.semilogy(self.time, self.buoPower, label='Thermal buoyancy')
            ax.semilogy(self.time, self.viscDiss, label='Viscous diss.')
            ax.legend(loc='best', frameon=False)
            ax.set_xlabel('Time')
            ax.set_ylabel('Power')
