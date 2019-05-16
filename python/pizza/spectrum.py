# -*- coding: utf-8 -*-
import copy
import os
import re
import matplotlib.pyplot as plt
import numpy as np
from .log import PizzaSetup
from .libpizza import scanDir, fast_read


class PizzaSpectrum(PizzaSetup):
    """
    This class can be used to read and display the spectra 'spec_#.TAG'
    or 'spec_avg.TAG' and the vorticity balance spectra 'vort_terms_avg.TAG'

    >>> # display the content of spec_1.TAG
    >>> # where TAG is the most recent file in the current directory
    >>> sp = PizzaSpectrum(ispec=1)
    >>> # stack the content of spec_ave.test_a to spec_ave.test_c
    >>> sp = PizzaSpectrum(tag='test_[a-c]', iplot=False)
    """

    def __init__(self, field='spec', datadir='.', iplot=True, ispec=None,
                 ave=True, tag=None, all=False):
        """
        :param field: a string to decide which spectra one wants to load
                      and plot
        :type field: str
        :param iplot: display the output plot when set to True
                      (default is True)
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
        :param all: when set to True, the complete time series is reconstructed
                    by stacking all the corresponding files from the working
                    directory (False by default)
        :type all: bool
        """

        if ispec is not None:
            ave = False

        self.ave = ave

        if field in ('vort_terms', 'forces', 'forcebal', 'force_bal'):
            self.name = 'vort_terms_avg'
            self.ave = True
            ispec = None
        else:
            if self.ave:
                self.name = 'spec_avg'
            else:
                self.name = 'spec_'

        if not all:
            if tag is not None:
                if ispec is not None:
                    self.name += '%i' % ispec
                pattern = os.path.join(datadir, '%s.%s' % (self.name, tag))
                files = scanDir(pattern)
                # Either the log.tag directly exists and the setup is easy
                # to obtain
                if os.path.exists(os.path.join(datadir, 'log.%s' % tag)):
                    PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                        nml='log.%s' % tag)
                # Or the tag is a bit more complicated and we need to find
                # the corresponding log file
                else:
                    mask = re.compile(r'%s/%s\.(.*)' % (datadir, self.name))
                    if mask.match(files[-1]):
                        ending = mask.search(files[-1]).groups(0)[0]
                        pattern = os.path.join(datadir, 'log.%s' % ending)
                        if os.path.exists(pattern):
                            PizzaSetup.__init__(self, datadir=datadir,
                                                quiet=True,
                                                nml='log.%s' % ending)

                # Sum the files that correspond to the tag
                mask = re.compile(r'%s\.(.*)' % self.name)
                for k, file in enumerate(files):
                    print('reading %s' % file)
                    tag = mask.search(file).groups(0)[0]
                    nml = PizzaSetup(nml='log.%s' % tag, datadir=datadir,
                                     quiet=True)
                    filename = file
                    if k == 0:
                        self.tstart = nml.start_time
                        self.tstop = nml.stop_time  # will be replaced later
                        data = fast_read(filename)
                    else:
                        if os.path.exists(filename):
                            tmp = fast_read(filename)
                            data = self.add(data, tmp, nml.stop_time,
                                            nml.start_time)

            else:
                if ispec is not None:
                    pattern = os.path.join(datadir, '%s%i*' % (self.name,
                                                               ispec))
                else:
                    pattern = os.path.join(datadir, '%s*' % self.name)
                files = scanDir(pattern)
                filename = files[-1]
                print('reading %s' % filename)
                # Determine the setup
                if ispec is not None:
                    mask = re.compile(r'%s%i\.(.*)' % (self.name, ispec))

                else:
                    mask = re.compile(r'%s\.(.*)' % self.name)
                ending = mask.search(files[-1]).groups(0)[0]
                if os.path.exists('log.%s' % ending):
                    try:
                        PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                            nml='log.%s' % ending)
                    except AttributeError:
                        pass

                data = fast_read(filename, skiplines=0)

        else:  # if all is requested
            pattern = os.path.join(datadir, '%s.*' % self.name)
            files = scanDir(pattern)

            # Determine the setup
            mask = re.compile(r'%s\.(.*)' % self.name)
            for k, file in enumerate(files):
                print('reading %s' % file)
                tag = mask.search(file).groups(0)[0]
                nml = PizzaSetup(nml='log.%s' % tag, datadir=datadir,
                                 quiet=True)
                filename = file
                if k == 0:
                    self.tstart = nml.start_time
                    self.tstop = nml.stop_time  # will be replaced later
                    data = fast_read(filename)
                else:
                    if os.path.exists(filename):
                        tmp = fast_read(filename)
                        data = self.add(data, tmp, nml.stop_time,
                                        nml.start_time)
            PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                nml='log.%s' % tag)

            # Determine the setup
            mask = re.compile(r'.*\.(.*)')
            ending = mask.search(files[-1]).groups(0)[0]
            if os.path.exists(os.path.join(datadir, 'log.%s' % ending)):
                PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % ending)

        self.index = data[:, 0]

        if self.name == 'vort_terms_avg':
            self.buo_mean = data[:, 1]
            self.buo_std = data[:, 2]
            self.cor_mean = data[:, 3]
            self.cor_std = data[:, 4]
            self.adv_mean = data[:, 5]
            self.adv_std = data[:, 6]
            self.domdt_mean = data[:, 7]
            self.domdt_std = data[:, 8]
            self.visc_mean = data[:, 9]
            self.visc_std = data[:, 10]
            self.pump_mean = data[:, 11]
            self.pump_std = data[:, 12]
            self.thwind_mean = data[:, 13]
            self.thwind_std = data[:, 14]
            self.iner_mean = data[:, 15]
            self.iner_std = data[:, 16]
            self.cia_mean = data[:, 17]
            self.cia_std = data[:, 18]
        else:
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

    def add(self, data, tmp, stop_time, start_time):
        """
        Clean way to stack data
        """
        out = copy.deepcopy(data)
        out[:, 0] = tmp[:, 0]

        nm_new = len(tmp[:, 0])
        nm_old = len(data[:, 0])

        fac_old = self.tstop-self.tstart
        fac_new = stop_time-start_time
        self.tstop = stop_time
        fac_tot = self.tstop-self.tstart

        if nm_new == nm_old:  # Same grid before and after
            if self.name == 'vort_terms_avg':
                for j in [1, 3, 5, 7, 9, 11, 13, 15, 17]:
                    out[:, j] = (fac_old*data[:, j]+fac_new*tmp[:, j])/fac_tot
                for j in [2, 4, 6, 8, 10, 12, 14, 16]:
                    out[:, j] = np.sqrt((fac_old*data[:, j]**2 +
                                         fac_new*tmp[:, j]**2) / fac_tot)
            else:
                if not self.ave:
                    out[:, 1:] = 0.5*(data[:, 1:]+tmp[:, 1:])
                else:
                    for j in [1, 3, 5]:
                        out[:, j] = (fac_old*data[:, j]+fac_new*tmp[:, j]) / \
                                     fac_tot
                    for j in [2, 4, 6]:
                        out[:, j] = np.sqrt((fac_old*data[:, j]**2 +
                                             fac_new*tmp[:, j]**2) / fac_tot)
        else:
            print('Not implemented yet ...')

        return out

    def plot(self):
        """
        Plotting function
        """
        if self.name == 'vort_terms_avg':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            sd = self.buo_std/np.sqrt(self.buo_mean)/2.
            ax.fill_between(self.index, np.sqrt(self.buo_mean)-sd,
                            np.sqrt(self.buo_mean)+sd, alpha=0.1)
            ax.plot(self.index, np.sqrt(self.buo_mean), label='Buoyancy')

            sd = self.cor_std/np.sqrt(self.cor_mean)/2.
            ax.fill_between(self.index, np.sqrt(self.cor_mean)-sd,
                            np.sqrt(self.cor_mean)+sd, alpha=0.1)
            ax.plot(self.index, np.sqrt(self.cor_mean), label='Coriolis')

            sd = self.iner_std/np.sqrt(self.iner_mean)/2.
            ax.fill_between(self.index, np.sqrt(self.iner_mean)-sd,
                            np.sqrt(self.iner_mean)+sd, alpha=0.1)
            ax.plot(self.index, np.sqrt(self.iner_mean), label='Inertia')

            sd = self.visc_std/np.sqrt(self.visc_mean)/2.
            ax.fill_between(self.index, np.sqrt(self.visc_mean)-sd,
                            np.sqrt(self.visc_mean)+sd, alpha=0.1)
            ax.plot(self.index, np.sqrt(self.visc_mean), label='Viscosity')

            sd = self.pump_std/np.sqrt(self.pump_mean)/2.
            ax.fill_between(self.index, np.sqrt(self.pump_mean)-sd,
                            np.sqrt(self.pump_mean)+sd, alpha=0.1)
            ax.plot(self.index, np.sqrt(self.pump_mean), label='Ekman pumping')

            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.set_xlabel('m')
            ax.set_ylabel('Forces')
            ax.set_xlim(1, self.index[-1])
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()

        else:
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
                    ax.loglog(self.index[1:]+1, self.enst[1:],
                              label='omega**2')
                ax.set_xlabel('m+1')
                ax.set_xlim(1, self.index[-1])
                ax.legend(loc='best', frameon=False)
                fig.tight_layout()
            else:
                fig = plt.figure()
                ax = fig.add_subplot(111)

                if abs(self.up2_mean[0]) > 0.:
                    ax.fill_between(self.index[1:]+1,
                                    self.us2_mean[1:]-self.us2_std[1:],
                                    self.us2_mean[1:]+self.us2_std[1:],
                                    alpha=0.1)
                    ax.plot(self.index[1:]+1, self.us2_mean[1:], label='us**2')

                    ax.fill_between(self.index+1, self.up2_mean-self.up2_std,
                                    self.up2_mean+self.up2_std, alpha=0.1)
                    ax.plot(self.index+1, self.up2_mean, label='up**2')

                    ax.fill_between(self.index+1, self.enst_mean-self.enst_std,
                                    self.enst_mean+self.enst_std, alpha=0.1)
                    ax.plot(self.index+1, self.enst_mean, label='omega**2')
                else:
                    ax.fill_between(self.index[1:]+1,
                                    self.us2_mean[1:]-self.us2_std[1:],
                                    self.us2_mean[1:]+self.us2_std[1:],
                                    alpha=0.1)
                    ax.plot(self.index[1:]+1, self.us2_mean[1:], label='us**2')

                    ax.fill_between(self.index[1:]+1,
                                    self.up2_mean[1:]-self.up2_std[1:],
                                    self.up2_mean[1:]+self.up2_std[1:],
                                    alpha=0.1)
                    ax.plot(self.index[1:]+1, self.up2_mean[1:], label='up**2')

                    ax.fill_between(self.index[1:]+1,
                                    self.enst_mean[1:]-self.enst_std[1:],
                                    self.enst_mean[1:]+self.enst_std[1:],
                                    alpha=0.1)
                    ax.plot(self.index[1:]+1, self.enst_mean[1:],
                            label='omega**2')

                ax.set_yscale('log')
                ax.set_xscale('log')
                ax.set_xlabel('m+1')
                ax.set_xlim(1, self.index[-1])
                ax.legend(loc='best', frameon=False)
                fig.tight_layout()


class Pizza2DSpectrum(PizzaSetup):
    """
    This class can be used to read and display the 2D spectra '2D_spec_avg.TAG'

    >>> # display the content of 2D_spec_avg.TAG
    >>> # where TAG is the most recent file in the current directory
    >>> sp = Pizza2DSpectrum()
    """

    def __init__(self,  datadir='.', iplot=False, tag=None, endian='l',
                 all=False):
        """
        :param iplot: display the output plot when set to True (default is
                      False)
        :type iplot: bool
        :param tag: file suffix (tag), if not specified the most recent one in
                    the current directory is chosen
        :type tag: str
        :param datadir: current working directory
        :type datadir: str
        :param endian: control the endianness of the file
        :type endian: str
        :param all: when set to True, the complete time series is reconstructed
                    by stacking all the corresponding files from the working
                    directory (False by default)
        :type all: bool
        """

        self.name = '2D_spec_avg'

        if not all:
            if tag is not None:
                pattern = os.path.join(datadir, '%s.%s' % (self.name, tag))
                files = scanDir(pattern)
                # Either the log.tag directly exists and the setup is
                # easy to obtain
                if os.path.exists(os.path.join(datadir, 'log.%s' % tag)):
                    PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                        nml='log.%s' % tag)
                # Or the tag is a bit more complicated and we need to find
                # the corresponding log file
                else:
                    mask = re.compile(r'%s/%s\.(.*)' % (datadir, self.name))
                    if mask.match(files[-1]):
                        ending = mask.search(files[-1]).groups(0)[0]
                        pattern = os.path.join(datadir, 'log.%s' % ending)
                        if os.path.exists(pattern):
                            PizzaSetup.__init__(self, datadir=datadir,
                                                quiet=True,
                                                nml='log.%s' % ending)

                # Sum the files that correspond to the tag
                mask = re.compile(r'%s\.(.*)' % self.name)
                for k, file in enumerate(files):
                    print('reading %s' % file)
                    tag = mask.search(file).groups(0)[0]
                    nml = PizzaSetup(nml='log.%s' % tag, datadir=datadir,
                                     quiet=True)
                    filename = file
                    if k == 0:
                        self.tstart = nml.start_time
                        self.tstop = nml.stop_time  # will be replaced later
                        data = self.read(filename, endian)
                    else:
                        if os.path.exists(filename):
                            tmp = self.read(filename, endian)
                            data = self.add(data, tmp, nml.stop_time,
                                            nml.start_time)
            else:
                pattern = os.path.join(datadir, '%s*' % self.name)
                files = scanDir(pattern)
                filename = files[-1]
                print('reading %s' % filename)
                # Determine the setup
                mask = re.compile(r'%s\.(.*)' % self.name)
                ending = mask.search(files[-1]).groups(0)[0]
                if os.path.exists('log.%s' % ending):
                    try:
                        PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                            nml='log.%s' % ending)
                    except AttributeError:
                        pass

                data = self.read(filename, endian)

        else:  # if all is requested
            pattern = os.path.join(datadir, '%s.*' % self.name)
            files = scanDir(pattern)

            # Determine the setup
            mask = re.compile(r'%s\.(.*)' % self.name)
            for k, file in enumerate(files):
                print('reading %s' % file)
                tag = mask.search(file).groups(0)[0]
                nml = PizzaSetup(nml='log.%s' % tag, datadir=datadir,
                                 quiet=True)
                filename = file
                if k == 0:
                    self.tstart = nml.start_time
                    self.tstop = nml.stop_time  # will be replaced afterwards
                    data = self.read(filename, endian)
                else:
                    if os.path.exists(filename):
                        tmp = self.read(filename, endian)
                        data = self.add(data, tmp, nml.stop_time,
                                        nml.start_time)
            PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                nml='log.%s' % tag)

            # Determine the setup
            mask = re.compile(r'.*\.(.*)')
            ending = mask.search(files[-1]).groups(0)[0]
            if os.path.exists(os.path.join(datadir, 'log.%s' % ending)):
                PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % ending)

        self.assemble(data)

        if iplot:
            self.plot()

    def add(self, data, tmp, stop_time, start_time):
        """
        Clean way to stack data
        """
        if data.shape[0] >= tmp.shape[0]:
            out = copy.deepcopy(data)
        else:
            out = copy.deepcopy(tmp)
        fac_old = self.tstop-self.tstart
        fac_new = stop_time-start_time
        self.tstop = stop_time
        fac_tot = self.tstop-self.tstart

        if data.shape[0] == tmp.shape[0]:
            if self.version == 1:
                out = (fac_old*data+fac_new*tmp)/fac_tot
            else:
                for j in [0, 2, 4]:
                    out[j, ...] = (fac_old*data[j, ...] +
                                   fac_new*tmp[j, ...]) / fac_tot
                for j in [1, 3, 5]:
                    out[j, ...] = np.sqrt((fac_old*data[j, ...]**2 +
                                           fac_new*tmp[j, ...]**2) / fac_tot)
        else:
            if tmp.shape > data.shape:
                for j in [0, 1, 2]:
                    out[2*j, ...] = (fac_old*data[j, ...] +
                                     fac_new*tmp[2*j, ...]) / fac_tot
                for j in [1, 3, 5]:
                    out[j, ...] = tmp[j, ...]
            else:
                for j in [0, 1, 2]:
                    out[2*j, ...] = (fac_old*data[2*j, ...] +
                                     fac_new*tmp[j, ...]) / fac_tot
                for j in [1, 3, 5]:
                    out[j, ...] = data[j, ...]

        return out

    def read(self, filename, endian='l'):
        """
        :param filename: name of the input file
        :type filename: str
        :param endian: endianness of the binary file
        :type endian: str
        """
        file = open(filename, 'rb')
        dt = np.dtype("i4, 6f8")
        self.version, params = np.fromfile(file, dtype=dt, count=1)[0]
        self.ra, self.pr, self.raxi, self.sc, self.ek, self.radratio = params
        dt = np.dtype("4i4")
        self.n_r_max, self.n_m_max, self.m_max, self.minc = \
            np.fromfile(file, dtype=dt, count=1)[0]

        dt = np.dtype("%if8" % self.n_r_max)
        self.radius = np.fromfile(file, dtype=dt, count=1)[0]
        self.idx2m = np.zeros(self.n_m_max)
        for i in range(self.n_m_max):
            self.idx2m[i] = i*self.minc

        dt = np.dtype("(%i,%i)f8" % (self.n_r_max, self.n_m_max))
        if self.version == 1:
            data = np.fromfile(file, dtype=dt, count=3)
        elif self.version == 2:
            data = np.fromfile(file, dtype=dt, count=6)

        file.close()

        return data

    def assemble(self, data):
        """
        Once the reading is over then assemble the object

        :param data: a numpy array that contains all the time-averaged fields
        :type data: numpy.ndarray
        """
        if data.shape[0] == 3:
            self.us2_mean = data[0, ...].T
            self.up2_mean = data[1, ...].T
            self.enst_mean = data[2, ...].T
        elif data.shape[0] == 6:
            self.us2_mean = data[0, ...].T
            self.us2_std = data[1, ...].T
            self.up2_mean = data[2, ...].T
            self.up2_std = data[3, ...].T
            self.enst_mean = data[4, ...].T
            self.enst_std = data[4, ...].T

        self.ekin_mean = self.us2_mean+self.up2_mean

        ind = self.us2_mean[:, 1:-1].argmax(axis=0)
        self.peaks = np.zeros_like(ind)
        for k, idx in enumerate(ind):
            self.peaks[k] = self.idx2m[idx]

    def plot_rad(self, r=1.0):
        """
        This routine allows to plot the spectra at a given radius

        :param r: the radius one wants to display
        :type r: float
        """

        idx = np.where(abs(self.radius-r) == abs(self.radius-r).min(), 1, 0)
        idx = np.nonzero(idx)[0][0]

        fig = plt.figure()
        ax = fig.add_subplot(111)

        # if hasattr(self, 'us2_std'):
        #     sd = self.us2_std[1:,idx]
        #     ax.fill_between(self.idx2m[1:]+1, self.us2_mean[1:,idx]-sd, \
        #                     self.us2_mean[1:,idx]+sd, alpha=0.1)
        ax.plot(self.idx2m[1:]+1, self.us2_mean[1:, idx], label='us2')
        # if hasattr(self, 'up2_std'):
        #     sd = self.up2_std[:,idx]
        #     ax.fill_between(self.idx2m+1, self.up2_mean[:,idx]-sd, \
        #                     self.up2_mean[:,idx]+sd, alpha=0.1)
        ax.plot(self.idx2m+1, self.up2_mean[:, idx], label='up2')
        ax.plot(self.idx2m+1, self.ekin_mean[:, idx], label='Ekin')
        # if hasattr(self, 'enst_std'):
        #     sd = self.enst_std[1:,idx]
        #     ax.fill_between(self.idx2m[1:], self.enst_mean[1:,idx]-sd, \
        #                     self.enst_mean[1:,idx]+sd, alpha=0.1)
        ax.plot(self.idx2m+1, self.enst_mean[:, idx], label='Enst')

        ax.set_xscale('log')
        ax.set_yscale('log')

        ax.set_xlabel('m+1')
        ax.set_ylabel('Spectra')

        ax.legend(loc='best', frameon=False)

    def plot(self, levels=17, cm='magma', cut=1., solid_contour=True,
             log_yscale=True):
        """
        Plotting function

        :param cm: name of the colormap ('jet', 'seismic', 'RdYlBu_r', etc.)
        :type cm: str
        :param levels: the number of levels used in the contour plot
        :type levels: int
        :param  cut: a coefficient to change the dynamics of the contour levels
        :param cut: float
        :param solid_contour: a boolean to decide whether one also wants the
                              solid contour lines
        :type solid_contour: bool
        :param log_yscale: a boolean to decide whether one wants a logarithmic
                           y-axis
        :type log_yscale: bool
        """

        vmax = cut*np.log10(self.us2_mean[1:, :]+1e-34).max()
        vmin = cut*np.log10(self.us2_mean[self.us2_mean > 1e-15]).min()
        levs = np.linspace(vmin, vmax, levels)

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        im = ax1.contourf(self.radius[1:-1], self.idx2m[1:],
                          np.log10(self.us2_mean[1:, 1:-1]+1e-20),
                          levs, extend='both', cmap=plt.get_cmap(cm))
        if solid_contour:
            ax1.contour(self.radius[1:-1], self.idx2m[1:],
                        np.log10(self.us2_mean[1:, 1:-1]+1e-20), levs,
                        extend='both', linestyles=['-'], colors=['k'],
                        linewidths=[0.5])
        ax1.plot(self.radius[1:-1], self.peaks, ls='--')
        ax1.set_title('us**2')
        if log_yscale:
            ax1.set_yscale('log')
        ax1.set_xlabel('Radius')
        ax1.set_ylabel('m')

        # vmax = cut*np.log10(self.enst[1:,:]+1e-34).max()
        # vmin = cut*np.log10(self.enst[self.enst>1e-15]).min()
        # levs = np.linspace(vmin, vmax, levels)
        # ax2 = fig.add_subplot(122, sharey=ax1, sharex=ax1)
        # im = ax2.contourf(self.radius, self.idx2m[1:],
        #                   np.log10(self.enst[1:]+1e-20),
        #                   levs, extend='both', cmap=plt.get_cmap(cm))
        # if solid_contour:
        #     ax2.contour(self.radius, self.idx2m[1:],
        #                 np.log10(self.enst[1:]+1e-20), levs,
        #                 extend='both', linestyles=['-'], colors=['k'],
        #                 linewidths=[0.5])
        # ax2.set_title(r'Enstrophy')
        # if log_yscale: ax2.set_yscale('log')
        # ax2.set_xlabel('Radius')
        #
        # plt.setp(ax2.get_yticklabels(), visible=False)
        # ax1.set_xlim(self.radius[-1], self.radius[0])

        fig.colorbar(im)
        fig.tight_layout()
