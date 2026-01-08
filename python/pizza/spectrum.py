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
                 ave=True, tag=None, all=False, quiet=False):
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
        :param quiet: when set to True, makes the output silent (default False)
        :type quiet: bool
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

        self._speclut = None # None by default: to be filled later

        if not all:
            if tag is not None:
                if ispec is not None:
                    self.name += f'{ispec}'
                pattern = os.path.join(datadir, f'{self.name}.{tag}')
                files = scanDir(pattern)
                # Either the log.tag directly exists and the setup is easy
                # to obtain
                if os.path.exists(os.path.join(datadir, f'log.{tag}')):
                    PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                        nml=f'log.{tag}')
                # Or the tag is a bit more complicated and we need to find
                # the corresponding log file
                else:
                    mask = re.compile(fr'{datadir}/{self.name}\.(.*)')
                    if mask.match(files[-1]):
                        ending = mask.search(files[-1]).groups(0)[0]
                        pattern = os.path.join(datadir, f'log.{ending}')
                        if os.path.exists(pattern):
                            PizzaSetup.__init__(self, datadir=datadir,
                                                quiet=True, nml=f'log.{ending}')

                # Sum the files that correspond to the tag
                mask = re.compile(fr'{self.name}\.(.*)')
                for k, file in enumerate(files):
                    if not quiet:
                        print(f'reading {file}')
                    tag = mask.search(file).groups(0)[0]
                    nml = PizzaSetup(nml=f'log.{tag}', datadir=datadir,
                                     quiet=True)
                    filename = file
                    data = fast_read(filename)

                    if k == 0:
                        self._speclut = SpecLookUpTable(data, self.name,
                                                        nml.start_time,
                                                        nml.stop_time)
                    else:
                        self._speclut += SpecLookUpTable(data, self.name,
                                                         nml.start_time,
                                                         nml.stop_time)
            else:
                if ispec is not None:
                    pattern = os.path.join(datadir, f'{self.name}{ispec}*')
                else:
                    pattern = os.path.join(datadir, f'{self.name}*')
                files = scanDir(pattern)
                filename = files[-1]
                if not quiet:
                    print(f'reading {filename}')
                # Determine the setup
                if ispec is not None:
                    mask = re.compile(fr'{self.name}{ispec}\.(.*)')

                else:
                    mask = re.compile(fr'{self.name}\.(.*)')
                ending = mask.search(files[-1]).groups(0)[0]
                if os.path.exists(f'log.{ending}'):
                    try:
                        PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                            nml=f'log.{ending}')
                    except AttributeError:
                        self.start_time = None
                        self.stop_time = None
                        pass

                if not hasattr(self, 'stop_time'):
                    self.stop_time = None
                data = fast_read(filename, skiplines=0)
                self._speclut = SpecLookUpTable(data, self.name, self.start_time,
                                                self.stop_time)

        else:  # if all is requested
            pattern = os.path.join(datadir, f'{self.name}.*')
            files = scanDir(pattern)
            for k, file in enumerate(files):
                if not quiet:
                    print(f'reading {file}')

                mask = re.compile(fr'{self.name}\.(.*)')
                ending = mask.search(files[-1]).groups(0)[0]
                nml = PizzaSetup(nml=f'log.{ending}', datadir=datadir, quiet=True)
                data = fast_read(file)
                if k == 0:
                    self._speclut = SpecLookUpTable(data, self.name,
                                                    nml.start_time,
                                                    nml.stop_time)
                else:
                    self._speclut += SpecLookUpTable(data, self.name,
                                                     nml.start_time,
                                                     nml.stop_time)
            if os.path.exists(os.path.join(datadir, f'log.{ending}')):
                PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml=f'log.{ending}')

        # Copy look-up table arguments into MagicSpectrum object
        if self._speclut is not None:
            for attr in self._speclut.__dict__:
                setattr(self, attr, self._speclut.__dict__[attr])

        if iplot:
            self.plot()

    def __add__(self, new):
        """
        Clean way to stack data
        """
        out = copy.deepcopy(new)

        out._speclut += self._speclut
        for attr in out._speclut.__dict__:
            setattr(out, attr, out._speclut.__dict__[attr])

        return out

    def plot(self):
        """
        Plotting function
        """
        if self.name == 'vort_terms_avg':
            fig, ax = plt.subplots()
            sd = self.buo_m_SD/np.sqrt(self.buo_m)/2.
            ax.fill_between(self.index, np.sqrt(self.buo_m)-sd,
                            np.sqrt(self.buo_m)+sd, alpha=0.1)
            ax.plot(self.index, np.sqrt(self.buo_m), label='Buoyancy')

            sd = self.cor_m_SD/np.sqrt(self.cor_m)/2.
            ax.fill_between(self.index, np.sqrt(self.cor_m)-sd,
                            np.sqrt(self.cor_m)+sd, alpha=0.1)
            ax.plot(self.index, np.sqrt(self.cor_m), label='Coriolis')

            sd = self.iner_m_SD/np.sqrt(self.iner_m)/2.
            ax.fill_between(self.index, np.sqrt(self.iner_m)-sd,
                            np.sqrt(self.iner_m)+sd, alpha=0.1)
            ax.plot(self.index, np.sqrt(self.iner_m), label='Inertia')

            sd = self.visc_m_SD/np.sqrt(self.visc_m)/2.
            ax.fill_between(self.index, np.sqrt(self.visc_m)-sd,
                            np.sqrt(self.visc_m)+sd, alpha=0.1)
            ax.plot(self.index, np.sqrt(self.visc_m), label='Viscosity')

            sd = self.pump_m_SD/np.sqrt(self.pump_m)/2.
            ax.fill_between(self.index, np.sqrt(self.pump_m)-sd,
                            np.sqrt(self.pump_m)+sd, alpha=0.1)
            ax.plot(self.index, np.sqrt(self.pump_m), label='Ekman pumping')

            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.set_xlabel('m')
            ax.set_ylabel('Forces')
            ax.set_xlim(1, self.index[-1])
            ax.legend(loc='best', frameon=False)
            fig.tight_layout()

        else:
            if not self.ave:
                fig, ax = plt.subplots()
                if abs(self.up2_m[0]) > 0.:
                    ax.loglog(self.index[1:]+1, self.us2_m[1:], label='us**2')
                    ax.loglog(self.index+1, self.up2_m, label='up**2')
                    ax.loglog(self.index+1, self.enst_m, label='omega**2')
                else:
                    ax.loglog(self.index[1:]+1, self.us2_m[1:], label='us**2')
                    ax.loglog(self.index[1:]+1, self.up2_m[1:], label='up**2')
                    ax.loglog(self.index[1:]+1, self.enst_m[1:],
                              label='omega**2')
                ax.set_xlabel('m+1')
                ax.set_xlim(1, self.index[-1])
                ax.legend(loc='best', frameon=False)
                fig.tight_layout()
                if abs(self.temp2_m).max() > 0 or abs(self.xi2_m).max() > 0:
                    fig1, ax1 = plt.subplots()
                    if abs(self.temp2_m).max() > 0:
                        ax1.loglog(self.index[1:], self.temp2_m[1:], label='T**2')
                    if abs(self.xi2_m).max() > 0:
                        ax1.loglog(self.index[1:], self.xi2_m[1:], label='xi**2')
                    ax1.set_xlabel('m')
                    ax1.set_xlim(1, self.index[-1])
                    ax1.legend(loc='best', frameon=False)
                    fig1.tight_layout()
            else:
                fig, ax = plt.subplots()

                if abs(self.up2_m[0]) > 0.:
                    ax.fill_between(self.index[1:]+1,
                                    self.us2_m[1:]-self.us2_m_SD[1:],
                                    self.us2_m[1:]+self.us2_m_SD[1:],
                                    alpha=0.1)
                    ax.plot(self.index[1:]+1, self.us2_m[1:], label='us**2')

                    ax.fill_between(self.index+1, self.up2_m-self.up2_m_SD,
                                    self.up2_m+self.up2_m_SD, alpha=0.1)
                    ax.plot(self.index+1, self.up2_m, label='up**2')

                    ax.fill_between(self.index+1, self.enst_m-self.enst_std,
                                    self.enst_m+self.enst_std, alpha=0.1)
                    ax.plot(self.index+1, self.enst_m, label='omega**2')
                else:
                    ax.fill_between(self.index[1:]+1,
                                    self.us2_m[1:]-self.us2_m_SD[1:],
                                    self.us2_m[1:]+self.us2_m_SD[1:],
                                    alpha=0.1)
                    ax.plot(self.index[1:]+1, self.us2_m[1:], label='us**2')

                    ax.fill_between(self.index[1:]+1,
                                    self.up2_m[1:]-self.up2_m_SD[1:],
                                    self.up2_m[1:]+self.up2_m_SD[1:],
                                    alpha=0.1)
                    ax.plot(self.index[1:]+1, self.up2_m[1:], label='up**2')

                    ax.fill_between(self.index[1:]+1,
                                    self.enst_m[1:]-self.enst_std[1:],
                                    self.enst_m[1:]+self.enst_std[1:],
                                    alpha=0.1)
                    ax.plot(self.index[1:]+1, self.enst_m[1:],
                            label='omega**2')

                ax.set_yscale('log')
                ax.set_xscale('log')
                ax.set_xlabel('m+1')
                ax.set_xlim(1, self.index[-1])
                ax.legend(loc='best', frameon=False)
                fig.tight_layout()

                if abs(self.temp2_m).max() > 0 or abs(self.xi2_m).max() > 0:
                    fig1, ax1 = plt.subplots()
                    if abs(self.temp2_m).max() > 0:
                        ax1.fill_between(self.index[1:],
                                         self.temp2_m[1:]-self.temp2_m_SD[1:],
                                         self.temp2_m[1:]+self.temp2_m_SD[1:],
                                         alpha=0.1)
                        ax1.plot(self.index[1:], self.temp2_m[1:], label='temp**2')
                    if abs(self.xi2_m).max() > 0:
                        ax1.fill_between(self.index[1:],
                                        self.xi2_m[1:]-self.xi2_m_SD[1:],
                                        self.xi2_m[1:]+self.xi2_m_SD[1:],
                                        alpha=0.1)
                        ax1.plot(self.index[1:], self.xi2_m[1:], label='xi**2')

                    ax1.set_yscale('log')
                    ax1.set_xscale('log')
                    ax1.set_xlabel('m')
                    ax1.set_xlim(1, self.index[-1])
                    ax1.legend(loc='best', frameon=False)
                    fig1.tight_layout()


class SpecLookUpTable:
    """
    The purpose of this class is to create a lookup table between the numpy
    array that comes from the reading of the spec files and the corresponding
    columns.
    """

    def __init__(self, data, name, tstart=None, tstop=None):
        """
        :param data: numpy array that contains the data
        :type data: numpy.ndarray
        :param name: name of the field
        :type name: str
        :param tstart: starting time that was used to compute the time average
        :type tstart: float
        :param tstop: stop time that was used to compute the time average
        :type tstop: float
        """

        self.name = name
        self.start_time = tstart
        self.stop_time = tstop

        self.index = data[:, 0]

        if self.name == 'vort_terms_avg':
            self.buo_m = data[:, 1]
            self.buo_m_SD = data[:, 2]
            self.cor_m = data[:, 3]
            self.cor_m_SD = data[:, 4]
            self.adv_m = data[:, 5]
            self.adv_std = data[:, 6]
            self.domdt_m = data[:, 7]
            self.domdt_m_SD = data[:, 8]
            self.visc_m = data[:, 9]
            self.visc_m_SD = data[:, 10]
            self.pump_m = data[:, 11]
            self.pump_m_SD = data[:, 12]
            self.thwind_m = data[:, 13]
            self.thwind_m_SD = data[:, 14]
            self.iner_m = data[:, 15]
            self.iner_m_SD = data[:, 16]
            self.cia_m = data[:, 17]
            self.cia_m_SD = data[:, 18]
        elif self.name == 'spec_avg':
            self.us2_m = data[:, 1]
            self.us2_m_SD = data[:, 2]
            self.up2_m = data[:, 3]
            self.up2_m_SD = data[:, 4]
            self.enst_m = data[:, 5]
            self.enst_std = data[:, 6]
            if data.shape[1] > 7:
                self.temp2_m = data[:, 7]
                self.temp2_m_SD = data[:, 8]
                self.xi2_m = data[:, 9]
                self.xi2_m_SD = data[:, 10]
            else:
                self.temp2_m = np.zeros_like(self.us2_m)
                self.temp2_m_SD = np.zeros_like(self.us2_m)
                self.xi2_m = np.zeros_like(self.us2_m)
                self.xi2_m_SD = np.zeros_like(self.us2_m)
            self.ekin_m = self.us2_m + self.up2_m
            self.ekin_m_SD = np.sqrt(self.us2_m_SD**2+self.up2_m_SD**2)
        else:
            self.us2_m = data[:, 1]
            self.up2_m = data[:, 2]
            self.enst_m = data[:, 3]
            if data.shape[1] > 4:
                self.temp2_m = data[:, 4]
                self.xi2_m = data[:, 5]
            else:
                self.temp2_m = np.zeros_like(self.us2_m)
                self.xi2_m = np.zeros_like(self.us2_m)
            self.ekin_m = self.us2_m + self.up2_m

    def __add__(self, new):
        """
        This is a python built-in method to stack two look-up tables.
        """

        out = copy.deepcopy(new)
        if self.start_time is not None:
            fac_old = self.stop_time-self.start_time
            out.start_time = self.start_time
        else:
            fac_old = 0.
        if new.stop_time is not None:
            fac_new = new.stop_time-new.start_time
            out.stop_time = new.stop_time
        else:
            fac_new = 0.
        if fac_old != 0 or fac_new != 0:
            fac_tot = fac_new+fac_old
        else:
            fac_tot = 1.

        idx_old_max = len(self.index)
        idx_new_max = len(new.index)

        if idx_old_max == idx_new_max:
            for attr in new.__dict__.keys():
                if attr not in ['index', 'name', 'start_time', 'stop_time']:
                    # Only stack if both new and old have the attribute available
                    if attr in self.__dict__:
                        # Standard deviation
                        if attr.endswith('SD'):
                            if abs(self.__dict__[attr]).max() > 0.:
                                out.__dict__[attr] = np.sqrt(( \
                                                  fac_new*new.__dict__[attr]**2 + \
                                                  fac_old*self.__dict__[attr]**2) / \
                                                  fac_tot)
                            else:
                                out.__dict__[attr] = new.__dict__[attr]
                        # Regular field
                        else:
                            if abs(self.__dict__[attr]).max() > 0.:
                                out.__dict__[attr] = (fac_new*new.__dict__[attr] + \
                                                      fac_old*self.__dict__[attr]) / \
                                                      fac_tot
                            else:
                                out.__dict__[attr] = new.__dict__[attr]
        else: # Different truncations
            idx_min = min(idx_old_max, idx_new_max)
            for attr in new.__dict__.keys():
                if attr not in ['index', 'name', 'start_time', 'stop_time']:
                    # Only stack if both new and old have the attribute available
                    if attr in self.__dict__:
                        # Standard deviation
                        if attr.endswith('SD'):
                            if abs(self.__dict__[attr]).max() > 0.:
                                out.__dict__[attr][:idx_min] = \
                                   np.sqrt((fac_new*new.__dict__[attr][:idx_min]**2+\
                                            fac_old*self.__dict__[attr][:idx_min]**2)\
                                            /fac_tot)
                            else:
                                out.__dict__[attr] = new.__dict__[attr]
                        # Regular field
                        else:
                            if abs(self.__dict__[attr]).max() > 0.:
                                out.__dict__[attr][:idx_min] = \
                                   (fac_new*new.__dict__[attr][:idx_min] + \
                                    fac_old*self.__dict__[attr][:idx_min]) / fac_tot
                            else:
                                out.__dict__[attr] = new.__dict__[attr]

        return out


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
                pattern = os.path.join(datadir, f'{self.name}.{tag}')
                files = scanDir(pattern)
                # Either the log.tag directly exists and the setup is
                # easy to obtain
                if os.path.exists(os.path.join(datadir, f'log.{tag}')):
                    PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                        nml=f'log.{tag}')
                # Or the tag is a bit more complicated and we need to find
                # the corresponding log file
                else:
                    mask = re.compile(fr'{datadir}/{self.name}\.(.*)')
                    if mask.match(files[-1]):
                        ending = mask.search(files[-1]).groups(0)[0]
                        pattern = os.path.join(datadir, f'log.{ending}')
                        if os.path.exists(pattern):
                            PizzaSetup.__init__(self, datadir=datadir,
                                                quiet=True, nml=f'log.{ending}')

                # Sum the files that correspond to the tag
                mask = re.compile(fr'{self.name}\.(.*)')
                for k, file in enumerate(files):
                    print(f'reading {file}')
                    tag = mask.search(file).groups(0)[0]
                    nml = PizzaSetup(nml=f'log.{tag}', datadir=datadir,
                                     quiet=True)
                    filename = file
                    if k == 0:
                        self.tstart = nml.start_time
                        self.tstop = nml.stop_time  # will be replaced later
                        data = self._read(filename, endian)
                    else:
                        if os.path.exists(filename):
                            tmp = self._read(filename, endian)
                            data = self.add(data, tmp, nml.stop_time,
                                            nml.start_time)
            else:
                pattern = os.path.join(datadir, f'{self.name}*')
                files = scanDir(pattern)
                filename = files[-1]
                print(f'reading {filename}')
                # Determine the setup
                mask = re.compile(fr'{self.name}\.(.*)')
                ending = mask.search(files[-1]).groups(0)[0]
                if os.path.exists(f'log.{ending}'):
                    try:
                        PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                            nml=f'log.{ending}')
                    except AttributeError:
                        pass

                data = self._read(filename, endian)

        else:  # if all is requested
            pattern = os.path.join(datadir, f'{self.name}.*')
            files = scanDir(pattern)

            # Determine the setup
            mask = re.compile(fr'{self.name}\.(.*)')
            for k, file in enumerate(files):
                print(f'reading {file}')
                tag = mask.search(file).groups(0)[0]
                nml = PizzaSetup(nml=f'log.{tag}', datadir=datadir,
                                 quiet=True)
                filename = file
                if k == 0:
                    self.tstart = nml.start_time
                    self.tstop = nml.stop_time  # will be replaced afterwards
                    data = self._read(filename, endian)
                else:
                    if os.path.exists(filename):
                        tmp = self._read(filename, endian)
                        data = self.add(data, tmp, nml.stop_time,
                                        nml.start_time)
            PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                nml=f'log.{tag}')

            # Determine the setup
            mask = re.compile(r'.*\.(.*)')
            ending = mask.search(files[-1]).groups(0)[0]
            if os.path.exists(os.path.join(datadir, f'log.{ending}')):
                PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml=f'log.{ending}')

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

    def _read(self, filename, endian='l'):
        """
        :param filename: name of the input file
        :type filename: str
        :param endian: endianness of the binary file
        :type endian: str
        """
        with open(filename, 'rb') as file:
            self.version = np.fromfile(file, dtype=np.int32, count=1)[0]
            self.ra, self.pr, self.raxi, self.sc, self.ek, self.radratio = \
                np.fromfile(file, dtype=np.float64, count=6)
            self.n_r_max, self.n_m_max, self.m_max, self.minc = \
                np.fromfile(file, dtype=np.int32, count=4)

            self.radius = np.fromfile(file, dtype=np.float64,
                                      count=self.n_r_max)
            self.idx2m = np.zeros(self.n_m_max)
            for i in range(self.n_m_max):
                self.idx2m[i] = i*self.minc

            dt = np.dtype("(%i,%i)f8" % (self.n_r_max, self.n_m_max))
            if self.version == 1:
                data = np.fromfile(file, dtype=dt, count=3)
            elif self.version == 2:
                data = np.fromfile(file, dtype=dt, count=6)

        return data

    def assemble(self, data):
        """
        Once the reading is over then assemble the object

        :param data: a numpy array that contains all the time-averaged fields
        :type data: numpy.ndarray
        """
        if data.shape[0] == 3:
            self.us2_m = data[0, ...].T
            self.up2_m = data[1, ...].T
            self.enst_m = data[2, ...].T
        elif data.shape[0] == 6:
            self.us2_m = data[0, ...].T
            self.us2_m_SD = data[1, ...].T
            self.up2_m = data[2, ...].T
            self.up2_m_SD = data[3, ...].T
            self.enst_m = data[4, ...].T
            self.enst_std = data[4, ...].T

        self.ekin_m = self.us2_m+self.up2_m

        ind = self.us2_m[:, 1:-1].argmax(axis=0)
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

        fig, ax = plt.subplots()
        # if hasattr(self, 'us2_m_SD'):
        #     sd = self.us2_m_SD[1:,idx]
        #     ax.fill_between(self.idx2m[1:]+1, self.us2_m[1:,idx]-sd, \
        #                     self.us2_m[1:,idx]+sd, alpha=0.1)
        ax.plot(self.idx2m[1:]+1, self.us2_m[1:, idx], label='us2')
        # if hasattr(self, 'up2_m_SD'):
        #     sd = self.up2_m_SD[:,idx]
        #     ax.fill_between(self.idx2m+1, self.up2_m[:,idx]-sd, \
        #                     self.up2_m[:,idx]+sd, alpha=0.1)
        ax.plot(self.idx2m+1, self.up2_m[:, idx], label='up2')
        ax.plot(self.idx2m+1, self.ekin_m[:, idx], label='Ekin')
        # if hasattr(self, 'enst_std'):
        #     sd = self.enst_std[1:,idx]
        #     ax.fill_between(self.idx2m[1:], self.enst_m[1:,idx]-sd, \
        #                     self.enst_m[1:,idx]+sd, alpha=0.1)
        ax.plot(self.idx2m+1, self.enst_m[:, idx], label='Enst')

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

        vmax = cut*np.log10(self.us2_m[1:, :]+1e-34).max()
        vmin = cut*np.log10(self.us2_m[self.us2_m > 1e-15]).min()
        levs = np.linspace(vmin, vmax, levels)

        fig, ax1 = plt.subplots()
        im = ax1.contourf(self.radius[1:-1], self.idx2m[1:],
                          np.log10(self.us2_m[1:, 1:-1]+1e-20),
                          levs, extend='both', cmap=plt.get_cmap(cm))
        if solid_contour:
            ax1.contour(self.radius[1:-1], self.idx2m[1:],
                        np.log10(self.us2_m[1:, 1:-1]+1e-20), levs,
                        extend='both', linestyles=['-'], colors=['k'],
                        linewidths=[0.5])
        ax1.plot(self.radius[1:-1], self.peaks, ls='--')
        ax1.set_title('us**2')
        if log_yscale:
            ax1.set_yscale('log')
        ax1.set_xlabel('Radius')
        ax1.set_ylabel('m')

        # vmax = cut*np.log10(self.enst_m[1:,:]+1e-34).max()
        # vmin = cut*np.log10(self.enst_m[self.enst_m>1e-15]).min()
        # levs = np.linspace(vmin, vmax, levels)
        # ax2 = fig.add_subplot(122, sharey=ax1, sharex=ax1)
        # im = ax2.contourf(self.radius, self.idx2m[1:],
        #                   np.log10(self.enst_m[1:]+1e-20),
        #                   levs, extend='both', cmap=plt.get_cmap(cm))
        # if solid_contour:
        #     ax2.contour(self.radius, self.idx2m[1:],
        #                 np.log10(self.enst_m[1:]+1e-20), levs,
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
