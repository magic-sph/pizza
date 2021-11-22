# -*- coding: utf-8 -*-
from .npfile import npfile
import numpy as np
import matplotlib.pyplot as plt
import copy
from .log import PizzaSetup
from .libpizza import scanDir, costf
from scipy.integrate import simps
import os
import re


def avg_std(time, y):
    """
    This routine computes the time-average and the standard deviation of
    an input array

    :param time: the time vector
    :type time: numpy.ndarray
    :param y: the quantity one wants to time average
    :type y: numpy.ndarray
    """

    tfac = 1./(time[-1]-time[0])
    ymean = tfac*simps(y, time, axis=0)
    ystd = np.sqrt(tfac*simps((y-ymean)**2, time, axis=0))

    return ymean, ystd


def interp_dct(arr, nr_new):
    """
    This routine interpolates an input array of size (n_snap, n_r_max)
    onto a new grid of size (n_snap, n_r_max_new) using Discrete Cosine
    Transforms
    """

    nr = arr.shape[-1]
    # Transform to Chebyshev space and renormalise for the new grid
    fhat = costf(arr)/np.sqrt((nr-1.)/(nr_new-1.))
    fout = np.zeros((arr.shape[0], nr_new), np.float64)
    # Set the first nr Cheb coefficients and then pad with zeros
    fout[:, :nr] = fhat
    # Finally bring back the array on the physical grid
    fout = costf(fout)

    return fout


class PizzaBalance(PizzaSetup):
    """
    This module is used to load and display the output files that
    contain the azimuthal force balance vphi_bal.TAG

    >>> bal = PizzaBalance(all=True)
    >>> bal.plot()
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
            pattern = os.path.join(datadir, 'vphi_bal.%s' % tag)
            files = scanDir(pattern)

            #  Either the log.tag directly exists and the setup is easy
            #  to obtain
            if os.path.exists(os.path.join(datadir, 'log.%s' % tag)):
                PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % tag)
            #  Or the tag is a bit more complicated and we need to find
            #  the corresponding log file
            else:
                mask = re.compile(r'%s/vphi_bal\.(.*)' % (datadir))
                if mask.match(files[-1]):
                    ending = mask.search(files[-1]).groups(0)[0]
                    pattern = os.path.join(datadir, 'log.%s' % ending)
                    if os.path.exists(pattern):
                        PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                            nml='log.%s' % ending)

            for k, file in enumerate(files):
                print('reading %s' % file)
                if k == 0:
                    self.radius, self.time, self.vp, self.dvpdt, \
                        self.rey_stress, self.ek_pump, self.visc = \
                        self.read(file, endian)
                else:
                    radius, time, vp, dvpdt, rey_stress, ek_pump, visc = \
                                              self.read(file, endian)
                    self.add(time, radius, vp, dvpdt, rey_stress, ek_pump,
                             visc)

        # If no tag is specified, the most recent is plotted
        elif not all:
            if len(logFiles) != 0:
                PizzaSetup.__init__(self, quiet=True, nml=logFiles[-1])
                name = 'vphi_bal.%s' % self.tag
                filename = os.path.join(datadir, name)
                print('reading %s' % filename)
                self.radius, self.time, self.vp, self.dvpdt, self.rey_stress, \
                    self.ek_pump, self.visc = self.read(filename, endian)
            else:
                dat = scanDir('vphi_bal.*')
                filename = dat[-1][1]
                print('reading %s' % filename)
                self.radius, self.time, self.vp, self.dvpdt, self.rey_stress, \
                    self.ek_pump, self.visc = self.read(filename, endian)

        # If no tag is specified but all=True, all the directory is plotted
        else:
            if len(logFiles) != 0:
                PizzaSetup.__init__(self, quiet=True, nml=logFiles[-1])
            pattern = os.path.join(datadir, 'vphi_bal.*')
            files = scanDir(pattern)
            for k, file in enumerate(files):
                print('reading %s' % file)
                if k == 0:
                    self.radius, self.time, self.vp, self.dvpdt, \
                        self.rey_stress, self.ek_pump, self.visc = \
                        self.read(file, endian)
                else:
                    radius, time, vp, dvpdt, rey_stress, ek_pump, visc = \
                        self.read(file, endian)
                    self.add(time, radius, vp, dvpdt, rey_stress, ek_pump,
                             visc)

        if iplot:
            self.plot()

    def __add__(self, new):
        """
        Built-in functions to sum two vphi force balances
        """
        out = copy.deepcopy(new)
        nr_new = new.vp.shape[-1]
        nr_max = self.vp.shape[-1]
        if nr_max != nr_new:
            self.vp = interp_dct(self.vp, nr_new)
            self.dvpdt = interp_dct(self.dvpdt, nr_new)
            self.rey_stress = interp_dct(self.rey_stress, nr_new)
            self.ek_pump = interp_dct(self.ek_pump, nr_new)
            self.visc = interp_dct(self.visc, nr_new)

            # We have a situation here:
            self.radius = radius

        if new.time[0] == self.time[-1]:
            out.time = np.concatenate((self.time, new.time[1:]), axis=0)
            out.vp = np.concatenate((self.vp, new.vp[1:, :]), axis=0)
            out.dvpdt = np.concatenate((self.dvpdt, new.dvpdt[1:, :]), axis=0)
            out.rey_stress = np.concatenate((self.rey_stress,
                                             new.rey_stress[1:, :]), axis=0)
            out.ek_pump = np.concatenate((self.ek_pump, new.ek_pump[1:, :]),
                                         axis=0)
            out.visc = np.concatenate((self.visc, new.visc[1:, :]), axis=0)
        else:
            out.time = np.concatenate((self.time, new.time), axis=0)
            out.vp = np.concatenate((self.vp, new.vp), axis=0)
            out.dvpdt = np.concatenate((self.dvpdt, new.dvpdt), axis=0)
            out.rey_stress = np.concatenate((self.rey_stress, new.rey_stress),
                                            axis=0)
            out.ek_pump = np.concatenate((self.ek_pump, new.ek_pump), axis=0)
            out.visc = np.concatenate((self.visc, new.visc), axis=0)

        return out

    def add(self, time, radius, vp, dvpdt, rey_stress, ek_pump, visc):
        nr_new = vp.shape[-1]
        nr_max = self.vp.shape[-1]
        if nr_max != nr_new:
            self.vp = interp_dct(self.vp, nr_new)
            self.dvpdt = interp_dct(self.dvpdt, nr_new)
            self.rey_stress = interp_dct(self.rey_stress, nr_new)
            self.ek_pump = interp_dct(self.ek_pump, nr_new)
            self.visc = interp_dct(self.visc, nr_new)

            self.radius = radius

        if time[0] == self.time[-1]:
            self.time = np.concatenate((self.time, time[1:]), axis=0)
            self.vp = np.concatenate((self.vp, vp[1:, :]), axis=0)
            self.dvpdt = np.concatenate((self.dvpdt, dvpdt[1:, :]), axis=0)
            self.rey_stress = np.concatenate((self.rey_stress,
                                              rey_stress[1:, :]), axis=0)
            self.ek_pump = np.concatenate((self.ek_pump, ek_pump[1:, :]),
                                          axis=0)
            self.visc = np.concatenate((self.visc, visc[1:, :]), axis=0)
        else:
            self.time = np.concatenate((self.time, time), axis=0)
            self.vp = np.concatenate((self.vp, vp), axis=0)
            self.dvpdt = np.concatenate((self.dvpdt, dvpdt), axis=0)
            self.rey_stress = np.concatenate((self.rey_stress, rey_stress),
                                             axis=0)
            self.ek_pump = np.concatenate((self.ek_pump, ek_pump), axis=0)
            self.visc = np.concatenate((self.visc, visc), axis=0)

    def read(self, filename, endian):
        """
        Read one vphi_bal.TAG file
        """

        # Since n_r_max is not in the file, one needs to get from the log file
        # This is right now the only way to properly stack vphi_bal files that
        # have different radial resolution
        tag = filename.split('vphi_bal.')[-1]
        if os.path.exists('log.%s' % tag):
            stp = PizzaSetup(nml='log.%s' % tag, quiet=True)
            n_r_max = stp.n_r_max
        else:
            n_r_max = self.n_r_max

        try:  # Old file format (with record markers)
            file = npfile(filename, endian=endian)
            self.ra, self.ek, self.pr, self.radratio, self.raxi, self.sc \
                = file.fort_read(np.float64)
            radius = file.fort_read(np.float64)
            time = file.fort_read(np.float64)
            vp = file.fort_read(np.float64)
            dvpdt = file.fort_read(np.float64)
            rey_stress = file.fort_read(np.float64)
            ek_pump = file.fort_read(np.float64)
            visc = file.fort_read(np.float64)
            while 1:
                try:
                    time = np.append(time, file.fort_read(np.float64))
                    vp = np.vstack((vp, file.fort_read(np.float64)))
                    dvpdt = np.vstack((dvpdt, file.fort_read(np.float64)))
                    rey_stress = np.vstack((rey_stress,
                                            file.fort_read(np.float64)))
                    ek_pump = np.vstack((ek_pump, file.fort_read(np.float64)))
                    visc = np.vstack((visc, file.fort_read(np.float64)))
                except TypeError:
                    break

            file.close()
        except:  # New file format (without record marker)
            file = open(filename, 'rb')
            data = np.fromfile(file, dtype=np.float64)
            self.ra, self.ek, self.pr, self.radratio, self.raxi, self.sc = \
                data[0:6]
            radius = data[6:7+n_r_max-1]
            data = data[7+n_r_max-1:]
            nsteps = len(data)//(5*n_r_max+1)

            data = data.reshape((nsteps, 5*n_r_max+1))
            time = data[:, 0]
            vp = data[:, 1:n_r_max+1]
            dvpdt = data[:, n_r_max+1:2*n_r_max+1]
            rey_stress = data[:, 2*n_r_max+1:3*n_r_max+1]
            ek_pump = data[:, 3*n_r_max+1:4*n_r_max+1]
            visc = data[:, 4*n_r_max+1:5*n_r_max+1]

            file.close()

        return radius, time, vp, dvpdt, rey_stress, ek_pump, visc

    def write(self, filename):
        """
        This routine writes a snap using stream

        :param filename: the name of the output file
        :type filename: str
        """

        out = open('%s' % filename, 'wb')
        x = np.array([self.ra, self.ek, self.pr, self.radratio,
                      self.raxi, self.sc], dtype=np.float64)
        x.tofile(out)
        self.radius.tofile(out)

        nsteps = len(self.time)

        for i in range(nsteps):
            x = np.array([self.time[i]], dtype=np.float64)
            x.tofile(out)
            self.vp[i, :].tofile(out)
            self.dvpdt[i, :].tofile(out)
            self.rey_stress[i, :].tofile(out)
            self.ek_pump[i, :].tofile(out)
            self.visc[i, :].tofile(out)

        out.close()

    def plot(self, nstep=1):
        """
        Display when ``iplot = True``

        :param nstep: time interval (by default 1)
        :type nstep: int
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)
        dvpdtm, dvpdtstd = avg_std(self.time[::nstep], self.dvpdt[::nstep])
        ax.fill_between(self.radius, dvpdtm-dvpdtstd, dvpdtm+dvpdtstd,
                        alpha=0.1)
        ax.plot(self.radius, dvpdtm, label='dvp/dt')

        reym, reystd = avg_std(self.time[::nstep], self.rey_stress[::nstep])
        ax.fill_between(self.radius, reym-reystd, reym+reystd, alpha=0.1)
        ax.plot(self.radius, reym, label='Reynolds stress')

        ekpm, ekpstd = avg_std(self.time[::nstep], self.ek_pump[::nstep])
        ax.fill_between(self.radius, ekpm-ekpstd, ekpm+ekpstd, alpha=0.1)
        ax.plot(self.radius, ekpm, label='Ekman pumping')

        vim, vistd = avg_std(self.time[::nstep], self.visc[::nstep])
        ax.fill_between(self.radius, vim-vistd, vim+vistd, alpha=0.1)
        ax.plot(self.radius, vim, label='Viscosity')
        tot = self.rey_stress[::nstep]+self.ek_pump[::nstep]+self.visc[::nstep]
        ax.plot(self.radius, tot.mean(axis=0), label='Total')

        ax.set_xlim(self.radius[-1], self.radius[0])
        ax.legend(loc='best', frameon=False)
        ax.set_xlabel('Radius')
        ax.set_ylabel('Forces')
        fig.tight_layout()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        vpm = self.vp[::nstep].mean(axis=0)
        vpstd = self.vp[::nstep].std(axis=0)
        vpm_p = vpm + vpstd
        vpm_m = vpm - vpstd
        # xs, ys = mlab.poly_between(self.radius, vpm_m, vpm_p)
        ax.fill_between(self.radius, vpm_m, vpm_p, alpha=0.1)
        ax.plot(self.radius, vpm)
        # ax.plot(self.radius, vm)
        ax.set_xlim(self.radius[-1], self.radius[0])
        ax.set_xlabel('time')
        ax.set_ylabel('up')
        fig.tight_layout()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        vmax = abs(self.vp[::nstep]).max()
        cs = np.linspace(-vmax, vmax, 65)
        ax.contourf(self.time[::nstep], self.radius, self.vp[::nstep].T, cs,
                    cmap=plt.get_cmap('seismic'), extend='both')

        ax.set_xlabel('Time')
        ax.set_ylabel('Radius')

        """
        omz = get_dr(self.vp*self.radius)/self.radius
        #omz_filt = np.zeros_like(omz)
        #b, a = butter(3, 0.1)
        #for i in range(omz.shape[0]):
            #omz_filt[i, :] = filtfilt(b, a, omz[i, :])
        #omz = omz_filt
        domz = get_dr(omz)
        beta = np.zeros_like(self.radius)
        beta[1:] = -self.radius[1:]/(self.radius[0]**2-self.radius[1:]**2)
        crit = 2./self.ek*beta-domz

        #im =  ax.contourf(self.time, self.radius, crit.T, cs,
                          #cmap=plt.get_cmap('seismic'), extend='both')
        im =  ax.contour(self.time, self.radius, crit.T, colors=['k'])
        """

        fig.tight_layout()


class PizzaVortBalance(PizzaSetup):

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

        self.name = 'vort_bal'

        if not all:
            if tag is not None:
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
                    self.tstop = nml.stop_time  # will be overwritten later
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

    def read(self, filename, endian='l'):
        """
        :param filename: name of the input file
        :type filename: str
        :param endian: endianness of the binary file
        :type endian: str
        :returns data: a numpy array that contains all the time-averaged fields
        :type data: numpy.ndarray
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
            data = np.fromfile(file,  dtype=dt, count=9)
        elif self.version == 2:
            data = np.fromfile(file,  dtype=dt, count=18)

        file.close()

        return data

    def assemble(self, data):
        """
        Once the reading is over then assemble the object

        :param data: a numpy array that contains all the time-averaged fields
        :type data: numpy.ndarray
        """
        if data.shape[0] == 9:
            self.buo_mean = data[0, ...].T
            self.cor_mean = data[1, ...].T
            self.adv_mean = data[2, ...].T
            self.domdt_mean = data[3, ...].T
            self.visc_mean = data[4, ...].T
            self.pump_mean = data[5, ...].T
            self.thwind_mean = data[6, ...].T
            self.iner_mean = data[7, ...].T
            self.cia_mean = data[8, ...].T
        elif data.shape[0] == 18:
            self.buo_mean = data[0, ...].T
            self.buo_std = data[1, ...].T
            self.cor_mean = data[2, ...].T
            self.cor_std = data[3, ...].T
            self.adv_mean = data[4, ...].T
            self.adv_std = data[5, ...].T
            self.domdt_mean = data[6, ...].T
            self.domdt_std = data[7, ...].T
            self.visc_mean = data[8, ...].T
            self.visc_std = data[9, ...].T
            self.pump_mean = data[10, ...].T
            self.pump_std = data[11, ...].T
            self.thwind_mean = data[12, ...].T
            self.thwind_std = data[13, ...].T
            self.iner_mean = data[14, ...].T
            self.iner_std = data[15, ...].T
            self.cia_mean = data[16, ...].T
            self.cia_std = data[17, ...].T

    def add(self, data, tmp, stop_time, start_time):
        """
        Clean way to stack data

        :returns out: a numpy array that contains all the time-averaged fields
        :type out: numpy.ndarray
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
                for j in [0, 2, 4, 6, 8, 10, 12, 14, 16]:
                    out[j, ...] = (fac_old*data[j, ...] +
                                   fac_new*tmp[j, ...])/fac_tot
                for j in [1, 3, 5, 7, 9, 11, 13, 15, 17]:
                    out[j, ...] = np.sqrt((fac_old*data[j, ...]**2 +
                                           fac_new*tmp[j, ...]**2) / fac_tot)
        else:
            if tmp.shape > data.shape:
                for j in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
                    out[2*j, ...] = (fac_old*data[j, ...] +
                                     fac_new*tmp[2*j, ...]) / fac_tot
                for j in [1, 3, 5, 7, 9, 11, 13, 15, 17]:
                    out[j, ...] = tmp[j, ...]
            else:
                for j in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
                    out[2*j, ...] = (fac_old*data[2*j, ...] +
                                     fac_new*tmp[j, ...]) / fac_tot
                for j in [1, 3, 5, 7, 9, 11, 13, 15, 17]:
                    out[j, ...] = data[j, ...]

        return out

    def plot_rad(self, r=1.0):
        """
        This routine allows to display the force balance at a given radius

        :param r: the radius one wants to display
        :type r: float
        """

        idx = np.where(abs(self.radius-r) == abs(self.radius-r).min(), 1, 0)
        idx = np.nonzero(idx)[0][0]

        fig = plt.figure()
        ax = fig.add_subplot(111)

        if hasattr(self, 'buo_std'):
            sd = self.buo_std[1:, idx]/np.sqrt(self.buo_mean[1:, idx])/2.
            ax.fill_between(self.idx2m[1:], np.sqrt(self.buo_mean[1:, idx])-sd,
                            np.sqrt(self.buo_mean[1:, idx])+sd, alpha=0.1)
        ax.plot(self.idx2m[1:], np.sqrt(self.buo_mean[1:, idx]),
                label='Buoyancy')

        if hasattr(self, 'cor_std'):
            sd = self.cor_std[1:, idx]/np.sqrt(self.cor_mean[1:, idx])/2.
            ax.fill_between(self.idx2m[1:], np.sqrt(self.cor_mean[1:, idx])-sd,
                            np.sqrt(self.cor_mean[1:, idx])+sd, alpha=0.1)
        ax.plot(self.idx2m[1:], np.sqrt(self.cor_mean[1:, idx]),
                label='Coriolis')

        if hasattr(self, 'iner_std'):
            sd = self.iner_std[1:, idx]/np.sqrt(self.iner_mean[1:, idx])/2.
            ax.fill_between(self.idx2m[1:],
                            np.sqrt(self.iner_mean[1:, idx])-sd,
                            np.sqrt(self.iner_mean[1:, idx])+sd, alpha=0.1)
        ax.plot(self.idx2m[1:], np.sqrt(self.iner_mean[1:, idx]),
                label='Inertia')

        if hasattr(self, 'visc_std'):
            sd = self.visc_std[1:, idx]/np.sqrt(self.visc_mean[1:, idx])/2.
            ax.fill_between(self.idx2m[1:],
                            np.sqrt(self.visc_mean[1:, idx])-sd,
                            np.sqrt(self.visc_mean[1:, idx])+sd, alpha=0.1)
        ax.plot(self.idx2m[1:], np.sqrt(self.visc_mean[1:, idx]),
                label='Viscosity')

        if hasattr(self, 'pump_std'):
            sd = self.pump_std[1:, idx]/np.sqrt(self.pump_mean[1:, idx])/2.
            ax.fill_between(self.idx2m[1:],
                            np.sqrt(self.pump_mean[1:, idx])-sd,
                            np.sqrt(self.pump_mean[1:, idx])+sd, alpha=0.1)
        ax.plot(self.idx2m[1:], np.sqrt(self.pump_mean[1:, idx]),
                label='Ekman pumping')

        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlabel('m')
        ax.set_ylabel('Forces')
        ax.set_xlim(1, self.idx2m[-1])
        ax.legend(loc='best', frameon=False)
        fig.tight_layout()

    def plot(self, levels=15, cm='magma', cut=1., solid_contour=True,
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

        vmax = cut*np.log10(self.iner_mean[1:, 1:-1]).max()
        vmin = cut*np.log10(self.iner_mean[1:, 1:-1]).min()
        levs = np.linspace(vmin, vmax, levels)

        fig = plt.figure(figsize=(15, 4))
        ax1 = fig.add_subplot(151)
        im = ax1.contourf(self.radius[1:-1], self.idx2m[1:],
                          np.log10(self.buo_mean[1:, 1:-1]),
                          levs, extend='both', cmap=plt.get_cmap(cm))
        if solid_contour:
            ax1.contour(self.radius[1:-1], self.idx2m[1:],
                        np.log10(self.buo_mean[1:, 1:-1]),
                        levs, extend='both', linestyles=['-'], colors=['k'],
                        linewidths=[0.5])
        ax1.set_title('Buoyancy')
        if log_yscale:
            ax1.set_yscale('log')
        ax1.set_xlabel('Radius')
        ax1.set_ylabel('m')

        ax2 = fig.add_subplot(152, sharey=ax1, sharex=ax1)
        im = ax2.contourf(self.radius[1:-1], self.idx2m[1:],
                          np.log10(self.cor_mean[1:, 1:-1]),
                          levs, extend='both', cmap=plt.get_cmap(cm))
        if solid_contour:
            ax2.contour(self.radius[1:-1], self.idx2m[1:],
                        np.log10(self.cor_mean[1:, 1:-1]),
                        levs, extend='both', linestyles=['-'], colors=['k'],
                        linewidths=[0.5])
        ax2.set_title('Coriolis')
        if log_yscale:
            ax2.set_yscale('log')
        ax2.set_xlabel('Radius')

        ax3 = fig.add_subplot(153, sharey=ax1, sharex=ax1)
        im = ax3.contourf(self.radius, self.idx2m[1:],
                          np.log10(self.iner_mean[1:, :]),
                          levs, extend='both', cmap=plt.get_cmap(cm))
        if solid_contour:
            ax3.contour(self.radius, self.idx2m[1:],
                        np.log10(self.iner_mean[1:, :]),
                        levs, extend='both', linestyles=['-'], colors=['k'],
                        linewidths=[0.5])
        ax3.set_title('Inertia')
        if log_yscale:
            ax3.set_yscale('log')
        ax3.set_xlabel('Radius')

        ax4 = fig.add_subplot(154, sharey=ax1, sharex=ax1)
        im = ax4.contourf(self.radius, self.idx2m[1:],
                          np.log10(self.visc_mean[1:, :]),
                          levs, extend='both', cmap=plt.get_cmap(cm))
        if solid_contour:
            ax4.contour(self.radius, self.idx2m[1:],
                        np.log10(self.visc_mean[1:, :]),
                        levs, extend='both', linestyles=['-'], colors=['k'],
                        linewidths=[0.5])
        ax4.set_title('Viscosity')
        if log_yscale:
            ax4.set_yscale('log')
        ax4.set_xlabel('Radius')

        ax5 = fig.add_subplot(155, sharey=ax1, sharex=ax1)
        im = ax5.contourf(self.radius[1:-1], self.idx2m[1:],
                          np.log10(self.pump_mean[1:, 1:-1]),
                          levs, extend='both', cmap=plt.get_cmap(cm))
        if solid_contour:
            ax5.contour(self.radius[1:-1], self.idx2m[1:],
                        np.log10(self.pump_mean[1:, 1:-1]),
                        levs, extend='both', linestyles=['-'], colors=['k'],
                        linewidths=[0.5])
        ax5.set_title('Ekman pumping')
        if log_yscale:
            ax5.set_yscale('log')
        ax5.set_xlabel('Radius')
        fig.colorbar(im)

        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.setp(ax3.get_yticklabels(), visible=False)
        plt.setp(ax4.get_yticklabels(), visible=False)
        plt.setp(ax5.get_yticklabels(), visible=False)

        ax1.set_xlim(self.radius[-1], self.radius[0])

        fig.tight_layout()


if __name__ == '__main__':
    r = PizzaVortBalance()
    plt.show()
