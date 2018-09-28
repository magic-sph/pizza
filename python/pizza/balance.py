from .npfile import *
import numpy as np
import matplotlib.pyplot as plt
import copy
from .log import PizzaSetup
from .libpizza import scanDir, costf
import os, re


def interp_dct(arr, nr_new):
    """
    This routine interpolates an input array of size (n_snap, n_r_max)
    onto a new grid of size (n_snap, n_r_max_new) using Discrete Cosine
    Transforms
    """

    nr = arr.shape[-1]
    # Transform to Chebyshev space and renormalise for the new grid
    fhat = costf(arr)/np.sqrt((nr-1.)/(nr_new-1.))
    fout = np.zeros((arr.shape[0], nr_new), 'Float64')
    # Set the first nr Cheb coefficients and then pad with zeros
    fout[:, :nr] = fhat
    # Finally bring back the array on the physical grid
    fout = costf(fout)

    return fout



class PizzaBalance(PizzaSetup):
    """
    This module is used to load and display the output files that
    contain the azimuthal force balance vphi_bal.TAG
    """

    def __init__(self, datadir='.', tag=None, endian='l', iplot=False, all=False):
        """
        :param datadir: working directory
        :type datadir: str
        :param tag: a specific trailing tag, default is None
        :type tag: str
        :param endian: endianness of the file ('B' or 'l')
        :type endian: str
        :param iplot: boolean to toggle the display (False by default)
        :type iplot: bool
        :param all: when set to True, the complete time series is reconstructed by
                    stacking all the corresponding files from the working directory
                    (False by default)
        :type all: bool
        """
        pattern = os.path.join(datadir, 'log.*')
        logFiles = scanDir(pattern)

        if tag is not None:
            pattern = os.path.join(datadir, 'vphi_bal.%s' % tag)
            files = scanDir(pattern)

            # Either the log.tag directly exists and the setup is easy to obtain
            if os.path.exists(os.path.join(datadir, 'log.%s' % tag)):
                PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % tag)
            # Or the tag is a bit more complicated and we need to find 
            # the corresponding log file
            else:
                mask = re.compile(r'vphi_bal\.(.*)')
                if mask.match(files[-1]):
                    ending = mask.search(files[-1]).groups(0)[0]
                    if logFiles.__contains__('log.%s' % ending):
                        PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                            nml='log.%s' % ending)

            for k, file in enumerate(files):
                print('reading %s' % file)
                if k == 0:
                    self.radius, self.time, self.vp, self.dvpdt, self.rey_stress, \
                               self.ek_pump, self.visc = self.read(file, endian)
                else:
                    radius, time, vp, dvpdt, rey_stress, ek_pump, visc = \
                                              self.read(file, endian)
                    self.add(time, radius, vp, dvpdt, rey_stress, ek_pump, visc)

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
                mot = 'vphi_bal.*'
                dat = [(os.stat(i).st_mtime, i) for i in glob.glob(mot)]
                dat.sort()
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
                    self.radius, self.time, self.vp, self.dvpdt, self.rey_stress, \
                                self.ek_pump, self.visc = self.read(file, endian)
                else:
                    radius, time, vp, dvpdt, rey_stress, ek_pump, visc = \
                                              self.read(file, endian)
                    self.add(time, radius, vp, dvpdt, rey_stress, ek_pump, visc)

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

            self.radius = radius

        if new.time[0] == self.time[-1]:
            out.time = np.concatenate((self.time, new.time[1:]), axis=0)
            out.vp = np.concatenate((self.vp, new.vp[1:, :]), axis=0)
            out.dvpdt = np.concatenate((self.dvpdt, new.dvpdt[1:, :]), axis=0)
            out.rey_stress = np.concatenate((self.rey_stress, new.rey_stress[1:, :]), axis=0)
            out.ek_pump = np.concatenate((self.ek_pump, new.ek_pump[1:, :]), axis=0)
            out.visc = np.concatenate((self.visc, new.visc[1:, :]), axis=0)
        else:
            out.time = np.concatenate((self.time, new.time), axis=0)
            out.vp = np.concatenate((self.vp, new.vp), axis=0)
            out.dvpdt = np.concatenate((self.dvpdt, new.dvpdt), axis=0)
            out.rey_stress = np.concatenate((self.rey_stress, new.rey_stress), axis=0)
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
            self.rey_stress = np.concatenate((self.rey_stress, rey_stress[1:, :]), axis=0)
            self.ek_pump = np.concatenate((self.ek_pump, ek_pump[1:, :]), axis=0)
            self.visc = np.concatenate((self.visc, visc[1:, :]), axis=0)
        else:
            self.time = np.concatenate((self.time, time), axis=0)
            self.vp = np.concatenate((self.vp, vp), axis=0)
            self.dvpdt = np.concatenate((self.dvpdt, dvpdt), axis=0)
            self.rey_stress = np.concatenate((self.rey_stress, rey_stress), axis=0)
            self.ek_pump = np.concatenate((self.ek_pump, ek_pump), axis=0)
            self.visc = np.concatenate((self.visc, visc), axis=0)

    def read(self, filename, endian):
        """
        Read one vphi_bal.TAG file
        """

        file = npfile(filename, endian=endian)

        self.ra, self.ek, self.pr, self.radratio, self.raxi, self.sc \
                 = file.fort_read('Float64')
        radius = file.fort_read('Float64')
        time = file.fort_read('Float64')
        vp = file.fort_read('Float64')
        dvpdt = file.fort_read('Float64')
        rey_stress = file.fort_read('Float64')
        ek_pump = file.fort_read('Float64')
        visc = file.fort_read('Float64')
        while 1:
            try:
                time = np.append(time, file.fort_read('Float64'))
                vp = np.vstack((vp, file.fort_read('Float64')))
                dvpdt = np.vstack((dvpdt, file.fort_read('Float64')))
                rey_stress = np.vstack((rey_stress, file.fort_read('Float64')))
                ek_pump = np.vstack((ek_pump, file.fort_read('Float64')))
                visc = np.vstack((visc, file.fort_read('Float64')))
            except TypeError:
                break

        file.close()

        return radius, time, vp, dvpdt, rey_stress, ek_pump, visc

    def plot(self):
        """
        Display when ``iplot = True``
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)
        dvpdtm = self.dvpdt.mean(axis=0)
        dvpdtstd = self.dvpdt.std(axis=0)
        ax.fill_between(self.radius, dvpdtm-dvpdtstd, dvpdtm+dvpdtstd, alpha=0.1)
        ax.plot(self.radius, dvpdtm, label='dvp/dt')

        reym = self.rey_stress.mean(axis=0)
        reystd = self.rey_stress.std(axis=0)
        ax.fill_between(self.radius, reym-reystd, reym+reystd, alpha=0.1)
        ax.plot(self.radius, reym, label='Reynolds stress')

        ekpm = self.ek_pump.mean(axis=0)
        ekpstd = self.ek_pump.std(axis=0)
        ax.fill_between(self.radius, ekpm-ekpstd, ekpm+ekpstd, alpha=0.1)
        ax.plot(self.radius, ekpm, label='Ekman pumping')
        
        vim = self.visc.mean(axis=0)
        vistd = self.visc.std(axis=0)
        ax.fill_between(self.radius, vim-vistd, vim+vistd, alpha=0.1)
        ax.plot(self.radius, vim, label='Viscosity')
        tot = self.rey_stress+self.ek_pump+self.visc
        ax.plot(self.radius, tot.mean(axis=0), label='Total')

        ax.set_xlim(self.radius[-1], self.radius[0])
        ax.legend(loc='best', frameon=False)
        ax.set_xlabel('Radius')
        ax.set_ylabel('Forces')
        fig.tight_layout()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        vpm = self.vp.mean(axis=0)
        vpstd = self.vp.std(axis=0)
        vpm_p = vpm + vpstd
        vpm_m = vpm - vpstd
        #xs, ys = mlab.poly_between(self.radius, vpm_m, vpm_p)
        ax.fill_between(self.radius, vpm_m, vpm_p, alpha=0.1)
        ax.plot(self.radius, vpm)
        #ax.plot(self.radius, vm)
        ax.set_xlim(self.radius[-1], self.radius[0])
        ax.set_xlabel('time')
        ax.set_ylabel('up')
        fig.tight_layout()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        vmax = abs(self.vp).max()
        cs = np.linspace(-vmax, vmax, 65)
        ax.contourf(self.time, self.radius, self.vp.T, cs,
                    cmap=plt.get_cmap('seismic'), extend='both')

        ax.set_xlabel('Time')
        ax.set_ylabel('Radius')
        fig.tight_layout()


class PizzaVortBalance(PizzaSetup):

    def __init__(self, datadir='.', tag=None, endian='l', iplot=True, all=False):
        """
        :param datadir: working directory
        :type datadir: str
        :param tag: a specific trailing tag, default is None
        :type tag: str
        :param endian: endianness of the file ('B' or 'l')
        :type endian: str
        :param iplot: boolean to toggle the display (False by default)
        :type iplot: bool
        :param all: when set to True, the complete time series is reconstructed by
                    stacking all the corresponding files from the working directory
                    (False by default)
        :type all: bool
        """
        pattern = os.path.join(datadir, 'log.*')
        logFiles = scanDir(pattern)

        if tag is not None:
            pattern = os.path.join(datadir, 'vort_bal.%s' % tag)
            files = scanDir(pattern)

            # Either the log.tag directly exists and the setup is easy to obtain
            if os.path.exists(os.path.join(datadir, 'log.%s' % tag)):
                PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % tag)
            # Or the tag is a bit more complicated and we need to find 
            # the corresponding log file
            else:
                mask = re.compile(r'vort_bal\.(.*)')
                if mask.match(files[-1]):
                    ending = mask.search(files[-1]).groups(0)[0]
                    if logFiles.__contains__('log.%s' % ending):
                        PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                            nml='log.%s' % ending)

            filename = files[-1]
            print('reading %s' % filename)
            self.read(filename, endian)

        else:
            if len(logFiles) != 0:
                PizzaSetup.__init__(self, quiet=True, nml=logFiles[-1])
                name = 'vort_bal.%s' % self.tag
                filename = os.path.join(datadir, name)
                print('reading %s' % filename)
                self.read(filename, endian)
            else:
                mot = 'vort_bal.*'
                dat = [(os.stat(i).st_mtime, i) for i in glob.glob(mot)]
                dat.sort()
                filename = dat[-1][1]
                print('reading %s' % filename)
                self.read(filename, endian)

        if iplot == True:
            self.plot()


    def read(self, filename, endian='l'):
        """
        :param filename: name of the input file
        :type filename: str
        :param endian: endianness of the binary file
        :type endian: str
        """
        file = open(filename, 'rb')
        dt = np.dtype("i4, 6Float64")
        self.version, params = np.fromfile(file, dtype=dt, count=1)[0]
        self.ra, self.pr, self.raxi, self.sc, self.ek, self.radratio = params
        dt = np.dtype("4i4")
        self.n_r_max, self.n_m_max, self.m_max, self.minc = np.fromfile(file, dtype=dt, count=1)[0]

        dt = np.dtype("%iFloat64" % self.n_r_max)
        self.radius = np.fromfile(file, dtype=dt, count=1)[0]
        self.idx2m = np.zeros(self.n_m_max)
        for i in range(self.n_m_max):
            self.idx2m[i] = i*self.minc

        dt = np.dtype("(%i,%i)Float64" % (self.n_r_max,self.n_m_max))
        self.buo = np.fromfile(file, dtype=dt, count=1)[0, :, :]
        self.buo = self.buo.T
        self.cor = np.fromfile(file, dtype=dt, count=1)[0, :, :]
        self.cor = self.cor.T
        self.adv = np.fromfile(file, dtype=dt, count=1)[0, :, :]
        self.adv = self.adv.T
        self.domdt = np.fromfile(file, dtype=dt, count=1)[0, :, :]
        self.domdt = self.domdt.T
        self.visc = np.fromfile(file, dtype=dt, count=1)[0, :, :]
        self.visc = self.visc.T
        self.pump = np.fromfile(file, dtype=dt, count=1)[0, :, :]
        self.pump = self.pump.T
        self.thwind = np.fromfile(file, dtype=dt, count=1)[0, :, :]
        self.thwind = self.thwind.T
        self.iner = np.fromfile(file, dtype=dt, count=1)[0, :, :]
        self.iner = self.iner.T
        self.cia = np.fromfile(file, dtype=dt, count=1)[0, :, :]
        self.cia = self.cia.T

        file.close()

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
        :param solid_contour: a boolean to decide whether one also wants the solid
                              contour lines
        :type solid_contour: bool
        :param log_yscale: a boolean to decide whether one wants a logarithmic
                           y-axis
        :type log_yscale: bool
        """

        vmax = cut*np.log10(self.iner[1:,1:-1]).max()
        vmin = cut*np.log10(self.iner[1:,1:-1]).min()
        levs = np.linspace(vmin, vmax, levels)

        fig = plt.figure(figsize=(15,4))
        ax1 = fig.add_subplot(151)
        im = ax1.contourf(self.radius[1:-1], self.idx2m[1:], np.log10(self.buo[1:,1:-1]),
                         levs, extend='both', cmap=plt.get_cmap(cm))
        if solid_contour:
            ax1.contour(self.radius[1:-1], self.idx2m[1:], np.log10(self.buo[1:,1:-1]),
                       levs, extend='both', linestyles=['-'], colors=['k'],
                       linewidths=[0.5])
        ax1.set_title('Buoyancy')
        if log_yscale: ax1.set_yscale('log')
        ax1.set_xlabel('Radius')
        ax1.set_ylabel('m')
        
        ax2 = fig.add_subplot(152, sharey=ax1, sharex=ax1)
        im = ax2.contourf(self.radius[1:-1], self.idx2m[1:], np.log10(self.cor[1:,1:-1]),
                         levs, extend='both', cmap=plt.get_cmap(cm))
        if solid_contour:
            ax2.contour(self.radius[1:-1], self.idx2m[1:], np.log10(self.cor[1:,1:-1]),
                       levs, extend='both', linestyles=['-'], colors=['k'],
                       linewidths=[0.5])
        ax2.set_title('Coriolis')
        if log_yscale: ax2.set_yscale('log')
        ax2.set_xlabel('Radius')

        ax3 = fig.add_subplot(153, sharey=ax1, sharex=ax1)
        im = ax3.contourf(self.radius, self.idx2m[1:], np.log10(self.iner[1:,:]),
                         levs, extend='both', cmap=plt.get_cmap(cm))
        if solid_contour:
            ax3.contour(self.radius, self.idx2m[1:], np.log10(self.iner[1:,:]),
                       levs, extend='both', linestyles=['-'], colors=['k'],
                       linewidths=[0.5])
        ax3.set_title('Inertia')
        if log_yscale: ax3.set_yscale('log')
        ax3.set_xlabel('Radius')

        ax4 = fig.add_subplot(154, sharey=ax1, sharex=ax1)
        im = ax4.contourf(self.radius, self.idx2m[1:], np.log10(self.visc[1:,:]),
                         levs, extend='both', cmap=plt.get_cmap(cm))
        if solid_contour:
            ax4.contour(self.radius, self.idx2m[1:], np.log10(self.visc[1:,:]),
                       levs, extend='both', linestyles=['-'], colors=['k'],
                       linewidths=[0.5])
        ax4.set_title('Viscosity')
        if log_yscale: ax4.set_yscale('log')
        ax4.set_xlabel('Radius')

        ax5 = fig.add_subplot(155, sharey=ax1, sharex=ax1)
        im = ax5.contourf(self.radius[1:-1], self.idx2m[1:], np.log10(self.pump[1:,1:-1]),
                         levs, extend='both', cmap=plt.get_cmap(cm))
        if solid_contour:
            ax5.contour(self.radius[1:-1], self.idx2m[1:], np.log10(self.pump[1:,1:-1]),
                       levs, extend='both', linestyles=['-'], colors=['k'],
                       linewidths=[0.5])
        ax5.set_title('Ekman pumping')
        if log_yscale: ax5.set_yscale('log')
        ax5.set_xlabel('Radius')
        fig.colorbar(im)

        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.setp(ax3.get_yticklabels(), visible=False)
        plt.setp(ax4.get_yticklabels(), visible=False)
        plt.setp(ax5.get_yticklabels(), visible=False)

        ax1.set_xlim(self.radius[-1], self.radius[0])

        fig.tight_layout()


if __name__ == '__main__':
   
   r = PizzaVortBal()
   plt.show()
