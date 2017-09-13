from .npfile import *
import numpy as np
import matplotlib.pyplot as plt
from .log import PizzaSetup


class PizzaBalance:

    def __init__(self, tag='None', endian='B', iplot=False):
        filename = 'vphi_bal.%s' % tag
        file = npfile(filename, endian=endian)

        self.ra, self.ek, self.pr, self.radratio, self.raxi, self.sc \
                 = file.fort_read('Float64')
        self.radius = file.fort_read('Float64')
        self.time = file.fort_read('Float64')
        self.vp = file.fort_read('Float64')
        self.dvpdt = file.fort_read('Float64')
        self.rey_stress = file.fort_read('Float64')
        self.ek_pump = file.fort_read('Float64')
        self.visc = file.fort_read('Float64')
        while 1:
            try:
                self.time = np.append(self.time, file.fort_read('Float64'))
                self.vp = np.vstack((self.vp, file.fort_read('Float64')))
                self.dvpdt = np.vstack((self.dvpdt, file.fort_read('Float64')))
                self.rey_stress = np.vstack((self.rey_stress, file.fort_read('Float64')))
                self.ek_pump = np.vstack((self.ek_pump, file.fort_read('Float64')))
                self.visc = np.vstack((self.visc, file.fort_read('Float64')))
            except TypeError:
                break

        file.close()

        if iplot:
            self.plot()

    def plot(self):
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

        fig = plt.figure()
        ax = fig.add_subplot(111)
        vmax = abs(self.vp).max()
        cs = np.linspace(-vmax, vmax, 65)
        ax.contourf(self.time, self.radius, self.vp.T, cs,
                    cmap=plt.get_cmap('seismic'), extend='both')

        ax.set_xlabel('Time')
        ax.set_ylabel('Radius')
