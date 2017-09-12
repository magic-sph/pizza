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
        self.dvpdt = file.fort_read('Float64')
        self.rey_stress = file.fort_read('Float64')
        self.ek_pump = file.fort_read('Float64')
        self.visc = file.fort_read('Float64')
        while 1:
            try:
                self.time = np.append(self.time, file.fort_read('Float64'))
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
        ax.plot(self.radius, self.dvpdt.mean(axis=0), label='dvp/dt')
        ax.plot(self.radius, self.rey_stress.mean(axis=0), label='Reynolds stress')
        ax.plot(self.radius, self.ek_pump.mean(axis=0), label='Ekman pumping')
        ax.plot(self.radius, self.visc.mean(axis=0), label='Viscosity')
        tot = self.rey_stress+self.ek_pump+self.visc
        ax.plot(self.radius, tot.mean(axis=0), label='Total')

        ax.legend(loc='best', frameon=False)
        ax.set_xlabel('Radius')
        ax.set_ylabel('Forces')
