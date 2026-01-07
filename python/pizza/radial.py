# -*- coding: utf-8 -*-
import copy
import os
import re
import matplotlib.pyplot as plt
import numpy as np
from .log import PizzaSetup
from .libpizza import fast_read, scanDir, get_dr


def plotSpan(ax, idx, s):
    for i in range(len(idx)):
        if i > 0:
            if i % 2:
                ax.axvspan(s[idx[i-1]], s[idx[i]], facecolor='0.75',
                           edgecolor='None')
            else:
                ax.axvspan(s[idx[i-1]], s[idx[i]], facecolor='0.97',
                           edgecolor='None')


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
        :param all: if all=True, then all the radial profiles from the
                    directory are stacked together and averaged by their
                    respective time span
        :type all: bool
        """

        name = 'radial_profiles'

        if not all:
            if tag is not None:
                pattern = os.path.join(datadir, f'{name}.{tag}')
                files = scanDir(pattern)

                # Either the log.tag directly exists and the setup is easy to
                # obtain
                if os.path.exists(os.path.join(datadir, f'log.{tag}')):
                    PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                        nml=f'log.{tag}')
                # Or the tag is a bit more complicated and we need to find
                # the corresponding log file
                else:
                    mask = re.compile(fr'{datadir}\/{name}\.(.*)')
                    if mask.match(files[-1]):
                        ending = mask.search(files[-1]).groups(0)[0]
                        pattern = os.path.join(datadir, f'log.{ending}')
                        if os.path.exists(pattern):
                            PizzaSetup.__init__(self, datadir=datadir,
                                                quiet=True, nml=f'log.{ending}')

                # Sum the files that correspond to the tag
                mask = re.compile(fr'{name}\.(.*)')
                for k, file in enumerate(files):
                    print(f'reading {file}')
                    tag = mask.search(file).groups(0)[0]
                    nml = PizzaSetup(nml=f'log.{tag}', datadir=datadir,
                                     quiet=True)
                    filename = file
                    if k == 0:
                        self.tstart = nml.start_time
                        self.tstop = nml.stop_time  # will be overwritten
                        data = fast_read(filename)
                    else:
                        if os.path.exists(filename):
                            tmp = fast_read(filename)
                            data = self.add(data, tmp, nml.stop_time,
                                            nml.start_time)

            else:  # if all
                pattern = os.path.join(datadir, f'{name}.*')
                files = scanDir(pattern)
                filename = files[-1]
                print(f'reading {filename}')
                # Determine the setup
                mask = re.compile(fr'{name}\.(.*)')
                ending = mask.search(files[-1]).groups(0)[0]
                if os.path.exists(os.path.join(datadir, f'log.{ending}')):
                    try:
                        PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                            nml=f'log.{ending}')
                    except AttributeError:
                        pass

                data = fast_read(filename, skiplines=0)
        else:
            pattern = os.path.join(datadir, f'{name}.*')
            files = scanDir(pattern)

            # Determine the setup
            mask = re.compile(fr'{name}\.(.*)')
            for k, file in enumerate(files):
                print(f'reading {file}')
                tag = mask.search(file).groups(0)[0]
                nml = PizzaSetup(nml=f'log.{tag}', datadir=datadir,
                                 quiet=True)
                filename = file
                if k == 0:
                    self.tstart = nml.start_time
                    self.tstop = nml.stop_time  # will be overwritten later
                    data = fast_read(filename)
                else:
                    if os.path.exists(filename):
                        tmp = fast_read(filename)
                        data = self.add(data, tmp, nml.stop_time,
                                        nml.start_time)
            PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                nml=f'log.{tag}')

        self.radius = data[:, 0]
        self.us2_r = data[:, 1]
        self.us2_r_SD = data[:, 2]
        self.up2_r = data[:, 3]
        self.up2_r_SD = data[:, 4]
        self.enst_r = data[:, 5]
        self.enst_r_SD = data[:, 6]
        self.uphi_r = data[:, 7]
        self.uphi_r_SD = data[:, 8]
        self.temp_r = data[:, 9]
        self.temp_r_SD = data[:, 10]
        if data.shape[-1] == 15:
            self.xi_r = data[:, 11]
            self.xi_r_SD = data[:, 12]
            self.nushell_r = data[:, 13]
            self.nushell_r_SD = data[:, 14]
        else:
            self.nushell_r = data[:, 11]
            self.nushell_r_SD = data[:, 12]
            self.xi_r = np.zeros_like(self.temp_r)
            self.xi_r = np.zeros_like(self.temp_r)

        self.vortz_r = 1./self.radius * get_dr(self.radius*self.uphi_r)
        h = np.sqrt(self.radius[0]**2-self.radius**2)
        self.qplan = np.zeros_like(h)
        self.qpot = np.zeros_like(h)
        self.qplan[1:] = 2./self.ek/h[1:]
        self.qpot[1:] = self.vortz_r[1:]/h[1:]+self.qplan[1:]

        self.qplan[0] = self.qplan[1]
        self.qpot[0] = self.qpot[1]

        if iplot:
            self.plot()

    def add(self, data, tmp, stop_time, start_time):
        """
        Clean way to stack data
        """
        out = copy.deepcopy(data)
        out[:, 0] = tmp[:, 0]

        nr_new = len(tmp[:, 0])
        nr_old = len(data[:, 0])

        fac_old = self.tstop-self.tstart
        fac_new = stop_time-start_time
        self.tstop = stop_time
        fac_tot = self.tstop-self.tstart

        if nr_new == nr_old:  # Same grid before and after
            for j in [1, 3, 5, 7, 9, 11]:
                out[:, j] = (fac_old*data[:, j]+fac_new*tmp[:, j])/fac_tot
            for j in [2, 4, 6, 8, 10, 12]:
                out[:, j] = np.sqrt((fac_old*data[:, j]**2 +
                                     fac_new*tmp[:, j]**2) / fac_tot)
        else:
            print('Not implemented yet ...')

        return out

    def plot(self):
        """
        Display the result when ``iplot=True``
        """

        fig, ax = plt.subplots()
        ax.fill_between(self.radius, self.temp_r-self.temp_r_SD,
                        self.temp_r+self.temp_r_SD, alpha=0.1)
        ax.plot(self.radius, self.temp_r)
        ax.set_xlabel('Radius')
        ax.set_ylabel('Temperature')
        ax.set_ylim(self.temp_r.min(), self.temp_r.max())
        ax.set_xlim(self.radius[-1], self.radius[0])
        fig.tight_layout()

        if abs(self.uphi_r).max() > 0.:
            fig, ax = plt.subplots()
            ax.fill_between(self.radius, self.uphi_r-self.uphi_r_SD,
                            self.uphi_r+self.uphi_r_SD, alpha=0.1)
            ax.plot(self.radius, self.uphi_r)
            ax.set_xlabel('Radius')
            ax.set_ylabel('Uphi')
            ax.set_xlim(self.radius[-1], self.radius[0])
            fig.tight_layout()

        fig, ax = plt.subplots()
        ax.fill_between(self.radius, self.us2_r-self.us2_r_SD,
                        self.us2_r+self.us2_r_SD, alpha=0.1)
        ax.plot(self.radius, self.us2_r, label='us**2')
        ax.fill_between(self.radius, self.up2_r-self.up2_r_SD,
                        self.up2_r+self.up2_r_SD, alpha=0.1)
        ax.plot(self.radius, self.up2_r, label='up**2')
        ax.fill_between(self.radius, self.enst_r-self.enst_r_SD,
                        self.enst_r+self.enst_r_SD, alpha=0.1)
        ax.plot(self.radius, self.enst_r, label='omega**2')
        ax.legend(loc='best', frameon=False)
        ax.set_xlabel('Radius')
        ax.set_xlim(self.radius[-1], self.radius[0])
        ax.set_yscale('log')
        fig.tight_layout()

        if self.l_non_rot == 'F' or not self.l_non_rot:
            fig, ax = plt.subplots()
            ax.plot(self.radius, self.qpot)
            ax.plot(self.radius, self.qplan)

            # Find the zeroes of uphi
            idx = np.array([], dtype=np.int8)
            for i in range(len(self.radius)):
                if i > 0:
                    if self.uphi_r[i] < 0. and self.uphi_r[i-1] > 0.:
                        idx = np.append(idx, i)
                    elif self.uphi_r[i] > 0. and self.uphi_r[i-1] < 0.:
                        idx = np.append(idx, i)

            plotSpan(ax, idx, self.radius)

            ax.set_xlim(self.radius[-1], self.radius[0])
            ax.set_ylim(0.95*self.qplan[-1], self.qplan[-1]*2)
            ax.set_xlabel('Radius')
            ax.set_ylabel('Potential vorticity')
            fig.tight_layout()

            # Determine jets width and Rhines scale
            jet_widths = np.zeros((len(idx)-1), dtype=np.float64)
            lamb_uzon = np.zeros_like(jet_widths)
            lamb_us = np.zeros_like(jet_widths)
            beta = np.zeros_like(self.radius)
            coord = np.zeros_like(jet_widths)
            beta[1:] = self.radius[1:]/(self.radius[0]**2-self.radius[1:]**2)
            beta[0] = beta[1]
            for k in range(len(idx)):
                if k > 0:
                    jet_widths[k-1] = self.radius[idx[k-1]]-self.radius[idx[k]]
                    vp_loc = abs(self.uphi_r[idx[k-1]:idx[k]+1]).mean()
                    vs_loc = np.sqrt(self.us2_r[idx[k-1]:idx[k]+1].mean())
                    beta_loc = abs(beta[idx[k-1]:idx[k]+1]).mean()
                    lamb_uzon[k-1] = 2.*np.pi*np.sqrt(vp_loc*self.ek/beta_loc)
                    lamb_us[k-1] = 2.*np.pi*np.sqrt(vs_loc*self.ek/beta_loc)
                    coord[k-1] = self.radius[idx[k-1]:idx[k]].mean()

            fig, ax = plt.subplots()
            ax.plot(coord, jet_widths, marker='o', label='Jet width')
            ax.plot(coord, lamb_uzon, marker='s', label='Rhines scale (Uzon)')
            ax.plot(coord, lamb_us, marker='s', label='Rhines scale (Us)')
            plotSpan(ax, idx, self.radius)

            ax.set_xlabel('Radius')
            ax.set_ylabel('Lengthscales')
            ax.set_xlim(self.radius[-1], self.radius[0])
            ax.legend(loc='best')
            fig.tight_layout()
