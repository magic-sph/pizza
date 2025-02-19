# -*- coding: utf-8 -*-
import os
import numpy as np
from .plotlib import equatContour
from .libpizza import spec_spat, scanDir, symmetrize
from .frame import Frame
import matplotlib.pyplot as plt


class PizzaMovie:

    def __init__(self, field='up', nvar='all', datadir='.', tag=None, step=1,
                 endian='l', levels=64, cm=None, deminc=True, png=False,
                 cut=1., bgcolor=None, dpi=100, normed=True,
                 rm_hor_avg=False, pcolor=False):

        if png:
            plt.ioff()
            if not os.path.exists('movie'):
                os.mkdir('movie')
        else:
            plt.ion()

        if field in ('om', 'vortz', 'vort', 'omega', 'Vorticity', 'Omega'):
            st = 'om'
            if cm is None:
                try:
                    import cmocean.cm as cmo
                    cm = cmo.curl
                except ModuleNotFoundError:
                    cm = 'Spectral_r'
        elif field in ('temperature', 'Temperature', 'temp', 'Temp', 't', 'T'):
            st = 'temp'
            if cm is None:
                try:
                    import cmocean.cm as cmo
                    cm = cmo.thermal
                except ModuleNotFoundError:
                    cm = 'magma'
            normed = False
        elif field in ('comp', 'Comp', 'composition', 'xi', 'Xi'):
            st = 'xi'
            if cm is None:
                try:
                    import cmocean.cm as cmo
                    cm = cmo.haline
                except ModuleNotFoundError:
                    cm = 'viridis'
            normed = False
        elif field in ('us', 'Us', 'ur', 'Ur', 'vs', 'Vs', 'Vr', 'vr'):
            st = 'us'
            if cm is None:
                cm = 'seismic'
        elif field in ('up', 'Up', 'uphi', 'Uphi', 'vp', 'Vp', 'Vphi', 'vphi'):
            st = 'up'
            if cm is None:
                cm = 'seismic'

        # Determine the file name pattern
        if tag is not None:
            pattern = 'frame_{}_*.{}'.format(st, tag)
        else:
            pattern = 'frame_{}_*'.format(st)
        pattern = os.path.join(datadir, pattern)

        # Assemble a list of files that match the pattern
        files = scanDir(pattern)

        # Only every steps
        files = files[::step]

        # Only the last nvar files
        if type(nvar) == int:
            if nvar <= len(files):
                files = files[-nvar:]

        for k, file in enumerate(files):
            if png:
                filename = 'movie/img{:05d}.png'.format(k)
                if os.path.exists(filename) and k != 0:
                    continue
            f = Frame(file, endian=endian)
            if st == 'temp':
                f.field_m[0, :] += f.tcond
            elif st == 'xi':
                f.field_m[0, :] += f.xicond
            dat = spec_spat(f.field_m, f.n_phi_max)
            if deminc:
                dat = symmetrize(dat, ms=f.minc)

            if rm_hor_avg:
                dat = dat - np.mean(dat, axis=0)

            if k == 0:
                self.ra = f.ra
                self.ek = f.ek
                self.pr = f.pr
                self.radratio = f.radratio
                self.sc = f.sc
                self.raxi = f.raxi
                self.n_r_max = f.n_r_max
                self.n_m_max = f.n_m_max
                self.m_max = f.m_max
                self.n_phi_max = f.n_phi_max
                self.minc = f.minc
                self.radius = f.radius
                self.tcond = f.tcond
                self.idx2m = f.idx2m

                if deminc:
                    phi = np.linspace(0., 2.*np.pi, dat.shape[0])
                else:
                    phi = np.linspace(0., 2.*np.pi/self.minc, dat.shape[0])

                if normed:
                    vmin = -max(abs(dat.max()), abs(dat.min()))
                    vmin = cut * vmin
                    vmax = -vmin
                else:
                    vmin = dat.min()
                    vmax = dat.max()

                print(vmin, vmax)
                cs = np.linspace(vmin, vmax, levels)

                fig, xx, yy = equatContour(dat, self.radius, minc=self.minc,
                                           levels=levels, cm=cm, vmin=vmin,
                                           vmax=vmax, deminc=deminc,
                                           cbar=False, pcolor=pcolor)

                man = plt.get_current_fig_manager()
                man.canvas.draw()
            else:
                if not png:
                    print(k)

                plt.cla()
                if normed:
                    vmin = -max(abs(dat.max()), abs(dat.min()))
                    vmin = cut * vmin
                    vmax = -vmin
                    cs = np.linspace(vmin, vmax, levels)

                ax = fig.axes[0]
                if pcolor:
                    ax.pcolormesh(xx, yy, dat, cmap=plt.get_cmap(cm),
                                  vmin=vmin, vmax=vmax, shading='gouraud',
                                  antialiased=True, rasterized=True)
                else:
                    ax.contourf(xx, yy, dat, cs, cmap=plt.get_cmap(cm),
                                extend='both')

                ax.plot(self.radius[0]*np.cos(phi), self.radius[0]*np.sin(phi),
                        'k-')
                ax.plot(self.radius[-1]*np.cos(phi),
                        self.radius[-1]*np.sin(phi), 'k-')
                if not deminc and self.minc > 1:
                    ax.plot(self.radius, np.zeros_like(self.radius), 'k-',
                            lw=1.5)
                    xa = self.radius[-1]*np.cos(2.*np.pi/self.minc)
                    ya = self.radius[-1]*np.sin(2.*np.pi/self.minc)
                    xb = self.radius[0]*np.cos(2.*np.pi/self.minc)
                    x = np.linspace(xa, xb, 32)
                    y = np.tan(2.*np.pi/self.minc)*(x-xa)+ya
                    ax.plot(x, y, 'k-', lw=1.5)
                    ax.plot(self.radius, np.zeros_like(self.radius), 'k-',
                            lw=1.5)
                ax.axis('off')

                if xx.min() < 0:
                    ax.set_xlim(1.01*xx.min(), 1.01*xx.max())
                elif xx.min() == 0.:
                    ax.set_xlim(xx.min()-0.01, 1.01*xx.max())
                else:
                    ax.set_xlim(0.99*xx.min(), 1.01*xx.max())
                if yy.min() < 0:
                    ax.set_ylim(1.01*yy.min(), 1.01*yy.max())
                elif yy.min() == 0.:
                    ax.set_ylim(yy.min()-0.01, 1.01*yy.max())
                else:
                    ax.set_ylim(0.99*yy.min(), 1.01*yy.max())

                man.canvas.draw()

            if png:
                filename = 'movie/img{:05d}.png'.format(k)
                print('write {}'.format(filename))
                # st = 'echo {}'.format(ivar) + ' > movie/imgmax'
                if bgcolor is not None:
                    fig.savefig(filename, facecolor=bgcolor, dpi=dpi)
                else:
                    fig.savefig(filename, dpi=dpi)
