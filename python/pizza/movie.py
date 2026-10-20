# -*- coding: utf-8 -*-
import os
import numpy as np
from .plotlib import equatContour
from .libpizza import spec_spat, scanDir, symmetrize
from .frame import Frame
import matplotlib.pyplot as plt
import matplotlib.animation as animation


class PizzaMovie:

    def __init__(self, field='up', nvar='all', datadir='.', tag=None, nstep=1,
                 endian='l', levels=64, cm=None, deminc=True, png=False,
                 cut=1., bgcolor=None, dpi=100, normed=True,
                 rm_hor_avg=False, pcolor=False):

        if png:
            if not os.path.exists('movie'):
                os.mkdir('movie')

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
            pattern = f'frame_{st}_*.{tag}'
        else:
            pattern = f'frame_{st}_*'
        pattern = os.path.join(datadir, pattern)

        # Assemble a list of files that match the pattern
        files = scanDir(pattern)

        # Only every steps
        files = files[::nstep]

        # Only the last nvar files
        if type(nvar) == int:
            if nvar <= len(files):
                files = files[-nvar:]
        else:
            nvar = len(files)

        f = Frame(files[0], endian=endian)
        if st == 'temp':
            f.field_m[0, :] += f.tcond
        elif st == 'xi':
            f.field_m[0, :] += f.xicond
        dat = spec_spat(f.field_m, f.n_phi_max)
        if deminc:
            dat = symmetrize(dat, ms=f.minc)
        if rm_hor_avg:
            dat = dat - np.mean(dat, axis=0)

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

        cs = np.linspace(vmin, vmax, levels)

        fig, xx, yy, im = equatContour(dat, self.radius, minc=self.minc,
                                       levels=levels, cm=cm, vmin=vmin,
                                       vmax=vmax, deminc=deminc,
                                       cbar=False, pcolor=pcolor)
        ax = fig.get_axes()[0]

        props = dict(boxstyle='round', facecolor='wheat')
        timestamp = ax.text(0.95, 0.95, f'{f.time:g}',
                            transform=ax.transAxes, ha='right', bbox=props)

        def animate_contourf(k):
            f = Frame(files[k], endian=endian)
            if st == 'temp':
                f.field_m[0, :] += f.tcond
            elif st == 'xi':
                f.field_m[0, :] += f.xicond
            dat = spec_spat(f.field_m, f.n_phi_max)
            if deminc:
                dat = symmetrize(dat, ms=f.minc)
            if rm_hor_avg:
                dat = dat - np.mean(dat, axis=0)

            for coll in ax.collections:
                coll.remove()

            im = ax.contourf(xx, yy, dat, cs, cmap=plt.get_cmap(cm), extend='both')
            timestamp.set_text(f'{f.time:g}')

            if png:
                filename = f'movie/img{k:05d}.png'
                print(f'write {filename}')
                if bgcolor is not None:
                    fig.savefig(filename, facecolor=bgcolor, dpi=dpi)
                else:
                    fig.savefig(filename, dpi=dpi)

            return im

        def animate_pcolor(k):
            f = Frame(files[k], endian=endian)
            if st == 'temp':
                f.field_m[0, :] += f.tcond
            elif st == 'xi':
                f.field_m[0, :] += f.xicond
            dat = spec_spat(f.field_m, f.n_phi_max)
            if deminc:
                dat = symmetrize(dat, ms=f.minc)
            if rm_hor_avg:
                dat = dat - np.mean(dat, axis=0)

            im.set_array(dat)
            timestamp.set_text(f'{f.time:g}')

            if png:
                filename = f'movie/img{k:05d}.png'
                print(f'write {filename}')
                if bgcolor is not None:
                    fig.savefig(filename, facecolor=bgcolor, dpi=dpi)
                else:
                    fig.savefig(filename, dpi=dpi)

            return im

        if pcolor:
            print(files)
            anim = animation.FuncAnimation(fig, animate_pcolor, frames=nvar)
        else:
            anim = animation.FuncAnimation(fig, animate_contourf, frames=nvar)
        if png:
            anim.save('movie/out.gif')
        else:
            plt.show()


