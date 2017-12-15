import os
from .npfile import *
import numpy as np
from .plotlib import equatContour
from .libpizza import spec_spat, symmetrize
import scipy.interpolate as inp

def my_interp2d(f, rad, radnew):
    r = rad[::-1]
    rnew = radnew[::-1]
    fnew = np.zeros_like(f)
    for i in range(f.shape[0]):
        val = f[i, ::-1]
        tckp = inp.splrep(r, val)
        fnew[i, :] = inp.splev(rnew, tckp)
    
    return fnew[:, ::-1]


class Frame:

    def __init__(self, filename, endian='l'):

        file = npfile(filename, endian=endian)

        version = file.fort_read('i4')
        self.time = file.fort_read('Float64')
        self.ra, self.ek, self.pr, self.radratio, self.sc, \
                      self.raxi = file.fort_read('Float64')
        self.n_r_max, self.n_m_max, self.m_max, self.minc, \
                      self.n_phi_max = file.fort_read('i4')

        self.radius = file.fort_read('Float64')
        self.tcond = file.fort_read('Float64')
        self.idx2m = file.fort_read('i4')

        self.field_m = file.fort_read('Complex64')
        self.field_m = self.field_m.reshape((self.n_r_max, self.n_m_max))
        self.field_m = self.field_m.T

        file.close()

class PizzaFields:

    def __init__(self, ivar=1, datadir='.', tag=None, endian='l'):

        filename = 'frame_temp_%i.%s' % (ivar, tag)
        filename = os.path.join(datadir, filename)
        f = Frame(filename, endian=endian)
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
        self.time = f.time

        self.temp_m = f.field_m
        self.temp_m[0, :] += self.tcond
        self.temp = spec_spat(self.temp_m, self.n_phi_max)

        filename = 'frame_us_%i.%s' % (ivar, tag)
        filename = os.path.join(datadir, filename)
        f = Frame(filename, endian=endian)
        self.us_m = f.field_m
        self.us = spec_spat(self.us_m, self.n_phi_max)

        filename = 'frame_up_%i.%s' % (ivar, tag)
        filename = os.path.join(datadir, filename)
        f = Frame(filename, endian=endian)
        self.uphi_m = f.field_m
        self.uphi = spec_spat(self.uphi_m, self.n_phi_max)

        filename = 'frame_om_%i.%s' % (ivar, tag)
        filename = os.path.join(datadir, filename)
        f = Frame(filename, endian=endian)
        self.vortz_m = f.field_m
        self.vortz = spec_spat(self.vortz_m, self.n_phi_max)

    def equat(self, field='vort', cm='seismic', levels=65, deminc=True,
              normed=True, vmax=None, vmin=None, normRad=False, stream=False,
              streamNorm='vel', streamDensity=1.5, cbar=True, label=None):

        if field in ('om', 'vortz', 'vort', 'omega', 'Vorticity', 'Omega'):
            data = self.vortz
        elif field in ('temperature', 'Temperature', 'temp', 'Temp', 't', 'T'):
            data = self.temp
        elif field in ('tfluct', 'tempfluct'):
            data = self.temp-self.temp_m[0,:]
        elif field in ('us', 'Us', 'ur', 'Ur', 'vs', 'Vs', 'Vr', 'vr'):
            data = self.us
        elif field in ('up', 'Up', 'uphi', 'Uphi', 'vp', 'Vp', 'Vphi', 'vphi'):
            data = self.uphi

        if deminc:
            data = symmetrize(data, ms=self.minc)

        self.fig, xx, yy = equatContour(data, self.radius, minc=self.minc, levels=levels,
                          cm=cm, deminc=deminc, normed=normed, vmax=vmax, vmin=vmin,
                          normRad=normRad, cbar=cbar, label=label)
        if stream:
            ax = self.fig.get_axes()[0]
            bbox = ax.get_position()
            # For now matplotlib polar plots do not handle theta_min
            # theta_max (set_xlim does not work)
            # As a trick I simply shift the polar Axes outside the domain
            # when minc /= 1
            if not deminc:
                a = bbox.bounds[0]
                b = bbox.bounds[1]
                c = bbox.bounds[2]
                d = bbox.bounds[3]
                if self.minc == 2:
                    bbox.bounds = (a, b-d, c, 2*d)
                elif self.minc == 4:
                    bbox.bounds = (a-c, b-d, 2*c, 2*d)
            ax1 = self.fig.add_axes(bbox, polar=True)
            ax1.axis('off')

            if deminc:
                #theta = np.linspace(-np.pi, np.pi, data.shape[0])
                theta = np.linspace(0., 2.*np.pi, data.shape[0])
            else:
                theta = np.linspace(0., 2*np.pi/self.minc, data.shape[0])

            rad = np.linspace(self.radius[0], self.radius[-1], data.shape[1])
            rr, ttheta = np.meshgrid(rad, theta)
            if deminc:
                u = symmetrize(self.us, self.minc)
                v = symmetrize(self.uphi,self.minc)
            else:
                u = self.us
                v = self.uphi
            v /= self.radius
            
            u = my_interp2d(u, self.radius, rad)
            v = my_interp2d(v, self.radius, rad)
            speed = np.sqrt(u**2+v**2)

            if streamNorm == 'vel':
                lw = 3.*speed.T/speed.max()
            else:
                lw = 3.*abs(data.T)/abs(self.vortz).max()
            ax1.streamplot(ttheta.T, rr.T, v.T, u.T, density=streamDensity,
                           linewidth=lw, color='k')
            ax1.set_ylim(0., rad.max())
