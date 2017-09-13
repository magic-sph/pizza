import os
from .npfile import *
import numpy as np
from .plotlib import equatContour
from .libpizza import spec_spat, symmetrize


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
              normed=True, vmax=None, vmin=None, normRad=False):

        if field in ('om', 'vortz', 'vort', 'omega', 'Vorticity', 'Omega'):
            data = self.vortz
        elif field in ('temperature', 'Temperature', 'temp', 'Temp', 't', 'T'):
            data = self.temp
        elif field in ('us', 'Us', 'ur', 'Ur', 'vs', 'Vs', 'Vr', 'vr'):
            data = self.us
        elif field in ('up', 'Up', 'uphi', 'Uphi', 'vp', 'Vp', 'Vphi', 'vphi'):
            data = self.uphi

        if deminc:
            data = symmetrize(data, ms=self.minc)

        equatContour(data, self.radius, minc=self.minc, levels=levels,
                     cm=cm, deminc=deminc, normed=normed, vmax=vmax, vmin=vmin,
                     normRad=normRad)
