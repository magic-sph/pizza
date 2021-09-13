# -*- coding: utf-8 -*-
import os
import re
from .npfile import npfile
import numpy as np
from .log import PizzaSetup
from .plotlib import equatContour
from .libpizza import spec_spat, symmetrize, scanDir
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
    """
    This module is used to read the 2-D snapshots
    """

    def __init__(self, filename, endian='l'):
        """
        :param filename: name of the file
        :type filename: str
        :param endian: endianness of the file ('B' or 'l')
        :type endian: str
        """
        try:

            file = npfile(filename, endian=endian)
            self.version = file.fort_read('i4')
            self.time = file.fort_read(np.float64)
            self.ra, self.ek, self.pr, self.radratio, self.sc, \
                self.raxi = file.fort_read(np.float64)
            self.n_r_max, self.n_m_max, self.m_max, self.minc, \
                self.n_phi_max = file.fort_read('i4')

            self.radius = file.fort_read(np.float64)
            self.tcond = file.fort_read(np.float64)
            self.idx2m = file.fort_read('i4')

            self.field_m = file.fort_read(np.complex128)
            self.field_m = self.field_m.reshape((self.n_r_max, self.n_m_max))
            self.field_m = self.field_m.T
            file.close()

        except:

            file = open(filename, 'rb')
            self.version = np.fromfile(file, dtype='i4', count=1)[0]
            self.time, self.ra, self.ek, self.pr, self.radratio, self.sc, \
                self.raxi = np.fromfile(file, dtype='7f8', count=1)[0]
            self.n_r_max, self.n_m_max, self.m_max, self.minc, \
                self.n_phi_max = np.fromfile(file, dtype='5i4', count=1)[0]

            self.radius = np.fromfile(file, dtype='%if8' % self.n_r_max,
                                      count=1)[0]
            self.tcond = np.fromfile(file, dtype='%if8' % self.n_r_max,
                                     count=1)[0]
            if self.version == 2:
                self.xicond = np.fromfile(file,
                                          dtype='%if8' % self.n_r_max,
                                          count=1)[0]

            self.idx2m = np.fromfile(file, dtype='%ii4' % self.n_m_max,
                                     count=1)[0]

            dt = np.dtype("(%i,%i)c16" % (self.n_r_max, self.n_m_max))
            self.field_m = np.fromfile(file, dtype=dt, count=1)[0]
            self.field_m = self.field_m.T
            file.close()

    def write(self, filename):
        """
        This routine writes a snap using stream

        :param filename: the name of the output file
        :type filename: str
        """

        out = open('%s' % filename, 'wb')
        x = np.array([1], dtype="i4")
        x.tofile(out)
        x = np.array([self.time, self.ra, self.ek, self.pr, self.radratio,
                      self.sc, self.raxi], dtype=np.float64)
        x.tofile(out)
        x = np.array([self.n_r_max, self.n_m_max, self.m_max, self.minc,
                      self.n_phi_max], dtype='i4')
        x.tofile(out)

        self.radius.tofile(out)
        self.tcond.tofile(out)
        self.idx2m.tofile(out)
        self.field_m.T.tofile(out)

        out.close()


class PizzaFields(PizzaSetup):
    """
    This module is used to read the frames, transform them back on the
    physical grid and display them.
    """

    def __init__(self, ivar=None, datadir='.', tag=None, endian='l',
                 verbose=False):
        """
        :param ivar: the number of the snapshot file
        :type ivar: int
        :param datadir: working directory
        :type datadir: str
        :param tag: extension TAG of the snapshot files
        :type tag: str
        :param endian: endianness of the file ('B' or 'l')
        :type endian: str
        :param verbose: a boolean to display some informations
        :type verbose: bool
        """

        filename = self.get_filename('frame_us', ivar, datadir, tag, verbose)
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
        if f.version == 2:
            self.xicond = f.xicond
        else:
            self.xicond = np.zeros_like(self.tcond)

        self.idx2m = f.idx2m
        self.time = f.time

        self.us_m = f.field_m
        self.us = spec_spat(self.us_m, self.n_phi_max)

        filename = self.get_filename('frame_up', ivar, datadir, tag, verbose)
        f = Frame(filename, endian=endian)
        self.uphi_m = f.field_m
        self.uphi = spec_spat(self.uphi_m, self.n_phi_max)

        filename = self.get_filename('frame_om', ivar, datadir, tag, verbose)
        f = Frame(filename, endian=endian)
        self.vortz_m = f.field_m
        self.vortz = spec_spat(self.vortz_m, self.n_phi_max)

        filename = self.get_filename('frame_temp', ivar, datadir, tag, verbose)
        if filename is not None:
            f = Frame(filename, endian=endian)
            self.temp_m = f.field_m
            self.temp_m[0, :] += self.tcond
            self.temp = spec_spat(self.temp_m, self.n_phi_max)

        filename = self.get_filename('frame_xi', ivar, datadir, tag, verbose)
        if filename is not None:
            f = Frame(filename, endian=endian)
            self.xi_m = f.field_m
            self.xi_m[0, :] += self.xicond
            self.xi = spec_spat(self.xi_m, self.n_phi_max)

    def get_filename(self, prefix, ivar, datadir, tag, verbose):
        """
        This routine determines the filename based on what is available
        in the current directory

        :param prefix: the file prefix ('frame_temp', 'frame_us',
                       'frame_om', ...)
        :type prefix: str
        :param ivar: the number of the snapshot file
        :type ivar: int
        :param datadir: working directory
        :type datadir: str
        :param tag: extension TAG of the snapshot files
        :type tag: str
        :param verbose: a boolean to display some informations
        :type verbose: bool
        """

        if tag is not None:
            if ivar is not None:
                file = '%s_%i.%s' % (prefix, ivar, tag)
                filename = os.path.join(datadir, file)
            else:
                files = scanDir('%s_*%s' % (prefix, tag))
                if len(files) != 0:
                    filename = os.path.join(datadir, files[-1])
                else:
                    return

            if os.path.exists('log.%s' % tag):
                PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % tag)
        else:
            if ivar is not None:
                files = scanDir('%s_%i.*' % (prefix, ivar))
                if len(files) != 0:
                    filename = os.path.join(datadir, files[-1])
                else:
                    return
            else:
                files = scanDir('%s_*' % prefix)
                if len(files) != 0:
                    filename = os.path.join(datadir, files[-1])
                else:
                    return
            # Determine the setup
            mask = re.compile(r'.*\.(.*)')
            ending = mask.search(files[-1]).groups(0)[0]
            if os.path.exists('log.%s' % ending):
                PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % ending)

        if not os.path.exists(filename):
            return

        if verbose:
            print('read %s' % filename)
        return filename

    def equat(self, field='vort', cm='seismic', levels=65, deminc=True,
              normed=True, vmax=None, vmin=None, normRad=False, stream=False,
              streamNorm='vel', streamDensity=1.5, cbar=True, label=None,
              streamColor='k'):
        """
        Display an equatorial planform of a scalar quantity

        :param field: the name of the input physical quantity you want to
                      display
        :type field: str
        :param normRad: when set to True, the contour levels are normalised
                        radius by radius (default is False)
        :type normRad: bool
        :param levels: the number of levels in the contour
        :type levels: int
        :param cm: name of the colormap ('jet', 'seismic', 'RdYlBu_r', etc.)
        :type cm: str
        :param tit: display the title of the figure when set to True
        :type tit: bool
        :param cbar: display the colorbar when set to True
        :type cbar: bool
        :param vmax: maximum value of the contour levels
        :type vmax: float
        :param vmin: minimum value of the contour levels
        :type vmin: float
        :param normed: when set to True, the colormap is centered around zero.
                       Default is True, except for entropy/temperature plots.
        :type normed: bool
        :param stream: a boolean to control the plotting of streamlines (False
                       by default)
        :type stream: bool
        :param streamNorm: scalar that is used to normalise the thickness of
                           the streamlines (by default 'vel')
        :type streamNorm: str
        :param streamDensity: control the number of streamlines
        :type streamDensity: float
        :param streamColor: color of the streamlines
        :type streamColor: str
        """

        if field in ('om', 'vortz', 'vort', 'omega', 'Vorticity', 'Omega'):
            data = self.vortz
        elif field in ('temperature', 'Temperature', 'temp', 'Temp', 't', 'T'):
            data = self.temp
        elif field in ('composition', 'Composition', 'xi', 'Xi', 'chem',
                       'Chem', 'comp', 'Comp'):
            data = self.xi
        elif field in ('tfluct', 'tempfluct'):
            data = self.temp-self.temp_m[0, :]
        elif field in ('xifluct', 'chemfluct', 'compfluct'):
            data = self.xi-self.xi_m[0, :]
        elif field in ('us', 'Us', 'ur', 'Ur', 'vs', 'Vs', 'Vr', 'vr'):
            data = self.us
        elif field in ('up', 'Up', 'uphi', 'Uphi', 'vp', 'Vp', 'Vphi', 'vphi'):
            data = self.uphi

        if deminc:
            data = symmetrize(data, ms=self.minc)

        self.fig, xx, yy = equatContour(data, self.radius, minc=self.minc,
                                        levels=levels, cm=cm, deminc=deminc,
                                        normed=normed, vmax=vmax, vmin=vmin,
                                        normRad=normRad, cbar=cbar,
                                        label=label)

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
                # theta = np.linspace(-np.pi, np.pi, data.shape[0])
                theta = np.linspace(0., 2.*np.pi, data.shape[0])
            else:
                theta = np.linspace(0., 2*np.pi/self.minc, data.shape[0])

            rad = np.linspace(self.radius[0], self.radius[-1], data.shape[1])
            rr, ttheta = np.meshgrid(rad, theta)
            if deminc:
                u = symmetrize(self.us, self.minc)
                v = symmetrize(self.uphi, self.minc)
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
                           linewidth=lw, color=streamColor)
            ax1.set_ylim(0., rad.max())
