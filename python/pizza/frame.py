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
    r = rad
    rnew = radnew
    fnew = np.zeros_like(f)
    for i in range(f.shape[0]):
        val = f[i, :]
        tckp = inp.splrep(r, val)
        fnew[i, :] = inp.splev(rnew, tckp)

    return fnew


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
            self.version = file.fort_read(np.int32)
            self.time = file.fort_read(np.float64)
            self.ra, self.ek, self.pr, self.radratio, self.sc, \
                self.raxi = file.fort_read(np.float64)
            self.n_r_max, self.n_m_max, self.m_max, self.minc, \
                self.n_phi_max = file.fort_read(np.int32)

            self.radius = file.fort_read(np.float64)
            self.tcond = file.fort_read(np.float64)
            if self.version >= 2:
                self.xicond = file.fort_read(np.float64)
            self.idx2m = file.fort_read(np.int32)

            self.field_m = file.fort_read(np.complex128)
            self.field_m = self.field_m.reshape((self.n_r_max, self.n_m_max))
            self.field_m = self.field_m.T
            file.close()

        except:

            with open(filename, 'rb') as file:
                self.version = np.fromfile(file, dtype=np.int32, count=1)[0]
                self.time, self.ra, self.ek, self.pr, self.radratio, self.sc, \
                    self.raxi = np.fromfile(file, dtype=np.float64, count=7)
                self.n_r_max, self.n_m_max, self.m_max, self.minc, \
                    self.n_phi_max = np.fromfile(file, dtype=np.int32, count=5)

                self.radius = np.fromfile(file, dtype=np.float64,
                                          count=self.n_r_max)
                self.tcond = np.fromfile(file, dtype=np.float64,
                                         count=self.n_r_max)
                if self.version == 2:
                    self.xicond = np.fromfile(file, dtype=np.float64,
                                              count=self.n_r_max)

                self.idx2m = np.fromfile(file, dtype=np.int32,
                                         count=self.n_m_max)

                dt = np.dtype("(%i,%i)c16" % (self.n_r_max, self.n_m_max))
                self.field_m = np.fromfile(file, dtype=dt, count=1)[0]
                self.field_m = self.field_m.T

    def write(self, filename):
        """
        This routine writes a snap using stream

        :param filename: the name of the output file
        :type filename: str
        """

        with open(filename, 'wb') as out:
            x = np.array([1], dtype=np.int32)
            x.tofile(out)
            x = np.array([self.time, self.ra, self.ek, self.pr, self.radratio,
                          self.sc, self.raxi], dtype=np.float64)
            x.tofile(out)
            x = np.array([self.n_r_max, self.n_m_max, self.m_max, self.minc,
                          self.n_phi_max], dtype=np.int32)
            x.tofile(out)

            self.radius.tofile(out)
            self.tcond.tofile(out)
            self.idx2m.tofile(out)
            self.field_m.T.tofile(out)

    def frame2vtk(self, filename, name):
        """
        This routine converts the snapshot as a vtk file.

        :parameter filename: name of the vts file (without the trailing vts)
        :type filename: str
        :parameter name: name of the field (to appear in paraview scalars)
        :type filename: str
        """
        try:
            from evtk.hl import structuredToVTK
            gridToVTK = structuredToVTK
        except:
            import evtk
            gridToVTK = evtk.hl.gridToVTK

        field = spec_spat(self.field_m, self.n_phi_max)
        field = symmetrize(field, self.minc)

        phi = np.linspace(0., 2.*np.pi, field.shape[0])
        X = np.zeros((1, field.shape[0], self.n_r_max), dtype=np.float32)
        Y = np.zeros_like(X)
        Z = np.zeros_like(X)
        rr, pphi = np.meshgrid(self.radius, phi)
        X[0, :, :] = rr*np.cos(pphi)
        Y[0, :, :] = rr*np.sin(pphi)

        dat = np.zeros_like(X)
        dat[0, ...] = field

        point_data = {}
        point_data[name] = dat

        gridToVTK(filename, X, Y, Z, pointData=point_data)
        print('Store {}.vts'.format(filename))



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

        filename = self.get_filename('frame_phase', ivar, datadir, tag, verbose)
        if filename is not None:
            f = Frame(filename, endian=endian)
            self.phase_m = f.field_m
            self.phase = spec_spat(self.phase_m, self.n_phi_max)

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
                file = f'{prefix}_{ivar}.{tag}'
                filename = os.path.join(datadir, file)
            else:
                files = scanDir(os.path.join(datadir, f'{prefix}_*{tag}'))
                if len(files) > 0:
                    filename = files[-1]
                else:
                    return

            if os.path.exists(f'log.{tag}'):
                PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml=f'log.{tag}')
        else:
            if ivar is not None:
                files = scanDir(os.path.join(datadir, f'{prefix}_{ivar}.*'))
                if len(files) > 0:
                    filename = files[-1]
                else:
                    return
            else:
                files = scanDir(os.path.join(datadir, f'{prefix}_*'))
                if len(files) > 0:
                    filename = files[-1]
                else:
                    return
            # Determine the setup
            mask = re.compile(r'.*\.(.*)')
            ending = mask.search(files[-1]).groups(0)[0]
            if os.path.exists(f'log.{ending}'):
                PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml=f'log.{ending}')

        if not os.path.exists(filename):
            return

        if verbose:
            print(f'read {filename}')

        return filename

    def equat(self, field='vort', cm=None, levels=65, deminc=True,
              normed=True, vmax=None, vmin=None, normRad=False, stream=False,
              streamNorm='vel', streamDensity=1.5, cbar=True, label=None,
              streamColor='k', pcolor=False, rasterized=True):
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
        :param pcolor: when set to True, use pcolormesh instead of contourf
        :type pcolor: bool
        :param rasterized: when set to True, the rasterization for vector graphics
                           is turned on
        :type rasterized: bool
        """

        if field in ('om', 'vortz', 'vort', 'omega', 'Vorticity', 'Omega'):
            data = self.vortz
            if cm is None:
                try:
                    import cmocean.cm as cmo
                    cm = cmo.curl
                except ModuleNotFoundError:
                    cm = 'Spectral_r'
        elif field in ('temperature', 'Temperature', 'temp', 'Temp', 't', 'T'):
            data = self.temp
            if cm is None:
                try:
                    import cmocean.cm as cmo
                    cm = cmo.thermal
                except ModuleNotFoundError:
                    cm = 'magma'
            normed = False
        elif field in ('phase', 'Phase', 'phi', 'Phi'):
            if cm is None:
                try:
                    import cmocean.cm as cmo
                    cm = cmo.tempo
                except ModuleNotFoundError:
                    cm = 'binary'
            data = self.phase
            normed = False
        elif field in ('composition', 'Composition', 'xi', 'Xi', 'chem',
                       'Chem', 'comp', 'Comp'):
            if cm is None:
                try:
                    import cmocean.cm as cmo
                    cm = cmo.haline
                except ModuleNotFoundError:
                    cm = 'viridis'
            data = self.xi
            normed = False
        elif field in ('tfluct', 'tempfluct'):
            data = self.temp-self.temp_m[0, :].real
            cm = 'PuOr'
        elif field in ('xifluct', 'chemfluct', 'compfluct'):
            data = self.xi-self.xi_m[0, :].real
            cm = 'PiYG'
        elif field in ('us', 'Us', 'ur', 'Ur', 'vs', 'Vs', 'Vr', 'vr'):
            data = self.us
            if cm is None:
                cm = 'seismic'
        elif field in ('up', 'Up', 'uphi', 'Uphi', 'vp', 'Vp', 'Vphi', 'vphi'):
            data = self.uphi
            if cm is None:
                cm = 'seismic'
        elif field in ('Ek', 'ek', 'u2', 'U2', 'v2', 'V2', 'ekin', 'Ekin'):
            data = 0.5*(self.us**2+self.uphi**2)
            normed = False
            if cm is None:
                cm = 'Blues'

        if deminc:
            data = symmetrize(data, ms=self.minc)

        self.fig, xx, yy, im  = equatContour(data, self.radius, minc=self.minc,
                                             levels=levels, cm=cm, deminc=deminc,
                                             normed=normed, vmax=vmax, vmin=vmin,
                                             normRad=normRad, cbar=cbar,
                                             label=label, pcolor=pcolor,
                                             rasterized=rasterized)

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

            rad = np.linspace(self.radius[-1], self.radius[0], data.shape[1])
            rr, ttheta = np.meshgrid(rad, theta)
            if deminc:
                u = symmetrize(self.us[:, ::-1], self.minc)
                v = symmetrize(self.uphi[:, ::-1], self.minc)
            else:
                u = self.us
                v = self.uphi

            if self.radius[-1] > 0:
                v /= self.radius
            else: # Full disk
                v[:, :-1] /= self.radius[:-1]

            u = my_interp2d(u, self.radius[::-1], rad)
            v = my_interp2d(v, self.radius[::-1], rad)
            speed = np.sqrt(u**2+v**2)

            if streamNorm == 'vel':
                lw = 5.*speed.T/speed.max()
            else:
                lw = 5.*abs(data.T)/abs(self.vortz).max()
            ax1.streamplot(ttheta.T, rr.T, v.T, u.T, density=streamDensity,
                           linewidth=lw, color=streamColor)
            ax1.set_ylim(0., rad.max())
