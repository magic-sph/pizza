# -*- coding: utf-8 -*-
import numpy as np
from .libpizza import scanDir, symmetrize
from .frame import my_interp2d
from .plotlib import equatContour, merContour, radialContour
from .log import PizzaSetup
import os
import re


class Frame3D:
    """
    This module is used to read 3-D snapshots
    """

    def __init__(self, filename, endian='l'):
        """
        :param filename: name of the file
        :type filename: str
        :param endian: endianness of the file ('B' or 'l')
        :type endian: str
        """

        file = open(filename, 'rb')

        self.version = np.fromfile(file, dtype='i4', count=1)[0]
        self.time, self.ra, self.ek, self.pr, self.radratio, self.sc, \
            self.raxi = np.fromfile(file, dtype='7Float64', count=1)[0]
        self.n_r_max_3D, self.l_max, self.m_max_3D, self.lm_max, \
            self.minc_3D, self.n_theta_max, self.n_phi_max_3D = \
            np.fromfile(file, dtype='7i4', count=1)[0]

        self.radius_3D = np.fromfile(file, dtype='%iFloat64' % self.n_r_max_3D,
                                     count=1)[0]
        self.tcond_3D = np.fromfile(file, dtype='%iFloat64' % self.n_r_max_3D,
                                    count=1)[0]
        self.theta = np.fromfile(file, dtype='%iFloat64' % self.n_theta_max,
                                 count=1)[0]

        dt = np.dtype("(%i,%i,%i)Float32" % (self.n_r_max_3D,
                                             self.n_theta_max,
                                             self.n_phi_max_3D))
        self.field = np.fromfile(file, dtype=dt, count=1)[0]
        self.field = self.field.T


class Pizza3DFields(PizzaSetup):
    """
    This module is used to read the 3-D frames
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

        filename = self.get_filename('frame_ur_3D', ivar, datadir, tag,
                                     verbose)
        f = Frame3D(filename, endian=endian)
        self.ra = f.ra
        self.ek = f.ek
        self.pr = f.pr
        self.radratio = f.radratio
        self.sc = f.sc
        self.raxi = f.raxi
        self.n_r_max_3D = f.n_r_max_3D
        self.l_max = f.l_max
        self.m_max_3D = f.m_max_3D
        self.lm_max = f.lm_max
        self.minc_3D = f.minc_3D
        self.n_theta_max = f.n_theta_max
        self.n_phi_max_3D = f.n_phi_max_3D
        self.radius_3D = f.radius_3D
        self.theta = f.theta
        self.tcond_3D = f.tcond_3D
        self.time = f.time

        self.ur = f.field

        filename = self.get_filename('frame_ut_3D', ivar, datadir, tag,
                                     verbose)
        f = Frame3D(filename, endian=endian)
        self.utheta = f.field

        filename = self.get_filename('frame_up_3D', ivar, datadir, tag,
                                     verbose)
        f = Frame3D(filename, endian=endian)
        self.uphi = f.field

        try:
            filename = self.get_filename('frame_temp_3D', ivar, datadir, tag,
                                         verbose)
            f = Frame3D(filename, endian=endian)
        	self.temp = f.field
        except TypeError:
            print('frame_temp_3D.%s heat fields not found' % tag)
            pass

        try:
            filename = self.get_filename('frame_br_3D', ivar, datadir, tag,
                                         verbose)
            f = Frame3D(filename, endian=endian)
            self.br = f.field

            filename = self.get_filename('frame_bt_3D', ivar, datadir, tag,
                                         verbose)
            f = Frame3D(filename, endian=endian)
            self.bt = f.field

            filename = self.get_filename('frame_bp_3D', ivar, datadir, tag,
                                         verbose)
            f = Frame3D(filename, endian=endian)
            self.bp = f.field
        except TypeError:
            print('frame_b*_3D.%s mag fields not found' % tag)
            pass

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

    def avg(self, field='up', cm='seismic', levels=65, normed=True, vmax=None,
            vmin=None, cbar=True, label=None):
        """
        Display an azimuthal average of a scalar quantity

        :param field: the name of the input physical quantity you want to
                      display
        :type field: str
        :param levels: the number of levels in the contour
        :type levels: int
        :param cm: name of the colormap ('jet', 'seismic', 'RdYlBu_r', etc.)
        :type cm: str
        :param cbar: display the colorbar when set to True
        :type cbar: bool
        :param vmax: maximum value of the contour levels
        :type vmax: float
        :param vmin: minimum value of the contour levels
        :type vmin: float
        :param normed: when set to True, the colormap is centered around zero.
                       Default is True, except for entropy/temperature plots.
        :type normed: bool
        """

        if field in ('temperature', 'Temperature', 'temp', 'Temp', 't', 'T'):
            data = self.temp.mean(axis=0)
        elif field in ('ur', 'Ur', 'Vr', 'vr'):
            data = self.ur.mean(axis=0)
        elif field in ('ut', 'Ut', 'utheta', 'Utheta', 'vt', 'Vt', 'Vtheta',
                       'vtheta'):
            data = self.utheta.mean(axis=0)
        elif field in ('up', 'Up', 'uphi', 'Uphi', 'vp', 'Vp', 'Vphi', 'vphi'):
            data = self.uphi.mean(axis=0)
        elif field in ('br', 'Br'):
            data = self.br.mean(axis=0)
        elif field in ('bt', 'Bt', 'bth', 'Bth', 'btheta', 'Btheta'):
            data = self.bt.mean(axis=0)
        elif field in ('bp', 'Bp', 'bph', 'Bph', 'bphi', 'Bphi'):
            data = self.bp.mean(axis=0)

        self.fig, xx, yy = merContour(data, self.radius_3D, colat=self.theta,
                                      label=label,
                                      levels=levels, cm=cm, normed=normed,
                                      vmax=vmax, vmin=vmin, cbar=cbar,
                                      tit=False)

    def slice(self, p=0.0, field='up', cm='seismic', levels=65, normed=True, vmax=None,
            vmin=None, cbar=True, label=None):
        """
        Display an meridional slice of a scalar quantity

        :param field: the name of the input physical quantity you want to
                      display
        :type field: str
        :param p: azimuth at which the cut is done
        :type p: float
        :param levels: the number of levels in the contour
        :type levels: int
        :param cm: name of the colormap ('jet', 'seismic', 'RdYlBu_r', etc.)
        :type cm: str
        :param cbar: display the colorbar when set to True
        :type cbar: bool
        :param vmax: maximum value of the contour levels
        :type vmax: float
        :param vmin: minimum value of the contour levels
        :type vmin: float
        :param normed: when set to True, the colormap is centered around zero.
                       Default is True, except for entropy/temperature plots.
        :type normed: bool
        """

        self.phi = np.linspace(self.theta[0]/2.,np.pi,self.n_phi_max_3D)
        phmax = self.phi[-1]
        pnorm = self.phi/phmax
        pmask = np.where(abs(pnorm-p) == abs(pnorm-p).min(), 1, 0)
        iphi = np.nonzero(pmask)[0][0]
        if field in ('temperature', 'Temperature', 'temp', 'Temp', 't', 'T'):
            data = self.temp[iphi, :, :]
        elif field in ('tfluct', 'tempfluct'):
            data = self.temp-self.temp.mean(axis=0)
            data = data[iphi, :, :]
        elif field in ('ur', 'Ur', 'Vr', 'vr'):
            data = self.ur[iphi, :, :]
        elif field in ('ut', 'Ut', 'utheta', 'Utheta', 'vt', 'Vt', 'Vtheta',
                       'vtheta'):
            data = self.utheta[iphi, :, :]
        elif field in ('up', 'Up', 'uphi', 'Uphi', 'vp', 'Vp', 'Vphi', 'vphi'):
            data = self.uphi[iphi, :, :]
        elif field in ('br', 'Br'):
            data = self.br[iphi, :, :]
        elif field in ('bt', 'Bt','bth', 'Bth', 'btheta', 'Btheta'):
            data = self.bt[iphi, :, :]
        elif field in ('bp', 'Bp', 'bph', 'Bph', 'bphi', 'Bphi'):
            data = self.bp[iphi, :, :]

        self.fig, xx, yy = merContour(data, self.radius_3D, colat=self.theta,
                                      label=label,
                                      levels=levels, cm=cm, normed=normed,
                                      vmax=vmax, vmin=vmin, cbar=cbar,
                                      tit=False)

    def equat(self, field='ur', cm='seismic', levels=65, deminc=True,
              normed=True, vmax=None, vmin=None, normRad=False, stream=False,
              streamNorm='vel', streamDensity=1.5, cbar=True, label=None,
              streamColor='k'):
        """
        Display an equatorial cut of a scalar quantity

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

        if field in ('temperature', 'Temperature', 'temp', 'Temp', 't', 'T'):
            data = self.temp[:, self.n_theta_max/2, :]
        elif field in ('tfluct', 'tempfluct'):
            data = self.temp-self.temp.mean(axis=0)
            data = data[:, self.n_theta_max/2, :]
        elif field in ('ur', 'Ur', 'Vr', 'vr'):
            data = self.ur[:, self.n_theta_max/2, :]
        elif field in ('ut', 'Ut', 'utheta', 'Utheta', 'vt', 'Vt', 'Vtheta',
                       'vtheta'):
            data = self.utheta[:, self.n_theta_max/2, :]
        elif field in ('up', 'Up', 'uphi', 'Uphi', 'vp', 'Vp', 'Vphi', 'vphi'):
            data = self.uphi[:, self.n_theta_max/2, :]
        elif field in ('br', 'Br'):
            data = self.br[:, self.n_theta_max/2, :]
        elif field in ('bt', 'Bt', 'bth', 'Bth', 'btheta', 'Btheta'):
            data = self.bt[:, self.n_theta_max/2, :]
        elif field in ('bp', 'Bp', 'bph', 'Bph', 'bphi', 'Bphi'):
            data = self.bp[:, self.n_theta_max/2, :]

        if deminc:
            data = symmetrize(data, ms=self.minc_3D)

        self.fig, xx, yy = equatContour(data, self.radius_3D,
                                        minc=self.minc_3D, levels=levels,
                                        cm=cm, deminc=deminc, normed=normed,
                                        vmax=vmax, vmin=vmin, normRad=normRad,
                                        cbar=cbar, label=label)

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
                phi = np.linspace(0., 2.*np.pi, data.shape[0])
            else:
                phi = np.linspace(0., 2*np.pi/self.minc_3D, data.shape[0])

            rad = np.linspace(self.radius_3D[0], self.radius_3D[-1],
                              data.shape[-1])
            rr, pphi = np.meshgrid(rad, phi)
            if deminc:
                u = symmetrize(self.ur[:, self.n_theta_max/2, :], self.minc_3D)
                v = symmetrize(self.uphi[:, self.n_theta_max/2, :],
                               self.minc_3D)
            else:
                u = self.ur[:, self.n_theta_max/2, :]
                v = self.uphi[:, self.n_theta_max/2, :]
            v /= self.radius_3D

            u = my_interp2d(u, self.radius_3D, rad)
            v = my_interp2d(v, self.radius_3D, rad)
            speed = np.sqrt(u**2+v**2)

            if streamNorm == 'vel':
                lw = 3.*speed.T/speed.max()
            ax1.streamplot(pphi.T, rr.T, v.T, u.T, density=streamDensity,
                           linewidth=lw, color=streamColor)
            ax1.set_ylim(0., rad.max())

    def surf(self, r=0.5, field='ur', cm='seismic', levels=65, deminc=True,
             normed=True, vmax=None, vmin=None, cbar=True, label=None,
             streamColor='k', lon_0=0, lines=False):
        """
        Display an radial cut of a scalar quantity

        :param field: the name of the input physical quantity you want to
                      display
        :type field: str
        :param r: radius at which the cut is done
        :type r: float
        :param lon_0: the central longitude
        :type lon_0: float
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
        :param lines: when set to True, overplot the solid contour levels
        :type lines: bool
        :param normed: when set to True, the colormap is centered around zero.
                       Default is True, except for entropy/temperature plots.
        :type normed: bool
        """

        rcmb = self.radius_3D[0]
        rnorm = self.radius_3D/rcmb
        mask = np.where(abs(rnorm-r) == abs(rnorm-r).min(), 1, 0)
        irad = np.nonzero(mask)[0][0]
        if field in ('temperature', 'Temperature', 'temp', 'Temp', 't', 'T'):
            data = self.temp[:, :, irad]
        elif field in ('tfluct', 'tempfluct'):
            data = self.temp-self.temp.mean(axis=0)
            data = data[:, :, irad]
        elif field in ('ur', 'Ur', 'Vr', 'vr'):
            data = self.ur[:, :, irad]
        elif field in ('ut', 'Ut', 'utheta', 'Utheta', 'vt', 'Vt', 'Vtheta',
                       'vtheta'):
            data = self.utheta[:, :, irad]
        elif field in ('up', 'Up', 'uphi', 'Uphi', 'vp', 'Vp', 'Vphi', 'vphi'):
            data = self.uphi[:, :, irad]
        elif field in ('br', 'Br'):
            data = self.br[:, :, irad]
        elif field in ('bt', 'Bt','bth', 'Bth', 'btheta', 'Btheta'):
            data = self.bt[:, :, irad]
        elif field in ('bp', 'Bp', 'bph', 'Bph', 'bphi', 'Bphi'):
            data = self.bp[:, :, irad]

        if deminc:
            data = symmetrize(data, ms=self.minc_3D)

        fig = radialContour(data, rad=r, label=label, lon_0=lon_0, vmax=vmax,
                            vmin=vmin, levels=levels, cm=cm, normed=normed,
                            cbar=cbar, tit=False, lines=lines)
