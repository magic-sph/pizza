# -*- coding: utf-8 -*-
import numpy as np
from .libpizza import scanDir
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

        filename = self.get_filename('frame_temp_3D', ivar, datadir, tag,
                                     verbose)
        f = Frame3D(filename, endian=endian)
        self.temp = f.field

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
