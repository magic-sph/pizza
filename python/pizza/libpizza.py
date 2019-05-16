# -*- coding: utf-8 -*-
from .npfile import npfile
import numpy as np
from scipy.fftpack import dct
import glob
import os
import re


def cc2real(f):

    return 2.*np.sum(abs(f[1:, :])**2, axis=0) + f[0, :].real*f[0, :].real


def costf(f, fac=True):
    """
    This routine transform an input array from real to Chebyshev space

    :param f: the input array
    :type f: numpy.ndarray
    :param fac: normalisation factor is used
    :type f: bool
    :returns: a transformed array
    :rtype: numpy.ndarray
    """
    # nr = f.shape[-1]
    if fac:
        norm = np.sqrt(0.5/(f.shape[-1]-1))
    else:
        norm = 1.
    # fbig = np.hstack((f[..., :], f[..., -2:0:-1]))
    # fbig = fbig.astype('complex256')
    # fhat = norm*np.fft.fft(fbig, axis=-1)[..., :nr]

    fhat = norm*dct(f, type=1, axis=-1)

    return fhat


def chebgrid(nr, a, b):
    """
    This function defines a Gauss-Lobatto grid from a to b.

    >>> r_icb = 0.5 ; r_cmb = 1.5; n_r_max=65
    >>> rr = chebgrid(n_r_max, r_icb, r_cmb)

    :param nr: number of radial grid points
    :type nr: int
    :param a: lower limit of the Gauss-Lobatto grid
    :type a: float
    :param b: upper limit of the Gauss-Lobatto grid
    :type b: float
    :returns: the Gauss-Lobatto grid
    :rtype: numpy.ndarray
    """
    rst = (a+b)/(b-a)
    rr = 0.5*(rst+np.cos(np.pi*(1.-np.arange(nr+1.)/nr)))*(b-a)
    return rr


def avgField(time, field, tstart=None, tstop=None, std=False):
    """
    This subroutine computes the time-average (and the std) of a time series

    >>> ts = PizzaTs(field='e_kin', iplot=False, all=True)
    >>> us2vg = avgField(ts.time, ts.us2, 0.35)
    >>> print(us2avg)

    :param time: time
    :type time: numpy.ndarray
    :param field: the time series of a given field
    :type field: numpy.ndarray
    :param tstart: the starting time of the averaging
    :type tstart: float
    :param tstart: the stopping time of the averaging
    :type tstart: float
    :param std: when set to True, the standard deviation is also calculated
    :type std: bool
    :returns: the time-averaged quantity
    :rtype: float
    """
    if tstart is not None:
        mask = np.where(abs(time-tstart) == min(abs(time-tstart)), 1, 0)
        ind = np.nonzero(mask)[0][0]
    else:  # the whole input array is taken!
        ind = 0
    if tstop is not None:
        mask = np.where(abs(time-tstop) == min(abs(time-tstop)), 1, 0)
        ind1 = np.nonzero(mask)[0][0]+1
    else:  # the whole input array is taken!
        ind1 = len(time)
    fac = 1./(time[ind1-1]-time[ind])
    avgField = fac*np.trapz(field[ind:ind1], time[ind:ind1])

    if std:
        stdField = np.sqrt(fac*np.trapz((field[ind:ind1]-avgField)**2,
                           time[ind:ind1]))
        return avgField, stdField
    else:
        return avgField


def get_dr(f):
    """
    This routine calculates the first radial derivative of a input array using
    Chebyshev recurrence relation.

    :param f: the input array
    :type f: numpy.ndarray
    :returns: the radial derivative of f
    :rtype: numpy.ndarray
    """
    Nr = f.shape[-1]
    fhat = costf(f)

    # eps = np.finfo(1.0e0).eps
    # valmin = 500. * eps*abs(fhat).max()

    df = np.zeros_like(fhat)
    df[..., -1] = 0.
    df[..., -2] = (Nr-1)*fhat[..., -1]

    for i in range(Nr-3, -1, -1):
        df[..., i] = df[..., i+2]+2.*(i+1)*fhat[..., i+1]

    df[..., :] = 2.*df[..., :]

    df = costf(df)

    return df


def intcheb(f):
    """
    This routine computes an integration of a function along radius using
    Chebyshev recurrence relation.

    :param f: the input array
    :type f: numpy.ndarray
    :returns: the integral of f between z1 and z2
    :rtype: float
    """
    nr = f.shape[-1]-1
    w2 = costf(f)
    w2[..., 0] *= 0.5
    w2[..., -1] *= 0.5

    if len(f.shape) == 1:
        int = 0.
        for i in range(0, nr+1, 2):
            int = int-1./(i**2-1)*w2[i]
    elif len(f.shape) == 2:
        int = np.zeros(f.shape[0], dtype=f.dtype)
        for m in range(f.shape[0]):
            int[m] = 0.
            for i in range(0, nr+1, 2):
                int[m] = int[m]-1./(i**2-1)*w2[m, i]

    # Be careful if a mapping is used this would be modified
    int *= np.sqrt(2./(f.shape[-1]-1))

    return int


def spat_spec(arr_grid, n_m_max):
    """
    This routine computes a spectral transform from a spatial represenation
    to a spectral representation.

    :param f: the input array
    :type f: numpy.ndarray
    :param n_m_max: the number of modes
    :type n_m_max: int
    :returns: an array in the spectral space
    :rtype: numpy.ndarray
    """
    n_phi = arr_grid.shape[0]
    return np.fft.fft(arr_grid, axis=0)[:n_m_max]/n_phi


def spec_spat(arr_M, n_phi_max):
    """
    This routine computes a spectral transform from a spectral represenation
    to a spatial representation.

    :param f: the input array
    :type f: numpy.ndarray
    :param n_phi_max: the number of azimuthal grid points
    :type n_phi_max: int
    :returns: an array in the physical space
    :rtype: numpy.ndarray
    """
    n_m = arr_M.shape[0]
    tmp = np.zeros((int(n_phi_max/2)+1, arr_M.shape[-1]), np.complex128)
    tmp[:n_m, :] = arr_M
    return np.fft.irfft(tmp, n=n_phi_max, axis=0)*n_phi_max


def scanDir(pattern, tfix=None):
    """
    This function sorts the files which match a given input pattern from the
    oldest to the most recent one (in the current working directory)

    >>> dat = scanDir('log.*')
    >>> print(log)

    :param pattern: a classical regexp pattern
    :type pattern: str
    :param tfix: in case you want to add only the files that are more recent
                 than   a certain date, use tfix (computer 1970 format!!)
    :type tfix: float
    :returns: a list of files that match the input pattern
    :rtype: list
    """
    dat = [(os.stat(i).st_mtime, i) for i in glob.glob(pattern)]
    dat.sort()
    if tfix is not None:
        out = []
        for i in dat:
            if i[0] > tfix:
                out.append(i[1])
    else:
        out = [i[1] for i in dat]
    return out


def symmetrize(data, ms, reversed=False):
    """
    Symmetrise an array which is defined only with an azimuthal symmetry
    minc=ms

    :param data: the input array
    :type data: numpy.ndarray
    :param ms: the azimuthal symmetry
    :type ms: int
    :param reversed: set to True, in case the array is reversed (i.e. n_phi
                     is the last column)
    :type reversed: bool
    :returns: an output array of dimension (data.shape[0]*ms+1)
    :rtype: numpy.ndarray
    """
    if reversed:
        nphi = data.shape[-1]*ms+1
        size = [nphi]
        size.insert(0, data.shape[-2])
        if len(data.shape) == 3:
            size.insert(0, data.shape[-3])
        out = np.zeros(size, dtype=data.dtype)
        for i in range(ms):
            out[..., i*data.shape[-1]:(i+1)*data.shape[-1]] = data
        out[..., -1] = out[..., 0]
    else:
        nphi = data.shape[0]*ms + 1
        size = [nphi]
        if len(data.shape) >= 2:
            size.append(data.shape[1])
        if len(data.shape) == 3:
            size.append(data.shape[2])
        out = np.zeros(size, dtype=data.dtype)
        for i in range(ms):
            out[i*data.shape[0]:(i+1)*data.shape[0], ...] = data
        out[-1, ...] = out[0, ...]

    return out


def fast_read(file, skiplines=0, binary=False, precision=np.float64):
    """
    This function reads an input ascii table
    (can read both formatted or unformatted fortran)

    >>> # Read 'e_kin.test', skip the first 10 lines
    >>> data = fast_read('e_kin.test', skiplines=10)

    :param file: name of the input file
    :type file: str
    :param skiplines: number of header lines to be skept during reading
    :type skiplines: int
    :param binary: when set to True, try to read an unformatted binary
                   Fortran file
                   (default is False)
    :type binary: bool
    :param precision: single (np.float32) or double precision (np.float64)
    :type precision: str
    :returns: an array[nlines, ncols] that contains the data of the ascii file
    :rtype: numpy.ndarray
    """
    if not binary:
        f = open(file, 'r')
        X = []
        for k, line in enumerate(f.readlines()):
            st = line.replace('D', 'E')
            if k >= skiplines:
                X.append(st.split())
        X = np.array(X, dtype=precision)
        f.close()
    else:
        f = npfile.npfile(file, endian='B')
        X = []
        while 1:
            try:
                X.append(f.fort_read(precision))
            except TypeError:
                break
        X = np.array(X, dtype=precision)
        f.close()

    return X

def getCpuTime(file):
    """
    This function calculates the CPU time from one given log file

    :param file: the log file you want to analyze
    :type file: file
    :returns: the total CPU time
    :rtype: float
    """
    threads = re.compile(r'[\s]*\![\s]*nThreads\:[\s]*(.*)')
    ranks = re.compile(r'[\s\w]*n_procs[\s\w]*=[\s]*(.*)')
    runTime = re.compile(r'[\s\!\w]*run time:[\s]*([0-9]*)\sh[\s]*([0-9]*)\sm[\s]*([0-9]*)\ss[\s]*([0-9]*).*')
    f = open(file, 'r')
    tab = f.readlines()
    nThreads = 1  # In case a pure MPI version is used
    nRanks = 1  # In case the old OpenMP version is used
    realTime = 0.
    for line in tab:
        if threads.match(line):
            nThreads = int(threads.search(line).groups()[0])
        elif ranks.match(line):
            nRanks = int(ranks.search(line).groups()[0])
        elif runTime.match(line):
            hours = int(runTime.search(line).groups()[0])
            min = int(runTime.search(line).groups()[1])
            sec = int(runTime.search(line).groups()[2])
            ms = int(runTime.search(line).groups()[3])
            realTime = hours+1./60*min+1./3600*sec+1./3.6e6*ms
    f.close()
    cpuTime = nThreads*nRanks*realTime

    return cpuTime


def getTotalRunTime(datadir='.'):
    """
    This function calculates the total CPU time of one run directory

    :param datadir: name of the directory
    :type datadir: str
    :returns: the total RUN time
    :rtype: float
    """
    logFiles = glob.glob(os.path.join(datadir, 'log.*'))
    totCpuTime = 0
    for file in logFiles:
        totCpuTime += getCpuTime(file)

    return totCpuTime
