from .npfile import *
import numpy as np
import glob
import os

def spat_spec(arr_grid, n_m_max):
    n_phi = arr_grid.shape[0]
    return np.fft.fft(arr_grid, axis=0)[:n_m_max]/n_phi

def spec_spat(arr_M, n_phi_max):
    n_m = arr_M.shape[0]
    tmp = np.zeros((int(n_phi_max/2)+1, arr_M.shape[-1]), 'Complex64')
    tmp[:n_m, :] = arr_M
    return np.fft.irfft(tmp, n=n_phi_max, axis=0)*n_phi_max

def scanDir(pattern, tfix=None):
    """
    This function sorts the files which match a given input pattern from the oldest
    to the most recent one (in the current working directory)
    
    >>> dat = scanDir('log.*')
    >>> print(log)

    :param pattern: a classical regexp pattern
    :type pattern: str
    :param tfix: in case you want to add only the files that are more recent than   
                 a certain date, use tfix (computer 1970 format!!)
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
    Symmetrise an array which is defined only with an azimuthal symmetry minc=ms

    :param data: the input array
    :type data: numpy.ndarray
    :param ms: the azimuthal symmetry
    :type ms: int
    :param reversed: set to True, in case the array is reversed (i.e. n_phi is the last column)
    :type reversed: bool
    :returns: an output array of dimension (data.shape[0]*ms+1)
    :rtype: numpy.ndarray
    """
    if reversed:
        nphi = data.shape[-1]*ms+1
        size = [nphi]
        size.insert(0,data.shape[-2])
        if len(data.shape) == 3:
            size.insert(0,data.shape[-3])
        out = np.zeros(size, dtype=data.dtype)
        for i in range(ms):
            out[..., i*data.shape[-1]:(i+1)*data.shape[-1]] = data
        out[..., -1] = out[..., 0]
    else:
        nphi = data.shape[0]*ms +1
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


def fast_read(file, skiplines=0, binary=False, precision='Float64'):
    """
    This function reads an input ascii table
    (can read both formatted or unformatted fortran)

    >>> # Read 'e_kin.test', skip the first 10 lines
    >>> data = fast_read('e_kin.test', skiplines=10)

    :param file: name of the input file
    :type file: str
    :param skiplines: number of header lines to be skept during reading
    :type skiplines: int
    :param binary: when set to True, try to read an unformatted binray Fortran file
                   (default is False)
    :type binary: bool
    :param precision: single ('Float32') or double precision ('Float64')
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
        f = npfile(file, endian='B')
        X = []
        while 1:
            try:
                X.append(f.fort_read(precision))
            except TypeError:
                break
        X = np.array(X, dtype=precision)
        f.close()
    return X
