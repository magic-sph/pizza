from scipy import optimize
from scipy.special import jv, yn, hankel1
from scipy.interpolate import interp1d
from .libpizza import chebgrid, intcheb, progressbar
from multiprocessing import Pool
import numpy as np
import time
import os, sys, uuid
import pickle
try:
    from scipy.integrate import simps
except:
    from scipy.integrate import simpson as simps
import matplotlib.pyplot as plt

# For m > msafety, the roots of the annulus can be approximated
# by only retaining the type I Bessel terms
msafety = 200
# This is the cut-off value (x*rin < safety*m) which allows to
# truncate the expansion
safety = 0.8

def remap(fold, rold, rnew):
    """
    This function remaps onto a new grid using splines
    """
    ro = rold[0]
    ri = rold[-1]

    if ro > ri:
        ir = interp1d(rold[::-1], fold[..., ::-1], axis=-1, kind='cubic')
    else:
        ir = interp1d(rold, fold, axis=-1, kind='cubic')
    fnew = ir(rnew)

    return fnew


def jnyn_zeros_guess(rin, rout, m, zeros_m_minus_one, nroots):
    """
    This subroutine determines the roots of the transcendental equation
    when Dirichlet boundary conditions are considered. Number of roots
    is specified by nroots. This makes use of the roots found for m-1
    to determine the initial brackets.

    :param rin: inner radius
    :type rin: float
    :param rout: outer radius
    :type rout: float
    :param zeros_m_minus_one: the zeros obtained at m-1
    :type zeros_m_minus_one: numpy.ndarray
    :param m: azimuthal wavenumber
    :type m: int
    :param nroots: nroots
    :type nroots: number of roots to be computed
    :returns: a numpy array which contains the roots
    :rtype: numpy.ndarray
    """

    if rin == 0.: # Full disk configuration
        def f(x):
            return jv(m, x)
    else:
        def f(x):
            if m > msafety and x*rin < safety*m:
                return -jv(m, x*rout)
            else:
                return (jv(m, x*rout)*yn(m, x*rin) - jv(m, x*rin)*yn(m, x*rout)) / \
                        abs(hankel1(m, x*rin))

    len0 = len(zeros_m_minus_one)

    import mpmath as mp

    if rin == 0.:
        def fmpmath(x):
            return mp.besselj(m, x)
    else:
        def fmpmath(x):
            if m > msafety and x*rin < safety*m:
                return -mp.besselj(m, x*rout)
            else:
                return (mp.besselj(m, x*rout)*mp.bessely(m, x*rin) -  \
                        mp.besselj(m, x*rin)*mp.bessely(m, x*rout)) / \
                        abs(mp.hankel1(m, x*rin))

    roots = np.zeros(nroots, np.float64)
    nroot = 0
    while nroot < nroots:
        if nroot < len0:
            start = zeros_m_minus_one[nroot]
        else:
            spacing = roots[nroot-1]-roots[nroot-2]
            start = roots[nroot-1]+0.5*spacing
        if nroot < len0-1:
            spacing = zeros_m_minus_one[nroot+1] - zeros_m_minus_one[nroot]
            stop = zeros_m_minus_one[nroot]+0.5*spacing
        else:
            spacing = roots[nroot-1]-roots[nroot-2]
            stop = roots[nroot-1]+1.5*spacing

        compute_mpmath = False
        try:
            sol = optimize.root_scalar(f, bracket=[start, stop],
                                       method='brentq')
            if sol.converged and abs(sol.root) > 1e-1:
                root = sol.root
            else:
                compute_mpmath = True
        except ValueError:
            compute_mpmath = True

        if compute_mpmath:
            root = mp.findroot(fmpmath, x0=[start, stop])
        roots[nroot] = root
        nroot += 1

    return roots


def jnyn_zeros(rin, rout, m, nroots):
    """
    This subroutine determines the roots of the transcendental equation
    when Dirichlet boundary conditions are considered. Number of roots
    is specified by nroots.

    :param rin: inner radius
    :type rin: float
    :param rout: outer radius
    :type rout: float
    :param m: azimuthal wavenumber
    :type m: int
    :param nroots: nroots
    :type nroots: number of roots to be computed
    :returns: a numpy array which contains the roots
    :rtype: numpy.ndarray
    """

    def f(x):
        if m > msafety and x*rin < safety*m:
            return -jv(m, x*rout)
        else:
            return (jv(m, x*rout)*yn(m, x*rin) - jv(m, x*rin)*yn(m, x*rout)) / \
                    abs(hankel1(m, x*rin))

    def full_annulus(x):
        return (jv(m, x*rout)*yn(m, x*rin) - jv(m, x*rin)*yn(m, x*rout)) / \
                abs(hankel1(m, x*rin))
    def simplified_annulus(x):
        return -jv(m, x*rout)

    if m <= 15:
        npoints = 1024
        x = np.linspace(0.1, 100, npoints)
    else:
        # Fit for the first root
        xmin = int(1./rout*m)
        # Fit for the fifth root
        xmax = int(1.06/rout*m+30)
        dx = 0.2
        npoints = int((xmax-xmin)/dx)
        # Spacing between two m for first roots is like
        # lambda(m+1)-lambda(m)\sim0.65
        x = np.linspace(xmin, xmax, npoints)
    if m > msafety:
        mask = (x*rin < safety*m)
        y = np.zeros_like(x)
        y[mask] = simplified_annulus(x[mask])
        y[~mask] = full_annulus(x[~mask])
    else:
        y = full_annulus(x)
    mask = np.isnan(y)
    improve_prec = False
    if len(x[mask]) > 0:
        improve_prec = True

    if improve_prec:
        print('Precision of scipy is not enough, using mpmath')

        import mpmath as mp

        def fmpmath(x):
            if m > msafety and x*rin < safety*m:
                return -mp.besselj(m, x*rout)
            else:
                return (mp.besselj(m, x*rout)*mp.bessely(m, x*rin) -  \
                        mp.besselj(m, x*rin)*mp.bessely(m, x*rout)) / \
                        abs(mp.hankel1(m, x*rin))

        xmax = x[mask][-1]
        xmin = x[mask][0]
        dx = 1.
        lenx = int((xmax-xmin)/dx)
        xprec = np.linspace(xmin, xmax, lenx)
        yprec = np.zeros_like(xprec)
        for i in range(len(xprec)):
            yprec[i] = fmpmath(xprec[i])
        y = np.hstack((yprec, y[~mask]))
        x = np.hstack((xprec, x[~mask]))
        npoints = len(x)
    else:
        xmax = x[-1]

    signs = np.sign(y)
    nroot = 0
    roots = np.zeros((nroots), np.float64)
    for i in range(npoints-1):
        if signs[i] + signs[i+1] == 0:
            if improve_prec and x[i] < 1.1*xmax:
                root = mp.findroot(fmpmath, x0=[x[i], x[i+1]])
            else:
                sol = optimize.root_scalar(f, bracket=[x[i], x[i+1]],
                                           method='brentq')
                if sol.converged and abs(sol.root) > 1e-1:
                    root = sol.root
                else:
                    root = mp.findroot(fmpmath, x0=[x[i], x[i+1]])
                # print(root)
            if nroot < nroots:
                roots[nroot] = root

            nroot += 1

    # Only works if previous spacing between 2 roots is not far from next
    while nroot < nroots:
        spacing = roots[nroot-1] - roots[nroot-2]
        start = roots[nroot-1] + 0.5*spacing
        stop = roots[nroot-1] + 1.5*spacing
        if improve_prec and start < 1.1*xmax:
            root = mp.findroot(fmpmath, x0=[start, stop])
        else:
            sol = optimize.root_scalar(f, bracket=[start, stop],
                                       method='brentq')
            root = sol.root
        roots[nroot] = root

        nroot += 1

    """
    print(roots)
    plt.plot(x, f(x))
    #plt.ylim(-10*f(x[-1]), 10*f(x[-1]))
    plt.plot(roots, f(roots), ls='None', marker='o', mfc='None')
    if improve_prec:
        print(xmin, xmax)
        xmp = mp.linspace(xmin, xmax, 128)
        yy = mp.linspace(0., 0., len(xmp))
        for i in range(len(xmp)):
            yy[i] = fmpmath(xmp[i])
        i = 0
        while roots[i] < xmax:
            plt.plot(roots[i], fmpmath(roots[i]), ls='None', marker='o',
                     mfc='None', mec='C2')
            i += 1
        plt.plot(xmp, yy)

    plt.show()

    #print(f(roots))
    """

    return roots


def jn_zeros(m, nroots):
    """
    This subroutine determines the roots of the transcendental equation
    when Dirichlet boundary conditions are considered in the full disk.

    :param m: azimuthal wavenumber
    :type m: int
    :param nroots: nroots
    :type nroots: number of roots to be computed
    :returns: a numpy array which contains the roots
    :rtype: numpy.ndarray
    """
    def f(x):
        return jv(m, x)

    if m < 10:
        npoints = 1024
        x = np.linspace(m, 100, npoints)
    else:
        npoints = 10*m
        x = np.linspace(m, 2*m, npoints)
    signs = np.sign(f(x))
    nroot = 0
    roots = np.zeros((nroots), np.float64)
    for i in range(npoints-1):
        if signs[i] + signs[i+1] == 0:
            sol = optimize.root_scalar(f, bracket=[x[i], x[i+1]], method='brentq')
            if nroot < nroots:
                roots[nroot] = sol.root
            nroot += 1

    # Only works if previous spacing between 2 roots is not far from next
    while nroot < nroots:
        spacing = roots[nroot-1] - roots[nroot-2]
        start = roots[nroot-1] + 0.5*spacing
        stop = roots[nroot-1] + 1.5*spacing
        sol = optimize.root_scalar(f, bracket=[start, stop], method='brentq')
        roots[nroot] = sol.root

        nroot += 1

    return roots



class HankelDisk:
    """
    This is an implementation of the discrete Hankel transform for a disk.
    This is based on ``The Discrete Hankel Transform'' by Natalie Baddour.
    """

    def __init__(self, rout, m, nroots, storage_dir=None, verbose=False):
        """
        :param rout: outer radius
        :type rout: float
        :param m: azimuthal wavenumber
        :type m: int
        :param nroots: nroots
        :type nroots: number of roots to be computed
        :param storage_dir: a directory when one case store the Hankel
                            transforms
        :type storage_dir: char
        :param verbose: a boolean to trigger some printouts of infos
        :type verbose: bool
        """
        compute = True
        store = False
        if storage_dir is not None:
            filename = f'Hankel_{rout:.2f}_{nroots:04d}_{m:04d}.pickle'
            file = os.path.join(storage_dir, filename)
            if os.path.exists(file):
                self._read(file, verbose)
                compute = False
            else:
                store = True

        if compute:
            self.rout = rout
            self.m = m
            self.nroots = nroots

            # Get the roots
            if storage_dir is not None:
                filename = \
                     f'Hankel_{rout:.2f}_{nroots:04d}_{m-1:04d}.pickle'
                file = os.path.join(storage_dir, filename)
                if os.path.exists(file):
                    with open(file, 'rb') as fi:
                        tmp = pickle.load(fi)
                        roots_m_minus_1, grid_m_minus_1, k_m_minus_1 = \
                            pickle.load(fi)
                    self.roots = jnyn_zeros_guess(0., self.rout, self.m,
                                                  roots_m_minus_1,
                                                  nroots+1)
                else:
                    self.roots = jn_zeros(self.m, nroots+1)
            else:
                self.roots = jn_zeros(self.m, nroots+1)

            self.roots_last = self.roots[-1]
            self.roots = self.roots[:-1]

            self.grid = self.roots * self.rout / self.roots_last

            self.kr = self.roots / self.rout

            self.compute_kernels()

            if store:
                self._store(storage_dir)

    def compute_kernels(self):
        """
        This routine is used to store the kernels employed in the Hankel
        transforms on the disk
        """

        self.jp1 = abs(jv(self.m+1, self.roots))
        jp = jv(self.m, self.roots[:, np.newaxis]@self.roots[np.newaxis, :] /
                self.roots_last)
        self.kernels = 2.*jp/(self.jp1[:, np.newaxis]@self.jp1[np.newaxis, :] *
                              self.roots_last)

    def _read(self, file, verbose):
        """
        This routine is used to read the inverse matrix and the roots

        :param file: the pickle file which contains the structure
        :type file: char
        :param verbose: a boolean to trigger some printouts of infos
        :type verbose: bool
        """
        with open(file, 'rb') as f:
            if verbose:
                print(f'Reading {file}')
            self.rout, self.m, self.nroots, self.roots_last = pickle.load(f)
            self.roots, self.grid, self.kr = pickle.load(f)
            self.jp1, self.kernels = pickle.load(f)

    def _store(self, dir):
        """
        This routine is used to store the inverse matrix and the roots

        :param dir: a directory when one case store the Hankel transforms
        :type dir: char
        """

        filename = f'Hankel_{self.rout:.2f}_{self.nroots:04d}_{self.m:04d}.pickle'
        file = os.path.join(dir, filename)

        with open(file, 'wb') as f:
            pickle.dump([self.rout, self.m, self.nroots, self.roots_last], f)
            pickle.dump([self.roots, self.grid, self.kr], f)
            pickle.dump([self.jp1, self.kernels], f)

    def HT(self, f):
        """
        This is the actual discrete Hankel transform from physical
        to spectral space.

        :param f: the function defined in physical space
        :type f: numpy.ndarray
        :returns: the function in spectral space
        :rtype: numpy.ndarray
        """
        fac = self.rout**2/self.roots_last
        fhat = self.jp1 * fac * np.matmul(self.kernels, f/self.jp1)

        return fhat

    def iHT(self, fhat):
        """
        This is the inverse discrete Hankel transform on the disk.

        :param f: the function defined in spectral space
        :type f: numpy.ndarray
        :returns: the function in physical space
        :rtype: numpy.ndarray
        """

        fac = self.roots_last/self.rout**2
        f = self.jp1 * fac * np.matmul(self.kernels, fhat/self.jp1)

        return f

    def spectra(self, fhat):
        """
        This function handles the computation of the spectra
        :param fhat: the function defined in spectral space
        :type fhat: numpy.ndarray
        :returns: the energy content as a function of the order
        :rtype: numpy.ndarray
        """
        jp1 = jv(self.m+1, self.roots)
        En = 2. / self.rout**2 * abs(fhat)**2 / jp1**2

        return En


class HankelAnnulus:

    def __init__(self, rin, rout, m, ngrid_points, storage_dir=None,
                 verbose=False):
        """
        :param rin: inner radius
        :type rin: float
        :param rout: outer radius
        :type rout: float
        :param m: azimuthal wavenumber
        :type m: int
        :param ngrid_points: number of grid_points
        :type ngrid_points: int
        :param storage_dir: a directory when one case store the Hankel
                            transforms
        :type storage_dir: char
        :param verbose: a boolean to trigger some printouts of infos
        :type verbose: bool
        """
        compute = True
        store = False
        if storage_dir is not None:
            nroots = ngrid_points // 2
            filename = f'WeberOrr_{rin/rout:.2f}_{nroots:05d}_{m:05d}.npz'
            file = os.path.join(storage_dir, filename)
            if os.path.exists(file):
                self._read(file, verbose)
                compute = False
            else:
                store = True

        if compute:
            self.rin = rin
            self.rout = rout
            self.m = m
            # We use a Gauss-Lobatto grid to ensure a proper estimate
            # of the integrals once the kernels have been evaluated
            self.grid = chebgrid(ngrid_points-1, rout, rin)
            # Number of wave numbers = half the number of grid points
            self.nroots = ngrid_points // 2

            # Get the roots
            if storage_dir is not None:
                filename = \
                     f'WeberOrr_{rin/rout:.2f}_{self.nroots:05d}_{m-1:05d}.npz'
                file = os.path.join(storage_dir, filename)
                if os.path.exists(file):
                    tmp = np.load(file)
                    roots_m_minus_1 = tmp['roots']
                    self.roots = jnyn_zeros_guess(rin, rout, m,
                                                  roots_m_minus_1,
                                                  self.nroots)
                else:
                    self.roots = jnyn_zeros(rin, rout, m, self.nroots)
            else:
                self.roots = jnyn_zeros(rin, rout, m, self.nroots)

            # No need to compute wavenumber beyond the max of kzero
            # Take one extra beyond nroots to make the last bin works fine
            self.kzero = jnyn_zeros(rin, rout, 0, self.nroots+1)

            self.kmax = self.kzero[-1]
            self.kr = self.roots[self.roots <= self.kmax]

            self.compute_kernels(verbose=verbose)

            if store:
                self._store(storage_dir)

    def _read(self, file, verbose):
        """
        This routine is used to read the kernels employed in the Weber-Orr
        transform

        :param file: the npz file which contains the structure
        :type file: char
        :param verbose: a boolean to trigger some printouts of infos
        :type verbose: bool
        """
        if verbose:
            print(f'Reading {file}')

        data = np.load(file)
        self.rin = data['rin']
        self.rout = data['rout']
        self.m = data['m']
        self.nroots = data['nroots']
        self.roots = data['roots']
        self.grid = data['grid']
        self.kr = data['kr']
        self.kernels = data['kernels']

    def _store(self, dir):
        """
        This routine is used to store the inverse matrix and the roots

        :param dir: a directory when one case store the Hankel transforms
        :type dir: char
        """

        filename = f'WeberOrr_{self.rin/self.rout:.2f}_{self.nroots:05d}_{self.m:05d}.npz'
        file = os.path.join(dir, filename)

        np.savez(file, rin=self.rin, rout=self.rout, m=self.m, nroots=self.nroots,
                 grid=self.grid, roots=self.roots, kr=self.kr, kernels=self.kernels)

    def compute_kernels(self, verbose=False):
        """
        This routine is used to store the kernels employed in the inverse
        Weber-Orr transform.

        :param verbose: display some messages if requested
        :type verbose: bool
        """

        out = np.outer(self.grid, self.kr)
        if self.m > msafety:
            mask_k = (self.kr*self.rin < safety*self.m)
            k_yn = self.kr[~mask_k]
            if verbose and len(k_yn) < len(self.kr):
                ndiscards = len(self.kr)-len(k_yn)
                print(f'For wavenumber m={self.m}, I approximate the first {ndiscards} modes')
            self.kernels = np.zeros_like(out)

            self.kernels[:, mask_k] = -jv(self.m, out[:, mask_k])
            self.kernels[:, ~mask_k] = (jv(self.m, out[:, ~mask_k]) * \
                                        yn(self.m, self.rin*k_yn)  -  \
                                        yn(self.m, out[:, ~mask_k]) * \
                                        jv(self.m, self.rin*k_yn)) /  \
                                        abs(hankel1(self.m, self.rin*k_yn))
        else:
            self.kernels = (jv(self.m, out)*yn(self.m, self.rin*self.kr)  -  \
                            yn(self.m, out)*jv(self.m, self.rin*self.kr)) /  \
                            abs(hankel1(self.m, self.rin*self.kr))

        self.kernels = self.kernels.T # Wavenumber first, grid second

    def hankel_mat(self):
        """
        This subroutine constructs the Hankel matrix from the zeroes and
        the collocation grid
        """
        self.A = np.zeros((len(self.grid), len(self.grid)), np.float64)
        fac = self.kr**2 * jv(self.m, self.kr*self.rout)**2 / \
              (jv(self.m, self.kr*self.rin)**2 -
               jv(self.m, self.kr*self.rout)**2)
        norm = abs(hankel1(self.m, self.kr*self.rin))
        for k, rk in enumerate(self.grid):
            kernel = jv(self.m, self.kr*rk)*yn(self.m, self.kr*self.rin) - \
                     jv(self.m, self.kr*self.rin)*yn(self.m, self.kr*rk)
            kernel *= norm
            self.A[k, :] = fac * kernel * np.pi**2 / 2.

    def spectra(self, fhat):
        """
        :param fhat: the function defined in spectral space
        :type fhat: numpy.ndarray
        :returns: the energy content as a function of the order
        :rtype: numpy.ndarray
        """
        vmin = jv(self.m, self.kr*self.rout)**2
        vmax = jv(self.m, self.kr*self.rin)**2

        # To prevent round-off errors
        fac = np.ones_like(self.kr)
        mask = ( vmin > 1e-10 )
        fac[mask] = vmax[mask] / (vmax[mask] - vmin[mask])
        #fac[~mask] = np.ones_like(self.kr[~mask])
        En = 0.5 * np.pi**2 * fac * self.kr**2 * abs(fhat)**2 * \
             abs(hankel1(self.m, self.kr*self.rout))**2

        return En

    def HT(self, f):
        """
        This is a Hankel transform based on integration

        :param f: the function defined in physical space
        :type f: numpy.ndarray
        :returns: the function in spectral space
        :rtype: numpy.ndarray
        """
        assert f.shape[-1] == self.kernels.shape[-1]
        fhat = intcheb(self.kernels * f * self.grid)

        return fhat

    def iHT(self, fhat):
        """
        :param fhat: the function defined in spectral space
        :type fhat: numpy.ndarray
        :returns: the function in physical space
        :rtype: numpy.ndarray
        """
        vmin = jv(self.m, self.kr*self.rout)**2
        vmax = jv(self.m, self.kr*self.rin)**2
        # To prevent round-off errors
        fac = np.ones_like(self.kr)
        mask = ( vmin > 1e-10 )
        fac[mask] = vmax[mask] / (vmax[mask] - vmin[mask])
        fac *= self.kr**2*abs(hankel1(self.m, self.kr*self.rout))**2
        f = np.sum(fac*fhat*self.kernels.T, axis=1) * np.pi**2/2

        return f


def hankel_spectra(f, r, weight, idx2m, nroots, storage_dir=None, file_save=None):
    """
    This routine computes the Hankel spectra on the full disk.
    It takes a complex array of dimension (n_m_max, n_r) as an input.

    :param f: a a complex input array (first dimension is m
    :type f: numpy.ndarray
    :param r: the radius
    :type f: numpy.ndarray
    :param weight: a weighting function (typically the height of the container)
    :type weight: numpy.ndarray
    :param idx2m: an array which convert the index to ms
    :type idx2m: numpy.ndarray
    :param nroots: the number of roots
    :type nroots: int
    :param storage_dir: the directory where the Hankel transform are stored
    :type storage_dir: char
    :param file_save: filename to save the computation once done
    :type file_save: char
    """
    rout = r.max()
    nms = f.shape[0]
    Ek = np.zeros((nms, nroots), np.float64)
    K = np.zeros_like(Ek)

    if file_save is not None and os.path.exists(file_save):
        with open(file_save, 'rb') as fi:
            K, Ek = pickle.load(fi)
    else:
        for idx in progressbar(range(nms)):
            m = idx2m[idx]
            if storage_dir is not None:
                h = HankelDisk(rout, m, nroots, storage_dir=storage_dir)
            else:
                h = HankelDisk(rout, m, nroots)

            fnew = remap(f[idx, :], r, h.grid)
            weight_new = remap(weight, r, h.grid)
            fhat = h.HT(fnew*np.sqrt(weight_new))
            if m == 0:
                Ek[idx, :] = np.pi * h.spectra(fhat)
            else:
                Ek[idx, :] = 2.*np.pi * h.spectra(fhat)
            K[idx, :] = h.kr

    if file_save is not None:
        with open(file_save, 'wb') as fi:
            pickle.dump([K, Ek], fi)

    return K, Ek

def globalize(func):
    def result(*args, **kwargs):
        return func(*args, **kwargs)
    result.__name__ = result.__qualname__ = uuid.uuid4().hex
    setattr(sys.modules[result.__module__], result.__name__, result)
    return result

def weber_orr_dot_product(f, g, r, weight, idx2m, nranks=1,
                          storage_dir=None, file_save=None):
    """
    This routine computes a spectra based on the dot product of two fields.

    :param f: a complex input array (first dimension is m, second dimension
              is r) or a list of complex-type arrays. In the latter case
              the dot-product is computed over both directions.
    :type f: numpy.ndarray
    :param g: a complex input array (first dimension is m, second dimension
              is r) or a list of complex-type arrays.
    :type g: numpy.ndarray
    :param r: the radius
    :type f: numpy.ndarray
    :param weight: a weighting function (typically the height of the container)
    :type weight: numpy.ndarray
    :param idx2m: an array which convert the index to ms
    :type idx2m: numpy.ndarray
    :param nranks: number of ranks to compute the spectra
    :type nranks: int
    :param storage_dir: the directory where the Hankel transform are stored
    :type storage_dir: char
    :param file_save: filename to save the computation once done
    :type file_save: char
    """
    try:
        assert type(f) == type(g)
    except AssertionError:
        print("Inputs 'f' and 'g' have to share the same type!")
        print('Either numpy.ndarray or list')

    rin = r.min()
    rout = r.max()
    if type(f) == np.ndarray:
        nms = f.shape[0]
    elif type(f) == list:
        nms = f[0].shape[0]
    nr = len(r)
    # This is to determine whether data is already on cheb grid
    grid = chebgrid(nr-1, r[0], r[-1])
    diff = abs(grid-r).max()
    if diff > 1e-6:
        spectral = False
        print('I will interpolate the data on a Gauss Lobatto grid')
    else:
        spectral = True
        print('Data is already on a Gauss Lobatto grid')

    if not spectral:
        # Interpolate on half the number of points
        rnew = chebgrid(nr//2, r[0], r[-1])
        wnew = remap(weight, r, rnew)
        if type(f) == np.ndarray:
            fnew = remap(f, r, rnew)
            gnew = remap(g, r, rnew)
        elif type(f) == list:
            fnew = []
            gnew = []
            for i in range(len(f)):
                fnew.append(remap(f[i], r, rnew))
                gnew.append(remap(g[i], r, rnew))
    else:
        fnew = f
        gnew = g
        rnew = r
        wnew = weight

    if file_save is not None and os.path.exists(file_save):
        dat = np.load(file_save)
        K = dat['K']
        trans = dat['trans']
    else:
        @globalize
        def _compute_nm(nm):
            m = idx2m[nm]
            h = HankelAnnulus(rin, rout, m, len(rnew), storage_dir=storage_dir)
            if type(fnew) == np.ndarray:
                fhat = h.HT(fnew[nm, :]*np.sqrt(wnew))
                ghat = h.HT(gnew[nm, :]*np.sqrt(wnew))
                tmp = fhat.conjugate()*ghat+fhat*ghat.conjugate()
            elif type(fnew) == list:
                tmp = 0
                for i in range(len(fnew)):
                    fhat = h.HT(fnew[i][nm, :]*np.sqrt(wnew))
                    ghat = h.HT(gnew[i][nm, :]*np.sqrt(wnew))
                    tmp += fhat.conjugate()*ghat+fhat*ghat.conjugate()
            trans = np.zeros_like(h.roots)
            fac = np.ones_like(h.kr)
            vmin = jv(m, h.kr*h.rout)**2
            vmax = jv(m, h.kr*h.rin)**2
            mask = ( vmin > 1e-10 )
            fac[mask] = vmax[mask]/(vmax[mask]-vmin[mask])
            trans[:len(h.kr)] = np.pi**3 * fac * h.kr**2 * tmp.real * \
                                abs(hankel1(m, h.kr*h.rout))**2
            if m == 0:
                trans *= 0.5
            k_nm = h.roots

            return k_nm, trans

        import time
        t1 = time.time()
        with Pool(processes=nranks) as pool:
            results = pool.map(_compute_nm, range(nms))
        t2 = time.time()
        print(f'Computing time: {t2-t1:.2e} s')

        K = np.vstack([r[0] for r in results])
        trans = np.vstack([r[1] for r in results])

        if file_save is not None:
            np.savez(file_save, K=K, trans=trans)

    return K, trans

def weber_orr_spectra(f, r, weight, idx2m, nranks=1,
                      storage_dir=None, file_save=None):
    """
    This routine computes the Weber-Orr spectra. It takes a complex array
    of dimension (n_m_max, nroots) as input.

    :param f: a complex input array (first dimension is m, second dimension
              is r)
    :type f: numpy.ndarray
    :param r: the radius
    :type f: numpy.ndarray
    :param weight: a weighting function (typically the height of the container)
    :type weight: numpy.ndarray
    :param idx2m: an array which convert the index to ms
    :type idx2m: numpy.ndarray
    :param nranks: number of ranks to compute the spectra
    :type nranks: int
    :param storage_dir: the directory where the Hankel transform are stored
    :type storage_dir: char
    :param file_save: filename to save the computation once done
    :type file_save: char
    """

    rin = r.min()
    rout = r.max()
    nms = f.shape[0]
    nr = len(r)
    # This is to determine whether data is already on cheb grid
    grid = chebgrid(nr-1, r[0], r[-1])
    diff = abs(grid-r).max()
    if diff > 1e-6:
        spectral = False
    else:
        spectral = True

    if not spectral:
        # Interpolate on half the number of points
        rnew = chebgrid(nr//2, r[0], r[-1])
        fnew = remap(f, r, rnew)
        wnew = remap(weight, r, rnew)
    else:
        fnew = f
        rnew = r
        wnew = weight

    if file_save is not None and os.path.exists(file_save):
        dat = np.load(file_save)
        K = dat['K']
        Ek = dat['Ek']
    else:
        @globalize
        def _compute_nm(nm):
            m = idx2m[nm]
            h = HankelAnnulus(rin, rout, m, len(rnew), storage_dir=storage_dir)
            fhat = h.HT(fnew[nm, :]*np.sqrt(wnew))
            Ek_nm = np.zeros_like(h.roots)
            Ek_nm[:len(h.kr)] = 2 * np.pi * h.spectra(fhat)
            if m == 0:
                Ek_nm *= 0.5
            k_nm = h.roots

            return k_nm, Ek_nm

        import time
        t1 = time.time()
        with Pool(processes=nranks) as pool:
            results = pool.map(_compute_nm, range(nms))
        t2 = time.time()
        print(f'Computing time: {t2-t1:.2e} s')

        K = np.vstack([r[0] for r in results])
        Ek = np.vstack([r[1] for r in results])

        if file_save is not None:
            np.savez(file_save, K=K, Ek=Ek)

    return K, Ek


if __name__ == '__main__':
    eta = 0.35
    ri = eta/(1.-eta)
    ro = 1./(1.-eta)
    m = 11

    ht = HankelAnnulus(ri, ro, m=11, ngrid_points=256)
    f = np.sin(13*np.pi * (ht.grid-ht.rin)) * np.exp(-3*ht.grid)
    t1 = time.time()
    fhat = ht.HT(f)
    t2 = time.time()
    print(f'Timing HT: {t2-t1:.2e}')

    fback = ht.iHT(fhat)
    fig, ax = plt.subplots()
    ax.plot(ht.grid, f, lw=2, label='orig')
    ax.plot(ht.grid, fback, label='iHT')
    err = abs(f-fback).max()
    print(f'Err: {err:2e}')
    ax.set_ylim(f.min(), f.max())
    fig.legend()
    plt.show()
