from scipy import optimize
from scipy.special import jv, yn, hankel1
from .libpizza import chebgrid, intcheb, progressbar
import numpy as np
import time
import os
import pickle
from scipy.integrate import simps
import matplotlib.pyplot as plt


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
    :type m: integer
    :param nroots: nroots
    :type nroots: number of roots to be computed
    :returns: a numpy array which contains the roots
    :rtype: numpy.ndarray
    """

    def f(x):
        return (jv(m, x*rout)*yn(m, x*rin)-jv(m, x*rin)*yn(m, x*rout)) / \
                abs(hankel1(m, x*rin))

    len0 = len(zeros_m_minus_one)

    import mpmath as mp

    def fmpmath(x):
        return (mp.besselj(m, x*rout)*mp.bessely(m, x*rin) -
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
    :type m: integer
    :param nroots: nroots
    :type nroots: number of roots to be computed
    :returns: a numpy array which contains the roots
    :rtype: numpy.ndarray
    """

    def f(x):
        return (jv(m, x*rout)*yn(m, x*rin)-jv(m, x*rin)*yn(m, x*rout)) / \
                abs(hankel1(m, x*rin))

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
    y = f(x)
    mask = np.isnan(y)
    improve_prec = False
    if len(x[mask]) > 0:
        improve_prec = True

    if improve_prec:
        print('Precision of scipy is not enough, using mpmath')

        import mpmath as mp

        def fmpmath(x):
            return (mp.besselj(m, x*rout)*mp.bessely(m, x*rin) -
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
        xmax = x[0]

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


class HankelAnnulus:

    def __init__(self, rin, rout, m, nroots, storage_dir=None, grid_type='GL',
                 verbose=False):
        """
        :param rin: inner radius
        :type rin: float
        :param rout: outer radius
        :type rout: float
        :param m: azimuthal wavenumber
        :type m: integer
        :param nroots: nroots
        :type nroots: number of roots to be computed
        :param storage_dir: a directory when one case store the Hankel
                            transforms
        :type storage_dir: char
        :param grid_type: a parameter to select the collocation grid: 'GL'
                          stands for Gauss-Lobatto and allows a more accurate
                          spectral integration. Possible alternatives are
                          'EQUI' for an equidistant grid, and 'ZER' for a grid
                          built upon zeros of the transcendental equation.
        :type grid_type: char
        :param verbose: a boolean to trigger some printouts of infos
        :type verbose: bool
        """
        self.grid_type = grid_type
        compute = True
        store = False
        if storage_dir is not None:
            filename = f'WeberOrr_{rin/rout:.2f}_{nroots:04d}_{m:04d}.pickle'
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
            self.nroots = nroots

            # Get the roots
            if storage_dir is not None:
                filename = \
                     f'WeberOrr_{rin/rout:.2f}_{nroots:04d}_{m-1:04d}.pickle'
                file = os.path.join(storage_dir, filename)
                if os.path.exists(file):
                    with open(file, 'rb') as fi:
                        tmp = pickle.load(fi)
                        roots_m_minus_1, grid_m_minus_1, k_m_minus_1 = \
                            pickle.load(fi)
                    self.roots = jnyn_zeros_guess(rin, rout, m,
                                                  roots_m_minus_1,
                                                  nroots+1)
                else:
                    self.roots = jnyn_zeros(rin, rout, m, nroots+1)
            else:
                self.roots = jnyn_zeros(rin, rout, m, nroots+1)
            roots_last = self.roots[-1]
            self.roots = self.roots[:-1]

            if self.grid_type == 'GL':
                # Gauss-Lobatto grid
                self.grid = chebgrid(nroots-1, rin, rout)
            elif self.grid_type == 'EQUI':
                # Equi-distant grid
                self.grid = np.linspace(rin, rout, nroots+2)
                self.grid = self.grid[1:-1]
            else:
                # Remap the roots between rin and rout
                self.grid = self.roots/roots_last * (rout-rin) + rin

            self.kr = self.roots  # / 2. / np.pi

            # compute the Hankel matrix
            # self.hankel_mat()

            # Inverse of Hankel matrix
            # self.Ainv = np.linalg.inv(self.A)
            # self.Ainv = np.zeros((self.nroots, self.nroots), np.float64)

            self.compute_kernels()

            if store:
                self._store(storage_dir)

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
            self.rin,  self.rout, self.m, self.nroots = pickle.load(f)
            self.roots, self.grid, self.kr = pickle.load(f)
            self.kernels = pickle.load(f)[0]

    def _store(self, dir):
        """
        This routine is used to store the inverse matrix and the roots

        :param dir: a directory when one case store the Hankel transforms
        :type dir: char
        """

        filename = f'WeberOrr_{self.rin/self.rout:.2f}_{self.nroots:04d}_{self.m:04d}.pickle'
        file = os.path.join(dir, filename)

        with open(file, 'wb') as f:
            pickle.dump([self.rin, self.rout, self.m, self.nroots], f)
            pickle.dump([self.roots, self.grid, self.kr], f)
            pickle.dump([self.kernels], f)

    def compute_kernels(self):
        """
        This routine is used to store the kernels employed in the inverse
        Weber-Orr transform.
        """
        self.kernels = np.zeros((len(self.grid), len(self.grid)), np.float64)

        def f(root, grid):
            return (jv(self.m, root*grid)*yn(self.m, root*self.rin) -
                    jv(self.m, root*self.rin)*yn(self.m, root*grid)) / \
                    abs(hankel1(self.m, root*self.rin))

        def fmpmath(root, grid):
            import mpmath as mp
            return (mp.besselj(self.m, root*grid)*mp.bessely(self.m, root*self.rin) -
                    mp.besselj(self.m, root*self.rin)*mp.bessely(self.m, root*grid)) / \
                    abs(mp.hankel1(self.m, root*self.rin))

        for k, root in enumerate(self.roots):
            dat = f(root, self.grid)
            mask = np.isnan(dat)
            if len(dat[mask]) > 0:
                vals = np.zeros(len(dat[mask]), np.float64)
                for i in range(len(self.grid[mask])):
                    gr = self.grid[mask][i]
                    val = fmpmath(root, gr)
                    vals[i] = val
                dat[mask] = vals
            self.kernels[k, :] = dat

    def hankel_mat(self):
        """
        This subroutine constructs the Hankel matrix from the zeroes and
        the collocation grid
        """
        self.A = np.zeros((len(self.grid), len(self.grid)), np.float64)
        fac = self.roots**2 * jv(self.m, self.roots*self.rout)**2 / \
              (jv(self.m, self.roots*self.rin)**2 - 
               jv(self.m, self.roots*self.rout)**2)
        norm = abs(hankel1(self.m, self.roots*self.rin))
        for k, rk in enumerate(self.grid):
            kernel = jv(self.m, self.roots*rk)*yn(self.m, self.roots*self.rin) - \
                     jv(self.m, self.roots*self.rin)*yn(self.m, self.roots*rk)
            kernel *= norm
            self.A[k, :] = fac * kernel * np.pi**2 / 2.

    def spectra(self, fhat):
        """
        :param fhat: the function defined in spectral space
        :type fhat: numpy.ndarray
        :returns: the energy content as a function of the order
        :rtype: numpy.ndarray
        """
        fac = np.zeros_like(self.roots)
        for k, root in enumerate(self.roots):
            vmin = jv(self.m, root*self.rout)**2
            vmax = jv(self.m, root*self.rin)**2

            if vmin > 1e-10:  # to prevent round-off errors
                fac[k] = root**2 * vmax * abs(hankel1(self.m, root*self.rout))**2 / \
                         (vmax - vmin)
            else:
                fac[k] = root**2 * abs(hankel1(self.m, root*self.rout))**2
        En = 0.5 * np.pi**2 * fac * abs(fhat)**2

        return En

    def HT(self, f):
        """
        :param f: the function defined in physical space
        :type f: numpy.ndarray
        :returns: the function in spectral space
        :rtype: numpy.ndarray
        """
        # if self.m <= 10 and hasattr(self, 'Ainv'):
        #     fhat = self.Ainv@f
        # else:
        fhat = self.HTsimps(f)
        # from scipy.sparse.linalg import gmres
        # fhat, info = gmres(self.A, f, x0=fhat_guess)

        return fhat

    def HTsimps(self, f):
        """
        This is a Hankel transform based on integration

        :param f: the function defined in physical space
        :type f: numpy.ndarray
        :returns: the function in spectral space
        :rtype: numpy.ndarray
        """
        fhat = np.zeros_like(f)
        if self.grid_type == 'GL':
            for k in range(self.nroots):
                fhat[k] = intcheb(self.kernels[k, :] * f * self.grid)
                # fhat[k] = simps(self.kernels[k, :] * f * self.grid, self.grid)
        else:
            for k in range(self.nroots):
                fhat[k] = simps(self.kernels[k, :] * f * self.grid, self.grid)

        return fhat

    def iHT(self, fhat):
        """
        :param fhat: the function defined in spectral space
        :type fhat: numpy.ndarray
        :returns: the function in physical space
        :rtype: numpy.ndarray
        """
        # if hasattr(self, 'A'):
        #     f = self.A@fhat
        # else:
        if self.m > 30:
            import sys
            sys.exit('Inverse transform is not accurate beyond m=30')
        else:
            f = self.iHTsum(fhat)

        return f

    def iHTsimps(self, fhat):
        """
        This is an inverse Hankel transform based on integrals (for debug only)

        :param fhat: the function defined in spectral space
        :type fhat: numpy.ndarray
        :returns: the function in physical space
        :rtype: numpy.ndarray
        """
        f = np.zeros_like(fhat)
        norm = abs(hankel1(self.m, self.roots*self.rin))
        for k, rk in enumerate(self.grid):
            kernel = jv(self.m, self.roots*rk)*yn(self.m, self.roots*self.rin) - \
                     jv(self.m, self.roots*self.rin)*yn(self.m, self.roots*rk)
            kernel /= norm
            f[k] = simps(kernel * fhat * self.roots, self.roots)

        return f

    def iHTsum(self, fhat):
        """
        :param fhat: the function defined in spectral space
        :type fhat: numpy.ndarray
        :returns: the function in physical space
        :rtype: numpy.ndarray
        """
        f = np.zeros_like(fhat)
        fac = self.roots**2 * jv(self.m, self.roots*self.rout)**2 / \
                (jv(self.m, self.roots*self.rin)**2 -
                 jv(self.m, self.roots*self.rout)**2)
        norm = abs(hankel1(self.m, self.roots*self.rin))
        for k, rk in enumerate(self.grid):
            kernel = jv(self.m, self.roots*rk)*yn(self.m, self.roots*self.rin) - \
                     jv(self.m, self.roots*self.rin)*yn(self.m, self.roots*rk)
            kernel *= norm
            f[k] = np.sum(fac*fhat*kernel) * np.pi**2/2.

        return f


def weber_orr_spectra(f, r, weight, idx2m, storage_dir=None, file_save=None):
    """
    This routine computes the Weber-Orr spectra. It takes a complex array
    of dimension (n_m_max, nroots) as input.

    :param f: a a complex input array (first dimension is m
    :type f: numpy.ndarray
    :param r: the radius
    :type f: numpy.ndarray
    :param weight: a weighting function (typically the height of the container)
    :type weight: numpy.ndarray
    :param idx2m: an array which convert the index to ms
    :type idx2m: numpy.ndarray
    :param storage_dir: the directory where the Hankel transform are stored
    :type storage_dir: char
    :param file_save: filename to save the computation once done
    :type file_save: char
    """
    rin = r.min()
    rout = r.max()
    nroots = len(r)
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
                h = HankelAnnulus(rin, rout, m, nroots,
                                  storage_dir=storage_dir)
            else:
                h = HankelAnnulus(rin, rout, m, nroots)
            fhat = h.HT(f[idx, :]*np.sqrt(weight))
            if m == 0:
                Ek[idx, :] = np.pi * h.spectra(fhat)
            else:
                Ek[idx, :] = 2.*np.pi * h.spectra(fhat)
            K[idx, :] = h.kr

    if file_save is not None:
        with open(file_save, 'wb') as fi:
            pickle.dump([K, Ek], fi)

    return K, Ek


if __name__ == '__main__':

    eta = 0.35
    ri = eta/(1.-eta)
    ro = 1./(1.-eta)

    nroots = 10
    m = 61

    def f(x):
        return jv(m, x*ro)*yn(m, x*ri)-jv(m, x*ri)*yn(m, x*ro)

    roots = np.zeros(768, np.float64)
    for m in range(768):
        ht = HankelAnnulus(ri, ro, m, 1)
        roots[m] = ht.roots[0]

    plt.plot(np.diff(roots))
    plt.show()

    ht = HankelAnnulus(ri, ro, m, nroots)

    x = np.linspace(0.66*m, 5.66*m, 10*m)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, f(x))
    ax.plot(ht.roots, np.zeros_like(ht.roots), ls='None',
            marker='o', mfc='None')
    ax.set_xlabel('x')
    ax.set_xlim(0.66*m, 5.66*m)
    ax.set_ylim(-f(5.66*m)*10, f(5.66*m)*10)
    fig.tight_layout()

    nroots = 2048
    m_max = 2049

    for m in range(m_max):
        ht = HankelAnnulus(ri, ro, m, nroots,
                           storage_dir='/home/gastine/hankel_mats')

    """
    #f = np.sin(13*np.pi * (ht.grid-ht.rin)/ht.rin) * np.exp(-3*ht.grid)
    f = np.sin(13*np.pi * (ht.grid-ht.rin)) * np.exp(-3*ht.grid)

    t1 = time.time()
    fhat = ht.HTsimps(f)
    fhat1 = ht.HT(f)
    t2 = time.time()
    print(f'Timing HT: {t2-t1:.2e}')


    t1 = time.time()
    fback = ht.iHTsimps(fhat1)
    t2 = time.time()
    print(f'Timing iHT (simps): {t2-t1:.2e}')
    t1 = time.time()
    fback1 = ht.iHTsum(fhat1)
    t2 = time.time()
    print(f'Timing iHT (series): {t2-t1:.2e}')
    t1 = time.time()
    fback2 = ht.iHT(fhat1)
    t2 = time.time()
    print(f'Timing iHT (dgem): {t2-t1:.2e}')

    plt.plot(ht.grid, f, lw=2, label='orig')
    plt.plot(ht.grid, fback, label='iHT integ')
    print(abs(f-fback).max())
    plt.plot(ht.grid, fback1, label='iHT sum')
    print(abs(f-fback1).max())
    plt.plot(ht.grid, fback2, label='iHT mat')
    print(abs(f-fback2).max())
    plt.ylim(f.min(), f.max())
    plt.legend()
    plt.show()

    plt.loglog(ht.kr, ht.spectra(fhat))
    plt.loglog(ht.kr, ht.spectra(fhat1))
    """
    plt.show()
