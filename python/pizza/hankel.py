from scipy import optimize
from scipy.special import jv, yn
import numpy as np
import time
import os
import pickle
from scipy.integrate import simps
import matplotlib.pyplot as plt


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
        return jv(m, x*rout)*yn(m, x*rin)-jv(m, x*rin)*yn(m, x*rout)

    if m < 10:
        npoints = 1024
        x = np.linspace(0.1, 100, npoints)
    else:
        npoints = 10*m
        # Spacing between two m for first roots is like lambda(m+1)-lambda(m)\sim0.65
        x = np.linspace(int(0.65*m), 2*m, npoints)
    signs = np.sign(f(x))
    nroot = 0
    roots = np.zeros((nroots), np.float64)
    for i in range(npoints-1):
        if signs[i] + signs[i+1] == 0:
            sol = optimize.root_scalar(f, bracket=[x[i], x[i+1]],
                                       method='brentq')
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

    #plt.plot(x, f(x))
    #plt.ylim(-10*f(x[-1]), 10*f(x[-1]))
    #plt.plot(roots, f(roots), ls='None', marker='o', mfc='None')
    #plt.show()

    #print(f(roots))

    return roots


class HankelAnnulus:

    def __init__(self, rin, rout, m, nroots, storage_dir=None):
        """
        :param rin: inner radius
        :type rin: float
        :param rout: outer radius
        :type rout: float
        :param m: azimuthal wavenumber
        :type m: integer
        :param nroots: nroots
        :type nroots: number of roots to be computed
        :param storage_dir: a directory when one case store the Hankel transforms
        :type storage_dir: char
        """

        compute = True
        store = False
        if storage_dir is not None:
            filename = f'Hankel_{rin/rout:.2f}_{nroots:04d}_{m:04d}.pickle'
            file = os.path.join(storage_dir, filename)
            if os.path.exists(file):
                self._read(file)
                compute = False
            else:
                store = True

        if compute:
            self.rin = rin
            self.rout = rout
            self.m = m
            self.nroots = nroots

            # Get the roots
            self.roots = jnyn_zeros(rin, rout, m, nroots+1)
            roots_last = self.roots[-1]
            self.roots = self.roots[:-1]
            # Remap the roots between rin and rout
            self.grid = self.roots/roots_last * (rout-rin) + rin
            #self.grid = np.linspace(rin, rout, nroots+2)
            #self.grid = self.grid[1:-1]

            self.kr = self.roots * (self.rout-self.rin) / 2. / np.pi

            # compute the Hankel matrix
            self.hankel_mat()

            # Inverse of Hankel matrix
            self.Ainv = np.linalg.inv(self.A)

            if store:
                self._store(storage_dir)

    def _read(self, file):
        """
        This routine is used to read the inverse matrix and the roots

        :param file: the pickle file which contains the structure
        :type file: char
        """
        with open(file, 'rb') as f:
            print(f'Reading {file}')
            self.rin,  self.rout, self.m, self.nroots = pickle.load(f)
            self.roots, self.grid, self.kr = pickle.load(f)
            self.Ainv = pickle.load(f)[0]

    def _store(self, dir):
        """
        This routine is used to store the inverse matrix and the roots

        :param dir: a directory when one case store the Hankel transforms
        :type dir: char
        """

        filename = f'Hankel_{self.rin/self.rout:.2f}_{self.nroots:04d}_{self.m:04d}.pickle'
        file = os.path.join(dir, filename)

        with open(file, 'wb') as f:
            pickle.dump([self.rin, self.rout, self.m, self.nroots], f)
            pickle.dump([self.roots, self.grid, self.kr], f)
            pickle.dump([self.Ainv], f)

    def hankel_mat(self):
        """
        This subroutine constructs the Hankel matrix from the zeroes and
        the collocation grid
        """
        self.A = np.zeros((len(self.grid), len(self.grid)), np.float64)
        for k, rk in enumerate(self.grid):
            kernel = jv(self.m, self.roots*rk)*yn(self.m, self.roots*self.rin) - \
                     jv(self.m, self.roots*self.rin)*yn(self.m, self.roots*rk)
            fac = self.roots**2 * jv(self.m, self.roots*self.rout)**2 / \
                  (jv(self.m, self.roots*self.rin)**2 - jv(self.m, self.roots*self.rout)**2)
            self.A[k, :] = fac * kernel * np.pi**2 / 2.

    def spectra(self, fhat):
        """
        :param fhat: the function defined in spectral space
        :type fhat: numpy.ndarray
        :returns: the energy content as a function of the order
        :rtype: numpy.ndarray
        """
        fac = self.roots**2 * jv(self.m, self.roots*self.rout)**2 / \
              (jv(self.m, self.roots*self.rin)**2 - jv(self.m, self.roots*self.rout)**2)
        En = 0.5 * np.pi**2 * fac * abs(fhat)**2

        return En

    def HT(self, f):
        """
        :param f: the function defined in physical space
        :type f: numpy.ndarray
        :returns: the function in spectral space
        :rtype: numpy.ndarray
        """
        if self.m <= 10:
            fhat = self.Ainv@f
        else:
            fhat = self.HTsimps(f)
        #from scipy.sparse.linalg import gmres
        #fhat, info = gmres(self.A, f, x0=fhat_guess)

        return fhat

    def HTsimps(self, f):
        """
        This is a Hankel transform based on integration (only for debug)

        :param f: the function defined in physical space
        :type f: numpy.ndarray
        :returns: the function in spectral space
        :rtype: numpy.ndarray
        """
        fhat = np.zeros_like(f)
        for k, root in enumerate(self.roots):
            kernel = jv(self.m, root*self.grid)*yn(self.m, root*self.rin) - \
                     jv(self.m, root*self.rin)*yn(self.m, root*self.grid)
            if k == 2:
                plt.plot(self.grid, kernel * f.real * self.grid)
                plt.show()
            fhat[k] = simps(kernel * f * self.grid, self.grid)

        return fhat

    def iHT(self, fhat):
        """
        :param fhat: the function defined in spectral space
        :type fhat: numpy.ndarray
        :returns: the function in physical space
        :rtype: numpy.ndarray
        """
        if hasattr(self, 'A'):
            f = self.A@fhat
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
        for k, rk in enumerate(self.grid):
            kernel = jv(self.m, self.roots*rk)*yn(self.m, self.roots*self.rin) - \
                     jv(self.m, self.roots*self.rin)*yn(self.m, self.roots*rk)
            kernel /= jv(m, self.roots*self.rin)**2 + yn(m, self.roots*self.rin)**2
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
        for k, rk in enumerate(self.grid):
            kernel = jv(self.m, self.roots*rk)*yn(self.m, self.roots*self.rin) - \
                     jv(self.m, self.roots*self.rin)*yn(self.m, self.roots*rk)
            fac = self.roots**2 * jv(self.m, self.roots*self.rout)**2 / \
                 (jv(self.m, self.roots*self.rin)**2 - jv(self.m, self.roots*self.rout)**2)
            f[k] = np.sum(fac*fhat*kernel) * np.pi**2/2.

        return f


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
    ax.plot(ht.roots, np.zeros_like(ht.roots), ls='None', marker='o', mfc='None')
    ax.set_xlabel('x')
    ax.set_xlim(0.66*m, 5.66*m)
    ax.set_ylim(-f(5.66*m)*10, f(5.66*m)*10)
    fig.tight_layout()

    nroots = 2048
    m_max = 2049

    for m in range(m_max):
        ht = HankelAnnulus(ri, ro, m, nroots, storage_dir='/home/gastine/hankel_mats')

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
