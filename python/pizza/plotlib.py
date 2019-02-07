# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt


def cut(dat, vmax=None, vmin=None):
    """
    This functions truncates the values of an input array that are beyond
    vmax or below vmin and replace them by vmax and vmin, respectively.

    >>> # Keep only values between -1e3 and 1e3
    >>> datNew = cut(dat, vmin=-1e3, vmax=1e3)

    :param dat: an input array
    :type dat: numpy.ndarray
    :param vmax: maximum upper bound
    :type vmax: float
    :param vmin: minimum lower bound
    :type vmin: float
    :returns: an array where the values >=vmax have been replaced by vmax
              and the values <=vmin have been replaced by vmin
    :rtype: numpy.ndarray
    """
    if vmax is not None:
        mask = np.where(dat >= vmax, 1, 0)
        dat = dat*(mask == 0) + vmax*(mask == 1)
    if vmin is not None:
        mask = np.where(dat <= vmin, 1, 0)
        dat = dat*(mask == 0) + vmin*(mask == 1)

    return dat


def equatContour(data, radius, minc=1, label=None, levels=65,
                 cm='seismic', normed=True, vmax=None, vmin=None, cbar=True,
                 tit=True, normRad=False, deminc=True, stream=False,
                 contourLines=False, fig=None, ax=None):
    """
    Plot the equatorial cut of a given field

    :param data: the input data (an array of size (nphi,nr))
    :type data: numpy.ndarray
    :param radius: the input radius
    :type radius: numpy.ndarray
    :param minc: azimuthal symmetry
    :type minc: int
    :param label: the name of the input physical quantity you want to
                  display
    :type label: str
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
    :param deminc: a logical to indicate if one wants do get rid of the
                   possible azimuthal symmetry
    :type deminc: bool
    :param contourLines: a switch to overplot the contour lines
    :type contourLines: bool
    :param fig: a pre-existing figure (if needed)
    :type fig: matplotlib.figure.Figure
    :param ax: a pre-existing axis
    :type ax: matplotlib.axes._subplots.AxesSubplot
    """

    nphi, ntheta = data.shape

    if deminc:
        phi = np.linspace(0., 2.*np.pi, nphi)
    else:
        phi = np.linspace(0., 2.*np.pi/minc, nphi)
    rr, pphi = np.meshgrid(radius, phi)
    xx = rr * np.cos(pphi)
    yy = rr * np.sin(pphi)

    if normRad:  # Normalise each radius
        maxS = np.sqrt(np.mean(data**2, axis=0))
        data[:, maxS != 0.] /= maxS[maxS != 0.]

    if fig is None and ax is None:
        if tit and label is not None:
            if cbar:
                fig = plt.figure(figsize=(6.5, 5.5))
                ax = fig.add_axes([0.01, 0.01, 0.76, 0.9])
            else:
                fig = plt.figure(figsize=(5, 5.5))
                ax = fig.add_axes([0.01, 0.01, 0.98, 0.9])
            ax.set_title(label, fontsize=24)
        else:
            if cbar:
                fig = plt.figure(figsize=(6.5, 5))
                ax = fig.add_axes([0.01, 0.01, 0.76, 0.98])
            else:
                fig = plt.figure(figsize=(5, 5))
                ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])

        if not deminc:
            if minc == 2:
                w, h = fig.get_size_inches()
                fig.set_size_inches((2*w, h))

    cmap = plt.get_cmap(cm)
    if vmax is not None or vmin is not None:
        normed = False
        cs = np.linspace(vmin, vmax, levels)
        im = ax.contourf(xx, yy, data, cs, cmap=cmap, extend='both')
        if contourLines:
            ax.contour(xx, yy, data, cs, colors=['k'], linewidths=[0.5],
                       extend='both')
    else:
        cs = levels
        im = ax.contourf(xx, yy, data, cs, cmap=cmap)
        if contourLines:
            ax.contour(xx, yy, data, cs, colors=['k'], linewidths=[0.5])
        # im = ax.pcolormesh(xx, yy, data, cmap=cmap, antialiased=True)

    ax.plot(radius[0]*np.cos(phi), radius[0]*np.sin(phi), 'k-', lw=1.5)
    ax.plot(radius[-1]*np.cos(phi), radius[-1]*np.sin(phi), 'k-', lw=1.5)

    if not deminc and minc > 1:
        ax.plot(radius, np.zeros_like(radius), 'k-', lw=1.5)
        xa = radius[-1]*np.cos(2.*np.pi/minc)
        ya = radius[-1]*np.sin(2.*np.pi/minc)
        xb = radius[0]*np.cos(2.*np.pi/minc)
        x = np.linspace(xa, xb, 32)
        y = np.tan(2.*np.pi/minc)*(x-xa)+ya
        ax.plot(x, y, 'k-', lw=1.5)
        ax.plot(radius, np.zeros_like(radius), 'k-', lw=1.5)

    if xx.min() < 0:
        ax.set_xlim(1.01*xx.min(), 1.01*xx.max())
    elif xx.min() == 0.:
        ax.set_xlim(xx.min()-0.01, 1.01*xx.max())
    else:
        ax.set_xlim(0.99*xx.min(), 1.01*xx.max())
    if yy.min() < 0:
        ax.set_ylim(1.01*yy.min(), 1.01*yy.max())
    elif yy.min() == 0.:
        ax.set_ylim(yy.min()-0.01, 1.01*yy.max())
    else:
        ax.set_ylim(0.99*yy.min(), 1.01*yy.max())
    ax.axis('off')

    # Add the colorbar at the right place
    pos = ax.get_position()
    l, b, w, h = pos.bounds
    if cbar:
        if tit and label is not None:
            cax = fig.add_axes([0.85, 0.46-0.7*h/2., 0.03, 0.7*h])
        else:
            cax = fig.add_axes([0.85, 0.5-0.7*h/2., 0.03, 0.7*h])
        fig.colorbar(im, cax=cax)

    # Normalise data
    if normed:
        im.set_clim(-max(abs(data.max()), abs(data.min())),
                    max(abs(data.max()), abs(data.min())))

    return fig, xx, yy


def spec2D(field_m, n_cheb_max=None):
    """
    This subroutine allows to display one 2-D spectrum in the Fourier/Chebyshev
    space

    :param field_m: the input array in Fourier/Chebyshev space
    :type field_m: numpy.ndarray
    :param n_cheb_max: the maximum Chebyshev degree (in case one does not want
                       to display up to n_r_max)
    :type n_cheb_max: int

    >>> f = PizzaFields(tag='test', ivar=1)
    >>> om_hat = costf(f.vortz_m)
    >>> spec2D(om_hat)
    """
    n_m_max, n_r_max = field_m.shape
    if n_cheb_max is None:
        n_cheb_max = n_r_max

    x = np.arange(n_m_max-1)+1
    x = np.log10(x)
    y = np.arange(n_cheb_max)+1
    y = np.log10(y)

    xx, yy = np.meshgrid(y, x)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    eps = 1e-15
    eps = 1e-15
    dat = np.log10(abs(field_m[1:, :n_cheb_max])+eps)
    im = ax.plot_surface(xx, yy, dat, cmap=plt.get_cmap('magma'),
                         antialiased=False, linewidth=0.)
    ax.set_ylabel(r'$\log_{10}(m)$')
    ax.set_xlabel(r'$\log_{10}(N_c+1)$')
    ax.set_zlabel(r'$\log_{10}(|\hat{f}_m|)$')
    vmax = dat.max()
    vmin = dat.min()
    ax.set_zlim(vmin+0.5, vmax+0.5)
    ax.view_init(30, 60)
    fig.colorbar(im)
    fig.tight_layout()
