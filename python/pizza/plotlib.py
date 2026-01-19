# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

mplMaj, mplMin, _ = np.int32(matplotlib.__version__.split('.'))

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


def equatContour(data, radius, minc=1, label=None, levels=64,
                 cm='seismic', normed=None, vmax=None, vmin=None, cbar=True,
                 title=True, normRad=False, deminc=True, bounds=True,
                 lines=False, linewidths=0.5, pcolor=False, rasterized=False,
                 fig=None, ax=None, shading='flat'):
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
    :param title: display the title of the figure when set to True
    :type title: bool
    :param cbar: display the colorbar when set to True
    :type cbar: bool
    :param vmax: maximum value of the contour levels
    :type vmax: float
    :param vmin: minimum value of the contour levels
    :type vmin: float
    :param normed: when set to True, the colormap is centered around zero.
                   Default is None, it tries to find it by itself.
    :type normed: bool
    :param deminc: a logical to indicate if one wants do get rid of the
                   possible azimuthal symmetry
    :type deminc: bool
    :param bounds: a boolean to determine if one wants to plot the limits
                   of the domain (True by default)
    :type bounds: bool
    :param lines: when set to True, over-plot solid lines to highlight
                  the limits between two adjacent contour levels
    :type lines: bool
    :param linewidths: the thickness of the solid lines, whenever plotted
    :type linewidths: float
    :param pcolor: when set to True, use pcolormesh instead of contourf
    :type pcolor: bool
    :param rasterized: when set to True, the rasterization for vector graphics
                       is turned on
    :type rasterized: bool
    :param fig: a pre-existing figure (if needed)
    :type fig: matplotlib.figure.Figure
    :param ax: a pre-existing axis
    :type ax: matplotlib.axes._subplots.AxesSubplot
    :param shading: shading if pcolormesh is employed (can be 'flat', 'gouraud')
    :type shading: str
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
        if title and label is not None:
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

    if normed is None:
        if abs(data.min()) < 1e-8:
            normed = False
        else:
            if data.max() > 0 and data.min() < 0:
                normed = True
            else:
                normed = False

        if cm == 'seismic' and not normed:
            cm = 'viridis'

    cmap = plt.get_cmap(cm)
    if vmax is not None or vmin is not None:
        normed = False
        cs = np.linspace(vmin, vmax, levels)
        if pcolor:
            if shading == 'flat':
                data = data[:-1, :-1]
            im = ax.pcolormesh(xx, yy, data, cmap=cmap, antialiased=True,
                               shading=shading, vmax=vmax, vmin=vmin,
                               rasterized=rasterized, edgecolors='face')
        else:
            im = ax.contourf(xx, yy, data, cs, cmap=cmap, extend='both')
        if lines:
            ax.contour(xx, yy, data, cs, colors=['k'], linewidths=linewidths,
                       extend='both', linestyles=['-'])
    else:
        if not normed:
            cs = levels
        else:
            vmax = max(abs(data.max()), abs(data.min()))
            vmin = -vmax
            cs = np.linspace(vmin, vmax, levels)
        if pcolor:
            if shading == 'flat':
                data = data[:-1, :-1]
            if normed:
                im = ax.pcolormesh(xx, yy, data, cmap=cmap, antialiased=True,
                                   shading=shading, vmax=vmax, vmin=vmin,
                                   rasterized=rasterized, edgecolors='face')
            else:
                im = ax.pcolormesh(xx, yy, data, cmap=cmap, antialiased=True,
                                   shading=shading, rasterized=rasterized,
                                   edgecolors='face')
        else:
            im = ax.contourf(xx, yy, data, cs, cmap=cmap, antialiased=True)
        if lines:
            ax.contour(xx, yy, data, cs, colors=['k'], linewidths=linewidths,
                       linestyles=['-'])

    if bounds:
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

    eps = 1e-4
    if xx.min() < -eps:
        ax.set_xlim(1.01*xx.min(), 1.01*xx.max())
    elif abs(xx.min()) < eps :
        ax.set_xlim(xx.min()-0.01, 1.01*xx.max())
    else:
        ax.set_xlim(0.99*xx.min(), 1.01*xx.max())
    if yy.min() < -eps:
        ax.set_ylim(1.01*yy.min(), 1.01*yy.max())
    elif abs(yy.min()) < eps:
        ax.set_ylim(yy.min()-0.01, 1.01*yy.max())
    else:
        ax.set_ylim(0.99*yy.min(), 1.01*yy.max())
    ax.axis('off')

    # Add the colorbar at the right place
    pos = ax.get_position()
    l, b, w, h = pos.bounds
    if cbar:
        if title and label is not None:
            cax = fig.add_axes([0.85, 0.46-0.7*h/2., 0.03, 0.7*h])
        else:
            cax = fig.add_axes([0.85, 0.5-0.7*h/2., 0.03, 0.7*h])
        mir = fig.colorbar(im, cax=cax)

    #To avoid white lines on pdfs
    if not pcolor:
        for c in ax.collections:
            c.set_edgecolor('face')
            if mplMaj >= 3 and mplMin < 9:
                if rasterized:
                    c.set_rasterized(True)

    return fig, xx, yy, im


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
