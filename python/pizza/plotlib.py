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
        mask = np.where(dat>=vmax, 1, 0)
        dat = dat*(mask == 0) + vmax*(mask == 1)
        normed = False
    if vmin is not None:
        mask = np.where(dat<=vmin, 1, 0)
        dat = dat*(mask == 0) + vmin*(mask == 1)
        normed = False
    return dat


def equatContour(data, radius, minc=1, label=None, levels=65,
                 cm='seismic', normed=True, vmax=None, vmin=None, cbar=True,
                 tit=True, normRad=False, deminc=True):
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
    """

    nphi, ntheta = data.shape

    if deminc:
        phi = np.linspace(0., 2.*np.pi, nphi)
    else:
        phi = np.linspace(0., 2.*np.pi/minc, nphi)
    rr, pphi = np.meshgrid(radius, phi)
    xx = rr * np.cos(pphi)
    yy = rr * np.sin(pphi)

    if normRad: # Normalise each radius
        maxS = np.sqrt(np.mean(data**2, axis=0))
        data[:, maxS!=0.] /= maxS[maxS!=0.]

    if tit and label is not None:
        if cbar:
            fig = plt.figure(figsize=(6.5,5.5))
            ax = fig.add_axes([0.01, 0.01, 0.76, 0.9])
        else:
            fig = plt.figure(figsize=(5,5.5))
            ax = fig.add_axes([0.01, 0.01, 0.98, 0.9])
        ax.set_title(label, fontsize=24)
    else:
        if cbar:
            fig = plt.figure(figsize=(6.5,5))
            ax = fig.add_axes([0.01, 0.01, 0.76, 0.98])
        else:
            fig = plt.figure(figsize=(5, 5))
            ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])

    cmap = plt.get_cmap(cm)
    if vmax is not None or vmin is not None:
        normed = False
        cs = np.linspace(vmin, vmax, levels)
        im = ax.contourf(xx, yy, data, cs, cmap=cmap, extend='both')
    else:
        cs = levels
        im = ax.contourf(xx, yy, data, cs, cmap=cmap)
        #im = ax.pcolormesh(xx, yy, data, cmap=cmap, antialiased=True)

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
        mir = fig.colorbar(im, cax=cax)

    # Normalise data 
    if normed:
        im.set_clim(-max(abs(data.max()), abs(data.min())),
                     max(abs(data.max()), abs(data.min())))

    return fig, xx, yy
