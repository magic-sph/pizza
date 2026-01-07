import numpy as np
import matplotlib.pyplot as plt
from .libpizza import scanDir
from .log import PizzaSetup
import os
import re


class PizzaTransfer(PizzaSetup):
    """
    This module is used to handle the output files energy_transfer.TAG
    """

    def __init__(self, datadir='.', tag=None, endian='l', iplot=True, cm='seismic'):
        """
        :param datadir: working directory
        :type datadir: str
        :param tag: a specific trailing tag, default is None
        :type tag: str
        :param endian: endianness of the file ('B' or 'l')
        :type endian: str
        :param iplot: boolean to toggle the display (False by default)
        :type iplot: bool
        :param cm: name of the colormap
        :type cm: str
        """
        pattern = os.path.join(datadir, 'log.*')
        logFiles = scanDir(pattern)

        if tag is not None:
            pattern = os.path.join(datadir, f'energy_transfer.{tag}')
            files = scanDir(pattern)

            #  Either the log.tag directly exists and the setup is easy
            #  to obtain
            if os.path.exists(os.path.join(datadir, f'log.{tag}')):
                PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml=f'log.{tag}')
            #  Or the tag is a bit more complicated and we need to find
            #  the corresponding log file
            else:
                mask = re.compile(fr'{datadir}/energy_transfer\.(.*)')
                if mask.match(files[-1]):
                    ending = mask.search(files[-1]).groups(0)[0]
                    pattern = os.path.join(datadir, f'log.{ending}')
                    if os.path.exists(pattern):
                        PizzaSetup.__init__(self, datadir=datadir, quiet=True,
                                            nml=f'log.{ending}')
            filename = files[-1]
        else:
            if len(logFiles) != 0:
                PizzaSetup.__init__(self, quiet=True, nml=logFiles[-1])
                name = f'energy_transfer.{self.tag}'
                filename = os.path.join(datadir, name)
                print(f'reading {filename}')
            else:
                dat = scanDir('energy_transfer.*')
                filename = dat[-1][1]
                print(f'reading {filename}')

        self._read(filename, endian)

        if iplot:
            self.plot(cm=cm)

    def _read(self, filename, endian):
        """
        Read one energy_transfer.TAG file

        :param filename: full name of the file
        :type filename: str
        :param endian: endianness of the file ('B' or 'l')
        :type endian: str
        """

        with open(filename, 'rb') as fi:
            version = np.fromfile(fi, dtype=np.int32, count=1)[0]
            self.ra, self.ek, self.pr, self.radratio, self.sc, self.raxi = \
               np.fromfile(fi, dtype=np.float64, count=6)
            self.n_r_max, self.n_m_max, self.m_max, self.minc, self.n_phi_max, \
               self.n_m_transfer = np.fromfile(fi, dtype=np.int32, count=6)
            self.idx2m = np.fromfile(fi, dtype=np.int32, count=self.n_m_max)

            self.trans = np.fromfile(fi, dtype=np.float64)
            self.trans = self.trans.reshape((self.n_m_transfer,
                                             self.n_m_transfer))

    def plot(self, cm='seismic'):
        """
        This routines displays the transfer function

        :param cm: name of the colormap
        :type cm: str
        """

        fig, ax = plt.subplots()
        vmax = abs(self.trans[1:,1:]).max()
        vmin = -vmax
        im = ax.pcolormesh(self.idx2m[:self.n_m_transfer]+1,
                           self.idx2m[:self.n_m_transfer]+1, self.trans,
                           cmap=plt.get_cmap(cm),
                           vmax=vmax, vmin=vmin)
        ax.set_xlabel("m+1")
        ax.set_ylabel("m'+1")
        fig.colorbar(im)
        fig.tight_layout()


if __name__ == '__main__':

    tt = PizzaTransfer()
    plt.show()
