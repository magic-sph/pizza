# -*- coding: utf-8 -*-
import os
import numpy as np
from .series import PizzaTs
from .libpizza import avgField


class PizzaAverages:
    """
    This class calculates the time-average properties from time series. It will
    store the input starting time in a small file named ``tInitAvg``, such that
    the next time you use it you don't need to give ``tstart`` again.
    """

    def __init__(self, tstart=None, tag=None, std=False, tstartHeat=None,
                 l_3D=True):
        """
        :param tstart: the starting time for averaging
        :type tstart: float
        :param tag: if you specify an input tag (generic regExp pattern),
                    the averaging process will only happen on the time series
                    that match this input pattern
        :type tag: str
        :type std: compute the standard deviation when set to True
        :type std: bool
        :param l_3D: when turned on, 3D files are used for diagnostic
        :type l_3D: bool
        """
        self.l_3D = l_3D

        if os.path.exists('tInitAvg') and tstart is None:
            with open('tInitAvg', 'r') as file:
                st = file.readline().strip('\n')
                tstart = float(st)
        elif tstart is not None:
            with open('tInitAvg', 'w') as file:
                file.write(f'{tstart}')

        if os.path.exists('tstartHeat') and tstartHeat is None:
            with open('tstartHeat', 'r') as file:
                st = file.readline().strip('\n')
                tstartHeat = float(st)
        else:
            if tstartHeat is None:
                tstartHeat = tstart
            with open('tstartHeat', 'w') as file:
                file.write(f'{tstartHeat}')

        print(tstart, tstartHeat)

        self.std = std

        # e_kin_3D file
        if self.l_3D:
            file = 'e_kin_3D'
        else:
            file = 'e_kin'
        ts = PizzaTs(field=file, all=True, tag=tag, iplot=False)
        mask = np.where(abs(ts.time-tstart) == min(abs(ts.time-tstart)), 1, 0)
        ind = np.nonzero(mask)[0][0]

        if self.std:
            self.us2_avg, self.us2_std = avgField(ts.time[ind:], ts.us2[ind:],
                                                  std=True)
            self.up2_avg, self.up2_std = avgField(ts.time[ind:], ts.up2[ind:],
                                                  std=True)
            if self.l_3D:
                self.uz2_avg, self.uz2_std = avgField(ts.time[ind:],
                                                      ts.uz2[ind:],
                                                      std=True)
            self.up2_axi_avg, self.up2_axi_std = avgField(ts.time[ind:],
                                                          ts.up2_axi[ind:],
                                                          std=True)
        else:
            self.us2_avg = avgField(ts.time[ind:], ts.us2[ind:])
            self.up2_avg = avgField(ts.time[ind:], ts.up2[ind:])
            if self.l_3D:
                self.uz2_avg = avgField(ts.time[ind:], ts.uz2[ind:])
            self.up2_axi_avg = avgField(ts.time[ind:], ts.up2_axi[ind:])

        self.tavg = ts.time[-1]-ts.time[ind]  # Averaging time

        self.ra = ts.ra
        self.pr = ts.pr
        self.ek = ts.ek
        self.radratio = ts.radratio

        # reynolds_3D file
        if self.l_3D:
            file = 'reynolds_3D'
        else:
            file = 'reynolds'
        ts2 = PizzaTs(field=file, all=True, iplot=False, tag=tag)
        mask = np.where(abs(ts2.time-tstart) == min(abs(ts2.time-tstart)),
                        1, 0)
        ind = np.nonzero(mask)[0][0]

        if self.std:
            self.rey_avg, self.rey_std = avgField(ts2.time[ind:],
                                                  ts2.rey[ind:], std=True)
            self.rey_zon_avg, self.rey_zon_std = avgField(ts2.time[ind:],
                                                          ts2.rey_zon[ind:],
                                                          std=True)
            self.rey_fluct_avg, self.rey_fluct_std \
                = avgField(ts2.time[ind:], ts2.rey_fluct[ind:],
                           std=True)
        else:
            self.rey_avg = avgField(ts2.time[ind:], ts2.rey[ind:])
            self.rey_zon_avg = avgField(ts2.time[ind:], ts2.rey_zon[ind:])
            self.rey_fluct_avg = avgField(ts2.time[ind:], ts2.rey_fluct[ind:])

        # heat.TAG files
        ts3 = PizzaTs(field='heat', all=True, iplot=False, tag=tag)
        mask = np.where(abs(ts3.time-tstartHeat) == min(abs(ts3.time-tstartHeat)), 1, 0)
        ind = np.nonzero(mask)[0][0]

        if self.std:
            self.nu_bot_avg, self.nu_bot_std = \
                avgField(ts3.time[ind:], ts3.botnuss[ind:], std=True)
            self.nu_top_avg, self.nu_top_std = \
                avgField(ts3.time[ind:], ts3.topnuss[ind:], std=True)
            self.nu_vol_avg, self.nu_vol_std = \
                avgField(ts3.time[ind:], ts3.volnuss[ind:], std=True)
            self.nu_shell_avg, self.nu_shell_std = \
                avgField(ts3.time[ind:], ts3.shellnuss[ind:], std=True)

        else:
            self.nu_bot_avg = avgField(ts3.time[ind:], ts3.botnuss[ind:])
            self.nu_top_avg = avgField(ts3.time[ind:], ts3.topnuss[ind:])
            self.nu_vol_avg = avgField(ts3.time[ind:], ts3.volnuss[ind:])
            self.nu_shell_avg = avgField(ts3.time[ind:], ts3.shellnuss[ind:])

        # power3D.TAG files
        if self.l_3D:
            file = 'power_3D'
        else:
            file = 'power'
        ts4 = PizzaTs(field=file, all=True, iplot=False, tag=tag)
        mask = np.where(abs(ts4.time-tstart) == min(abs(ts4.time-tstart)),
                        1, 0)
        ind = np.nonzero(mask)[0][0]

        if self.std:
            self.power_avg, self.power_std = \
                avgField(ts4.time[ind:], ts4.buoPower[ind:], std=True)
            self.visc_avg, self.visc_std = \
                avgField(ts4.time[ind:], ts4.viscDiss[ind:], std=True)
        else:
            self.power_avg = avgField(ts4.time[ind:], ts4.buoPower[ind:])
            self.visc_avg = avgField(ts4.time[ind:], ts4.viscDiss[ind:])

        # length_scales.TAG files
        ts5 = PizzaTs(field='length_scales', all=True, iplot=False, tag=tag)
        mask = np.where(abs(ts5.time-tstart) == min(abs(ts5.time-tstart)),
                        1, 0)
        ind = np.nonzero(mask)[0][0]

        if self.std:
            self.lpeak_avg, self.lpeak_std = \
                avgField(ts5.time[ind:], ts5.lus_peak[ind:], std=True)
            self.lint_avg, self.lint_std = \
                avgField(ts5.time[ind:], ts5.lint[ind:], std=True)
            self.ldiss_avg, self.ldiss_std = \
                avgField(ts5.time[ind:], ts5.ldiss[ind:], std=True)
        else:
            self.lpeak_avg = avgField(ts5.time[ind:], ts5.lus_peak[ind:])
            self.lint_avg = avgField(ts5.time[ind:], ts5.lint[ind:])
            self.ldiss_avg = avgField(ts5.time[ind:], ts5.ldiss[ind:])

    def __str__(self):
        """
        Formatted output
        """
        st_std = ''
        if self.l_3D:
            st = f'{self.ra:.3e}{self.ek:9.2e}{self.pr:9.2e}{self.radratio:9.2e}'
            st += f'{self.us2_avg:12.5e}{self.up2_avg:12.5e}{self.uz2_avg:12.5e}'
            st += f'{self.up2_axi_avg:12.5e}'
        else:
            st = f'{self.ra:.3e}{self.pr:9.2e}{self.radratio:9.2e}'
            st += f'{self.us2_avg:12.5e}{self.up2_avg:12.5e}{self.up2_axi_avg:12.5e}'
        if self.std:
            st_std = f'{self.us2_std:12.5e}{self.up2_std:12.5e}'
            st_std += f'{self.uz2_std:12.5e}{self.up2_axi_avg:12.5e}'

        st += f'{self.rey_avg:10.3e}{self.rey_zon_avg:10.3e}'
        st += f'{self.rey_fluct_avg:10.3e}'

        if self.std:
            st_std += f'{self.rey_std:10.3e}{self.rey_zon_std:10.3e}'
            st_std += f'{self.rey_fluct_std:10.3e}'
        st += f'{self.nu_top_avg:10.3e}{self.nu_bot_avg:10.3e}'
        st += f'{self.nu_vol_avg:10.3e}{self.nu_shell_avg:10.3e}'
        if self.std:
            st_std += f'{self.nu_top_std:10.3e}{self.nu_bot_std:10.3e}'
            st_std += f'{self.nu_vol_std:10.3e}{self.nu_shell_std:10.3e}'

        st += f'{self.power_avg:12.5e}{self.visc_avg:12.5e}'
        if self.std:
            st_std += f'{self.power_std:12.5e}{self.visc_std:12.5e}'

        st += f'{self.lpeak_avg:10.3e}{self.lint_avg:10.3e}{self.ldiss_avg:10.3e}'
        if self.std:
            st_std += f'{self.lpeak_std:10.3e}{self.lint_std:10.3e}'
            st_std += f'{self.ldiss_std:10.3e}'
        st += st_std
        st += '\n'

        return st
