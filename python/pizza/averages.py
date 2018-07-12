#-*- coding: utf-8 -*-
import os, re
import numpy as np
from .series import PizzaTs
import glob
from .libpizza import avgField

class PizzaAverages:
    """
    This class calculates the time-average properties from time series. It will
    store the input starting time in a small file named ``tInitAvg``, such that
    the next time you use it you don't need to give ``tstart`` again.
    """

    def __init__(self, tstart=None, tag=None, std=False):
        """
        :param tstart: the starting time for averaging
        :type tstart: float
        :param tag: if you specify an input tag (generic regExp pattern), 
                    the averaging process will only happen on the time series 
                    that match this input pattern
        :type tag: str
        :type std: compute the standard deviation when set to True
        :type std: bool
        """

        if os.path.exists('tInitAvg') and tstart is None:
            file = open('tInitAvg', 'r')
            st = file.readline().strip('\n')
            tstart = float(st)
            file.close()
        elif tstart is not None:
            file = open('tInitAvg', 'w')
            file.write('%f' % tstart)
            file.close()

        self.std = std

        # e_kin_3D file
        ts = PizzaTs(field='e_kin_3D', all=True, tag=tag, iplot=False)
        mask = np.where(abs(ts.time-tstart) == min(abs(ts.time-tstart)), 1, 0)
        ind = np.nonzero(mask)[0][0]

        if self.std:
            self.us2_avg, self.us2_std = avgField(ts.time[ind:], ts.us2[ind:],
                                                  std=True)
            self.up2_avg, self.up2_std = avgField(ts.time[ind:], ts.up2[ind:],
                                                  std=True)
            self.uz2_avg, self.uz2_std = avgField(ts.time[ind:], ts.uz2[ind:],
                                                  std=True)
            self.up2_axi_avg, self.up2_axi_std = avgField(ts.time[ind:], 
                                                 ts.up2_axi[ind:], std=True)
        else:
            self.us2_avg = avgField(ts.time[ind:], ts.us2[ind:])
            self.up2_avg = avgField(ts.time[ind:], ts.up2[ind:])
            self.uz2_avg = avgField(ts.time[ind:], ts.uz2[ind:])
            self.up2_axi_avg = avgField(ts.time[ind:], ts.up2_axi[ind:])

        self.tavg = ts.time[-1]-ts.time[ind] # Averaging time

        self.ra = ts.ra
        self.pr = ts.pr
        self.ek = ts.ek
        self.radratio = ts.radratio

        # reynolds_3D file
        ts2 = PizzaTs(field='reynolds_3D', all=True, iplot=False, tag=tag)
        mask = np.where(abs(ts2.time-tstart) == min(abs(ts2.time-tstart)), 1, 0)
        ind = np.nonzero(mask)[0][0]

        if self.std:
            self.rey_avg, self.rey_std = avgField(ts2.time[ind:], ts2.rey[ind:],
                                                  std=True)
            self.rey_zon_avg, self.rey_zon_std = avgField(ts2.time[ind:],
                                                ts2.rey_zon[ind:], std=True)
            self.rey_fluct_avg, self.rey_fluct_std = avgField(ts2.time[ind:],
                                                ts2.rey_fluct[ind:], std=True)
        else:
            self.rey_avg = avgField(ts2.time[ind:], ts2.rey[ind:])
            self.rey_zon_avg = avgField(ts2.time[ind:], ts2.rey_zon[ind:])
            self.rey_fluct_avg = avgField(ts2.time[ind:], ts2.rey_fluct[ind:])

        # heat.TAG files
        ts3 = PizzaTs(field='heat', all=True, iplot=False, tag=tag)
        mask = np.where(abs(ts3.time-tstart) == min(abs(ts3.time-tstart)), 1, 0)
        ind = np.nonzero(mask)[0][0]

        if self.std:
            self.nu_bot_avg, self.nu_bot_std = avgField(ts3.time[ind:], 
                                               ts3.botnuss[ind:], std=True)
            self.nu_top_avg, self.nu_top_std = avgField(ts3.time[ind:], 
                                               ts3.topnuss[ind:], std=True)
            self.nu_vol_avg, self.nu_vol_std = avgField(ts3.time[ind:], 
                                               ts3.volnuss[ind:], std=True)
            self.nu_shell_avg, self.nu_shell_std = avgField(ts3.time[ind:], 
                                               ts3.shellnuss[ind:], std=True)

        else:
            self.nu_bot_avg = avgField(ts3.time[ind:], ts3.botnuss[ind:])
            self.nu_top_avg = avgField(ts3.time[ind:], ts3.topnuss[ind:])
            self.nu_vol_avg = avgField(ts3.time[ind:], ts3.volnuss[ind:])
            self.nu_shell_avg = avgField(ts3.time[ind:], ts3.shellnuss[ind:])

        # power3D.TAG files
        ts4 = PizzaTs(field='power_3D', all=True, iplot=False, tag=tag)
        mask = np.where(abs(ts4.time-tstart) == min(abs(ts4.time-tstart)), 1, 0)
        ind = np.nonzero(mask)[0][0]

        if self.std:
            self.power_avg, self.power_std = avgField(ts4.time[ind:], 
                                             ts4.buoPower[ind:], std=True)
            self.visc_avg, self.visc_std = avgField(ts4.time[ind:], 
                                           ts4.viscDiss[ind:], std=True)
        else:
            self.power_avg = avgField(ts4.time[ind:], ts4.buoPower[ind:])
            self.visc_avg = avgField(ts4.time[ind:], ts4.viscDiss[ind:])

        # length_scales.TAG files
        ts5 = PizzaTs(field='length_scales', all=True, iplot=False, tag=tag)
        mask = np.where(abs(ts5.time-tstart) == min(abs(ts5.time-tstart)), 1, 0)
        ind = np.nonzero(mask)[0][0]

        if self.std:
            self.lpeak_avg, self.lpeak_std = avgField(ts5.time[ind:], 
                                                 ts5.lus_peak[ind:], std=True)
            self.lint_avg, self.lint_std = avgField(ts5.time[ind:], 
                                                    ts5.lint[ind:], std=True)
            self.ldiss_avg, self.ldiss_std = avgField(ts5.time[ind:], 
                                                      ts5.ldiss[ind:], std=True)
        else:
            self.lpeak_avg = avgField(ts5.time[ind:], ts5.lus_peak[ind:])
            self.lint_avg = avgField(ts5.time[ind:], ts5.lint[ind:])
            self.ldiss_avg = avgField(ts5.time[ind:], ts5.ldiss[ind:])


    def __str__(self):
        """
        Formatted output
        """
        st_std=''
        st = '%.3e%9.2e%9.2e%9.2e%12.5e%12.5e%12.5e%12.5e' % \
              (self.ra, self.ek, self.pr, self.radratio, self.us2_avg, \
               self.up2_avg, self.uz2_avg, self.up2_axi_avg)
        if self.std:
            st_std = '%12.5e%12.5e%12.5e%12.5e' % \
                     (self.us2_std, self.up2_std, self.uz2_std, self.up2_axi_avg)
             
        st +='%10.3e%10.3e%10.3e' % (self.rey_avg, self.rey_zon_avg, self.rey_fluct_avg)

        if self.std:
            st_std +='%10.3e%10.3e%10.3e' % (self.rey_std, self.rey_zon_std,
                                          self.rey_fluct_std)
        st += '%10.3e%10.3e%10.3e%10.3e' % (self.nu_top_avg, self.nu_bot_avg,\
                                            self.nu_vol_avg, self.nu_shell_avg)
        if self.std:
            st_std += '%10.3e%10.3e%10.3e%10.3e' % (self.nu_top_std, self.nu_bot_std,
                                                    self.nu_vol_std, self.nu_shell_std)

        st += '%12.5e%12.5e' % (self.power_avg, self.visc_avg)
        if self.std:
            st_std += '%12.5e%12.5e' % (self.power_std, self.visc_std)

        st += '%10.3e%10.3e%10.3e' % (self.lpeak_avg, self.lint_avg, self.ldiss_avg)

        if self.std:
            st_std += '%10.3e%10.3e%10.3e' % (self.lpeak_std, self.lint_std,
                                           self.ldiss_std)
        st += st_std
        st += '\n'

        return st
