# -*- coding: utf-8 -*-
import os
import json
from types import SimpleNamespace
from collections import OrderedDict
import numpy as np
from .libpizza import scanDir, avgField
from .series import PizzaTs
from .log import PizzaSetup
from .spectrum import PizzaSpectrum
from .radial import PizzaRadial

default_model = os.path.join(os.environ['HOME'], 'pizza/python/pizza/model.json')

INDENT = 3
SPACE = " "
NEWLINE = "\n"
def to_json(o, level=0):
    """
    This is a manual JSON serializer function that makes the outputs look a little
    bit better than default.
    """
    ret = ""
    if isinstance(o, dict):
        ret += "{" + NEWLINE
        comma = ""
        for k, v in o.items():
            ret += comma
            comma = ",\n"
            ret += SPACE * INDENT * (level + 1)
            ret += '"' + str(k) + '":' + SPACE
            ret += to_json(v, level + 1)

        ret += NEWLINE + SPACE * INDENT * level + "}"
    elif isinstance(o, str):
        ret += '"' + o + '"'
    elif isinstance(o, list):
        ret += "[" + ",".join([to_json(e, level + 1) for e in o]) + "]"
    # Tuples are interpreted as lists
    elif isinstance(o, tuple):
        ret += "[" + ",".join(to_json(e, level + 1) for e in o) + "]"
    elif isinstance(o, bool):
        ret += "true" if o else "false"
    elif isinstance(o, int):
        ret += str(o)
    elif isinstance(o, float):
        if abs(o) > 1e2 or abs(o) < 1e-2:
            ret += f'{o:.8e}'
        else:
            ret += f'{o:.8g}'
    elif isinstance(o, np.ndarray) and np.issubdtype(o.dtype, np.integer):
        ret += "[" + ','.join(map(str, o.flatten().tolist())) + "]"
    elif isinstance(o, np.ndarray) and np.issubdtype(o.dtype, np.bool_):
        ret += "[" + ','.join(map(lambda x: f'"{x}"', o.flatten().tolist())) + "]"
    elif isinstance(o, np.ndarray) and np.issubdtype(o.dtype, np.str_):
        ret += "[" + ','.join(map(lambda x: f'"{x}"', o.flatten().tolist())) + "]"
    elif isinstance(o, np.ndarray) and np.issubdtype(o.dtype, np.inexact):
        ret += "[" + ','.join(map(lambda x: f'{x:.8e}', o.flatten().tolist())) + "]"
    elif o is None:
        ret += 'null'
    else:
        raise TypeError(f"Unknown type '{type(o)!s}' for json serialization")

    return ret


class PizzaAvg:
    """
    This class computes the time-average properties from time series, spectra
    and radial profiles. It will store the input starting time in a small file
    named ``tInitAvg``, such that the next time you use it you don't need to
    provide ``tstart`` again. By default, the outputs are stored in a
    fully documented JSON file named avg.json: this is split into several
    categories, namely numerical parameters, physical parameters, time averaged
    scalar quantities, time averaged spectra and time-averaged radial profiles.
    The quantities stored in the JSON file are entirely controlled by an input
    model file wich enlists the different quantities of interest.
    A default example file named ``model.json`` is provided in
    ``$MAGIC_HOME/python/magic``, and an example of how to build a dedicated
    one is provided below.

    >>> # Average from t=2.11
    >>> a = PizzaAvg(tstart=2.11)
    >>> # Average only the files that match the pattern N0m2[a-c]
    >>> a = PizzaAvg(tstart=2.11, tag='N0m2[a-c]')
    >>> print(a) # print the formatted output
    >>> # Custom JSON model to select averages
    >>> json_model = { 'phys_params': ['ek', 'ra'],
                       'time_series': { 'heat': ['topnuss', 'botnuss'],
                                        'e_kin': ['us2', 'up2'] },
                       'spectra': ['us2_m'],
                       'radial_profiles': ['uphi_r', 'temp_r']}
                     }
    >>> # Compute the selected averages in the dirctory mydir
    >>> a = PizzaAvg(datadir='mydir', model=json_model)
    """

    def __init__(self, tag=None, tstart=None, tstartHeat=None, std=False,
                 model=default_model, datadir='.', write=True,
                 tinit_file='tInitAvg'):
        """
        :param tag: if you specify an input tag (generic regExp pattern),
                    the averaging process will only happen on the time series
                    that match this input pattern
        :type tag: str
        :param tstart: the starting time for averaging
        :type tstart: float
        :param tstartHeat: the starting time for averaging heat transfer
                           (in case this is different)
        :type tstartHeat: float
        :param datadir: working directory
        :type datadir: str
        :param std: compute the standard deviation when set to True
        :type std: bool
        :param write: write the outputs in a JSON file
        :type write: bool
        :param model: this is the path of a JSON file which defines
                      which fields will be handled in the time-averaging
                      process. This can be any python attributes wich is
                      defined in PizzaTs, PizzaSpectrum or PizzaRadial.
        :type model: str
        :param tinit_file: name of the file which contains the start time
        :type tinit_file: str
        """

        if not os.path.exists(datadir):
            print(f'Directory "{datadir}" has not been found')
            return

        tInitFile = os.path.join(datadir, tinit_file)
        if os.path.exists(tInitFile) and tstart is None:
            with open(tInitFile, 'r') as f:
                st = f.readline().strip('\n')
                tstart = float(st)
        elif tstart is not None:
            with open(tInitFile, 'w') as f:
                f.write(f'{tstart}')

        if os.path.exists('tstartHeat') and tstartHeat is None:
            with open('tstartHeat', 'r') as file:
                st = file.readline().strip('\n')
                tstartHeat = float(st)
        else:
            if tstartHeat is None:
                tstartHeat = tstart
            with open('tstartHeat', 'w') as file:
                file.write(f'{tstartHeat}')

        if type(model) is str:
            with open(model, 'r') as f:
                params = json.load(f)
        else:  # This is directly a json dict
            params = model

        pattern = os.path.join(datadir, 'log.*')
        logs = scanDir(pattern)

        self.lut = OrderedDict()
        if 'phys_params' in params:
            self.lut['phys_params'] = {}
        if 'num_params' in params:
            self.lut['num_params'] = {}

        # First grab the requested control parameters
        if len(logs) > 0:
            stp = PizzaSetup(nml=logs[-1], quiet=True)
            if 'phys_params' in self.lut:
                for p in params['phys_params']:
                    if hasattr(stp, p):
                        self.lut['phys_params'][p] = getattr(stp, p)
                        setattr(self, p, getattr(stp, p))
            if 'num_params' in self.lut:
                for p in params['num_params']:
                    if hasattr(stp, p):
                        self.lut['num_params'][p] = getattr(stp, p)
                        setattr(self, p, getattr(stp, p))

        # Handle time series
        self.lut['time_series'] = {}
        for k, key in enumerate(params['time_series'].keys()):
            ts = PizzaTs(field=key, datadir=datadir, all=True, tag=tag,
                         iplot=False)
            if hasattr(ts, 'time'):  # Manage to read file

                if key == 'heat':
                    ind = np.argmin(abs(ts.time-tstartHeat))
                else:
                    ind = np.argmin(abs(ts.time-tstart))

                if k == 0:
                    self.lut['time_series']['tavg'] = ts.time[-1]-ts.time[ind]
                    self.lut['time_series']['total_time'] = ts.time[-1]-ts.time[0]
                    self.tavg = ts.time[-1]-ts.time[ind]
                    self.trun = ts.time[-1]-ts.time[0]

                for field in params['time_series'][key]:
                    if hasattr(ts, field):
                        if std and field != 'dt':
                            xmean, xstd = avgField(ts.time[ind:],
                                                   ts.__dict__[field][ind:],
                                                   std=True)
                            self.lut['time_series'][field+'_av'] = xmean
                            self.lut['time_series'][field+'_sd'] = xstd
                            setattr(self, field+'_av', xmean)
                            setattr(self, field+'_sd', xstd)
                        else:
                            xmean = avgField(ts.time[ind:],
                                             ts.__dict__[field][ind:])
                            self.lut['time_series'][field+'_av'] = xmean
                            setattr(self, field+'_av', xmean)
            else:  # If parameters is absent then put it to -1
                for field in params['time_series'][key]:
                    self.lut['time_series'][field+'_av'] = -1
                    setattr(self, field+'_av', -1)
                    if std:
                        self.lut['time_series'][field+'_sd'] = -1
                        setattr(self, field+'_sd', -1)

        # Get tags involved in averaging for spectra and radial profiles
        tags = self.get_tags(datadir, tstart)

        # Handle spectra
        self.lut['spectra'] = {}
        # Determine whether file exists
        if len(tags) > 0:
            file_exists = True
            # If only one tag is retained but averaged file does not exist yet
            if len(tags) == 1:
                file = os.path.join(datadir, 'spec_avg.' + tags[-1])
                if not os.path.exists(file):
                    file_exists = False
        else:
            file_exists = False

        if len(params['spectra']) > 0:
            if file_exists:
                for k, tag in enumerate(tags):
                    if k == 0:
                        sp = PizzaSpectrum(datadir=datadir, iplot=False,
                                           tag=tag, quiet=True)
                    else:
                        sp += PizzaSpectrum(datadir=datadir, iplot=False,
                                            tag=tag, quiet=True)
                if hasattr(sp, 'index'):  # Manage to read file
                    self.lut['spectra']['index'] = sp.index
                    for field in params['spectra']:
                        if hasattr(sp, field):
                            self.lut['spectra'][field+'_av'] = sp.__dict__[field]
                            if std and hasattr(sp, field + '_SD'):
                                self.lut['spectra'][field+'_sd'] = \
                                    sp.__dict__[field + '_SD']
            else:  # Set parameters to -1
                self.lut['spectra']['index'] = -1 * np.ones(32)
                for field in params['spectra']:
                    self.lut['spectra'][field+'_av'] = -1 * np.ones(32)
                    if std:
                        self.lut['spectra'][field+'_sd'] = -1 * np.ones(32)

        # Handle radial profiles
        self.lut['radial_profiles'] = {}
        if len(params['spectra']) > 0:
            if file_exists:
                for k, tag in enumerate(tags):
                    if k == 0:
                        rr = PizzaRadial(datadir=datadir, iplot=False,
                                         tag=tag, quiet=True)
                    else:
                        rr += PizzaRadial(datadir=datadir, iplot=False,
                                          tag=tag, quiet=True)
                if hasattr(rr, 'radius'):  # Manage to read file
                    self.lut['radial_profiles']['radius'] = rr.radius
                    for field in params['radial_profiles']:
                        if hasattr(rr, field):
                            self.lut['radial_profiles'][field+'_av'] = \
                                rr.__dict__[field]
                            #setattr(self, field+'_av', rr.__dict__[field])
                            if std and hasattr(rr, field + '_SD'):
                                self.lut['radial_profiles'][field+'_sd'] = \
                                    rr.__dict__[field + '_SD']
                                #setattr(self, field+'_sd', rr.__dict__[field+'_SD'])
            else:  # Set parameters to -1
                self.lut['radial_profiles']['radius'] = -1 * np.ones(33)
                for field in params['radial_profiles']:
                    self.lut['radial_profiles'][field+'_av'] = -1 * np.ones(33)
                    #setattr(self, field+'R_av', rr.__dict__[field])
                    if std:
                        self.lut['radial_profiles'][field+'_sd'] = -1 * np.ones(33)
                        #setattr(self, field+'R_sd', rr.__dict__[field+'_SD'])

        # Write a json file
        if write:
            self.write_json(datadir)

    def write_header(self):
        """
        Write header in case an ascii output is requested
        """
        st = '#'

        if 'phys_params' in self.lut:
            for key in self.lut['phys_params']:
                st += key.rjust(max(10, len(key)+1))
        for key in self.lut['time_series']:
            st += key.rjust(max(16, len(key)+1))
        if 'num_params' in self.lut:
            for key in self.lut['num_params']:
                par = self.lut['num_params'][key]
                if type(par) is not str and type(par) is not bool:
                    st += key.rjust(len(key)+1)

        return st

    def __str__(self):
        """
        Formatted output
        """
        st = ' '

        if 'phys_params' in self.lut:
            for par in self.lut['phys_params'].values():
                st += f'{par:10.3e}'
        for par in self.lut['time_series'].values():
            st += f'{par:16.8e}'
        if 'num_params' in self.lut:
            for par in self.lut['num_params'].values():
                if type(par) is not str and type(par) is not bool:
                    st += f' {par:g}'

        return st

    def get_tags(self, datadir, tstart):
        """
        This routine returns a list of tags which have been produced
        for t>tstart.

        :param datadir: working directory
        :type datadir: str
        :param tstart: starting averaging time
        :type tstart: float
        :returns: a list of tags
        :rtype: list
        """
        logFiles = scanDir(os.path.join(datadir, 'log.*'))
        tags = []
        for lg in logFiles:
            nml = PizzaSetup(nml=lg, quiet=True)
            if hasattr(nml, 'start_time'):
                if nml.start_time > tstart:
                    tags.append(nml.tag)

        return tags

    def write_json(self, datadir='.'):
        """
        This function writes the averages as a simple JSON file stored in the
        directory 'avg.json'

        :param datadir: working directory
        :type datadir: str
        """
        with open(os.path.join(datadir, 'avg.json'), 'w') as f:
            st = to_json(self.lut)
            f.write(st)


class PizzaAvgStack:
    """
    This class is made to go through a list of directories and gather all the
    local avg files to compute a global output which summarises all the
    outputs in one single file

    >>> # Simply go through the directories listed in "runs.txt" and produce a
    >>> # local file named "my_results.json"
    >>> st = PizzaAvgStack(dirList="runs.txt", filename="my_results.json")
    >>> # Now also recompute each individual avg.json file in each directory
    >>> st = PizzaAvgStack(dirList="runs.txt", filename="my_results.json",
                           recompute=True, module="path_to_model.json",
                           std=True)
    """

    def __init__(self, dirList='runs.txt', filename='avg.json', datadir='.',
                 recompute=False, model=default_model, std=False,
                 readonly=False):
        """
        :param dirList: the filename of the list of directories, lines starting
                        with a hash are omitted by the reader
        :type dirList: str
        :param filename:
        :param recompute: recompute each individual average file in each single
                          directory when set to True
        :type recompute: bool
        :param datadir: working directory
        :type datadir: str
        :param std: compute the standard deviation when set to True
        :type std: bool
        :param model: this is the path of a JSON file which defines which
                      fields will be handled in the time-averaging process.
                      This can be any python attributes wich is defined in
                      PizzaTs, PizzaSpectrum or PizzaRadial.
        :type model: str
        :param readonly: when set to True, simply read the summary json file
        :param readonly: bool
        """

        if readonly:  # Only read the json file
            with open(os.path.join(datadir, filename), 'r') as f:
                self.lut = json.load(f)
        else:
            with open(os.path.join(datadir, dirList), 'r') as f:
                dirs = f.readlines()

            startdir = os.getcwd()
            for k, dir in enumerate(dirs):
                if not dir.startswith('#'):
                    print(f"In dir {dir.rstrip('\n')}")
                    os.chdir(os.path.join(datadir, dir.rstrip('\n')))
                    if recompute:
                        PizzaAvg(model=model, std=std)
                    if not hasattr(self, 'lut'):
                        self.lut = self.load()
                        keys = ['spectra', 'radial_profiles']
                        for key in keys:
                            if key in self.lut.keys():
                                for key1 in self.lut[key].keys():
                                    self.lut[key][key1] = \
                                        [np.atleast_1d(self.lut[key][key1])]
                    else:
                        lut = self.load()
                        self.merge(lut)

                    os.chdir(startdir)

            with open(os.path.join(datadir, filename), 'w') as f:
                st = to_json(self.lut)
                f.write(st)

        self.simple_namespace()

    def __add__(self, new):
        """
        This routine is used to add two PizzaAvgStack objects.

        :param new: the lookup table which needs to be added
        :type new: PizzaAvgStack
        """
        keys = ['phys_params', 'num_params', 'time_series']
        for key in keys:
            if key in new.lut.keys() and key in self.lut.keys():
                for key1 in self.lut[key].keys():
                    if key1 in new.lut[key].keys():
                        arr = np.atleast_1d(self.lut[key][key1])
                        arr1 = np.atleast_1d(new.lut[key][key1])
                        self.lut[key][key1] = np.concatenate((arr, arr1))

        keys = ['spectra', 'radial_profiles']
        for key in keys:
            if key in new.lut.keys() and key in self.lut.keys():
                for key1 in self.lut[key].keys():
                    if key1 in new.lut[key].keys():
                        for lst in new.lut[key][key1]:
                            arr1 = np.atleast_1d(lst)
                            self.lut[key][key1].append(arr1)

        self.simple_namespace()

        return self

    def simple_namespace(self):
        """
        This routine creates a simpler namespace from the lookup table
        """
        keys = ['phys_params', 'num_params', 'time_series']
        for key in keys:
            if key in self.lut.keys():
                for key1 in self.lut[key].keys():
                    setattr(self, key1, np.atleast_1d(self.lut[key][key1]))

        keys = ['spectra', 'radial_profiles']
        for key in keys:
            if key in self.lut.keys():
                for key1 in self.lut[key].keys():
                    g = []
                    for ff in self.lut[key][key1]:
                        g.append(np.atleast_1d(ff))
                    setattr(self, key1, g)

    def merge(self, newlut):
        """
        This routine is used to merge two lookup tables.

        :param newlut: the lookup table which needs to be added
        :type newlut: dict
        """
        keys = ['phys_params', 'num_params', 'time_series']
        for key in keys:
            if key in newlut.keys() and key in self.lut.keys():
                for key1 in self.lut[key].keys():
                    if key1 in newlut[key].keys():
                        arr = np.atleast_1d(self.lut[key][key1])
                        arr1 = np.atleast_1d(newlut[key][key1])
                        self.lut[key][key1] = np.concatenate((arr, arr1))

        keys = ['spectra', 'radial_profiles']
        for key in keys:
            if key in newlut.keys() and key in self.lut.keys():
                for key1 in self.lut[key].keys():
                    if key1 in newlut[key].keys():
                        arr1 = np.atleast_1d(newlut[key][key1])
                        self.lut[key][key1].append(arr1)
                    else: # Fill with -1 dummy arrays to make sure everything has the right size
                        arr1 = -1 * np.ones(32)
                        self.lut[key][key1].append(arr1)

    def load(self):
        """
        This routine reads the avg.json file stored in one local directory
        which contains a numerical simulation. It returns the corresponding
        lookup table.

        :returns: a lookup table
        :rtype: dict
        """
        with open('avg.json', 'r') as f:
            lut = json.load(f)

        return lut


if __name__ == '__main__':
    a = PizzaAvg(std=True)
    st = a.write_header()
    print(st)
    print(a)
