#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
import astropy.time as asttime
import astropy.coordinates as astcoords
import astropy.units as astunits
from astropy.utils.exceptions import AstropyUserWarning
#import astropy.constants as astconst
import matplotlib
matplotlib.use('agg')    # Otherwise keeps hanging when using CERES; GUIs work, but plt.show() won't; matplotlib.interactive(True) doesn't show
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import time
import datetime
import operator
import copy
import random
import warnings
from scipy.optimize import curve_fit            # Curve-fit Gauss takes the longest: 20ms
# import scipy.interpolate as inter             # Not used anymore
import scipy
import scipy.signal
try:
    from tqdm import tqdm
except:
    print('Error: tqdm could not be loaded. Did you activate the Anaconda hiflex environment?')
    exit(1)
import tkcanvas as tkc
import json
# detect python version
if sys.version_info[0] < 3:     # Python 2
    import Tkinter as Tk
    from collections import OrderedDict as dict
    import urllib2
else:                           # Python 3
    import tkinter as Tk
    import urllib
    raw_input = input
if sys.version_info[0] == 3 and sys.version_info[1] >= 5:
    try:
        from deepCR import deepCR   # Only available for python 3.5
    except:
        print('Error: deepCR could not be loaded. Did you activate the Anaconda hiflex environment?')
        deepCR = None
else:
    deepCR = None
import plot_img_spec
import psutil
import ephem
import math
import multiprocessing
import subprocess
# Necessary because of https://github.com/astropy/astropy/issues/9427       # commented out again on 20/5/2020 after "WARNING: IERSStaleWarning: IERS_Auto predictive values are older than 15.0 days but downloading the latest table did not find newer values [astropy.utils.iers.iers]"
try:
    import astropy
except:
    print('Error: astropy could not be loaded. Did you activate the Anaconda hiflex environment?')
    exit(1)
from astropy.utils.iers import conf as iers_conf 
#iers_conf.iers_auto_url = 'https://astroconda.org/aux/astropy_mirror/iers_a_1/finals2000A.all'         # is too old as of May 2020
iers_conf.iers_auto_url = 'https://datacenter.iers.org/data/9/finals2000A.all'                          # untested
iers_conf.iers_auto_url = 'ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all'                 # worked on 25 May 2020
iers_conf.auto_max_age = None 
success = False
for ii in range(5):
    if sys.version_info[0] < 3:     # Python 2
        try:
            import barycorrpy
            success = True
            break
        except (urllib2.URLError, ValueError, astropy.utils.iers.iers.IERSRangeError) as e:
            print('Warn: Cannot import barrycorrpy. Will try {0} more times. Error: {1}, Reason: {2}'.format(4-ii, e, e.reason))
            try:
                print('Warn: Cannot import barrycorrpy. Will try {0} more times. Error: {1}, Reason: {2}'.format(4-ii, e, e.reason))
            except:
                print('Warn: Cannot import barrycorrpy. Will try {0} more times. Error: {1}'.format(4-ii, e))
    else:
        try:
            import barycorrpy
            success = True
            break
        except (urllib.error.URLError, ValueError, astropy.utils.iers.iers.IERSRangeError) as e:
            print('Warn: Cannot import barrycorrpy. Will try {0} more times. Error: {1}, Reason: {2}'.format(4-ii, e, e.reason))
if not success:
    print('Error: barrycorrpy could not be loaded. It needs an active internet connection in order to download the IERS_B file. This failure will lead to a crash of the program later!'+os.linesep)
import glob
import pickle

""" only needed for BJD calculation using jplephem and BVC calculation from CERES pipeline
import jplephem                     # jplehem from the pip install jplehem
import de423
# SSEphem package http://www.cv.nrao.edu/~rfisher/Python/py_solar_system.html coupled to SOFA http://www.iausofa.org/
# to clean up: mv SSEphem/ ssephem_update.py old/
baryc_dir = os.path.dirname(os.path.abspath(__file__))+'/SSEphem/'
sys.path.append(baryc_dir)
if not os.path.isfile(baryc_dir+'man_jplephem.py'):
    os.system('cd {0} && ln jplephem.py man_jplephem.py'.format(baryc_dir))
    time.sleep(1)
    if not os.path.isfile(baryc_dir+'man_jplephem.py'):
        print('Please run in a different termina the following line and then press Enter')
        print('cd {0} && ln jplephem.py man_jplephem.py'.format(baryc_dir))
        raw_input('')
import man_jplephem                 # jplephem in the SSEphem folder
ephemeris = 'DEc403'        # To use the JPL DE403 ephemerides, https://en.wikipedia.org/wiki/Jet_Propulsion_Laboratory_Development_Ephemeris
"""
tqdm.monitor_interval = 0   #On the virtual machine at NARIT the code raises an exception otherwise

calimages = dict()  # dictionary for all calibration images used by create_image_general and read_file_calibration
# oldorder = 12     used for correlate_UI
beg_ts = time.time()
GLOBALutils, correlation, lowess = None, None, None

class Constants:                    # from CERES
    "Here I declare the constants I will use in the different functions"
    "G,c: the gravity and speed of light constants in SI units; mearth and mmoon: the earth and moon masses in kg"
    G,c = 6.673E-11,2.99792458E8
    mearth, mmoon = 5.9736E24,7.349E22 
    "the mass of the planets,moon and sun in earth masses, the ennumeration is sun, moon, mercury,venus,mars,jupiter,saturn, uranus,neptune,pluto"    
    mplanets = np.array([332981.787,0.0123000371,0.05528,0.815,0.10745,317.83,95.19,14.536,17.147,0.0021])
    "conversion from degrees to radians, and from hours to degrees, and from degrees to hours"
    degtorad = np.pi/180.0
    radtodeg = 180.0/np.pi
    HtoDeg = 360.0/24.0
    DegtoH = 1.0/HtoDeg
    "Req: equatorial radius    f: polar flatenning , w angular speedof earth in rad/s"
    Req = 6378136.6 #in m
    f = 1.0/298.256420
    w = 7.2921158554E-5 
    DaystoYear = 1.0/365.256363004
    MJD0 = 2400000.5

def get_timing(text=''):
    global beg_ts
    end_ts = time.time()
    print(text+"elapsed time: %f" % (end_ts - beg_ts))
    beg_ts = time.time()
    
def time_usage(func):
    def wrapper(*args, **kwargs):
        beg_ts = time.time()
        retval = func(*args, **kwargs)
        end_ts = time.time()
        print("elapsed time: %f" % (end_ts - beg_ts))
        return retval
    return wrapper
    
def logger(message, show=True, printarrayformat=[], printarray=[], logfile='logfile'):
    """
    Saves the status information to a logfile
    :param message: Text to log
    :param show: if False than it will be logged only to the logfile but not printed on the screen
    :param printarrayformat: Format of the columns in printarray, e.g. '%3.1f', '%2.2i'
    :param printarray: Array with values to log
    :param logfile: filename in which the information will be written
    """
    if show:
        print(message)
    with open(logfile, 'a') as file:
        if multiprocessing.current_process().name == 'MainProcess': pid = os.getpid()
        else:                                                       pid = '{0}-{1}'.format(os.getppid(), os.getpid())
        file.write('{1} - {2} - {3}{0}'.format( os.linesep, time.strftime("%Y%m%d%H%M%S", time.localtime()), pid, message ))
        if printarrayformat != [] and printarray != []:
            for line in printarray:
                text = ''
                for i,printformat in enumerate(printarrayformat):
                    #print printformat,line[i]
                    text += printformat%line[i] + '\t'
                file.write(text[:-1]+os.linesep)
                if show:
                    print(text)
    if message.find('Error') == 0:
        print('\t-> exiting')
        if 'params' in locals() or 'params' in globals(): 
            log_params(params) 
        exit(1)

def log_params(params):
    """
    formats the dictionary to be saved in the logfile
    """
    # Finding the python files
    text = ''
    list_of_files = glob.glob('{0}/*.py'.format(os.path.realpath(__file__).rsplit(os.sep,1)[0]))
    list_of_files = sorted(list_of_files, key=os.path.getmtime, reverse=True)
    for fname in list_of_files:
        filedata1 = read_text_file(fname, no_empty_lines=False)
        filedata2 = read_text_file(fname, no_empty_lines=True)
        text += os.linesep+'    {1}  {2}  {3} {4}  {0}'.format( fname, 
                        datetime.datetime.utcfromtimestamp(os.path.getmtime(fname)).strftime('%Y-%m-%dT%H:%M:%S'),
                        '%10.1i'%os.stat(fname).st_size, '%6.1i'%len(filedata1), '%6.1i'%len(filedata2) )
    logger('python files with unix timestamp of modification, size in bytes, number of lines, and number of non-empty lines: '+text, show=False, logfile='logfile_params')
    logger('Info: Using procedures file {0}'.format( os.path.realpath(__file__) ), show=False)
    # Log the versions of the packages
    # maybe later
    # Logging the parameters
    paramstxt = dict()
    for key in params.keys():
        if type(params[key]) in [list, np.ndarray]:
            paramstxt[key] = str(params[key])
        else:
            paramstxt[key] = params[key]
    logger('params: '+json.dumps(paramstxt, sort_keys = False, indent = 4), show=False, logfile='logfile_params')
    # Old way as real json, but not very compact
    # logger('params: '+json.dumps(params, sort_keys = False, indent = 4), show=False, logfile='logfile_params')

def read_parameterfile(textfile):
    # load text file (remove all white spaces)
    if not os.path.exists(textfile):
        logger('Error: The parameterfile {0} does not exist.'.format(textfile))
    """ # Extra steps to check that the user didn't make mistakes in the file:
    data = read_text_file(textfile, no_empty_lines=True)
    data = convert_readfile(data, [str,str], delimiter='=', replaces=[' ', '\t',['\\',' ']], ignorelines=['#'])         # this replaces too many spaces
    for line in data:
        if len(line) != 2:
            logger(('Error: The line in file {0} containing the entries {1} has the wrong format. Expected was "parameter = value(s)" .'+\
                    'Please check if the "=" sign is there or if a comment "#" is missing.').format(textfile, line))"""
    try:
        keys, values = np.genfromtxt(textfile, dtype=str, comments='#', delimiter='=', filling_values='', autostrip=True, unpack=True)
    except ValueError as error:
        print(error)
        logger(('Error: One line (see previous output (empty lines are missing in the line number counting)) in file {0} has the wrong format. Expected was "parameter = value(s)" .'+\
                'Please check if the "=" sign is there or if a comment "#" is missing.').format(textfile))
        # raise                 # just this to show the error
        # raise ValueError      # Don't do this, you'll lose the stack trace!
    if len(keys) == 0 or len(values) == 0:
        logger('Error: No values found when reading the parameterfile {0}.'.format(textfile))
    """# This bit was necesary before adding autostrip and unpack
    keys, values = data[:, 0], data[:, 1]
    for k in range(len(keys)):
        keys[k] = keys[k].replace(' ', '')
        values[k] = values[k].replace(' ', '')"""
    textparams = dict(zip(keys, values))
    
    return textparams

def read_cmdparams():
    
    cmdargs = sys.argv[1:]
    cmdparams = dict()
    for arg in cmdargs:
        if '=' in arg:
            key, value = arg.replace(' ', '').split('=')
            cmdparams[key] = value
            # From Neil, helpful to make the below list of formating shorter
            #VARIABLES = dict(filepath=str, plotpath=str, savepath=str, files=list,
            #     order_direction=str, xblur=float, yblur=float,
            #     fitmin=int, fitmax=int, width=int, minpixelsinorder=int)
            #if key in VARIABLES:
            #    # try to convert to required type (using VARIABLES dict)
            #    try:
            #        if VARIABLES[key] in [list, np.ndarray]:
            #            ckwargs[key] = VARIABLES[key](value.split(','))
            #        else:
            #            ckwargs[key] = VARIABLES[key](value)
            #    except ValueError:
            #        emsg = [key, str(VARIABLES[key])]
            #        print('Command line input not understood for' +
            #              'argument {0} must be {1}'.format(*emsg))
        elif arg in ['nocheck', 'prepare'] or 0 <= arg.find('nocheck') <= 2 or 0 <= arg.find('nogui') <= 2:
            continue        # Do nothing, just prevent the warning below
        else:
            logger('Warn: I dont know how to handle command line argument: {0}'.format(arg))
    
    return cmdparams

def textfileargs(params, textfile=None):
    """
    Imports parameters from text file located at textfile, parameters given in the command line and the parameters given with the programm code
    Parameters in the programm code will be overwritten by parameters in the text file and parameters in the text file will be overwritten by parameters from the command line
    :param params: dictionary, default parameter dictionary (the default values)
    :param textfile: string, location of config file
    :return params: dictionary with added parameters and updated default parameters 
    
    Please note, the parameters are only lists or single entries. No numpy arrays allowed
    """
    
    # Set up standard parameters for 
    params['in_shift'] = -0
    params['started_from_p3'] = 'False'
    
    emsg = 'Error in config file: '
    
    textparams = read_parameterfile(textfile)
    params.update(textparams)
    
    cmdparams = read_cmdparams()
    params.update(cmdparams)
    
    if 'configfile_fitsfiles' in params.keys():
        if os.path.exists(params['configfile_fitsfiles']):
            textparams = read_parameterfile(params['configfile_fitsfiles'])
            params.update(textparams)
    
    params.update(cmdparams)        # Repeat, in case textparams overwrites cmdparams
    
    # Add standard parameters, if they are missing, more at the end
    params['logging_arc_line_identification_residuals_hist'] = params.get('logging_arc_line_identification_residuals_hist', 'arc_line_identification_residuals_hist.png')
    params['logging_em_lines_bisector'] = params.get('logging_em_lines_bisector', 'wavelength_solution_emmission_lines_bisector.png')
    params['logging_blaze_spec'] = params.get('logging_blaze_spec', 'blaze_spectrum.pdf')
    # Added 20201022:
    params['logging_crossdispersion_shift'] = params.get('logging_crossdispersion_shift', 'crossdispersions_shift.txt')
    # Added 20201210:
    params['logging_resolution_form'] = params.get('logging_resolution_form', 'wavelength_solution_resolution_on_detector.png')
    
    list_txt = ['reference_catalog', 'use_catalog_lines', 'raw_data_file_endings', 'raw_data_mid_exposure_keys', 'raw_data_paths', 'raw_data_object_name_keys', 'cosmic_ray_settings']
    list_int = ['arcshift_range', 'order_offset', 'px_offset', 'polynom_order_traces', 'polynom_order_intertraces',
             'bin_search_apertures', 'bin_adjust_apertures', 'polynom_bck', 'dataset_rv_analysis']
    list_float = ['opt_px_range', 'px_offset_order', 'background_width_multiplier', 'sigmaclip_spectrum', 'traces_searchlimit_brightness']
    list_abs = ['arcshift_range']
    ints = ['polynom_order_apertures', 'rotate_frame']
    floats = ['max_good_value', 'catalog_file_wavelength_muliplier', 'extraction_width_multiplier', 'arcextraction_width_multiplier',
              'resolution_offset_pct', 'diff_pxs', 'maxshift_orders', 'wavelength_scale_resolution', 'width_percentile', 'raw_data_timezone_cor',
              'altitude', 'latitude', 'longitude', 'in_shift', 'extraction_shift', 'extraction_min_ident_part_of_trace_percent',
              'max_cores_used_pct', 'pxshift_between_wavesolutions']
    bools = ['flip_frame', 'update_width_orders', 'GUI', 'started_from_p3']
    text_selection = ['arcshift_side', 'extraction_precision']
    results = ['path_extraction', 'path_extraction_single', 'logging_path', 'path_reduced', 'path_rv_ceres', 'path_rv_terra', 'path_rv_serval',
               'path_harpsformat', 'master_blaze_spec_norm_filename']     #, 'configfile_fitsfiles' (excluded, as should be handled as conf.txt
    full_filenames = ['badpx_mask_filename', 'original_master_traces_filename', 'original_master_wavelensolution_filename',
                'configfile_fitsfiles', 'raw_data_file_list', 'terra_jar_file']   # deal with full filenames -> nothing to do
    texts = ['editor', 'extracted_bitpix', 'site', 'object_file', 'raw_data_imtyp_keyword', 
             'raw_data_imtyp_bias', 'raw_data_imtyp_dark', 'raw_data_imtyp_flat', 'raw_data_imtyp_trace1', 'raw_data_imtyp_blaze', 'raw_data_imtyp_trace2',
             'raw_data_exptim_keyword', 'raw_data_dateobs_keyword', 'raw_data_timezone_cor', 'blazercor_function']      # -> nothing to do
    paths, loggings = [], []
    #list_raw = []              # not necessary anymore
    for entry in params.keys():                         # make the live easier by adding entries automatically to the lists above
        if entry not in list_int and (entry.find('subframe') >= 0):                                         # add to the list
            list_int.append(entry)
        if entry not in list_txt and (entry.find('_rawfiles') > 0 or entry.find('calibs_') > -1):           # add to the list
            list_txt.append(entry)
        #if entry not in list_raw and (entry.find('_rawfiles') > 0):
        #    list_raw.append(entry)
        if entry not in paths    and (entry.find('path_') != -1 or entry.find('_path') != -1):
            paths.append(entry)
        if entry not in results  and ((entry.find('master_') == 0 or entry.find('background_') == 0) and entry.find('_filename') > 0):      # Put these files also into the result path
            results.append(entry)
        if entry not in loggings and (entry.find('logging_') == 0 and entry.find('logging_path') == -1):
            loggings.append(entry)
        if entry not in texts    and (entry.find('raw_data_imtyp_') == 0 or (entry.find('raw_data_') == 0 and entry.find('_keyword') > 0)):
            texts.append(entry)
        if entry not in list_txt and (type(params[entry]) == str):
            if params[entry].find(',') != -1:
                list_txt.append(entry)
    all_parameters = list(np.unique( list_txt + list_int + list_float + list_abs + ints + floats + bools + text_selection + results + full_filenames + paths + loggings + texts ))       # + list_raw
    trues =     ['yes', 'true', '1']
    falses =    ['no', 'false', '0']
    lefts =     ['left','l']
    rights =    ['right','r']
    centers =   ['center','centre','c']
    standards = ['standard']
    linfits =   ['linfit']
    undeclared_params = ''
    # Important settings first, as other things depend on them:
    for entry in paths:
        if entry in params.keys() and entry not in list_txt:
            params[entry] = (params[entry]+os.sep).replace(os.sep+os.sep, os.sep)      # Add a / at the end in case the user didn't
    for entry in results:
        if entry in params.keys():                                            # deal with result filenames/folders -> add result_path
            params[entry] = params['result_path'] + params[entry]
    # All parameters
    for entry in params.keys():
        if entry in list_txt+list_int+list_float:       # Split the lists
            temp = params[entry]
            for i in ['[', ']']:                        # remove the brakets
                temp = temp.replace(i,'')
            temp = [x.strip() for x in temp.split(',')] # strip the entries from leading or ending spaces
            for i in range(len(temp))[::-1]:            # remove empty entries
                if temp[i] == '':
                    del temp[i]
            params[entry] = temp
        for [list_no, funct, functxt] in [ [list_int, int, 'intergers'], [list_float, float, 'floats'], [list_abs, abs, 'numbers'] ]:
            if entry in list_no:        # deal with list of integers, floats, abs, numpy array (int and float could be also done with dtype, but then needs to be converted back to list
                for i in range(len(params[entry])):
                    try:
                        params[entry][i] = funct(params[entry][i])
                    except:
                        logger(emsg + 'Parameter "{0}" (value of "{1}") must be a list of {2}. Error occured with index {3}.'.format(entry, params[entry], functxt, i))
        for [ints_floats, funct, functxt] in [ [ints, int, 'intergers'], [floats, float, 'floats'] ]:
            if entry in ints_floats:    # deal with integers or floats
                try:
                    params[entry] = funct(params[entry])
                except:
                    logger(emsg + 'Parameter "{0}" (value of "{1}") must be a {2}'.format(entry, params[entry], functxt))
        if entry in bools:
            if params[entry].lower() in trues:      params[entry] = True
            elif params[entry].lower() in falses:   params[entry] = False
            else:
                logger(emsg + 'Parameter "{0}" (value of "{1}") must be any of the following: {2}'.format(entry, params[entry], trues+falses))
        if entry in text_selection:
            if entry == 'arcshift_side':
                if params[entry].lower() in lefts:          params[entry] = -1
                elif params[entry].lower() in rights:       params[entry] = 1
                elif params[entry].lower() in centers:      params[entry] = 0
                else:
                    logger(emsg + 'Parameter "{0}" (value of "{1}") must be one of the following: {2}'.format(entr, params[entry], lefts+rights+centers))
            if entry == 'extraction_precision':
                if params[entry].lower() in standards:   params[entry] = standards[0]
                elif params[entry].lower() in linfits:   params[entry] = linfits[0]
                else:
                    logger(emsg + 'Parameter "{0}" (value of "{1}") must be one of the following: {2}'.format(entr, params[entry], standards+linfits))
        if entry in paths:                                              # Create folders
            if type(params[entry]).__name__ in ['list']:    
                for ii in range(len(params[entry])):
                    params[entry][ii] = (params[entry][ii]+os.sep).replace(os.sep+os.sep, os.sep)      # Add a / at the end in case the user didn't
            if entry in ['raw_data_paths', 'path_ceres', 'path_serval']:         # raw_data_paths is list, hence has to be checked before
                continue
            if params[entry].lower() not in ['na'+os.sep, params['result_path']+'na'+os.sep, params['result_path']+os.sep+'na'+os.sep]:
                make_directory(params[entry])               # Create folders, if necessary
        #if entry in list_raw:                                           # deal with lists of raw data filenames -> add path
        #    for i in range(len(params[entry])):
        #        params[entry][i] = params['raw_data_path'] + params[entry][i]              # Not necessary anymore
        if entry in loggings:                                           # deal with logging filenames/folders -> add logging_path
            params[entry] = params['logging_path'] + params[entry]
        if entry not in all_parameters:
            undeclared_params += '{0}, '.format(entry)
    overdeclared_params = ''
    for entry in all_parameters:
        if entry not in params.keys():
            overdeclared_params += '{0}, '.format(entry)
    if len(overdeclared_params) > 0:
        logger('Warn: The following parameters are expected, but do not appear in the configuration files. : {1}{0}\t!!! The missing parameters might cause crashes later !!!'.format(os.linesep, overdeclared_params[:-2]))
    if len(undeclared_params) > 0:
        logger('Warn: The following parameters appear in the configuration files, but the programm did not expect them: {0}'.format(undeclared_params[:-2]))
     
    # Use standard parameters, if not given; more at the beginning, e.g. when modification of paths is necessary
    params['use_cores'] = int(multiprocessing.cpu_count()*params.get('max_cores_used_pct',80)/100.)     # Use 80% of the CPU cores
    params['dataset_rv_analysis'] = params.get('dataset_rv_analysis', [5, 6])
    params['pxshift_between_wavesolutions'] = params.get('pxshift_between_wavesolutions', 0)
    params['cosmic_ray_settings'] =  params.get('cosmic_ray_settings', ['deepCR', 'ACS-WFC-F606W-2-32', '0.999'])
    params['traces_searchlimit_brightness'] =  params.get('traces_searchlimit_brightness', 50)
    return params

def make_directory(directory, errormsg='Warn: Cannot create directory {0}'):
    if not os.path.exists(directory):
        try:                                                    # Create folders, if necessary
            os.makedirs(directory)
        except:
            logger(errormsg.format(directory))

def log_returncode(code, explain=''):
    """
    Logs a message if a subprocess fails
    """
    if code != 0 and code != None:
        logger('Warn: Subprocess returned with error code {0}. {1}'.format(code, explain))
    return code

#def update_calibration_memory(key, value):         # Not needed anymore with global calimages in hiflex.py and in procedures here
#    """
#    Add new information to the global variable calimages, which is not accessable from other python files
#    :param key: string, key for the dictionary
#    :param value: string, number, array, or anything: value for the dictionary
#    """
#    global calimages
#    calimages[key] = copy.deepcopy(value)

def round_sig(value, sig=5):
    """
    Rounds to a significant number, e.g. to have 5 digits
    :param value: number
    :param sig: int, number of dgits allowed
    :return rounded: float, rounded to the number of sig digits
    """
    if np.isnan(value) or np.isinf(value):
        return value
    if value == 0.0:
        return value
    rounded = round(value, sig-int(np.floor(np.log10(abs(value))))-1)
    
    return rounded

def sort_for_multiproc_map(inputlist, number_cores):
    """
    multiprocessing.Pool.map will give the first part of the list to the first subprocess, second part to second subprocess and so on.
    So in order to get data reduced in a more linear way, resorting the list will be helpful
    :param inputlist: list of anthing
    :param number_cores: integer
    :return outlist: same length as inputlist, only resorted (e.g. number_cores=3: inputlist=[1,'a',3,4,5,6,7,8] -> outlist=[1,4,7,'a',5,8,3,6]
    """
    return inputlist
    # !!! Doesn't work as expected
    outlist = []
    for ii in range(number_cores): 
        outlist.append([])       # [[]]*number will create empty arrays that are all the same, so if one is filled, also the others are filled
    # Resort
    for ii in range(len(inputlist)): 
        outlist[ii%number_cores].append(inputlist[ii])
    # Combine into one list again
    for entry in outlist[1:]: 
        outlist[0] += entry
    outlist = outlist[0]
    
    return outlist

def read_text_file(filename, no_empty_lines=False, warn_missing_file=True):
    """
    Read a textfile and put it into a list, one entry per line
    Empty lines means they consist only of tabs and spaces
    """
    text = []
    if os.path.isfile(filename) == True:
        text1 = open(filename,'r').readlines()
        for line in text1:
            line = line.replace(os.linesep, '')
            linetemp = line.replace('\t', '')
            linetemp = linetemp.replace(' ', '')
            if ( line == '' or linetemp == '') and no_empty_lines:
                continue
            text.append(line)
    elif warn_missing_file:
        logger('Warn: File {0} does not exist, assuming empty file'.format(filename))    
    return text

def add_text_file(text, filename):
    """
    Adds a line or lines of text to a file without checking if the line already exists
    """
    with open(filename, 'a') as file:
        file.write(text+os.linesep)

def add_text_to_file(text, filename, warn_missing_file=True):
    """
    If the text is not yet in the file, the text is added to the file
    """
    oldtext = read_text_file(filename, warn_missing_file=warn_missing_file)
    exists = False
    for line in oldtext:
        if line.find(text) != -1:
            exists = True
            break
    if exists == False:
        add_text_file(text, filename)
        
def convert_readfile(input_list, textformats, delimiter='\t', replaces=[], ignorelines=[], expand_input=False, shorten_input=False, ignore_badlines=False, replacewithnan=False):
    """
    Can be used convert a read table into entries with the correct format. E.g integers, floats
        Ignories the lines which have less entries than entries in textformats
        replacements will be done before 
    :param input_list: 1d list or array of strings from a read file
    :param textformats: 1d list or array of formats, e.g. [str, str, int, float, float]. 
                        Individual entries can consist of a list of entries, then the conversion will be done in the order given, e.g. [str, [float, int], [float, '%Y-%m-%dT%H:%M:%S.%f'], float]
                        If a string is given this will be used to convert it into a datetime object
                        If a float, int, or str is run on a datetime object, then the datetime object will be converted into a number
    :param delimiter: string, used to split each line into the eleements
    :param replaces: 1d list or array of strings and/or lists of 2 strings, contains the strings which should be replaced by '' (if string) or by the second entry (if list)
    :param ignorelines: List of strings and/or lists. A list within ignorelines needs to consist of a string and the maximum position this string can be found.
                        If the string is found before the position, the entry of the input list is ignored
    :param expand_input: bool, if the line in the input_list contains less elements than textformats the line will be filled up with ''s
    :param shorten_input: bool, if the  line in the input_list contains more elements than textformats the line will be shortened to len(textformats)
    :param ignore_badlines: bool, if True will raise Warnings, if False will raise Errors (which will terminate code)
    :param replacewithnan: bool, if True and conversion with textformats is not possible will replace with NaN, otherwise will raise Warning/Error
    :retrun result_list: 2d list with numbers or strings, formated acording to textformats
    """
    # Make the textformats, replaces, and ignorelines consistent
    for index in range(len(textformats)):
        if type(textformats[index]) != list:
            textformats[index] = [textformats[index]]     # convert the single entries into list, e.g. make a list of lists
    for index in range(len(ignorelines)):
        error = False
        if type(ignorelines[index]) == str:     # only one entry
            ignorelines[index] = [ignorelines[index], 1E10]
        elif type(ignorelines[index]) == list:  # list with two entries: string, and position after which the string will be ignored
            if len(ignorelines[index]) != 2:
                error = True
        else:
            error = True
        if error:
            logger('Error: Programming error: ignorelines which where hand over to convert_readfile are wrong. '+\
                   'It has to be a list consisting of strings and/or 2-entry lists of string and integer. Please check ignorelines: {0}'.format(ignorelines))
    for index in range(len(replaces)):
        error = False
        if type(replaces[index]) == str:     # only one entry
            replaces[index] = [replaces[index], '']
        elif type(replaces[index]) == list:  # list with two entries: search-string and replace-string
            if len(replaces[index]) != 2:
                error = True
        else:
            error = True
        if error:
            logger('Error: Programming error: replaces which where hand over to convert_readfile are wrong. '+\
                   'It has to be a list consisting of strings and/or 2-entry lists of strings. Please check replaces: {0}'.format(replaces))
    # Convert the text
    error = {False:'Error:',True:'Warn:'}[ignore_badlines]
    result_list = []
    for entry in input_list:
        notuse = False
        for ignore in ignorelines:
            if entry.find(ignore[0]) <= ignore[1] and entry.find(ignore[0]) != -1:
                notuse = True
                break
        if notuse:
            continue
        for replce in replaces:
            entry = entry.replace(replce[0],replce[1])
        entry = entry.split(delimiter)
        if len(entry) < len(textformats):           # add extra fields, if not enough
            if expand_input:
                entry += [''] * ( len(textformats) - len(entry) )
            else:
                continue
        for index in range(len(textformats)):
            for textformat in textformats[index]:
                if type(textformat) == type:
                    if type(entry[index]) == datetime.datetime:
                        epoch = datetime.datetime.utcfromtimestamp(0)
                        entry[index] = (entry[index] - epoch).total_seconds()         # (obsdate - epoch) is a timedelta
                    try:
                        entry[index] = textformat(entry[index])
                    except:
                        if replacewithnan:
                            entry[index] = np.nan
                        else:
                            logger(error+' Cannot convert {0} into a {1}. The problem happens in line {2} at index {3}.'.format(entry[index], textformat, entry, index))
                elif type(textformat) == str:
                    try:
                        entry[index] = datetime.datetime.strptime(entry[index], textformat)
                    except:
                        logger(error+' Cannot convert {0} into a datetime object using the format {1}. The problem happens in list line {2} at index {3}.'.format(entry[index], textformat, entry, index))
                else:
                    logger('Error: Programming error: textformats which where hand over to convert_readfile are wrong. It has to be a type or a string')
        if shorten_input and len(entry) > len(textformats):
            del entry[len(textformats):]
        result_list.append(entry)
    return result_list

def read_badpx_mask(params, imageshape): 
    """
    Reads the bad-pixel-mask and applies the corrections to it
    :param params: Dictionary with all the parameters: 
                badpx_mask_filename is used in order to get the bad px mask. If it doesn't exist, than all pixels are fine
                calibs is used in order to if a subframe of the image is used. In this case subframe is used to find the area used
    :param imageshape: tuple of two integers, gives the shape of the image to which the badpixel mask should be applied
    :return badpx_mask: 2d numpy array of the bad pixel mask
    """
    filename = params['badpx_mask_filename']
    if os.path.isfile(filename) == True:
        badpx_mask = np.array(fits.getdata(filename), dtype=np.float64)
        badpx_mask = rotate_flip_frame(badpx_mask, params )
        if badpx_mask.shape != imageshape:
            if subframe in params.keys():
                subframe = params['subframe']
                if len(subframe) >= 4:
                    badpx_mask = badpx_mask[subframe[2]: subframe[0]+subframe[2], subframe[3]: subframe[1]+subframe[3]]
        if badpx_mask.shape != imageshape:
            logger('Warn: The loaded badpixel mask {0} has a different format than the image file: {1} vs {2} (image file)'.format(filename, badpx_mask.shape, imageshape))
        else:
            logger('Info: badpixel mask loaded: {0}'.format(filename))
    else:
        logger('Warn: No badpixel mask found ({0}). Assuming no bad pixels'.format(filename))
        badpx_mask = np.ones(imageshape)
        
    return badpx_mask

"""def read_background(params, filename):               # not used anymore as not really useful
    "#""
    Reads the background map and applies the corrections to it
    :param params: Dictionary with all the parameters
    :param filename: path and name of the background file
    :return background: 2d numpy array of the background, normalised to 1s exposure time
    "#""
    imtype = 'background'
    params['calibs'] = params['calibs_read']
    if '{0}_calibs_read'.format(imtype) in params.keys():
        params['calibs'] = params['{0}_calibs_read'.format(imtype)]
    im, im_head = read_file_calibration(params, filename)
    exptime = im_head[params['raw_data_exptim_keyword']]
    background = im/exptime
    
    return background"""

def warn_images_not_same(ims, names):
    """
    Check that two images have the same size
    :param ims: list or array of arrays
    :param names: list of strings, same length as ims
    """
    if len(ims) <= 1:
        return
    if len(ims) != len(names):
        logger('Warn: len(ims) != len(names), that seems like a coding error.')
    problems = ''
    for i in range(len(ims)-1):
        for j in range(i,len(ims)):
            if ims[i].shape != ims[j].shape:
                problems += '\tImage: {0} ({2}) and Image {1} ({3})\n'.format(names[i], names[j], ims[i].shape, ims[j].shape)
    if problems != '':
        logger('Error: The following images have not the same size, but should have. '+\
               'This is most likely caused by a missing "subframe" in one of the parameters "calib*". Please check.\n {0}'.format(problems[:-1]) )
    return

def read_file_calibration(params, filename, level=0, realrun=True):
    """
    Reads the filename and applies the calibrations as given in params['calibs']
    This works also if the files are stored with the following header: BITPIX = 16, BZERO = 32768
    :param params: Dictionary with all the parameters
    :param filename: path and name of the file
    :return: 2D numpy array of the file, and the header of the file read
    """
    global calimages
    im, im_head = read_fits_file(filename, realrun=realrun)
    if realrun:
        ims = im.shape
        while len(ims) > 2:
            if ims[0] == 1:                 # e.g. MRES from NARIT
                im = im[0,:,:]
            elif ims[-1] == 1:
                im = im[:,:,0]
            else:
                logger(('Error: The file is stored in a multi-demensional array, which I do not know how to handle. The size of the image is {0}. '+\
                        '\tThis requires a small adjustment to the code in procedure read_file_calibration.').format(ims))
            ims = im.shape
    if params['raw_data_exptim_keyword'] in im_head.keys():
        exptime = im_head[params['raw_data_exptim_keyword']]        #saved as float
    else:
        exptime = 0.0
        logger('Warn: Cannot find the raw_data_exptim_keyword = {0} in the header. Assuming 0 seconds.'.format(params['raw_data_exptim_keyword'] ))
    if realrun:
        logger('Info: image loaded: {0}'.format(filename))
        im = rotate_flip_frame(im, params )
        filename_s = filename[max(0,len(filename)-(80-25-1)):]               # Shorten the filename so it fits into the header, 21 failed, 25 worked, between not tested
        im_head['HIERARCH HiFLEx orig'] = filename_s
        filename_nopath = filename_s.rsplit(os.sep,1)[-1]
        im_head['HIERARCH HiFLEx orid'] = filename_nopath
        calimages['{0}_saturated'.format(filename_s)] = np.where( im > params['max_good_value'] )         # Find the saturated pixel in the original image, consists of (x,y) with x and y are arrays
    for entry in params['calibs']:
        if entry == '':
            continue
        logtxt, headtxt = [], []
        if entry.lower().find('subframe') > -1 and realrun:
            subframe = copy.deepcopy(params[entry])
            if len(subframe) >= 4:
                # First: extend the image, if necessary
                if subframe[2] < 0:
                    im = np.vstack(( np.zeros((-subframe[2], im.shape[1])), im ))
                    subframe[2] = 0
                if subframe[3] < 0:
                    im = np.hstack(( np.zeros((im.shape[0], -subframe[3])), im ))
                    subframe[3] = 0
                if im.shape[0] < subframe[0]+subframe[2]:
                    im = np.vstack(( im, np.zeros((subframe[0]+subframe[2]-im.shape[0], im.shape[1])) ))
                if im.shape[1] < subframe[1]+subframe[3]:
                    im = np.hstack(( im, np.zeros((im.shape[0], subframe[1]+subframe[3]-im.shape[1])) ))
                # Last: cut the image
                if im.shape != (subframe[0],subframe[1]):                   # only apply subframe if the file doesn't have the size already
                    im = im[subframe[2]: subframe[0]+subframe[2], subframe[3]: subframe[1]+subframe[3]]
            ims = im.shape
            logger('Info: {1}: subframe applied: {0} ({2})'.format(entry, level, params[entry]))
            calimages['{0}_saturated'.format(filename_s)] = np.where( im > params['max_good_value'] )         # Find the saturated pixel in the original image, do it again if the image was cut
            im_head['HIERARCH HiFLEx redu{0}a'.format(level)] = 'Subframe: {0}'.format(entry)
        elif entry.lower().find('bias') > -1:
            if entry not in calimages:
                 create_image_general(params, entry, level=level+1)
            if realrun:
                warn_images_not_same([im, calimages[entry]], [filename,entry])
                if np.percentile(calimages[entry], 90) > 2000 or np.percentile(calimages[entry], 10) < -100:
                    logger('Warn: The bias ({0}) has unphysical values: 10%-percentile = {1}, 90%-percentile = {2}'.format(entry, np.percentile(calimages[entry], 10), np.percentile(calimages[entry], 90)))
                im = im - calimages[entry]
                logtxt, headtxt = ['bias correction applied'], ['redu{0}b'.format(level), 'Bias']
        elif entry.lower().find('dark') > -1:
            if entry == 'dark':             #only add exposure time, if just dark is given
                entry = entry+str(exptime)
            if entry not in calimages:
                 create_image_general(params, entry, level=level+1)
            if realrun:
                warn_images_not_same([im, calimages[entry]], [filename,entry])
                if np.percentile(calimages[entry], 90) > 2000 or np.percentile(calimages[entry], 10) < -100:
                    logger('Warn: The dark ({0}) has unphysical values: 10%-percentile = {1}, 90%-percentile = {2}'.format(entry, np.percentile(calimages[entry], 10), np.percentile(calimages[entry], 90)))
                im = im - calimages[entry]
                logtxt, headtxt = ['dark correction applied'], ['redu{0}c'.format(level), 'Dark']
        elif entry.lower().find('rflat') > -1:
            if entry not in calimages:
                 create_image_general(params, entry, level=level+1)
            if realrun:
                warn_images_not_same([im, calimages[entry]], [filename,entry])
                im = im / (calimages[entry] / np.median(calimages[entry]) )
                logtxt, headtxt = ['flat correction applied with normalised flat (rflat)'], ['redu{0}d'.format(level), 'Flat']
        elif entry.lower().find('badpx_mask') > -1 and realrun:
            if entry not in calimages:
                calimages[entry] = read_badpx_mask(params, ims)
            warn_images_not_same([im, calimages[entry]], [filename,entry])
            im = im * calimages[entry]
            nonzeroind = np.nonzero(1-calimages[entry])       #find the indexes where the badpx mask is zero
            for i in range(len(nonzeroind[0])):         #find a suitable amound of surrounding pixels
                for j in range(1,5):                    # try to find it in the surrounding 8, 24, 48, 80, 120 pixel
                    section = im[max(0,nonzeroind[0][i]-j) : min(ims[0],nonzeroind[0][i]+j+1) , max(0,nonzeroind[1][i]-j) : min(ims[1],nonzeroind[1][i]+j+1)]
                    section = section[section!=0]           #only use areas <> 0 (exclude bad pixel
                    if len(section) >= 5:
                        break
                if len(section) == 0:
                    logger('Warn: cannot replace bad pixel ({0}, {1}) with surrounding area in {2}'.format(nonzeroind[0][i],nonzeroind[1][i],filename))
                else:
                    im[nonzeroind[0][i],nonzeroind[1][i]] = np.median(section)  #replace bad px with the median of each surrounding area
            logger('Info: {1}: badpx correction applied: {0}'.format(entry, level))
            im_head['HIERARCH HiFLEx redu{0}f'.format(level)] = 'Bad-pixel-mask: {0}'.format(entry)
        elif entry.lower().find('cosmic_ray') > -1 and realrun:
            cr_setting = params['cosmic_ray_settings']
            if cr_setting[0].lower() == 'deepcr' and deepCR is not None:
                if multiprocessing.current_process().name == 'MainProcess':     n_jobs = params['use_cores']
                else:                                                           n_jobs = 1
                mdl = deepCR(mask=cr_setting[1], inpaint=cr_setting[1], device="CPU")     # mdl could be initialised earlier
                #logger('Step: Cosmic ray removal')
                mask, im_clean = mdl.clean(im, threshold=float(cr_setting[2]), segment=True, n_jobs=n_jobs)            # best threshold is highest value that generate mask covering full extent of CR; choose threshold by visualizing outputs.
                pos = np.where(mask)
                #for ii in range(pos[0].shape[0]):
                #    print(filename, ii, pos[0][ii],pos[1][ii])
                fname = filename.rsplit(os.sep,1)
                if mask.any() and os.path.exists(params['path_reduced']) and params['path_reduced'].lower() != 'na'+os.sep:      # at least one entry is True/1
                    im_head['Comment'] = 'In this file {0} cosmic rays were marked.'.format(np.sum(mask)) +\
                                         ' The first image shows the cleaned image, the second the mask and the third the uncleaned image'
                    save_im_fits(params, [im_clean, mask, im], im_head,  params['path_reduced']+fname[-1].replace('.fit','_cr.fit'))
                im = im_clean
                logger('Info: {1}: cosmic ray correction applied: {0}'.format(entry, level))
                im_head['HIERARCH HiFLEx redu{0}g'.format(level)] = 'Cosmic ray correction: {0}'.format(entry)
                im_head['HIERARCH HiFLEx deepCR'] = '{0} pixel from cosmic rays'.format(np.sum(mask))
        elif entry.lower().find('background') > -1 and realrun:         # was localbackground before 20200108, this search covers *background
            if 'sci_trace' in calimages.keys() and 'cal_trace' in calimages.keys():
                logger('Step: Performing the background fit')
                sci_tr_poly, xlows, xhighs, widths = calimages['sci_trace']
                cal_tr_poly, axlows, axhighs, awidths = calimages['cal_trace']
                bck_px_sci = find_bck_px(im, sci_tr_poly, xlows, xhighs, widths, params['background_width_multiplier'][0])
                bck_px_cal = find_bck_px(im, cal_tr_poly, axlows, axhighs, awidths, params['background_width_multiplier'][1])
                bck_px = bck_px_sci * bck_px_cal        # mask with 0 for all the traces and 1 for background data
                bad_values = ( im*bck_px > np.percentile(im[bck_px==1],95) )        # exclude brightest data in the background data (emission lines or borders of the traces
                bck_px[bad_values] = 0
                
                # Some deviation at the red side with not many lines, computational heavy, binning in dispersion direction made it faster
                bck_im = fit_2d_image(im, params['polynom_bck'][1], params['polynom_bck'][0], w=bck_px, binning=[int(im.shape[0]/200.),1])
                #plot_img_spec.plot_image(im, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'orig')
                #plot_img_spec.plot_image(bck_px, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'bck_px')
                #plot_img_spec.plot_image(bck_im, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'bck_im')
                #plot_img_spec.plot_image(im-bck_im, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'diff')
                plot_img_spec.plot_image((im - bck_im)*bck_px, \
                                        [params['logging_path']+'background_subtracted-'+filename.rsplit(os.sep,1)[-1]+'.png'],\
                                         1, False, [0.05,0.95,0.95,0.05], 'difference between image and background fit')
                im = im - bck_im
                bck_px[ bck_px==0 ] = np.nan
                bck_noise_std, bck_noise_var = measure_background_noise(im * bck_px)
                if np.isnan(bck_noise_var):
                    bck_noise_var = -1
                logger('Info: {1}: background correction applied: {0}'.format(entry, level))
                im_head['HIERARCH HiFLEx redu{0}f'.format(level)] = 'Background: {0}'.format(entry)
                im_head['HIERARCH HiFLEx BCKNOISE'] = round(bck_noise_std,8)
                im_head['HIERARCH HiFLEx BNOISVAR'] = (round(bck_noise_var,8), 'Variation of noise through image')             # Background noise variation can be very high, because some light of the traces remains
            else:
                logger('Warn: Could not apply the calibration step {0} because the science and/or calibration traces are not yet known.'.format(entry))
        elif entry.lower().find('combine_sum') > -1 or entry.lower().find('combine_mean') > -1 or entry.lower().find('normalise') > -1:
            'nothing to do, as for a different step'
        elif realrun:
            logger('Warn: do not know what to do with this correction: {0}'.format(entry))
        if len(logtxt) > 0 and len(headtxt) > 0 and realrun:
            #print "np.where(np.isnan(calimages[entry])), calimages[entry].shape", np.where(np.isnan(calimages[entry])), calimages[entry].shape, np.median(calimages[entry]), np.nanmedian(calimages[entry]), np.mean(calimages[entry]), np.nanstd(calimages[entry], ddof=1, axis=None), np.where( np.isinf(calimages[entry]) )
            im_median, im_std = int(round(np.nanmedian(calimages[entry], axis=None))), int(round(np.nanstd(calimages[entry], ddof=1, axis=None)))
            logger('Info: {4}: {3}: {0} (median={1}, std={2})'.format(entry, im_median, im_std, logtxt[0], level))
            im_head['HIERARCH HiFLEx '+headtxt[0]] = '{3}: {0}, median={1}, std={2}'.format(entry, im_median, im_std, headtxt[1])
    #logger('Info: image loaded and processed: {0}'.format(filename))
    if os.path.exists(params['path_reduced']) and params['path_reduced'].lower() != 'na'+os.sep and realrun:       # Save the reduced image
        fname = filename.rsplit(os.sep,1)
        save_im_fits(params, im, im_head,  params['path_reduced']+fname[-1])
    return im, im_head

def create_image_general(params, imtype, level=0, realrun=True):
    """
    Reads or creates the imtype file. If the key and file for master_<imtype>_filename exists the file is read, otherwise the file is created by combining the <imtype>_rawfiles
    :param params: Dictionary with all the parameters
    :param imtype: type of images, e.g. flat, dark5.0, bias
    :param level: Sometimes the calibration images has to be created before it can be applied, to destinguish between this and the current image, the level can be increased by one
    :return im: 2d array of the combined file
    :return im_head: header of the last read file
    """
    global calimages
    mem = psutil.virtual_memory()                   # svmem(total=33221091328, available=28485840896, percent=14.3, used=4202041344, free=25513508864, active=..., inactive=..., buffers=..., cached=.., shared=...)
    loaded = False
    if 'master_{0}_filename'.format(imtype) in params.keys():
        if params['master_{0}_filename'.format(imtype)] != '' and os.path.isfile(params['master_{0}_filename'.format(imtype)]) == True:
            if realrun:
                logger('Info: Using existing {0}: {1}'.format(imtype,params['master_{0}_filename'.format(imtype)]))
            params['calibs'] = params['calibs_read']
            if '{0}_calibs_read'.format(imtype) in params.keys():
                params['calibs'] = params['{0}_calibs_read'.format(imtype)]
            im, im_head = read_file_calibration(params, params['master_{0}_filename'.format(imtype)], level=level, realrun=realrun)
            loaded = True
    if loaded == False:
        if '{0}_calibs_create'.format(imtype) not in params.keys():
            if 'standard_calibs_create' not in params.keys() and realrun:
                logger('Error: Missing entry in the configuration file. Neigther "{0}_calibs_create" nor "standard_calibs_create" is given. Please update the configuration file(s).'.format(imtype))
            params['{0}_calibs_create'.format(imtype)] = params['standard_calibs_create']
        for i in range(len(params['{0}_calibs_create'.format(imtype)])):                                                                    # make it safe from different user input
            params['{0}_calibs_create'.format(imtype)][i] = params['{0}_calibs_create'.format(imtype)][i].lower()
            params['{0}_calibs_create'.format(imtype)][i] = params['{0}_calibs_create'.format(imtype)][i].replace('normaliz', 'normalis')
            params['{0}_calibs_create'.format(imtype)][i] = params['{0}_calibs_create'.format(imtype)][i].replace('normalisation', 'normalise')
        im, med_fluxes, std_fluxes = None, [], []
        if '{0}_rawfiles'.format(imtype) not in params.keys() and realrun:
            logger('Error: The list of raw files for image type {0} is not defined in the configuration. Please check the configuration files.'.format(imtype))
        if len(params['{0}_rawfiles'.format(imtype)]) == 0 and realrun:
            logger('Error: The list of raw files for image type {0} is empty. Please check the configuration files.'.format(imtype))
        num_imgs = len(params['{0}_rawfiles'.format(imtype)])                # how many images are expected
        header_updates = np.zeros((num_imgs,2))
        for im_index, imf in enumerate(params['{0}_rawfiles'.format(imtype)]):                   # Only works for maximum 40 images on neils machine
            params['calibs'] = params['{0}_calibs_create'.format(imtype)]       # get's overwritten when other files are being read
            img, im_head = read_file_calibration(params, imf, level=level, realrun=realrun)
            if realrun:
                im_head, obsdate_midexp, obsdate_mid_float, jd_midexp = get_obsdate(params, im_head)    # unix_timestamp of mid exposure time
                header_updates[im_index,:] = [im_head['HIERARCH HiFLEx EXPOSURE'], obsdate_mid_float]         # !!! Improve this calculation and write in the header so it can be used later by get_obsdate 
                med_flux = np.median(img, axis=None)
                med_fluxes.append(med_flux)
                std_fluxes.append(np.std(img, axis=None, ddof=1))
                if 'normalise' in params['{0}_calibs_create'.format(imtype)]:
                    img = img/(med_flux+0.0)
                if im is None:                                                # Initiate the array with correct precission to avoid swapping
                    if num_imgs * np.prod(img.shape) * 2 > mem[1] * 0.49:
                        prec = np.float16
                        logger('Warn: The ammount of pictures will most likely cause swapping, which might cause the system to become unresponsive.')
                    elif num_imgs * np.prod(img.shape) * 4 > mem[1] * 0.49:
                        prec = np.float16  
                    elif num_imgs * np.prod(img.shape) * 8 > mem[1] * 0.49:
                        prec = np.float32  
                    else:
                        prec = np.float64
                    #logger('Test: Precision {0} is used'.format(prec))
                    im = np.zeros( (num_imgs, img.shape[0], img.shape[1]) ).astype(prec)
                    #print 'Created:', im.dtype, im.itemsize, im.nbytes, sys.getsizeof(im), im.nbytes/7979408000., mem
                im[im_index,:,:] = img
                #print im.dtype, im.itemsize, im.nbytes, sys.getsizeof(im), im.nbytes/7979408000.
        if realrun:
            for i in range(len(med_fluxes)):
                im_head['HIERARCH HiFLEx NORM_{0}'.format(i)] = (med_fluxes[i], 'Median flux in image {0}'.format(i))
            for i in range(len(std_fluxes)):
                im_head['HIERARCH HiFLEx STDV_{0}'.format(i)] = (round(std_fluxes[i],5), 'Stdev of flux')
            if 'combine_mean' in params['{0}_calibs_create'.format(imtype)]:
                im = combine_sum(im)/(len(im)+0.0)
                im_head['HIERARCH HiFLEx redu07'] = 'Average of {0} images'.format(len(params['{0}_rawfiles'.format(imtype)]))
                exposure_time = np.mean(header_updates[:,0])                     # Average of the exposure times
            elif 'combine_sum' in params['{0}_calibs_create'.format(imtype)]:
                im = combine_sum(im)
                im_head['HIERARCH HiFLEx redu07'] = 'Sum of {0} images'.format(len(params['{0}_rawfiles'.format(imtype)]))
                exposure_time = np.sum(header_updates[:,0])                      # Sum of the exposure times
            else:           # Median combine
                im = combine_median(im)
                #imt = combine_median(im)
                #aa = np.where(np.isinf(imt))
                #for ii in range(aa[0].shape):
                #    print aa[0],aa[1], im[:,aa[0],aa[1]]
                im_head['HIERARCH HiFLEx redu07'] = 'Median of {0} images'.format(len(params['{0}_rawfiles'.format(imtype)]))
                exposure_time = np.median(header_updates[:,0])                   # Median of the exposure times
            if 'normalise' in params['{0}_calibs_create'.format(imtype)]:
                norm_factor = np.median(med_fluxes)
                im = im * norm_factor
                im_head['HIERARCH HiFLEx NORM_MED'] = norm_factor
            im_head['HIERARCH HiFLEx MID_'+params['raw_data_dateobs_keyword']] = datetime.datetime.utcfromtimestamp(np.median(header_updates[:,1])).strftime('%Y-%m-%dT%H:%M:%S.%f')
            im_head['HIERARCH HiFLEx '+params['raw_data_exptim_keyword']] = exposure_time
            first, last = np.argmin(header_updates[:,1]), np.argmax(header_updates[:,1])
            first = header_updates[first,1]-header_updates[first,0]/2.
            last  = header_updates[last, 1]+header_updates[last, 0]/2.
            im_head['HIERARCH HiFLEx BEGIN FIRST'] = datetime.datetime.utcfromtimestamp(first).strftime('%Y-%m-%dT%H:%M:%S.%f')
            im_head['HIERARCH HiFLEx END LAST']    = datetime.datetime.utcfromtimestamp(last ).strftime('%Y-%m-%dT%H:%M:%S.%f')
            im_head['HIERARCH HiFLEx EXP_RANGE']   = (last - first, 'sec, from BEGIN to END')
            filename = params['master_{0}_filename'.format(imtype)]
            filename_s = filename[max(0,len(filename)-(80-25-1)):]               # Shorten the filename so it fits into the header, 21 failed, 25 worked, between not tested
            im_head['HIERARCH HiFLEx orig'] = filename_s
            filename_nopath = filename_s.rsplit(os.sep,1)[-1]
            im_head['HIERARCH HiFLEx orid'] = filename_nopath
            if 'master_{0}_filename'.format(imtype) in params.keys():
                if params['master_{0}_filename'.format(imtype)] != '':
                    save_im_fits(params, im, im_head,  params['master_{0}_filename'.format(imtype)])
    
    if realrun:    
        calimages[imtype] = im
        calimages['{0}_head'.format(imtype)] = im_head
    return im, im_head

def get_minimum_data_type(arr, allow_unsigned=True):
    """
    Finds the smallest necessary data type (in terms of memory) to handle the data in an array by checking, if numbers can be intergers or not and then searching for the required precession
    :param arr: list or array of numbers
    :param allow_unsigned: an unsigned data type might be less storage dependent, but not all programs will be able to handle it
    :return arr: array of numbers, using the smallest possible data type
    """
    arr = np.array(arr)
    if np.sum(np.abs(arr - arr.astype(int)), axis=None) == 0.:
        if np.min(arr, axis=None) >= 0 and allow_unsigned == True:          # Unsigned integers
            if np.max(arr, axis=None) <= 255:
                arr.dtype = np.uint8
            elif np.max(arr, axis=None) <= 65535:
                arr.dtype = np.uint16
            elif np.max(arr, axis=None) <= 4294967295:
                arr.dtype = np.uint32
            elif np.max(arr, axis=None) <= 18446744073709551615:
                arr.dtype = np.uint64
        else:
            if np.max(np.abs(arr), axis=None) <= 127:
                arr.dtype = np.int8
            elif np.max(np.abs(arr), axis=None) <= 32767:
                arr.dtype = np.int16
            elif np.max(np.abs(arr), axis=None) <= 2147483647:
                arr.dtype = np.int32
            elif np.max(np.abs(arr), axis=None) <= 9223372036854775807:
                arr.dtype = np.int64
    else:
        if (arr == arr.astype(np.float16) ).all():      # that doesn't work, as float16(3.14) == 3.140625, but I can't think about a better way at the moment
            arr.dtype = np.float16
        elif (arr == arr.astype(np.float32) ).all():
            arr.dtype = np.float32
    return arr

def save_obj(obj, name ):
    try:
        with open(name + '.pkl', 'wb') as f:
            pickle.dump(obj, f, 0)
    except:
        logger('Warn: Cannot save {0}.'.format(name + '.pkl'))

def load_obj(name ):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

def read_fits_file(filename, dtype=np.float64, realrun=True):
    if os.path.isfile(filename) == False:
        logger('Error: File {0} is missing.'.format(filename))
    hdu = fits.open(filename)
    founddata = False
    for ii in range(len(hdu)):      # Use all headers but only the first data set
        if ii == 0:
            im_head = hdu[ii].header
        else:
            im_head.extend(hdu[ii].header)
        if len(hdu[ii].shape) >= 2 and not founddata and realrun:
            im = np.array(hdu[ii].data, dtype=dtype)
            founddata = True        
    if not founddata and realrun:
        logger('Error: Found no useful data in file {0} . This probably requires a modification in the code. header.info() gives:\n{1}'.format(filename, hdu.info() ))
    if not realrun:       # just open the header
        im = np.array([])
        
    return im, im_head

def save_im_fits(params, im, im_head, filename):
    """
    Saves an image to a fits file                       # This procedure can possibly be combined with save_multispec
    :param params: Dictionary with all the parameters
    :param im: 2d array to be written into fits file
    :param im_head: header of the 2d array
    :param filename: filename in which to save the fits file
    """
    if len(filename.rsplit(os.sep,1)) == 1:     # no path is in the filename
        logger('Warn: no folder to save {0} was given, using the current folder ({1}).'.format( filename, os.getcwd() ))
    elif not os.path.exists(filename.rsplit(os.sep,1)[0]):
        logger('Error: Folder to save {0} does not exists.'.format(filename))
    if type(im).__name__ == 'list':
        for ii in range(len(im)):
            im[ii] = rotate_flip_frame(im[ii], params, invert=True)
        im = np.array(im)
    elif len(im.shape) == 2:
        im = rotate_flip_frame(im, params, invert=True)
    #im = get_minimum_data_type(im, allow_unsigned=False)   #That doen't make the fits files smaller for uint16, int16, or float32. Additionally, plotting a file with uint16 or int16 with ds9 or gaia doesn't show the data correctly
    fits.writeto(filename, im, im_head, overwrite=True)
    logger('Info: image saved: {0}'.format(filename))

def read_fits_width(filename):
    """
    Reads fits file containign a table, which defines the individual orders by parameters of a polynomial fit)
    :param filename: string, location and file name of the file containing the polynomial fits
    fits file should look like the following:
     Order  cen_center cen_poly5 ... cen_poly0 left_center left_poly5 ... left_poly0 right_center right_poly5 ... right_poly0   low_x   high_x   width_left   width_right   width_gauss  
    float64  float64    float64  ...  float64    float64     float64  ...   float64   float64     float64     ...   float64    float64 float64    float64       float64       float64    
    ------- ---------- --------- ... --------- ----------- ---------- ... ---------- ------------ ----------- ... ------------ ------- ------- ------------- ------------- --------------
        0.0      601.0 -5.54e-14 ... 22.010051       601.0 6.3123e-13 ... 21.4594466        601.0 1.06323e-12 ... 26.154773181     0.0   569.0 1.68522324742  2.8547112791   1.3453208386
    ...
       54.0      601.0  1.04e-14 ... 499.13391       601.0 4.8593e-14 ... 497.680562        601.0 -5.3207e-13 ... 500.67363664   379.0  1069.0 1.45191600068 1.42653756074  1.16066630044
    from which the polynomials are calculated (y for cen, left, and right):
        y_poly5*(x-y_center)**4 + y_poly4*(x-y_center)**(3) + ... + y_poly1*(x-y_center) + y_poly0
    
    #fits file should look like the following (old version):
    #Order          6                  5          ...       0        low_x   high_x  width_left width_right width_gauss
    #float64      float64            float64       ...    float64    float64 float64    float64     float64     float64
    #------- ------------------ ------------------ ... ------------- ------- ------- ---------- ----------- -----------
    #    1.0  7.80380819018e-20 -7.28455089523e-16 ... 1086.73647399     0.0  3028.0       13.5        11.2      2.234
    #    2.0  2.09138850722e-19 -2.01802656266e-15 ... 1123.38227429     0.0  3082.0       11.5        12.7      2.546
    #...
    #where 6, 5, 4, 3, 2, 1, 0 are the polynomial powers in p
    #    i.e. p where:
    #        p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
    :return: array of the parameters for the polynomial fit
    :return: minimum and maximum position of the trace in dispersion direction
    :return: array with the width of the trace in cross-dispersion direction, giving the width to the left, to the right, and the Gaussian width
    """
    # convert to astropy table
    if os.path.isfile(filename) == False:
        logger('Error: File {0} is missing.'.format(filename))
    atable = Table.read(filename)
    orders = np.array(atable['Order'], dtype=int)
    
    """if 'cen_poly0' not in atable.colnames:                                # Old way of storing the data, can be removed at a few months after July 2018
        pfits, xlows, xhighs, widths = [], [], [], np.array([])
        for order in orders:
            pfit = []
            for p in atable.colnames:
                if p not in ['Order', 'low_x', 'high_x', 'width_left', 'width_right', 'width_gaus', 'width_gauss']:  
                    pfit.append(atable[p][order])        
            pfits.append(pfit)
            xlows.append(int(atable['low_x'][order]))
            xhighs.append(int(atable['high_x'][order]))
            try:
                widths.append([atable['width_left'][order], atable['width_right'][order], atable['width_gauss'][order]])
            except:
                widths = []
        logger('Info: Master file to trace orders read: {0}'.format(filename))
        return np.array(pfits), np.array(xlows), np.array(xhighs), np.array(widths)"""
    
    # new way of storing the data
    pfits, widths = [], []
    for order in orders:
        cen, left, right = [], [], []
        for p in atable.colnames:
            if p.find('cen_poly') == 0 or p.find('cen_center') == 0:
                cen.append(atable[p][order])
            if p.find('left_poly') == 0 or p.find('left_center') == 0:
                left.append(atable[p][order])
            if p.find('right_poly') == 0 or p.find('right_center') == 0:
                right.append(atable[p][order])
        if len(left) == 0:                          # only 2d pfits
            pfits.append(cen)
        else:                                       # 3d pfits
            pfits.append([cen, left, right])
    xlows = atable['low_x']
    xhighs = atable['high_x']
    if 'width_left' in atable.colnames:
        widths = np.vstack([atable['width_left'], atable['width_right'], atable['width_gauss']]).T
        
    logger('Info: Master file to trace orders read: {0}'.format(filename))
    return np.array(pfits), np.array(xlows, dtype=int), np.array(xhighs, dtype=int), widths

def save_fits_width(pfits, xlows, xhighs, widths, filename):
    """
    save the fits to file, for use in extracting future spectra
    :param pfits: 2d array or 3d array of floats. If 2d: length same as number of orders, width same as order of polynomial + 1
                                                  If 3d: 3 entries of the 2d: for center, left, right
                      Each last dimension is the output of np.polyval, with a refit of the center as zerost entry
                      i.e. p where:
                      p[1]*(x-p[0])**(N-1) + p[2]*(x-p[0])**(N-2) + ... + p[N-2]*(x-p[0]) + p[N-1]
    :param xlows: list, length same as number of orders, the lowest x pixel
                   (wavelength direction) used in each order
    :param xhighs: list, length same as number of orders, the highest x pixel
                   (wavelength direction) used in each order
    :param widths: 2d list, length same as number of orders, each entry contains left border, right border, and Gaussian width of the lines, as estimated in the master flat
    :param filename: string, location and file name of the file containing the polynomial fits
    """
    if len(filename.rsplit(os.sep,1)) == 1:     # no path is in the filename
        logger('Warn: no folder to save {0} was given, using the current folder ({1}).'.format( filename, os.getcwd() ))
    elif not os.path.exists(filename.rsplit(os.sep,1)[0]):
        logger('Error: Folder to save {0} does not exists.'.format(filename))
        
    data = []
    temp = filename.rsplit('.',1)
    if temp[1] not in ['fit', 'fits']:
        logger('Warn: the filename does not end in fits: {0} . The proper ending will be added to save the file without problems, but this might cause issues when loading data. Please check your parameters.'.format(filename))
        filename += '.fits'
    
    # What is the maximum numbers of orders used for the polynomial fit
    pfitss = pfits.shape
    porder = pfitss[-1] - 1         # -1 in order to remove the central pixel
    
    # Loop around each order and plot fits
    for pp, pf in enumerate(pfits):
        row = [pp]
        if len(pfitss) == 3:        # 3d array
            for pfs in pf:
                row += list(pfs)
        else:                       # old 2d array
            row += list(pf)
        row += list([xlows[pp], xhighs[pp]])
        if len(widths) > 0:
            row += list(widths[pp])
        data.append(row)
    data = np.array(data)
    # make the header
    cols = ['Order']
    cols += ['cen_center'] + ['cen_poly'+str(i) for i in range(porder)[::-1]]
    if len(pfitss) == 3:        # 3d array
        cols += ['left_center'] + ['left_poly'+str(i) for i in range(porder)[::-1]]
        cols += ['right_center'] + ['right_poly'+str(i) for i in range(porder)[::-1]]
    cols += ['low_x', 'high_x']
    if len(widths) > 0:
        cols += ['width_left', 'width_right', 'width_gauss']

    # convert to astropy table
    atable = Table()
    for c, col in enumerate(cols):
        atable[str(col)] = data[:, c]
    atable.write(filename, overwrite=True)
    logger('Info: master file to trace orders written: {0}'.format(filename))

def save_wavelength_solution_to_fits(wave_sol_dict, filename):
    """
    save the wavelength solution in a file for later use
    :param wavelength_solution: 2d array of floats, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
                      (output of np.polyval)
                      i.e. p where:
                      p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
    :param wavelength_solution_arclines: 2d array of floats with one line for each order. Each order contains the wavelengths of the identified reference lines, 
                                        sorted by the brightest to faintest (in spectrum from which the solution is derived) and 0 to make into an array
    :param filename: string, location and file name of the file
    """
    if len(filename.rsplit(os.sep,1)) == 1:     # no path is in the filename
        logger('Warn: no folder to save {0} was given, using the current folder ({1}).'.format( filename, os.getcwd() ))
    elif not os.path.exists(filename.rsplit(os.sep,1)[0]):
        logger('Error: Folder to save {0} does not exists.'.format(filename))
    
    filename = filename.replace('.fits','').replace('.fit','')+'.fits'
    wavelength_solution = wave_sol_dict['wavesol']
    wavelength_solution_arclines = wave_sol_dict['reflines']
    line_stats = wave_sol_dict.get('linestat', None)
    
    data = np.hstack(( wavelength_solution, wavelength_solution_arclines ))
    if line_stats is not None:  data = np.hstack(( data, line_stats ))
    
    cols = ['real_order', 'central_px']
    cols += list(range(len(wavelength_solution[0]) -2))[::-1]
    cols += list([ 'arclin' + str(i) for i in range(len(wavelength_solution_arclines[0])) ])
    if line_stats is not None:  cols += ['width_med', 'width_std']
    
    # convert to astropy table
    atable = Table()
    for c, col in enumerate(cols):
        atable[str(col)] = data[:, c]
    atable.write(filename, overwrite=True)
    logger('Info: Wavelength solution written: {0}'.format(filename))

def read_wavelength_solution_from_fits(filename):
    """
    Reads the wavelength solution from a file
    :param filename: string, location and file name of the file
    fits file should look like the following:
    real_order central_px               2 ...             0 arclin1 ... arclin345 width_med width_std
       float64    float64         float64 ...       float64  loat64 ...   float64   float64   float64
    ---------- ---------- --------------- ... ------------- ------- ... --------- --------- ---------
        -101.0     2150.0 -7.82221887e-07 ... 5.6899272e+03 5678.91 ...  6543.987     1.234     0.456
        -100.0     2150.0 -7.39532842e-07 ... 5.7469279e+03 6542.10 ...   6677.32     1.534     0.257
    ...
    where n, ..., 3, 2, 1, 0 are the polynomial powers in p
        i.e. p where:
            p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
    :return wave_sol_dict: dictionary, which contains the wavelength_solution, wavelength_solution_arclines, and other information  
    :return wavelength_solution: 2d array of floats, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    :return wavelength_solution_arclines: 2d array of floats, same length as number of orders, each line contains the wavelength of the used reference lines. 
                                          Each order must have the same number of entries, if less reference lines are available in the order, the array is filled with 0.
    """
    if os.path.isfile(filename) == False:
        logger('Error: The file {0} with the wavelength solution does not exist'.format(filename)) 
    # convert to astropy table
    atable = Table.read(filename)
    wave_sol_dict = dict()    
    for p in atable.colnames:
        data = np.expand_dims(np.array(atable[p]), axis=1)        # expand_dims to make into a column vector, e.g. (9,1)
        if p.startswith('arclin'):      field = 'reflines'
        elif p.startswith('width_'):    field = 'linestat'
        else:                           field = 'wavesol'
        if field not in wave_sol_dict.keys():
            wave_sol_dict[field] = data
        else:
            wave_sol_dict[field] = np.hstack(( wave_sol_dict[field], data ))
    logger('Info: Master file for wavelength solution read: {0}'.format(filename))

    return wave_sol_dict

def save_multispec(data, fname, im_head, bitpix='-32'):
    """
    Saves a n-d array into a fits file
    :param data: n-d array of data or list of arrays
    :param fname: filename in which the data should be written. The filename will end in .fits
    :param im_head: header, which should be written to the fits file
    :param bitpix: Precission in which the data should be stored
    """
    if len(fname.rsplit(os.sep,1)) == 1:
        logger('Warn: No folder is given for the file {0}. File will be stored in the current working directory.'.format(fname))
    elif not os.path.exists(fname.rsplit(os.sep,1)[0]):
        logger('Error: Folder to save {0} does not exists.'.format(fname))
        
    fname = fname.replace('.fits', '')
    fname = fname.replace('.fit', '')
    if bitpix in ['-32', -32]:
        # relative difference is less than 1E-6 of the wavelength/flux value compared to float64, only needs half the size
        hdu = fits.PrimaryHDU(np.array(data).astype(np.float32), header=im_head)
    elif bitpix in ['-64', -64]:
        hdu = fits.PrimaryHDU(np.array(data).astype(np.float64), header=im_head)
    else:
        logger('Warn: The bitpix parameter ({0}) is unknown, using the data suggested one ({1})'.format(bitpix,hdu.header['BITPIX']))
        hdu = fits.PrimaryHDU(np.array(data), header=im_head)
    for index in range(hdu.header['NAXIS']):                # NAXISj doesn't has a Comment and Serval fails without it
        hdu.header['NAXIS{0}'.format(index+1)] = ( hdu.header['NAXIS{0}'.format(index+1)], 'length of data axis {0}'.format(index+1) )
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=AstropyUserWarning)
        hdu.writeto(fname+'.fits', overwrite=True)    # as .fits # clobber=True was used earlier
    
def rotate_flip_frame(im, params, invert=False):
    """
    Rotates and the flips the images in order to have traces along the right axis. Before saving, the images are processed in the opposite direction (invert=True) for consistency
    :param im: 2d numpy array with the image
    :param params: Dictionary with all the parameters. rotate_frame and flip_frame are required
    :param invert: If True,then the image will be first flipped and afterward rotated backwards
    :return im: processed image
    """
    #plot_img_spec.plot_image(im, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'orig')
    if invert == True:
        if params['flip_frame'] == True:
            im = np.fliplr(im)
        im = np.rot90(im, k=-int(params['rotate_frame']/90), axes=(0,1))
    else:
        im = np.rot90(im, k=int(params['rotate_frame']/90), axes=(0,1))
        if params['flip_frame'] == True:
            im = np.fliplr(im)
    #plot_img_spec.plot_image(im, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'result')

    return im
    
def bin_im(im, binxy, method='median'):
    """
    :param im: 2d numpy array
    :param binxy: list of integers (entries), contains the number of pixel which should be binned into one in x and y direction
    :return: 2d numpy arrays of the image, the number of elements which are not NaN, and of the standard deviation of each combined pixel
    """
    def bin_im_function(data, method, axis=None):
        if method == 'median':  return np.nanmedian(data, axis=axis)    # median binning
        elif method == 'mean':  return np.nanmean(data, axis=axis)    # averge binning
        elif method == 'sum':   return np.nansum(data, axis=axis)    # sum binning
        elif method == 'min':   return np.nanmin(data, axis=axis)    # min binning
        elif method == 'max':   return np.nanmax(data, axis=axis)    # max binning
        else:   logger('Error: method {0} does not exisit in bin_im_function. This is a programming bug.'.format(method))
    
    [binx, biny] = binxy
    ims = im.shape
    if binx <1 or biny <1 or (binx == 1 and biny == 1):
        logger('Warn: no binning possible: {0},{1}'.format(binx,biny))
        return im, im*0+1, im*0
    nim, gim, sim = [], [], []
    if binx > 1 and biny > 1:
        xx = np.arange(0, ims[0]+binx-1, binx)
        yy = np.arange(0, ims[1]+biny-1, biny)
        xx[-1] = ims[0]
        yy[-1] = ims[1]
        data = np.zeros(( xx.shape[0], yy.shape[0], 3 ))
        data.fill(np.nan)
        for ii in range(xx.shape[0]-1):
            for jj in range(yy.shape[0]-1):
                temdata = im[ xx[ii]:xx[ii+1], yy[jj]:yy[jj+1] ]
                nonan = np.sum( ~np.isnan(temdata) )
                if nonan >= 0.9 * np.prod(temdata.shape):    # only if at least 90% of data points available
                    data[ii,jj,0] = bin_im_function(temdata, method, axis=None)     # binning
                data[ii,jj,1] = nonan                        # number of the elements not nan
                if nonan > 1:
                    data[ii,jj,2] = np.nanstd(temdata, ddof=1)                      # standard deviation
                elif nonan == 1:
                    data[ii,jj,2] = 0
                
        """x_range = list(range(int((ims[0]+binx-1)/binx)))
        y_range = list(range(int((ims[1]+biny-1)/biny)))
        data = np.zeros(( len(x_range), len(y_range), 3 ))
        data.fill(np.nan)
        for ii in x_range:
            iline = im[ii*binx:min(ims[0],(ii+1)*binx),:]
            for jj in y_range:
                temdata = iline[:,jj*biny:min(ims[1],(jj+1)*biny)]
                if np.sum( ~np.isnan(temdata) ) >= 0.9 * np.prod(temdata.shape):    # only of at least 10% of data points available
                    data[ii,jj,0] = bin_im_function(temdata, method, axis=None)     # binning
                data[ii,jj,1] = np.sum( ~np.isnan(temdata) )                        # number of the elements not nan
                if data[ii,jj,1] > 1:
                    data[ii,jj,2] = np.nanstd(temdata, ddof=1)                      # standard deviation
                elif data[ii,jj,1] == 1:
                    data[ii,jj,2] = 0
                 
                #nline.append(sum(percentile_list(list((im[i*binx:(i+1)*binx,j*biny:(j+1)*biny]).flat),.2)))   #sum after 80% clipping -> more flux -> more lines find with gauss of height 50
        """
        nim = data[:,:,0]
        gim = data[:,:,1].astype(int)
        sim = data[:,:,2]
    elif binx > 1:
        for i in range(int((ims[0]+binx-1)/binx)):
            temdata = im[i*binx:min(ims[0],(i+1)*binx),:]
            nim.append(bin_im_function(temdata, method, axis=0))    # binning
            gim.append(np.sum( ~np.isnan(temdata), axis=0))  # number of the elements not nan
            sim.append(np.nanstd(temdata, ddof=min(1,temdata.shape[0]-1), axis=0))      # standard deviation, if only one column of pixels left, then std with ddof=0
    elif biny > 1:
        for i in range(int((ims[1]+biny-1)/biny)):
            temdata = im[:,i*biny:min(ims[1],(i+1)*biny)]
            nim.append(bin_im_function(temdata, method, axis=1))    # binning
            gim.append(np.sum( ~np.isnan(temdata), axis=1))  # number of the elements not nan
            sim.append(np.nanstd(temdata, ddof=min(1,temdata.shape[1]-1), axis=1))      # standard deviation
        nim = np.transpose(np.array(nim))
    #logger('Info: image binned')
    return np.array(nim), np.array(gim), np.array(sim)

def scale_image_plot(im, scale_type):
    """
    Applies a function to a 2d array, e.g. log
    :param im: 2d array, e.g. the image
    :param scale_type: Function, which should be applied to the {im}. So far the following functions are available:
                log10: decade logarithm
    :return im_plot: 2d image, to which the function was applied
    """
    im = np.array(im)
    if scale_type == 'log10':
        min_im = np.min(im)
        if min_im <= 0:
            add_value = 1-min_im
            im_plot = np.log10(im+add_value) - np.log10(add_value)
        elif min_im < 1:
            add_value = 1-min_im
            im_plot = np.log10(im+add_value)
        else:
            im_plot = np.log10(im)
    else:
        logger('Warn: The scaletype to scale the image is not defined: {}'.format(scale_type))
        im_plot = im
    
    return im_plot

def combine_median(ims):
    """
    Combines the 3d array of images by calculating the median (pixelwise).
    :param ims: 3d numpy array of the images which should be combined using the median
    :return im: 2d array of float, the median of the images (pixelwise)
    """
    im = np.median(ims, axis=0).astype(np.float64)
    #x,y = 1908,1572
    #x,y = 1699,1234
    #print ims[:,x,y],sorted(ims[:,x,y]),im[x,y]
    #exit(100)
    logger('Info: median combined {0} files'.format(ims.shape[0]))
    return im

def combine_sum(ims):           #no sigmacliping, as the sum needs to add always the same number of pixels
    """
    Combines the 3d array of images by calculating the sum (pixelwise). The maximum/minimum values are of each pixel are excluded
    :param ims: 3d numpy array of the images which should be added
    :return im: 2d array of float, the sum of the images (pixelwise)
    """
    if ims.shape[0] > 3:
        logger('Info: cleaning away maximum and minumum value before adding the stack')
        ims = np.sort(ims, axis=0)
        ims = ims[1:-1,:,:]
    im = np.sum(ims, axis=0).astype(np.float64)
    logger('Info: added {0} files together'.format(ims.shape[0]))
    return im

def percentile_list(data, prcentl):
    """
    Removes the highest and lowest prcentl percent of the data
    :param data: 1-d aray of data
    :param prcentl: integer or float, gives the number (in fractions, e.g. percent/100) of how much data to remove
    :return data: sorted list without the highest and smallest values, and without NaNs
    """
    if prcentl < 0 or prcentl > 0.5:
        print('Warn: Percentile must be in the range [0,0.5]')
    data = data[~np.isnan(data)]
    data = sorted(data)
    length = len(data)
    return data[int(round(prcentl*length)):int(round((1-prcentl)*length))]

def oneD_gauss(x, amplitude, x0=0, sigma=0, offset=0):
    """
    Calculates the Gauss for an array of data
    :param x: 1d aray of the x-axis
    :param amplitude: amplitude of the Gauss curve; or array of the parameters: [a,x0,sigma,b]
    :param x0: center of the Gauss curve
    :param sigma: width of the Gauss
    :param offset: offset of the Gauss in y
    :return: float/1d array, values of the Gauss
    """
    if x0==0 and sigma==0 and offset==0:     # Necessary in order for the fit to work: curve_fit will send individual parameters
        [amplitude, x0, sigma, offset] = amplitude
    if sigma <= 0:
        return np.repeat([np.NaN],len(x))
    return amplitude * np.exp(-(x-x0)**2/(2.*sigma**2))+offset

def oneD_blended_gauss(x, parameters, p01=0, p02=0, p03=0, p10=np.nan, p11=np.nan, p12=np.nan, p13=np.nan):
    """
    Calculates the Gauss for an array of data
    :param x: 1d aray of the x-axis
    :param parameters: 2d list or array of floats, first dimension is the number of Gauss to fit, second dimension is the 4: the parameters for a*exp(-(x-x0)**2/(2*sigma**2))+b:
                                    a: amplitude of the Gauss curve
                                    x0: center of the Gauss curve
                                    sigma: width of the Gauss
                                    b: offset of the Gauss in y
    :return: float/1d array, values of the Gauss
    """
    if type(parameters) is not np.ndarray and type(parameters) is not list:     # Necessary in order for the fit to work: curve_fit will send individual parameters
        parameters = [[ parameters, p01, p02, p03 ]]
        if not np.isnan(p10):
            parameters.append([ p10, p11, p12, p13 ])
    parameters = np.array(parameters)
    parameters.shape = ( int(np.prod(parameters.shape)/4), 4 )                       # curve_fit will create a 1d array
    if np.min(parameters[:,2]) <= 0:       #check that the parameters give useful results
        return np.repeat([np.NaN],len(x))
    result = np.zeros(len(x))
    for i in range(len(parameters)):
        result += oneD_gauss(x, parameters[i])
    return result

def twoD_Gaussian(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):          # Changed from "x, y, amplitude,..." to  "(x, y), amplitude,..." on 20190528 after curve fit failed. Does it depend on the scipy version for curve fit? (x,y) doesn't work with python3
    """
    Calculates the Gauss in 2 dimensions
    :param (x,y)/xy: lists of the x- and y- values for which the Gauss should be calculated
    :param amplitude: Amplitude of the Gauss function
    :param x0: center of the Gauss in x direction
    :param y0: center of the Gauss in y direction
    :param sigma_x: Gaussian width on x direction
    :param sigma_y: Gaussian width on y direction
    :param theta: rotation of the Gauss compared to the x-axis
    :param offset: zero level of the Gauss
    """
    [x,y] = xy      # for python 3
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2.*sigma_x**2) + (np.sin(theta)**2)/(2.*sigma_y**2)
    b = -(np.sin(2*theta))/(4.*sigma_x**2) + (np.sin(2*theta))/(4.*sigma_y**2)
    c = (np.sin(theta)**2)/(2.*sigma_x**2) + (np.cos(theta)**2)/(2.*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

def polyfit_adjust_order(xarr, yarr, p_orders, w=None):
    """
    Finds the polynomial with the highest posible order to fit the data
    """
    xarr = np.array(xarr)
    yarr = np.array(yarr)
    poly = np.array([np.mean(yarr)])            # if the order = 0 is failing, then report at least the average (which is a polynom of order 0)
    if yarr.shape[0] > 1:
        if np.std(yarr) < 1E-10:
            return poly
    for order in range(p_orders+1)[::-1]:       # Try smaller number of polynomial orders, if there is a RankWarning
            with warnings.catch_warnings():
                warnings.filterwarnings('error')
                try:
                    poly = np.polyfit(xarr, yarr, order, w=w)
                    break
                except np.RankWarning:
                    poly = poly
    return poly

def polynomial_value_2d(xx, yy, xord, yord, poly2d_params):
    """
    Calculates the Dot product between the 2d-coordinates and the parameters of the polynomial fit, e.g. gives the result of fit
    :param xx: array of floats, flattened results of a meshgrid of the x-coordinates
    :param yy: array of floats, same dimension of xx, flattened results of a meshgrid of the y-coordinates
    :params xord, yord: number of orders for the polynomial to be used in each direction, order = 0 means constant value, same conversion as np.polyfit
    :param poly2d_params: Result of fit of a 2d polynomial to the data
    :return: 
    """
    a = polynomial_array_2d(xx, yy, xord, yord)
    return np.dot(a,poly2d_params)

def polynomial_array_2d(xx, yy, xord, yord):
    """
    Calculates an array with the possible combination for xx**i times yy**i
    :param xx: array of floats, flattened results of a meshgrid of the x-coordinates (or just a 1d array, if not datapoints in an 2d-array)
    :param yy: array of floats, same dimension of xx, flattened results of a meshgrid of the y-coordinates
    :params xord, yord: number of orders for the polynomial to be used in each direction, order = 0 means constant value, same conversion as np.polyfit
    :return: Arry with the dimension of xx or yy in one axis and number possible combinations in the other direction (i.e. 16 for xord=3 and yord=3, or 25 for 4,4)
    """
    
    ims = np.insert(xx.shape, 0, (xord+1)*(yord+1)).astype(int)
    a = np.ones(ims)                                    # first entry is 1 -> constant term
    index = 1
    for k in range(1, (xord+1)*(yord+1) ):             # All combinations, starting with the lowest order ones in order
        for j in range(yord+1):
            for i in range(xord+1):
                if i+j == k:
                    data = xx**i * yy**j
                    #if k==3:
                    #    print 'xx**{0} * yy**{1}'.format(i,j), data
                    if len(ims) == 1:
                        a[index] = data
                    else:
                        a[index,:] = data
                    index += 1
    #if a.nbytes > 0.1E9:# and a.nbytes > array_polyfit_2d.nbytes:
    #    array_polyfit_2d = copy.deepcopy(a).T
    return a.T

def polynomial_fit_2d_norm(xx, yy, zz, xord, yord, w=[], w_range=1E50):
    """
    Fits a 2d poynormial to the data.
      Before the fit the axes are normalised in order to weight both axis the same way
      The weights are 
    :params xx,yy,zz: all data needs to have the same entries: x-value, y-value, and z-value
    :params xord, yord: number of orders for the polynomial to be used in each direction, order = 0 means constant value, same conversion as np.polyfit (was order=1 before 20171207)
    :params w: weights of the individual points. Weights are created by copying the value the respective times, depending on the weight
    :params w_range: break down the weights to that range. The smaller the faster the calculations, but also the less the weights are taken into account
    :returns coeff: 1d array with the coefficients of the 2d fit
    """
    ww = w
    good_values = ( ~np.isnan(xx) & ~np.isnan(yy) & ~np.isnan(zz) )                    # remove the nan values
    if np.sum(good_values) < len(good_values):          # log, if there were nan-values
        logger('Warn: NaN-values in the input values for the 2d-fit. This probably means that something went wrong before. '+\
               'Of the available {0} values, {1} in x, {2} in y, and {3} entries in z-direction are NaN'.format(len(good_values), 
                    np.sum(np.isnan(xx),axis=None), np.sum(np.isnan(yy),axis=None), np.sum(np.isnan(zz),axis=None)))
    xx = xx[good_values]                    # remove the nan values
    yy = yy[good_values]                    # remove the nan values
    zz = zz[good_values]                    # remove the nan values
    # Get rid of the data points with weight == 0
    if len(ww) > 0:
        ww = ww[good_values]                    # remove the nan values
        if np.min(ww, axis=None) < 0:
            logger('Warn: The weights for the 2d polynomial fit contains values below 0. This is probably a programming mistake.')
        good_values = (ww > 1E-10)
        xx = xx[good_values]
        yy = yy[good_values]
        zz = zz[good_values]
        ww = ww[good_values]
    # Normalise the data so that the maximum is -1 or +1, make sure that a proper float division is done
    x_normf = np.max(np.abs(xx))+0.0
    y_normf = np.max(np.abs(yy))+0.0
    A = polynomial_array_2d(xx, yy, xord, yord)
    #print(sys.getsizeof(A), A.nbytes, A.shape)
    norm_f = polynomial_array_2d(x_normf, y_normf, xord, yord).T            # get the normalisiation for each possible combination for xx**i times yy**i
    for i,norm in enumerate(norm_f):     # Normalise the array
        A[:,i] = A[:,i]/(norm+0.0)
    # Add the weights -> more entries if the weight is higher (crude solution)
    if len(ww) > 0:
        if min(ww)/(max(ww)+0.0) < 0.5:
            # Normalise so that the smallest weight has a value of 1
            #ww = ((ww/min(ww))**2 * 10).astype(int)       # Normalise so that the smallest weight has a value of 10, but this increaese the array too massively
            Azz = np.hstack((A, np.expand_dims(zz, axis=1) ))
            ww = ww/(min(ww)+0.0)
            if max(ww) > w_range:
                ww = ww / ( max(ww) / (w_range+0.0) )             # scale by w_range
            mem = psutil.virtual_memory()                   # svmem(total=33221091328, available=28485840896, percent=14.3, used=4202041344, free=25513508864, active=..., inactive=..., buffers=..., cached=.., shared=...)
            if np.sum(ww) * (A.shape[1] + 1) * 8 > mem[1] * 0.4:        # if more than 40% of the available memory needs to be used for Azz_n (subprocedure needs memory as well), 8 bytes per entry
                ww = (ww / ( np.sum(ww) * (A.shape[1] + 1.) * 8. / (mem[1] * 0.4) )).astype(int)      # ignore the lowest weight entries
            else:
                ww = ww.astype(int)
            #print "ww = 0 / > 0 :", len(ww[(ww == 0)]), len(ww[(ww != 0)])
            ww[ww == 0] = 1                             # even the lowest weights should apear at least once
            wws = np.sum(ww)
            #if wws > 1E6:          # so quick now, that not necessary anymore
            #    logger('Info: To use the correct weights for the fit of a 2d polynomial an array with size ({0}, {1}) will be created.'.format(wws, Azz.shape[1]))
            Azz_n = np.empty((wws, Azz.shape[1]))
            start = 0
            for i in range(ww.shape[0])[::-1]:
                Azz_n[start:start+ww[i], :] = Azz[[i],:].repeat(ww[i],axis=0)
                start += ww[i]
                #if ww[i] > 1:  # old and slow
                #    A = np.concatenate( (A, A[[i],:].repeat(ww[i]-1,axis=0)) , axis=0)
                #    zz = np.concatenate( (zz, zz[[i]].repeat(ww[i]-1,axis=0)) , axis=0)
            if start != Azz_n.shape[0]:
                logger('Warn: A difference occured in procedure polynomial_fit_2d_norm, which should nor occur. Please inform the Author. Values: {0}, {1}'.format(start, Azz_n.shape[0]))
            A = Azz_n[:,:-1]
            zz = Azz_n[:,-1]
            del Azz
    if np.prod(A.shape) > 1E9:
        logger('Step: Performing the fit, that might take a bit')
    coeff, res, rank, s = np.linalg.lstsq(A, zz, rcond=1E-50)        # for wavelengths fit
    #fit = polynomial_value_2d(xx, yy, xord, yord, coeff)             # calculates the fit
    #print np.array([xx,yy,zz, fit]).T                 to see the residuals between zz and the fit
    
    return coeff/(norm_f+0.0)

def fit_2d_image(im, xord, yord, w=[], binning=[]):
    """
    Returns the fitted image
    :params im: 2d numpy array, i.e. the image
    :params xord, yord: number of orders for the polynomial to be used in each direction, yord is along dispersion axis
    :params w: weights of the individual points.
    :params binning: List of two integers or empty list: Binning before the fit increases the calculation speed.
                    Binning should be done only in x-axis (dispersion), because binning in one direction is much faster and to avoid lossing areas to fit.
    """
    """print(time.time(),'Begin ori')
    ims = im.shape
    cen = [ ims[0]/2. , ims[1]/2. ]
    x, y = np.meshgrid(np.arange(ims[1])-cen[1], np.arange(ims[0])-cen[0] )       # reverse order, otherwise the flatten doesn't create the right values for z
    xx, yy, zz = x.flatten(), y.flatten(), im.flatten()
    if len(w) > 0:
        w1 = w.flatten()

    poly2d_params = polynomial_fit_2d_norm(xx, yy, zz, xord, yord, w=w1)       #x,y,z,order_x,order_y, w=weights
    im_fito = polynomial_value_2d(xx, yy, xord, yord, poly2d_params)
    #print np.array([xx,yy,zz, im_fit]).T
    im_fito.shape = ims
    imo = copy.copy(im)
    wo = copy.copy(w)
    print(time.time(),'End old')"""
    ims = im.shape
    if len(binning) > 0:         # Fit on a binned image
        binx = binning[0]
        biny = binning[1]            # don't bin in cross-dispersion direction (checked in image)
        im, dummy,dummy = bin_im(im, [binx,biny], method='median')
        w,  dummy,dummy = bin_im(w,  [binx,biny], method='min')
    cen = [ ims[0]/2. , ims[1]/2. ]
    x, y = np.meshgrid(np.arange(0, ims[1], biny)-cen[1], np.arange(0,ims[0],binx)-cen[0] )       # reverse order, otherwise the flatten doesn't create the right values for z
    xx, yy, zz = x.flatten(), y.flatten(), im.flatten()
    if len(w) > 0:
        w = w.flatten()

    poly2d_params = polynomial_fit_2d_norm(xx, yy, zz, xord, yord, w=w)       #x,y,z,order_x,order_y, w=weights
    # Full coordinates
    x, y = np.meshgrid(np.arange(ims[1])-cen[1], np.arange(ims[0])-cen[0] )       # reverse order, otherwise the flatten doesn't create the right values for z
    xx, yy = x.flatten(), y.flatten()
    im_fit = polynomial_value_2d(xx, yy, xord, yord, poly2d_params)
    #print np.array([xx,yy,zz, im_fit]).T
    im_fit.shape = ims

    #plot_img_spec.plot_image(im, [''], 1, True, [0.05,0.95,0.95,0.05], 'orig')
    #plot_img_spec.plot_image(im_fit, [''], 1, True, [0.05,0.95,0.95,0.05], 'fit')
    #plot_img_spec.plot_image((imo-im_fit)*wo, [''], 1, True, [0.05,0.95,0.95,0.05], 'residuals')
    
    return im_fit

def sigmaclip(data, x=[], nc=1, sub=[], ll=3., lu=3., repeats=5, exclude_absorption=False, x_cen=np.nan):
    """
    It ignores absorption lines
    :param data: 1d array with the data
    :param x: 1d array with the x-values of the data, can be empty
    :param nc: integer, orders of the polynomial
    :param sub: List or array of bool (or ones and zeros), True for entries that should be used to calculated the polynomial and the the standard deviation
    :param ll: float or integer, Data off by this sigma are rejected on the lower side of the fit
    :param lu: float or integer, Data off by this sigma are rejected on the higher side of the fit
    :param repeats: redo the fit with the (cleaned) data how many times? 
    :param exclude_absorption: bool, if true than absorption lines are ignored and noise (std) is calculaded mirroring the difference values around 0
    :param x_cen: float, pressumed centre of the values in x
    :return goodvalues: 1d array with True/False. The values which are inside the limits are True
    :return poly: parameters of the last polynomial fit
    """
    poly = np.repeat([0], nc)
    if len(data) == 0:
        logger('Warn: empty array for sigma clipping.')
        return np.array([]), poly
    if len(np.array(data).shape) > 1:
        logger('Warn: Data array for sigma clipping is not 1-dimensional ({0}). This is a programming error and will most likely lead to a crash.'.format(np.array(data).shape))
    if len(x) != len(data):
        if len(x) > 0:
            logger('Warn: Got different sized arrays for sigma clipping (y: {0}, x: {1}). This is a programming error and should not happen.'.format(len(data), len(x)))
        x = np.arange(len(data))                        # Create a numbered list
    if len(data) > 1:
        if np.std(data) < 1E-10:
            return np.ones(len(data)), np.array([0]*nc+[data[0]])
    if np.isnan(x_cen):
        x_cen = (np.nanmin(x) + np.nanmax(x))/2.        # central wavelength
    if len(sub) != len(data):
        if len(sub) > 0:
            logger('Warn: Got different sized arrays for the good values (y: {0}, sub: {1}). This is a programming error and should not happen.'.format(len(data), len(sub)))
        sub = np.ones(len(data),dtype=bool)             # Use all data
    good_data = ~np.isnan(data)
    good_values = np.ones(len(data),dtype=bool)         # Use all data
    for ii in range(repeats):
        if np.sum(good_values*good_data*sub) <= nc:
            break
        poly = polyfit_adjust_order(x[good_values*good_data*sub]-x_cen, data[good_values*good_data*sub], nc )        # Fit the data
        res = data - np.polyval(poly, x-x_cen)              # Residuals for the whole array, >0 for emission lines
        if exclude_absorption:
            res[np.isnan(res)] = -1E6                       # To avoid res >= 0 giving invalid value warning, they are excluded in dev and will give False in good_values
            IU = (res >= 0)                                 # Positive residuals, use only non-absorption lines
            subarr = res[good_values*sub*IU]
            dev = np.nanstd((subarr,-subarr), ddof=nc+1)    # (subarr,-subarr) as all negative values are ignored -> stdev will be 2 low. Only correct if scatter is similar to both sides
        else:
            dev = np.nanstd(res[good_values*sub], ddof=nc+1)
        #print ii, dev, sum(IU), sum(sub), sum(good_values), sum(good_values*sub*IU)
        good_values = ( (res < lu*dev) & (res > -ll*dev) )  # Sigmacliping, but allow already clipped data back in
    
    return good_values, poly

def width_order(y, oldcenter, maxFWHM):
    """
    Determines the width and center of a line by checking when the line merges into the next line or reaches the basic level
    :param y: 1d array
    :param oldcenter: initial position of the center (e.g. maximum)
    :return: integers, half-width of the order, centerposition of the order, left border of line, right border of the line
    """
    smoothdiff = 2      # 1 means neighboring pixels
    # find minimum to either side of the oldcenter
    leftmin, dummy = min(enumerate(y[max(0,int(round(oldcenter-maxFWHM))):max(1,int(round(oldcenter)))]), key=operator.itemgetter(1))
    rightmin, dummy = min(enumerate(y[min(len(y)-2,int(round(oldcenter+1))):min(len(y)-1,int(round(oldcenter+maxFWHM+1)))]), key=operator.itemgetter(1))
    leftmin += max(0,int(round(oldcenter-maxFWHM)))            # first minimum left of oldcenter
    rightmin += min(len(y)-2,int(round(oldcenter+1)))               # first minimum right of oldcenter
    linemax = max(y[max(0,int(round(oldcenter))-1):min(len(y),int(round(oldcenter))+2)])
    linemin = min(y[leftmin:rightmin+1])
    
    #find half width 30% width in this range
    range30pct = np.where(y[leftmin:rightmin+1] >= 0.3*(linemax-linemin)+linemin)[0] + leftmin
    
    # find end of the line: below 2% threashold
    range2pct = np.where(y[leftmin:rightmin+1] >= 0.02*(linemax-linemin)+linemin)[0] + leftmin
    
    # find end of the line: changes in flux are small
    diff = np.abs(y[0:len(y)-1] - y[1:len(y)])
    diffsmooth = np.zeros(len(diff))
    for i in range(smoothdiff):
        diffsmooth = diffsmooth[0:len(diff)]+diff[0:len(diff)]
        #print diffsmooth[:5].astype(int),diffsmooth[-5:].astype(int)
        diff = diff[1:]
    rangediff = np.where(diffsmooth[leftmin:rightmin+1] >= 0.01*linemax)[0] + leftmin
    left,right = oldcenter,oldcenter
    if len(rangediff) != 0 and len(range2pct) != 0:
        left = max(min(rangediff),min(range2pct))
        right = min(max(rangediff)+smoothdiff,max(range2pct))
    elif len(rangediff) != 0:
        left = min(rangediff)
        right = max(rangediff)+smoothdiff
    elif len(range2pct) != 0:
        left = min(range2pct)
        right = max(range2pct)
    #if right - left < 2:
    #    print 'leftmin,rightmin,linemax,linemin,range30pct,range2pct,rangediff,y[leftmin:rightmin+1].astype(int),left,right',leftmin,rightmin,linemax,linemin,range30pct,range2pct,rangediff,y[leftmin:rightmin+1].astype(int),left,right
    return int(round((right-left+.1)/2)), int(round((right+left+.1)/2)), leftmin, rightmin

def centroid_order(x, y, center, width, significance=3, bordersub_fine=True, blended_gauss=False, bugfix=False):
    """
    Fits the center of a line
    :param x: 1d array of px
    :param y: 1d array of flux
    :param center: expected center of the line
    :param width: width of the line
    :return popt: list of 4 floats, parameter of the fit of the Gauss: amplitude of the Gauss curve, center of the Gauss, width of the Gauss, offset of the Gauss in y
    """
    x, y = np.array(x), np.array(y)
    width = abs(width)
    notnan = ~np.isnan(y)
    x = x[notnan]
    y = y[notnan]
    if len(y) < 7:
        return [0,0,0,0]
    rangey= max(y)-min(y)
    pos_max = y.argmax()
    if abs(x[pos_max] - center) > width/2.:              # changed on 20180525 from (x[pos_max] - center) < width/2
        center = x[pos_max]
    p0=[rangey,center,width,min(y)]
    bounds=((0.2*rangey, center-max(1.2,0.3*width), 0.1*width, min(y)-0.2*rangey), (5.*rangey, center+max(1.2,0.3*width), 2.*width, min(y)+0.2*rangey))
    if blended_gauss:
        p0 = p0 + p0
        bounds = (bounds[0] + bounds[0], bounds[1] + bounds[1])
    stdlin = np.std(y, ddof=1)
    if bordersub_fine:
        border_subs = [ [0,0], [0,1], [1,0], [1,1], [2,1], [1,2], [2,2] ]   # number of pixels to be removed from the begin/end of the area for the Gaussian fit
    else:
        border_subs = [ [0,0], [2,2] ]              # Testing all the ones above can take too long
    significant = False
    popt_test, signif_test = [], []
    for border_sub in border_subs:    
        range_data = list(range(border_sub[0],len(x)-border_sub[1]))
        if len(range_data) <= 6:
            break
        if not blended_gauss:
            #get_timing('{0} '.format(border_sub))
            try:
                popt,pcov = curve_fit(oneD_gauss,x[range_data],y[range_data],p0=p0, bounds=bounds)            #a, x0, sigma, b: a*np.exp(-(x-x0)**2/(2*sigma**2))+b
            except:
                # print 'curve fit failed'
                continue
            #get_timing('after ')
            stdfit = np.std(oneD_gauss(x[range_data],popt)-y[range_data], ddof=len(popt))           # if the average is 0, than this gives the stray of the residuals
            if stdlin/(stdfit+0.0) >= significance or popt[0]/(stdfit+0.0) >= significance:
                significant = True
                popt_test.append(popt)
                signif_test.append(stdfit)
                if len(signif_test) >= 3 and np.std(signif_test, ddof=1) < 0.1*np.median(signif_test):  # Several non-scattering data points
                    break
                #plot_img_spec.plot_spectra(np.array([x,x]),np.array([y,oneD_gauss(x,popt)]),['data','fit'], ['savepaths'], True, [0.06,0.92,0.95,0.06, 1.0,1.01], 
                #                           'Significant fit, stdlin={0}, stdGauss={1}, height={2}, needed significance={3}'.format(stdlin,stdfit,popt[0],significance))
            elif bugfix:
                print('Gauss not significant', p0,bounds, stdfit, border_sub)
                plot_img_spec.plot_spectra(np.array([x,x]),np.array([y,oneD_gauss(x,popt)]),['data','fit'], ['savepaths'], True, [0.06,0.92,0.95,0.06, 1.0,1.01], 
                                           'No significant fit, stdlin={0}, stdGauss={1}, height={2}, needed significance={3}'.format(stdlin,stdfit,popt[0],significance))
        else:
            try:
                popts,pcov = curve_fit(oneD_blended_gauss,x[range_data],y[range_data],p0=p0, bounds=bounds)            #a, x0, sigma, b: a*np.exp(-(x-x0)**2/(2*sigma**2))+b
            except:
                # print 'curve fit failed'
                continue
            stdfit = np.std( oneD_blended_gauss( x[range_data], popts) - y[range_data], ddof=len(popts) )
            if stdlin/(stdfit+0.0) >= significance or popts[0]/(stdfit+0.0) >= significance:
                # not real: sigma1/sigma2 > 2 or sigma2/sigma1 > 2 -> 1d fit
                # if x0s are within one px (diff within 10% of max) and sigmas more or less the same: -> 1d fit
                # check that each gauss has a significant height, otherwise -> 1d
                # allow more shift ub the bounds, at the moment double lines a few pixels apart are not resolved
                plot_img_spec.plot_spectra(np.array([x,x,x,x]),np.array([y,oneD_blended_gauss( x, popts), oneD_gauss( x, popts[:4]), oneD_gauss( x, popts[4:])]),
                                           ['data','fit','fit1','fit2'], ['savepaths'], True, [0.1,0.9,0.9,0.1, 1.0,1.01], 
                                           'Significant fit, stdlin={0}, stdGauss={1}, popts={2}'.format(stdlin,stdfit,popts))
                significant = True
                popts.shape = ( (int(len(popts)/4), 4) )
                for popt in popts:
                    print(popt)
                    oneD_gauss( x[range_data], popt )
                popt_test.append(popts)
                signif_test.append(stdfit)
            elif bugfix:
                print('Gauss not significant', p0,bounds, stdfit, border_sub)
                plot_img_spec.plot_spectra(np.array([x,x]),np.array([y,oneD_blended_gauss( x, popts)]),['data','fit'], ['savepaths'], True, [0.06,0.92,0.95,0.06, 1.0,1.01], 
                                           'No significant fit, stdlin={0}, stdGauss={1}, height={2}, needed significance={3}'.format(stdlin,stdfit,popt[0],significance))
            
    #print 'stdfit,stdlin', stdfit,stdlin,stdlin/stdfit < significance, popt[0]/stdfit < significance, popt
    if not significant:   # fit is not significant
        return  np.array([0,0,0,0])
    else:
        best_pos = np.argmin(signif_test)
        popt = popt_test[best_pos]
    return popt

def find_border_pctl(data, border_pctl=50, full_px=False):
    """
    Calculates the position where the signal falls below a certain threashold. If enough data points are available, then a linear fit is done in the transition area to get sub-pixel precision
    :param data: 1d list or array of floats, contains a Gauss-like signal
    :param border_pctl: integer or float, defines the percentage of the signal at which the cut should be done
    :param full_px: bool, if only full pixels are required the calculation can be sped up
    :return result[0]: float (int if full_px), position on the left of the maximum where the signal in data falls below border_pctl percent
                                              (or the full pixel before the flux goes above border_pctl)
    :return result[1]: float (int if full_px), position on the right of the maximum where the signal in data falls below
                                               (or the full pixel before the flux goes above border_pctl)
    """
    if len(data) < 2:
        return 0,len(data)
    data = copy.copy(data)                                              # Otherwise will be modified also in the routine calling find_border_pctl
    data -= min(data)                                                   # Set min(data) to zero
    ran_pctl = max(data)*border_pctl/100.
    higher = np.where( data > ran_pctl)[0]                              # Find the indexes above the percentile value, can have gaps
    higher = [ max(1, min(higher)) , min(len(data)-2, max(higher)) ]    # To allow a 2px fit, exclude the outermost pixels; otherwise it contains the first and last index over border_pctl
    lower = np.where( data < ran_pctl)[0]
    lower1 = lower[ lower < higher[0] ]                                 # The left side lower than the percentile value
    if len(lower1) < 1:                                                 # happens if the line is cut before it reaches the lower level
        lower1 = [0]
    lower2 = lower[ lower > higher[1] ]                                 # The righ side lower than the percentile value
    if len(lower2) < 1:                                                 # happens if the line is cut before it reaches the lower level
        lower2 = [len(data)-1]
    if full_px:
        return max(lower1), min(lower2)                                 # Returns full pixel to speed things up
    ran1 = list(range( max(lower1), higher[0]+1 ))                      # This should be a 2px (or more) transition around the transition
    ran2 = list(range( higher[1], min(lower2)+1 ))                      # This should be a 2px (or more) transition around the transition
    result = []
    for ran in [ran1, ran2]:
        if np.std(data[ran], ddof=1) == 0.0:                            # same flux values in all pixel -> fit is not possible
            result.append(np.median(ran))
        else:
            pfit = np.polyfit( data[ran], ran, 1)                       # linear fit px over data (corrected by minimum)
            fit_pos = np.polyval(pfit, ran_pctl)
            if min(ran)-1 <= fit_pos <= max(ran)+1:
                result.append( fit_pos )
            else:
                result.append( np.nan )
    #if np.isnan(result[0]).any():
    #    print result, ran1, ran2, data[ran1], data[ran2], 'ran_pctl', ran_pctl, border_pctl, len(data), data
    #else:
    #    print result
    return result[0], result[1]
    
def find_center(imslice, oldcenter, x, maxFWHM, border_pctl=0, border_pctl_poly=75, significance=3.0, polymax=False, bugfix=False, dummy=[]):
    """
    Determines the center of a gauss
    :param imslice: 1d array to search for the center
    :param oldcenter: integer, initial guess of the center (e.g. maximum)
    :param x: integer, index of the 1d array in the overal image, only used for error investigation
    :param maxFWHM: integer, estimated width of the order
    :param significance: float, Significance of the Gaussian fit compared to a simple average
    :param polymax: bool, fits the centre also with a 3rd order polynom to find the maximum in this
    :return centerfit: float, center from the centroid fit
    :return width: float, half width of the order at the base, set to 0 in case the fit wasn't good
    :return leftmin: integer, left border of the line
    :return rightmin: integer, right border of the line
    """
    imslice = copy.copy(imslice)        # Otherwise the data in the original image might be modified
    maxFWHM = max(3, maxFWHM)           # To have enough pixel (in order to fit the Gauss)
    oldcenter = max(0, oldcenter)       # The oldcenter needs to be in the data
    oldcenter = min(len(imslice)-1,oldcenter)
    width, center, leftmin, rightmin = width_order(imslice,oldcenter,maxFWHM)
    cen_poly = np.nan
    if width > maxFWHM or width == 0:                                        
        width = 0
        centerfit = center
    else:
        while rightmin-leftmin < max(5,2*maxFWHM):             # To have enough pixel to fit the Gauss, 2 times the maximum FWHM
            leftmin = max(0, leftmin-1)
            rightmin = min(len(imslice)-1, rightmin+1)
        leftminori, rightminori = copy.copy(leftmin), copy.copy(rightmin)
        sub_slice = imslice[leftmin:rightmin+1]
        if len(range(leftmin,rightmin+1)) != len(sub_slice):
            print('Error: not same length', x, oldcenter, len(imslice), leftmin, rightmin)
        for border in [0,1]:                # Try to fit different widths 
            popt = centroid_order(list(range(leftmin+border,rightmin+1-border)), sub_slice[border:len(sub_slice)-border], center, width, significance=significance, bugfix=bugfix)
            if popt[2] > 0:                 # run the second fit only when the first fails
                break 
        centerfit = popt[1]
        width = popt[2]
        
        # Fit the centre with a polynomial as well, using 75% of the flux
        if centerfit > 0 and width > 0 and polymax:
            cen_poly, poly,l,r = fit_poly_center(sub_slice, border_pctl=border_pctl_poly)
            cen_poly += leftmin

        # Redefine the borders of the order
        if border_pctl > 0:
            left, right = find_border_pctl(sub_slice, border_pctl=border_pctl)          # left, right in terms of sub_slices pixels
            rightmin = {True:rightmin, False:right+leftmin}[np.isnan(right)]            # leftmin as this is the zero-point in sub_slice
            leftmin  = {True:leftmin,  False:left+leftmin}[np.isnan(left)]
    
    """datax, datay, label = [], [], []
    datax.append(np.arange(leftminori,rightminori+0.1,0.1))
    datay.append(oneD_gauss(datax[-1],dummy) )
    label.append('old gauss fit')
    datax.append(np.arange(leftminori,rightminori+0.1,0.1))
    datay.append(oneD_gauss(datax[-1],popt) )
    label.append('gauss fit')
    datax.append(np.arange(l+leftminori, r+leftminori+.1, 0.1))
    datay.append(np.polyval(poly, datax[-1]-leftminori) )
    label.append('poly fit')
    datax.append(list(range(leftminori,rightminori+1)))
    datay.append(sub_slice)
    label.append('data')
    plot_img_spec.plot_points(datax, datay,label,'path',show=True, x_title='Pixel', y_title='Flux', title='' )"""

    return centerfit, cen_poly, width, leftmin,rightmin

def fit_poly_center(imslice, x=[], border_pctl=75):
    """
    Fits a third order polynomial against a non-gaussian curve to find the maximum
    :param imslice: list or array of floats, non-gaussian line 
    :param x: list or array of floats, x values to imslice
    :return cen_poly: float of the maximum of the polynomial fitted to imslice
    :return left, right: borders of the area used to fit the polynomial
    """
    nc = 3      # Third order polynomial to fit an asymetric parabola
    cen_poly, poly = np.nan, np.array([np.nan]*4)
    if len(x) != len(imslice):
        x = np.arange(len(imslice))
    data = np.vstack((x,imslice))
    data = data[:, data[0,:].argsort() ]            # sort by x
    left, right = find_border_pctl(data[1,:], border_pctl=border_pctl, full_px=True)     # The pixels before and after the flux is above 75% of max
    if right - left >= 3:
        x = np.arange(data.shape[1])
        mask = ( ( x >= left ) & ( x <= right ) )
        poly = np.polyfit(data[0,mask], data[1,mask], max(2, min(np.sum(mask)-1, nc)) )
        left, right = np.nanmin(data[0,mask]), np.nanmax(data[0,mask])
        crit = [ x.real for x in np.poly1d(poly).deriv().r if x.imag == 0 and left < x.real < right ]   # Critical points of the curve
        if len(crit) > 0:
            y_crit = np.polyval(poly, crit)
            cen_poly = crit[y_crit.argmax()]
    
    """ with fraction of pixel
    if not np.isnan(left+right):
        x = np.arange(data.shape[1])
        left, right = int(left), int(np.ceil(right))
        mask = ( ( x >= left ) & ( x <= right ) )
        if np.sum(mask) >= 3:
            poly = np.polyfit(data[0,mask], data[1,mask], max(2, min(np.sum(mask)-1, nc)) )
            left, right = np.nanmin(data[0,mask]), np.nanmax(data[0,mask])
            crit = [ x.real for x in np.poly1d(poly).deriv().r if x.imag == 0 and left < x.real < right ]   # Critical points of the curve
            if len(crit) > 0:
                y_crit = np.polyval(poly, crit)
                cen_poly = crit[y_crit.argmax()]"""
                
    return cen_poly, poly, left, right

def measure_noise(spectra, p_order=16, semi_window=10):
    """
    Measures the noise in a spectrum by fitting a high order polynomial against the data and deriving the residuals
    !!! Possible improvement: Take into account the wavelength solution -> use a semi_window that only few lines are covered
    :param spectra: 2d array of floats, spectrum
    :param p_order: integer, order of the polynomial
    :param semi_window: integer, the noise is calculated using a sample of +/- semi_window pixels around the current value; tests were done with semi_window=10
    :return noise: 2d array of floats, same dimension as spectra
    """
    while 2*semi_window+1 - p_order < 6:           # make sure that the values make sense
        if semi_window < 10:
            semi_window += 2
        else:
            p_order -= 1
    specs = spectra.shape
    noise = []
    for order in range(specs[0]):
        xarr = np.arange(specs[1])
        spec = spectra[order]
        nonnan = ~np.isnan(spec)
        spec_fit = scipy.signal.savgol_filter(spec[nonnan], 2*semi_window+1, p_order)            # Fit the flux
        res = spec[nonnan] - spec_fit
        noise_order = np.repeat([np.NaN], specs[1])
        for index, index_f in enumerate(xarr[nonnan]):
            x_s = xarr[ max(0,index_f-semi_window):min(specs[1],index_f+semi_window) + 1 ]        # range on full list
            nonnan_s = ~np.isnan(spec[x_s])
            x_s = x_s[nonnan_s]
            left = np.where(x_s == index_f)[0][0]
            res_sub = res[ max(0,index-left):min(len(res),index+(len(x_s)-left)) + 1 ]
            res_sub = percentile_list(res_sub, 1./len(res))         # remove the highest and lowest outlier
            noise_order[index_f] = np.std(res_sub, ddof=1)        # ddof=1 and not ddof=p_order+1, as the latter one doesn't take into account that there are real physical features (absorption lines)
        # !!! Maybe remove the outliers, as a cosmics affect the fit and then the noise in this area is much higher
        noise_order = scipy.signal.medfilt(noise_order, 2*semi_window+1)
        noise.append(noise_order)
        #print len(spec_fit), sum(nonnan), specs, len(xarr)
        #if order>45:
        #    plot_img_spec.plot_points([xarr, xarr[nonnan], xarr[nonnan], xarr, xarr],[spec, spec_fit, res, noise_nofilt, noise_order],['data','fit', 'res', 'noise', 'noisefiltered'],'path',show=True, x_title='Pixel', y_title='Flux', title='' )
        
    return np.array(noise)
    
    
    """# Old, before August 2019 - fits the whole order and then uses part of the fit to calculate stdev
    specs = spectra.shape
    cen_px = int(specs[1]/2)                 # integer of the central pixel
    noise = []
    for order in range(specs[0]):
        xarr = np.arange(specs[1])
        spec = spectra[order]
        nonnan = ~np.isnan(spec)
        noise_order = np.repeat([np.NaN], specs[1])
        if len(xarr[nonnan]) > p_order:                                     # Only do the fit if enough data is available
            poly = polyfit_adjust_order(xarr[nonnan]-cen_px, spec[nonnan], p_order)
            res = spec - np.polyval(poly, xarr-cen_px)                         # Residuals between spectrum and polynomial fit
            for index in range(semi_range, specs[1]-semi_range):
                res_sub = res[ max(0,index-semi_range) : min(specs[1],index+semi_range) + 1]
                res_sub = res_sub[~np.isnan(res_sub)]
                if len(res_sub) >= p_order+2:
                    noise_order[index] = np.std(res_sub, ddof=p_order+1)        # Noise will be high at areas wih many absorption lines
        noise.append(noise_order)
        
    return np.array(noise)"""

def estimate_width(im):
    """
    Finding a guess of the Gaussian width of the traces in an image by using a set of random positions
    :param im: 2d array of float, image with the orders
    :return width: float, median of the measured Gaussian widths
    """
    ims = im.shape
    widths = []
    logger('Info: Estimate the width of the bright features')
    maxtests = max( 200, int(np.prod(ims)/50) )
    for i in range(maxtests):     # Use at least 200 positions in order to be sure to find at least 10
        x = random.randint(int(0.05*ims[0]), int(0.95*ims[0]))
        y = random.randint(int(min(50,0.05*ims[1])), int(max(ims[1]-50,0.95*ims[1])))
        ys = max(0,y-30)
        yr = list(range(ys, min(ims[1],y+31)))
        data = im[x,yr]
        pos_max = data.argmax() + ys       # Maximum in total y (from im)
        widths1 = []
        for w1 in range(2,10,1):
            for w2 in range(1,4,1):
                yy = im[x, max(0,pos_max-w1*w2):min(ims[1],pos_max+w1*w2+1)]
                if yy.shape[0] < 5:
                    continue
                xx = np.arange(yy.shape[0])
                popt = centroid_order(xx, yy, int(yy.shape[0]/2), w1)    # x,y,pos,width ; result: a,x0,sigma,b in a*np.exp(-(x-x0)**2/(2*sigma**2))+b
                
                """title = '{}, {}, {}, {} - {}, {}, {}, {}'.format(x,y,len(yr),pos_max, popt[0],popt[1],popt[2],popt[3])
                label_datapoins = '{}, '.format(popt[2])
                gauss = popt
                label_gauss = 'fit'
                shifts = popt[1]
                label_centroids = 'position'
                plot_gauss_data_center(np.array([xx]), np.array([yy]), [label_datapoins], np.array([gauss]), [label_gauss], np.array([shifts]), [label_centroids], filename='', title=title)"""
                if popt[2] != 0:
                    widths1.append(popt[2])
        if len(widths1) < 1:
            continue
        widths.append(np.median(widths1))
        #print i, x,y, np.median(widths1), np.std(widths1, ddof=1)
        if len(widths) >= 50:
            break
    if len(widths) < 10:
        logger('Warn: Width of the traces could not be determinded, only {0} traces found when testing {1} positions. Assuming 1 pixel as Gaussian width.'.format(len(widths), maxtests ))
        return 1.
    width = np.median(widths)
    logger('Info: The median Gaussian width of the orders when testing {0} positions is {1} pixel.'.format( len(widths), round(width,2) ))
    return width

def find_trace_orders(params, im, imageshape):
    """
    :param params: Dictionary with all the parameters. 'binx', 'biny', and 'polynom_order_apertures' are required
    :param im: 2d numpy array
    return: array of the parameters for the polynomial, orders ordered from left
    return: minimum and maximum position of the trace in dispersion direction
    """
    binx = params['binx']
    biny = params['biny']
    maxFWHM = max(1, int(round(estimate_width(im)*2.35482*1.5)))   # The average Gaussian width, transformed into a FWHM and extendet, as this is the maximum FWHM
    maxshift = max(0.5,0.17*binx)                    # Only 0.15px shift per pixel along one order
    ims = im.shape
    cen_px = int(ims[0]*binx/2)             # central pixel to base the fit on
    #breakvalue = np.percentile(im,15)       # The brightes pixel of one order needs to be brighter than this value, otherwise it won't be identified (lower values -> darker orders are found)
    breakvalue = params['traces_searchlimit_brightness']
    im_orig = copy.deepcopy(im)
    im_traces = np.zeros(ims)               # Masks the orders in order to avoid overlapping traces
    traces = []
    trace_pos = [list(range(ims[0]))]     # fill array with index of slice, later for each order found the y values calculated from the polyfit are added -- this entry is not needed
    #bright10 = np.percentile(im, 5, axis=None)      # background
    #searchnumbers = np.sum(im >= breakvalue+bright10, axis=None)     # At most so many searches
    searchnumbers = np.sum( (im >= breakvalue), axis=None)     # At most so many searches
    searchnumbers = searchnumbers *0.1 / maxFWHM       # FWHM pixel are bright, 90% along the length of an order is brighter than seachnumbers -> less area to search
    for dummy in tqdm(range(int(searchnumbers)), desc='Searching for traces'):      # 2600 is necessary for the tight orders using Clark's lens and a low diffraction prism
        pos_max = np.unravel_index(im.argmax(), ims)
        if im_orig[pos_max] <= breakvalue:
            break
        for i in range(max(0,int(pos_max[0]-100/binx)), min(ims[0],int(pos_max[0]+100/binx+1)) ):
            im[i,max(0,pos_max[1]-maxFWHM): min(ims[1],pos_max[1]+maxFWHM+1)] = breakvalue              #don't find this and sourinding values again
        if pos_max[0] < 100./binx or pos_max[0] > ims[0]-100./binx:                                   # The brightest point of an order shouldn't be close to the border
            continue                # the maximum is too close to the border
        #if np.max(im_traces[pos_max[0],max(0,pos_max[1]-maxFWHM/3):min(ims[1],pos_max[1]+maxFWHM/3+1)]) != 0:        # old, which could result in too wide exclution area
        if np.max(im_traces[pos_max[0],max(0,pos_max[1]-2):min(ims[1],pos_max[1]+3)]) != 0:
            #print 'Too close to order', pos_max
            continue                # the maximum was found too close to an already traced order
        #print 'breakvalue, im_orig[pos_max], pos_max', breakvalue, im_orig[pos_max], pos_max
        expected_positions = []
        # expected_positions = np.average(np.array(trace_pos)[1:,:], axis=0, weights=1/(np.array(trace_pos)[1:,pos_max[0]]-pos_max[1])**2)    # Expected position, without an offset; doesn't work due to bad data
        diff_ord = np.array(trace_pos)[1:,pos_max[0]]-pos_max[1]     # Take care for nans
        #print trace_pos
        #print 'diff_ord, pos_max, np.array(trace_pos)[1:,pos_max[0]]', diff_ord, pos_max, np.array(trace_pos)[1:,pos_max[0]]
        diff_ord[np.isnan(diff_ord)] = -1E10    # small weight
        for px in range(ims[0]):
            data_px = np.array(trace_pos)[1:,px]
            notnan = ~np.isnan(data_px)
            data_px = data_px[notnan]
            if len(data_px) == 0:               # can not calculate the expected position for one px -> don't use expected positions
                expected_positions = []
                break
            diff = diff_ord[notnan]
            expected_positions.append( np.average(data_px - diff, weights=1/diff**2) )
            #if expected_positions[-1] > 1000:
            #    print 'expected_positions[-1], data_px, diff, data_px - diff, 1/diff**2', expected_positions[-1], data_px, diff, data_px - diff, 1/diff**2
        #print expected_positions
        positions = np.array([]).reshape(0,3)
        widths = []
        oldcenter, lastadd = pos_max[1], pos_max[0]
        order_overlap = False
        last_trustworth_position, no_center = 0, 0
        #print pos_max  # Test !
        for i in range(pos_max[0],-1,-1):               # check positions upwards
            center, cen_poly, width, leftmin,rightmin = find_center(im_orig[i,:], int(round(oldcenter)), i, maxFWHM, significance=3.5)       # significance=4.0 tested as useful for HARPS, EXOhSPEC
            #if pos_max[1] >= 2900 and pos_max[1] <= 320:            # bugfixing
            #    #if not ( width != 0 and ( abs(center-oldcenter)<maxshift or (positions.shape[0]==0 and abs(center-oldcenter)<maxshift*2) ) ):
            #    find_center(im_orig[i,:], int(round(oldcenter)), i, maxFWHM, significance=3.0, bugfix=True)
            #    #print(pos_max, i, center, width, leftmin,rightmin, oldcenter, abs(center-oldcenter), maxshift, len(positions))
            if width != 0 and ( abs(center-oldcenter)<maxshift or (positions.shape[0]==0 and abs(center-oldcenter)<maxshift*2) ):        #first entry can be further off
                if im_traces[i,min(ims[1]-1,int(center))] != 0 or im_traces[i,min(ims[1]-1,int(center+1))] != 0:        # this order shouldn't cross another order
                    #print i, "order_overlap"   # Test !
                    order_overlap = True
                    break
                positions = np.vstack([positions, [i, center, 0]])
                widths.append([center-leftmin,rightmin-center, width])
                oldcenter, lastadd, last_trustworth_position, no_center  = center, i, len(positions)-1, 0
            else:#if expected_positions != []:
                #if width !=0:
                #    print i, abs(center-oldcenter), maxshift   # Test !
                no_center += 1
                if abs(oldcenter - 0.) < maxshift/2. or abs(oldcenter - ims[1]) < maxshift/2.:  # if the trace leaves shortly the CCD
                    no_center -= 0.6
                if no_center >= max(5, int(40/binx), 1*maxFWHM):         # stop searching if too many fits are unseccessful, as otherwise the fit might drift off
                    break
                center1, cen_poly, width, leftmin,rightmin = find_center(im_orig[i,:], int(round(oldcenter)), i, maxFWHM, significance=3.0, bugfix=False)
                if width != 0 and abs(center1-oldcenter) < maxshift:
                    positions = np.vstack([positions, [i, center1, 1]])
                if lastadd-i == 5 and expected_positions != []:        #add entries every 10 empty line in order to avoid the fit of the trace going off
                    positions = np.vstack([positions, [i, oldcenter, 1]])
                    lastadd = i
                if positions.shape[0] > 0:
                    if expected_positions != []:                        # use the solution of the other traces to keep following this trace
                        oldcenter = positions[last_trustworth_position][1] - expected_positions[int(positions[last_trustworth_position,0])] + expected_positions[max(0,i-1)]
                    else:                                               # use the last values of this trace to keep following it
                        good_values = np.where( (positions[:,0] <= i+30) & (positions[:,2] == 0) )[0]
                        if len(good_values) > 5:
                            poly = np.polyfit(positions[good_values,0], positions[good_values,1], 1)
                            oldcenter = np.polyval(poly, i-1)
                    if abs(oldcenter - positions[last_trustworth_position,1]) > maxFWHM * 1.5:    # Stop, if the shift is in danger of touching the next order
                        #print 'too far off', oldcenter, positions[last_trustworth_position,1], no_center
                        break
        if positions.shape[0] != 0:
            oldcenter = positions[0][1]
        else:
            oldcenter = pos_max[1]
        lastadd = pos_max[0]
        last_trustworth_position, no_center = 0, 0
        for i in range(pos_max[0]+1,ims[0]):               # check positions downwards
            center, cen_poly, width, leftmin,rightmin = find_center(im_orig[i,:], int(round(oldcenter)), i, maxFWHM, significance=3.5)       # significance=4.0 tested as useful for HARPS, EXOhSPEC
            #if pos_max[1] >= 2900 and pos_max[1] <= 320:            # bugfixing
            #    #if not ( width != 0 and ( abs(center-oldcenter)<maxshift or (positions.shape[0]==0and abs(center-oldcenter)<maxshift*2) ) ):
            #    find_center(im_orig[i,:], int(round(oldcenter)), i, maxFWHM, significance=3.0, bugfix=True)
            #    #print(pos_max, i, center, width, leftmin,rightmin, oldcenter, abs(center-oldcenter), maxshift, len(positions), len(expected_positions), last_trustworth_position)
            if width != 0 and ( abs(center-oldcenter)<maxshift or (positions.shape[0]==0 and abs(center-oldcenter)<maxshift*2) ):
                if im_traces[i,min(ims[1]-1,int(center))] != 0 or im_traces[i,min(ims[1]-1,int(center+1))] != 0 or order_overlap == True:
                    #print i, "order_overlap"   # Test !
                    order_overlap = True
                    break
                positions = np.vstack([positions, [i, center, 0]])
                widths.append([center-leftmin,rightmin-center, width])
                oldcenter, lastadd, last_trustworth_position, no_center  = center, i, len(positions)-1, 0
            else:#if expected_positions != []:
                #if width !=0:
                #    print i, abs(center-oldcenter), maxshift   # Test !
                no_center += 1
                if abs(oldcenter - 0) < maxshift/2. or abs(oldcenter - ims[1]) < maxshift/2.:  # if the trace leaves shortly the CCD
                    no_center -= 0.6
                if no_center >= max(5, int(40/binx), 1*maxFWHM):         # stop searching if too many fits are unseccessful, as otherwise the fit might drift off
                    break
                center1, cen_poly, width, leftmin,rightmin = find_center(im_orig[i,:], int(round(oldcenter)), i, maxFWHM, significance=3.0, bugfix=False)
                if width != 0 and abs(center1-oldcenter) < maxshift:
                    positions = np.vstack([positions, [i, center1, 1]])
                if i-lastadd == 5 and expected_positions != []:
                    positions = np.vstack([positions, [i, oldcenter, 1]])
                    lastadd = i
                if positions.shape[0] > 0:
                    if expected_positions != []:
                        oldcenter = positions[last_trustworth_position,1] - expected_positions[int(positions[last_trustworth_position,0])] + expected_positions[min(ims[0]-1,i+1)]
                    else:
                        good_values = np.where( (positions[:,0] >= i-30) & (positions[:,2] == 0) )[0]
                        if len(good_values) > 5:
                            poly = np.polyfit(positions[good_values,0], positions[good_values,1], 1)
                            oldcenter = np.polyval(poly, i+1)
                    if abs(oldcenter - positions[last_trustworth_position,1]) > maxFWHM * 1.5:    # Stop, if the shift is in danger of touching the next order
                        break
        if len(positions) < ims[0]/4.:                                   # first this check and afterwards the check for order_overlap == True, as otherwise the order_overlap floats the output
            if len(positions) > ims[0]/10. and order_overlap == False:
                logger('Warn: the order around {0}, {1}, was skipped, as only {2} centroids along the order have been found'.format(pos_max[0]*binx,pos_max[1]*biny, len(positions)), show=False)
            continue
        if order_overlap == True:
            if len(positions) < ims[0]/2.:
                logger('Warn: the order around {0}, {1}, was skipped, as it overlapped with another order'.format(pos_max[0]*binx,pos_max[1]*biny))
                continue
            positions = positions[:-1,:]
            del widths[-1]
            logger('Warn: the order around {0}, {1}, was shortened, as it otherwise started overlapping with another order'.format(pos_max[0]*binx,pos_max[1]*biny))
        positions = np.array(sorted(positions, key=operator.itemgetter(0)))
        width = [np.nanmean(percentile_list(np.array(widths)[:,0],0.1)),
                 np.nanmean(percentile_list(np.array(widths)[:,1],0.1)),
                 np.nanmean(percentile_list(np.array(widths)[:,2],0.1))]      #average left, right, and gaussian width
        #print positions
        # Fit to the original image
        pfs = np.polyfit(positions[:,0]*binx - cen_px, positions[:,1]*biny, params['polynom_order_apertures'])
        # Remove values which were added only to follow the curve
        positions = positions[positions[:,2]==0,0:2]
        # Calculate the values from the fitted curve for each slices
        yarr = np.polyval(pfs, positions[:,0]*binx - cen_px)/(biny+0.0)
        trace_pos_add = []      # add the information where the trace is located
        for i in range(ims[0]):
            j = np.where(positions[:,0]==i)[0]
            if len(j) == 1:
                trace_pos_add.append(yarr[j[0]])
            else:
                trace_pos_add.append(np.nan)
        trace_pos.append(trace_pos_add)
        yarr = (np.round(yarr)).astype(int)   # if using this, than also change 'np.array(trace_pos[0])' to 'positions[:,0]'
        #y = (np.round(np.polyval(pfs, np.array(trace_pos[0])*binx)/biny+i)).astype(int)
        for i in range(-int(width[2]*2),int(width[2]*2)+1):         # fills im and im_traces with values in order to avoide searching orders there
            goodpos = (yarr+i>=0) & (yarr+i<=ims[1]-1)
            im_traces[positions[goodpos,0].astype(int),yarr[goodpos]+i] = 0.2
            im[positions[goodpos,0].astype(int),yarr[goodpos]+i] = breakvalue
        goodpos = (yarr>=0) & (yarr<=ims[1]-1)
        im_traces[positions[goodpos,0].astype(int),yarr[goodpos]] = 1
        traces.append( [cen_px] + list(pfs) + [min(positions[:,0])*binx, min(imageshape[0],(max(positions[:,0])+1)*binx-1), np.polyval(pfs, cen_px)] )
        logger('Step: order found at central pixel: {0} / {3} (frame/trace). The trace was identified between Pixels {1} and {2}'.format(round(np.polyval(pfs, cen_px),1),
                                                round(traces[-1][-3]), round(traces[-1][-2]), round(np.polyval(pfs, np.mean(traces[-1][-3:-1]) ),1) ))
    # sort
    traces = np.array(traces)
    cen_px2 = np.mean([ np.max(traces[:,-3]), np.min(traces[:,-2]) ])      # Center of the area covered by all orders
    for order in range(traces.shape[0]):
        pfs = traces[order, 0:-3]
        if cen_px2 >= traces[order, -3] and cen_px2 <= traces[order, -2]:
            traces[order, -1] = np.polyval(pfs[1:], cen_px2 - pfs[0])     # Position in Cross-dispersion direction
        else:
            logger('Warn: Please check that all orders have been identified correctly')
            if abs(cen_px2 - traces[order, -3]) < abs(cen_px2 - traces[order, -2]):     # Use the closest value for which the trace is still defined
                traces[order, -1] = np.polyval(pfs[1:], traces[order, -3] - pfs[0])
            else:
                traces[order, -1] = np.polyval(pfs[1:], traces[order, -2] - pfs[0])
    traces = np.array(sorted(traces, key=operator.itemgetter(-1)))
    #plot_img_spec.plot_image(im+im_traces*np.max(im), ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'cleared')
    if traces.shape[0] < 2:
        logger('Error: Only found {0} traces. Please check that the binned image ({1}) looks as expected.'.format(traces.shape[0], params['logging_trace1_binned'])+\
               'Please check that the right image was used to search for the traces (e.g. {0}) '.format(params['logging_traces_im_binned'])+\
               'After you selected different file(s) for parameter "trace1_rawfiles", please run the following command before restarting the pipeline'+\
               '\nrm {0}'.format( params['master_trace1_filename'] ))
    logger('Info: {0} orders found and traced'.format(traces.shape[0]))
    return traces[:,0:-3], traces[:,-3].astype(int), traces[:,-2].astype(int)            # polyfits, xlows, xhighs

def adjust_trace_order(kwargs):
    """
    This is the routine which does the work for adjust_trace_orders to allow multiprocessing. It performs the step for one order.
    """
    [order, kwargs] = kwargs
    trace_pos = kwargs['trace_pos']
    maxFWHM = kwargs['maxFWHM']
    im = kwargs['im']
    params = kwargs['params']
    maxshift = kwargs['maxshift']
    ims = kwargs['ims']
    binx = kwargs['binx']
    cen_px = kwargs['cen_px']
    biny = kwargs['biny']
    imus = kwargs['imus']
    
    if True:
        lastadd = trace_pos[0,0]
        positions, widths_o, shifts = [], [], []
        for j,i in enumerate(trace_pos[0,:].astype(int)):           # each pixel, is is the position in the array and i the real position on the CCD
            oldcenter = trace_pos[order+1,j]                        # +1 because first entry are pixels
            if oldcenter > -maxFWHM:                                # center is expected to be in the image
                cen_gauss, cen_poly, width, leftmin,rightmin = find_center(im[i,:], int(round(oldcenter)), i, maxFWHM, border_pctl=params['width_percentile'], significance=3.0, polymax=True)
                #logger('{0} {1} {2} {3}'.format(order,j,cen_gauss,cen_poly), show=False, printarrayformat=[], printarray=[], logfile='gauss_vs_poly.txt')
                center = cen_poly
            else:
                center, width, lastadd = 0,0, i
            #if width == 0 and oldcenter > -maxFWHM:
            #if width != 0 and abs(center-oldcenter)>=3:
            #    print center,oldcenter,i, leftmin,rightmin
            if width != 0 and abs(center-oldcenter)<maxshift:
                #print 'maxFWHM,order,j,i,leftmin,center, rightmin, width',maxFWHM,order,j,i,leftmin,center, rightmin, width
                if leftmin < center-width*2.35482*3 or leftmin > center:        # unreasonable small position (3xFWHM)
                    leftmin = np.nan
                if rightmin > center+width*2.35482*3 or rightmin < center:      # unreasonable large position (3xFWHM)
                    rightmin = np.nan
                #leftmin = max(leftmin, center-width*2.35482*3)          # anything bigger than 
                #rightmin = min(rightmin, center+width*2.35482*3)
                positions.append([i, center, 0, leftmin, rightmin])                 # leftmin and rightmin are positions
                widths_o.append([i, center-leftmin, rightmin-center, width])        # center-leftmin, rightmin-center are widths
                lastadd = i
                shifts.append(center-oldcenter)
            elif oldcenter > -maxFWHM and i-lastadd == 10:          # If inside the image add the center after every 10 data points
                positions.append([i, oldcenter, 1, np.nan, np.nan])
                lastadd = i
            #else:
            #    print center, oldcenter, maxshift, maxFWHM
                
        #print len(positions), 
        if widths_o == [] or len(shifts) < ims[0]/4.:                # Not enough data for that order
            #print 'len(positions),len(shifts)', len(positions),len(shifts)
            return [], [], [], [], [], [], []
        positions = np.array(positions)
        widths_o = np.array(widths_o)
        # Fit the center
        pfs = np.polyfit(positions[:,0]*binx - cen_px, positions[:,1]*biny, params['polynom_order_apertures'])
        #centerfit.append(pfs)
        # Fit the borders
        good_vaules = ~np.isnan(positions[:,3])
        if sum(good_vaules) > len(good_vaules)*0.05:        # Position available at least 5% of positions
            leftfit_o =  np.polyfit(positions[good_vaules,0]*binx - cen_px, positions[good_vaules,3]*biny, params['polynom_order_apertures'])
            widthleft = np.nanmean(percentile_list(widths_o[:,1],0.1))          # Average of the widths except the 10% highest and lowest
        else:                                                   # left side at FWHM of line
            widthleft = np.nanmean(percentile_list(widths_o[:,3],0.1)) * 2.35482
            leftfit_o = copy.deepcopy(pfs)
            leftfit_o[-1] -= widthleft
        good_vaules = ~np.isnan(positions[:,4])
        if sum(good_vaules) > len(good_vaules)*0.05:        # Position available at least 5% of positions
            rightfit_o = np.polyfit(positions[good_vaules,0]*binx - cen_px, positions[good_vaules,4]*biny, params['polynom_order_apertures'])
            widthright = np.nanmean(percentile_list(widths_o[:,2],0.1))         # Average of the widths except the 10% highest and lowest
        else:                                                   # left side at FWHM of line
            widthright = np.nanmean(percentile_list(widths_o[:,3],0.1)) * 2.35482
            rightfit_o = copy.deepcopy(pfs)
            rightfit_o[-1] += widthright
        centerfit_o = [[cen_px] + list(pfs), [cen_px] + list(leftfit_o), [cen_px] + list(rightfit_o)]
        center = np.polyval(pfs, widths_o[:,0]*binx - cen_px)
        widths_o = [widthleft, widthright, np.nanmean(percentile_list(widths_o[:,3],0.1))]      #average left, right, and gaussian width
        xlows_o = max(0,min(positions[:,0])*binx-5)
        xhighs_o = min(imus[0],(max(positions[:,0])+1)*binx+5)
        avg_shifts_o = np.mean(shifts)

    return centerfit_o, leftfit_o, rightfit_o, xlows_o, xhighs_o, widths_o, avg_shifts_o

def adjust_trace_orders(params, im, im_unbinned, pfits, xlows, xhighs):
    """
    Re-traces the position of the orders
    :param params: Dictionary with all the parameters. 'binx', 'biny', and 'polynom_order_apertures' are required
    :param im: 2d numpy array of floats, binned image
    :param im_unbinned: original image without binning
    :param pfits: list, length same as number of orders, polynomial values
                      (output of np.polyval)
                      i.e. p where:
                      p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
    :param xlows: list, length same as number of orders, the lowest x pixel, minimum is 0
                   (wavelength direction) used in each order
    :param xhighs: list, length same as number of orders, the highest x pixel, maximum is the number of pixel+1 -> range(xlows,xhighs) can cover the whole CCD
                   (wavelength direction) used in each order
    :return centerfit: 3d numpy array of floats. First dimensions has the same length as number of orders.
                        The second dimension is dividing between center of the trace and the left and right border of the trace.
                        Third deminsion contains the central pixel and the polynomial values
    :return xlows: 1d array of integers, length same as number of orders, minimum position of the trace in dispersion direction
    :return xhighs: 1d array of integers, length same as number of orders, maximum position of the trace in dispersion direction
    :return widths: 2d array of floats, length same as number of orders, width of the trace in cross-dispersion direction. 
                    Each line contains the giving the width to the left, to the right, and the Gaussian width
    """
    binx = params['binx']
    biny = params['biny']
    ims = im.shape
    imus = im_unbinned.shape
    cen_px = int(ims[0]*binx/2)
    if len(pfits.shape) == 3:     # in case a finalised trace will be used (e.g. cp master_traces_sci.fits logging/master_traces_sci_searched.fits
        pfits = pfits[:,0,:]
    # maxFWHM = 25/biny for old lens
    # maxFWHM = 15/biny+2            # !! Determine in values of space between orders
    maxFWHM = max(1, int(round(estimate_width(im)*2.35482*1.5)))   # The average Gaussian width, transformed into a FWHM and extendet, as this is the maximum FWHM
    if len(pfits.shape) > 1:                # move at maximum by 20% of the minimum space between the orders, but at old binning in cross-dispersion direction is allowed
        mask = ( (xlows-10 <= np.percentile(xlows, 80.)) & (xhighs+10 > np.percentile(xhighs, 20.)) )   # use the orders, which cover most of the CCD
        maxshift = max(params['bin_search_apertures'][1], int(min(np.abs(pfits[mask,-1][1:] - pfits[mask,-1][:-1] ))/5) )
    else:
        maxshift = params['bin_search_apertures'][1]    # old binning in cross-dispersion direction
    trace_pos = [list(range(int(min(xlows)),int(max(xhighs)),binx))]      #x-values as first entry in the array (pixel in the unbinned CCD
    for order, pfit in enumerate(pfits):
        trace_pos.append(np.polyval(pfit[1:], trace_pos[0]-pfit[0])/(biny+0.0))       # fit of the trace in binned frame
        trace_pos[-1][:int((xlows[order]-min(xlows))/binx)] = -1e10     # do not try on bad data
        trace_pos[-1][int((xhighs[order]-min(xlows)+1)/binx):] = -1e10  # do not try on bad data
    trace_pos = np.array(trace_pos)
    trace_pos[0,:] = (trace_pos[0,:]/binx).astype(int)
    trace_pos = trace_pos[:, trace_pos[0,:] < ims[0] ]              # Due to binning, entries outside the binned images can appear in trace_pos[0,:]
    #traces = []
    centerfit, leftfit, rightfit, xlows, xhighs, widths, avg_shifts = [], [], [], [], [], [], []
    kwargs = dict(trace_pos=trace_pos, maxFWHM=maxFWHM, im=im, params=params, maxshift=maxshift, ims=ims, binx=binx, cen_px=cen_px, biny=biny, imus=imus)
    orders = list(range(trace_pos.shape[0]-1))
    logger('Step: Adjusting orders')
    if params['use_cores'] > 1:
        orders2 = [[order, kwargs] for order in orders]
        p = multiprocessing.Pool(params['use_cores'])
        results = p.map(adjust_trace_order, orders2)
        p.terminate()
    else:
        results = []
        for order in tqdm(orders, desc='Adjusting orders'):
            results.append( adjust_trace_order([order, kwargs]) )
        
    for result in results:
        centerfit_o, leftfit_o, rightfit_o, xlows_o, xhighs_o, widths_o, avg_shifts_o = result
        if len(centerfit_o) > 0:
            centerfit.append(centerfit_o)
            leftfit.append(leftfit_o)
            rightfit.append(rightfit_o)
            xlows.append(xlows_o)
            xhighs.append(xhighs_o)
            widths.append(widths_o)
            avg_shifts.append(avg_shifts_o)
    """
    for order in tqdm(orders, desc='Adjusting orders'):
        centerfit_o, leftfit_o, rightfit_o, xlows_o, xhighs_o, widths_o, avg_shifts_o = adjust_trace_order([order, kwargs])
        if len(centerfit_o) > 0:
            centerfit.append(centerfit_o)
            leftfit.append(leftfit_o)
            rightfit.append(rightfit_o)
            xlows.append(xlows_o)
            xhighs.append(xhighs_o)
            widths.append(widths_o)
            avg_shifts.append(avg_shifts_o)"""

    if len(widths) == 0:
        logger('Error: When adjusting the traces, all traces were rejected due to poor fit. '+\
               'Please check the binned image {0}: Is the orientation and the binning right? Are the orders covering at least half of the CCD (in dispersion correction)'.format(params['logging_trace1_binned'])+\
               'Please check that the right image was used to search for the traces (e.g. {0}) '.format(params['logging_traces_im_binned'])+\
               'After you selected different file(s) for parameter "trace1_rawfiles", please run the following command before restarting the pipeline'+\
               '\nrm {0} {1}'.format( params['logging_traces_binned'], params['master_trace1_filename'] ))
    logger( ('Info: traces of the {0} apertures adjusted. The average shift of the individual apertures was between '+\
             '{1} and {2} pixel between the searching of the traces and this solution. '+\
             'The maximum allowed shift was {3} pixel.').format(len(centerfit), np.round(np.min(avg_shifts),1), np.round(np.max(avg_shifts),1), np.round(maxshift,1) ))
    
    return np.array(centerfit), np.array(xlows).astype(int), np.array(xhighs).astype(int), np.array(widths)

def adjust_width_orders(center, left, right, w_mult):
    """
    Adjusts the width of the orders by stretching the distance between the border of the trace and the centre of the trace
    :param center, left, right: 1d arrays of floats, center of order and borders of the order
    :param w_mult: list of float, multiplicator in order to adjust the width of the order to the left and right
    :return left, right: 1d arrays of floats, borders of the order
    """
    diff = center - left                # equivalent to w_mult == 1
    left = center - diff * w_mult[0]
    diff = right - center               # equivalent to w_mult == 1
    right = center + diff * w_mult[1]
    return left, right

def asign_bad_px_value(input_value, badpx_mask, saturated_mask, data_range):
    """
    :param input_value: already asigned value for the bad pixel mask
    :param badpx_mask: n-d array (normally 2d array) of int or float, with 0 to mark the bad pixels
    :param im: n-d array (normally 2d array) of int or float, image from which the extraction will be done
    :param data_range: number, list, or array of integers, to mark the extracted area
    """
    if np.nanmin(badpx_mask[data_range]) == 0:
        input_value = 0.2
    elif np.nanmin(saturated_mask[data_range]) == 0:
        input_value = 0.1
    
    return input_value

def no_tqdm(input, desc=''):
    """
    In order to switch between tqdm on and off, this is needed
    """
    return input

def extract_order(kwarg):
    """
    Extract the spectra for one order
    """
    [pp, kwargs] = kwarg
    xarr = kwargs['xarr']
    image = kwargs['image']
    badpx_mask = kwargs['badpx_mask']
    saturated_mask = kwargs['saturated_mask']
    pfits = kwargs['pfits']
    xlows = kwargs['xlows']
    xhighs = kwargs['xhighs']
    widths = kwargs['widths']
    w_mult = kwargs['w_mult']
    offset = kwargs['offset']
    var = kwargs['var']
    maxpx = image.shape[1]-1                    # End of the image in cross-dispersion direction
    #for pp in plot_tqdm(range(pfits.shape[0]), desc='Extract Spectrum'):        # pp is order
    if True:
        #print('\t Order {0}'.format(pp + 1))
        # find nearest y pixel
        yarr = np.polyval(pfits[pp, 0, 1:], xarr-pfits[pp, 0, 0])+offset
        xarr1 = np.arange(xlows[pp], xhighs[pp])
        ospecs = np.repeat([np.NaN], image.shape[0])
        ogood_px_mask = np.repeat([0.0], image.shape[0])            # Array needs to be float in order to save 0.1, 0.2, ... correctly
        #w = widths[pp]
        # getting the lower and upper bounds of the trace, doing everything in np, to speed up things
        #old: lowers = yarr - w[2]*w_mult
        #old: uppers = yarr + w[2]*w_mult
        lowers = np.polyval(pfits[pp, 1, 1:], xarr-pfits[pp, 1, 0]) + offset 
        uppers = np.polyval(pfits[pp, 2, 1:], xarr-pfits[pp, 2, 0]) + offset
        #print( 2580, pp, np.sum(np.isnan(lowers)), np.sum(np.isnan(uppers)), np.sum(uppers<lowers) )
        lowers, uppers = adjust_width_orders(yarr, lowers, uppers, [w_mult, w_mult])          # Adjust the width of the orders
        #print( 2582, pp, np.sum(np.isnan(lowers)), np.sum(np.isnan(uppers)), np.sum(uppers<lowers) )
        widths_o = uppers - lowers      # can be below 0, but shouldn't be!! e.g. 20201104/con4_50mu_bc
        #print pp,np.nanmean(uppers-lowers), np.nanmedian(uppers-lowers), np.nansum(uppers-lowers)
        if w_mult == 0.0:                                                       # only the central full pixel
            widths_o = widths_o * 0 + 1.
            order_in_frame = ((yarr[xarr1] <= maxpx) & (yarr[xarr1] >= 0) )     # Due to order shifting the order might be shifted to the ouside of the images, ignore these values
            yarr = np.round(yarr).astype(int)
            for xrow in xarr1[order_in_frame]:
                ssum = image[xrow, yarr[xrow]]                      # Verified 20190122
                good_px = asign_bad_px_value(1, badpx_mask, saturated_mask, (xrow, yarr[xrow]) )
                ospecs[xrow] = ssum
                ogood_px_mask[xrow] = good_px  
        elif var == 'linfit':         # new solution, better than the one below, tested with many values in excel
            lowersf = (np.ceil(lowers+.5)).astype(int)       #Full lowest pixel of the order
            uppersf = (np.floor(uppers-.5)).astype(int)      #Full highest pixel of the order
            lowers[lowers<-0.5] = -0.5
            uppers[uppers>maxpx+0.5] = maxpx+0.5
            lowersf[lowersf<0] = 0
            uppersf[uppersf>maxpx] = maxpx
            order_in_frame = ((lowersf[xarr1] <= maxpx) & (uppersf[xarr1] >= 0) & (lowers[xarr1] <= uppers[xarr1]) )     # Due to order shifting the order might be shifted to the ouside of the images, ignore these values, ignore also areas where the width of the order is kind of negative (however, when extracting narrow orders, this will lead to a problem
            # Loop around x-px, this is necessary, as the number of full pixels between lowersf and uppersf varies 
            for xrow in xarr1[order_in_frame]:
                good_px = 1
                ssum = 0
                if lowersf[xrow] <= uppersf[xrow]:
                    ssum = np.sum(image[xrow, lowersf[xrow]:uppersf[xrow]+1])
                    good_px = asign_bad_px_value(good_px, badpx_mask, saturated_mask, (xrow, range(lowersf[xrow],uppersf[xrow]+1)) )
                # Get fractions of pixel from a polynom, if they are not outside the frame borders (spline actes weired between points, when there is an outlier)
                #if w[2]*w_mult < 0.5:             # There will be no full pixel, which messes up everything
                if uppers[xrow] - lowers[xrow] < 1.0:       # There will be no full pixel, which messes up everything
                    fracpixparam = [[yarr[xrow],'b']]
                else:
                    fracpixparam = [[lowers[xrow], 'l'] , [uppers[xrow],'u']]
                for [borderpos, pos] in fracpixparam:
                    x = np.arange(max(0,np.floor(borderpos-1)), min(maxpx,np.ceil(borderpos+1))+1, dtype=int)
                    if len(x) <= 1:
                        continue
                    good_px = asign_bad_px_value(good_px, badpx_mask, saturated_mask, (xrow, x) )
                    y = image[xrow, x]
                    weight = 1./(np.abs(x-borderpos)+0.1)**2                   # weight the values so that the data around the fraction of the pixel is used most
                    p = polyfit_adjust_order(x, y, max(1,len(x)-3), w=weight)
                    poly = np.poly1d(p)
                    polyint = poly.integ()
                    if pos == 'l':
                        if lowersf[xrow] != 0:
                            ssum += polyint(lowersf[xrow]-0.5) - polyint(lowers[xrow])
                    elif pos == 'u':
                        if uppersf[xrow] != maxpx:
                            ssum += polyint(uppers[xrow]) - polyint(uppersf[xrow]+0.5)
                    elif pos == 'b':
                        ssum += polyint(uppers[xrow]) - polyint(lowers[xrow])
                    else:
                        print('Programming error around line 2200')
                ospecs[xrow] = ssum
                ogood_px_mask[xrow] = good_px
        else:           # old solution
            # +.5/-.5 because if light between pixels 3.5 and 4.5 should be extracted, that would be exactly all light in pixel 4; if light between 3 and 5: 0.5*3 + full px 4 + 0.5*5; between 2.6 and 5.4: 0.9*3 + full px 4 + 0.9*5
            lowersf = (np.ceil (lowers+.5)).astype(int)                         # Full pixels, round up
            uppersf = (np.floor(uppers-.5)).astype(int)                         # Full pixels, round down
            lowersf[lowersf<0] = 0                                              # If the order is at the boundary, only extract until the boundary
            uppersf[uppersf>maxpx] = maxpx                                      # If the order is at the boundary, only extract until the boundary
            order_in_frame = ((lowersf[xarr1] <= maxpx) & (uppersf[xarr1] >= 0) & (lowers[xarr1] <= uppers[xarr1]))     # Due to order shifting the order might be shifted to the ouside of the images, ignore these values
            lowersr = (np.round(lowers)).astype(int)                            # left pixel, of which a fraction should be extracted
            uppersr = (np.round(uppers)).astype(int)                            # right pixel, of which a fraction should be extracted
            # Loop around x-px, this is necessary, as the number of full pixels between lowersf and uppersf varies 
            for xrow in xarr1[order_in_frame]:
                good_px = 1
                ssum = 0
                if lowersf[xrow] <= uppersf[xrow]:
                    ssum = np.sum(image[xrow, lowersf[xrow]:uppersf[xrow]+1])
                    good_px = asign_bad_px_value(good_px, badpx_mask, saturated_mask, (xrow, range(lowersf[xrow],uppersf[xrow]+1)) )
                    # Tested on 15/3/19 for individual entries and ranges of entries: image[xrow, lowersf[xrow]:uppersf[xrow]+1] is the same as data_range = (xrow, range(lowersf[xrow],uppersf[xrow]+1)) ; image[data_range]
                    """if min(badpx_mask[xrow, lowersf[xrow]:uppersf[xrow]+1]) == 0:                       # Bad pixel in extracted data
                        good_px = 0.2
                    elif max(image[xrow, lowersf[xrow]:uppersf[xrow]+1]) >= params['max_good_value']:   # Saturated pixel in extracted data
                        good_px = 0.1"""
                #print 'xrow, ssum, ...', xrow, ssum, lowers[xrow], uppers[xrow], lowersf[xrow], uppersf[xrow]+1, 
                # Get fractions of pixel, if they are not outside the frame borders
                if lowersr[xrow] >= 0:
                    ssum += image[xrow, lowersr[xrow]]*( (1.5 - lowers[xrow]%1)%1 )     # lowers[xrow]%1=[0.5,0.0,0.6,0.4] -> [1.0,1.5,0.9,1.1] -> [.0,.5,.9,.1]
                    good_px = asign_bad_px_value(good_px, badpx_mask, saturated_mask, (xrow, lowersr[xrow]) )
                #print 'l',ssum, lowersr[xrow],(1.5 - lowers[xrow]%1)%1,
                if uppersr[xrow] <= maxpx:      # maxpx is shape-1
                    ssum += image[xrow, uppersr[xrow]]*( (0.5 + uppers[xrow]%1)%1 )     # uppers[xrow]%1)=[0.5,0.0,0.6,0.4] -> [1.0,0.5,1.1,0.9] -> [.0,.5,.1,.9]
                    good_px = asign_bad_px_value(good_px, badpx_mask, saturated_mask, (xrow, uppersr[xrow]) )
                #print 'u', ssum, uppersr[xrow],(0.5 + uppers[xrow]%1)%1
                if lowersr[xrow] == uppersr[xrow] and lowersr[xrow] >= 0 and uppersr[xrow] <= maxpx:
                    ssum += image[xrow, lowersr[xrow]]*( - 1 )     # remove the one added too much when using the fractions on the same pixel
                    #print 'f', ssum, uppersr[xrow],-1
                ospecs[xrow] = ssum
                ogood_px_mask[xrow] = good_px
                #if pp == 57 and xrow == 1350:
                #    print pp,xrow, image[xrow, lowersf[xrow]:uppersf[xrow]+1], image[xrow, lowersr[xrow]], (0.5-lowers[xrow]%1), image[xrow, uppersr[xrow]], (uppers[xrow]%1-.5), ssum
    
    return ospecs, ogood_px_mask, widths_o
     
def extract_orders(params, image, pfits, xlows, xhighs, widths, w_mult, offset=0, var='standard', plot_tqdm=True, header=dict(), plot_traces=False ):
    """
    Extract the spectra from each of the orders and return in a numpy array
    :param params: Dictionary with all the parameters. 'maxgood_value' is required
    :param image: numpy 2D array, containing the CCD image
    :param pfits: list, length same as number of orders, polynomial values
                      (output of np.polyval)
                      i.e. p where:
                      p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
    :param xlows: list, length same as number of orders, the lowest x pixel
                   (wavelength direction) used in each order
    :param xhighs: list, length same as number of orders, the highest x pixel
                   (wavelength direction) used in each order
    :param widths: (not needed anymore) 2d list, length same as number of orders, each entry contains left border, right border, and Gaussian width of the lines, as estimated in the master flat,
                   not needed anymore because of the new format of pfits, giving the lower and upper border of the trace
    :param w_mult: the space between the lower and upper boundary compared to the center of the traces is adjusted by w_mult
                   (old: w_mult * widths (Gaussian) defines the range on either side of the trace to be extracted)
                   if 0 or 0.0 then only the full central pixel will be extracted
    :param offset: allows to shift the spectrum
    :param var: can be 'standard' or 'linfit'.
        'standard': extraction of the fraction of the border pixel
        'linfit': fit a polynomial around the border of the extraction width to extract with higher precission
    :return spec: array,  length same as number of orders, each containing
                  a 1D numpy array with the spectrum of each order in
    :return good_px_mask: same format as spec, contains 1 for good data and values between 0 and 1 for bad data (saturated, bad pixel, ...)
    :return widths_m: same format as spec, contains the extracted with for each spectral pixel

    """
    if plot_tqdm:
        plot_tqdm = tqdm
    else:
        plot_tqdm = no_tqdm
    logger('Step: extracting spectra', show=False)
    # Prepare bad-pixel map
    badpx_mask_name = 'badpx_mask'
    if badpx_mask_name not in calimages:
        calimages[badpx_mask_name] = read_badpx_mask(params, image.shape)
    badpx_mask = calimages[badpx_mask_name]
    # Prepare saturated map
    saturated_mask = ( image < params['max_good_value'] )       # 1 when not saturated, analogue to bad-pixel
    filename = header.get('HIERARCH HiFLEx orig', header.get('HiFLEx orig', ''))
    if plot_traces and params['logging_traces_im'].find('*') != -1 :
        fname_short = header.get('HIERARCH HiFLEx orid', header.get('HiFLEx orid', '')).replace('.fits','').replace('.fit','')
        fname = '{0}_{1}.{2}'.format( params['logging_traces_im'].replace('*', fname_short) )
        plot_traces_over_image(image, fname, pfits, xlows, xhighs, widths, offset=offset, w_mult=w_mult)
    if '{0}_saturated'.format(filename) in calimages.keys():
        saturated_mask = saturated_mask*0 + 1
        x, y = calimages['{0}_saturated'.format(filename)]
        saturated_mask[x,y] = 0
    spec, good_px_mask, widths_m = [], [], []
    xarr = range(0, image.shape[0])
    #maxpx = image.shape[1]-1                    # End of the image in cross-dispersion direction
    kwargs = dict( xarr=xarr, image=image, badpx_mask=badpx_mask, saturated_mask=saturated_mask, 
                   pfits=pfits, xlows=xlows, xhighs=xhighs, widths=widths, w_mult=w_mult, offset=offset, var=var)
    orders = list(range(pfits.shape[0]))
    if w_mult != 0 and params['use_cores'] > 1 and multiprocessing.current_process().name == 'MainProcess':     #var == 'linfit':
        orders2 = [[pp, kwargs] for pp in orders]
        p = multiprocessing.Pool(params['use_cores'])
        results = p.map(extract_order, orders2)
        p.terminate()
    else:
        results = []
        for pp in plot_tqdm(orders, desc='Extract Spectrum'):        # pp is order
            results.append( extract_order([pp, kwargs]) )
    for result in results:
        ospecs, ogood_px_mask, widths_o = result
        spec.append(ospecs)
        good_px_mask.append(ogood_px_mask)
        widths_m.append(widths_o)
    
    return np.array(spec), np.array(good_px_mask), np.array(widths_m)

def shift_orders_multicore(kwarg):
        [order, params, im, sci_tr_poly, xlows, xhighs, steps, oldwidths, in_shift] = kwarg
        xarr = range(max(0,int(xlows[order])), min(im.shape[0],int(xhighs[order])), steps)
        yarr = np.polyval(sci_tr_poly[order, 0, 1:], xarr-sci_tr_poly[order, 0, 0])+in_shift          #apply input shift
        widths, shifts, shift_map = [], [], []
        for px,oldcenter in enumerate(yarr):        #each pixel
            #if px>6000 or px%100 == 0:
            #    print(order, px, len(xarr), yarr.shape, im.shape)
            #    print(xarr[px])
            cen_gauss, cen_poly, width, leftmin,rightmin = find_center(im[xarr[px],:], int(round(oldcenter)), xarr[px], 2*params['maxshift_orders'], border_pctl=0, significance=3.0, polymax=True)
            if width*2.35482*3 < 2*params['maxshift_orders']:
                cen_gauss, cen_poly, width, leftmin,rightmin = find_center(im[xarr[px],:], int(round(cen_gauss)), xarr[px], width*2.35482*3, border_pctl=params['width_percentile'], significance=3.0, polymax=True)
            #if width != 0:
            #    print 'order,px, cen_gauss, cen_poly, width, leftmin,rightmin',order,px,center, width, leftmin,rightmin
            #center = cen_gauss         # Old, but this might not cover the centre well
            center = cen_poly
            if width != 0 and abs(center-oldcenter) < params['maxshift_orders']:
                widths.append([center-leftmin, rightmin-center, width])
                #shifts.append(center-oldcenter)
                #shift_map[max(0,int(xarr[px]-steps/2)):min(im.shape[0],int(xarr[px]+steps/2))+1,order] = center-oldcenter
                shifts.append(center-oldcenter)
                shift_map.append([ max(0,int(xarr[px]-steps/2)), min(im.shape[0],int(xarr[px]+steps/2))+1, order ])
        if len(widths) == 0:
            width = oldwidths[order]
            problem_order = order
            #twidths.append(oldwidths[order])
            #problem_order.append(order)
        else:
            widths = np.array(widths)
            width = np.zeros(3)*np.nan
            for ii in range(3):
                if np.sum(~np.isnan(widths[:,ii])) > 0:                         # left and right might be np.nan all the time
                    width[ii] = np.nanmean(percentile_list(widths[:,ii],0.1))   # average left, right, and gaussian width
            problem_order = None
        return width, problem_order, np.array(shifts), shift_map

def shift_orders(im, params, sci_tr_poly, xlows, xhighs, oldwidths, in_shift = 0):
    """
    Determines the shift in spacial direction between the current flat {im} and an available solution. This is done by fitting the center again
    :param im: 2d numpy array with the flat lines
    :param params: Dictionary with all the parameters. 'maxshift_orders' is required
    :param sci_tr_poly: list, length same as number of orders, polynomial values
                      (output of np.polyval)
                      i.e. p where:
                      p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
    :param xlows: list, length same as number of orders, the lowest x pixel
                   (wavelength direction) used in each order
    :param xhighs: list, length same as number of orders, the highest x pixel
                   (wavelength direction) used in each order
    :param oldwidths: 2d list, length same as number of orders, each entry contains left border, right border, and Gaussian width of the lines, estimated for the previous solution
    :param in_shift: offset to the positions giving in sci_tr_poly
    :return shift+in_shift: absolute shift
    :return twidths: widths of the orders, one line per order, containing the average size to the left, right, and the gaussian width
    :return shift_map: 2d image with number of pixels in wavelength direction and number of orders in y. It's 0 where no shift could be determined
    :return shift_error: Standard deviation of the shifts.
    """
    steps = 10                              # is worse than binning, only every step will be checked
    shifts, twidths = [], []
    problem_order = []
    shift_map = np.zeros((im.shape[0],sci_tr_poly.shape[0]))
    logger('Step: Searching for shifts of the orders', show=False) 
    kwargs = []
    for order in range(sci_tr_poly.shape[0]):
        kwargs.append([order, params, im, sci_tr_poly, xlows, xhighs, steps, oldwidths, in_shift])
    if params['use_cores'] > 1 and multiprocessing.current_process().name == 'MainProcess':
        logger('Step: Searching for shifts of the orders (multicore).')
        p = multiprocessing.Pool(params['use_cores'])
        data = p.map(shift_orders_multicore, kwargs)
        p.terminate()
    else:
        data = []
        for kwarg in tqdm(kwargs, desc='Searching for shifts of the orders'):        # pp is order
            data.append( shift_orders_multicore(kwarg) )
    for entry in data:          # width, problem_order, np.array(shifts), shift_map
        shifts += list(entry[2])        # center-oldcenter
        twidths.append(entry[0])
        if entry[1] is not None:
            problem_order.append(entry[1])
        for ii, entry3 in enumerate(entry[3]):
            shift_map[entry3[0]:entry3[1], entry3[2] ] = entry[2][ii]  
    if len(shifts) > 0:
        shift = np.mean(percentile_list(np.array(shifts),0.1))
        shift_error = np.std(percentile_list(np.array(shifts),0.1), ddof=1)
    else:
        shift = 0
    twidths = np.array(twidths)
    printarrayformat = ['%1.1i', '%4.2f', '%4.2f']
    printarray = np.array([ range(sci_tr_poly.shape[0]), np.array(oldwidths)[:,2], twidths[:,2] ]).T
    problem_text = ''
    if len(problem_order) > 0:
        problem_text = ' The shift in the following orders could not be measured: {0}'.format(problem_order)
        if len(problem_order) > 0.5*sci_tr_poly.shape[0]:
            problem_text += '. No shift could be measured for {0} out of {1} orders.'.format(len(problem_order), sci_tr_poly.shape[0])
            shift_error = -1
    logger('Info: The shift of the traces was determinded to {0} +- {2}, including the inputshift of {1}.{3}'.format(round(shift+in_shift,2), in_shift, round(shift_error,2), problem_text ))
    #logger('Info: The gaussian width of the orders changed: order\told\tnew', printarrayformat=printarrayformat, printarray=printarray)
    return shift+in_shift, twidths, shift_map, shift_error

def find_bck_px(im, pfits, xlows, xhighs, widths, w_mult):
    """
    Creates the map of the pixel which can be used for background correction
    - Maybe at some point it might be worth to add the traces of the second fiber in order to allow background correction with the science file
    :param im: 2d array, which needs to have the same size as the images. The content is irrelevant
    :param pfits: list, length same as number of orders, polynomial values
                      (output of np.polyval)
                      i.e. p where:
                      p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
    :param xlows: list, length same as number of orders, the lowest x pixel
                   (wavelength direction) used in each order
    :param xhighs: list, length same as number of orders, the highest x pixel
                   (wavelength direction) used in each order
    :param w_mult: w_mult * widths (Gaussian) defines the range on either side of the trace to be extracted
    """
    im = im*0+1
    ims = im.shape
    #xarr = range(0,ims[0])
    # pfits, widths = update_tr_poly_width_multiplicate(pfits, widths, [w_mult, w_mult], xlows-10, xhighs+10) # not suitable as calculation has to be done again
    for pp in range(pfits.shape[0]):
        xarr = range(max(0,xlows[pp]-10),min(ims[0],xhighs[pp]+10))             # extend the area slightly
        yarr = np.polyval(pfits[pp, 0, 1:], xarr-pfits[pp, 0, 0])
        left = np.polyval(pfits[pp, 1, 1:], xarr-pfits[pp, 1, 0])
        right = np.polyval(pfits[pp, 2, 1:], xarr-pfits[pp, 2, 0])
        left, right = adjust_width_orders(yarr, left, right, [w_mult, w_mult])              # Adjust width
        right += 1
        left[ left < 0] = 0
        right[ right > ims[1] ] = ims[1]
        left = np.floor(left).astype(int)
        right = np.ceil(right).astype(int)
        for i in range(len(xarr)):
            im[xarr[i], left[i]:right[i]] = 0
    return im

def measure_background_noise(im):
    """
    Measures the background noise in an image by binning it into 10x10 sized areas and measuring the median and std of the noise in these 100px small areas
    :param im: 2d array of float, the traces are already removed and set to np.nan
    :return median(sim): float, median of the standard deviations from parts of the images
    :return std(sim): float, variation of the standard deviations from parts of the images
    """
    #if np.nanmedian(im) >= 10:          # Images with dark correction or images with zero bias level might be set to completely NaN otherwise
    #    im[im == 0.0] = np.nan                                                              # Replace zeros with NaNs
    #    im[im == 0] = np.nan                                                              # Replace zeros with NaNs
    #plot_img_spec.plot_image(im, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'orig')
    bim, gim, sim = bin_im(im, [10, 10])          # gim: number of eleements, sim: standard deviation
    #plot_img_spec.plot_image(bim, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'bined_im')
    #plot_img_spec.plot_image(gim, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'bined_datapoints')
    #plot_img_spec.plot_image(sim, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'bined_std')
    mask = ~np.isnan(bim)
    while np.sum(bim[mask] - np.percentile(bim[mask],30) >= 3 * np.nanstd(bim[mask], ddof=1)) > 0:  # Excluded brightest areas (MRES: unidentiefied orders define the noise otherwise
        #plot_img_spec.plot_image(mask, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'mask')
        bim_0 = copy.copy(bim)
        bim_0[np.isnan(bim)] = -1E10            # otherwise the line below will produce "RuntimeWarning: invalid value encountered"
        mask[bim_0 - np.percentile(bim[mask],30) >= 3 * np.nanstd(bim[mask], ddof=1)] = False
    sim[~mask] = np.nan
    gim[~mask] = 0
    for threshold in range(90,30,-5):
        if np.sum( gim >= threshold ) >= gim.shape[0] * gim.shape[1] * 0.10:              # 10% of the binned data is available
            break
    sim = sim[ gim >= threshold ]                  # only use data where 90% of the data is defined
    return np.nanmedian(sim), np.nanstd(sim, ddof=1)
    
    """ first idea
    :return std_full: float, std of the full data set. Might be too high, if a background gradient is present
    
    good_values = ~np.isnan(im)
    im = im[good_values]
    ims = im.shape
    print ims, ' <- should be onedimensional'
    std_full = np.std(im, ddof=1)                            # Standard deviation of the whole data
    im = im[:(ims/100)*100].reshape(100, ims/100)          # reshape in order to get a better value for standard deviation
    std_part = np.std(im, ddof=1, axis=1)
    print std_part.shape ' <- should not be 100. If 100, then axis=0'
    return std_full, np.median(std_part), np.std(std_part, ddof=1) """

def find_shift_images_2d(im, im_ref, shift_range):
    """
    Correlates two images against each other
    :param im:
    :param im_ref:
    :param shift_range: 1d list or array with 4 intergers. Gives the range of how to shift the image, first in x, followed by y
    :return popt: 1d array of 7 floats. The parameters to fit twoD_Gaussian
    """
    ims = im.shape
    if ims != im_ref.shape:
        logger('Error: Image and comparison image do not have the same dimension: {0} versus {1}'.format(ims, im_ref.shape))
    data = []
    shiftxs = range(shift_range[2],shift_range[3]+1)
    shiftys = range(shift_range[0],shift_range[1]+1)
    for shiftx in shiftxs:
        for shifty in shiftys:
            posx1 = 0      + max(0, shiftx)
            posx2 = ims[0] - max(0,-shiftx)
            posy1 = 0      + max(0, shifty)
            posy2 = ims[1] - max(0,-shifty)
            data.append([shiftx, shifty, posx2-posx1, posy2-posy1, np.nan])
    data = np.array(data)
    minxrange = np.min(data[:,2])   # Smallest image
    minyrange = np.min(data[:,3])
    fluxdiff = []
    for index, [shiftx, shifty, rangex, rangey, dummy] in enumerate(data):
        posx1 = int( 0      + max(0, shiftx) )
        posx2 = int( ims[0] - max(0,-shiftx) )
        posy1 = int( 0      + max(0, shifty) )
        posy2 = int( ims[1] - max(0,-shifty) )
        porx1 = int( 0      + max(0,-shiftx) )
        porx2 = int( ims[0] - max(0, shiftx) )
        pory1 = int( 0      + max(0,-shifty) )
        pory2 = int( ims[1] - max(0, shifty) )
        diffim = im[posx1:posx2, posy1:posy2] - im_ref[porx1:porx2, pory1:pory2]
        # print('should be the same;', diffim.shape, rangex, rangey)        # Yes, it is :)
        fluxes = []
        ims2 = diffim.shape
        extrax = int(rangex - minxrange)
        extray = int(rangey - minyrange)
        for ii in range(extrax + 1):
            for jj in range(extray + 1):      # For a small shift, the resulting image is bigger than for a big shift -> avoid impact total flux
                fluxes.append(np.sum(np.abs(diffim[ii:ims2[0]-(extrax-ii), jj:ims2[0]-(extray-jj)])))
        data[index,4] = -np.median(fluxes)              # using a negative flux, so the minimum becomes a maximum, which only can be handled by the script in order to fit a gaussian
        #print data[index,0:4], data[index,4]
    pos_max = data[:,4].argmax()
    initial_guess = [np.max(data[:,4])-np.min(data[:,4]), data[pos_max,0], data[pos_max,1], 2, 2, 0, np.min(data[:,4])] # amplitude, xo, yo, sigma_x, sigma_y, theta, offset
    x, y = np.meshgrid(shiftxs, shiftys)
    popt, pcov = curve_fit(twoD_Gaussian, (x,y), data[:,4], p0=initial_guess)
    
    return popt

def find_shift_images(params, im, im_ref, sci_traces, w_mult, cal_tr_poly, extract=True, im_head=dict()):
    """
    Finds the shift between two images by cross corellating both. The images are substracted and the total flux [sum(abs())] is calculated. At the best position, a minimum will be reached
    :param im: 2d array with the image, for which the shift should be calculated
    :param im_ref: 2d array with the reference image
    :return shift: shift of the image (float)
    """
    if params['extraction_shift'] == 0.0:
        return 0.0, im_head
        
    shift = None
    fname = im_head.get('HiFLEx orid', im_head.get('HIERARCH HiFLEx orid', '') )
    if len(fname) > 3 and os.path.isfile(params['logging_crossdispersion_shift']):
        data_str = read_textfile_remove_double_entries(params['logging_crossdispersion_shift'], [float, float, float, float, str], delimiter='\t', replaces=['\n'], equal_indexes=[4], warn='')
        index = ( fname == data_str[:,4] )
        if np.sum(index) == 1:
            [shift, width, min_shifts, max_shifts] = data_str[index][0,0:4].astype(float)
            logger('Info: The shift between this frame and the reference frame in Cross-Dispersion direction was read from {4}: shift = {0} px, gaussian width = {1} px. A shift between {2} and {3} pixel was tested.'.format(shift, width, min_shifts, max_shifts, params['logging_crossdispersion_shift'] ))
            
    if shift is None:
        [sci_tr_poly, xlows, xhighs, widths] = sci_traces
        # Find the maximum shift by using the space between science and calibration orders
        for start_order in range(int(sci_tr_poly.shape[0]/3)):             # If orders are to close in the red, then no shift can be determinded -> might ignore the first third of the red orders
            shifts = []
            if sci_tr_poly.shape[0] == 1:   # only one order, e.g. when using live_extraction
                shift = 10 * max(widths[:,2])
            else:                           # normal case
                if params['arcshift_side'] != 0:
                    shifts.append( np.abs(sci_tr_poly[start_order  :  ,0,-1] - cal_tr_poly[start_order  :  ,0,-1]) )    # not when science and calibration at the same position
                shifts.append( np.abs(sci_tr_poly[start_order+1:  ,0,-1] - cal_tr_poly[start_order  :-1,0,-1]) )        # calibration fiber could be left or right of science fiber
                shifts.append( np.abs(sci_tr_poly[start_order  :-1,0,-1] - cal_tr_poly[start_order+1:  ,0,-1]) )        # calibration fiber could be left or right of science fiber
                if sci_tr_poly.shape[0] >= 4:           # If enough data, then do a polyfit
                    for i in range(len(shifts)):
                        poly = np.polyfit(range(len(shifts[i])), shifts[i], 2)
                        shifts[i] = np.polyval(poly, range(len(shifts[i])) )
                for i in range(len(shifts)):
                    shifts[i] = np.min(shifts[i])
                shift = min(shifts)
            shift = min(shift, 20 * max(widths[:,2]))
            if shift >= 3:              
                break
        shift = min(shift, params['extraction_shift'] *3.0)     # Added 20201027, 3 as division by 2 in next line and a bit wider is necessary for Gaussian fit
        range_shifts = [-int(shift/2),int(shift/2)]
        ims1 = im.shape[1]
        #shifts = range(min(range_shifts),max(range_shifts)+1)
        fluxdiff, shifts = [], []
        oldcen, olderr = np.nan, np.nan
        if extract:                                                     # takes a bit longer
            logger('Step: checking if the current image is shifted compared to the reference frame in which the traces were searched (in cross-dispersion direction)')
        no_orders = len(xlows)-start_order
        step = max( 1, int(no_orders/20) )              # only use every 20ths order
        mask = list(range(start_order,no_orders,step))            # only use part of the orders to speed up calculations
        for shift in range(max(np.abs(range_shifts))+1):
            for pm in [-1, +1]:
                if extract:         # extract and find the maximum flux
                    spec, good_px_mask, extr_width = extract_orders(params, im, sci_tr_poly[mask,:,:], xlows[mask], xhighs[mask], widths[mask,:], w_mult, pm*shift, plot_tqdm=False)       
                    fluxes = np.nansum(spec, axis=1)
                    fluxes[np.isnan(fluxes)] = np.nanmin(fluxes)          # replace nans with the minimum flux to avoid problems caused at the borders at the image
                else:               # Use the difference of the 2 images (old way)
                    fluxes = []
                    diffim = im[:,max(0,pm*shift):ims1-max(0,-pm*shift)] - im_ref[:,max(0,-pm*shift):ims1-max(0,pm*shift)]      # Tested that the minus signs are at the right place
                    ims2 = diffim.shape[1]
                    temp = max(np.abs(range_shifts).astype(int))-abs(shift)
                    for i in range(temp+1):         # For a small shift, the resulting image is bigger than for a big shift -> avoid impact total flux
                        #print shift, i, diffim.shape, diffim[:,i:ims2-(temp-i)].shape, np.sum(np.abs(diffim[:,i:ims2-(temp-i)]))
                        fluxes.append(-np.sum(np.abs(diffim[:,i:ims2-(temp-i)])))                   # using a negative flux, so the minimum becomes a maximum, which only can be handled by the script in order to fit a gaussian
                fluxdiff.append(np.median(fluxes))      # Changed from mean to median on 20190515
                shifts.append(pm*shift)
                #print 'shift, pm, fluxdiff[-1], fluxes', shift, pm, fluxdiff[-1], fluxes
                if shift == 0:      # For 0 shift + and - will be the same
                    break
            if len(fluxdiff) >= 7 or abs( shift - abs(min(range_shifts)) ) < 0.001:       # Run at least after the last shift
                popt = centroid_order(shifts,fluxdiff,shifts[np.argmax(fluxdiff)],max(shifts)-min(shifts))    #a,x0,sigma,b in a*np.exp(-(x-x0)**2/(2*sigma**2))+b
                #print shifts,fluxdiff,shifts[np.argmax(fluxdiff)],max(shifts)-min(shifts)
                if abs(oldcen - popt[1]) < 0.05 and abs(olderr - popt[2]) < 0.2 and popt[2] < 0.3 * (max(range_shifts) - min(range_shifts)):        # Stop early, if good enough precission
                    break
                oldcen, olderr = popt[1], popt[2]

        if popt[2] > 0.5 * (max(range_shifts) - min(range_shifts)) or popt[1] < min(range_shifts) or popt[1] > max(range_shifts) or abs(popt[1]) > params['extraction_shift']:
            comment = (' The calculated values (shift = {0} px, width = {1} px) were not very precise or'+\
                       ' outside of the allowed range (Parameter extraction_shift: {4},'+\
                       ' borders due to next order [{2},{3}]) .').format(round(popt[1],2), round(popt[2],2), 
                                            min(range_shifts), max(range_shifts), params['extraction_shift'] )
            shift, width = 0.0, 0.0
        elif not extract:                                                                               # keep everything if not extracting
            shift, width, comment = round(popt[1],3), round(popt[2],3), ''
        else:                                                                                           # Find the place of the maximum flux (and redo the gaussian fit)
            logger('Step: Do the finetuning of the center and find the maximum', show=False)
            #nshifts, fluxes = [], []
            nshifts, fluxes = shifts, fluxdiff
            fshifts = np.linspace(popt[1] - 0.1*popt[2], popt[1] + 0.1*popt[2], 7)                  # Search in the central area again
            for dummy in range(5):                                                                  # repeat several times in order to get closer and closer to the maximum
                for fshift in fshifts:
                    if np.min( np.abs( np.array(nshifts) - fshift ) ) < 0.005:                  # Don't do the fit again, if it was done before
                        continue
                    spec, good_px_mask, extr_width = extract_orders(params, im, sci_tr_poly[mask,:,:], xlows[mask], xhighs[mask], widths[mask,:], w_mult, fshift, plot_tqdm=False)
                    fluxes.append( np.median( np.nansum(spec, axis=1) ) )
                    nshifts.append(fshift)
                sort_arr = np.argsort(fluxes)
                best_shifts = np.array(nshifts)[sort_arr][-3:]               # best 3 values
                #print shifts, fluxes, best_shifts
                fshifts = np.linspace( min(best_shifts), max(best_shifts), 7)
            # Do the gaussian fit again, although changes are less than 1%
            popt = centroid_order(shifts,fluxdiff,shifts[np.argmax(fluxdiff)],max(shifts)-min(shifts))    #a,x0,sigma,b in a*np.exp(-(x-x0)**2/(2*sigma**2))+b
            shift, width = round(popt[1],3), round(popt[2],3)
            shift_poly, poly,l,r = fit_poly_center(fluxdiff, x=shifts)
            shift_poly = round(shift_poly, 3)
            mshift = round( best_shifts[-1], 3)                 # find the maximum
            comment = ' The center of the gauss was shifted by {0} px, the center of a third order polynomial was shifted by {2} px, and the maximum flux found at a shift of {1} px.'.format(shift, mshift, shift_poly)
            im_head['HIERARCH HiFLEx CD_SHIFT_GAUSS']  = (shift, 'Shift of the Gauss-centre in C-D [px]')        # (value, comment)
            if not np.isnan(shift_poly):
                im_head['HIERARCH HiFLEx CD_SHIFT_POLY']   = (shift_poly, 'Shift of the polynom maximum in C-D [px]')        # (value, comment)
            im_head['HIERARCH HiFLEx CD_SHIFT_MAXFLUX']  = (mshift, 'Shift of the maximum flux in C-D [px]')        # (value, comment)
            #shift = np.mean([shift,mshift]     # if not symetrical traces then, a big shift
            if not np.isnan(shift_poly):
                shift = shift_poly
            elif width == 0 and min(shifts) < mshift < max(shifts):      # If Gaussian fit didn't work, or if shift is close to borders, then mshift might be dogy
                shift = mshift                                      # use the maximum flux as shift
        min_shifts, max_shifts = min(shifts), max(shifts)
        logger('Info: The shift between this frame and the reference frame in Cross-Dispersion direction is {0} px. The gaussian width is {1} px. A shift between {2} and {3} pixel was tested.{4}'.format(shift, width, min_shifts, max_shifts, comment ))
        fname = im_head.get('HiFLEx orid', im_head.get('HIERARCH HiFLEx orid', '') )
        if (shift !=0 or width !=0) and len(fname) > 3:
            add_text_to_file('{0}\t{1}\t{2}\t{3}\t{4}'.format(shift, width, round(min_shifts,1), round(max_shifts,1), fname), params['logging_crossdispersion_shift'] ) 
    
    im_head['HIERARCH HiFLEx CD_SHIFT']  = (shift, 'Applied Shift in Cross-dispersion [px]')        # (value, comment)
    im_head['HIERARCH HiFLEx CD_S_WDTH'] = (width, 'Width of shift in Cross-dispersion [px]')        # (value, comment)
    im_head['HIERARCH HiFLEx CD_S_MIN']  = (min_shifts, 'Minimum tested shift in Cross-disp. [px]')        # (value, comment)
    im_head['HIERARCH HiFLEx CD_S_MAX']  = (max_shifts, 'Maximum tested shift in Cross-disp. [px]')        # (value, comment)
    
    return shift, im_head

def arc_shift_multicore(kwargs):
    [shift, params, im, pfits, xlows, xhighs, widths, wmult, plot_tqdm] = kwargs
    arc_spec, good_px_mask, extr_width = extract_orders(params, im, pfits, xlows, xhighs, widths, 0, shift, plot_tqdm=False)          # w_mult == 0: only central pixel
    flux = np.nansum(arc_spec, axis=1)
    flux[np.isnan(flux)] = np.nanmin(flux)          # replace nans with the minimum flux
    
    return flux

def arc_shift(params, im, pfits, xlows, xhighs, widths):
    """
    Determines the difference between the science orders and arc orders in cross-dispersion direction
    :param params: Dictionary with all the paramerters. arcshift_range (is made positive when reading the configuration), arcshift_side, and arcextraction_width_multiplier will be used
    :param im: 2d image in which the arc orders are iluminated
    :param pfits: list, length same as number of orders, polynomial values
                      (output of np.polyval)
                      i.e. p where:
                      p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
    :param xlows: list, length same as number of orders, the lowest x pixel
                   (wavelength direction) used in each order
    :param xhighs: list, length same as number of orders, the highest x pixel
                   (wavelength direction) used in each order
    :param widths: 2d list, length same as number of orders, each entry contains left border, right border, and Gaussian width of the lines, as estimated in the master flat
    :return shifts: 1d array with the shift for each order
    """
    if params['arcshift_side'] == 0:
        return np.array(xlows)*0        # All shifts are 0
    sigma = 2.5
    # w_mult = params['arcextraction_width_multiplier']
    w_mult = 0 # In which pixel row is the highest flux
    orders = np.arange(pfits.shape[0])
    cen_pos, cen_pos_diff = [], []
    arcshift_range = np.array(params['arcshift_range'])
    ims = im.shape
    if len(orders) >= 10:                    # Determine the search area automatically
        # Search everything between the maximum difference between consecutive orders, excluding the width of the orders
        xarr = np.arange(ims[0])
        cen_pos = np.zeros((len(orders), ims[0])) * np.nan                                          # Initiate a Nan-array for all px and orders
        for order in orders:
            mask = ( (xarr >= xlows[order]) & (xarr <= xhighs[order]) )                             # Covered by the order
            cen_pos[order, mask] = np.polyval(pfits[order, 0, 1:], xarr[mask]-pfits[order, 0, 0])   # The traces
        cen_pos_diff = np.nanmedian( np.abs(cen_pos[1:,:] - cen_pos[:-1,:]) , axis=1 )              # The difference between the traces into one number
        cen_pos_diff_medfilt = scipy.signal.medfilt(cen_pos_diff, 5)                                # Takes into account the NaNs
        # Test everything between the 4 times Gauss width and next Order-4*width (some values will be ignored, when fitting the gauss):   # factor 2 to avoid the science fiber
        for i in range(2):                                              # To run a second time in order to account for the widths[:,2]
            arcshift_range = np.array([ 0                               + np.median(widths[:,2])*w_mult*2, \
                                        np.nanmax(cen_pos_diff_medfilt) - np.median(widths[:,2])*w_mult*2  \
                                     ]).astype(int)
            if arcshift_range[1] - arcshift_range[0] >= 15:                                 # have at least 15 pixel
                break
            widths[:,2] /= 1.5                                                              # otherwise use a smaller width for the moment
        cen_pos_diff = np.insert(cen_pos_diff, 0, cen_pos_diff[0])                          # create the same length array
    #print cen_pos, cen_pos_diff
    arcshifts = params['arcshift_side'] * np.abs(arcshift_range)                            # Translate into the right shifts
    arcshifts = np.arange(min(arcshifts),max(arcshifts)+1, 1)
    kwargs = []
    for shift in arcshifts:
        kwargs.append([shift, params, im, pfits, xlows, xhighs, widths, w_mult, False])
    if params['use_cores'] > 1 and multiprocessing.current_process().name == 'MainProcess':
        logger('Step: Search for the shift of the calibration traces, compared to the science traces (multicore).')
        p = multiprocessing.Pool(params['use_cores'])
        fluxes = p.map(arc_shift_multicore, kwargs)
        p.terminate()
    else:
        fluxes = []
        for kwarg in tqdm(kwargs, desc='Search for the shift of the calibration traces, compared to the science traces'):        # pp is order
            fluxes.append( arc_shift_multicore(kwarg) )
    fluxes = np.array(fluxes)
    # Find the center in each order: where is the flux the highest
    gauss, goodvalues = [], []
    label_datapoins, label_gauss, label_poly, label_centroids = [],[],[],[]
    w_mult = max(1, w_mult)
    for order in orders:
        if len(cen_pos_diff) > 0:                          # Automatic determination of the search area
            for i in range(5):                          # If the cen_pos_diff is NaN for this order then check neightboring orders
                orderi = min(len(cen_pos_diff)-1,order+i)
                if not np.isnan(cen_pos_diff[orderi]):
                    break
                orderi = max(0,order-i)
                if not np.isnan(cen_pos_diff[orderi]):
                    break
            goodpos = [0]
            w_mult_test = w_mult*2./0.9
            while len(arcshifts[goodpos]) < widths[order,2]*w_mult_test:           # Enough data to have at least one extra order between the 2 orders
                w_mult_test *= 0.9
                goodpos = ( (np.abs(arcshifts) > widths[order,2]*w_mult_test) & (np.abs(arcshifts) < cen_pos_diff[orderi]-widths[order,2]*w_mult_test) )   # factor 2 to avoid the science fiber
                #print( arcshifts[goodpos], widths[order,2]*w_mult_test, cen_pos_diff[orderi]-widths[order,2]*w_mult_test )
                if w_mult_test < 0.1 * w_mult*2:
                    goodpos = ( (np.abs(arcshifts) >= min(params['arcshift_side']*np.array(params['arcshift_range']))) & (np.abs(arcshifts) <= max(params['arcshift_side']*np.array(params['arcshift_range']))) )
                    break
            x, y = arcshifts[goodpos], fluxes[goodpos,order]
            #print order,orderi,x
        else:
            x, y = arcshifts, fluxes[:,order]
        popts = []
        for centeroffset in range(0,int(len(x)/2-widths[order,2]),3):                # Test for different positions of the line
            for sign in [-1,+1]:
                width_g = min(15,max(x)-min(x))
                #data_range = range( max( 0, len(x)/2+sign*centeroffset-width_g ), min( len(x), len(x)/2+sign*centeroffset+width_g+1 ) )
                #print np.mean(x)+sign*centeroffset, width_g, data_range
                popt = centroid_order(x, y, np.mean(x)+sign*centeroffset, width_g )    #input: x,y,center, width; result: a,x0,sigma,b in a*np.exp(-(x-x0)**2/(2*sigma**2))+b
                #print 'fit',order, sign*centeroffset, np.mean(x)+sign*centeroffset, data_range[0], data_range[-1], width_g,popt 
                if np.all(popt != [0,0,0,0]):
                    diff = np.sum(np.abs(oneD_gauss(x,popt) - y))           # Used only to determine the best fit for this order
                    popts.append([diff, list(popt)])
                    #print order,diff, popt, np.mean(x)+sign*centeroffset, min(15,max(x)-min(x))
                if centeroffset == 0:       # independent of sign
                    break
        if len(popts) != 0:
            goodvalues.append(True)
            popts = sorted(popts, key=operator.itemgetter(0))               # Minimum deviation from the fit
            popt = popts[0][1]
            mask = ( np.abs(x-popt[1]) < 5*popt[2] )                        # Within +- 5* width
            cen_poly, poly,l,r = fit_poly_center(y[mask], x[mask], border_pctl=50)
            if np.isnan(cen_poly):       #cen_poly - gauss[-1][1] < 3:         
                cen_poly = popt[1]
                poly = [np.nan]*4
            gauss.append(popt+[cen_poly,l,r]+list(poly))      #shift=popt[1]
        else:
            goodvalues.append(False)
            gauss.append([np.nan]*11)
        label_datapoins.append('extracted flux in order {0}'.format('%2.2i'%order))
        label_gauss.append('fitted gauss in order {0}'.format('%2.2i'%order))
        label_poly.append('fitted 3rd order polynomial\n in order {0}'.format('%2.2i'%order))
        label_centroids.append('final center in order {0}'.format('%2.2i'%order))
    
    gauss = np.array(gauss)
    if (orders[goodvalues]).shape[0] < 3:
        shift = round(np.mean(params['arcshift_side']*np.array(params['arcshift_range'])),2)
        logger('Warn: A shift to the calibration orders could not be measured. The average value of the input parameters "arcshift_range" are used instead: {0}'.format(shift))
        shifts = np.repeat([shift], orders.shape[0])
        poly, diff, min_gauss, max_gauss, used_orders = [0], [0], 0, 0, 0
    else:
        # Fit a polynomial to the center
        goodval, poly = sigmaclip(gauss[:,4], x=orders, nc=2, sub=goodvalues, ll=sigma, lu=sigma, repeats=5, x_cen=0)
        shifts = np.round(np.polyval(poly, orders),2)
        diff = shifts[goodvalues*goodval] - gauss[goodvalues*goodval,1]
        min_gauss, max_gauss = round(min(gauss[goodvalues&goodval,2]),2), round(max(gauss[goodvalues&goodval,2]),2)
        used_orders = np.sum(goodval)
    # Log the results
    arcshifts = np.repeat([arcshifts],len(orders),axis=0)
    title = 'Determining the shift of the calibration traces'
    plot_gauss_data_center(arcshifts, fluxes.T, label_datapoins, gauss[:,:4], label_gauss, shifts, label_centroids, params['logging_find_arc_traces'], title=title, poly=gauss[:,5:], label_poly=label_poly)
    logger(('Info: The shift between science traces and calibration orders is between {0} and {1} px. '+\
            'The average residuals to the 2nd order polynomial is {6} px (standard deviation of {4} px), {5} orders were used. '+\
            'The gaussian width of the arc lines in spacial direction is between {2} and {3} px').format(min(shifts), max(shifts), 
                    min_gauss, max_gauss, round(np.std(np.abs(diff), ddof=len(poly)),2), used_orders, round(np.mean(np.abs(diff)),2) ))
    logger('Info: The parameters of the polynomial are: {0}'.format(poly), show=False)
    return shifts

def identify_emission_lines_single_order(kwargs):
        [ order, im, im_short, im_badpx, im_short_badpx, maxFWHM ] = kwargs
        lines = []
        ims = im.shape
        xarr = np.arange(ims[1])
        yarr = im[order,:]
        notnans = ~np.isnan(yarr)
        if np.sum(notnans) < 10:         # At the border of the chip the arc traces might be outside of the CCD
            return lines
        yarr1 = yarr[notnans]
        xarr1 = xarr[notnans]
        notline_pos, pfit = sigmaclip(yarr1, x=xarr1, nc=12, sub=[], ll=10, lu=2.2, repeats=20)  #orders, sigma low, sigma high ; tested that 2.5 is too high for some lines (UNe), If lines are not identified correctly, the problem isn't here, probably
        line_pos = ~notline_pos
        for i in range(1,len(line_pos)-1):
            if (line_pos[i-1:i+2]==[False, True, False]).all():     #exclude lines over only one px
                line_pos[i] = False
        #print 'xarr1[line_pos], pfit',xarr1[line_pos], pfit
        already_found = yarr*0
        small_yarr = yarr1[line_pos]
        for i in range(len(small_yarr)):
            pos_max = np.argmax(small_yarr)
            small_yarr[pos_max] = min(small_yarr)-1     #do not find the line again
            pos_real = xarr1[line_pos][pos_max]
            if already_found[pos_real] == 1:
                continue
            # get the data to which a gauss should be fitted
            range_arr = range(max(0,pos_real-maxFWHM*2),min(ims[1],pos_real+maxFWHM*2+1))
            # check if line is saturated and if so use the short exposure time
            y_data = yarr[range_arr]
            yarr_use = copy.copy(yarr)
            if im_short is not None and im_badpx is not None and not (im == im_short).all():
                #print order, pos_max, np.sum(im_badpx[order,range_arr]==0), np.sum(im_badpx[order,range_arr]==0.1), np.sum(im_badpx[order,range_arr]==0.2), np.sum(im_badpx[order,range_arr]==1), np.sum(im_badpx[order,range_arr]!=-1)
                if 0.1 in im_badpx[order,range_arr]:
                    if im_short_badpx is not None:
                        if 0.1 in im_short_badpx[order,range_arr]:          # ignore saturated lines in the short exposure
                            continue
                    y_data = im_short[order,range_arr]
                    yarr_use = copy.copy(im_short[order,:])
                    #if np.sum(already_found[range_arr]) == 0:           # Only print once in the area
                    #    print('Used the short exposure time to find the centroid of the line in order {1} @ px {0}'.format(pos_real, order))
            # Fit a double Gauss 
            #popts = oneD_blended_gauss
            # fit the gauss, pos_real is the expected center
            if True:       # old way until v1.0.0 (including)
                popt = centroid_order(xarr[range_arr],y_data,pos_real,maxFWHM, significance=2, blended_gauss=False)    #a,x0,sigma,b in a*np.exp(-(x-x0)**2/(2*sigma**2))+b
                centre, width, height = popt[1], popt[2], popt[0]
            else:           # using polynomial of order 3, that is not better, especially if FWHM is just few pixel
                popt = centroid_order(xarr[range_arr],y_data,pos_real,maxFWHM, significance=2, blended_gauss=False)
                cen_gauss, cen_poly, width, leftmin,rightmin = find_center(yarr_use, pos_real, order, FWHM, border_pctl=10, border_pctl_poly=50, significance=2.0, polymax=True, dummy=popt)
                if not np.isnan(cen_poly):      centre = cen_poly
                elif not np.isnan(cen_gauss):   centre = cen_gauss
                else:                           centre, width = 0,0
                height = np.max(yarr_use[max(0, int(round(centre-width))):min(ims[1],int(round(centre+width))+1)])
            #if order == 2 and pos_real>100 and pos_real<13150:
            #    print xarr[range_arr],np.round(y_data,0).astype(int), pos_real, popt, cen_gauss, cen_poly, width, leftmin,rightmin
            #if order == 3:
            #    print("exit 3117")
            #    exit(100)
            if centre == 0 or width > maxFWHM or int(round(centre)) >= ims[1] or int(round(centre)) < 0:
                #print 'not used',order, pos_real,popt
                already_found[max(0,pos_real-1):min(len(already_found),pos_real+1)+1] = 1
                continue
            if already_found[int(round(centre))] == 1:
                #print 'already found',order, pos_real,popt
                continue
            already_found[max(0,int(centre-width*3)):min(len(yarr),int(centre+width*3))+1] = 1
            
            lines.append([order, centre, width, height ])
        return lines

def identify_emission_lines(params, im, im_short=None, im_badpx=None, im_short_badpx=None):
    """
    Identifies the lines in a spectrum by searching for the significant outliers in a polynomial fit to the data and subsequent fitting of Gaussian profiles to this positions
    :param im: 2d array with the extracted spectra
    :param im_short: 2d array with the extracted spectra in which the saturated lines of {im} should be identified
    :param im_badpx: The bad-pixel-mask for the extracted spectrum {im}. Used to identify the lines with saturated pixels
    :return lines: 2d array with one line for each identified line, sorted by order and amplitude of the line. For each line the following informaiton is given:
                    order, pixel, width of the line, and height of the line
    """
    ims = im.shape
    lines = []
    maxFWHM = max(1, int(round(estimate_width(im)*2.35482*1.5)))
    kwargs = []
    for order in range(ims[0]):
        kwargs.append([ order, im, im_short, im_badpx, im_short_badpx, maxFWHM ])
           
    if params['use_cores'] > 1 and multiprocessing.current_process().name == 'MainProcess':
       logger('Identify lines in the emission spectrum (multicore)')
       p = multiprocessing.Pool(params['use_cores'])
       lines_unsorted = p.map(identify_emission_lines_single_order, kwargs)
       p.terminate()
       p.join()
       for entry in lines_unsorted:
            lines += entry
    else:
       for kwarg in tqdm(kwargs, desc='Identify lines in the emission spectrum'):
            lines += identify_emission_lines_single_order(kwarg)

    lines = np.array(lines)
    # Remove the lines which are too wide/narrow
    good_values, pfit = sigmaclip(lines[:,2], nc=0, ll=3, lu=3)  #orders, sigma low, sigma high; x doesn't matter for nc=0
    lines = lines[good_values, :]
    
    return lines

def read_reference_catalog(params, filenames=None, wavelength_muliplier=None, arc_lines=None):
    """
    Reads the reference catalogue from a file and extracts the lines which should be used
    :param filename: text, filename of the catalogue file. The file needs to consist of 3 columns: wavelength of the line, 
                            strength/heigt of the line [can be empty or followed by text (as from the NIST database)], and name of the line
    :param wavelength_muliplier: float, if the resolution is not given in Angstrom, than it can be converted into Angstrom with this factor
    :param arc_lines: list of text, subsection of the arc lines, which should be extracted
    :return reference_catalog: 2d array with one entry for each line. Each entry contains the wavelength, the intensity of the line, and the index in the catalogue
    :return reference_names: list with same length as reference_catalog, name of each line
    """
    if filenames is None:                   filenames = params['reference_catalog']
    if type(filenames).__name__ == 'str':   filenames = [filenames]
    if wavelength_muliplier is None:        wavelength_muliplier = params['catalog_file_wavelength_muliplier']
    if arc_lines is None:                   arc_lines = params['use_catalog_lines']
    for i in range(len(arc_lines)):
        arc_lines[i] = arc_lines[i].replace(' ','')
    reference_lines_dict = dict(reference_catalog=[], reference_names=[])
    for ii, filename in enumerate(filenames):
        reference_catalog, reference_names = [], []
        refernce_files, refernce_files_full = find_file_in_allfolders(filename, [params['result_path']] + params['raw_data_paths'] + [os.path.realpath(__file__).rsplit(os.sep,1)[0]+os.sep])
        if len(refernce_files) == 0:
            logger('Error: file for the reference coordinates does not exist. Checked: {0}'.format(refernce_files_full) )
        for refernce_file in refernce_files:
            with open(refernce_file, 'r') as file:
                for line in file:
                    line = line[:-1].split('\t')
                    if len(line) < 3:
                        continue
                    if line[2].replace(' ','') not in arc_lines:        # if Ar IV
                        continue
                    try:
                        line[0] = float(line[0])*wavelength_muliplier   # wavelength in Angstrom
                    except:
                        logger('Warn: Wavelength {0} of line {1} in file {2} cannot be transformed to float'.format(line[0], line, refernce_file))
                        continue                                        # wavelength can't be transformed to float
                    for i in range(len(line[1])+1)[::-1]:               # Make intensity into a number, intensity can contain text at the end
                        try:
                            line[1] = float(line[1][:i])                # It's a number
                            break
                        except:
                            continue                                    # Try without the last character
                    if i == 0:
                        line[1] = 1.
                    if line[1] < 0:
                        logger('Warn: Line intensity cannot be smaller than 0.0, check line {0} in {1}'.format(line, refernce_file))
                        line[1] = 1.
                    if line[1] == 0:
                        line[1] = 0.1
                    reference_names.append(line[2])
                    line = line[0:2]
                    if reference_names[-1].find('Ar I') == 0:
                        line[1] *= 50
                    line.append(len(reference_catalog))
                    reference_catalog.append(line)   # wavelength in Angstrom, relative intensity of the line, index of line in reference_names
            if len(reference_catalog) > 0:
                break
        if len(reference_catalog) == 0:
            logger('Error: no reference lines found in {0} for the requested lines {1}'.format(refernce_files, arc_lines))
        reference_catalog = np.array(reference_catalog)
        arcs = reference_catalog.shape
        # Remove the faintest lines, if too many lines are in the catalogue
        """if arcs[0] > 100000:
            breakvalue = np.percentile(reference_catalog[:,1],min(90.,arcs[0]/120.))
            keep = np.logical_or(reference_catalog[:,1] >= breakvalue , reference_catalog[:,1] == 1)
            for i in range(arcs[0])[::-1]:
                if keep[i] == False:
                    del reference_names[i]                      # remove the names
            reference_catalog = reference_catalog[keep,:]       # remove the wavelengths, intensities
            logger('Info The faintest {0} of {1} entries in the arc reference file {2} will not be used '.format(arcs[0]-reference_catalog.shape[0], arcs[0], refernce_file ))"""
        reference_lines_dict['reference_catalog'].append(reference_catalog)
        reference_lines_dict['reference_names'].append(reference_names)
        
    return reference_lines_dict

def remove_orders_from_calimages(params, calimages, keep_orders):
    """
    In order to remove some orders from the extracted data
    """
    todo_full = ['sci_trace', 'cal_trace', 'flat_spec_norm']
    todo_start = ['wave_sol_dict_', 'wave_sols_']
    todo_end = []
    todo_not = ['wavesol2d']
        
    def remove_orders_from_array(arr, keep_orders, entry):
        if type(arr).__name__ == 'ndarray':
            # arr.shape[0] == keep_orders.shape[0] doesn't work, as keep_orders was created with where
            if np.max(keep_orders) < arr.shape[0]:
                if len(arr.shape) == 1:
                    arr = arr[keep_orders]
                elif len(arr.shape) >= 2:
                    arr = arr[keep_orders,:]
            elif len(arr.shape) >= 2:
                if np.max(keep_orders) < arr.shape[1]:
                    if len(arr.shape) == 2:
                        arr = arr[:,keep_orders]
                    elif len(arr.shape) >= 3:
                        arr = arr[:,keep_orders,:]
            else:
                logger(('Error: Something went wrong with the numbers of orders: highest order to keep: '+\
                        '{0}, wave_sol[{1}] has shape {2}. This should not have happened. '+\
                        'Please consult Ronny').format(np.max(keep_orders), entry, arr.shape))
        else:
            logger('Not yet ndarry type: {1} has {0}'.format(type(calimages[entry]).__name__, entry))
        return arr
    
    def circle_through(data_set, keep_orders, entry):
        if type(data_set).__name__ == 'ndarray':
            data_set = remove_orders_from_array(data_set, keep_orders, entry)
        elif type(data_set).__name__ == 'list':
            for ii in range(len(data_set)):
                 data_set[ii] = circle_through(data_set[ii], keep_orders, entry)
        elif type(data_set).__name__ in ['dict', 'OrderedDict']:
            for entry2 in data_set.keys():
                if entry2 not in todo_not:
                    data_set[entry2] = circle_through(data_set[entry2], keep_orders, entry+' -> '+entry2)
        elif type(data_set).__name__ in ['str', 'float']:
            'do nothing'
        else:
            logger('Warn: Unexpected data type {0} for {1}. Please consult Ronny'.format(type(data_set).__name__, entry))
            
        return data_set
    
    for entry in calimages.keys():
        process = False
        if entry in todo_full:      process = True
        else:
            for todo_part in todo_start:
                if entry.startswith(todo_part):     process = True
            for todo_part in todo_end:
                if entry.endswith(todo_part):       process = True
        if not process:
            continue
            
        calimages[entry] = circle_through(calimages[entry], keep_orders, entry)
    
    return calimages

def shift_wavelength_solution(params, aspectra, wave_sol_dict, reference_lines_dict, xlows, xhighs, obsdate_float, jd_midexp, sci_tr_poly, cal_tr_poly, objname, maxshift=20.0, in_shift=0, fib='cal', fine=False, im_head=dict()):
    """
    Determines the pixelshift between the current arc lines and the wavelength solution
    Two ways for a wavelength shift can happen:
        1: The CCD moves -> all ThAr lines are moved by the same pixel distance -> this was implemented from the beginning -> This implemented
        2: The light hits the grating at a different position and passes through different are of the lens -> Pixelshift depends on wavelength, pixel
            maybe replace by a cross correlation between arc spectra: x1, y1 from wavelength solution, x2, y2 from the current file -> y2'(x) = y1(x+dx)*a+b so that y2 and y2' match best
        3: Re-fit the wavelength solution each time -> replace the lines that have been identified in this aspectra
            if not enough lines available -> linear shift (method 1)
    :param aspectra: spectrum of the reference orders
    :param wave_sol_dict: dictionary with the wavelength solution and the catalog lines and other information
    #:param wavelength_solution: 2d array of floats, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    #:param wavelength_solution_arclines: 2d array of floats with one line for each order. Each order contains the wavelengths of the identified reference lines, 
                                        sorted by the brightest to faintest (in spectrum from which the solution is derived) and 0 to make into an array
    :param reference_catalog: 2d array of floats with one entry for each line. Each entry contains the wavelength and the intensity of the line
    :param reference_names: list of strings with same length as reference_catalog, name of each line
    :param xlows: list of floats, length same as number of orders, the lowest x pixel (wavelength direction) used in each order
    :param xhighs: list of floats, length same as number of orders, the highest x pixel (wavelength direction) used in each order
    :param in_shift: integer, gives an offset of where to expect the lines (positive -> means line in aspectra is to the right, compared to wavelength solution)
    :param fib: for which fiber is wavelength solution
    :return wavelength_solution: 2d array of floats, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    """
    sigma = 3.5     # before 20200708
    sigma = 2       # test on 20200708 -> no difference
    sigma = 3.0
    precision = 0.0001      # in pixel; 0.001 is a milli-pixel
    ratio_lines_identified = 0.15       # if more than ratio_lines_identified of the checked_arc_lines has been identified, then a sigma_clipping will be applied. If less than this number of lines remain after sigma clipping, then the assumption is, that the calibration fiber wasn't used and therefore no wavelength shift is applied
    #print(aspectra, wave_sol_dict, xlows, xhighs, obsdate_float, jd_midexp, sci_tr_poly, cal_tr_poly, objname, maxshift, in_shift, fib, fine)
    
    wavelength_solution = wave_sol_dict['wavesol']
    wavelength_solution_arclines = wave_sol_dict['reflines']
    reference_catalog, reference_names = reference_lines_dict['reference_catalog'][-1], reference_lines_dict['reference_names'][-1]
    #logger('Step: Finding the wavelength shift for this exposure')
    if np.max(wavelength_solution[:,-1]) < 1000:                 # pseudo solution, wavelength of the central pixel is < 1000 Angstrom
        return copy.deepcopy(wavelength_solution), 0., im_head               # No shift is necessary
    if np.nansum(np.abs(sci_tr_poly - cal_tr_poly)) == 0.0 and np.nansum(aspectra) == 0:         # science and calibration traces are at the same position and it's not the calibration spectrum
        # Not necessary anymore, as interpolate the solution later
        # wavelength_solution_shift, shift, im_head = shift_wavelength_solution_times(params, wavelength_solution, obsdate_float, jd_midexp, objname, im_head)
        return copy.deepcopy(wavelength_solution), 0., im_head
    
    if 'linestat' in wave_sol_dict.keys():
        FWHM = wave_sol_dict['linestat']            # per order: [median, std]
    else:
        FWHM = np.array([[3.5,0]]*wavelength_solution.shape[0])       # maybe seach for width instead?
    ratio_lines_identified = 0.15       # if more than ratio_lines_identified of the checked_arc_lines has been identified, then a sigma_clipping will be applied. If less than this number of lines remain after sigma clipping, then the assumption is, that the calibration fiber wasn't used and therefore no wavelength shift is applied
    in_shift_int = int(round(in_shift))
    #print 'input_shift', in_shift, in_shift_int
    
    aspectra = np.array(aspectra)
    ass = aspectra.shape
    shifts = []
    checked_arc_lines = 0
    xarr = np.arange(ass[1])                                    # array of the x-positions
    for order_index in range(wavelength_solution.shape[0]):
        if FWHM[order_index][0] == 0 or wavelength_solution_arclines[order_index][0] == 0.0:           # No emission line was identified, or no line from catalogue was identified
            continue
        warr = np.polyval(wavelength_solution[order_index,2:], xarr-wavelength_solution[order_index,1])     # Wavelength of each pixel in the order
        for arcline in wavelength_solution_arclines[order_index]:
            if arcline == 0.0:
                continue
            ref_line_index = np.argmin(np.abs( reference_catalog[:,0] - arcline ))
            if abs( reference_catalog[ref_line_index,0] - arcline ) > 0.05:           # check that it is in the reference catalog, allowing for uncertainty; catches also zeros used to fill the array
                continue                                                                # it's not in the reference catalog (is this necessary???)
            #diff = np.abs(warr-reference_catalog[ref_line_index,0])                     # Diff in wavelengths, Before 20200720
            diff = np.abs(warr-arcline)                                                 # Diff in wavelengths
            if min(diff) <= wavelength_solution[order_index,-2]*maxshift:               # Distance should be less than 2.5 px
                checked_arc_lines += 1
                pos = np.argmin(diff)                                               # Position of the line in the array
                range_arr = list(range( max(0,pos+in_shift_int-int(maxshift)), min(ass[1],pos+in_shift_int+int(maxshift)+1) ))     # Range to search for the arc line
                #print range_arr, xarr, aspectra.shape, order_index
                #get_timing('{0}'.format(arcline))
                popt = centroid_order(xarr[range_arr],aspectra[order_index,range_arr], pos+in_shift_int, FWHM[order_index,0]*3, significance=sigma, bordersub_fine=False)    #a,x0,sigma,b in a*np.exp(-(x-x0)**2/(2*sigma**2))+b
                #print popt, range_arr
                #get_timing('{0}'.format(popt))
                if popt[1] == 0 or popt[2] > FWHM[order_index,0]*1.5:
                    #if order_index == 50:
                    #    print 'not used',order_index, pos,popt, reference_catalog[ref_line_index,0]
                    #    plot_img_spec.plot_points([xarr[range_arr]],[aspectra[order_index,range_arr]],[str(pos)],'path',show=True, x_title='Pixel', y_title='Flux')
                    continue
                #if order_index == 50:
                #    print 'used',order_index, pos,popt, reference_catalog[ref_line_index,0]
                #    x = xarr[range_arr]
                #    plot_img_spec.plot_spectra(np.array([x, x]),np.array([aspectra[order_index,range_arr], oneD_gauss(x,popt)]),['data','fit'], ['savepaths'], True, [0.06,0.92,0.95,0.06, 1.0,1.01])
                pos_fine = invert_polyval(wavelength_solution[order_index,2:], xarr, arcline, poly0=wavelength_solution[order_index,1])
                """ before 20200811
                xarr_fine = np.arange(xarr[pos]-2, xarr[pos]+2+precision, precision)       # fine check of where the reference line is lockated 1/10000ths of a pixel, ?? interpolate
                warr_fine = np.polyval(wavelength_solution[order_index,2:], xarr_fine-wavelength_solution[order_index,1])       # Wavelength of the fine array
                #diff = np.abs(warr_fine-reference_catalog[ref_line_index,0])       # Before 20200720
                diff = np.abs(warr_fine-arcline)
                pos_fine = xarr_fine[np.argmin(diff)]"""
                if np.abs(popt[1] - (pos_fine+in_shift)) <= maxshift:                                   # Maximal 2.5 px shift
                    # Index of order, Index of reference line, wavelength of reference line, px of arc line, \n wavelength of arc line, \n height of arc line, sigma/width of gauss, original pixel position of reference line, residuals in px, position in px, wavelength from the wavelength solution
                    shifts.append([ order_index, ref_line_index, reference_catalog[ref_line_index,0], popt[1], \
                                    np.polyval(wavelength_solution[order_index,2:], popt[1]-wavelength_solution[order_index,1]), \
                                    popt[0], popt[2], pos_fine, popt[1]-pos_fine, popt[1], arcline])
    # shifts: order_index, ref_line_index, ref_wavelength, gauss_centre (3), wavelength at gauss_centre, gauss_height, gauss_width, position where catalogue line is located to 1/10000ths of a pixel
    shifts = np.array(shifts)
        
    if shifts.shape[0] >= ratio_lines_identified * checked_arc_lines and shifts.shape[0] != 0:    
        # Using the wavelenght shift, but this is wrong
        #good_values, pfit = sigmaclip(shifts[:,3]-shifts[:,2], nc=0, ll=3, lu=3, repeats=20)   # x doesn't matter for nc=0
        # Using the pixel shift
        good_values, pfit = sigmaclip(shifts[:,8], nc=0, ll=3, lu=3, repeats=20)    # x doesn't matter for nc=0
        shifts = shifts[good_values,:]
    
    if fine == True:
        return shifts
    elif maxshift > 5. and shifts.shape[0] >= ratio_lines_identified * checked_arc_lines and shifts.shape[0] != 0:
        shift_med = np.median( (shifts[:,8]) )      # Better: check within bins where the most matches are found
        shifts = shift_wavelength_solution(params, aspectra, wave_sol_dict, reference_lines_dict, 
                                           xlows, xhighs, obsdate_float, jd_midexp, sci_tr_poly, cal_tr_poly, objname, maxshift=5.5, in_shift=shift_med, fib=fib, fine=True)
    shift_med, shift_std, width_avg, width_std, textsource = in_shift, 0, 0, 0 , ''
    wavelength_solution_new, px_shift = copy.deepcopy(wavelength_solution), True
    
    if shifts.shape[0] >= ratio_lines_identified * checked_arc_lines and shifts.shape[0] != 0:                          # Only if enough lines have been detected
        # Log the results
        objname_txt = objname.rsplit(os.sep,1)[-1]      # If the name contains a / there will be a problem
        fname = '{0}_{1}.{2}'.format(params['logging_em_lines_gauss_width_form'].rsplit('.',1)[0], objname_txt, params['logging_em_lines_gauss_width_form'].rsplit('.',1)[-1])
        plot_wavelength_solution_width_emmission_lines(fname, ass, shifts, [0,3,6]) # order at 0, pixel at 3, gauss at 6
        fname = '{0}_{1}.{2}'.format(params['logging_arc_line_identification_residuals_hist'].rsplit('.',1)[0], objname_txt, params['logging_arc_line_identification_residuals_hist'].rsplit('.',1)[-1])
        plot_hist_residuals_wavesol(fname, shifts, [0,3,4,8] )
        fname = '{0}_{1}.{2}'.format(params['logging_em_lines_bisector'].rsplit('.',1)[0], objname_txt, params['logging_em_lines_bisector'].rsplit('.',1)[-1])
        bisector_measurements_emission_lines(fname, aspectra, shifts, [0,3,6])         # order, pixel, gauss width
        
        y_title='Residual (O-C) [Pixel]'
        width_avg, width_std = np.mean(shifts[:,6]), np.std(shifts[:,6], ddof=1)
        shift_med, shift_std = np.median( (shifts[:,8]) ), np.std( (shifts[:,8]), ddof=1)
        logtext = 'Median and Standard deviation of the offset of all {0} refitted lines is {1} +- {2} px.'.format(shifts.shape[0], round(shift_med,4), round(shift_std,4))
        textsource = ' Only {0} lines available, therefore used the median and standard deviation.'.format(shifts.shape[0])
        
        if shifts.shape[0] > 200:
            popts = []
            for bins in np.linspace(shifts.shape[0]/30, shifts.shape[0]/7, 11, dtype=int):
                hist_y, hist_x = np.histogram( shifts[:,8], bins )
                hist_x2 = (hist_x[1:]+hist_x[:-1])/2.           # center of the bin
                popt = centroid_order(hist_x2, hist_y, np.median(hist_x2), (hist_x2[-1]-hist_x2[0])/2., significance=3, bordersub_fine=False)    #a,x0,sigma,b in a*np.exp(-(x-x0)**2/(2*sigma**2))+b
                if popt[2] != 0:
                    popts.append(popt)
            if len(popts) > 3:
                popt_med = np.median(popts, axis=0)
                shift_med, shift_std = popt_med[1], popt_med[2]
                logtext += (' The Gaussian centre and width of the distribution of the refitted lines gives: {0} +- {1} px '+\
                            '(the values depend on the binning, standard deviation between {2} different bins is {3} px). ').format(round(shift_med,4), 
                                                                round(shift_std,4), len(popts), round_sig(np.std(popts, axis=0, ddof=1)[1],3))
                textsource = ' Used the Gaussian center and width of the distribution.'
                #im_head['HIERARCH HiFLEx D_SHIFT_GAUSS'] = (round(shift_med,5), 'Offset from Gauss distribution [px]')
                #im_head['HIERARCH HiFLEx D_SHIFT_G_ERR'] = (round(shift_std,5), 'Offset from Gauss distribution [px]')
        
        """ This doesn't make it better, in case of bifurcated fiber the scatter between the fiber increases by factor of 3
        if shifts.shape[0] >= 500:      # refitting the wavelength solution if at least 500 lines are available
            arc_lines_wavelength = shifts[:,[0, 9, 10, 0,0,0,0,0,0, 6]]         # order, pixel, wavelength, ..., width for fit_wavelengths_solution_2d and transform_wavelength_solution_2d_to_n1d
            cen_pxs = arc_lines_wavelength[:,1] * np.nan
            for order_index in range(wavelength_solution.shape[0]):
                cen_pxs[ (order_index == arc_lines_wavelength[:,0]) ] = wavelength_solution[order_index, 1]
            order_offset = wavelength_solution[0,0]
            polynom_order_trace, polynom_order_intertrace = params['polynom_order_traces'][-1], params['polynom_order_intertraces'][-1]
            poly2d_params, new_waves = fit_wavelengths_solution_2d(arc_lines_wavelength, cen_pxs, order_offset, polynom_order_trace, polynom_order_intertrace)
            shifts[:,8] = shifts[:,2] - new_waves               # new differences in Angstrom
            cen_pxs = wavelength_solution[:, 1]
            wavelength_solution_new, wavelength_solution_arclines, line_stats = transform_wavelength_solution_2d_to_n1d(ass[0], ass[1], 
                                    polynom_order_trace, polynom_order_intertrace, poly2d_params, order_offset, cen_pxs, arc_lines_wavelength)
            
            xarr = np.arange(ass[1])
            ori_px = 0.5*ass[1]
            oo = int(0.5*ass[0])
            cenwave = np.polyval(wavelength_solution[oo,2:], ori_px - wavelength_solution[oo,1])      # central wavelength
            shifts_px, shifts_wave = [], []
            for order in range(int(0.5*ass[0])-2, int(0.5*ass[0])+3):     #range(int(0.1*ass[0]),int(0.9*ass[0])):     # the shift is wavelength dependent
                waveold = np.polyval(wavelength_solution[order,2:], ori_px - wavelength_solution[order,1])      # central wavelength
                new_px = invert_polyval(wavelength_solution_new[order,2:], xarr, waveold, poly0=wavelength_solution_new[order,1])
                shifts_px.append(new_px - ori_px)
                shifts_wave.append(np.polyval(wavelength_solution_new[order,2:], ori_px - wavelength_solution_new[order,1]) - waveold)
            shift_med = np.median(shifts_px)
            textsource = ' Refitted the wavelength solution and using the central wavelength of the central orders to give a single shift.'
            logtext += ' The central wavelength of the central order ({0}) of the refitted the wavelength solution is shifted by {1} +- {2} px or {3} +- {4} Angstrom.'.format(round_sig(cenwave,7), 
                            round(shift_med,4), round(np.std(shifts_px,ddof=1),4), round(np.median(shifts_wave),4), round(np.std(shifts_wave,ddof=1),4))
            y_title='Residual (O-C) [Angstrom]'
            px_shift = False
            im_head['HIERARCH HiFLEx D_SHIFT_CENWAVE'] = (round(shift_med,5), 'Offset from spectrum center [px]')
            ## Revert to better value
            #shift_med = im_head['HIERARCH HiFLEx D_SHIFT_GAUSS']
        #"""
        
        # More plotting        
        x, w, y, l = [], [], [], []
        for order in range(ass[0]):
            inorder = ( shifts[:,0]==order )
            x.append(shifts[inorder,3])    # Pixel
            w.append(shifts[inorder,4])    # Wavelength
            y.append(shifts[inorder,8])    # Residuals
            l.append('{0}'.format(order))
        text = 'Residuals to the wavelength solution: used {0} out of {1} lines'.format(shifts.shape[0], checked_arc_lines)
        fname_px   =   '{0}_{1}_px.{2}'.format(params['logging_arc_line_identification_residuals'].rsplit('.',1)[0], objname_txt, params['logging_arc_line_identification_residuals'].rsplit('.',1)[-1])
        fname_wave = '{0}_{1}_wave.{2}'.format(params['logging_arc_line_identification_residuals'].rsplit('.',1)[0], objname_txt, params['logging_arc_line_identification_residuals'].rsplit('.',1)[-1])
        plot_img_spec.plot_points(x,y,l,[fname_px],   show=False, title=text, x_title='Pixel', 
                                  y_title=y_title, marker=['o','s','*','P','^','v','>','<','x'])
        plot_img_spec.plot_points(w,y,l,[fname_wave], show=False, title=text, x_title='Wavelength [Angstrom]', 
                                  y_title=y_title, marker=['o','s','*','P','^','v','>','<','x'])
    
        logger('Info: '+logtext+' The resulting shift between {2}-fiber and the emission line spectrum in {3} of {0} +- {1} px will be used'.format(round(shift_med,4), 
                        round(shift_std,4), fib, objname))
        # Save the shift for later use
        if params['extract_wavecal']:            # This values are not correct, as later the shift has to be applied
            add_text_to_file('{0}\t{1}\t{2}\t{3}'.format(jd_midexp, round(shift_med-in_shift,4), round(shift_std,4), fib), params['master_wavelengths_shift_filename'] )    # Not the shift between the two solutions
        params['force_DT_SHIFT'] = False
    else:       # not enough lines identified
        params['force_DT_SHIFT'] = True         # If wavelength solution only in calibration fiber, but not in this file use stored information
    
    logtext = 'The overal shift between wavelength solution and the current calibration spectrum'
    if not params['extract_wavecal']:         # science and calibration traces are at the same position and it's not the calibration spectrum
        # In case of pixel shift available -> linear interpolation of pixel shift
        if fib == 'cal':
            wavelength_solution_dummy, shift_stored, im_head = shift_wavelength_solution_times(params, wavelength_solution_new, obsdate_float, jd_midexp, objname, im_head)
            shift_med += shift_stored
            logtext = 'The overal shift to be applied to the wavelength solution in the science fiber for the current calibration spectrum'
    #if params['two_solutions']:
    #    text = '(from {0}) '.format( params['wavelength_solution_type'].replace('cal', 'calibration').replace('sci','science') )
    #else:
    #    text = ''
    d_shift_kms = round(shift_med*np.median(wavelength_solution[:,-2]/wavelength_solution[:,-1])*Constants.c/1000.,4)
    logger(('Info: '+logtext+' of file {9} at JD {10} is {0} +- {1} px ({8} km/s).{11} {2} reference lines have been used, {7} reference lines have been tested.'+\
            ' The calibration lines have a Gaussian width of {3} +- {4} px, which corresponds to a FWHM of {5} +- {6} px.').format(round(shift_med,4), round(shift_std,4), 
                        shifts.shape[0], round(width_avg,3), round(width_std,3), round(width_avg*2.35482,3), round(width_std*2.35482,3), checked_arc_lines,
                        d_shift_kms, objname, round(jd_midexp,5), textsource ))
    im_head['HIERARCH HiFLEx D_SHIFT'] = (round(shift_med,5), 'Offset in dispersion direction [px]')        # Includes DT_SHIFT
    im_head['HIERARCH HiFLEx D_SHIFT_ERR'] = (round(shift_std,5), 'Uncertainty of the offset [px]')
    im_head['HIERARCH HiFLEx D_SHIFT_NUMBER_LINES'] = (shifts.shape[0], 'out of {0} calibration lines'.format(checked_arc_lines))
    im_head['HIERARCH HiFLEx D_WIDTH'] = (round(width_avg,2), 'Gaussian width of the calibration lines [px]')
    im_head['HIERARCH HiFLEx D_WIDTH_ERR'] = (round(width_std,5), 'Uncertainty of the width [px]')
    im_head['HIERARCH HiFLEx D_SHIFT_KMS'] = (d_shift_kms, 'Offset in dispersion direction [km/s]')
    if len(shifts) >= ratio_lines_identified * checked_arc_lines and len(shifts) != 0:                          # Statistics only if enough lines were detected
        statistics_arc_reference_lines(shifts, [0,1,6,2], reference_names, wavelength_solution, xlows, xhighs, show=False)
    
    if px_shift:
        #wavelength_solution_new[:,1] -= shift_avg                       # shift the central pixel, - sign is right, tested before 19/9/2018
        wavelength_solution_new[:,1] += shift_med                       # shift the central pixel, + sign is right, tested on 19/9/2018
    
    if False:           # The changes below are not necessary, shift_med/shift_avg will contain the input shift
        shift_avg -= in_shift           # To separate the input shift from the rest
        logger('Info: Corrected for input shift. The shift of the currect spectrum is {0}'.format(round(shift_avg,4) ))
    #print 'return_shift', shift_avg
    return wavelength_solution_new, shift_med, im_head
 
def read_textfile_remove_double_entries(fname, dataformats, delimiter='\t', replaces=['\n'], equal_indexes=[0], warn='Warn: no data in file {0}'):
    """
    Reads all the information from params['master_wavelengths_shift_filename'], removes the older entries
    :param fname: string, filename
    :param dataformat: list of data types, e.g [float, int, str], must contain
    :param delimiter: string, how are coloumns are separated
    :param replaces: sting, what data should be replaced
    :equal_indexes: list of integers, what indexes (coloumns) should be compared 
    :return data: 2d array of type string, with each entry containing the last version of the line in fname
    """
    if len(dataformats) < max(equal_indexes)+1:
        logger('Error: Coding error when calling read_textfile_remove_double_entries(): equal_indexes ({0}) contains indexes higher than the length of the dataformats ({1})'.format(equal_indexes, dataformats))
    data = read_text_file(fname, no_empty_lines=True)
    data = convert_readfile(data, dataformats, delimiter=delimiter, replaces=['\n'])       # jd_midexp, shift_avg, shift_std, can contain duplicate jd_midexp (last one is the reliable one)
    if len(data) == 0:
        if len(warn) > 3:
            logger(warn.format(fname))
        data = np.array([])
    else:
        data = np.array(data)                           # Becomes an array of strings
        if len(equal_indexes) > 0:
            # Remove the double entries, e.g. if better data was added later
            data_no_doubles = []
            for i in range(data.shape[0])[::-1]:
                index = equal_indexes[0]
                equal_arr = ( data[i:,index] == data[i,index] ) 
                for index in equal_indexes[1:]:
                    equal_arr &= ( data[i:,index] == data[i,index] ) 
                if np.sum( equal_arr ) == 1:    # Find only exactly this entry
                    data_no_doubles.append(data[i,:])
            data = np.array(data_no_doubles)                     # Only contains the latest entries, newest entry first
    
    return data
   
def shift_wavelength_solution_times(params, wavelength_solution, obsdate_float, jd_midexp, objname, im_head):
    """
    In case no calibration spectrum was taken at the same time as the science spectra use the stored information to aply a shift in the lines
    :param wavelength_solution: 2d array of floats, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    
    :return wavelength_solution: 2d array of floats, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    :return shift_avg: Float, shift from science to calibration fiber
    """
    warn = 'Warn: No pixel-shift for the wavelength solution is available. Please re-run prepare_file_list.py and asign "w" to the emission line spectra.'
    all_shifts_str = read_textfile_remove_double_entries(params['master_wavelengths_shift_filename'], [float, float, float, str], equal_indexes=[0,3], warn=warn) # Find only exactly this entry with same jd_midexp and same fiber
    # all_shifts_str = 2d array of type string, with each entry containing JD of the mid exposure, Shift to the wavelength solution, Standard deviation of the shift, fiber type
    if len(all_shifts_str) == 0:
        logger('Warn: No pixel-shift for the wavelength solution is available. Please re-run file_assignment.py and asign "ws" to the emission line spectra.')
        return copy.deepcopy(wavelength_solution), 0., im_head               # No shift is possible
    
    # select the right data, depending on params['two_solutions']
    hours = 1.
    """ before 20200813
    if not params['two_solutions']:                                 # single soution: use all lines
        data = all_shifts_str[:,0:3].astype(float)
        data_bef = data[ (data[:,0] <= jd_midexp), : ]
        data_aft = data[ (data[:,0] >= jd_midexp), : ]
        if data_bef.shape[0] > 1:
            goodtime = data_bef[:,0] >= np.max(data_bef[:,0]) - hours / 24.      # only use data taken within 1 hour of the closest data for the fit
            data_bef = data_bef[goodtime ,:]
        if data_aft.shape[0] > 1:
            goodtime = data_aft[:,0] <= np.min(data_aft[:,0]) + hours / 24.      # only use data taken within 1 hour of the closest data for the fit
            data_aft = data_aft[goodtime ,:]
        if data_bef.shape[0] == 0:
            shift_avg = np.median(data_aft[:,1])
            text = 'Median of the datapoints in {0} after the exposure ( {1} to {2} ).'.format( params['master_wavelengths_shift_filename'], np.min(data_aft[:,0]), np.max(data_aft[:,0]) )
        elif data_aft.shape[0] == 0:
            shift_avg = np.median(data_bef[:,1])
            text = 'Median of the datapoints in {0} before the exposure ( {1} to {2} ).'.format( params['master_wavelengths_shift_filename'], np.min(data_bef[:,0]), np.max(data_bef[:,0]) )
        else:
            data = np.vstack(( data_bef, data_aft ))
            p = np.polyfit(data[:,0], data[:,1], 1)
            shift_avg = np.polyval(p, jd_midexp)
            text = 'Linear fit of the datapoints in {0} around the exposure ( {1} to {2} ).'.format( params['master_wavelengths_shift_filename'], np.min(data[:,0]), np.max(data[:,0]) )
    else:                                                           # two wavelength solutions -> compare the two"""
    if True:
        all_shifts_cal = all_shifts_str[all_shifts_str[:,3]=='cal',0:3].astype(float)       # jd_midexp, shift_avg, shift_std, fiber
        all_shifts_sci = all_shifts_str[all_shifts_str[:,3]=='sci',0:3].astype(float)
        modes = ['cal','sci']
        if len(all_shifts_cal)*len(all_shifts_sci) == 0:
            logger('Warn: Offset between the fibers cannot be determined as the offset is only available for one fiber. '+\
                    'Please re-run file_assignment.py and asign "wc" and "ws" to the emission line spectra of the calibration and science fiber, respectively.')
            #return copy.deepcopy(wavelength_solution), 0., im_head               # No shift is possible
            if len(all_shifts_sci) == 0 and params['force_DT_SHIFT']:           # Do it only if no lines in the cal are available
                modes = ['cal']
            else:
                return copy.deepcopy(wavelength_solution), 0., im_head               # No shift is possible
        data = dict()
        if 'cal' in modes:
            data['cal_bef'] = all_shifts_cal[ (all_shifts_cal[:,0] <= jd_midexp), : ]
            data['cal_aft'] = all_shifts_cal[ (all_shifts_cal[:,0] >= jd_midexp), : ]
        if 'sci' in modes:
            data['sci_bef'] = all_shifts_sci[ (all_shifts_sci[:,0] <= jd_midexp), : ]
            data['sci_aft'] = all_shifts_sci[ (all_shifts_sci[:,0] >= jd_midexp), : ]
        for entry in data.keys():
            if data[entry].shape[0] > 1:
                if entry.endswith('_bef'):
                    goodtime = data[entry][:,0] >= np.max(data[entry][:,0]) - hours / 24.      # only use data taken within 1 hour of the closest data for the fit
                if entry.endswith('_aft'):
                    goodtime = data[entry][:,0] <= np.min(data[entry][:,0]) + hours / 24.      # only use data taken within 1 hour of the closest data for the fit
                data[entry] = data[entry][goodtime ,:]
        for mode in modes:
            if data[mode+'_bef'].shape[0] == 0:
                data[mode+'_jd'] = np.median(data[mode+'_aft'][:,1])
                data[mode+'_text'] = 'Median after the exposure: {0} - {1} '.format( round(np.min(data[mode+'_aft'][:,0]),4), round(np.max(data[mode+'_aft'][:,0]),4) )
            elif data[mode+'_aft'].shape[0] == 0:
                data[mode+'_jd'] = np.median(data[mode+'_bef'][:,1])
                data[mode+'_text'] = 'Median before the exposure: {0} - {1} '.format( round(np.min(data[mode+'_bef'][:,0]),4), round(np.max(data[mode+'_bef'][:,0]),4) )
            else:
                data[mode] = np.vstack(( data[mode+'_bef'], data[mode+'_aft'] ))
                data[mode+'_p'] = np.polyfit(data[mode][:,0], data[mode][:,1], 1)
                data[mode+'_jd'] = np.polyval(data[mode+'_p'], jd_midexp)
                data[mode+'_text'] = 'Linear fit around the exposure: {0} - {1} '.format( round(np.min(data[mode][:,0]),4), round(np.max(data[mode][:,0]),4) )
        if 'cal' in modes and 'sci' in modes:
            shift_avg = (data['sci_jd'] - data['cal_jd'])                 # Thought about no - on 20200612
            text = 'It is the difference between science ({1}) and calibration ({2}) datapoints in {0}.'.format( params['master_wavelengths_shift_filename'], data['sci_text'], data['cal_text'] )
        else:       # only cal data
            shift_avg = 0 - data['cal_jd']
            text = 'It is the difference between calibration ({1}) datapoints in {0}.'.format( params['master_wavelengths_shift_filename'], data['cal_text'] )
        
    logger('Info: The wavelength solutions of both fibers drifted to a relative offset at the current file {0} (center of exposure is {1}, JD = {2}) of {3} px ({4} km/s). {5}'.format(\
                        objname, datetime.datetime.utcfromtimestamp(obsdate_float).strftime('%Y-%m-%d %H:%M:%S'), round(jd_midexp,5),
                        round(shift_avg,4), round(shift_avg*np.median(wavelength_solution[:,-2]/wavelength_solution[:,-1])*Constants.c/1000.,4),
                        text ) )
    im_head['HIERARCH HiFLEx DT_SHIFT'] = (round(shift_avg,4), 'Applied shift in dispersion direction [px]')
    #im_head['HIERARCH HiFLEx DT_SHIFT1'] = (round(shifts[0,1],4), 'Shift before @ {0}'.format(round(shifts[0,0],5) ))
    #im_head['HIERARCH HiFLEx DT_SHIFT2'] = (round(shifts[-1,1],4), 'Shift after  @ {0}'.format(round(shifts[-1,0],5) ))
    im_head['HIERARCH HiFLEx DT_SHIFT_KMS'] = (round(shift_avg*np.median(wavelength_solution[:,-2]/wavelength_solution[:,-1])*Constants.c/1000.,4), 'Applied shift in dispersion d. [km/s]')

    wavelength_solution_new = copy.deepcopy(wavelength_solution)
    #wavelength_solution_new[:,1] -= shift_avg                       # shift the central pixel, - sign is right, tested before 19/9/2018
    wavelength_solution_new[:,1] += shift_avg                       # shift the central pixel, + sign is right, tested on 19/9/2018
    
    return wavelength_solution_new, shift_avg, im_head

def create_wavelengths_from_solution(params, wavelength_solution, spectra, wave_sols_sci=None, jd_mid=0, shift=0):
    """
    Converts the wavelength solution into a 2d array with the wavelength for each pixel and order
    If the wavelength solutions of science fiber exist, use this value
    :param wavelength_solution: 2d array of floats, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    :param spectra: 2d array of floats, spectrum, only needed for the shape of the wavelengths array
    :param wave_sols_sci: list of wavelength solutions, each entry contains a list of entries: array of the wavelength solution, array of the reference lines, jd, name
    :param jd_mid: float, JD of the Midexposure to create the wavelength array for
    :param shift: float, input shift. By how many pixel need the solution be shifted
    :return wavelengths: 2d array of floats, same dimensions as spectra, contains the wavelength for each pixel
    :return wavelength_solution_recreated: only when wave_sols_sci is given, returns a wavelength_solution
    """
    if wave_sols_sci is None:           # Original, no interpolation of different solutions
        wavelength_solution_lst = [copy.deepcopy(wavelength_solution)]     # copy, otherwise it will be changed in wavesol[:,1] += 
        jd_weight_shift = np.array([[0,1,0]])
    else:
        warn = 'Warn: No pixel-shift for the wavelength solution is available. Please re-run prepare_file_list.py and asign "w" to the emission line spectra.'
        all_shifts_str = read_textfile_remove_double_entries(params['master_wavelengths_shift_filename'], [float, float, float, str], equal_indexes=[0,3], warn=warn)
        all_shifts_sci = all_shifts_str[all_shifts_str[:,3]=='sci',0:3].astype(float)       # jd_midexp, shift_avg, shift_std (, fiber)
        # Find the values closest to the opservation (before and after)
        jd_weight_shift = np.zeros((len(wave_sols_sci), 3))
        wavelength_solution_lst, name = [], []
        for ii, entry in enumerate(wave_sols_sci):
            wavelength_solution_lst.append(entry[0])
            jd_weight_shift[ii,0] = entry[2]
            name.append(entry[3])
            if all_shifts_sci.shape[0] > 0:
                diff = np.abs( all_shifts_sci[:,0] - entry[2] )
                if np.min(diff) < 1E-4:
                    jd_weight_shift[ii,2] = all_shifts_sci[np.argmin(diff), 1]              # only necessary for bifurcated fibers, for single fibers this will be corrected by the average value
        name = np.array(name)   # defining beforehand might cut the string
        jd_weight_shift[:,1] = jd_mid-jd_weight_shift[:,0]
        closestp = ( jd_weight_shift[:,1] >= 0 )      # after
        closestn = ( jd_weight_shift[:,1] <= 0 )      # before
        if not closestp.any():          # no solution after -> use the ones from before as dummy
            closestp = closestn
        elif not closestn.any():        # no solution before -> use the ones from after as dummy
            closestn = closestp
        jd_weight_shift[:,1] = np.abs(jd_weight_shift[:,1])
        closest = ( ( ( jd_weight_shift[:,1] <= np.min(jd_weight_shift[closestp,1] + 1./24) ) & closestp ) | ( ( jd_weight_shift[:,1] <= np.min(jd_weight_shift[closestn,1] + 1./24) ) & closestn ) )
        jd_weight_shift = jd_weight_shift[closest,:]
        jd_weight_shift[:,1] = 1./(jd_weight_shift[:,1]+5./24/60)                 # everything within 5 minutes is good, then decrease the weight
        wavelength_solution_lst = np.array(wavelength_solution_lst)[closest]
        name = name[closest]
        result = []
        if not params['two_solutions']:
            avg_waveshift = np.average( jd_weight_shift[:,2], weights=jd_weight_shift[:,1] )
            jd_weight_shift[:,2] -= avg_waveshift
        for ii in range(jd_weight_shift.shape[0]):
            result.append(list(jd_weight_shift[ii,0:3])+ [name[ii] ])
        printarrayformat = ['\t\t%4.4f', '%3.1f', '%4.4f', '%s']
        textshift = ''
        if shift != 0:
            textshift = ' (The input shift was {0} px)'.format(round(shift,4))
        logger('Info: Used the following wavelength solutions of the science fibers and their corresponding weights to'+\
               ' calculate the wavelength array for JD (midexposure) of {0} days{1}: '.format(round(jd_mid,4), textshift), show=False)
        logger('JD (mid)\tweight\tshift\tname', printarrayformat=printarrayformat, printarray=result, show=False)
        order_offset = wavelength_solution[0,0]
        
    wavelengths = np.zeros(( len(wavelength_solution_lst), spectra.shape[0], spectra.shape[1] ))
    xarr = np.arange(np.array(spectra).shape[1])
    for jj, wavesol in enumerate(wavelength_solution_lst):
        wavesol[:,1] += shift-jd_weight_shift[jj,2]      # Add shift (e.g. between Sci and Cal fiber and between this ThAr and Cal solution
        for ii, wls in enumerate(wavesol):
            wavelengths[jj,ii,:] = np.polyval(wls[2:], xarr - wls[1])
        #print( wavelengths[jj, 34, int(spectra.shape[1]/2)] )
    wavelengths = np.average( wavelengths, axis=0, weights=jd_weight_shift[:,1] )
    
    if wave_sols_sci is None:
        return wavelengths, wavelength_solution
    else:
        wavelength_solution_recreated = create_wavelength_solution_from_array(params, wavelengths, order_offset)
        
        return wavelengths, wavelength_solution_recreated

def create_wavelength_solution_from_array(params, wavelengths, order_offset):
    """
    Creates a solition from the wavelength array
    """
    polynom_order_trace = params['polynom_order_traces'][-1]
    wavelength_solution = np.zeros((wavelengths.shape[0], 2+polynom_order_trace+1))
    xarr = np.arange(0, wavelengths.shape[1])
    cen_px = int(wavelengths.shape[1]/2.)
    wavelength_solution[:,0] = order_offset+np.arange(wavelengths.shape[0])
    wavelength_solution[:,1] = cen_px
    for order in range(wavelengths.shape[0]):
        polyfit = np.polyfit(xarr-cen_px, wavelengths[order,:], polynom_order_trace)      #lambda from px
        wavelength_solution[order,2:] = polyfit
    
    return wavelength_solution

def find_px_from_wave(wavelength, wl_array, px_array):
    """
    Find the closest pixel for a given wavelength, something like a procedure that does invert_polyval()
    :param wavelength: float, a single wavelength 
    :param wl_array: 1d-array or list of floats, containing the wavelengths
    :param px_array: 1d-array or list of floats, containing the pixels
    :return px_array[index]: pixel closest to the wavelength
    :return diff: difference between the given wavelength and the wavelength of the pixel
    """
    wl_array = np.array(wl_array)
    index = np.argmin(np.abs( wl_array - wavelength ))
    diff =  wl_array[index] - wavelength
    
    return px_array[index], diff

def invert_polyval(poly, xx, yy, poly0=0):
    """
    :param poly: parameters of a polynomial
    :param xx: 1d array
    :param yy: float, find the x-value to this value
    :param poly0: 
    """
    yfit = np.polyval(poly, xx-poly0)
    index = np.argmin(np.abs( yy - yfit ))
    xx_sub = np.linspace(xx[max(0,index-3)], xx[min(xx.shape[0]-1,index+3)], 100)        # just get 100 values
    yfit_sub = np.polyval(poly, xx_sub-poly0)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', np.RankWarning)
        poly_inv = np.polyfit(yfit_sub, xx_sub, poly.shape[0]-1)        # fit x(y), might cause rank warning
    xx_inv = np.polyval(poly_inv, yy)

    return xx_inv

def find_adjust_trace_orders(params, im_sflat, im_sflat_head):
    """
    Procedure to search and adjust the traces. This was originally a standalone python script.
    :param params: Dictionary with all the parameters.
    :param im_sflat: 2d numpy array of floats, image of the science fiber ilumintated with a continuum source
    :param im_sflat_head: Header of im_sflat
    :return polyfits: 2d numpy array of floats, length same as number of orders, each line contains the polynomial values
    :return xlows: 1d array of integers, length same as number of orders, minimum position of the trace in dispersion direction
    :return xhighs: 1d array of integers, length same as number of orders, maximum position of the trace in dispersion direction
    :return widths: 2d array of floats, length same as number of orders, width of the trace in cross-dispersion direction. 
                    Each line contains the giving the width to the left, to the right, and the Gaussian width
    """
    logger('Step: Tracing the orders')
    # Find the traces in a binned file
    if os.path.isfile(params['logging_traces_binned']) == True:
        logger('Info: Using existing order file: {0}'.format(params['logging_traces_binned']))
        polyfits, xlows, xhighs, dummy = read_fits_width(params['logging_traces_binned'])
    else:
        # make flat smaller: combine px in spatial direction
        params['bin_search_apertures'] = adjust_binning_UI(im_sflat, params['bin_search_apertures'], userinput=params['GUI'])
        params['binx'], params['biny'] = params['bin_search_apertures']             # Required for find_trace_orders        
        sim_sflat, dummy, dummy = bin_im(im_sflat, params['bin_search_apertures'], method='mean' )
        save_im_fits(params, sim_sflat, im_sflat_head, params['logging_trace1_binned'])
     
        # Search for orders in the small image
        polyfits, xlows, xhighs = find_trace_orders(params, sim_sflat, im_sflat.shape)
        # save parameters of the polynoms into a fitsfile (from Neil)
        save_fits_width(polyfits, xlows, xhighs, [], params['logging_traces_binned'])
        plot_traces_over_image(im_sflat, params['logging_traces_im_binned'], polyfits, xlows, xhighs)
    if params['GUI']:
        logger('Step: Allowing user to remove orders')
        fmask, dummy, dummy = remove_adjust_orders_UI( im_sflat, polyfits, xlows, xhighs, userinput=params['GUI'], do_rm=True)
        polyfits, xlows, xhighs = polyfits[fmask], xlows[fmask], xhighs[fmask]
        plot_traces_over_image(im_sflat, params['logging_traces_im_binned'], polyfits, xlows, xhighs)
    # retrace orders in the original image to finetune the orders
    params['bin_adjust_apertures'] = adjust_binning_UI(im_sflat, params['bin_adjust_apertures'], userinput=params['GUI'])
    params['binx'], params['biny'] = params['bin_adjust_apertures']
    sim_sflat, dummy, dummy = bin_im(im_sflat, params['bin_adjust_apertures'])        # Not saved
    polyfits, xlows, xhighs, widths = adjust_trace_orders(params, sim_sflat, im_sflat, polyfits, xlows, xhighs)
    if params['GUI']:
        logger('Step: Allowing user to remove orders and to adjust the extraction width')
        fmask, polyfits, widths = remove_adjust_orders_UI( im_sflat, polyfits, xlows, xhighs, widths, userinput=params['GUI'], do_rm=True, do_adj=True)
        polyfits, xlows, xhighs, widths = polyfits[fmask], xlows[fmask], xhighs[fmask], widths[fmask]      

    return polyfits, xlows, xhighs, widths

def prepare_measure_background_noise(params, im, sci_tr_poly, xlows, xhighs, widths, cal_tr_poly, axlows, axhighs, awidths):
    """
    Removes the area with traces from an images and then measures the noise (std) in the remaining image
    Part of it is the same as read_file_calibration, if localbackground is not part of the calibrations
    :return bck_noise_std: median of the noise
    :return bck_noise_var: variance of the noise in the image (std of noise in different areas)
    """
    bck_px_sci = find_bck_px(im, sci_tr_poly, xlows, xhighs, widths, params['background_width_multiplier'][0])
    bck_px_cal = find_bck_px(im, cal_tr_poly, axlows, axhighs, awidths, params['background_width_multiplier'][1])
    bck_px = bck_px_sci * bck_px_cal
    bad_values = ( im*bck_px > np.percentile(im[bck_px==1],95) )
    bck_px[ bck_px==0 ] = np.nan
    bck_px[bad_values] = np.nan
    bck_noise_std, bck_noise_var = measure_background_noise(im * bck_px)
    return bck_noise_std, bck_noise_var

def combine_photonnoise_readnoise(spectra, bck_noise):
    """
    Combine the photon noise and background noise. Used for the extracted spectra
    Better would be flux/gain under the sqrt, but should be close to 1 anyway
    :param background_noise: background noise over the extracted area
    """
    espectra = []
    pos_spectra = copy.copy(spectra)
    pos_spectra[np.isnan(pos_spectra)] = 0          # otherwise the line below creates a invalid value error
    pos_spectra[pos_spectra < 0] = 0                # avoid errors
    for order in range(spectra.shape[0]):
        espectra.append( np.sqrt( pos_spectra[order,:] + bck_noise[order,:]**2 ) )
    espectra = np.array(espectra)
    espectra[np.isnan(spectra)] = np.nan            # revert replacing NaNs with 0
    return espectra

def correct_blaze(flat_spec_norm, minflux=0.1):
    """
    Do the blaze correction (e.g. divide through the normalised flat)
    :param spectra: 2d array of floats with the ectracted spectrum
    :param flat_spec_norm: 3d array of floats, extracted and normalised flat as flat_spec_norm[1]
    :param minflux: float, exclude the data with less flux than this value relative to the medium flux
    """
    # Fitting the blaze with some function? or smoothing it over serveral pixels?
    flat_spec = copy.copy(flat_spec_norm)
    fss = flat_spec.shape
    low_flux = np.zeros(fss, dtype=bool)
    # Use global value for minflux
    #with np.errstate(invalid='ignore'):
    #    low_flux = ( flat_spec < minflux )          # only use pixels with enough flux, (e.g. flat_spec_norm[1] < 0.1 means flux needs to be at least 10% of median flux)
    # Replace single values which might have a bit more flux than the minflux, but are surrounded by pixels with not enough flux
    for order in range(fss[0]):                                                 # each order
        # Find the low flux in each order
        with np.errstate(invalid='ignore'):
            low_flux[order,:] = ( flat_spec[order,:] < minflux*np.nanmax(flat_spec[order,:]) )          # only use pixels with enough flux, (e.g. flat_spec_norm[1] < 0.1 means flux needs to be at least 10% of median flux)
        for i in range(1, int(fss[1]/10)):
            with np.errstate(invalid='ignore'):
                nonnan = np.sum(~np.isnan(flat_spec[order,:i+1]))
            if low_flux[order,i] and np.sum(low_flux[order,:i+1]) < nonnan:
                #print('a', order, i)
                if nonnan < 20 or np.sum(low_flux[order,:i+1]) > nonnan*0.9:     # 90% data points have low flux or beginning of array
                    #print('c', order, i)
                    low_flux[order,:i+1] = True
            nonnan = np.sum(~np.isnan(flat_spec[order,-i-1:]))
            if low_flux[order,-i-1] and np.sum(low_flux[order,-i-1:]) < nonnan:
                #print('b', order, i)
                if nonnan < 20 or np.sum(low_flux[order,-i-1:]) > nonnan*0.9:     # 90% data points have low flux or end of array
                    #print('d', order, i)
                    low_flux[order,-i-1:] = True
    flat_spec[low_flux] = np.nan
    return flat_spec

def clip_noise(spectra, maxnoise=10, noisedataset=4, correctdatasets=[5,6]):
    """
    :param spectra: 3d array of floats
    :param maxnoise: float, exclude the data which has higher noise than this value, relative to the median noise
    :param dataset: integer, which data to use for noise clipping. (dataset=4: error of the flat corrected spectrum (residuals to a polynomial) )
    :param correctdatasets: list/array of integers, which data to correct
    """
    spectra = copy.deepcopy(spectra)
    with np.errstate(invalid='ignore'):
        high_noise = (spectra[noisedataset] > maxnoise * np.nanmedian(spectra[noisedataset]) )
    # Smooth the results to have fixed borders
    for i in correctdatasets:
        spectra[i, high_noise] = np.nan
    #print np.nanmedian(spectra[noisedataset]), np.sum(high_noise), np.sum(high_noise, axis=1).shape, np.sum(high_noise, axis=1)
    return spectra

def get_julian_datetime(date):      # from https://stackoverflow.com/questions/31142181/calculating-julian-date-in-python
    """
    Convert a datetime object into julian float.
    Args:
        date: datetime-object of date in question

    Returns: float - Julian calculated datetime.
    Raises: 
        TypeError : Incorrect parameter type
        ValueError: Date out of range of equation
    """

    # Ensure correct format
    if not isinstance(date, datetime.datetime):
        raise TypeError('Invalid type for parameter "date" - expecting datetime')
    elif date.year < 1801 or date.year > 2099:
        raise ValueError('Datetime must be between year 1801 and 2099')

    # Perform the calculation
    julian_datetime =   367 * date.year \
                      - int((7 * (date.year + int((date.month + 9) / 12.0))) / 4.0) \
                      + int((275 * date.month) / 9.0) \
                      + date.day + 1721013.5 \
                      + (date.hour + date.minute / 60.0 + date.second / math.pow(60,2)) / 24.0 \
                      - 0.5 * math.copysign(1, 100 * date.year + date.month - 190002.5) + 0.5

    return julian_datetime

def get_obsdate(params, im_head):
    """
    Get the observation date and time using the header and exposure time
    !!! Necessary improvement: check what header keywords are already available to not overwrite information created by create_image_general()
    :return obsdate_midexp: datetime.datetime object, containing the the (light collection) center of the observation (in UTC)
    :return obsdate_mid_float: float, unix timestamp of the mid observation time (in UTC)
    :return jd_midexp: float, Julian data of mid-exposure
    """
    obsformats = ['%Y-%m-%dT%H:%M:%S.%f','%Y-%m-%dT%H:%M:%S']           # Put into the parameters?
    exp_fraction_keys = params['raw_data_mid_exposure_keys']
    # Get the obsdate
    obsdate = -1
    for header_key in [ params['raw_data_dateobs_keyword'] ]:       # 'HiFLEx MID_'+params['raw_data_dateobs_keyword'] as it's mid exposure time
        if header_key in im_head.keys():
            obsdate = im_head[header_key]
    if obsdate == -1:
        logger('Warn: Cannot find the raw_data_dateobs_keyword = {0} in the header. Assuming 2000-01-01T00:00:00. This will cause wrong barycentric corrections.'.format(params['raw_data_dateobs_keyword'] ))
        obsdate = '2000-01-01T00:00:00'
    for obsformat in obsformats:
        try:
            obsdate = datetime.datetime.strptime(obsdate, obsformat)      # datetime object
        except:
            continue
        break               # Found a working time
    obsdate -= datetime.timedelta(0, 3600*params['raw_data_timezone_cor'])         # days, seconds, then other fields -> add the time zone difference, raw_data_timezone_cor: + for east of UTC, - for west of UTC -> subtract it
    
    fraction = 0.5                                          # mid exposure is half of the exposure time
    for exp_fraction in exp_fraction_keys:
        if exp_fraction in im_head.keys():                  # Weight of the exposure in header
            fraction = im_head[exp_fraction]                   # replace the fraction of half the exposure time
            break
    # Get the exposure time    
    exposure_time = -1
    for header_key in [ 'HiFLEx '+params['raw_data_exptim_keyword'], params['raw_data_exptim_keyword'] ]:
        if header_key in im_head.keys():
            exposure_time = im_head[header_key]
    if exposure_time == -1:
        logger('Warn: Cannot find the raw_data_exptim_keyword = {0} in the header. Assuming 0 seconds.'.format(params['raw_data_exptim_keyword'] ))
        exposure_time = 0
        
    obsdate_midexp = obsdate + datetime.timedelta(0, fraction*exposure_time)      # days, seconds, then other fields.
    jd_midexp = get_julian_datetime(obsdate_midexp)
    jd_begin = get_julian_datetime(obsdate)
    im_head['HIERARCH HiFLEx DATE-OBS']   = (obsdate.strftime("%Y-%m-%dT%H:%M:%S.%f"), 'UTC, Begin of expose')
    im_head['HIERARCH HiFLEx DATE-MID']   = (obsdate_midexp.strftime("%Y-%m-%dT%H:%M:%S.%f"), 'UTC, Middle of expose')
    im_head['HIERARCH HiFLEx EXPOSURE']   = (exposure_time, 'Exposure time in s')
    im_head['HIERARCH HiFLEx EXP_FRAC']   = (fraction, 'Normalised mean exposure time')
    im_head['HIERARCH HiFLEx JD']         = (round(jd_midexp,6), 'mid-exposure Julian date')             # MJD = JD - 2400000.5 from http://www.csgnetwork.com/julianmodifdateconv.html
    im_head['HIERARCH HiFLEx MJD']        = (round(jd_midexp - Constants.MJD0,6), 'mid-exposure modified JD')        # round 5 -> precision is 1 second, timing is not more precise
    im_head['HIERARCH HiFLEx JD_START']   = (round(jd_begin,6), 'JD at start of exposure')
    im_head['HIERARCH HiFLEx MJD_START']  = (round(jd_begin - Constants.MJD0,6), 'modified JD at begin of exposure')

    epoch = datetime.datetime.utcfromtimestamp(0)
    obsdate_mid_float = (obsdate_midexp - epoch).total_seconds()                                   # (obsdate - epoch) is a timedelta
    return im_head, obsdate_midexp, obsdate_mid_float, jd_midexp
    #return obsdate_midexp, obsdate_mid_float, exposure_time, obsdate, fraction, jd_midexp
    
def extraction_wavelengthcal(params, im, im_name, im_head, sci_traces, cal_traces, wave_sol_dict, reference_lines_dict, im_trace, objname):
    """
    Extracts the wavelength calibration
    
    """
    shift = 0
    [sci_tr_poly, xlows, xhighs, widths] = sci_traces
    [cal_tr_poly, axlows, axhighs, awidths] = cal_traces
    im_head, obsdate_midexp, obsdate_mid_float, jd_midexp = get_obsdate(params, im_head)               # in UTC, mid of the exposure
    if params['arcshift_side'] == 0:                       # single fibre spectrograph
        shift, im_head = find_shift_images(params, im, im_trace, sci_traces, 1, cal_tr_poly, im_head=im_head)     # w_mult=1 so that the same area is covered as for the find traces
        #shift = 0
        #logger('Warn: Line 3469: Searching for shift is switched off in the code')
    if im_name.endswith('_wavecal'):
        aspectra, agood_px_mask, extr_width = extract_orders(params, im, cal_tr_poly, axlows, axhighs, awidths, params['arcextraction_width_multiplier'], offset=shift, var='standard', plot_tqdm=False, header=im_head, plot_traces=True)
        fib = 'cal'
    elif im_name.endswith('_wavesci'):   
        aspectra, agood_px_mask, extr_width = extract_orders(params, im, sci_tr_poly, xlows, xhighs, widths, params['arcextraction_width_multiplier'], offset=shift, var='standard', plot_tqdm=False, header=im_head, plot_traces=True)
        fib = 'sci'
    else:
        logger('Error: The filename does not end as expected: {0} . It should end with _wavecal or _wavesci. This is probably a programming error.'.format(im_name))
    wavelength_solution_shift, shift, im_head = shift_wavelength_solution(params, aspectra, wave_sol_dict, reference_lines_dict, 
                                                              xlows, xhighs, obsdate_mid_float, jd_midexp, sci_tr_poly, cal_tr_poly, objname, fib=fib, im_head=im_head)   # This is only for a shift of the pixel, but not for the shift of RV
    wavelengths, dummy = create_wavelengths_from_solution(params, wavelength_solution_shift, aspectra)
    im_head['Comment'] = 'File contains a 3d array with the following data in the form [data type, order, pixel]:'
    im_head['Comment'] = ' 0: wavelength for each order and pixel in barycentric coordinates'
    im_head['Comment'] = ' 1: spectrum of the emission line lamp'
    im_head['Comment'] = ' 2: Mask with good areas of the spectrum: 0.1=saturated_px, 0.2=badpx'
    ceres_spec = [wavelengths, aspectra, agood_px_mask]
    save_multispec(ceres_spec, params['path_extraction']+im_name, im_head, bitpix=params['extracted_bitpix'])
    wave_sol_dict_new = copy.deepcopy(wave_sol_dict)
    wave_sol_dict_new['wavesol'] = wavelength_solution_shift
    save_wavelength_solution_to_fits( wave_sol_dict_new, params['path_extraction']+im_name.replace('_wavecal', '_wavesol_{0}_cal').replace('_wavesci', '_wavesol_{0}_sci').format(im_head['HiFLEx JD']) )
    

def get_possible_object_names(filename, header, header_keywords, replacements=['_arc','arc', '_thar','thar', '_une','une', 'extract_combine']):
    """
    Analyses the filename in order to find the possible Name of the Object (for exampled stored in parameter object_file
    The filename is subsequently stripped from the "_" or "-" separated entings
    :param filename: string, filename with removed path, file ending, and \n
    :param replacements: list of strings with entries to be removed from the filename
    return obnames: list of strings with possible object names
    """
    obnames = []
    for header_keyword in header_keywords:                  # Check the header first, if the object is stored there
        if header_keyword in header.keys():
            new = header[header_keyword].replace(' ','')
            if new not in obnames and len(new) >= 2:
                obnames.append(new)
    first_entry = filename.replace('-','_').split('_')      # most likely object is without any _ and -
    if first_entry[0] not in obnames:
        obnames.append(first_entry[0])
    obname = filename + '-'                   # Have at least one run
    while obname.find('_') != -1 or obname.find('-') != -1:
        for splitter in ['-', '_']:
            obname = obname.rsplit(splitter, 1)     # remove stuff from the right
            if len(obname) == 1:            # Splitter isn't part of filename anymore
                obname = obname[0]          # remove the list from it
                continue
            for i in range(5):              # Try to extract different information from the filename
                if i == 0:
                    obnametemp = obname[1]                      # check what will be stripped away in case the filename is like "unuseful_objectname-rest
                if (( i == 1 or i == 3 ) and obname[0].find('_') == -1) or (( i == 2 or i == 3 ) and obname[0].find('-') == -1):
                    continue                                    # if nothing has to be replaced, then running the step would be a waste of time
                if i == 1:
                    obnametemp = obname[0].replace('_','')      # in case the object name has a "_" in the name, but the object in the object_file doesn't
                if i == 2:
                    obnametemp = obname[0].replace('-','')      # in case the object name has a "-" in the name, but the object in the object_file doesn't
                if i == 3:
                    obnametemp = obname[0].replace('-','').replace('_','')  # in case the object name has both a "-" and "_" in the name, but the object in the object_file doesn't
                if i == 4:
                    obnametemp = obname[0]                      # use the front bit of the object name without any modification
                for rplc in replacements:               # Ignore everything that has to do with calibration
                    posi = obnametemp.lower().find(rplc)
                    if posi + len(rplc) == len(obnametemp) and posi != -1:           # the searchtext is at the end of the filename
                        obnametemp = obnametemp[:posi]              # get rid of the Arc in the filename
                    if posi == 0:
                        obnametemp = obnametemp[len(rplc):]
                    if i == 0 and len(obnametemp) < 5:
                        continue
                    if obnametemp not in obnames and obnametemp.lower() not in replacements:
                        obnames.append(obnametemp)
                if i == 4:                            # Use the stripped filename as new filename   
                    obname = obnametemp

    return obnames

def find_file_in_allfolders(filename, folders=[]):
    """
    Returns a list of places where to find the filename, both alone and checked in all folders
    :param filename: string: filename, can contain a path
    :param folders: list of strings: paths in which the filename and the filename without its path should be searched
    :return pathfiles: list of strings: All places in which the file was found
    """
    filename_no_path = filename.rsplit(os.sep,1)[-1]
    pathfiles, pathfiles_full = [], []
    for folder in [''] + folders:      # Check also the result and raw data paths for object names
        for fname in [filename, filename_no_path]:
            totest = folder+fname
            if totest not in pathfiles_full:
                pathfiles_full.append( totest )
                if os.path.isfile(totest):
                    pathfiles.append( totest )
    
    return pathfiles, pathfiles_full
    
def create_blaze_norm(params, im_trace1, sci_traces, cal_traces, wave_sol_dict, reference_lines_dict):
    """
    Creates the file for the blaze correction
    """
    fit_poly_orders = 15
    minflux = 0.002          # The blaze below 
    [sci_tr_poly, xlows, xhighs, widths] = sci_traces
    [cal_tr_poly, axlows, axhighs, awidths] = cal_traces
    im_blazecor, im_head = create_image_general(params, 'blazecor')
    logger('Step: Create the normalised blaze for the night')
    im_head, obsdate_midexp, obsdate_mid_float, jd_midexp = get_obsdate(params, im_head)
    shift, im_head = find_shift_images(params, im_blazecor, im_trace1, sci_traces, 0, cal_tr_poly, extract=True, im_head=im_head)
    #shift = 0  # for test
    flat_spec, good_px_mask, extr_width = extract_orders(params, im_blazecor, sci_tr_poly, xlows, xhighs, widths, params['extraction_width_multiplier'], var='linfit', offset=shift)
    #flat_spec, good_px_mask, extr_width = extract_orders(params, im_blazecor, sci_tr_poly, xlows, xhighs, widths, params['extraction_width_multiplier'], offset=shift) # for test
    med_flux = np.nanmedian(flat_spec)
    flat_spec_norm = flat_spec/med_flux
    flat_spec_norm_cor = correct_blaze(flat_spec_norm, minflux=0.01)         # Ignore all areas where the flux is 0.1% of median flux
    if params['blazercor_function'].find('poly') == 0:
        try:
            fit_poly_orders = int(params['blazercor_function'].replace('poly',''))
        except:
            'keep the standard'
    blaze_fit = fit_blazefunction(params['logging_blaze_spec'], flat_spec_norm, fit_poly_orders, minflux)
    if np.nansum(np.abs(sci_tr_poly - cal_tr_poly)) == 0.0:                   # science and calibration traces are at the same position
        blazecor_spec, agood_px_mask = flat_spec*0, copy.copy(good_px_mask)
    else:
        blazecor_spec, agood_px_mask, extr_width = extract_orders(params, im_blazecor, cal_tr_poly, axlows, axhighs, awidths, params['arcextraction_width_multiplier'], offset=shift)
    wavelength_solution_shift, shift, im_head = shift_wavelength_solution(params, blazecor_spec, wave_sol_dict, reference_lines_dict,
                                            xlows, xhighs, obsdate_mid_float, jd_midexp, sci_tr_poly, cal_tr_poly, params['master_blaze_spec_norm_filename'], im_head=im_head)
    wavelengths, dummy = create_wavelengths_from_solution(params, wavelength_solution_shift, blazecor_spec)
    im_head['Comment'] = 'File contains a 3d array with the following data in the form [data type, order, pixel]:'
    im_head['Comment'] = ' 0: wavelength for each order and pixel in barycentric coordinates'
    im_head['Comment'] = ' 1: extracted spectrum normalised by the global median flux ({0} ADU)'.format(int(med_flux))
    im_head['Comment'] = ' 2: 1, with signal below 1% of the flux removed'
    im_head['Comment'] = ' 3: fit of 1 with a polynomial of order {0}'.format(fit_poly_orders)
    im_head['Comment'] = ' 4: spectrum of the emission line lamp'
    save_multispec([wavelengths, flat_spec_norm, flat_spec_norm_cor, blaze_fit, blazecor_spec], params['master_blaze_spec_norm_filename'], im_head)

def fit_blazefunction(fname, spec, fit_poly_orders, minflux=0.01):
    """
    Fits a polynomial to the blaze function
    :param spec: 2d array with the spectra
    :param fit_poly_orders: integer, number of orders to be used for the polynomial fit
    """
    xarr = np.arange(spec.shape[1], dtype=int)
    spec_fit = copy.deepcopy(spec)
    titel_f = 'Extracted spectra of order {0} from the blaze file and fit with a polynomial of order {1}'
    with PdfPages(fname) as pdf:
        for order in range(spec.shape[0]):
            goodvalues = ~np.isnan(spec[order,:])
            poly = np.polyfit(xarr[goodvalues], spec[order, goodvalues], fit_poly_orders)
            fit = np.polyval(poly,xarr[goodvalues])
            fit[ fit < minflux*np.max(fit) ] = np.nan       # remove all that is 0.1% of maximum flux
            spec_fit[order, goodvalues] = fit
            fig, frame = plt.subplots(1, 1)
            title = titel_f.format(order, fit_poly_orders)
            plot_img_spec.plot_points([xarr,xarr],[spec[order,:],spec_fit[order,:]],['data','fit'],'', show=False, title=title, x_title='x [px]', 
                                      y_title='Flux (normalsied to global median) [ADU]', marker=['o',''], linestyle=["",'-'], frame=frame, return_frame=True)
    
            fig.set_size_inches(16.2, 10)
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()
    
        # We can also set the file's metadata via the PdfPages object:
        d = pdf.infodict()
        d['Title'] = 'Blaze function '
        d['Author'] = 'HiFLEx pipeline'
        d['Subject'] = ''
        d['Keywords'] = 'Blaze function HiFLEx pipeline'
        d['CreationDate'] = datetime.datetime.today()
        d['ModDate'] = datetime.datetime.today()
    
    return spec_fit

def read_create_spec(params, fname, im, im_head, trace_def, wmult, offset):
    """
    For the wavelength solution: reads the extracted spectrum or creates it
    """
    if os.path.isfile(fname):
        logger('Info: extracted emission spectrum already exist: {0}'.format(fname))
        spec = np.array(fits.getdata(fname))
        good_px_mask = spec[1,:,:]
        emission_spec = spec[0,:,:]
    else:           # Create the extracted spectrum
        [tr_poly, xlows, xhighs, widths] = trace_def
        emission_spec, good_px_mask, extr_width = extract_orders(params, im, tr_poly, xlows, xhighs, widths, wmult, offset=offset, header=im_head, plot_traces=True)
        spec = np.array([emission_spec, good_px_mask])
        im_head['Comment'] = 'File contains a 3d array with the following data in the form [data type, order, pixel]:'
        im_head['Comment'] = ' 0: spectrum of the emission line lamp'
        im_head['Comment'] = ' 1: Mask with good areas of the spectrum: 0.1=saturated_px, 0.2=badpx'
        save_multispec(spec, fname, im_head)
    
    return emission_spec, good_px_mask

def extraction_steps(params, im, im_name, im_head, sci_traces, cal_traces, wave_sol_dict_cal, reference_lines_dict, flat_spec_norm, im_trace):
    """
    Extracts the spectra and stores it in a fits file
    
    """
    [sci_tr_poly, xlows, xhighs, widths] = sci_traces
    [cal_tr_poly, axlows, axhighs, awidths] = cal_traces
    if 'HiFLEx BCKNOISE' not in im_head.keys():        # if not already done because background is in parameters
        bck_noise_std, bck_noise_var = prepare_measure_background_noise(params, im, sci_tr_poly, xlows, xhighs, widths, cal_tr_poly, axlows, axhighs, awidths)
        if np.isnan(bck_noise_var):
            bck_noise_var = -1
        im_head['HIERARCH HiFLEx BCKNOISE'] = (round_sig(bck_noise_std,4), 'Background noise in ADU')
        im_head['HIERARCH HiFLEx BNOISVAR'] = (round_sig(bck_noise_var,4), 'Variation of noise through image')             # Background noise variation can be very high, because some light of the traces remains
    if im_head['HiFLEx BCKNOISE'] <= 0 or np.isnan(im_head['HiFLEx BCKNOISE']):
        logger('Warn: Measured an unphysical background noise in the data: {0}. Set the noise to 1'.format(im_head['HiFLEx BCKNOISE']))
        im_head['HIERARCH HiFLEx BNOISVAR'] = (1., '1, because of unphysical measurement')
    im_head, obsdate_midexp, obsdate_mid_float, jd_midexp = get_obsdate(params, im_head)               # in UTC, mid of the exposure
    
    im_name = im_name.replace('.fits','').replace('.fit','')                # to be sure the file ending was removed
    # not anymore: Object name needs to be split by '_', while numbering or exposure time needs to be split with '-'
    obname = im_name.replace('\n','').split(os.sep)    # get rid of the path
    im_head['HIERARCH HiFLEx NAME'] = (obname[-1], 'original filename')       # To know later what was the original filename
    #not necessary anymore: obname = obname.split('-')  # remove the numbering and exposure time from the filename
    #not necessary anymore: obname = obname[0]              # contains only the name, e.g. ArturArc, SunArc
    obnames = get_possible_object_names(obname[-1], im_head, params['raw_data_object_name_keys'])
    
    # Change the path to the object_file to the result_path, if necessary
    object_file = params['object_file']         # needed later
    object_files, object_files_full = find_file_in_allfolders(object_file, [params['result_path']] + params['raw_data_paths'])      # Check also the result and raw data paths for object names
    ra2, dec2, pmra, pmdec, epoch = 0., 0., 0., 0., 2000.
    for object_file in object_files:        # Will be run in any case, even if empty
        ra2, dec2, epoch, pmra, pmdec, obnames, dummy = getcoords_from_file(obnames, 0, filen=object_file, warn_notfound=False)        # mjd=0 because because not using CERES to calculated BCV
        if ra2 !=0 or dec2 != 0:                                           # Found the objects -> obnames becomes the single entry which matched entry in params['object_file']
            break
    if ra2 ==0 and dec2 == 0:
        logger('Warn: Reference coordinates files do not exist or object can not be found within them. Checked: {0}'.format(object_files_full))
    im_head['HIERARCH HiFLEx OBJNAME'] = (obnames[0], 'Used object name')
    # Get the baycentric velocity
    params, bcvel_baryc, mephem, im_head = get_barycent_cor(params, im_head, obnames[0], ra2, dec2, epoch, pmra, pmdec, obsdate_midexp)
    
    shift, im_head = find_shift_images(params, im, im_trace, sci_traces, 0, cal_tr_poly, im_head=im_head)     # w_mult=1 so that the same area is covered as for the find traces, w_mult=0 so that only the central pixel is extracted
    spectra, good_px_mask, extr_width = extract_orders(params, im, sci_tr_poly, xlows, xhighs, widths, params['extraction_width_multiplier'], offset=shift, header=im_head, var=params['extraction_precision'], plot_traces=True)
    orders = range(spectra.shape[0])
    if np.nansum(np.abs(sci_tr_poly - cal_tr_poly)) == 0.0:                                                 # science and calibration traces are at the same position
        aspectra, agood_px_mask = spectra*0, copy.copy(good_px_mask)
    else:
        aspectra, agood_px_mask, aextr_width = extract_orders(params, im, cal_tr_poly, axlows, axhighs, awidths, params['arcextraction_width_multiplier'], offset=shift, var='standard', plot_tqdm=False, header=im_head)
    wavelength_solution_shift, shift, im_head = shift_wavelength_solution(params, aspectra, wave_sol_dict_cal, reference_lines_dict, 
                                                              xlows, xhighs, obsdate_mid_float, jd_midexp, sci_tr_poly, cal_tr_poly, im_name, im_head=im_head)   # This is only for a shift of the pixel, but not for the shift of RV
    wavelengths, wavelength_solution_shift = create_wavelengths_from_solution(params, wavelength_solution_shift, spectra, wave_sols_sci=calimages['wave_sols_sci'], jd_mid=jd_midexp, shift=shift)
    espectra = combine_photonnoise_readnoise(spectra, im_head['HiFLEx BCKNOISE'] * np.sqrt(np.abs(extr_width)) )     # widths[:,2] is gaussian width
    if params['blazercor_function'].find('poly') == 0:  index = 3
    elif params['blazercor_function'] == 'flux':        index = 2
    else:       index = 2           # 1: extracted flat, 2: low flux removed, 3: fitted with a high order polynomial
    im_head['HIERARCH HiFLEx BLAZECOR'] = (params['blazercor_function'], 'Parameter blazercor_function')
    im_head['HIERARCH HiFLEx BLAZECOR INDEX'] = (index, 'Dataset for blaze correction')
    fspectra = spectra/(flat_spec_norm[index]+0.0)        # 1: extracted flat, 2: low flux removed
    # Doing a wavelength shift for the flat_spec_norm is probably not necessay, as it's only few pixel
    measure_noise_orders = 16
    measure_noise_semiwindow = 10                   # in pixel
    efspectra = measure_noise(fspectra, p_order=measure_noise_orders, semi_window=measure_noise_semiwindow)             # Noise will be high at areas wih absorption lines, takes about a minute for 75o6000px
    cspectra, noise_cont = normalise_continuum(fspectra, wavelengths, nc=6, semi_window=measure_noise_semiwindow, nc_noise=measure_noise_orders)  # takes about 3 (22) minutes for 75o6000px for Tungsten (for noise)

    # normalise_continuum measures the noise different than measure_noise
    if len(params['sigmaclip_spectrum']) == 3:        # sigmaclipping
        removed = 0
        for order in range(cspectra.shape[0]):
            good_values, dummy = sigmaclip(cspectra[order,:], x=wavelengths[order,:], nc=int(params['sigmaclip_spectrum'][0]),
                                           sub=[], ll=params['sigmaclip_spectrum'][1], lu=params['sigmaclip_spectrum'][2],
                                           repeats=3, exclude_absorption=True)       # if ll is too low (10) large absorption lines might be cliped (Halpha)
            cspectra[order,~good_values] = np.nan
            removed += np.sum(~good_values)
        logger('Info: Sigmaclipping removed {0} data points for file {1} (normalised spectrum).'.format(removed, im_name))
        im_head['HIERARCH HiFLEx SIGMA_CLEARED']  = (removed, 'points removed from normalised spectrum')
        
    # Correct wavelength by barycentric velocity
    wavelengths_bary = wavelengths * (1 + bcvel_baryc/(Constants.c/1E3) )

    do_RV = True
    for no_RV_name in params['no_RV_names']:
        if im_name.lower().find(no_RV_name) in [0,1,2,3,4,5]:
            do_RV = False
            break

    im_head_bluefirst = copy.copy(im_head)
    im_head = add_specinfo_head(spectra, spectra, noise_cont, extr_width, im_head)
    im_head_bluefirst = add_specinfo_head(spectra[::-1,:], spectra[::-1,:], noise_cont[::-1,:], extr_width[::-1,:], im_head_bluefirst)
    wspec = good_px_mask* flat_spec_norm[1]
    doo = dict(spec=True, blaze=True, weight=True, norm=True, blue=True, harps=do_RV)

    save_single_files(params, obnames, im_name, im_head, im_head_bluefirst, spectra, fspectra, cspectra, wspec, wavelengths, wavelength_solution_shift, doo)

    if params['wavelength_scale_resolution'] > 0.0 and np.max(wave_sol_dict_cal['wavesol'][:,-1]) > 100:
        # Create a linearised solution for the input spectrum and the continuum corrected spectrum
        #logger('Step: Linearising the spectrum (commented out)')
        #wavelenghts_lin, spectrum_lin = linearise_wavelength_spec(params, wavelength_solution_shift, spectra, method='sum', weight=espectra)
        #save_multispec([wavelenghts_lin,spectrum_lin], params['path_extraction_single']+im_name+'_lin', im_head)
        if True:
            wavearray = wavelengths_bary
        else:
            wavearray = wavelengths
        wavelenghts_lin, spectrum_lin = linearise_wavelength_spec( params, wavearray, cspectra, method='weight', weight=((flat_spec_norm[index]+0.0)/espectra)**2 ) # highest throughput and lowest error should have the highest weight
        save_multispec([wavelenghts_lin,spectrum_lin], params['path_extraction_single']+im_name+'_lin_cont', im_head)
        with open(params['path_extraction_single']+im_name+'_lin_cont.csv', 'w') as file:
            for ii in range(wavelenghts_lin.shape[0]):
                if ~np.isnan(wavelenghts_lin[ii]) and ~np.isnan(spectrum_lin[ii]):
                    file.write('{1},{2}{0}'.format(os.linesep, round(wavelenghts_lin[ii],4), round(spectrum_lin[ii],5) ))
        
    # For easier plotting
    add_text_to_file(params['path_extraction']+im_name+'.fits', 'plot_files.lst')
    # Checked that telluric lines (e.g 7600 Angstom) are at the same position when using wavelength solution without barycentric correction and 
    #   that stellar lines are overlaying when using wavelength solution with barycentric correction (MRES 20190801)
    ceres_spec = np.array([wavelengths_bary, spectra, espectra, fspectra, efspectra, cspectra, noise_cont, good_px_mask, aspectra, wavelengths])        
    im_head['Comment'] = 'File contains a 3d array with the following data in the form [data type, order, pixel]:'
    im_head['Comment'] = ' 0: wavelength for each order and pixel in barycentric coordinates'
    im_head['Comment'] = ' 1: extracted spectrum'
    im_head['Comment'] = ' 2: measure of error (photon noise, read noise)'
    im_head['Comment'] = ' 3: flat corrected spectrum'
    im_head['Comment'] = ' 4: error of the flat corrected spectrum (residuals to a {0} order polynomial)'.format(measure_noise_orders)
    im_head['Comment'] = ' 5: continuum normalised spectrum'
    im_head['Comment'] = ' 6: error in continuum (fit to residuals of {0} order polynomial)'.format(measure_noise_orders)
    im_head['Comment'] = ' 7: Mask with good areas of the spectrum: 0.1=saturated_px, 0.2=badpx'
    im_head['Comment'] = ' 8: spectrum of the emission line lamp'
    im_head['Comment'] = ' 9: wavelength for each order and pixel without barycentric correction'
    save_multispec(ceres_spec, params['path_extraction']+im_name, im_head, bitpix=params['extracted_bitpix'])
    logger('Info: Finished extraction of {0}'.format(im_name))

    return obnames[0]

def save_single_files(params, obnames, im_name, im_head, im_head_bluefirst, spectra, fspectra, cspectra, wspec, wavelengths, wavelength_solution_shift, doo):    
    """
    Save the spectrum as single files
    The wavelengths won't be barycentric corrected
    """
    im_head_wave, im_head_weight = copy.copy(im_head), copy.copy(im_head)
    im_head_iraf_format = copy.copy(im_head)
    im_head_iraf_format_bluefirst = copy.copy(im_head_bluefirst)
    orders = list(range(spectra.shape[0]))
    
    ## Save in a easier way
    # Single files
    im_head_wave['Comment'] = 'Contains the wavelength per order and exctracted pixel for file {0}'.format(im_name+'_extr/_blaze/_norm')
    im_head_weight['Comment'] = 'Contains the weights per order and exctracted pixel for file {0}'.format(im_name+'_extr/_blaze/_norm')
    im_head_iraf_format = wavelength_solution_iraf(params, im_head_iraf_format, wavelengths, wavelength_solution_shift, norder=params['polynom_order_traces'][-1]+2)
    im_head_iraf_format_bluefirst = wavelength_solution_iraf(params, im_head_iraf_format_bluefirst, wavelengths[::-1,:], wavelength_solution_shift[::-1,:], norder=params['polynom_order_traces'][-1]+2)
    if doo.get('spec',False):
        save_multispec(spectra,     params['path_extraction_single']+im_name+'_extr', im_head_iraf_format, bitpix=params['extracted_bitpix'])
    if doo.get('blaze',False):
        save_multispec(fspectra,    params['path_extraction_single']+im_name+'_blaze', im_head_iraf_format, bitpix=params['extracted_bitpix'])
    if doo.get('norm',False):
        save_multispec(cspectra,    params['path_extraction_single']+im_name+'_norm', im_head_iraf_format, bitpix=params['extracted_bitpix'])
    if doo.get('spec',False) or doo.get('blaze',False):
        save_multispec(wavelengths, params['path_extraction_single']+im_name+'_wave', im_head_wave, bitpix=params['extracted_bitpix'])
    if (doo.get('spec',False) or doo.get('blaze',False)) and doo.get('weight',False):
        save_multispec(wspec,       params['path_extraction_single']+im_name+'_weight', im_head_weight, bitpix=params['extracted_bitpix'])  # Weight shouldn't be the espectra, the flat provides a smoother function
    if doo.get('spec',False) and doo.get('blue',False):
        save_multispec(spectra[::-1,:],  params['path_extraction_single']+im_name+'_extr_bluefirst',  im_head_iraf_format_bluefirst, bitpix=params['extracted_bitpix'])
    if doo.get('blaze',False) and doo.get('blue',False):
        save_multispec(fspectra[::-1,:], params['path_extraction_single']+im_name+'_blaze_bluefirst', im_head_iraf_format_bluefirst, bitpix=params['extracted_bitpix'])
    if doo.get('norm',False) and doo.get('blue',False):
        save_multispec(cspectra[::-1,:], params['path_extraction_single']+im_name+'_norm_bluefirst', im_head_iraf_format_bluefirst, bitpix=params['extracted_bitpix'])
    
    if doo.get('harps',False):
        # Harps format
        im_head_harps_format = copy.copy(im_head_bluefirst)
        im_head_harps_format = wavelength_solution_harps(params, im_head_harps_format, wavelengths[::-1,:])        # [::-1,:] -> Blue orders first      # 20190509: wavelengths instead of wavelengths_bary
        serval_keys = []
        serval_keys.append(['INSTRUME', 'HIFLEX',                                                 'added for Serval'])
        serval_keys.append(['EXPTIME',  im_head_harps_format.get('HiFLEx EXPOSURE',0),                    'Exposure time, for Serval'])
        serval_keys.append(['DATE-OBS', im_head_harps_format.get('HiFLEx DATE-OBS',0),                    'UT start, for Serval'])
        serval_keys.append(['MJD-OBS',  im_head_harps_format.get('HiFLEx MJD_START',0),                   'MJD start ({0})'.format(im_head_harps_format.get('HiFLEx DATE-OBS',0)) ])
        if 'HiFLEx BJDTDB' in im_head.keys():
            serval_keys.append(['HIERARCH ESO DRS BJD',     im_head_harps_format['HiFLEx BJDTDB'],        'Barycentric Julian Day'])      # DRS produces them without leap seconds, e.g. 68.2s earlier at 2015
        if 'HiFLEx RA' in im_head_harps_format.keys():
            serval_keys.append(['RA',       im_head_harps_format['HiFLEx RA'],                            'RA start, for Serval'])
        if 'HiFLEx DEC' in im_head_harps_format.keys():
            serval_keys.append(['DEC',      im_head_harps_format['HiFLEx DEC'],                           'DEC start, for Serval'])
        if 'HiFLEx BCV' in im_head_harps_format.keys():
            serval_keys.append(['HIERARCH ESO DRS BERV',     im_head_harps_format['HiFLEx BCV'],          'Barycentric Earth Radial Velocity'])
            serval_keys.append(['HIERARCH ESO DRS BERVMX',   im_head_harps_format['HiFLEx BCV MAX'],      'Maximum BERV'])
            serval_keys.append(['HIERARCH ESO DRS BERVMN',   im_head_harps_format['HiFLEx BCV MIN'],      'Minimum BERV'])
        serval_keys.append(['HIERARCH ESO DPR TECH',         'ECHELLE ',        'Observation technique'])
        serval_keys.append(['HIERARCH ESO INS MODE',         'HiFLEx',          'Instrument mode used.'])
        serval_keys.append(['HIERARCH ESO DRS CAL LOC NBO',  spectra.shape[0],  'nb orders localised'])
        serval_keys.append(['HIERARCH ESO OBS TARG NAME',    obnames[0],        'OB target name'])
        serval_keys.append(['OBJECT',                        obnames[0],        'OB target name'])
        serval_keys.append(['HIERARCH ESO INS DET1 TMMEAN',  im_head_harps_format.get('HiFLEx EXP_FRAC',0.5),    'Normalised mean exposure time'])
        serval_keys.append(['HIERARCH ESO INS DET2 TMMEAN',  im_head_harps_format.get('HiFLEx EXP_FRAC',0.5),    'Normalised mean exposure time'])
        #serval_keys.append(['',         '',        ''])
        for order in orders:
            serval_keys.append([ 'HIERARCH ESO DRS SPE EXT SN{0}'.format(order), im_head.get('HiFLEx SN_order{0}'.format('%2.2i'%order),10), 'S_N order center{0}'.format(order) ])
        for [newkey, value, comment] in serval_keys:
            if newkey not in im_head_harps_format.keys():
                im_head_harps_format[newkey] = (value, comment)
        if 'COMMENT' in im_head_harps_format.keys():
            del im_head_harps_format['COMMENT']                 # Serval can't read comments
        for entry in im_head_harps_format.keys():
            #print "key, value, comment",(entry, im_head_harps_format[entry], im_head_harps_format.comments[entry])
            if im_head_harps_format.comments[entry] == '':
                im_head_harps_format.comments[entry] = '/'      # Serval can't read header keywords without comment, for NAXISj this needs to be done in save_multispec
        fname = params['path_harpsformat']+obnames[0].lower()+os.sep+'HARPS.{0}_e2ds_A.fits'.format(im_head_harps_format.get('HiFLEx DATE-OBS','dummytime_  ')[:-3])
        if len(fname.rsplit(os.sep,1)) == 1:     # no path is in the filename
            logger('Warn: no folder to save {0} was given, using the current folder ({1}).'.format( fname, os.getcwd() ))
        else:
            make_directory(fname.rsplit(os.sep,1)[0], errormsg='Error: Folder to save {0} does not exists and cannot be created.'.format(fname))
        #for entry in im_head_harps_format.keys():
        #    if entry.find('AXIS') != -1:
        #        print "key, value, comment",(entry, im_head_harps_format[entry], im_head_harps_format.comments[entry]) 
        save_multispec(spectra[::-1,:], fname, im_head_harps_format, bitpix=params['extracted_bitpix'])                 # [::-1,:] -> Blue orders first
    
    
def save_spec_csv(spec, wavelengths, good_px_mask, fname):
    """
    Save the spectra in a csv file to be compatible with Terra. The files need to have the form:
    ONUM,WAVELENGTH (in Angstrongs!),FLUX (arbitrary units, normalization to order each mean recommended).
    Note: all orders must of of the SAME LENGTH
    [email from Guillem Anglada, 18/10/2018 12:25
    """
    if len(fname.rsplit(os.sep,1)) == 1:     # no path is in the filename
        logger('Warn: no folder to save {0} was given, using the current folder ({1}).'.format( fname, os.getcwd() ))
    else:
        make_directory(fname.rsplit(os.sep,1)[0], errormsg='Error: Folder to save {0} does not exists and cannot be created.'.format(fname))
        
    specs = spec.shape
    spec_cor = spec * good_px_mask
    spec_cor[np.isnan(spec_cor)] = 0            # replace NaNs with 0
    #spec_cor[spec_cor < 0] = 0                  # replace negative values (test 20190701 after terra RVs are NaN - > doesn't solve the NaNs => problems was a full order without data)
    wave = copy.copy(wavelengths)
    wave[np.isnan(wave)] = 0
    fname = fname.replace('.csv','') + '.csv'
    with open(fname, 'w') as file:
        for order in range(specs[0]):
            for px in range(specs[1]):
                file.write('{0},{1},{2}\n'.format(order, wave[order,px], spec_cor[order,px]) )
    #print order,px

def wavelength_solution_harps(params, head, wavelengths):
    """
    Save the wavelength solution in the same way as the ESO DRS 
    Don't use the barycentric corrected wavelengths
    :return head: Modified header
    """
    ws = wavelengths.shape
    deg = max(params['polynom_order_traces'])
    head['HIERARCH ESO DRS CAL TH DEG LL'] = (deg, 'degre polyn fit ll(x,order)')
    x = np.arange(ws[1])
    for order in range(ws[0]):
        poly = np.polyfit( x, wavelengths[order,:], deg)
        #print poly, wavelengths[order,0], wavelengths[order,-1], max(np.abs(np.polyval(poly, x) - wavelengths[order,:]))
        for i in range(deg+1):
           head['HIERARCH ESO DRS CAL TH COEFF LL{0}'.format( i + order * ( deg + 1 ) )] = (poly[-(i+1)], 'coeff ll(x,order)')
    
    return head

def wavelength_solution_iraf(params, head, wavelengths, wavelength_solution, norder=10):
    """
    Save the wavelength solution in the same way as IRAF does
    Sources: chapter 5 in http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?specwcs
    Tested with: https://github.com/kgullikson88/General/blob/master/readmultispec.py
    :param norder: Degree of the fitting polynomials
    :return head: Modified header
    """
    ws = wavelengths.shape
    deg = max(params['polynom_order_traces'])
    head['WAT0_001'] = 'system=multispec'
    head['WAT1_001'] = 'wtype=multispec label=Wavelength units=Angstroms'
    text = 'wtype=multispec'
    for order in range(ws[0]):
        x_fit = np.arange(ws[1])+1
        pmin = 1                    # maybe replace this by the first non-nan pixel
        pmax = ws[1]  # ! without +1  # maybe replace this by the last non-nan pixel
        pmiddle = (pmax + pmin) / 2.
        prange = pmax - pmin
        x_fit = (np.arange(ws[1],dtype=float) + 1. - pmiddle) / (prange / 2.)
        y_fit = wavelengths[order,:]        # needs to be changed, if only non-nan flux is used
        avg_dwave = np.nanmean(wavelengths[order,1:] - wavelengths[order,:-1])
        cheb_pol = np.polynomial.chebyshev.chebfit(x_fit, y_fit, norder)
        cheb_pol_txt = np.str(cheb_pol)[1:-1].replace('\n',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ')
        text += ' spec{0} = "{1} {2} 2 {3} {4} {5} 0. {0}. {0}. 1. 0. 1 {6} {7} {8} {9}"'.format(
                            order+1, order, int(wavelength_solution[order,0]), round_sig(wavelengths[order,0],6), round_sig(avg_dwave,6), ws[1], 
                            len(cheb_pol), float(pmin), float(pmax), cheb_pol_txt)
        # 1/1: order (arbitary)
        # 2/2: wavelength_solution[order,0] : real order
        # 3: 2 for non-linear solution
        # 4/3: wavelength at first pixel
        # 5/4: average dispersion
        # 6/5: number of valid pixel
        # 7: 0. as doppler correction
        # 8/0: extraction aperture limit (float, as in iraf created spectrum: Keck-Spectra/3218_706r_wc.fits)
        # 9/0: extraction aperture limit (float, as in iraf created spectrum)
        # 10: 1. as weight
        # 11: 0. as zero point offset
        # 12: 1 for chebychev dispersion function
        # 13/6: order of the polynomial: len(cheb_pol) = norder + 1
        # 14/7: 1. as lowest pixel
        # 15/8: highest good pixel
        # 16-end/9: coefficients from the chebychev fit
    
    split = 68
    for i in range(int(len(text)/split)):
        head['WAT2_%3.3i'%(i+1)] = text[i*split:(i+1)*split]
    head['WAT2_%3.3i'%(i+2)] = text[(i+1)*split:]               # remaining text
    if i+2 > 999:
        logger('Warn: IRAF wavelength solution is too long. Try to use a smaller degree of freedom to fit the wavelength solution (polynom_order_traces)')
    return head
    
def wavelength_air_to_vacuum(wave):         # from CERES, naming and structure changed
    """
    Converts the the Air wavelength into vacuum wavelength, using the material in 
    https://www.as.utexas.edu/~hebe/apogee/docs/air_vacuum.pdf
    :param wave: nd array of floats, containing the air wavelengths per pixel
    :return wave_new: nd array of floats, same format as wave, containing the vacuum wavelengths per pixel
    """
    wave_prev = wave.copy()
    while(True):                                # Iterative 
        wave_new = n_Edlen(wave_prev) * wave
        if (np.nanmax(np.absolute(wave_new - wave_prev), axis=None) < 1e-10):
            break
        wave_prev = wave_new
    return wave_new

def wavelength_vacuum_to_air(wave):         # from CERES, naming changed
    return (wave / n_Edlen(wave))

def n_Edlen(l):
    """
    Calculates the dispersion index, depending on the wavelength l
    :params l: nd array of floats, containing the air wavelengths
    return n: nd array of floats, same format as l, distpersion index for each wavelength
    """
    s = 1e4 / l
    s2 = s*s
    # The following line "(number-s2)" or "(s2-number)" isn't calculated when the pipeline creates the python 2 environment and this is called through "run_ceres_multicore" -> "wavelength_air_to_vacuum"; however it works fine when called through "run_ceres" -> "wavelength_air_to_vacuum"
    n = 1 + 8.34213e-5 + 2.406030e-2/(130-s2) + 1.5997e-4/(38.9-s2)    # applied from https://www.as.utexas.edu/~hebe/apogee/docs/air_vacuum.pdf same formular as in CERES
    return n

def plot_traces_over_image(im, fname, pfits, xlows, xhighs, widths=[], w_mult=1, offset=0, mask=[], frame=None, return_frame=False, color=['r','b','g'], showtext=True, imscale=None):
    """
    Plot the found traces over the CCD image in order to allow the detection of mistakes made by the automatic steps
    :param im: 2d array of floats with the image
    :param fname: String, file to which should be saved
    :param pfits: list of floats, length same as number of orders, polynomial values
    :param xlows: list, length same as number of orders, the lowest x pixel (wavelength direction) used in each order
    :param xhighs: list, length same as number of orders, the highest x pixel (wavelength direction) used in each order
    :param widths: 2d list, length same as number of orders, each entry contains left border, right border, and Gaussian width of the lines, as estimated in the master flat
    :param w_mult: w_mult * widths (Gaussian) defines the range on either side of the trace to be extracted
    :param mask: list of Boolean, length same as number of orders, only entries True are plotted
    :param frame: if the results should be plotted in a GUI, then a frame already exists
    :param return_frame: if the results should be plotted in a GUI, then the frame needs to returned
    :return frame: if return_frame=True the frame is returned
    """
    leftfit, rightfit = [], []
    if len(pfits.shape) == 3:                           # new 3d array
        leftfit = pfits[:,1,:]
        rightfit = pfits[:,2,:]
        pfits = pfits[:,0,:]
    if widths == []:
        widths = np.repeat([[0,0,0]], len(pfits), axis=0)
    if mask == []:
        mask = np.repeat([True], len(pfits), axis=0)
    ims = im.shape
    im = scale_image_plot(im,'log10')
    colorbar = False
    if frame == None:
        fig, frame = plt.subplots(1, 1)
        colorbar = True
    title=''
    if showtext:
        title = 'Plot the traces in the image (log10 of image).'
        if np.mean(widths, axis=None) == 0 and np.std(widths, axis=None, ddof=1) == 0:
            title += ' The marked width (dashed lines) are shown for an extraction width multiplier of {0}.'.format(w_mult)
    if imscale is not None: pctile = imscale
    else:                   pctile = 0.2        # 1 was too high for MRES data binning 
    plot_img_spec.plot_image(im, [], pctile=pctile, show=False, adjust=[0.05,0.95,0.95,0.05], title=title, return_frame=True, frame=frame, autotranspose=False, colorbar=colorbar)
    colors = color*len(pfits)
    for pp, pf in enumerate(pfits):
            if mask[pp] == False:
                continue
            xarr = np.arange(xlows[pp], xhighs[pp], 1)
            yarr = np.polyval(pf[1:], xarr-pf[0]) + offset
            if len(leftfit) == 0:                        # old 2d array, use gaussian width (single value) as width
                yarrl = yarr - widths[pp][2]*w_mult
                yarrr = yarr + widths[pp][2]*w_mult
            else:
                #yarrl = np.polyval(leftfit[pp, 1:], xarr - leftfit[pp, 0]) - widths[pp][0]*(w_mult-1)
                #yarrr = np.polyval(rightfit[pp, 1:], xarr- rightfit[pp, 0]) + widths[pp][1]*(w_mult-1)
                yarrl = np.polyval(leftfit[pp, 1:], xarr - leftfit[pp, 0]) + offset
                yarrr = np.polyval(rightfit[pp, 1:], xarr- rightfit[pp, 0]) + offset
                yarrl, yarrr = adjust_width_orders(yarr, yarrl, yarrr, [w_mult, w_mult])              # Adjust width
            # Take into account the boundaries of the image
            yarr[yarr < 0] = 0
            yarr[yarr > ims[1]] = ims[1]
            yarrl[yarrl < 0] = 0
            yarrl[yarrl > ims[1]] = ims[1]
            yarrr[yarrr < 0] = 0
            yarrr[yarrr > ims[1]] = ims[1]
            mid_pos = int(len(xarr)/2)
            ymid = yarr[mid_pos]
            xmid = xarr[mid_pos]
            #if pp == 24:
            #    print 'xarr[0],xarr[-1],yarr[0],yarr[-1],yarrl[0],yarrl[-1],yarrr[0],yarrr[-1],pf',xarr[0],xarr[-1],yarr[0],yarr[-1],yarrl[0],yarrl[-1],yarrr[0],yarrr[-1],pf

            frame.plot(yarr,  xarr, color=colors[pp], linewidth=1)
            frame.plot(yarrl, xarr, color=colors[pp], linewidth=1, linestyle='dashed')
            frame.plot(yarrr, xarr, color=colors[pp], linewidth=1, linestyle='dashed')
            if showtext:
                frame.text(ymid, xmid, 'Order{0}'.format(pp),
                       horizontalalignment='center', verticalalignment='center',
                       rotation=90, color=colors[pp], zorder=5)
    if return_frame:
        return frame
    fig.set_size_inches(42, 42)        
    plt.savefig(fname, bbox_inches='tight')
    plt.close()

def scale_data(y, x=[], mode='gauss'):
    """
    Scales the data to be between 0 and 1
    :param y: list or array of floats or ints with the data
    if mode=='gauss' or 'poly3', then should be absorption line
    """
    if type(y).__name__ != 'ndarray':
        y = np.array(y)
    if len(x) != y.shape[0]:
        x = np.arange(len(y))
    elif type(x).__name__ != 'ndarray':
        x = np.array(x)
    maxy = max(y)               # Continuum
    miny = min(y)               # could be line
    rangey = maxy - miny
    popt = [0,0,0,0]
    if mode == 'gauss':
        # Try to get the minimum of the line with a fit (if undersampled)
        center = int(len(y)/2)
        width = max(x)-min(x)
        p0 = [rangey, x[center], 0.1*width, min(-y)]        # -y to make a positive Gauss
        bounds=((0.2*rangey, x[int(len(y)/4)], 0.01*width, min(-y)-0.2*rangey), (5.*rangey, x[int(len(y)*3/4)], 0.5*width, min(-y)+0.2*rangey))
        #print x,y,p0,bounds
        try:
            popt,pcov = curve_fit(oneD_gauss,x,-y,p0=p0, bounds=bounds)            #a, x0, sigma, b: a*np.exp(-(x-x0)**2/(2*sigma**2))+b
            #maxy = -popt[3]        # use not in case the line is not in the form of a Gauss
            #print miny, maxy,rangey, popt
            diff = np.abs( ( oneD_gauss( x, popt ) - y ) / y )
            if np.percentile(diff, 75) < 0.1:
                maxy = -popt[3]
                miny = -( popt[0] + popt[3] )
                rangey = popt[0]
        except:
            print('curve fit failed')
    elif mode == 'min':
        True        # Already all confirmed
    elif mode == 'poly3':
        cen_poly, popt,l,r = fit_poly_center(-y, x)
        if not np.isnan(cen_poly):
            miny = -np.polyval(popt, cen_poly)
            rangey = maxy - miny
        #else:
        #    print rangey, maxy, miny, cen_poly, popt, x, y
        #    exit()
    else:
        logger('Warn: Mode {0} is not implemented for procedure "scale_data"'.format(mode))
    y_scaled = (y - miny)*1./rangey
    
    return y_scaled, popt

def linear_interpolation(x, y, x_search):
    """
    Linear equation to interpolate. Similar to y_x=f(x_search); f=scipy.interpolate.interp1d(x,y,kind='linear'), but at least for single numbers it needs only 20% of time
    :param x, y: list or array of floats
    :param x_search: list or array of floats (one entry less than x,y) or float
    :return y_x: array or float (depending on type of x_seach), the interpolated x_search value.
    """
    numer = False
    if type(y).__name__ != 'ndarray':
        y = np.array(y)
    if type(x).__name__ != 'ndarray':
        x = np.array(x)
    if type(x_search).__name__ != 'ndarray':
        if type(x_search).__name__ == 'list':
            x_search = np.array(x_search)
        else:                   # It's a number
            number = True
            x_search = np.array([x_search])
    y_x = (y[:-1] - y[1:]) / (x[:-1] - x[1:])*(x_search-x[1:]) + y[1:]
    if number:
        y_x = y_x[0]
    
    return y_x

def create_grid_data(y, x=[], ysteps=20):
    """
    Grids the data into ysteps. First the data is fitted by a cubic spline
        (https://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html),
        afterwards the y-values interpolated to the gridded values
    :param y: list or array of floats or ints with the data, should be sorted
    :param ysteps: integer, number of steps for y
    :return data_g[0,:], data_g[1,:]: 1d arrays of the curve with values only at y
    """
    if len(x) != len(y):
        x = np.arange(len(y))
    # Fit a cubic spline
    if len(y) <= 3: kind = 'linear'
    else:           kind = 'cubic'
    f = scipy.interpolate.interp1d(x, y, kind=kind)
    xf = np.linspace( np.min(x), np.max(x), len(y)*100, endpoint=True ) 
    data = np.vstack(( xf,f(xf) ))
    """ Do it manually instead of just using a cubic spline: using a stepwise polynomial of grade 2
    start = x[0]
    data = np.zeros(0)
    data.shape = (3,0)
    while True:     # Go through the data left to right, and fit polynomials to short areas
        mask = ( ( x >= start) & (x <= start+xrang ) )
        multi = 1
        # Make sure we have enough data to fit, at least half of xrang is covered
        while start+xrang*multi <= max(x) and x[mask][-1]-x[mask][0] <= xrang/2.:
            multi *= 1.1                        # extend the range
            mask = ( ( x >= start) & (x <= start+xrang*multi ) )
        # Fit the data
        polx = np.polyfit( x[mask], y[mask], 2 )        #min( 2, len(np.unique((x[mask]*1E6).astype(int)))-1 )
        # Create a small stepped aray and calculate the y values for it
        xf = np.arange( x[mask][0], x[mask][-1]+xrang/1E3, xrang/1E3 )      
        y_f = np.polyval( polx, xf )
        dist = (xf-np.median(x[mask]))**2
        w = 1 - dist/np.max(dist)
        data = np.hstack(( data, np.vstack((xf,y_f,w)) ))
        
        #bfit = scipy.interpolate.splrep( x[mask],y[mask],k=min(3,sum(mask)-1) )       # Get the B-sline of the flux
        #y_f = scipy.interpolate.splev(x[mask],bfit)            # Apply the B-spline -> this will create the same points as the original points
        
        #"#""
        fig, frame = plt.subplots(1,1)
        frame.plot(x,y)
        frame.plot(x[mask],y[mask])
        frame.plot(xf, y_f)
        plt.show()#"#""
        
        # Prepare the next step
        if start+xrang*multi >= max(x):
            break
        start = start + 0.25*xrang*multi        # Any value
        start = x[x>=start][0]                  # The value of the first x value
    #data_ori = data
    # Combine the polynomial curves
    data = data[:,data[0,:].argsort()]
    step = int(10**(max(1,np.log10(data.shape[1])-2.5)))
    len_data = int(data.shape[1]/step)
    data = data[:,:len_data*step]
    data.shape = (3, len_data, step)            # Reorder the data to use np
    x = np.average( data[0,:,:], weights=data[2,:,:], axis=1 )
    y = np.average( data[1,:,:], weights=data[2,:,:], axis=1 )
    data = np.vstack((x,y))"""
    miny = np.nanmin(data[1,:])
    maxy = np.nanmax(data[1,:])
    y_values = np.linspace(miny, maxy, ysteps+1, endpoint=True)  # last entry is not necessary but for range
    #y_values = np.linspace(miny, maxy, ysteps, endpoint=False)  # without last entry
    
    # Find the y values
    data_g = np.zeros(0)
    data_g.shape = (2,0)
    for yy in y_values:
        diffsign = np.sign( data[1,:] - yy )
        signchange = (np.roll(diffsign, 1) - diffsign) != 0
        signchange[0] = False               # As last and first element are checked by np.roll
        signchangepos = np.where(signchange)[0]
        for ii in signchangepos:
            x_interp = linear_interpolation(data[1,ii-1:ii+1], data[0,ii-1:ii+1], yy)   # Interpolate, takes 2.5E-5s
            #f = scipy.interpolate.interp1d(data[1,ii-1:ii+1], data[0,ii-1:ii+1], kind='linear')     # together with next step takes 1.5E-4s
            #x_interp2= f(yy)                                                           # result is the same as with linear_interpolation
            data_g = np.hstack(( data_g, np.vstack((x_interp, yy)) ))

    data_g = data_g[:,data_g[0,:].argsort()]
    
    
    """fig, frame = plt.subplots(1,2)
    #frame[0].plot(x,y)
    #frame[0].plot(data_ori[0,:],data_ori[1,:])
    frame[0].plot(data[0,:],data[1,:])
    frame[0].plot(data_g[0,:],data_g[1,:])
    #frame[1].plot(y,x)
    #frame[1].plot(yy,xxf)
    plt.show()"""
        
    return data_g[0,:], data_g[1,:]

def find_overlap(stats):
    """
    Find the overlap-ordering of data with different data ranges
    :param stats: 2d array: min and max as columns, indexes as third column is optional
    :return indexes: list of list of integers: order of the indexes, with the best index starting, and then continuing
    stats: 0=min, 1=max, 2=global index, 3=index, 4=diff, 5=sorted, 6=mask
    """
    if stats.shape[0] == 0:
        return [[]]
    if stats.shape[1] == 2:
        stats = np.vstack(( stats.T, np.arange(stats.shape[0]) )).T                  # Add global index, if necessary
    stats = np.vstack(( stats.T, np.arange(stats.shape[0]) , stats[:,1] - stats[:,0], np.argsort(stats[:,1] - stats[:,0]), np.ones(stats.shape[0]) )).T     # Add index, range, sort by longest range, mask
    index = int(stats[-1,5])                            # index of the longest range
    min_val, max_val = stats[index,0], stats[index,1]    # first entry: min and max
    stats[index,6] = 0                                  # don't use the first entry again
    indexes = [[ int(stats[index,2]) ]]
    for step in np.arange(np.sum(stats[:,6])):
        sub_stats = stats[stats[:,6].astype(bool), :]   # Only entries which haven't been done yet
        sub_stats[:,4] = 0
        completely_coverd = ( ( sub_stats[:,0] >= min_val ) & ( sub_stats[:,1] <= max_val ) )
        if np.sum(completely_coverd) >= 1:              # Use the data points which have the biggest range
            sub_stats[completely_coverd,4] = sub_stats[completely_coverd,1] - sub_stats[completely_coverd,0]
            sub_stats[~completely_coverd,4] = 0
        else:
            sub_stats[:,4] = np.min( np.vstack(( sub_stats[:,1], [max_val]*sub_stats.shape[0])), axis=0 ) - np.max( np.vstack(( sub_stats[:,0], [min_val]*sub_stats.shape[0])), axis=0 )
            if np.sum( (sub_stats[:,4] > 0) ) == 0:
                indexes += find_overlap(sub_stats[:,:3]) # get the indexes for the sub_array
                break
        sub_stats[:,5] = np.argsort(sub_stats[:,4])
        sub_index = int(sub_stats[-1,5])                # Which position in the sub_array has the highest overlap
        index = int(sub_stats[sub_index,3])             # Which position in the stats (local) array
        min_val, max_val = np.min([min_val, stats[index,0]]), np.max([max_val, stats[index,1]])
        stats[index,6] = 0 
        indexes[0].append(int(stats[index,2]))          # Which position in the global stats
        #print( stats[index,0:3], stats[index,3:] , sub_stats[sub_index,0:3])
        
    return indexes

def calculate_bisector_multiple(data, data_x):
    """
    Combines the bisector measurement of different slices
    :param data: list of lists/arrays of floats, e.g. the flux of an absorption line. The feature needs to have lower flux than the continuum
    :param data_x: list of lists/arrays of floats, e.g. pixel or wavelength or velocity
    :return bisecs: 2d array (x,y as rows): All datapoints from the bisector mesurements, sorted by y
    :return bisec_fit: 2d array (x,y as rows): Fitted bisectors using a high order polynomial, sorted by y
    :return lines: 2d array (x,y as rows): The absorption lines, sorted by x
    """
    if len(data) == 0:
        return np.zeros(0), np.zeros(0), np.zeros(0)             # Empty array
    if type(data[0]).__name__ not in ['list','ndarray']:
        data = [data]
        data_x = [data_x]
    # Fold the data into a one-D curve by rescaling and fitting the centre
    bisecs_sub, line_data = [], []
    bisecs_stat = np.zeros(0)
    bisecs_stat.shape = (3,0)
    for ii in range(len(data)):
        bisec, norm_line = calculate_bisector(data[ii], data_x[ii])        # 2d array with x and y
        if bisec.shape[1] > 3:
            bisecs_sub.append(bisec)
            bisecs_stat = np.hstack(( bisecs_stat, np.vstack(( np.min(bisec[1,:]), np.max(bisec[1,:]), bisecs_stat.shape[1] )) ))     # min, max, index
            line_data.append( norm_line )
    # Combine the bisectors from each measurement using overlapping areas
    #print('find overlapping indexes')
    if bisecs_stat.shape[1] == 0:
        return np.zeros(0), np.zeros(0), np.zeros(0)             # Empty array
    indexes = find_overlap(bisecs_stat.T)
    #print('finished', indexes)
    lines = np.zeros(0)
    lines.shape = (2,0)
    bisecs = np.zeros(0)
    bisecs.shape = (3,0)
    for ii in range(len(indexes)):
        line_part = line_data[indexes[ii][0]]
        bisecs_part = np.vstack(( bisecs_sub[indexes[ii][0]], np.arange(bisecs_sub[indexes[ii][0]].shape[1]) ))
        for index in indexes[ii][1:]:
            for jj in range(bisecs_part.shape[1]-1):                # Avoid having the same y values
                if bisecs_part[1,jj] == bisecs_part[1,jj+1]:
                    bisecs_part[1,jj+1] *= (1+1E-9)
                    #print jj, bisecs_part[1,jj:jj+2]
            #f = scipy.interpolate.UnivariateSpline(bisecs_part[1,:], bisecs_part[0,:], s=0)       # y becomes x in this setting, this is no linear interpolation and does bad things
            f = scipy.interpolate.interp1d(bisecs_part[1,:], bisecs_part[0,:])                      # y becomes x in this setting
            bisec_new = bisecs_sub[index]
            mask = ( ( bisec_new[1,:] >= np.min(bisecs_part[1,:]) ) & ( bisec_new[1,:] <= np.max(bisecs_part[1,:]) ) )        # Only use values in the old array
            offset = np.median( f(bisec_new[1,mask]) - bisec_new[0,mask] )                            # calculate the difference in x
            """print ii, index, offset, np.median( bisecs_part[0,:]) - np.median( bisec_new[0,:] )
            fig, frame = plt.subplots(1,1)
            y_temp = np.arange(np.min(bisecs_part[1,:]), np.max(bisecs_part[1,:]), 0.001)
            frame.plot( f(y_temp), y_temp, label='Spline')
            frame.plot(bisecs[0,:], bisecs[1,:], label='bisecs')
            frame.plot(bisecs_part[0,:], bisecs_part[1,:], label='bisecs_part')
            frame.plot(bisec_new[0,:], bisec_new[1,:], label='new before')"""
            bisec_new[0,:] += offset
            """frame.plot(bisec_new[0,:], bisec_new[1,:], label='new after')
            frame.plot(line_part[0,:], line_part[1,:], label='line_part')
            frame.plot(line_data[index][0,:], line_data[index][1,:], label='line before')
            frame.legend(loc='upper right')
            plt.show()         # This shows that the result with higher scatter make sense"""
            bisecs_part = np.hstack(( bisecs_part, np.vstack(( bisec_new, np.arange(bisec_new.shape[1])+bisecs_part.shape[1] )) ))
            bisecs_part = bisecs_part[:,np.argsort(bisecs_part[1,:])]
            line_data[index][0,:] += offset
            line_part = np.hstack(( line_part, line_data[index] ))
        if bisecs.shape[1] > 0:
            #offset = np.median( bisecs[0,:] - bisecs_part[0,:] )                       # Before 20200816
            offset = np.median( bisecs[0,:] ) - np.median( bisecs_part[0,:] )
            bisecs_part[0,:] += offset
            line_part[0,:] += offset
        bisecs = np.hstack(( bisecs, bisecs_part ))
        lines = np.hstack(( lines, line_part ))
    bisecs = bisecs[0:2,np.argsort(bisecs[1,:])]        # sort by y
    lines = lines[:,np.argsort(lines[0,:])]        # sort by x
    # Fit a polynom, find the unvertainties
    nc = max(0,min(8,min(np.unique(bisecs[1,:]).shape[0], np.unique(bisecs[0,:]).shape[0])-2))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', np.RankWarning)
        poly = np.polyfit(bisecs[1,:], bisecs[0,:], nc)
    y = np.linspace(np.min(bisecs[1,:]), np.max(bisecs[1,:]), min(bisecs.shape[1]-2, max(10,int(0.1*bisecs.shape[1]))) )
    x = np.polyval(poly, y)
    med_x = np.median(x)
    diff_fit = np.polyval(poly, bisecs[1,:]) - bisecs[0,:]      # difference from fit
    maskl, maskr = (diff_fit < 0), (diff_fit > 0)               # left and right side
    x -= med_x
    bisecs[0,:] -= med_x
    lines[0,:] -= med_x
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', np.RankWarning)
        polyl = np.polyfit(bisecs[1,maskl], bisecs[0,maskl], nc)
        polyr = np.polyfit(bisecs[1,maskr], bisecs[0,maskr], nc)
    xl = np.polyval(polyl, y)
    xr = np.polyval(polyr, y)
    bisec_fit = np.vstack(( x, y, xl, xr ))
    
    """fig, frame = plt.subplots(1,1)
    frame.plot(bisecs[0,:], bisecs[1,:])
    frame.plot(bisec[0,:], bisec[1,:])
    plt.show()"""
    
    return bisecs, bisec_fit, lines
    
def calculate_bisector(y, x):
    """
    Calculates the bisector in an array of data given by interpolating the normalised data into a grid of y-values
    :params y: 1d list or array of floats: y-values of an absorption line
    :params x: 1d list or array of floats: x-values for the y values
    :return bisec: 2d array with x and y of the bisector (y=[1,:])
    :return (x, ys): 2d array with x and y of the normalised lines
    """
    x, y = np.array(x), np.array(y)
    # Rescale to between 1 (continuum) and 0 (absortion)
    ys, popt = scale_data(y, x, mode='poly3')       #mode='gauss')
    # Create a grid of the y data
    diffy = np.abs(ys[:-1]-ys[1:])          # Spacing between consecutive data points
    diffy = diffy[diffy > 0]                # consecutive values could have 0 at continuum level
    ysteps = int(10./np.nanpercentile(diffy, 90))        # 5 steps per step
    xg, yg = create_grid_data(ys, x, ysteps=ysteps)
    
    # Find for each yy the points left and right of the absorption line
    indmin = np.argmin(yg)
    x = x - xg[indmin]
    xg -= xg[indmin]
    mask = ( ( yg >= 0.02 ) & ( yg <= 0.95 ) )
    yys = np.unique(yg[mask])
    bisec = np.vstack(( yys*np.nan, yys ))
    #print y,x,ys, min(ys), max(ys), yg,xg,ysteps, xrang, indmin, xg[indmin], popt
    #print min(ys), max(ys), ysteps, xrang, indmin, xg[indmin], popt
    for ii in range(bisec.shape[1]):
        index = np.where( yg == bisec[1,ii] )[0]        # which indexes have the correct y value
        maskl = (index <= indmin)
        maskr = (index >= indmin)
        if np.sum(maskl) == 0 or np.sum(maskr) == 0:
            continue
        xl = xg[index[maskl][-1]]
        xr = xg[index[maskr][0]]
        bisec[0,ii] = 0.5*(xr+xl)     # Centre, corrected by the Gaussian centre
    """print bisec
    fig, frame = plt.subplots(1,1)
    frame.plot(x,ys)
    frame.plot(xg,yg)
    frame.plot(bisec[0,:], bisec[1,:])
    plt.show()#"""
    
    return bisec[:,~np.isnan(bisec[0,:])], np.vstack(( x, ys ))

def bisector_measurements_orders(im, fname, pfits, xlows, xhighs, widths=[]):
    """
    Measures the bisector at different areas of the image and plots the data
    """
    no = len(pfits)
    no2 = int(no/2)
    no4 = int(no/4)
    if no >= 20:
        orders_to_check = [[0,1,2], [no4-1,no4,no4+1], [no2-1,no2,no2+1], [3*no4-1,3*no4,3*no4+1], [no-3,no-2,no-1]]    # Modify by user, if wanted
        #orders_to_check = [[0], [no4], [no2], [3*no4], [no-1]]
    elif no >= 8:
        orders_to_check = [[0,1,2], [no2-1,no2,no2+1], [no-3,no-2,no-1]]    # Modify by user, if wanted
    elif no >= 3:
        orders_to_check = [[0], [no2], [no-1]]
    else:
        orders_to_check = [[0], [no-1]]
    number_of_positions = 7                         # Modify by user, if wanted
    number_of_pixel = 20                            # Modify by user, if wanted
    x_label = 'Position relative to centre [px], (order cross-sections (blue) were rescaled to 10% in x)'
    y_label = 'Flux (normalised)'
    xmin = min(xlows[ np.array(orders_to_check).flatten() ])
    xmax = max(xhighs[ np.array(orders_to_check).flatten() ])
    if number_of_positions == 1:
        xposis = np.array( [0.5*(xmax-number_of_pixel + xmin)], dtype=int )
    else:
        xposis = np.linspace(xmin, xmax-number_of_pixel, number_of_positions, endpoint=True, dtype=int)
    ims = im.shape
    leftfit, rightfit = [], []
    if len(pfits.shape) == 3:                           # new 3d array
        leftfit = pfits[:,1,:]
        rightfit = pfits[:,2,:]
        pfits = pfits[:,0,:]
        
    x_range = np.zeros(0)
    lx, ly = len(xposis), len(orders_to_check)
    fig, frame = plt.subplots(ly, lx, sharex=True, sharey=True, figsize=(2*lx, 2*ly) )
    plt.subplots_adjust(left=max(0.04,0.19-0.02*lx), right=min(0.98,0.90+0.01*lx), top=min(0.96,0.92+.0067*ly), bottom=max(0.04,0.17-0.023*ly), wspace=0., hspace=0.)
    for ii, order_range in enumerate(orders_to_check):
        for jj, posi in enumerate(xposis):
            imslices = []
            imslices_x = []
            for order in order_range:
                xarr = np.arange( max(posi,xlows[order]), min(posi+number_of_pixel,xhighs[order]), 1)
                if len(xarr) == 0:
                    continue
                yarr = np.polyval(pfits[order,1:], xarr-pfits[order,0])
                if len(widths) > 0:
                    yarrl = yarr-widths[order,0]
                    yarrr = yarr+widths[order,1]
                else:
                    yarrl = np.polyval(leftfit[order, 1:], xarr - leftfit[order, 0])
                    yarrr = np.polyval(rightfit[order, 1:], xarr- rightfit[order, 0])
                    yarrl, yarrr = adjust_width_orders(yarr, yarrl, yarrr, [1.5, 1.5])              # Adjust width, double it to more likely catch the whole trace, but too big and problem with MRES red orders
                for kk,x in enumerate(xarr):
                    imslices.append( -1*im[x,max(0,int(yarrl[kk])):min(ims[1],int(yarrr[kk]+.99))+1] )    # negative to create absorption line
                    imslices_x.append( range(len(imslices[-1])) )
            bisecs, bisec_fit, lines = calculate_bisector_multiple(imslices, imslices_x)
            # Plot
            if lx == 1 and ly == 1: frame_temp = frame
            elif ly == 1: frame_temp = frame[jj]
            elif lx == 1:          frame_temp = frame[ii]
            else:                           frame_temp = frame[ii,jj]
            if len(bisecs.shape) != 1:                                              # only when data
                #frame_temp.plot(bisecs[0,:], 1-bisecs[1,:], color='lightgray')   # 1-y to revese the creation of absorption line
                x = np.hstack(( bisec_fit[2,:], bisec_fit[3,:]))
                y = np.hstack(( 1-bisec_fit[1,:], 1-bisec_fit[1,:] ))
                sort_index = np.argsort(y)
                frame_temp.plot(x[sort_index], y[sort_index], color='lightgray', linewidth=3)
                frame_temp.plot(0.1*lines[0,:], 1-lines[1,:], color='blue', marker='o', markersize=1, linestyle='')       #decrese the x-scaling by factor 10
                frame_temp.plot(bisec_fit[0,:], 1-bisec_fit[1,:], color='k')   # 1-y to revese the creation of absorption line
                x_range = np.hstack(( x_range, bisec_fit[0,:] ))
            if ii == 0:
                secax = frame_temp.twiny()
                secax.axes.xaxis.set_ticklabels([])     # Disable the tick label
                secax.tick_params( axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
                secax.set_xlabel('Pixels: {0} - {1}'.format(posi, posi+number_of_pixel-1), fontsize=10)
                secax.xaxis.set_label_position("top")
            if jj == len(xposis)-1:
                secax = frame_temp.twinx()
                secax.axes.yaxis.set_ticklabels([])     # Disable the tick label
                secax.tick_params( axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
                secax.set_ylabel('Orders: {0}'.format(str(order_range)[1:-1]), fontsize=10)
                secax.yaxis.set_label_position("right")
            frame_temp.tick_params( axis='both', which='both', bottom=True, top=True, right=True, left=True, direction='inout' )           # Have the ticks in all frames
    
    if lx == 1 and ly == 1: f1, f2 = frame, frame
    elif ly == 1:   f1, f2 = frame[int(jj/2)] , frame[0]
    elif lx == 1:   f1, f2 = frame[-1], frame[int(ii/2)]
    else:           f1, f2 = frame[-1,int(jj/2)], frame[int(ii/2),0]
    f1.set_xlabel(x_label, fontsize=13)
    f2.set_ylabel(y_label, fontsize=13)
    x_min, x_max = np.percentile(x_range,2), np.percentile(x_range,98)
    x_range = x_max - x_min
    f1.set_xlim(x_min-0.2*x_range,x_max+0.2*x_range)
    f1.set_ylim(-0.1,1.1)
    plt.savefig(fname, bbox_inches='tight')
    plt.close()

def get_subframe(frame_list, nr_h, nr_v, ih, iv):
    if type(frame_list).__name__ == 'ndarray':
        if len(frame_list.shape) == 2:
            subframe = frame_list[iv,ih]
            fx = frame_list[-1,int(nr_h/2)]
            fy = frame_list[int(nr_v/2),0]
        elif nr_v == 1:
            subframe = frame_list[ih]
            fx = frame_list[int(nr_h/2)]
            fy = frame_list[0]
        elif nr_h == 1:
            subframe = frame_list[ivh]
            fy = frame_list[-1]
            fy = frame_list[int(nr_v/2)]
    else:
        subframe = frame_list
        fx = frame_list
        fy = frame_list

    return subframe, fx, fy

def bisector_measurements_emission_lines(fname, spec, data, indexes):     #im, fname, pfits, xlows, xhighs, widths=[]):
    """
    Measures the bisector for the emission lines used for the wavelength solution and plot it
    """
    ord_split = 5                       # Modify by user, if wanted
    px_split = 5                        # Modify by user, if wanted
    
    [oo, xx, gg] = indexes
    
    min_orders, max_orders = np.min(data[:,oo]), np.max(data[:,oo])+1E-6
    min_px, max_px = np.min(data[:,xx]), np.max(data[:,xx])+1E-6
    range_orders = (max_orders - min_orders + 0.0)/ord_split
    range_px = (max_px - min_px + 0.0)/px_split
    max_width = int(np.ceil(np.max(data[:,gg])*2.35482*2))        # 2*FWHM to either side
    
    x_range = np.zeros(0)
    x_label = 'Position relative to centre [px], (order cross-sections (blue) were rescaled to 10% in x)'
    y_label = 'Flux (normalised)'
    fig, frame = plt.subplots(px_split, ord_split, sharex=True, sharey=True, figsize=(2*px_split, 2*ord_split) )
    plt.subplots_adjust(left=max(0.04,0.19-0.02*px_split), right=min(0.98,0.90+0.01*px_split), top=min(0.96,0.92+.0067*ord_split), bottom=max(0.04,0.17-0.023*ord_split), wspace=0., hspace=0.)
    for ii in range(ord_split):
        for jj in range(px_split):
            mino, maxo = round(ii*range_orders+min_orders,0), round((ii+1)*range_orders+min_orders,0)
            minx, maxx = jj*range_px+min_px, (jj+1)*range_px+min_px
            good_data = ( (data[:,oo] >= mino) & (data[:,oo] < maxo) & (data[:,xx] >= minx) & (data[:,xx] < maxx) )
            imslices = []
            imslices_x = []
            for entry in data[good_data,:]:
                pos = int(round(entry[xx]))
                #imslices_x.append( range( max(0,pos-max_width), min(pos+max_width+1,spec.shape[1]) ) )
                #imslices.append( spec[int(entry[oo]),imslices_x[-1]] )
                spec_sub = spec[int(entry[oo]), max(0,pos-max_width): min(pos+max_width+1,spec.shape[1])]
                spec_sub = spec_sub[~np.isnan(spec_sub)]
                imslices.append( -1*spec_sub )
                imslices_x.append( range(spec_sub.shape[0]) )
                
            bisecs, bisec_fit, lines = calculate_bisector_multiple(imslices, imslices_x)
            # Plot
            if px_split == 1 and ord_split == 1: frame_temp = frame
            elif ord_split == 1:    frame_temp = frame[jj]
            elif px_split == 1:     frame_temp = frame[ii]
            else:                   frame_temp = frame[ii,jj]
            if len(bisecs.shape) != 1:                                              # only when data
                x = np.hstack(( bisec_fit[2,:], bisec_fit[3,:]))
                y = np.hstack(( 1-bisec_fit[1,:], 1-bisec_fit[1,:] ))
                sort_index = np.argsort(y)
                frame_temp.plot(x[sort_index], y[sort_index], color='lightgray', linewidth=3)
                frame_temp.plot(0.1*lines[0,:], 1-lines[1,:], color='blue', marker='o', markersize=1, linestyle='')       #decrese the x-scaling by factor 10
                frame_temp.plot(bisec_fit[0,:], 1-bisec_fit[1,:], color='k')   # 1-y to revese the creation of absorption line
                frame_temp.text( 0, 0, '{0} lines'.format(len(imslices)), 
                                    horizontalalignment='center', verticalalignment='bottom', rotation=0, color='r', zorder=5 )
                x_range = np.hstack(( x_range, bisec_fit[0,:] ))
            if ii == 0:
                secax = frame_temp.twiny()
                secax.axes.xaxis.set_ticklabels([])     # Disable the tick label
                secax.tick_params( axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
                secax.set_xlabel('Pixels: {0} - {1}'.format(int(round(minx)), int(round(maxx))), fontsize=10)
                secax.xaxis.set_label_position("top")
            if jj == px_split-1:
                secax = frame_temp.twinx()
                secax.axes.yaxis.set_ticklabels([])     # Disable the tick label
                secax.tick_params( axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
                secax.set_ylabel('Orders: {0} - {1}'.format(int(mino),int(maxo)), fontsize=10)
                secax.yaxis.set_label_position("right")
            frame_temp.tick_params( axis='both', which='both', bottom=True, top=True, right=True, left=True, direction='inout' )           # Have the ticks in all frames
    
    if px_split == 1 and ord_split == 1: f1, f2 = frame, frame
    elif ord_split == 1:    f1, f2 = frame[int(jj/2)] , frame[0]
    elif px_split == 1:     f1, f2 = frame[-1], frame[int(ii/2)]
    else:                   f1, f2 = frame[-1,int(jj/2)], frame[int(ii/2),0]
    f1.set_xlabel(x_label, fontsize=13)
    f2.set_ylabel(y_label, fontsize=13)
    x_min, x_max = np.percentile(x_range,2), np.percentile(x_range,98)
    x_range = x_max - x_min
    f1.set_xlim(min(x_min-0.2*x_range, -0.1*max_width/2),max(x_max+0.2*x_range, 0.1*max_width/2))
    f1.set_ylim(-0.1,1.1)
    plt.savefig(fname, bbox_inches='tight')
    plt.close()
    
    #print('exit 4936')
    #os.system('notify-send -t 2000 "done"')
    #exit()
    
def plot_wavelength_solution_form(fname, traces_def, wavelength_solution):
    """
    Creates a map of the wavelength solution
    :param fname: Filename to which the image is saved
    :param xlows: list, length same as number of orders, the lowest x pixel (wavelength direction) used in each order
    :param xhighs: list, length same as number of orders, the highest x pixel (wavelength direction) used in each order
    :param wavelength_solution: 2d array of floats, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    """
    step = 25
    [tr_poly, xlows, xhighs, widths] = traces_def
    im = np.zeros([int(max(xhighs)/step), len(xlows)])+np.nan
    text_order = ''
    for order in range(len(wavelength_solution)):
        xarr = np.arange(xlows[order], xhighs[order], step)
        yarr = np.polyval(wavelength_solution[order,2:], xarr-wavelength_solution[order,1])
        diff_wave = yarr[1:] - yarr[:-1]
        im[(xarr[:-1]/step).astype(int), order ] = diff_wave        # im[{floats} will work
        diff_wave2 = diff_wave[1:] - diff_wave[:-1]
        if np.sum(diff_wave2>1E-3) > 10 and np.sum(diff_wave2<-1E-3) > 10:     # There should'nt be any stationary points (10 points are ok, a change of less than 1E-3 as well)
            text_order += '{0},'.format(order)
    if text_order != '':
        logger('Warn: The wavelength solution for oder(s) {0} contain(s) at least one stationary point. Please check {1}'.format(text_order[:-1], fname))
    colorbar = True  
    title = 'Plot the wavelength difference between every {0} pixel (Angstrom/{0}px)'.format(step)
    axis_name = ['Order', 'Position in Dispersion direction [{0} px]'.format(step), 'wavelength difference [Angstrom/{0}px]'.format(step)]
    plot_img_spec.plot_image(im, [fname], pctile=0, show=False, adjust=[0.05,0.95,0.95,0.05], title=title, autotranspose=False, colorbar=colorbar, axis_name=axis_name, size=[16,10])

def plot_wavelength_solution_width_emmission_lines(fname, specs, data, indexes, step=20, title='Gaussian width'):
    """
    Creates a map of the Gaussian width of the emmission lines
    :param fname: Filename to which the image is saved
    :param specs: 1d list with two integers: shape of the extracted spectrum, first entry gives the number of orders and the second the number of pixel
    :param data: 2d array of floats with the data to plot
    :param indexes: list of int, position of the relevant columns: order, wavelength, gauss
    """
    [oo, xx, gg] = indexes
    
    gauss_ima = np.empty((specs[0], int(specs[1]/step)+1)) * np.nan
    for order in range(specs[0]):
        for i,px in enumerate(range(0, specs[1], step)):
            good_values = (data[:,oo] == order) & (data[:,xx] >= px) & (data[:,xx] < px+step)
            if np.sum(good_values) > 0:
                gauss_ima[order,i] = np.median(data[good_values,gg])
    gauss_ima[np.isnan(gauss_ima)] = np.nanmax(gauss_ima)+1
    title_f = 'Plot of the {1} every {0} pixel (white shows areas with no data available)'.format(step, title)
    axis_name = ['Position in Dispersion direction [{0} px]'.format(step), 'Order', '{0} of the emission lines [px] (white shows areas with no data available)'.format(title)]
    plot_img_spec.plot_image(gauss_ima, [fname], pctile=0, show=False, adjust=[0.05,0.95,0.95,0.05], title=title_f, autotranspose=False, colorbar=True, axis_name=axis_name, size=[16,10])

def plot_overlapping_orders(order, x_full, y1_full, y2_full=None, labels=[]):
    """
    :param order: integer, first index in x_full, y1_full, y2_full
    :param x_full: 2d array of floats: first index for order, second index for pixel; normally wavelength or pixel
    :param y1_full: same format and size as x_full, flux
    :param y2_full: same format and size as x_full, flux of a different exposure time
    :param labels: list of strings, label for y1_full and y2_full
    """
    x_data = [x_full[order,:]]
    y_data = [y1_full[order,:]]
    color = ['tab:blue']
    if y2_full is not None:
        x_data.append(x_full[order,:])
        y_data.append(y2_full[order,:])
        color.append('tab:orange')
    
    label = copy.deepcopy(labels)
    if np.min(x_full) > 1000:      # real wavelength solution, not pseudo solution
        if order < x_full.shape[0]-1:
            x_data.insert(0, x_full[order+1,:])
            y_data.insert(0, y1_full[order+1,:])
            label.insert(0, 'next order')
            color.insert(0, 'lightgrey')
        if order > 0:
            x_data.insert(0, x_full[order-1,:])
            y_data.insert(0, y1_full[order-1,:])
            label.insert(0, 'prev. order')
            color.insert(0, 'silver')
    return x_data, y_data, label, color

def adjust_data_log(data):
    """
    Calculates the log10 of the data, adds offsets before to avoid negative values
    input data: ndarray or list
    return data: ndarray: log10 of input data
    """
    if np.nanmin(data) < 1:
        data += 1 - np.nanmin(data)       # minus to cancel negative minimum
    data = np.log10(data)
    
    return data

def plot_wavelength_solution_spectrum(params, spec1, spec2, fname, wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, plot_log=False):
    """
    Creates a pdf with one page for each ThAr spectrum and adds the identified and further catalogue lines to the plots
    :param spec1: 2d array of floats, extracted spectrum of the long exposed arc
    :param spec2: 2d array of floats, extracted spectrum of the short exposed arc
    :param fname: string, Filename to which the image is saved
    :param wavelength_solution: 2d array of floats, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    :param wavelength_solution_arclines: 2d array of floats with one line for each order. Each order contains the wavelengths of the identified reference lines, 
                                        sorted by the brightest to faintest (in spectrum from which the solution is derived) and 0 to make into an array
    :param reference_catalog: 2d array or list of floats with one entry for each line. Each entry contains the wavelength and the intensity of the line (set to 1 if not provided by the catalog file), and the index in the file
    :param reference_names: list of strings, same length as reference_catalog. Names of the reference lines
    """
    max_reflines = 30
    plot_pos = [0.01, 0.035, 0.04]      # begin line, end line, beginn of text; all relative to the the flux at the position and proportional to the range of the flux
    
    labels = ['long exp','short exp']
    x_title = 'Wavelength [Angstroms]'
    y_title = 'extracted flux [ADU]'
    titel_f = 'Order {0}, real Order {1}\nmarked (proportional to line strength) are the identified reference lines (red, {2} lines) and\na subset of the omitted reference lines (green, showing the brightes {3} out of {4} lines)'
    
    if plot_log:
        y_title = 'log10 ( {0} )'.format(y_title)
        spec1 = adjust_data_log(spec1)
        spec2 = adjust_data_log(spec2)
    reference_catalog = np.array(sorted(reference_catalog, key=operator.itemgetter(1), reverse=True ))          # Sort by intensity in reverse order
    wavelengths, dummy = create_wavelengths_from_solution(params, wavelength_solution, spec1)
    # multiple pdf pages from https://matplotlib.org/examples/pylab_examples/multipage_pdf.html
    with PdfPages(fname) as pdf:
        for order in range(wavelength_solution.shape[0]):
            x_range = wavelengths[ order, ~np.isnan(spec1[order,:]) ]
            if len(x_range) < 10:                                       # the reference trace cod fall completely out of the CCD
                continue
            good_values = np.where((reference_catalog[:,0] >= min(x_range)) & (reference_catalog[:,0] <= max(x_range)) )[0]     # lines in the right wavelength range
            reference_lines = reference_catalog[good_values,:]                                                                  # possible reference lines in the order
            
            num_ident = len( np.where(np.array(wavelength_solution_arclines[order]) > 0)[0] )                                            # np.where only works on arrays
            #if num_ident == 0:                                          # Add a dummy value in case no lines are available for this order, as otherwise it will crash ; !! Not necessary after making everything into an numpy array?
            #    wavelength_solution_arclines[order].append(-1E5)
            num_notident = len(reference_lines)
            
            title = titel_f.format(order, int(wavelength_solution[order,0]), num_ident, min(max_reflines,num_notident-num_ident), num_notident-num_ident)
            fig, frame = plt.subplots(1, 1)
            
            # Create a the plotting data
            x_data, y_data, label, color = plot_overlapping_orders(order, wavelengths, spec1, spec2, labels)        # create the overlapping pixel
            plot_img_spec.plot_points(x_data, y_data, label, [], show=False, adjust=[0.12,0.88,0.85,0.12, 1.0,1.01], title=title, 
                                      return_frame=True, frame=frame, x_title=x_title, y_title=y_title, linestyle="-", marker="",
                                      color=color)
            axes = plt.gca()
            #ymin, ymax = axes.get_ylim()
            #y_range = ymax - ymin
            #xmin, xmax = axes.get_xlim()        # wavelengths contains also values for nans
            # Redifine the axis using only the current order
            pctl = 5./len(spec1[order,:])*100                   # Exclude the five highest and five lowest points
            xmin, xmax = np.nanmin(wavelengths[order,:]), np.nanmax(wavelengths[order,:])
            ymin1, ymax1 = np.nanpercentile(spec1[order,:], pctl), np.nanpercentile(spec1[order,:], 100-pctl)
            ymin2, ymax2 = np.nanpercentile(spec2[order,:], pctl), np.nanpercentile(spec2[order,:], 100-pctl)
            ymin, ymax = min(ymin1, ymin2), max(ymax1, ymax2)
            ymin, ymax = min(ymin1, ymin2), max(ymax1, ymax2)
            y_range = ymax - ymin
            x_range = xmax - xmin
            
            if len(reference_lines) > 0:
                y_scale = (plot_pos[1] - plot_pos[0]) / (max(reference_lines[:,1])+0.0)
                num_notident = 1
                for color in ['g', 'r']:        # plot the green lines before the red ones
                    for refline in reference_lines:
                        if np.min(np.abs(refline[0] - wavelength_solution_arclines[order])) < 1E-2:
                            if color != 'r':
                                continue        # don't plot a matching line when it's not red lines to be plotted
                        else:
                            if color == 'r':
                                continue        # don't plot the non matching lines when red lines to be plotted
                            if num_notident > max_reflines:
                                break           # stop plot non-matching lines
                            num_notident += 1   
                        x_position = np.argmin(np.abs(refline[0] - wavelengths[order,:]))
                        if np.isnan(spec1[order, x_position]):
                            continue
                        y_position = np.nanmax(spec1[order, max(0,x_position-2):min(len(spec1[order,:]),x_position+2)])
                        ymax = max(ymax, y_position+0.23*y_range)
                        frame.plot( [refline[0],refline[0]], [y_position+plot_pos[0]*y_range,y_position+(plot_pos[0]+y_scale*refline[1])*y_range], color=color )
                        frame.text( refline[0], y_position+plot_pos[2]*y_range, '{0} - {1}'.format(round(refline[0],4),reference_names[int(refline[2])]), 
                                    horizontalalignment='center', verticalalignment='bottom', rotation=90, color=color, zorder=5 )
                    
            axes.set_ylim(ymin,ymax)
            axes.set_xlim(xmin-0.01*x_range,xmax+0.01*x_range)
            fig.set_size_inches(16.2, 10)
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()
    
        # We can also set the file's metadata via the PdfPages object:
        d = pdf.infodict()
        d['Title'] = 'Spectral atlas'
        d['Author'] = 'HiFLEx pipeline'
        d['Subject'] = ''
        d['Keywords'] = 'Spectra atlas HiFLEx pipeline'
        d['CreationDate'] = datetime.datetime.today()
        d['ModDate'] = datetime.datetime.today()
    
def plot_wavelength_solution_image(im, fname, traces_def, wavelength_solution, wavelength_solution_arclines, reference_catalog_full):
    """
    Plots the identified and used reference lines and all available reference lines to the image
    :param im: Image of the emission lines
    :param fname: Filename to which the image is saved
    :param pfits: list of floats, length same as number of orders, polynomial values
    :param xlows: list, length same as number of orders, the lowest x pixel (wavelength direction) used in each order
    :param xhighs: list, length same as number of orders, the highest x pixel (wavelength direction) used in each order
    :param wavelength_solution: 2d array of floats, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    :param wavelength_solution_arclines: 2d array of floats with one line for each order. Each order contains the wavelengths of the identified reference lines, 
                                        sorted by the brightest to faintest (in spectrum from which the solution is derived) and 0 to make into an array
    :param reference_catalog: 2d array of floats with one entry for each line. Each entry contains the wavelength and the intensity of the line
    """
    logger('Step: logging results to {0}'.format(fname))
    [pfits, xlows, xhighs, widths] = traces_def
    # Prepare the refernece catalog to show only limit number of lines
    #wavelength_solution_arclines_flattened = [item for sublist in wavelength_solution_arclines for item in sublist]
    wavelength_solution_arclines_flattened = wavelength_solution_arclines.reshape( ( np.prod(wavelength_solution_arclines.shape), ) )
    orders = list(range(pfits.shape[0]))
    minmaxwav = np.zeros(( len(orders), 2 ))
    for order in orders:                                # need to check all orders, as with bed wavelength solution, the maximum can end up beeing smaller than the minimum
        minmaxwav[order,0] = np.polyval(wavelength_solution[order,2:], xlows[order]-wavelength_solution[order,1])                   # minimum wavelength in the order
        minmaxwav[order,1] = np.polyval(wavelength_solution[order,2:], xhighs[order]-wavelength_solution[order,1])                  # maximum wavelength on the chip
    min_wav = np.min(minmaxwav[:,0])
    max_wav = np.max(minmaxwav[:,1])
    d_wav = 0.02*(max_wav - min_wav) 
    good_values = np.where( (wavelength_solution_arclines_flattened >= min_wav) & (wavelength_solution_arclines_flattened <= max_wav))[0]
    wavelength_solution_arclines_flattened = wavelength_solution_arclines_flattened[good_values]                        # exclude the entries only added to fill the array
    wavelength_solution_arclines_uniq = np.unique(wavelength_solution_arclines_flattened)                               # unique lines 
    good_values = np.where((reference_catalog_full[:,0] >= min_wav -d_wav) & (reference_catalog_full[:,0] <= max_wav+d_wav) )[0]     # lines in the right wavelength regime
    reference_catalog_full = reference_catalog_full[good_values,:]                                                      # lines in the right wavelength regime
    wsaf = np.repeat( [wavelength_solution_arclines_uniq], len(reference_catalog_full), axis=0)        # array with wavelengths from the solution, duplicated for each reference line. shape = (num_reflines, num_solution)
    rcf = np.repeat( [reference_catalog_full[:,0]], len(wavelength_solution_arclines_uniq), axis=0 ).T # array with wavelengths from the reference, duplicated for each identified line. shape = (num_reflines, num_solution)
    diff = np.abs(rcf - wsaf)
    reference_catalog = reference_catalog_full[np.min(diff, axis=1) > 0.0001,:]                                         # Non-identified lines, even when rounding problems occur
    reference_catalog_sorted = np.array(sorted(reference_catalog, key=operator.itemgetter(1), reverse=True ))           # Sort the reference catalogue (only non-identified) by brightness of lines
    
    # Plot the data
    im = scale_image_plot(im,'log10')
    ims = im.shape
    fig, frame = plt.subplots(1, 1)
    colorbar = True  
    title = 'Plot of the identified and used reference lines (red, {2} lines) and a subset of the omitted reference lines (green, showing the brightes {0} out of {1} lines) to the image (log10 of image)'
    plot_img_spec.plot_image(im, [], pctile=1, show=False, adjust=[0.05,0.95,0.95,0.05], title='', return_frame=True, frame=frame, autotranspose=False, colorbar=colorbar)
    axes = plt.gca()
    ymin, ymax = axes.get_ylim()
    xmin, xmax = axes.get_xlim()
    num_non_ident = 0
    for pp in orders:
            xarr_fine = np.arange(xlows[pp], xhighs[pp], 0.2)
            yarr_fine_p5 = np.polyval(pfits[pp, 0, 1:], xarr_fine-pfits[pp, 0, 0]) + 5
            yarr_fine_m5 = np.polyval(pfits[pp, 0, 1:], xarr_fine-pfits[pp, 0, 0]) - 5
            # wavelength_solution_arclines contains the catalog wavelength -> make them into px
            xarr_wave = np.polyval(wavelength_solution[pp,2:], xarr_fine-wavelength_solution[pp,1])
            min_wave, max_wave = min(xarr_wave), max(xarr_wave)
            for source in [1,2]:
                if source == 1:
                    source_data = reference_catalog_sorted[:,0]
                    color = 'g'
                elif source == 2:
                    source_data = wavelength_solution_arclines[pp]
                    color = 'r'
                good_values = np.where( (source_data >= min_wave) & (source_data <= max_wave))[0]
                source_data = source_data[good_values]                                                      # only have the data for this order
                if source == 1:
                    source_data = source_data[:min(70,len(source_data))]                                        # only maximum 70 lines per order
                    num_non_ident += len(source_data)
                for arc_line_wave in source_data:
                    min_diff = (np.abs(arc_line_wave - xarr_wave)).argmin()
                    xarc = xarr_fine[min_diff]
                    yarc_p5 = yarr_fine_p5[min_diff]
                    yarc_m5 = yarr_fine_m5[min_diff]
                    #print pp, arc_line_wave, xarc, yarc        # Can be used to get the positions of the lines on the CCD
                    frame.plot( [yarc_m5,yarc_p5],[xarc,xarc],color=color )
                    frame.text(yarc_p5, xarc, '{0}'.format(round(arc_line_wave,2)), horizontalalignment='left', verticalalignment='center', rotation=0, color=color, zorder=5)
    axes.set_ylim(ymin,ymax)
    axes.set_xlim(xmin,xmax)
    frame.set_title(title.format(num_non_ident, len(reference_catalog_full), len(wavelength_solution_arclines_flattened) ), fontsize=16)
    fig.set_size_inches(42, 42)
    plt.savefig(fname, bbox_inches='tight')
    plt.close()

def plot_hist_residuals_wavesol(fname, data, indexes):
    """
    Plots a histogram of the data and a Gaussian fit to the data
    :param data: 2d darry of floats
    :param indexes: indexes in the data array for order, pixel, wavelength, residuals
    """
    splits = 3          # split pixel and orders in 3 parts, e.g. 9 different areas in total
    bins = 20           # How many bins
    [oo, xx, ww, dd] = indexes
    min_orders, max_orders = np.min(data[:,oo]), np.max(data[:,oo])+1E-6
    min_px, max_px = np.min(data[:,xx]), np.max(data[:,xx])+1E-6
    if data.shape[0] < 9*50:
        splits = 1
    range_orders = (max_orders - min_orders + 0.0)/splits
    range_px = (max_px - min_px + 0.0)/splits
    
    ran_hist = ( np.min(data[:,dd]), np.max(data[:,dd]) )
    data_x, data_y, label = [], [], []
    for ii in range(splits):            # for orders
        for jj in range(splits):        # for pixel
            mino, maxo = round(ii*range_orders+min_orders,0), round((ii+1)*range_orders+min_orders,0)
            minx, maxx = jj*range_px+min_px, (jj+1)*range_px+min_px
            good_data = ( (data[:,oo] >= mino) & (data[:,oo] < maxo) & (data[:,xx] >= minx) & (data[:,xx] < maxx) )
            label.append('Orders {0}-{1}\nPixel {2}-{3}'.format(int(mino), int(maxo), int(round(minx,0)), int(round(maxx,0))))
            hist_y, hist_x = np.histogram( data[good_data,dd], bins=bins, range=ran_hist )
            data_x.append(hist_x)
            data_y.append(hist_y)
            hist_x2 = (hist_x[1:]+hist_x[:-1])/2.           # center of the bin
            popt = centroid_order(hist_x2, hist_y, np.median(hist_x2), (hist_x2[-1]-hist_x2[0])/2., significance=3, bordersub_fine=False)    #a,x0,sigma,b in a*np.exp(-(x-x0)**2/(2*sigma**2))+b
            hist_x2_fine = np.linspace(np.min(hist_x), np.max(hist_x), 100)
            data_x.append(hist_x2_fine)
            data_y.append(oneD_gauss(hist_x2_fine,popt))
            label.append('centre {0}\nwidth {1}'.format(round(popt[1],4), round(popt[2],4)))
    text = 'Histogram of the residuals for {0} lines (width of 0.0 means a Gaussian fit was not possible)'.format(data.shape[0])
    plot_img_spec.plot_points(data_x, data_y, label, [fname], show=False, adjust=[0.12,0.88,0.85,0.12, 1.0,1.01], 
                              return_frame=False, title=text, x_title='Difference [px]', y_title='Count per bin', linestyle="-", marker="")       

def bck_px_UI(params, im_orig, pfits, xlows, xhighs, widths, w_mult, userinput=True):    # not used anymore
    fig, frame = plt.subplots(1, 1)
    def plot(frame, im_orig, pfits, xlows, xhighs, widths, w_mult):
        frame.clear()
        try:
            w_mult = float(gui3.data['width_multiplier'])
        except:
            w_mult = w_mult
        title = ('Determine area for the background, width_multiplier = {0})'.format(w_mult))
        im_bck_px = bck_px(im_orig, pfits, xlows, xhighs, widths, w_mult)
        plot_img_spec.plot_image(im_orig*im_bck_px, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], title, return_frame=True, frame=frame, colorbar=False)
        #frame.set_title(title)
    # get kwargs
    pkwargs = dict(frame=frame, im_orig = im_orig, pfits = pfits, xlows = xlows, xhighs = xhighs, widths = widths, w_mult = w_mult)
    # run initial update plot function
    plot(**pkwargs)
    
    # define valid_function
    # input is one variable (the string input)
    # return is either:
    #   True and values
    # or
    #   False and error message
    def vfunc(xs):
        try:
            w_mult = float(xs)
            return True, w_mult
        except:
            return False, ('Error, input must consist of a float')
    # define widgets
    widgets = dict()
    widgets['width_multiplier'] = dict(label='Multiplier to the measured width',
                                comment='The area of the measured width\nof the trace times this value\non either side of the trace\nwill be excluded from the\nbackground map' ,
                                kind='TextEntry', minval=0, maxval=None,
                                fmt=str, start=str(w_mult), valid_function=vfunc,
                                width=10)

    widgets['accept'] = dict(label='Accept Map', kind='ExitButton', position=Tk.BOTTOM)
    widgets['update'] = dict(label='Update', kind='UpdatePlot', position=Tk.BOTTOM)

    wprops = dict(orientation='v', position=Tk.RIGHT)

    gui3 = tkc.TkCanvas(figure=fig, ax=frame, func=plot, kwargs=pkwargs,
                        title='HiFlEx: Determine area for the background', widgets=widgets,
                        widgetprops=wprops)
    
    gui3.master.mainloop()
    
    params['width_multiplier'] = float(gui3.data['width_multiplier'])

    im_bck_px = bck_px(im_orig, pfits, xlows, xhighs, widths, params['width_multiplier'])
    #gui3.destroy 
    plt.close()
    return im_bck_px, params

def create_new_wavelength_UI( params, cal_l_spec, cal_s_spec, arc_lines_px, reference_lines_dict, adjust=[0.12,0.88,0.85,0.12, 1.0,1.01] ):   
    """
    :param arc_lines_px: numpy arry of floats: order, pixel, width, height of the line
    """
    reference_catalog, reference_names = reference_lines_dict['reference_catalog'][0], reference_lines_dict['reference_names'][0]
    px_to_wave_txt = read_text_file(params['px_to_wavelength_file'], no_empty_lines=True, warn_missing_file=False)              # list of strings
    if len(px_to_wave_txt) == 0:
        px_to_wave_txt = read_text_file(params['logging_path']+'tmp_'+params['px_to_wavelength_file'], no_empty_lines=True, warn_missing_file=False)              # list of strings
    px_to_wave = np.array( convert_readfile(px_to_wave_txt, [int, int, float, float], delimiter='\t', replaces=['\n',''], shorten_input=True, replacewithnan=True ))     # order, real order, px, wave

    nr_entries = len(px_to_wave)
    if nr_entries > 0:
        order_offset = np.nan
        if np.prod( np.isnan(px_to_wave[:,1]) ) == 0:
            order_offset = int(np.nanmedian(px_to_wave[:,1] - px_to_wave[:,0]))
        px_to_wave = np.hstack(( px_to_wave, np.zeros((nr_entries,2))*np.nan, np.expand_dims(np.arange(nr_entries), axis=1), np.zeros((nr_entries,1))*np.nan  ))     # order, real order, px, wave, width, height of line, index, nan
        order = -1E6
        for entry in arc_lines_px:
            if entry[0] != order:       # So this step is only done when neccessary 
                order = entry[0]
                px_to_wave_sub = px_to_wave[ px_to_wave[:,0] == order, : ]
            diff = np.abs(px_to_wave_sub[:,2] - entry[1])
            posi = np.where(diff < 2)[0]
            if len(posi) >= 1:
                index = int(px_to_wave_sub[posi[0], 6])      # first entry as highest signal, not min(diff)
                px_to_wave[index,4] = entry[2]
                px_to_wave[index,5] = entry[3]
            else:
                px_to_wave = np.vstack(( px_to_wave, [order, order+order_offset, entry[1], np.nan, entry[2], entry[3], len(px_to_wave), np.nan ] ))
    else:
        order_offset = np.nan
        tmp = arc_lines_px[:,0]*np.nan
        px_to_wave = np.vstack(( arc_lines_px[:,0], tmp, arc_lines_px[:,1], tmp, arc_lines_px[:,2], arc_lines_px[:,3], tmp, tmp )).T
    
    px_to_wave[np.isnan(px_to_wave[:,4]),4] = -1
    px_to_wave[np.isnan(px_to_wave[:,5]),5] = -1
    px_to_wave[:,6] = np.nan
    order = int(len(cal_l_spec)/2)
    
    fig, ax = plt.subplots(3, 1, gridspec_kw={'height_ratios': [2, 5, 2]})
    gui3 = tkc.TkCanvasGrid(figure=None, ax=None, func=None, title='HiFLEx: Test', kwargs=dict(), widgets=dict(), widgetprops=dict() )      # Only to get the display size
    dpi=100
    fig.set_dpi(dpi)
    fig.set_size_inches( (int(gui3.screen_w_h[0]) - 400)/dpi, (int(gui3.screen_w_h[1])-90)/dpi, forward=True)     # Make the plot as big as possible
    tkc.TkCanvasGrid.end(gui3)
    plt.subplots_adjust(left=adjust[0], right=adjust[1], top=adjust[2], bottom=adjust[3])
    pkwargs = dict()
    
    def mark_line_in_GUI(px):
        gui3.repopulate()
        px_to_wave = pkwargs['px_to_wave']
        px_to_wave_sub = px_to_wave[ px_to_wave[:,0] == pkwargs['order'], : ]
        match = np.where( np.abs(px_to_wave_sub[:,2] - px) < 20)[0]
        if len(match) == 0:
            unmark_line_in_GUI()
        else:
            index = match[0]
            for widget in gui3.ws.keys():
                if widget.find('wave_') > -1:
                    try:
                        number = int(widget.replace('wave_',''))
                    except:
                        continue
                    if number != index:
                        gui3.ws[widget].config(state=Tk.DISABLED)
        
    def unmark_line_in_GUI():
        gui3.repopulate()
    
    def onclick(event):
        #print('clicked',event.xdata, event.ydata, event.inaxes)
        if event.xdata is None or event.ydata is None:
            unmark_line_in_GUI()
        else:
            px = event.xdata
            mark_line_in_GUI(px)

    def plot(ax, cal_l_spec, cal_s_spec, px_to_wave, order, **pkwarks):
        max_reflines = 30
        for ii in range(3):
            ax[ii].clear()
            ordii = order+ii-1
            if ordii < 0 or ordii >= len(cal_l_spec):         # outside of useful range
                continue
            y1 = adjust_data_log( cal_l_spec[ordii,:] )
            y2 = adjust_data_log( cal_s_spec[ordii,:] )
            x = list(range(len(y1)))
            title   = None
            labels  = []
            x_title = ''
            y_title = 'log ( Flux [ADU] )'
            if ii == 0:
                title = 'Plot of the spectrum for previous ({0}), actual ({1}), and next order ({2})'.format(order-1, order, order+1)
                labels = ['long\nexp', 'short\nexp']
            if ii == 2:
                x_title = 'Dispersion direction [px]'
            ax[ii] = plot_img_spec.plot_points([x,x], [y1,y2], labels, [], show=False, adjust=[0.05,0.99,0.97,0.05, 0.92,1.2], title='', return_frame=True, frame=ax[ii], x_title=x_title, y_title=y_title, linestyle="-", marker="")
            
            pctl = 2./len(y1)*100                   # Exclude the two highest and two lowest points
            xmin, xmax = np.nanmin(x), np.nanmax(x)
            ymin1, ymax1 = np.nanpercentile(y1, pctl), np.nanpercentile(y1, 100-pctl)
            ymin2, ymax2 = np.nanpercentile(y2, pctl), np.nanpercentile(y2, 100-pctl)
            ymin, ymax = min(ymin1, ymin2), max(ymax1, ymax2)
            y_range = ymax - ymin
            x_range = xmax - xmin
            ax[ii].set_ylim([ymin-y_range*0.01, ymax+y_range*0.01])
            ax[ii].set_xlim([xmin-x_range*0.02, xmax+x_range*0.02])
            title = {0:'Previous', 1:'', 2:'Next'}[ii]
            ax[ii].set_title('{0} Order: {1}'.format(title, ordii), fontsize=11)
            
            px_to_wave_sub = px_to_wave[ px_to_wave[:,0] == ordii, : ]
            if ii == 1:
                goodentries = 0
            else:
                goodentries = 1E6           # to only plot the identified lines
            for jj, entry in enumerate(px_to_wave_sub):
                if entry[3] is np.nan or str(entry[3]) == 'nan':
                    if goodentries > max_reflines:
                        continue
                    text = '{0}'.format(round(entry[2],1))
                    color = 'g'
                    goodentries += 1
                else:
                    if ii == 1:
                        text = '{0} - {1}'.format(round(entry[2],1), round(entry[3],4))     # px and wave
                    else:
                        text = '{0}'.format(round(entry[3],4))                              # wave
                    color = 'r'
                ax[ii].plot( [entry[2],entry[2]], [ymin,ymin+y_range*0.01], color=color )
                ax[ii].text( entry[2], ymin+y_range*0.015, text, fontsize=8, 
                                    horizontalalignment='center', verticalalignment='bottom', rotation=90, color=color, zorder=5 )

            x_f = np.arange(xmin-x_range*0.02, xmax+x_range*0.02, 0.01)
            w_fa, w_fl = [], []
            if 'poly_lin_{0}'.format(ordii) in pkwargs.keys() and ii == 1:
                    w_fl = np.polyval(pkwargs['poly_lin_{0}'.format(ordii)], x_f)
            if 'wavelength_solution' in pkwargs.keys():       # when wavelength solution is not available
                    wavelength_solution = pkwargs['wavelength_solution']
                    w_fa = np.polyval(wavelength_solution[ordii,2:], x_f-wavelength_solution[ordii,1])
            if len(reference_catalog) == 0:
                    continue
            pos_index = 0
            ytop = ymax+y_range*0.01
            if np.sum(~np.isnan(px_to_wave_sub[:,3])) == 0 and px_to_wave_sub.shape[0] > 0:
                px_to_wave_sub[0,3] = -10000                    # Dummy wavelength
            for kk, w_f in enumerate([w_fl, w_fa]):
                if len(w_f) < 10:
                    continue
                inorder = (reference_catalog[:,0] >= min(w_f)) & (reference_catalog[:,0] <= max(w_f))
                reference_catalog_sub = np.array(sorted(reference_catalog[inorder,:], key=operator.itemgetter(1), reverse=True ))
                if len(reference_catalog_sub) == 0:
                    continue
                y_scale = y_range / (max(reference_catalog_sub[:,1])+0.0)
                pos_index += 1
                num_notident = 1
                for color in {0:['b', 'r'], 1:['g', 'r']}[kk]:        # plot the green/blue lines before the red ones
                    for jj,refline in enumerate(reference_catalog_sub):
                        if np.nanmin(np.abs(refline[0] - px_to_wave_sub[:,3])) < 1E-2:
                            if color != 'r':
                                continue        # don't plot a matching line when it's not red lines to be plotted
                        else:
                            if color == 'r':
                                continue        # don't plot the non matching lines when red lines to be plotted
                            if num_notident > max_reflines:
                                break           # stop plot non-matching lines
                            num_notident += 1   
                        index = np.argmin(np.abs(refline[0] - w_f))
                        x_pos = x_f[index]
                        #if np.isnan(spec1[order, x_position]):
                        #    continue
                        #y_position = np.nanmax(spec1[order, max(0,x_position-2):min(len(spec1[order,:]),x_position+2)])
                        #??? = max(y_pos, y_position+0.23*y_range)
                        text = '{0} - {1}'.format( round(refline[0],2+2*kk), reference_names[int(refline[2])] )
                        if kk == 0:                                         # Text below marker
                            y1 = [ymax, ymax+y_range*0.005+y_scale*0.03*refline[1]]
                            y2 = ymax
                            align = 'top'
                            ytop = max(y1[1], ytop)
                        else:                                               # Text above marker
                            y_pos = ymax - y_range*0.47                        # where is the second line
                            y1 = [y_pos-y_scale*0.03*refline[1]-y_range*0.005, y_pos]
                            y2 = y_pos+y_range*0.01
                            align = 'bottom'
                        ax[ii].plot( [x_pos,x_pos], y1, color=color )
                        ax[ii].text( x_pos, y2, text, fontsize=8, horizontalalignment='center', verticalalignment=align, rotation=90, color=color, zorder=5 )
            ax[ii].set_ylim([ymin-y_range*0.01, ytop])
    
    def calculate_linear_solution(pkwargs, order=None):
        px_to_wave = pkwargs['px_to_wave']
        if order is None:
            orders = np.unique(px_to_wave[:,0]).astype(int)
        else:
            orders = [order]
        for order in orders:
            inorder = px_to_wave[:,0] == order
            px_to_wave_sub = px_to_wave[ inorder, : ]
            good_values = ~np.isnan(px_to_wave_sub[:,3])
            if np.sum(good_values) <= 1:
                px_to_wave[inorder,7] = np.nan
                if 'poly_lin_{0}'.format(order) in pkwargs.keys():
                    del pkwargs['poly_lin_{0}'.format(order)]
                continue
            poly = np.polyfit( px_to_wave_sub[good_values,2], px_to_wave_sub[good_values,3], 1)     # 1 for linear fit
            px_to_wave[inorder,7] = np.polyval(poly, px_to_wave[inorder,2])
            pkwargs['poly_lin_{0}'.format(order)] = poly
        pkwargs['px_to_wave'] = px_to_wave
        return pkwargs
    
    def calculate_wavesolution_calc(px_to_wave, cal_l_spec):
        px_to_wave_sub = px_to_wave[ ~np.isnan(px_to_wave[:,3]), : ]

        wavelength_solution, wavelength_solution_arclines = fit_basic_wavelength_solution(params, px_to_wave_sub, cal_l_spec, 'GUI')
        
        for order in range(len(wavelength_solution)):
            values_order = px_to_wave[:,0] == order
            px_to_wave_sub = px_to_wave[ values_order, : ]
            px_to_wave[ values_order, 6 ] = np.polyval(wavelength_solution[order,2:], px_to_wave[ values_order, 2 ]-wavelength_solution[order,1])
        return px_to_wave, wavelength_solution, wavelength_solution_arclines
    
    def calculate_wavesolution():
        update_order()                      # To make sure the data is read again
        px_to_wave = pkwargs['px_to_wave']
        if np.sum(np.isnan(px_to_wave[:,1])) > 0:
            gui3.prompt('The order offset is not yet defined.\n\nPlease insert the correct number.')
            return
        pkwargs['wavelength_solution_info'] = True
        px_to_wave, wavelength_solution, wavelength_solution_arclines = calculate_wavesolution_calc(px_to_wave, cal_l_spec)
        
        pkwargs['px_to_wave']                   = px_to_wave
        pkwargs['wavelength_solution']          = wavelength_solution
        pkwargs['wavelength_solution_arclines'] = wavelength_solution_arclines
        # Replot
        """gui3.funkwargs['notclose'] = True
        tkc.TkCanvasGrid.end(gui3)          # will update also the plot"""
        if pkwargs['wavelength_solution_info']:
            gui3.prompt_info(calimages['wavelength_solution_result_text'], width=500)
            pkwargs['wavelength_solution_info'] = False
        update_order()                      # To update the wavelengths
        
    def remove_widgets():
        variables = [gui3.validation, gui3.fmts, gui3.entries, gui3.funcs, gui3.onclickvalue, gui3.data]
        for entry in pkwargs['widgets_change'].keys():
            if entry in pkwargs.keys():
                del pkwargs[entry]
            gui3.ws[entry].destroy()
            for variable in variables:
                if entry in variable.keys():
                    del variable[entry]
    
    def add_widgets(pkwargs):
        order = pkwargs['order']
        px_to_wave = pkwargs['px_to_wave']
        px_to_wave_sub = px_to_wave[ px_to_wave[:,0] == order, : ]
        wid_sub = dict()
        offset = 4
        ii = 0              # if the order doesn't contain any entries
        for ii, entry in enumerate(px_to_wave_sub):
            pkwargs['delete_{0}'.format(ii)]   = False
            wid_sub['delete_{0}'.format(ii)]   = dict(label=None,  kind='CheckBox', start=pkwargs['delete_{0}'.format(ii)], row=ii+offset, column=0)
            wid_sub['px_pos_{0}'.format(ii)] = dict(label='%6.1f'%entry[2], kind='Label', row=ii+offset, column=1, rowspan=1, columnspan=1, orientation=Tk.W)
            wave = {True:'', False:str(round(entry[3],6))}[entry[3] is np.nan or str(entry[3]) == 'nan']
            pkwargs['wave_{0}'.format(ii)] = wave
            wid_sub['wave_{0}'.format(ii)] = dict(kind='TextEntry', fmt=str, start=pkwargs['wave_{0}'.format(ii)], width=10, row=ii+offset, column=2, columnspan=2)
            text = '{0}  {1}'.format('%4.1f'%entry[4], '%5.1i'%entry[5])
            wid_sub['lineprop_{0}'.format(ii)] = dict(label=text, kind='Label', row=ii+offset, column=4, rowspan=1, columnspan=2, orientation=Tk.W)
            text = ''
            if np.isnan(entry[6]) == False and entry[6] > 1000 and entry[6] < 1E6:
                # Add button to copy the value
                text += '%8.3f'%entry[6]
            text += ' - '
            if np.isnan(entry[7]) == False and entry[7] > 1000 and entry[7] < 1E6:
                text += '%8.3f'%entry[7]
            wid_sub['wave_sol_{0}'.format(ii)] = dict(label=text, kind='Label', row=ii+offset, column=6, rowspan=1, columnspan=2, orientation=Tk.W)
        for jj in range(3):
            pkwargs['px_pos_extra_{0}'.format(jj)] = ''
            wid_sub['px_pos_extra_{0}'.format(jj)] = dict(kind='TextEntry', fmt=str, start=pkwargs['px_pos_extra_{0}'.format(jj)], width=6, row=ii+1+jj+offset, column=1, columnspan=1)
            pkwargs['wave_extra_{0}'.format(jj)] = ''
            wid_sub['wave_extra_{0}'.format(jj)] = dict(kind='TextEntry', fmt=str, start=pkwargs['wave_extra_{0}'.format(jj)], width=10, row=ii+1+jj+offset, column=2, columnspan=2)
        
        pkwargs['widgets_change'] = wid_sub
        return pkwargs, wid_sub
    
    def read_data(gui3):
        oldorder = gui3.funkwargs['oldorder']
        px_to_wave = gui3.funkwargs['px_to_wave']
        params['tmp_polynom_order_traces'] = [gui3.funkwargs['disporder']]
        params['tmp_polynom_order_intertraces'] = [gui3.funkwargs['crossdisporder']]
        keep = np.ones(( len(px_to_wave) )).astype(bool)
        # Read the results
        indexes = np.where(px_to_wave[:,0] == oldorder)[0]
        for ii,index in enumerate(indexes):
            #if 'wave_{0}'.format(ii) not in gui3.data.keys():       # can be empty when adding the extra stuff
            #    continue
            if gui3.data['delete_{0}'.format(ii)]:                  # delete
                keep[index] = False
            data = gui3.data['wave_{0}'.format(ii)]
            good_data = False
            if len(data) > 0:
                try:
                    data = float(data)
                    good_data = True
                except:
                    print('Warn: Can not convert entry into a number: {0}'.format(data))
            if good_data:
                px_to_wave[index,3] = data
            else:
                px_to_wave[index,3] = np.nan
        px_to_wave = px_to_wave[keep,:]
        for jj in range(3):
            good_data = False
            px = gui3.data['px_pos_extra_{0}'.format(jj)]
            wave = gui3.data['wave_extra_{0}'.format(jj)]
            if len(px) > 0 and len(wave) > 0:
                try:
                    px = float(px)
                    good_data = True
                except:
                    print('Warn: Can not convert entry into a number: {0}'.format(data))
                try:
                    wave = float(wave)
                except:
                    good_data = False
                    print('Warn: Can not convert entry into a number: {0}'.format(data))
            if good_data:
                px_to_wave = np.vstack(( px_to_wave, [oldorder, np.nan, px, wave, -1, -1, np.nan, np.nan ] ))
        if type(gui3.funkwargs['order_offset']).__name__ == 'int':
            px_to_wave[:,1] = px_to_wave[:,0] + gui3.funkwargs['order_offset']
        save_px_to_wave(px_to_wave, params['logging_path']+'tmp_'+params['px_to_wavelength_file'])
        gui3.funkwargs['px_to_wave'] = px_to_wave
        gui3.funkwargs = calculate_linear_solution(gui3.funkwargs, order=oldorder)
        
        #print 'new', gui3.data['order'], 'old', gui3.funkwargs['oldorder'], gui3.funkwargs['px_to_wave'][:,3]
    
    def update_order():
        tkc.TkCanvasGrid.update(gui3)               # Update order, get the new wavelengths
        read_data(gui3)              # read the new wavelengths, necessary here, as switched off at closing the gui
        remove_widgets()
        gui3.funkwargs['oldorder'] = copy.copy(gui3.data['order'])
        gui3.funkwargs, wid_sub = add_widgets(gui3.funkwargs)
        gui3.widgets = wid_sub
        gui3.add_widgets(parent=gui3.scrollable_frame)
        tkc.TkCanvasGrid.update(gui3)               # Update to get the plot with the new fitted wavelengths
    
    def save_px_to_wave(px_to_wave, fname):
        with open(fname,'w') as file:
            for entry in px_to_wave:
                entry = list(entry)
                for i in [1,3]:
                    if entry[i] is np.nan or str(entry[i]) == 'nan':
                        entry[i] = ''
                    elif i == 1:
                        entry[i] = int(entry[i])
                file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t\n'.format( int(entry[0]), entry[1], round(entry[2],2), entry[3], entry[4], entry[5] ))
    # define valid_function
    # input is one variable (the string input)
    # return is either:
    #   True and values
    # or
    #   False and error message
    def vfunc_int(xs):
        try:
            value = int(xs)
            return True, value
        except:
            return False, ('Error, input must be integer\n')
    def vfunc_float(xs):
        try:
            value = float(xs)
            return True, value
        except:
            return False, ('Error, input must be float\n')
    
    # get kwargs
    pkwargs = dict(ax=ax, cal_l_spec=cal_l_spec, cal_s_spec=cal_s_spec, px_to_wave=px_to_wave, order=order, order_offset=order_offset)
    pkwargs['wavelength_solution_info'] = False
    if np.sum(np.isnan(px_to_wave[:,1])) == 0:
        px_to_wave, wavelength_solution, wavelength_solution_arclines = calculate_wavesolution_calc(px_to_wave, cal_l_spec)
        pkwargs['wavelength_solution'] = wavelength_solution
        pkwargs['px_to_wave']          = px_to_wave
        pkwargs = calculate_linear_solution(pkwargs)
    # run initial update plot function
    plot(**pkwargs)
    # define widgets
    widgets = dict()
    widgets['ordertext']  = dict(label='Order:', kind='Label', row=0, column=0, rowspan=1, columnspan=4, orientation=Tk.W)
    pkwargs['order'] = order
    pkwargs['oldorder'] = order
    widgets['order'] = dict(kind='TextEntry', fmt=int, start=pkwargs['order'], valid_function=vfunc_int, minval=-1E-9, maxval=len(cal_l_spec)-1+1E-9, width=5, row=0, column=4, columnspan=1)
    widgets['update'] = dict(label='Update', kind='CommandButton', command=update_order, row=0, column=5, columnspan=2, width=9)
    widgets['accept'] = dict(label='Accept', kind='ExitButton', row=0, column=7, width=4)     # Replace by new fuction, as issue will occur if use changed the order but didn't updata
    widgets['order_offsettext']  = dict(label='Order offset between arbitary\nnumbering (starting at 0)\nand real physical order:', kind='Label', row=1, column=0, rowspan=1, columnspan=4, orientation=Tk.W)
    pkwargs['order_offset'] = {True:pkwargs['order_offset'], False:''}[pkwargs['order_offset'] is not np.nan]
    widgets['order_offset'] = dict(kind='TextEntry', fmt=int, start=pkwargs['order_offset'], valid_function=vfunc_int, minval=-1000, maxval=1000, width=5, row=1, column=4, columnspan=1)
    widgets['calculate_wavesolution'] = dict(label='Calculate\nWavelength\nsolution', kind='CommandButton', command=calculate_wavesolution, row=1, column=5, columnspan=2, width=9)
    widgets['disptxt']   = dict(label='Polynom-order for\ndispersion direction:', kind='Label', row=2, column=0, rowspan=1, columnspan=3, orientation=Tk.W)
    pkwargs['disporder'] = max(params.get('tmp_polynom_order_traces' ,params['polynom_order_traces']))
    widgets['disporder'] = dict(kind='TextEntry', fmt=int, start=pkwargs['disporder'], valid_function=vfunc_int, minval=1-1E-9, maxval=100, width=4, row=2, column=3, columnspan=1)
    widgets['crossdisptxt']   = dict(label='Polynom-order for cross-\ndispersion direction:', kind='Label', row=2, column=4, rowspan=1, columnspan=3, orientation=Tk.W)
    pkwargs['crossdisporder'] = max(params.get('tmp_polynom_order_intertraces' ,params['polynom_order_intertraces']))
    widgets['crossdisporder'] = dict(kind='TextEntry', fmt=int, start=pkwargs['crossdisporder'], valid_function=vfunc_int, minval=1-1E-9, maxval=100, width=4, row=2, column=7, columnspan=1)
    widgets['txtdelete']   = dict(label='Del-\nete', kind='Label', row=3, column=0, rowspan=1, columnspan=1, orientation=Tk.W)
    widgets['txtpx_pos']   = dict(label='Pixel', kind='Label', row=3, column=1, rowspan=1, columnspan=1, orientation=Tk.W)
    widgets['txtwave']     = dict(label='Wavelength', kind='Label', row=3, column=2, rowspan=1, columnspan=2, orientation=Tk.W)
    widgets['txtlineprop'] = dict(label='Line width\n+ height', kind='Label', row=3, column=4, rowspan=1, columnspan=2, orientation=Tk.W)   # Can't be float, as otherwise deleting an entry wouldn't work
    widgets['txtwave_sol'] = dict(label='Wavelenght-solution\nglobal - linear', kind='Label', row=3, column=6, rowspan=1, columnspan=2, orientation=Tk.W)   # Can't be float, as otherwise deleting an entry wouldn't work
    pkwargs, wid_sub = add_widgets(pkwargs)
    widgets.update(wid_sub)
    
    wprops = dict(fullscreen=False )
    gui3 = tkc.TkCanvasGrid(figure=fig, ax=ax, func=plot, title='HiFLEx: Create a new wavelength solution', 
                            kwargs=pkwargs, widgets=widgets, widgetprops=wprops )
    cid = fig.canvas.callbacks.connect('button_press_event', onclick)       # It works with GUI
    fig.set_size_inches( (int(gui3.screen_w_h[0]) - gui3.width_GUI-50)/dpi, (int(gui3.screen_w_h[1])-90)/dpi, forward=True)     # Make the plot as big as possible
        
    gui3.master.mainloop()
        
    read_data(gui3)             # read the new wavelengths only when closing the window
    plt.close()
    
    px_to_wave = gui3.funkwargs['px_to_wave']
    if os.path.isfile(params['px_to_wavelength_file']):
        os.system('mv {0} old_{0}'.format(params['px_to_wavelength_file']))
    save_px_to_wave(px_to_wave, params['px_to_wavelength_file'])
    del params['tmp_polynom_order_traces']
    del params['tmp_polynom_order_intertraces']
    
    if np.sum(np.isnan(px_to_wave[:,1])) == 0:
        px_to_wave, wavelength_solution, wavelength_solution_arclines = calculate_wavesolution_calc(px_to_wave, cal_l_spec)
    else:
        logger('Error: Not able to find a wavelength solution with the information provided by the user. Please rerun the program.')
    
    return dict(wavesol=wavelength_solution, reflines=wavelength_solution_arclines)

def correlate_px_wave_result_UI(im, arc_lines_wavelength, reference_catalog, arc_lines_px, reference_names, wavelength_solution, wavelength_solution_arclines, adjust=[0.12,0.88,0.85,0.12, 1.0,1.01]):
    ims = im.shape
            
    fig, frame = plt.subplots(1, 1)
    plt.subplots_adjust(left=adjust[0], right=adjust[1], top=adjust[2], bottom=adjust[3])
    
    order = 10
    high_limit = 98
    arc_stretch = 0.5
    
    def plot(frame, im, order, high_limit, arc_stretch):
        x_range, y_range = copy.copy(frame.get_xlim()), copy.copy(frame.get_ylim())
        frame.clear()
        """try:
            order = gui3.data['order']
        except:
            order = order
        try:
            high_limit = gui3.data['high_limit']
        except:
            high_limit = high_limit
        try:
            arc_stretch = gui3.data['arc_stretch']
        except:
            arc_stretch = arc_stretch"""
        
        #print 'order, arc_setting', order, arc_setting
        title = ('Order {0}: identified emission lines: blue [px],\n'+\
                 'catalog lines: red (0.1px precission, name and wavelength at the botom),\n'+\
                 'corellated lines: green (px, wavelength at top)').format(order)
        
        # Plot the extracted arc spectrum
        xarr = np.arange(len(im[order,:]))
        yarr = im[order,:]
        notnans = ~np.isnan(yarr)
        yarr = yarr[notnans]
        xarr = xarr[notnans]
        yarr_max = max(0,np.percentile(yarr,high_limit))
        yarr[yarr > yarr_max] = yarr_max
        frame = plot_img_spec.plot_points([xarr], [yarr], [], [], show=False, adjust=adjust, title=title, 
                                      return_frame=True, frame=frame, x_title='x [px]', y_title='flux [ADU]', linestyle="-", marker="")
        
        # Plot the correlated lines in the extracted spectrum
        arc_line_order = arc_lines_px[(arc_lines_px[:,0] == order),:]
        xarr = np.vstack(( arc_line_order[:, 1], arc_line_order[:, 1])).T
        yplot = np.repeat( [[min(yarr), max(yarr)]], len(arc_line_order), axis=0)
        frame = plot_img_spec.plot_points(xarr, yplot, [], [], show=False, adjust=adjust, title=title, 
                                      return_frame=True, frame=frame, x_title='x [px]', y_title='flux [ADU]', linestyle="-", marker="", color='b')

        # Plot the catalog lines using the fit
        xarr_catalog = np.arange(-100,ims[1]+100, 0.1)
        xarr_wave = np.polyval(wavelength_solution[order,2:], xarr_catalog-wavelength_solution[order,1])
        inorder = ( reference_catalog[:,0] >= min(xarr_wave) ) & ( reference_catalog[:,0] <= max(xarr_wave) )
        reference_catalog_sub = reference_catalog[inorder, :]
        reference_names_sub = np.array(reference_names)[inorder]
        arc_stretch *= (max(yarr)-min(yarr)) / max(reference_catalog_sub[:,1])      # Rescale the arc lines so that arc_stretch=1 fills it just once
        for line_index in range(reference_catalog_sub.shape[0]):
            diff = abs(reference_catalog_sub[line_index,0] - xarr_wave)
            xarc = xarr_catalog[np.argmin(diff)]
            frame.plot( [xarc,xarc],[min(yarr),min(yarr)+reference_catalog_sub[line_index,1]*arc_stretch], color='r' )
            frame.text(xarc, min(yarr)+reference_catalog_sub[line_index,1]*arc_stretch, '{0} {1}'.format(reference_names_sub[line_index],reference_catalog_sub[line_index,0]),
                            horizontalalignment='center', verticalalignment='bottom', rotation=90, color='k', zorder=5)
           
        # Plot the identified catalog lines
        arc_line_order = arc_lines_wavelength[(arc_lines_wavelength[:,0] == order),:]
        xarr = np.vstack(( arc_line_order[:, 1], arc_line_order[:, 1])).T
        yplot = np.repeat( [[min(yarr), max(yarr)]], len(arc_line_order), axis=0)
        frame = plot_img_spec.plot_points(xarr, yplot, [], [], show=False, adjust=adjust, title=title, 
                                      return_frame=True, frame=frame, x_title='x [px]', y_title='flux [ADU]', linestyle="-", marker="", color='g')
        for arc_line in arc_line_order:
            frame.text(arc_line[1], max(yarr), r'{0} $\pm$ {1}'.format(arc_line[2], round(arc_line[3],3) ), horizontalalignment='center', verticalalignment='top', rotation=90, color='k', zorder=5)

    # get kwargs
    pkwargs = dict(frame=frame, im=im, order=order, high_limit=high_limit, arc_stretch=arc_stretch)
    # run initial update plot function
    plot(**pkwargs)
    
    # define valid_function
    # input is one variable (the string input)
    # return is either:
    #   True and values
    # or
    #   False and error message
    def vfunc_int(xs):
        try:
            value = int(xs)
            return True, value
        except:
            return False, ('Error, input must be integer')
    def vfunc_float(xs):
        try:
            value = float(xs)
            return True, value
        except:
            return False, ('Error, input must be float')
            
    # define widgets
    widgets = dict()
    widgets['order'] = dict(label='Order',
                                comment='which order to show' ,
                                kind='TextEntry', minval=0, maxval=ims[0]-1,
                                fmt=int, start=order, valid_function=vfunc_int,
                                width=10)
    widgets['high_limit'] = dict(label='High limit [%]',
                                comment='Reject highest pixels' ,
                                kind='TextEntry', minval=1, maxval=105,
                                fmt=float, start=high_limit, valid_function=vfunc_float,
                                width=10)
    widgets['arc_stretch'] = dict(label='Scale of\ncatalogue lines',
                                comment='float value' ,
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=float, start=arc_stretch, valid_function=vfunc_float,
                                width=10)
    widgets['accept'] = dict(label='Close', kind='ExitButton', position=Tk.BOTTOM)
    widgets['update'] = dict(label='Update', kind='UpdatePlot', position=Tk.BOTTOM)

    wprops = dict(orientation='v', position=Tk.RIGHT)

    gui3 = tkc.TkCanvas(figure=fig, ax=frame, func=plot, kwargs=pkwargs,
                        title='HiFlEx: Check the identified and correlated emission lines', widgets=widgets,
                        widgetprops=wprops)
    gui3.master.mainloop()
    plt.close()
    

def correlate_UI(im, order, arc_settings, reference_catalog, reference_names, adjust=[0.07,0.93,0.94,0.06, 1.0,1.01]):      # not used anymore
    global oldorder
    oldorder = order
    
    fig, frame = plt.subplots(1, 1)
    plt.subplots_adjust(left=adjust[0], right=adjust[1], top=adjust[2], bottom=adjust[3])
    
    def plot(frame, im):
        global oldorder
        x_range, y_range = copy.copy(frame.get_xlim()), copy.copy(frame.get_ylim())
        frame.clear()
        try:
            order = gui3.data['order']
        except:
            order = oldorder
        try:
            high_limit = gui3.data['high_limit']
        except:
            high_limit = 98
        if order < 0 or order > im.shape[0]-1:
            print('Warn: order outside the allowed area: 0...{0}'.format(im.shape[0]-1))
            order = oldorder
        if order != oldorder:
            arc_setting = arc_settings[order]
            gui3.ws[2].delete(0, 200)
            gui3.ws[2].insert(0,str(arc_setting[0]))
            gui3.ws[3].delete(0, 200)
            gui3.ws[3].insert(0,'{0}:{1}'.format(arc_setting[1], arc_setting[2]))
            gui3.ws[4].delete(0, 200)
            gui3.ws[4].insert(0,str(arc_setting[3]))
        else:
            try:
                arc_settings[order,0] = gui3.data['resolution']
                arc_settings[order,1],arc_settings[order,2] = gui3.data['zeropoint']
                arc_settings[order,3] = gui3.data['arc_stretch']
            except:
                False
            arc_setting = arc_settings[order]
        #print 'order, arc_setting', order, arc_setting
        title = 'Order {}'.format(order)
        #print zeropoint[0],resolution,zeropoint[1]
        xarr_arc = (reference_catalog[:,0]-arc_setting[1]+0.0)/arc_setting[0]+arc_setting[2]
        xarr = np.arange(len(im[order,:]))
        yarr = im[order,:]
        notnans = ~np.isnan(yarr)
        yarr = yarr[notnans]
        xarr = xarr[notnans]
        yarr_max = max(0,np.percentile(yarr,high_limit))
        yarr[yarr > yarr_max] = yarr_max
        frame.plot(xarr, yarr)
        for i, xarc in enumerate(xarr_arc):
            if xarc > xarr[0]-100 and xarc < xarr[-1]+100:
                frame.plot( [xarc,xarc],[min(yarr),min(yarr)+reference_catalog[i,1]*arc_setting[3]],color='r' )
                frame.text(xarc, min(yarr)+reference_catalog[i,1]*arc_setting[3], '{0} {1}'.format(reference_names[i],reference_catalog[i,0]),
                       horizontalalignment='center', verticalalignment='center',
                       rotation=90, color='k', zorder=5) 
        #plot_ranges = min()
        #if x_range != (0.0, 1.0) and y_range != (0.0, 1.0):
        #    plt.axis([x_range[0], x_range[1], y_range[0], y_range[1]])
        #else:
        #    dx = (plot_ranges[1] - plot_ranges[0])*0.01
        #    dy = (plot_ranges[3] - plot_ranges[2])*0.01
        #    plt.axis([plot_ranges[0]-dx,plot_ranges[1]+dx, plot_ranges[2]-dy, plot_ranges[3]+dy])
        frame.set_xlabel('x [px]', fontsize=14)
        frame.set_ylabel('flux [ADU]', fontsize=14)
        frame.set_title(title, fontsize=16)
        #frame.legend(loc='upper left', bbox_to_anchor=(adjust[4], adjust[5]))
        oldorder = order
        return arc_settings
    # get kwargs
    pkwargs = dict(frame=frame, im = im)
    # run initial update plot function
    arc_settings = plot(**pkwargs)
    
    # define valid_function
    # input is one variable (the string input)
    # return is either:
    #   True and values
    # or
    #   False and error message
    def vfunc_colon(xs):
        try:
            new_xs = xs.split(':')
            xs = []
            for entry in new_xs:
                xs.append(float(entry))
            return True, xs
        except:
            return False, ('Error, input must consist of floats \n '
                           'separated by colon')
    def vfunc_int(xs):
        try:
            value = int(xs)
            return True, value
        except:
            return False, ('Error, input must be integer')
    def vfunc_float(xs):
        try:
            value = float(xs)
            return True, value
        except:
            return False, ('Error, input must be float')
            
    # define widgets
    widgets = dict()
    widgets['high_limit'] = dict(label='High limit',
                                comment='Reject highest pixels' ,
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=str, start=str(98), valid_function=vfunc_float,
                                width=10)
    widgets['order'] = dict(label='Order',
                                comment='which order to show' ,
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=str, start=str(oldorder), valid_function=vfunc_int,
                                width=10)
    widgets['resolution'] = dict(label='Resolution',
                                comment='in A/px' ,
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=str, start=str(arc_settings[oldorder,0]), valid_function=vfunc_float,
                                width=10)
    widgets['zeropoint'] = dict(label='zeropoint',
                                comment='wavelength@px' ,
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=str, start='{0}:{1}'.format(arc_settings[oldorder,1], arc_settings[oldorder,2]), valid_function=vfunc_colon,
                                width=15)
    widgets['arc_stretch'] = dict(label='Scale of Arc lines',
                                comment='float value' ,
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=str, start=str(arc_settings[oldorder,3]), valid_function=vfunc_float,
                                width=10)
    widgets['accept'] = dict(label='Save and Close', kind='ExitButton', position=Tk.BOTTOM)
    widgets['update'] = dict(label='Update', kind='UpdatePlot', position=Tk.BOTTOM)

    wprops = dict(orientation='v', position=Tk.RIGHT)

    gui3 = tkc.TkCanvas(figure=fig, ax=frame, func=plot, kwargs=pkwargs,
                        title='HiFlEx: Plot spectral orders', widgets=widgets,
                        widgetprops=wprops)
    
    gui3.master.mainloop()
    plt.close()
    
    #print 'arc_settings',arc_settings
    return arc_settings

def create_pseudo_wavelength_solution(number_orders):
    """
    Creates a pseudo wavelength solution with one Angstrom per pixel
    :param number_orders: interger, number of orders for which the solution should be created
    
    """
    wavelength_solution, wavelength_solution_arclines = [], []
    for order in range(number_orders):
        wavelength_solution.append([order, 0., 1., 1. ])
        wavelength_solution_arclines.append([0])     #The wavelength of the reference lines
    return dict(wavesol=np.array(wavelength_solution), reflines=np.array(wavelength_solution_arclines) )

def fit_wavelengths_solution_2d(arc_lines_wavelength, cen_pxs, order_offset, polynom_order_trace, polynom_order_intertrace):
        x = arc_lines_wavelength[:,1]-cen_pxs
        y = 1.0/(arc_lines_wavelength[:,0]+order_offset)
        
        """weight = arc_lines_wavelength[:,4] * arc_lines_wavelength[:,5] * arc_lines_wavelength[:,6] * abs(arc_lines_wavelength[:,1])**(polynom_order_trace-1)
        median = []
        for order in np.unique(arc_lines_wavelength[:,0]):
            median.append( np.max( weight[arc_lines_wavelength[:,0] == order] ))
        divisor = np.median(median)
        if divisor == 0:
            divisor = 1.
        median = np.array(median)/(divisor+0.0)
        median[median == 0] = 1.
        for i, order in enumerate(np.unique(arc_lines_wavelength[:,0])):
            weight[arc_lines_wavelength[:,0] == order] /= (median[i]+0.0)"""
        weight = []         # No weight fits the data better
        poly2d_params = polynomial_fit_2d_norm(x, y, arc_lines_wavelength[:,2], polynom_order_trace, polynom_order_intertrace, w=weight)
        new_waves = polynomial_value_2d(x, y, polynom_order_trace, polynom_order_intertrace, poly2d_params)
        
        return poly2d_params, new_waves

def transform_wavelength_solution_2d_to_n1d(number_orders, number_pixel, polynom_order_trace, polynom_order_intertrace, poly2d_params, order_offset, cen_px, arc_lines_wavelength):
    """
    Transform the 2d wavelength solution into the old n-times 1d wavelength solution
    """
    max_number_reflines = 0             # will be needed later in order to save the data correctly into a fits file
    for order in range(number_orders):
        inorder = ( arc_lines_wavelength[:,0]==order )
        max_number_reflines = max(max_number_reflines, np.sum(inorder) )
    
    wavelength_solution = np.zeros((number_orders, 2+polynom_order_trace+1))
    wavelength_solution_arclines = np.zeros((number_orders, max_number_reflines))
    line_stats = np.zeros((number_orders, 2))
    xarr = np.arange(0, number_pixel, 0.1)
    wavelength_solution[:,0] = order_offset+np.arange(number_orders)
    wavelength_solution[:,1] = cen_px
    for order in range(number_orders):
        yarr = polynomial_value_2d(xarr-cen_px[order], 1.0/(order+order_offset), polynom_order_trace, polynom_order_intertrace, poly2d_params)
        polyfit = np.polyfit(xarr-cen_px[order], yarr, polynom_order_trace)      #lambda from px
        wavelength_solution[order,2:] = polyfit
        inorder = (arc_lines_wavelength[:,0]==order)
        if np.sum(inorder) > 0:
            line_width = arc_lines_wavelength[inorder,9]                        # Width of the line
            # Get the real wavelength at the line positions instead of the catalogue lines:
            refline_wavereal = polynomial_value_2d(arc_lines_wavelength[inorder,1]-cen_px[order], 1.0/(order+order_offset), polynom_order_trace, polynom_order_intertrace, poly2d_params)
            #wavelength_solution_arclines[order,:np.sum(inorder)] = arc_lines_wavelength[inorder,2]       # Wavelength in reference catalogue
            wavelength_solution_arclines[order,:np.sum(inorder)] = refline_wavereal                    # Added on 20200720
            line_stats[order,0] = np.median(line_width)
            if np.sum(inorder) > 1:
                line_stats[order,1] = np.std(line_width, ddof=1)
    
    return wavelength_solution, wavelength_solution_arclines, line_stats

def compare_wavelength_solution_to_emission_lines(kwargs):
    [ resdiff, pxdiff, pxdifford, orderdiff, wavelength_solution_2, arc_lines_px, ignoreorders, specs, reference_catalog, max_diffs ] = kwargs
    matching_lines = np.empty(( len(arc_lines_px)*10, 4 ))
    index_matching = 0
    data_available = [0, 0, 0]
    # Caluculate the wavelength of all arclines with the current settings
    for order_arcline in np.arange(min(arc_lines_px[:,0]), max(arc_lines_px[:,0])+1, dtype=int):        # The new orders
        if order_arcline+orderdiff < 0 or order_arcline+orderdiff > len(wavelength_solution_2)-1 or order_arcline in ignoreorders:    # ignore orders that are not covered by orig solution
            # No solution is available for this shifted order
            continue
        arc_lines_order_px = arc_lines_px[arc_lines_px[:,0] == order_arcline,1]     # array with px position of the identified lines for this order
        # Only use lines for which the wavelength solution is defined (assuming that the whole CCD  was used, otherwise needs old xlows and old xhighs)
        arc_lines_order_px = arc_lines_order_px[ (arc_lines_order_px - pxdiff - pxdifford * (order_arcline+orderdiff) > 0) & \
                                                 (arc_lines_order_px - pxdiff - pxdifford * (order_arcline+orderdiff) < specs[1]) ]      # order_arcline+orderdiff is correct
        # Get a normalisierungs parameter
        px_range = [max(0, 0 - pxdiff - pxdifford * (order_arcline+orderdiff)), min(specs[1], specs[1] - pxdiff - pxdifford * (order_arcline+orderdiff))]
        ## Using pixel for data available ignores the different number of lines available in different orders (100 lines @ order 70; 40 @ 0)
        data_available[0] += px_range[1] - px_range[0]    # number of pixel covered
        wave_range = [ np.polyval(wavelength_solution_2[order_arcline+orderdiff, 2:], px_range[0] - wavelength_solution_2[order_arcline+orderdiff, 1]),
                       np.polyval(wavelength_solution_2[order_arcline+orderdiff, 2:], px_range[1] - wavelength_solution_2[order_arcline+orderdiff, 1]) ]
        covered_reference = (reference_catalog[:,0] >= wave_range[0]) & (reference_catalog[:,0] <= wave_range[1])
        data_available[1] += np.sum(covered_reference)
        data_available[2] += np.sum(reference_catalog[covered_reference,1])
        # Caluculate the wavelength of all arclines with the current settings
        arc_lines_order_wl = np.polyval(wavelength_solution_2[order_arcline+orderdiff, 2:], arc_lines_order_px - wavelength_solution_2[order_arcline+orderdiff, 1])
        # Get the offset to the closest lines
        for i, arc_line_wave in enumerate(arc_lines_order_wl):
            diff_catalog = reference_catalog[:,0] - arc_line_wave
            good_values = (abs(diff_catalog) <= max_diffs[order_arcline+orderdiff])
            diff_catalog = diff_catalog[good_values]
            if index_matching+len(diff_catalog) > matching_lines.shape[0]:                          # increase the array size
                matching_lines = np.vstack([ matching_lines, np.empty( matching_lines.shape ) ])    # double the array length
            matching_lines[index_matching:index_matching+len(diff_catalog), :] = \
                            np.vstack([ np.repeat([order_arcline],len(diff_catalog)), np.repeat([arc_line_wave],len(diff_catalog)), 
                                        diff_catalog, reference_catalog[good_values,1] ]).T
            index_matching += len(diff_catalog)
            #if order_arcline == 36:
            #    print resdiff, pxdiff, pxdifford, orderdiff, order_arcline, arc_lines_order_px[i], arc_line_wave, reference_catalog[good_values,0], reference_catalog[good_values,1]
            #print order_arcline, arc_line_wave, diff_catalog
            #for entry in diff_catalog:
            #    matching_lines.append([order_arcline, arc_line_wave, entry])
    if index_matching == 0:
        return []
    #matching_lines = np.array(matching_lines)
    matching_lines = matching_lines[:index_matching,:]                                  # get rid of unfilled array
    #best_matching: orderdiff, pxdiff, pxdifford, resdiff, number of identified reference lines, 
    #               stdev of the wavelength difference between current solution and reference wavelength, sum of the flux of the reference lines, 
    #               number of covered pixels, number of covered reference lines, sum of covered flux
    return [orderdiff, pxdiff, pxdifford, resdiff, len(matching_lines), np.std(matching_lines[:,2], ddof=1), np.sum(matching_lines[:,3])] + data_available


def adjust_wavelength_solution(params, spectrum, arc_lines_px, wavelength_solution_ori, wavelength_solution_arclines_ori, reference_lines_dict, traces_def, show_res=False, search_order_offset=False):
    """
    :param wavelength_solution_ori: 2d array of floats, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    :param arc_lines_px: 2d array with one line for each identified line, sorted by order and amplitude of the line. For each line the following informaiton is given:
                    order, pixel, width of the line, and height of the line
    :param wavelength_solution_arclines_ori: 2d array of floats with one line for each order. Each order contains the wavelengths of the identified reference lines, 
                                        sorted by the brightest to faintest (in spectrum from which the solution is derived) and 0 to make into an array
                                        not used, only to return the correct values when no solution is found
    :param reference_catalog: 2d array of floats with one entry for each line. Each entry contains the wavelength and the intensity of the line
    :param reference_names: list of strings with same length as reference_catalog, name of each line
    :return wavelength_solution: 2d array of floats, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    :return wavelength_solution_arclines: 2d array of floats with one line for each order. Each order contains the wavelengths of the identified reference lines, 
                                        sorted by the brightest to faintest (in spectrum from which the solution is derived).
                                        To create the same number of lines for each order, the array is filled with 0
    """
    ignoreorders = [] #[0,1,2,3,4,5,6,7,8]
    iter_break = 2     # Normally the script stops already at step 4-8, as no further improvements are possible
    steps = 20           # Standard: 20
    sigma = 3.5         # Standard: 3.5, the smaller the better the solution, but lines at the corners of the images might be rejected
    only_one_line = True
        
    [tr_poly, xlows, xhighs, widths] = traces_def
    reference_catalog_full, reference_names = reference_lines_dict['reference_catalog'][0], reference_lines_dict['reference_names'][0]
    specs = spectrum.shape
    orders = np.arange(specs[0])
    if len(arc_lines_px) <= 10:
        logger('Warn: no arc lines available -> creating a pseudo solution (1 step per px)')
        wavelength_solution, wavelength_solution_arclines = create_pseudo_wavelength_solution(specs[0])
        return wavelength_solution, wavelength_solution_arclines
    #reference_catalog = reference_catalog_full[ reference_catalog_full[:,1] >= np.percentile(reference_catalog_full[:,1], 50) ,:]       # only use the brightest lines at the beginning
    reference_catalog = copy.deepcopy(reference_catalog_full)
    cen_resolution = []
    for wls in wavelength_solution_ori:
        cen_resolution.append( abs( np.polyval(wls[2:], int(specs[1]/2) + 1 - wls[1]) - np.polyval(wls[2:], specs[1]/2 - 1 - wls[1]) ) / 2. )    #resolution around the middle of the chip
        #print np.polyval(wls[2:], 0 - wls[1]), np.polyval(wls[2:], specs[1] - wls[1])
    max_diffs = np.array(cen_resolution) * params['px_offset'][2]  # Assign lines only to resolution * step size of px_offset
    orderdiffs = range(min(params['order_offset'][0:2]),max(params['order_offset'][0:2])+1)
    pxdiffs = list(range(min(params['px_offset'][0:2]),max(params['px_offset'][0:2])+1,params['px_offset'][2]))
    #if not ( min(params['px_offset'][0:2]) < 0 and min(params['px_offset'][0:2]) > 0 ):                                     # if it's from -x to +y the user wants only the above pxdiffs
    #    pxdiffs += range(-1*max(params['px_offset'][0:2]),-1*min(params['px_offset'][0:2])+1,params['px_offset'][2])        # if it's from +x to +y, then the user wants the range in +/-
    pxdiffs = list(set(pxdiffs))            # remove duplicates
    pxdiffords = np.arange(min(params['px_offset_order'][0:2]),max(params['px_offset_order'][0:2])+params['px_offset_order'][2],params['px_offset_order'][2])
    resolution_offset = params['resolution_offset_pct']/100.
    res_steps = 11
    if resolution_offset == 0:
        res_steps = 1
    resdiffs = np.linspace(1-resolution_offset, 1+resolution_offset, res_steps)
    best_matching = []
    #logger('Step: Compare old wavelength solution with current arc lines')
    # Rescaling the original wavelength solution, calculating the wavelength for the identied lines with this solution
    kwargs = []
    for resdiff in resdiffs:
        wavelength_solution_1 = copy.deepcopy(wavelength_solution_ori)
        wavelength_solution_1[:,-2] *= resdiff                                                      # Scale the resolution
        for pxdiff in pxdiffs:
            for pxdifford in pxdiffords:
                wavelength_solution_2 = copy.deepcopy(wavelength_solution_1)
                wavelength_solution_2[:,1] += pxdiff + pxdifford * np.arange(len(wavelength_solution_2))        # that is equivalent to changing the px position in arc_lines_order_px
                for orderdiff in orderdiffs:
                    kwargs.append([ resdiff, pxdiff, pxdifford, orderdiff, wavelength_solution_2, arc_lines_px, ignoreorders, specs, reference_catalog, max_diffs ])
                    
    if params['use_cores'] > 1 and multiprocessing.current_process().name == 'MainProcess':
       logger('Compare old wavelength solution with current emission lines (multicore)')
       p = multiprocessing.Pool(params['use_cores'])
       best_matching = p.map(compare_wavelength_solution_to_emission_lines, kwargs)
       p.terminate()
       p.join()
    else:
       for kwarg in tqdm(kwargs, desc='Compare old wavelength solution with current emission lines'):
            best_matching.append( compare_wavelength_solution_to_emission_lines(kwarg) )
    for ii in range(len(best_matching))[::-1]:
        if len(best_matching[ii]) == 0:
            del best_matching[ii]
    if len(best_matching) == 0:
        if len(wavelength_solution_ori) == specs[0]:
            logger('Warn: No matching configuration of the lines in the emission line spectrum with the old wavelength solution found. Therefore the old solution will be used')
            return np.array(wavelength_solution_ori), wavelength_solution_arclines_ori
        else:
            logger('Warn: No matching configuration of the lines in the emission line spectrum with the old wavelength solution found. Additionally the number of orders in the original solution and this setting do not match -> creating a pseudo solution (1 step per px)')
            wavelength_solution, wavelength_solution_arclines = create_pseudo_wavelength_solution(specs[0])
            return wavelength_solution, wavelength_solution_arclines
    best_matching = np.array(best_matching, dtype=float)
    for i in [7,8,9]:
        best_matching[:,i] /= (np.median(best_matching[:,i])+0.0)                       # normalise
    
    """for y,yt in [[4,'lines'], [5,'variation'], [6,'total flux']]:
        for x,xt in [[0,'orddiff'], [1,'pxdiff'], [2,'pxdifford'], [3,'resdiff']]:
            plot_img_spec.plot_points([best_matching[:,x],best_matching[:,x],best_matching[:,x]],
                                      [best_matching[:,y]/best_matching[:,8],best_matching[:,y]/best_matching[:,9],best_matching[:,y]],
                                      ['norm1','norm2','ori'],'path',show=True, marker=['s','*','<'], x_title=xt, y_title=yt)"""
    
    best_matching = best_matching[np.argsort(best_matching[:,4]/best_matching[:,8], axis=0),:]     #sort by number of identified arclines (scaled by available lines) -> best at the end
    #print best_matching[-10:,:]
    best_matching = best_matching[int(best_matching.shape[0]*8/9):,:]     #only 25% of the best values

    best_matching = best_matching[np.argsort(best_matching[:,6]/best_matching[:,9], axis=0),:]     # sort by the sum of the intensity of the arc lines (scaled by available flux) -> best at end
    #print best_matching[-10:,:]
    best_matching = best_matching[int(best_matching.shape[0]*8/9):,:]     #only 25% of the best values
    #best_matching = best_matching[np.argsort(best_matching[:,5], axis=0),:]     #sort by standard deviation of the shift of the arc lines (don't scale!) -> best at the front =>> !! that isn't a good selection at all
    #print best_matching[:10]
    
    # Use this solution to identify lines and fit polynomials in each order
    orderdiff, pxdiff, pxdifford, resdiff = int(np.median(best_matching[:,0])), int(np.median(best_matching[:,1])), np.median(best_matching[:,2]), np.median(best_matching[:,3])
    logger('Info: To match the most lines in the emission line spectrum with the old wavelength solution, a shift of {0} orders, a multiplier to the resolution of {1}, a shift of {2} px, and a shift of {3} px per order needs to be applied. {4} lines were identified. The deviation is {5} Angstrom.'.format(orderdiff, round_sig(resdiff,3), pxdiff, round(pxdifford,2), int(best_matching[0,4]), round(best_matching[0,5],4) ))
    ## assign lines in the arcline with lines from the reference catalog
    med_arc_width = np.nanmedian(arc_lines_px[:,2])             # median width of the arc lines
    arc_lines_wavelength = []       #order, central pixel of the arc line, wavelength of the assigned reference line, diff between both solutions, height of the arc line, intensity of the reference line, 1/width of the arc line
    wavelength_solution_2 = copy.copy(wavelength_solution_ori)
    #wavelength_solution_2[:,0]  -= orderdiff                                                   # don't apply this. The central wavelength stays the same, hence also the real order should stay the same
    wavelength_solution_2[:,-2] *= resdiff                                                      # Scale the resolution
    wavelength_solution_2[:,1] += pxdiff + pxdifford * np.arange(len(wavelength_solution_1))   # Shift the central pixel
    #orders_prev = np.arange(min(arc_lines_px[:,0]) ,max(arc_lines_px[:,0])+1, dtype=int)
    orders_prev = np.arange(len(wavelength_solution_ori)) - orderdiff                           # The old orders are getting new indexes, but that is just a change in numbers. above this is reached by order_arcline+orderdiff < 0
    orders_prev = orders_prev[ (orders_prev >= 0) & (orders_prev <= len(xlows)-1) ]             # no solution is needed for orders outside the new area
    for order_arcline in orders_prev:               # use the orders covered by the previous wavelength solution
        if order_arcline in ignoreorders:
            continue
        arc_lines_order_px = arc_lines_px[ (arc_lines_px[:,0] == order_arcline),:]    #[:,1] is px position, [:,2] is width of lines, [:,3] is height of the line
        # Only use lines for which the wavelength solution is defined
        #print arc_lines_order_px[ (arc_lines_order_px[:,1] - pxdiff - pxdifford * (order_arcline+orderdiff) <= 0), :]
        arc_lines_order_px = arc_lines_order_px[ (arc_lines_order_px[:,1] - pxdiff - pxdifford * (order_arcline+orderdiff) > 0) & \
                                                 (arc_lines_order_px[:,1] - pxdiff - pxdifford * (order_arcline+orderdiff) < specs[1]), : ]
        arc_lines_order_wl = np.polyval(wavelength_solution_2[order_arcline+orderdiff, 2:], arc_lines_order_px[:,1] - wavelength_solution_2[order_arcline+orderdiff, 1])
        for i, arc_line_wave in enumerate(arc_lines_order_wl):
            diff_catalog = reference_catalog[:,0] - arc_line_wave
            good_values = (abs(diff_catalog) <= max_diffs[order_arcline+orderdiff])
            #if order_arcline == 36:
            #    print order_arcline, arc_lines_order_px[i,1], arc_line_wave, reference_catalog[good_values,0], reference_catalog[good_values,1]
            for entry in reference_catalog[good_values,:]:
                arc_lines_wavelength.append([order_arcline, arc_lines_order_px[i,1], entry[0],0, arc_lines_order_px[i,3], entry[1], 1./(np.abs(arc_lines_order_px[i,2]-med_arc_width)+med_arc_width)])
                #if order_arcline == 36:            #testing
                #    print arc_lines_wavelength[-1]          #testing
                #print arc_lines_wavelength[-1]
    
    ## Fit polynoms in and between orders
    order_offset = wavelength_solution_ori[0][0]+orderdiff      # checked that it is "+orderdiff"
    opt_px_range = (1 - np.linspace(params['opt_px_range'][0], params['opt_px_range'][1], int(steps*5/10) )) * int(specs[1]/2)                         # The last steps should be done with full image
    #opt_px_range = np.concatenate(( opt_px_range, np.repeat([opt_px_range[-1]],steps-len(opt_px_range)) ))                                  # Add the missing entries to the array
    opt_px_range = np.concatenate(( opt_px_range, (1 - np.linspace(params['opt_px_range'][1], 10, steps-len(opt_px_range) )) * int(specs[1]/2) ))    # Add the missing entries to the array
    polynom_order_traces = np.linspace(params['polynom_order_traces'][0], params['polynom_order_traces'][1], int(steps*5/10), dtype=int)         # The last steps should be done with full number of orders
    polynom_order_traces = np.concatenate(( polynom_order_traces, np.repeat([polynom_order_traces[-1]],steps-len(polynom_order_traces)) ))
    polynom_order_intertraces = params['polynom_order_intertraces']
    if len(polynom_order_intertraces) == 1:
        polynom_order_intertraces.append(params['polynom_order_intertraces'][0])
    polynom_order_intertraces = np.linspace(polynom_order_intertraces[0], polynom_order_intertraces[1], int(steps*5/10), dtype=int)              # The last steps should be done with full number of orders
    polynom_order_intertraces = np.concatenate(( polynom_order_intertraces, np.repeat([polynom_order_intertraces[-1]],steps-len(polynom_order_intertraces)) ))
    max_diff_pxs = np.linspace(params['px_offset'][2], params['diff_pxs'], int(steps*5/10))
    max_diff_pxs = np.concatenate(( max_diff_pxs, np.linspace(params['diff_pxs'], np.nanmedian(arc_lines_px[:,2]), steps-len(max_diff_pxs)) ))      # In the last diff pixel shouldn't be bigger than the median width of the lines
    ref_catalog_brigthness = np.linspace(90, 0, int(steps*5/10))
    ref_catalog_brigthness = np.concatenate(( ref_catalog_brigthness, np.repeat([ref_catalog_brigthness[-1]],steps-len(ref_catalog_brigthness)) ))
    # Find the central pixels, but take into account that the number of orders could have changed, this ignores order_offset for the moment
    poly = np.polyfit(range(len(wavelength_solution_ori)), wavelength_solution_ori[:, 1], 2)         # central pixels from previous solution
    cen_px = np.polyval(poly, orders)
    cen_px = np.repeat( [int(specs[1]/2)], len(orders) )
    arc_lines_wavelength = np.array(arc_lines_wavelength)
    old_arc_lines_wavelength = arc_lines_wavelength
    #logdata = np.zeros((len(arc_lines_px)+1, 2+steps*iter_break*2))
    #logdata[1:,0:2] = arc_lines_px[:,0:2]
    old_px_range = opt_px_range[0]
    for step in tqdm(range(steps+1), desc='Finding the new wavelength solution'):
        if step < steps:                # Normal steps
            max_diff_px = max_diff_pxs[step]
            px_range = opt_px_range[step]
            polynom_order_trace = polynom_order_traces[step]
            polynom_order_intertrace = polynom_order_intertraces[step]
            good_values = (arc_lines_px[:,1] >= px_range+pxdiff) & (arc_lines_px[:,1] <= specs[1]-px_range+pxdiff) & (arc_lines_px[:,0] not in ignoreorders)        #+pxdiff to start with the lines closest to the center of the old wavelength solution   !!! add cen_px, but keep in mind that opt_px_range = 1 should include the whole order
            arc_lines_px_step = arc_lines_px[good_values,:]                             # Subset of the emission lines in the spectrum to which the fitting is done in this step
            weights1, weights3 = arc_lines_px_step[:,3], 1./(abs(arc_lines_px_step[:,2]-med_arc_width)+med_arc_width)     # height of the arc line, 1/width of the arc line
            reference_catalog = reference_catalog_full[ reference_catalog_full[:,1] >= np.percentile(reference_catalog_full[:,1], ref_catalog_brigthness[step]) ,:]   # use only the brightest lines
            #print reference_catalog.shape, reference_catalog_full.shape
            #print(step,max_diff_px,px_range,polynom_order_trace,polynom_order_intertrace,len(arc_lines_px_step), len(good_values), reference_catalog.shape)
        else:                           # final step: using the 
            iter_break = 1
            reference_catalog, reference_names = reference_lines_dict['reference_catalog'][-1], reference_lines_dict['reference_names'][-1]
        for iteration in range(iter_break):
                cen_pxs = arc_lines_wavelength[:,0] * np.nan
                for i,order in enumerate(orders):
                    cen_pxs[ (order == arc_lines_wavelength[:,0]) ] = cen_px[i]
                x = arc_lines_wavelength[:,1] - cen_pxs
                y = 1.0/(arc_lines_wavelength[:,0]+order_offset)
                # Get the weight of each line, normalised by the mdeian value: height of the arc line + log10(intensity of the reference line) + 1/width of the arc line
                #medians = [np.median(arc_lines_wavelength[:,4]), np.median(arc_lines_wavelength[:,5]), np.median(arc_lines_wavelength[:,6])]
                #for i in range(len(medians)):
                #    if medians[i] == 0:
                #        medians[i] = np.mean(arc_lines_wavelength[:,i+4])
                #    if medians[i] == 0:
                #        medians[i] = 1
                #weight = (arc_lines_wavelength[:,4]/medians[0]+arc_lines_wavelength[:,5]/medians[1]+arc_lines_wavelength[:,6]/medians[2])**2
                # Modify the weight so that each order has the same weight
                weight = arc_lines_wavelength[:,4] * arc_lines_wavelength[:,5] * arc_lines_wavelength[:,6] * np.abs(x)**2   #(polynom_order_trace-1)
                weight_scale = []
                for order in np.unique(arc_lines_wavelength[:,0]):
                    weight_scale.append( np.sum( weight[arc_lines_wavelength[:,0] == order] ))
                divisor = np.median(weight_scale)
                if divisor == 0:
                    divisor = 1.
                weight_scale = np.array(weight_scale)/(divisor+0.0)
                weight_scale[weight_scale > 10] = 10
                weight_scale[weight_scale < 0.1] = 0.1                              # Otherwise that might scale a badly covered order too much
                for i, order in enumerate(np.unique(arc_lines_wavelength[:,0])):
                    weight[arc_lines_wavelength[:,0] == order] /= (weight_scale[i]+0.0)
                weight=[]
                #good_values = (arc_lines_wavelength[:, 0] == 35)        #testing
                #print(weight[good_values], medians)                      #testing
                poly2d_params = polynomial_fit_2d_norm(x, y, arc_lines_wavelength[:,2], polynom_order_trace, polynom_order_intertrace, w=weight, w_range=1E5)    # Fit against the available/identified lines
                # Convert the pixel positions of the found lines into wavelengths
                cen_pxs = arc_lines_px_step[:,1] * np.nan
                for i,order in enumerate(orders):
                    cen_pxs[ (order == arc_lines_px_step[:,0]) ] = cen_px[i]
                x2 = arc_lines_px_step[:,1] - cen_pxs                         # The pixel position of all lines in this step
                y2 = 1.0/(arc_lines_px_step[:,0]+order_offset)
                arc_line_wave = polynomial_value_2d(x2, y2, polynom_order_trace, polynom_order_intertrace, poly2d_params)
                #print(len(x2),len(x), px_range, pxdiff, np.sum(arc_lines_px_step[:,0]==26))
                #print(np.vstack([y2,x2,arc_line_wave]).T[arc_lines_px_step[:,0] == 26,:])
                ## Resolution around the line
                arc_line_res = np.abs(polynomial_value_2d(x2+1, y2, polynom_order_trace, polynom_order_intertrace, poly2d_params) - \
                                      polynomial_value_2d(x2-1, y2, polynom_order_trace, polynom_order_intertrace, poly2d_params) )/2.
                # Avoid running off the wavelength solution with high orders -> fit linear solution and use it for the outer areas
                if old_px_range - px_range > 10:    # only use linear solution when more than 10 px more
                    #print("use linear values")
                    good_values_l = ( arc_lines_px_step[:,1] <= old_px_range+pxdiff + min(100,4*(old_px_range - px_range)) ) & \
                                    ( arc_lines_px_step[:,1] >= old_px_range+pxdiff) & ( ~np.isnan(x2) ) & ( ~np.isnan(y2) )                # lines on the left, covered by the old solution
                    good_values_r = ( arc_lines_px_step[:,1] >= specs[1]-old_px_range+pxdiff - min(100,4*(old_px_range - px_range)) ) & \
                                    ( arc_lines_px_step[:,1] <= specs[1]-old_px_range+pxdiff ) & ( ~np.isnan(x2) ) & ( ~np.isnan(y2) )      # lines on the right, covered by the old solution
                    #print(px_range, old_px_range, np.min(x2[good_values_l]), np.max(x2[good_values_l]), np.min(x2[good_values_r]), np.max(x2[good_values_r]), 'bla', \
                    #      np.min(arc_lines_px_step[good_values_l,1]), np.max(arc_lines_px_step[good_values_l,1]), np.min(arc_lines_px_step[good_values_r,1]), np.max(arc_lines_px_step[good_values_r,1]), 'ble', \
                    #      np.min(cen_pxs[good_values_l]), np.max(cen_pxs[good_values_l]), np.min(cen_pxs[good_values_r]), np.max(cen_pxs[good_values_r]), 'bli', \
                    #      old_px_range+pxdiff, old_px_range+pxdiff + 2*(old_px_range - px_range), specs[1]-old_px_range+pxdiff - 2*(old_px_range - px_range), specs[1]-old_px_range+pxdiff, 'blu', \
                    #      specs[1],old_px_range,pxdiff,2*(old_px_range - px_range), 'blo', len(x2),len(good_values_l))
                    if np.sum(good_values_l) > 10:
                        poly2d_params_l = polynomial_fit_2d_norm(x2[good_values_l], y2[good_values_l], arc_line_wave[good_values_l], 1, polynom_order_intertrace)    # linear fit in the outsides of the order
                        good_values_l = ( arc_lines_px_step[:,1] < old_px_range+pxdiff )   # lines on the left, only covered by the new solution
                        #print(np.vstack([ arc_lines_px_step[good_values_l,1], arc_lines_px_step[good_values_l,0], arc_line_wave[good_values_l], 
                        #                  arc_line_wave[good_values_l] - polynomial_value_2d(x2[good_values_l], y2[good_values_l], 1, polynom_order_intertrace, poly2d_params_l) ]).T)
                        arc_line_wave[good_values_l] = polynomial_value_2d(x2[good_values_l], y2[good_values_l], 1, polynom_order_intertrace, poly2d_params_l)      # replace the wavelength with the linear fit
                    if np.sum(good_values_r) > 10:
                        poly2d_params_r = polynomial_fit_2d_norm(x2[good_values_r], y2[good_values_r], arc_line_wave[good_values_r], 1, polynom_order_intertrace)    # linear fit in the outsides of the order
                        good_values_r = ( arc_lines_px_step[:,1] > specs[1]-old_px_range+pxdiff )   # lines on the right, only covered by the new solution
                        #print(np.vstack([ arc_lines_px_step[good_values_r,1], arc_lines_px_step[good_values_r,0], arc_line_wave[good_values_r], 
                        #                  arc_line_wave[good_values_r] - polynomial_value_2d(x2[good_values_r], y2[good_values_r], 1, polynom_order_intertrace, poly2d_params_r) ]).T )
                        arc_line_wave[good_values_r] = polynomial_value_2d(x2[good_values_r], y2[good_values_r], 1, polynom_order_intertrace, poly2d_params_r)      # replace the wavelength with the linear fit
                    old_px_range = opt_px_range[step]
                    
                ## Resolution in the central area
                y_lin = y2[ sorted( np.unique(y2, return_index=True)[1] ) ]      # np.unique sorts the values -> get the indexes and sort them in order to avoid mixing up different lines
                x_lin1 = np.repeat( np.percentile(x, 25), len(y_lin) )          # for each order with found emission lines get the 25 percentile value of the pixel position
                x_lin2 = np.repeat( np.percentile(x, 75), len(y_lin) )
                res_cen = (polynomial_value_2d(x_lin2, y_lin, polynom_order_trace, polynom_order_intertrace, poly2d_params) - \
                           polynomial_value_2d(x_lin1, y_lin, polynom_order_trace, polynom_order_intertrace, poly2d_params) ) / (x_lin2 - x_lin1)       # Resolution for all orders in the center (25 to 75 percentile)
                # Identify the lines again, asign the reference line with the highest intensity
                orders_prev = list(orders_prev) + [min(orders_prev)-1, max(orders_prev)+1]      # Add only one order at the time
                arc_lines_wavelength = []
                oldorder = -1E5
                for i in range(arc_line_wave.shape[0]):
                    if arc_lines_px_step[i,0] not in orders_prev:                       # Ignore orders which are too far off
                        continue
                    if arc_lines_px_step[i,0] != oldorder:                              # Each order the reference lines can be used again
                        assigned_reflines = []
                        oldorder = arc_lines_px_step[i,0]
                    # Ignore lines, for which the resolution is unphysical
                    index = np.where( np.abs( 1./(arc_lines_px_step[i,0]+order_offset) - y_lin ) < 0.00001 )[0]      # Find the index of the order in res_cen, orders might be missing
                    if len(index) != 1:
                        logger('Error: that should not have happened, check code around line 4700. Index contains {0} entries'.format(len(index)))
                    index = index[0]                                                    # make the list into an integer
                    if abs(arc_line_res[i] - res_cen[index]) > res_cen[index] * 0.25:   # The variation in one order should be less than 15% (exohspec), 25% MRES
                        #1#print "ignore the line", arc_line_res[i], res_cen[int(arc_lines_px_step[i,0])], arc_lines_px_step[i,0], arc_lines_px_step[i,1]
                        #if arc_lines_px_step[i,0] == 26:
                        #    print arc_lines_px_step[i,0], arc_lines_px_step[i,1], arc_line_wave[i], "resolution too far off"
                        continue                # Ignore lines which differ a lot, hopefully other lines will fix that
                    # Wavelength difference to the catalog lines
                    diff_catalog = reference_catalog[:,0] - arc_line_wave[i]    # Diff between catalog wavelength and fitted wavelength
                    good_pos = np.where( np.abs(diff_catalog) <= max_diff_px*arc_line_res[i] )[0]                   # Lines which are close to the fitted wavelength
                    #print i, arc_line_wave[i], len(good_pos), max_diff_px, arc_line_res[i], max_diff_px*arc_line_res[i]
                    # Remove already identified lines, can't be done earlier, as the index in good_pos is important
                    if only_one_line == True:
                        good_pos = np.array( [j for j in good_pos if j not in assigned_reflines] )                  # Only keep the good_pos, which are not in assigned_reflines
                    if len(good_pos) > 1 and only_one_line == True:
                        pos_max = np.argmax(reference_catalog[good_pos,1] * np.sqrt(np.abs(diff_catalog[good_pos])) )          # Get the cataloge line with highes intensitiy or the one which is closer
                        good_pos = [good_pos[pos_max]]
                    for j in good_pos:
                        #logdata[pos_logdata, 2+(step*iter_break+iteration)*2+1] = reference_catalog[j,0]
                        weight2 = np.log10(reference_catalog[j,1])     # log10(intensity of the reference line)
                        arc_lines_wavelength.append([ arc_lines_px_step[i,0], arc_lines_px_step[i,1], reference_catalog[j,0], diff_catalog[j], 
                                                      weights1[i], weight2, weights3[i], arc_line_res[i], j, arc_lines_px_step[i,2] ])
                        """arc_lines_wavelength:
                        0: order
                        1: pixel
                        2: wavelength from reference cataloge
                        3: wavelength difference between fit and reference catalogue
                        4: weight: height of the line (in ThAr spectrum)
                        5: weight: log10 of intesity of the reference line
                        6: weight: 1/width of the line (in ThAr spectrum), the lines with width closest to the median are weighted most
                        7: resolution at this place (1/2 of 2-pixel (one pixel))
                        8: index in reference catalogue
                        9: width of the line (in ThAr spectrum)
                        """
                        assigned_reflines.append(j)                                                         # Don't use the line again
                        #if arc_lines_px_step[i,0] == 43:            #testing
                        #    print arc_lines_wavelength[-1]          #testing
                    #if arc_lines_px_step[i,0] == 26:
                    #    print arc_lines_px_step[i,0], arc_lines_px_step[i,1], arc_line_wave[i], len(good_pos), min(np.abs(diff_catalog)), len(arc_lines_wavelength), reference_catalog[j,0]
                    
                arc_lines_wavelength = np.array(arc_lines_wavelength)   
                # arc_lines_wavelength: order, px in spectrum, wavelength in reference catalog, difference in wavelength between wavelength solution and reference line, 
                #                       3x weights, 1px resolution at that position, index in refence catalog, width of the line
                if arc_lines_wavelength.shape[0] < 10:
                    logger('Error: Something went wrong, only {0} lines have been identified. '.format(arc_lines_wavelength.shape[0])+\
                            'Please compare the positon of orders in the current folder with the folder from which the previous wavelength solution was used. '+\
                            'If necessary adjust the parameters order_offset, px_offset, px_offset_order, and/or resolution_offset_pct. '+\
                            'Decreasing the number of orders in the polynomial fit (parameters polynom_order_traces and polynom_order_intertraces) might also help to fix the issue.')
                goodvalues, p = sigmaclip(arc_lines_wavelength[:,3], nc=0, ll=sigma, lu=sigma, repeats=20)    # poly_order=0 to sigma clip around the average
                arc_lines_wavelength = arc_lines_wavelength[goodvalues,:]
                cen_px, p_cen_px = find_real_center_wavelength_solution(order_offset, orders, [], cen_px, polynom_order_trace, polynom_order_intertrace, poly2d_params)
                
                #print max_diff_px, arc_lines_wavelength.shape[0],old_arc_lines_wavelength.shape[0],arc_lines_wavelength.shape[0] == old_arc_lines_wavelength.shape[0],np.mean(arc_lines_wavelength[:,3]), np.mean(old_arc_lines_wavelength[:]),np.mean(arc_lines_wavelength[:,3]) == np.mean(old_arc_lines_wavelength[:]), np.min(arc_lines_wavelength[:,3]), np.max(arc_lines_wavelength[:,3])
                if arc_lines_wavelength.shape[0] == old_arc_lines_wavelength.shape[0]:
                    if np.mean(arc_lines_wavelength[:,3]) == np.mean(old_arc_lines_wavelength[:]):
                        break
                old_arc_lines_wavelength = copy.deepcopy(arc_lines_wavelength[:,3])
                
                """#1# # Only for printing the progress: Plotting each step
                print step, iteration, max_diff_px, px_range, polynom_order_trace, polynom_order_intertrace, np.std(arc_lines_wavelength[:,3], ddof=(polynom_order_trace+polynom_order_intertrace+1))
                max_number_reflines = 0             # will be needed later in order to save the data correctly into a fits file
                for order in orders:
                    max_number_reflines = max(max_number_reflines, len(arc_lines_wavelength[(arc_lines_wavelength[:,0]==order),2]))
                # Transform the wavelength solution into the old wavelength solution
                wavelength_solution, wavelength_solution_arclines = [], []
                xarr = np.arange(0,specs[1], 0.1)
                for i,order in enumerat(orders):
                    yarr = polynomial_value_2d(xarr-cen_px[i], 1.0/(order+order_offset), polynom_order_trace, polynom_order_intertrace, poly2d_params)
                    polyfit = np.polyfit(xarr-cen_px[i], yarr, polynom_order_trace)      #lambda from px
                    solution = [order_offset+order, cen_px[i] ]
                    solution = np.append(solution, polyfit, axis=0)
                    wavelength_solution.append(solution)
                    reflines = list(arc_lines_wavelength[(arc_lines_wavelength[:,0]==order),2])     # Wavelength
                    for i in range(max_number_reflines - len(reflines)):                            # Fill up with zeros in order to create an array
                        reflines.append(0)
                    wavelength_solution_arclines.append(reflines)     #The wavelength of the reference lines
                    
                    #print order,polynomial_value_2d(0, 1.0/(order+order_offset), polynom_order_trace, polynom_order_intertrace, poly2d_params), solution
                    #diff_fit = yarr-np.polyval(polyfit, xarr-cen_px[i])        #<1e-10
                    #print [len(diff_fit), np.std(diff_fit, ddof=polynom_order_trace+1),min(diff_fit),max(diff_fit)]
                wavelength_solution = np.array(wavelength_solution)
                wavelength_solution_arclines = np.array(wavelength_solution_arclines)
                plot_wavelength_solution_spectrum(params, spectrum, spectrum, params['logging_arc_line_identification_spectrum'], wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names)
                raw_input('enter')"""
                
                #print arc_lines_wavelength.shape, poly2d_params
                """ # Only for printing the progress
                x,y,l=[],[],[]
                for order in orders:
                    x.append(arc_lines_wavelength[(arc_lines_wavelength[:,0]==order),1])
                    y.append(arc_lines_wavelength[(arc_lines_wavelength[:,0]==order),3])
                    l.append('{0}'.format(order))
                if iteration == 100:
                    plot_img_spec.plot_points(x,y,l,'path',show=True, x_title='Pixel', y_title='Residual')          """
    
    #for line in logdata:
    #    os.system('echo "{0}" >> logdata'.format("\t".join(line.astype(str))) )
    # Get the real order numbers using the grating equation: wavelength propotional to 1/(real_order) -> y(order)=(order_offset_new+order)*central_wavelength should have the smallest slope
    order_offset_old = int(order_offset)
    cenwave = polynomial_value_2d(orders*0, 1.0/(orders+order_offset_old), polynom_order_trace, polynom_order_intertrace, poly2d_params)    #old order_offset
    if search_order_offset:
        order_offset = find_order_offset(orders, cenwave)
        if 0.0 in order_offset+orders:            # Definitively the wrong order offset
            order_offset = order_offset_old
    p_real_cent = np.polyfit(1.0/(order_offset+orders), cenwave, 1)
    diff_real_cent = cenwave-np.polyval(p_real_cent, 1.0/(order_offset+orders))
    cen_px, p_cen_px = find_real_center_wavelength_solution(order_offset, orders, cenwave, cen_px, polynom_order_trace, polynom_order_intertrace, poly2d_params)
    # Redo the 2d polynomial fit (it should also re-create the arc_lines_wavelength table, but most times the difference is 0)
    cen_pxs = arc_lines_wavelength[:,1] * np.nan
    for i,order in enumerate(orders):
        cen_pxs[ (order == arc_lines_wavelength[:,0]) ] = cen_px[i]     
        
    # Fit the solution again
    poly2d_params, new_waves = fit_wavelengths_solution_2d(arc_lines_wavelength, cen_pxs, order_offset, polynom_order_trace, polynom_order_intertrace)
    arc_lines_wavelength[:,3] = arc_lines_wavelength[:,2] - new_waves               # new differences
    
    if len(arc_lines_wavelength) <= polynom_order_trace+polynom_order_intertrace:
        logger('Error: Not enough lines are remaining for the fit. Please check the parameters "order_offset", "px_offset", and "px_offset_order".')
    
    """# Testing
    xarr = np.repeat( [np.arange(max( arc_lines_wavelength[:,1] ))], len(orders), axis=0)
    yarr = []
    for i,order in enumerate(orders):
        wave = polynomial_value_2d(xarr[0]-cen_px[i], 1.0/(order+order_offset), polynom_order_trace, polynom_order_intertrace, poly2d_params)
        p_lin = np.polyfit(xarr[0]-cen_px[i], wave, 1)
        yarr.append( wave - np.polyval(p_lin, xarr[0]-cen_px[i]) +i*5 )         # Without linear part, offset between orders, can't compare identified lines
        #yarr.append( wave )                                                     # Complete solution, to compare with identified lines
    yarr = np.array(yarr)
    xp = list(xarr) + [arc_lines_wavelength[:,1]]
    yp = list(yarr) + [arc_lines_wavelength[:,2]]
    lp = list(np.array(orders,dtype=str)) + ['lines']
    plot_img_spec.plot_points(xp ,yp , lp,'dummy',show=True, title='', x_title='Pixel', y_title='wavelength [Angstrom]')
    good_values = ( np.abs( (xarr[0,:]/50).astype(int) - xarr[0,:].astype(float)/50.) < 1E-5 )       # get every 50th pxel
    l = xarr[0,good_values]                                         # x-values will be the lable
    xarr = np.repeat([orders], len(l), axis=0)                     # orders will be the new xarr
    yarr = yarr[:,good_values].T
    plot_img_spec.plot_points(xarr ,yarr,l,'dummy',show=True, title='', x_title='Orders', y_title='wavelength [Angstrom]')"""
    
    # Log the results
    # print np.vstack([ arc_lines_wavelength[:,2], arc_lines_wavelength[:,3], arc_lines_wavelength[:,-3], arc_lines_wavelength[:,2]/arc_lines_wavelength[:,3] ]).T
    #np.set_printoptions(threshold=np.inf)
    #print np.sort(arc_lines_wavelength[:,3])
    #np.set_printoptions(threshold=1000)
    #print np.histogram(arc_lines_wavelength[:,3], 20), arc_lines_wavelength.shape
    
    # Remove outliers
    len_orig = arc_lines_wavelength.shape[0]
    logger('Info: {0} before final sigma clipping'.format(len_orig), show=False)
    good_values = ( np.abs(arc_lines_wavelength[:,3]) < 100 )               # Initialise
    for ii in range(20):                                                    # make this graphical!!
        # Remove the most scattered lines (highest 0.5% and lowest 0.5%
        """todel = np.argsort(arc_lines_wavelength[:,3])           # diff in Ang between catalog wavelength and solution
        todel = np.delete(todel, np.arange(int(len(todel)*0.5/100), int(len(todel)*99.5/100)), axis=0)
        arc_lines_wavelength = np.delete(arc_lines_wavelength, todel, axis=0)
        cen_pxs = np.delete(cen_pxs, todel, axis=0)
        print '1percent', len(todel)"""
        if np.std(arc_lines_wavelength[good_values,3]) > 0.01:
            sigma = 2.5
        else:
            sigma = 2.8
        std_diff = np.std(arc_lines_wavelength[good_values,3], ddof=(polynom_order_trace+polynom_order_intertrace+1))             # real standard deviation, asume values +-1 to convice yourself, no abs
        good_values = ( (-sigma*std_diff <= arc_lines_wavelength[:,3]) & (arc_lines_wavelength[:,3] <= sigma*std_diff) )                                        # only works if average is 0
        # Fit the solution again, using only the good values
        poly2d_params, new_waves = fit_wavelengths_solution_2d(arc_lines_wavelength[good_values,:], cen_pxs[good_values], order_offset, polynom_order_trace, polynom_order_intertrace)
        # Calculate the wavelengths again, using all values
        new_waves = polynomial_value_2d(arc_lines_wavelength[:,1]-cen_pxs, 1.0/(arc_lines_wavelength[:,0]+order_offset), polynom_order_trace, polynom_order_intertrace, poly2d_params)
        arc_lines_wavelength[:,3] = arc_lines_wavelength[:,2] - new_waves               # new differences
        if sum(good_values) < len_orig*2/3:
            break
    arc_lines_wavelength = arc_lines_wavelength[good_values,:]
    cen_pxs = cen_pxs[good_values]
            
    #print np.histogram(arc_lines_wavelength[:,3], 20), arc_lines_wavelength.shape
    avg_diff_fit = np.mean(np.abs(arc_lines_wavelength[:,3]))           # Diff between catalog wavelength and fitted wavelength
    std_diff_fit = np.std(arc_lines_wavelength[:,3], ddof=(polynom_order_trace+polynom_order_intertrace+1))             # real standard deviation, asume values +-1 to convice yourself
    # Resolution from the precission of the fit
    std_R_fit = 1.0/np.std(arc_lines_wavelength[:,3]/(arc_lines_wavelength[:,2]+0.0), ddof=(polynom_order_trace+polynom_order_intertrace+1))    # lambda/d_lambda
    # Resolution using the Gaussian width of the arc lines
    R_gauss     = arc_lines_wavelength[:,2]/(1.*arc_lines_wavelength[:,9]*arc_lines_wavelength[:,7])    # lambda/d_lambda ; -1=9 is in Gauss width in px, -3=7 is resolution
    arc_lines_wavelength = np.hstack(( arc_lines_wavelength, (R_gauss/2.35482).reshape(( R_gauss.shape[0], 1 )) ))      # Resolution
    R_gauss     = arc_lines_wavelength[:,2]/(1.*arc_lines_wavelength[:,9]*arc_lines_wavelength[:,7])    # lambda/d_lambda ; -1=9 is in Gauss width in px, -3=7 is resolution
    good_values, p = sigmaclip(R_gauss, nc=0, ll=3, lu=3)   #y, order of polynom (if 0 than average -> x doesn't matter), sigma, sigma
    R_gauss     = R_gauss[good_values]
    avg_R_gauss = np.mean(R_gauss)
    std_R_gauss = np.std(R_gauss, ddof=1)
    avg_R_fwhm, std_R_fwhm = avg_R_gauss/2.35482, std_R_gauss/2.35482
    # 2px Resolution (using only the identified arc lines
    R_2px     = arc_lines_wavelength[:,2]/(2.*arc_lines_wavelength[:,7])    # lambda/d_lambda
    avg_R_2px = np.mean(R_2px)
    std_R_2px = np.std(R_2px, ddof=1) 
    text = 'Info: used {0} lines. The standard deviation (using {8} degrees of freedom) of the residuals between the lines and the fit is {1} Angstrom. '+\
                  'The FWHM of the emission lines results in an R = {3} +- {4}. The 2-pixel resolution (around the identified lines) is R = {5} +- {6}. '+\
                  'The deviation to the line fit converts into a resolution R = {2}. The average of the abs of the residuals is {7} Angstrom. '
    logger(text.format(arc_lines_wavelength.shape[0], round_sig(std_diff_fit,3), int(round_sig(std_R_fit,4)), int(avg_R_fwhm), int(std_R_fwhm),
                       int(avg_R_2px), int(std_R_2px), round_sig(avg_diff_fit,3), polynom_order_trace+polynom_order_intertrace ))
    p_cen_px = np.round(p_cen_px,3)
    text = 'Info: A 2d polynom fit with {0} orders in dispersion direction (along the traces) and {1} orders in cross-dispersion direction was used. '+\
                  'With this solution, the offset between aperture and real orders is {2}. To fulfil the grating equation the central pixel of the individual orders needs to be {5} + {6}*order + {7}*order**2. '+\
                  'With this values the standard deviation of the residuals between the central wavelengths and the grating equation is {3} Angstrom. Using the original solution gives an offset of {4}.'
    logger(text.format(polynom_order_trace, polynom_order_intertrace, int(order_offset), round_sig(np.std(diff_real_cent, ddof=len(p_real_cent)),3), order_offset_old, p_cen_px[2], p_cen_px[1], p_cen_px[0] ))
    
    # Create a wavelength solution
    wavelength_solution2d = np.array( [polynom_order_trace, polynom_order_intertrace, np.mean(cen_px), order_offset] + list(poly2d_params) )
    text = 'Info: Wavelenght solution in 2d (for pixel and order at the same time) is [No of orders1, No of orders2, mean central pixel, '+\
                  'offset to real order, parameters of the 2d polynomial(px^0*m^0, px^1*m^0, px^0*m^1, px^2*m^0, px^1*m^1, px^0*m^2, ....)]: \n'+\
                  '{0}'
    logger(text.format(wavelength_solution2d), show=False)
    # Log into text file and make it ready to copy to pixel_to_wavelength
    printarrayformat = ['%1.1i', '%1.1i', '%3.1f', '%9.4f', '%6.4f']
    printarray = copy.deepcopy(arc_lines_wavelength[:,[0,0,1,2,3]])
    printarray[:,1] += int(order_offset)
    logger('order\treal_o\tpixel\twavelength\tdwavel', show=False, printarrayformat=printarrayformat, printarray=printarray, logfile=params['logging_identified_arc_lines'])
    
    # Create an image of the gaussian widths of the identified lines
    plot_wavelength_solution_width_emmission_lines(params['logging_resolution_form'], specs, arc_lines_wavelength, [0,1,10], title='Resolution')   # order, wavelength, gauss width
    plot_wavelength_solution_width_emmission_lines(params['logging_em_lines_gauss_width_form'], specs, arc_lines_wavelength, [0,1,9])   # order, wavelength, gauss width
    plot_hist_residuals_wavesol(params['logging_arc_line_identification_residuals_hist'], arc_lines_wavelength, [0,1,2,3] )             # order, pixel, wavelength, residuals
    bisector_measurements_emission_lines(params['logging_em_lines_bisector'], spectrum, arc_lines_wavelength, [0,1,9])         # order, pixel, width
        
    x, w, y, l = [], [], [], []
    for order in orders:
        inorder = ( arc_lines_wavelength[:,0]==order )
        x.append(arc_lines_wavelength[inorder,1])    # Pixel
        w.append(arc_lines_wavelength[inorder,2])    # Wavelength
        y.append(arc_lines_wavelength[inorder,3])    # Residuals
        l.append('{0}'.format(order))
    text = ('Residuals of the identificated arc lines: identified {0} lines when using a 2d polynom fit with {1} (dispersion) '+\
            'and {2} (cross-dispersion) orders').format(arc_lines_wavelength.shape[0], polynom_order_trace, polynom_order_intertrace)
    fname_px   = params['logging_arc_line_identification_residuals'].rsplit('.',1)[0] + '_px.'   + params['logging_arc_line_identification_residuals'].rsplit('.',1)[-1]
    fname_wave = params['logging_arc_line_identification_residuals'].rsplit('.',1)[0] + '_wave.' + params['logging_arc_line_identification_residuals'].rsplit('.',1)[-1]
    plot_img_spec.plot_points(x,y,l,[fname_px],   show=False, title=text, x_title='Pixel', 
                              y_title='Residual (O-C) [Angstrom]', marker=['o','s','*','P','^','v','>','<','x'])
    plot_img_spec.plot_points(w,y,l,[fname_wave], show=False, title=text, x_title='Wavelength [Angstrom]', 
                              y_title='Residual (O-C) [Angstrom]', marker=['o','s','*','P','^','v','>','<','x'])

    # Transform the wavelength solution into the old wavelength solution
    wavelength_solution, wavelength_solution_arclines, line_stats = transform_wavelength_solution_2d_to_n1d(specs[0], specs[1], 
                                polynom_order_trace, polynom_order_intertrace, poly2d_params, order_offset, cen_px, arc_lines_wavelength)
    
    statistics_arc_reference_lines(arc_lines_wavelength, [0,8,9,2], reference_names, wavelength_solution, xlows, xhighs)
        
    if std_diff_fit > 2*max(wavelength_solution[:,-2]) or std_diff_fit < 1E-8:        # if residuals are bigger than 1px or unreasonable small
        plot_wavelength_solution_spectrum(params, spectrum, spectrum, params['logging_arc_line_identification_spectrum'], wavelength_solution, 
                                          wavelength_solution_arclines, reference_catalog, reference_names, plot_log=True)
        if 'master_cal2_l_filename' in params.keys():
            text = '{0} in the current folder'.format(params['master_cal2_l_filename'])
        else:
            text = 'files listed in the parameter cal2_l_rawfiles'
        logger(('Error: The wavelength solution seems wrong. Please check the parameters "order_offset", "px_offset", and "px_offset_order".' + \
               '\n\t\tIt might be useful to compare the file with the emission lines ({0}) and ' + \
               'the folder with the previous wavelength solution (see parameter "original_master_wavelensolution_filename")' +\
               '\n\t\tThe results of the identification can be seen in {1}.').format(text, params['logging_arc_line_identification_spectrum']))
    
    # See the results
    if show_res:
        correlate_px_wave_result_UI(spectrum, arc_lines_wavelength, reference_catalog, arc_lines_px, reference_names, wavelength_solution, wavelength_solution_arclines, adjust=[0.07,0.93,0.94,0.06, 1.0,1.01])
    
    return dict(wavesol=wavelength_solution, wavesol2d=wavelength_solution2d , reflines=wavelength_solution_arclines, linestat=line_stats)

def statistics_arc_reference_lines(data, positions, reference_names, wavelength_solution, xlows, xhighs, show=True):
    """
    positions: array with the columns in data for
                        - order
                        - index of reference line
                        - width of the data
                        - wavelength of the reference line
    """

    uniq_refnames = np.unique(reference_names)      # get the types of emission lines, e.g. Th I, Ne II
    uniq_orders = np.unique(data[:,positions[0]])   # get all available orders

    result = []
    for order in uniq_orders:
        index_o = (data[:,positions[0]] == order)
        data_o = data[index_o,:]                    # only data for the current order
        for refname in uniq_refnames:
            index_r = []
            for i,entry in enumerate(data_o[:,positions[1]].astype(int)):
                if reference_names[entry] == refname:
                    index_r.append(i)               # only indexes for the current name of the emission line
            if index_r != []:
                order = int(order)
                cenwave = wavelength_solution[order,-1]     # is the same as cenwave = np.polyval(wavelength_solution[order,2:], 0)
                minwave = np.polyval(wavelength_solution[order,2:], xlows[order]-wavelength_solution[order,1])
                maxwave = np.polyval(wavelength_solution[order,2:], xhighs[order]-wavelength_solution[order,1])
                min_refline = min(data_o[index_r,positions[3]])
                max_refline = max(data_o[index_r,positions[3]])
                range_refline = max(data_o[:,positions[3]])-min(data_o[:,positions[3]])
                std = np.nan
                if len(data_o[index_r,positions[2]]) > 1:                   # Message if np.std on array with one entry: "python2.7/site-packages/numpy/core/_methods.py:135: RuntimeWarning: Degrees of freedom <= 0 for slice"
                    std = np.std(data_o[index_r,positions[2]], ddof=1)      # ... and "python2.7/site-packages/numpy/core/_methods.py:127: RuntimeWarning: invalid value encountered in double_scalars"
                result.append([ order, cenwave, minwave, maxwave, maxwave-minwave, wavelength_solution[order,-2], refname, len(data_o[index_r,positions[2]]), 
                                np.average(data_o[index_r,positions[2]]), std, min_refline, max_refline, range_refline ])
                #print result[-1]
    # Adding the global information for each type of line (ThI, ArII, ...)
    for refname in uniq_refnames:
        index_r = []
        for i,entry in enumerate(data[:,positions[1]].astype(int)):
            if reference_names[entry] == refname:
                index_r.append(i)
        if index_r != []:
            min_refline = min(data[index_r,positions[3]])
            max_refline = max(data[index_r,positions[3]])
            result.append([ -1, -1, -1, -1, -1, -1, refname, len(data[index_r,positions[2]]), 
                            np.average(data[index_r,positions[2]]), np.std(data[index_r,positions[2]], ddof=1), min_refline, max_refline, max_refline-min_refline ])
    printarrayformat = ['%1.1i', '%3.1f', '%3.1f', '%3.1f', '%3.1f', '%5.3f', '%s', '%1.1i', '%4.2f', '\t%4.2f', '\t%4.1f', '\t%4.1f', '\t%4.1f']
    logger('\t\tapert\tcenwave\tminwave\tmaxwave\tranwave\tAng/px\tname\tnumber\tgausswidth_avg\tgausswidth_std\tmin_refline\tmax_refline\trange_reflines_whole_order',
           printarrayformat=printarrayformat, printarray=result, show=show)

def find_order_offset(orders, cenwave):
    """
    Uses the 2d wavelength solution to determine the real order offset between the [orders] and the grating equation.
    Get the real order numbers using the grating equation: wavelength is propotional to 1/(real_order) -> y(order)=(order_offset_new+order)*central_wavelength should have the smallest slope
    :param orders: 1d array of integers, normally from 0 to number of orders minus one
    :param cenwave: 1d array of floats, same length as orders. The central wavelength of each order
    :return order_offset: integer, offset so that orders fulfil the grating equation
    """
    slope_real_order = []
    for order_offset in range(-200,201):
        if np.min(np.abs(order_offset+orders)) < 1:     # Order 0 shouldn't be covered by Echelle Spectrograph
            continue
        y = (order_offset+orders)*cenwave   #(order_offset+order)*central_wavelength
        p = np.polyfit(orders, y, 1)
        slope_real_order.append( [abs(p[0]), order_offset ] + list(p) )
    slope_real_order = np.array(slope_real_order)
    order_offset = slope_real_order[ np.argmin(slope_real_order[:,0]), 1]
    #print slope_real_order[ np.argmin(slope_real_order[:,0]), :]
    return order_offset
    
def find_real_center_wavelength_solution(order_offset, orders, cenwave, cen_px, polynom_order_trace, polynom_order_intertrace, poly2d_params):
    """
    Find the real central pixel of each order using the wavelength solution
    :param order_offset: integer, offset so that orders fulfil the grating equation best
    :param orders: 1d array of integers, normally from 0 to number of orders minus one
    :param cenwave: 1d array of floats or empty, same length as orders. The central wavelength of each order
    :param cen_px: 1d array of floats, the original central pixels
    :param polynom_order_trace: integer, number of orders in dispersion direction to fit the 2d polynomial fit for the wavelength solution
    :param polynom_order_intertrace: integer, number of orders in cross-dispersion direction to fit the 2d polynomial fit for the wavelength solution
    :param poly2d_params: resulting paramers of the 2d polynomial fit for the wavelength solution
    :return cen_px: 1d array of floats, the new central pixels in order to fulfil the grating equation better
    """
    if len(cenwave) != len(orders):
        cenwave = polynomial_value_2d(orders*0, 1.0/(orders+order_offset), polynom_order_trace, polynom_order_intertrace, poly2d_params)    #old order_offset
    p_real_cent = np.polyfit(1.0/(order_offset+orders), cenwave, 1)         # cenwave = p[0] * 1/n + cenwave_offset
    real_cent = p_real_cent[0] * 1.0/(order_offset+orders)                  # real central wavelength for each order ignoring the offset
    # Find the pixel to real_cent
    xarr = np.arange(-1, 1.01, 0.05)                                      # only allow a shift of 50 px at a time !! decreased to 1 as sometimes it really run off, when it shouldn't
    for i, order in enumerate(orders):
        wave = polynomial_value_2d(xarr, 1.0/(order+order_offset), polynom_order_trace, polynom_order_intertrace, poly2d_params)       # xarr is around 0 -> cen_px doesn't need to be applied
        pos = np.argmin( np.abs( wave - real_cent[i] ) )                    # the index with the best matching central wavelength
        cen_px[i] += xarr[pos]
    cen_px_nofit = copy.deepcopy(cen_px)
    # fit a 2d polynomial against cen_px -> in case of runaway values and this should be physical
    poly = np.polyfit(orders, cen_px, 2)
    #print cen_px-np.polyval(poly, orders)
    cen_px = np.polyval(poly, orders)
    #print np.vstack([real_cent, cenwave, cen_px, cen_px_nofit]).T[ [0,1,21,41,61,81,-2,-1], :], p_real_cent, poly
    return cen_px, poly

"""def read_fit_wavelength_solution(params, filename, spec):
    "#""
    Reads the file with pixel - wavelength corellation and fits a 2d array against it to create a rough wavelength solution
    :params filename: string to a textfile with the following data: order (starting at 0), real order, pixel, wavelength. If one of the information is missing, this line will be skipped.
    "#""
    arc_lines_wavelength = []
    with open(filename, 'r') as file:
        for line in file:
            line = line[:-1].split('\t')
            if len(line) < 4:
                continue
            if line[0]!='' and line[1]!='' and line[2]!='' and line[3]!='':
                # get the order, real order, pixel, and wavelength
                try:
                    arc_lines_wavelength.append([ int(line[0]), int(line[1]), float(line[2]), float(line[3]),0 ])
                except:
                    print( 'Problems to convert to int/float:', line )
    if len(arc_lines_wavelength) == 0:
        logger(('Error: No useful information was found in {0}. Please make sure the entries in the file contain of the tab-separated '+\
                '(exactly one tab between each column) values: order (starting at 0), real order, pixel, wavelength').format(filename))
    arc_lines_wavelength = np.array(arc_lines_wavelength)
    
    wavelength_solution, wavelength_solution_arclines = fit_basic_wavelength_solution(params, arc_lines_wavelength, spec, filename)
    
    return wavelength_solution, wavelength_solution_arclines"""
    
def fit_basic_wavelength_solution(params, arc_lines_wavelength, spec, filename):
    if 0.0 in arc_lines_wavelength[:,1]:
        logger('Error: The second coloumn (real order) in {0} was not set correctly as it contains a zero. '+\
               'Please use the (estimated) real order (from grating equation).'.format(filename))
    orig_lines = arc_lines_wavelength.shape[0]
    
    order_offset = arc_lines_wavelength[:,1] - arc_lines_wavelength[:,0]
    orders = np.arange(spec.shape[0], dtype=int)
    if np.std(order_offset, ddof=1) != 0:
        logger('Error: There is an inconsistency between coloumn 1 (order/aperture) and coloumn 2 (real/physical order) in {0}. '+\
               'Please check that the difference between these two columns is always the same.'.format(filename))
    order_offset = order_offset[0]
    
    polynom_order_trace = max(2, max(params.get('tmp_polynom_order_traces', params['polynom_order_traces'])) )                  # Try to read the tmp values and if they don't exist, use the standard values
    polynom_order_intertrace = max(1, max(params.get('tmp_polynom_order_intertraces', params['polynom_order_intertraces'])) )
    
    cen_px = np.repeat([spec.shape[1]/2.], len(orders))       # cen_px needs to be float, cen_px defines the zeropoint of np.polyfit
    cen_px_old = copy.copy(cen_px)
    deleted_lines = ''
    for dummy in range(20):
        cen_pxs = arc_lines_wavelength[:,0] * np.nan
        for i,order in enumerate(orders):
            cen_pxs[ (order == arc_lines_wavelength[:,0]) ] = cen_px[i]
        for i in range(arc_lines_wavelength.shape[0] - polynom_order_trace*polynom_order_intertrace):
            x = arc_lines_wavelength[:,2]-cen_pxs
            y = 1.0/(arc_lines_wavelength[:,0]+order_offset)
            weight = arc_lines_wavelength[:,0]*0+1                                                  # No weights
            poly2d_params = polynomial_fit_2d_norm(x, y, arc_lines_wavelength[:,3], polynom_order_trace, polynom_order_intertrace, w=weight)
            arc_line_wave_fit = polynomial_value_2d(x, y, polynom_order_trace, polynom_order_intertrace, poly2d_params)
            arc_line_res = np.abs(arc_line_wave_fit - arc_lines_wavelength[:,3])
            if max(arc_line_res) < 1:       # only good data left, probably should be less than 1 Angstrom
                break
            deleted_lines += '{0}, '.format( arc_lines_wavelength[arc_line_res.argmax(),0:4] )
            arc_lines_wavelength = np.delete(arc_lines_wavelength, arc_line_res.argmax(), 0)        # Delete the least matching line
            cen_pxs = np.delete(cen_pxs, arc_line_res.argmax())                                     # Apply the same correction
        if arc_lines_wavelength.shape[0] <= polynom_order_trace*polynom_order_intertrace:
            logger(('Error: The solution seems unphysical, at least {0} lines should be used (degrees of freedom). '+\
                    'Please add more lines to {1}, or decrease polynom_order_traces or polynom_order_intertraces.').format(polynom_order_trace*polynom_order_intertrace+1, filename))
        
        # Find the new order_offset
        order_offset_old = int(order_offset)
        cenwave = polynomial_value_2d(orders*0, 1.0/(orders+order_offset_old), polynom_order_trace, polynom_order_intertrace, poly2d_params)    #old order_offset at cen_px
        order_offset = find_order_offset(orders, cenwave)
        if order_offset_old != order_offset:            # 
            text1 = 'The new real order offset will be used'
            if np.min(np.abs( (orders+order_offset) )) < 10:             # (orders+order_offset) should never be 0
                text1 = 'The new order offset seems unphysical and will not be used.'
            logger('Warn: The real orders are different from the ones given in {0} (column 2). The old order offset was {2}, the new order offset is {1}. {3}'\
                    .format(filename, order_offset, order_offset_old, text1))
            if np.min(np.abs( (orders+order_offset) )) < 10:             # (orders+order_offset) should never be 0
                order_offset = order_offset_old
        else:
            # The real center from the fit, ignoring the offset
            cen_px, p_cen_px = find_real_center_wavelength_solution(order_offset, orders, cenwave, cen_px, polynom_order_trace, polynom_order_intertrace, poly2d_params)
        if np.sum(np.abs(cen_px - cen_px_old)) < 1. and order_offset_old == order_offset:      # The total improvement of the central pixels is less than 1px
            break
        cen_px_old = copy.copy(cen_px)
    # Redo the fit
    cen_px, p_cen_px = find_real_center_wavelength_solution(order_offset, orders, cenwave, cen_px, polynom_order_trace, polynom_order_intertrace, poly2d_params)
    cen_pxs = arc_lines_wavelength[:,0] * np.nan
    for i,order in enumerate(orders):
        cen_pxs[ (order == arc_lines_wavelength[:,0]) ] = cen_px[i]
    x = arc_lines_wavelength[:,2]-cen_pxs
    y = 1.0/(arc_lines_wavelength[:,0]+order_offset)
    weight = arc_lines_wavelength[:,0]*0+1                                                  # No weights
    poly2d_params = polynomial_fit_2d_norm(x, y, arc_lines_wavelength[:,3], polynom_order_trace, polynom_order_intertrace, w=weight)
    arc_line_wave_fit = polynomial_value_2d(x, y, polynom_order_trace, polynom_order_intertrace, poly2d_params)
    arc_line_res = np.abs(arc_line_wave_fit - arc_lines_wavelength[:,3])
    p_cen_px = np.round(p_cen_px,3)
    text = ('Info: The standard deviation of the residual of the fit to the manual wavelength solution is {1} Angstrom (average is {0} Angstrom). '+\
            'Only input data, for which the residuals were less than 1 Angstrom have been used. '+\
            '{2} of {5} lines have been used to calculate the solution for a 2d polynom fit with {3} orders in dispersion direction (along the traces) '+\
            'and {4} orders in cross-dispersion direction. To fulfil the grating equation the central pixel of the individual orders needs to be '+\
            '{6} + {7}*order + {8}*order**2.').format( round(np.mean(arc_line_res),4), 
                            round(np.std(arc_line_res, ddof=polynom_order_trace+polynom_order_intertrace+1),4), arc_lines_wavelength.shape[0], 
                            polynom_order_trace, polynom_order_intertrace, orig_lines, p_cen_px[2], p_cen_px[1], p_cen_px[0] )
    logger(text)
    calimages['wavelength_solution_result_text'] = text
    if len(deleted_lines) > 0:
        logger('Warn: The following lines from {0} have not been used: {1}'.format(filename, deleted_lines) )
    
    """# Testing
    xarr = np.repeat( [np.arange(max( arc_lines_wavelength[:,2] ))], len(orders), axis=0)
    yarr = []
    for i,order in enumerate(orders):
        wave = polynomial_value_2d(xarr[0]-cen_px[i], 1.0/(order+order_offset), polynom_order_trace, polynom_order_intertrace, poly2d_params)
        p_lin = np.polyfit(xarr[0]-cen_px[i], wave, 1)
        #yarr.append( wave - np.polyval(p_lin, xarr[0]-cen_px[i]) +i*5 )         # Without linear part, offset between orders, can't compare identified lines
        yarr.append( wave )                                                     # Complete solution, to compare with identified lines
    yarr = np.array(yarr)
    xp = list(xarr) + [arc_lines_wavelength[:,2]]
    yp = list(yarr) + [arc_lines_wavelength[:,3]]
    lp = list(np.array(orders,dtype=str)) + ['lines']
    plot_img_spec.plot_points(xp ,yp , lp,'dummy',show=True, title='', x_title='Pixel', y_title='wavelength [Angstrom]')
    good_values = ( np.abs( (xarr[0,:]/50).astype(int) - xarr[0,:].astype(float)/50.) < 1E-5 )      # get every 50th pxel
    l = xarr[0, good_values]                                            # x-values will be the lable
    xarr = np.repeat([orders], len(l), axis=0)                          # orders will be the new xarr
    yarr = yarr[:,good_values].T
    plot_img_spec.plot_points(xarr ,yarr,l,'dummy',show=True, title='', x_title='Orders', y_title='wavelength [Angstrom]')"""
    
    # Transform the wavelength solution into the old wavelength solution
    wavelength_solution, wavelength_solution_arclines = [], []
    xarr = np.arange(0,max(arc_lines_wavelength[:,2]), 0.1)
    max_number_reflines = 0             # will be needed later in order to save the data correctly into a fits file
    for order in orders:
        max_number_reflines = max(max_number_reflines, len(arc_lines_wavelength[(arc_lines_wavelength[:,0]==order),3]))
    for i,order in enumerate(orders):
        yarr = polynomial_value_2d(xarr-cen_px[i], 1.0/(order+order_offset), polynom_order_trace, polynom_order_intertrace, poly2d_params)
        polyfit = np.polyfit(xarr-cen_px[i], yarr, polynom_order_trace)      #lambda from px
        solution = [order+order_offset, cen_px[i] ] + list(polyfit)
        wavelength_solution.append(solution)
        reflines = list(arc_lines_wavelength[(arc_lines_wavelength[:,0]==order),3])     # Wavelength
        for i in range(max_number_reflines - len(reflines)):                            # Fill up with zeros in order to create an array
            reflines.append(0)
        wavelength_solution_arclines.append(reflines)     #The wavelength of the reference lines
        
    wavelength_solution = np.array(wavelength_solution)
    
    return wavelength_solution, np.array(wavelength_solution_arclines)


def adjust_binning_UI(im1, binxy, userinput=True):
    """
    Adjusts the binning
    """
    if not userinput:
        return binxy
    
    #sim_sflat, dummy, dummy = bin_im(im_sflat, params['bin_search_apertures'] )
    # set up plot
    fig, frame = plt.subplots(1, 1)
    fig.set_size_inches(10, 7.5, forward=True)

    # get kwargs
    binx, biny = binxy
    pkwargs = dict(frame=frame, im1=im1, binx=binx, biny=biny)

    # define update plot function
    def plot(frame, im1, binx, biny):
        frame.clear()
        title = ('Adjusting the binning')
        im_bin, dummy, dummy = bin_im(im1, [binx,biny] )
        frame = plot_img_spec.plot_image(im_bin, 'dummy_filename', pctile=0, show=False, adjust=[0.07,0.95,0.95,0.07], title=title, return_frame=True, frame=frame, autotranspose=False, colorbar=False, axis_name=['Cross-dispersion axis [px]','Dispersion axis [px]','flux [ADU]'])
    
    # run initial update plot function
    plot(**pkwargs)

    # define valid_function
    # input is one variable (the string input)
    # return is either:
    #   True and values
    # or
    #   False and error message
    def vfunc_int(xs):
        try:
            value = int(xs)
            return True, value
        except:
            return False, ('Error, input must be integer')

    # define widgets
    widgets = dict()
    widgets['binx'] = dict(label='Binning in\nDispersion axis',
                           kind='TextEntry', minval=None, maxval=None,
                           fmt=str, start=binxy[0], valid_function=vfunc_int,
                           width=10)
    widgets['biny'] = dict(label='Binning in\nCross-dispersion axis',
                           kind='TextEntry', minval=None, maxval=None,
                           fmt=str, start=binxy[1], valid_function=vfunc_int,
                           width=10)
    widgets['accept'] = dict(label='Accept Binning',
                             kind='ExitButton',
                             position=Tk.BOTTOM)
    widgets['update'] = dict(label='Update', kind='UpdatePlot',
                             position=Tk.BOTTOM)
                             
    wprops = dict(orientation='v', position=Tk.RIGHT)

    gui3 = tkc.TkCanvas(figure=fig, ax=frame, func=plot, kwargs=pkwargs,
                        title='HiFlEx: Adjusting the binning', widgets=widgets,
                        widgetprops=wprops)
    gui3.master.mainloop()
    
    #binxy = pkwargs['binxy']       # This doesn't work, pkwargs are only updated within the plot function
    binxy = [ gui3.data['binx'], gui3.data['biny'] ]
    
    plt.close()
    return binxy

def get_leftwidth_rightwidth_tr_poly(centre, left, right):
    w_left  = np.median(centre-left)
    w_right = np.median(right-centre)
    return [w_left, w_right]

def update_tr_poly_width_multiplicate(tr_poly, widths, w_mults, xlows, xhighs):
    """
    Adjusts the polinomial fits for the trace of the order after the extraction width changed.
    Also modifies the widths
    :params w_mults: 1d-list of two floats or 2d-array/list of floats. In the first case gives the multiplicator to left an right borders of the trace
                                    In the second case gives for each order the multiplicator to left an right borders of the trace. Shape needs to be (orders,2)
    """
    w_mults = np.array(w_mults)
    if len(w_mults.shape) == 1:
        w_mults = np.repeat([w_mults],len(tr_poly), axis=0)
    else:
        if w_mults.shape[0] != len(tr_poly) and w_mults.shape[1] == len(tr_poly) and w_mults.shape[0] == 2:       # Transpose if necessay
            w_mults = w_mults.T
        if w_mults.shape[0] != len(tr_poly) or w_mults.shape[1] != 2:
            logger('Error: Expected w_mults.shape = ({0}, {1}) or ({1}, {0}), but got {2}. This is probably a programming error.'.format( len(tr_poly), 2, w_mults.shape ))
    for order in range(len(tr_poly)):
        xarr = np.arange(xlows[order], xhighs[order], 1)
        yarr =  np.polyval(tr_poly[order, 0, 1:], xarr - tr_poly[order, 0, 0])
        yarrl = np.polyval(tr_poly[order, 1, 1:], xarr - tr_poly[order, 1, 0])
        yarrr = np.polyval(tr_poly[order, 2, 1:], xarr - tr_poly[order, 2, 0])
        yarrl, yarrr = adjust_width_orders(yarr, yarrl, yarrr, w_mults[order,:])             # Adjust width
        tr_poly[order,1,1:] = np.polyfit( xarr - tr_poly[order, 1, 0], yarrl, len(tr_poly[order, 1, 1:])-1 )
        tr_poly[order,2,1:] = np.polyfit( xarr - tr_poly[order, 2, 0], yarrr, len(tr_poly[order, 2, 1:])-1 )
        widths[order,0:2] = get_leftwidth_rightwidth_tr_poly(yarr, yarrl, yarrr)
    
    return tr_poly, widths

def remove_orders(pfits, rm_orders):
    """
    Remove orders by masking them
    :param pfits:
    :param rm_orders: list of integers, contains the orders to be removed
    return mask: array of bool with the same length as original orders, True for the orders to keep
    """
    mask = np.repeat([True], len(pfits))
    if type(rm_orders) is not list:
        return mask
    # Remove orders (if rm_orders is populated)
    for r in range(len(pfits)):
        if r in rm_orders:
            mask[r] = False
    return mask

def remove_adjust_orders_UI(im1, pfits, xlows, xhighs, widths=[], shift=0, userinput=True, do_rm=False, do_adj=False, do_shft=False):
    """
    Removes orders from the pfits array in a GUI, allows to adjust the width of th extracted area
    return fmask: array of bool with the same length as original orders, True for the orders to keep
    return pfits: same format as pfits, with adjusted parameters for the polynomial for the left and right border
    return widths: same format as widths, with adjusted left and right stop of the trace
    """
    if not userinput or (not do_rm and not do_adj and not do_shft):
        return remove_orders(pfits, []), pfits, widths         # No orders to remove, pfits stays the same, widths stays the same

    # convert to numpy arrays
    pfits = np.array(pfits)
    xlows, xhighs = np.array(xlows), np.array(xhighs)

    # set up plot
    fig, frame = plt.subplots(1, 1)
    fig.set_size_inches(10, 7.5, forward=True)
    rm_orders = []
    w_mult = 1
    # get kwargs
    pkwargs = dict(frame=frame, im1=im1, pfits=pfits, xlows=xlows, xhighs=xhighs, widths=widths,
                   w_mult=w_mult, rm_orders=rm_orders, shift=shift)

    # define update plot function
    def plot(frame, im1, pfits, xlows, xhighs, widths, w_mult, rm_orders, shift):
        frame.clear()
        title = ''
        if do_adj:
            title += 'Defining the width of the traces.\n'
        if do_shft:
            title += 'Finding the shift of the traces.\n'
        if do_rm:
            title += 'Removing bad orders (Largest order number = {0})\n'.format(len(pfits)-1)
        mask = remove_orders(pfits, rm_orders)
        pfits_shift = copy.deepcopy(pfits)
        if len(pfits_shift.shape) == 3:
            pfits_shift[:,:,-1] += shift        # shift all traces
        else:
            pfits_shift[:,-1] += shift          # shift all traces
        frame = plot_traces_over_image(im1, 'dummy_filename', pfits_shift, xlows, xhighs, widths, w_mult=w_mult, mask=mask, frame=frame, return_frame=True)
        frame.set_title(title[:-1])

    # run initial update plot function
    plot(**pkwargs)

    # define valid_function
    # input is one variable (the string input)
    # return is either:
    #   True and values
    # or
    #   False and error message
    def vfunc(xs):
        try:
            xwhole = xs.replace(',', ' ')
            new_xs = xwhole.split()
            xs = []
            for nxs in new_xs:
                xs.append(int(nxs))
            return True, xs
        except:
            return False, ('Error, input must consist of integers \n '
                           'separated by commas or white spaces')
    def vfunc_float(xs):
        try:
            value = float(xs)
            return True, value
        except:
            return False, ('Error, input must be float')
    # define widgets
    widgets = dict()
    if do_rm:
        widgets['rm_orders'] = dict(label='Select orders to remove',
                                comment='Enter all order numbers to remove \n'
                                        'separated by a whitespace or comma \n'
                                        'to undo just delete the entered '
                                        'number',
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=str, start=" ", valid_function=vfunc,
                                width=40)
    if do_adj:
        widgets['w_mult'] = dict(label='Multiplier for the\nwidth of the traces',
                            comment='If the results are not as wished,\n'
                                    'a modification of the parameter "width_percentile"\n'
                                    'might help. To do this\n'
                                    'Cancel the script with CTRL+C in the terminal\n'
                                    'and then restart',
                            kind='TextEntry', minval=None, maxval=None,
                            fmt=float, start=1.0, valid_function=vfunc_float,
                            width=10)
    if do_shft:
        widgets['shift'] = dict(label='Shift the traces by:',
                            #comment='',
                            kind='TextEntry', minval=None, maxval=None,
                            fmt=float, start=0.0, valid_function=vfunc_float,
                            width=10)           
    widgets['accept'] = dict(label='Accept', kind='ExitButton',
                             position=Tk.BOTTOM)
    widgets['update'] = dict(label='Update', kind='UpdatePlot',
                             position=Tk.BOTTOM)
                             
    wprops = dict(orientation='v', position=Tk.RIGHT)

    gui3 = tkc.TkCanvas(figure=fig, ax=frame, func=plot, kwargs=pkwargs,
                        title='HiFlEx: Locating orders', widgets=widgets,
                        widgetprops=wprops)
    gui3.master.mainloop()

    if 'rm_orders' in gui3.data:
        rm_orders = gui3.data['rm_orders']
    fmask = remove_orders(pfits, rm_orders)
    
    if 'w_mult' in gui3.data and len(pfits.shape) == 3:
        w_mult = gui3.data['w_mult']
        pfits, widths = update_tr_poly_width_multiplicate(pfits, widths, [w_mult, w_mult], xlows, xhighs)
        
    if 'shift' in gui3.data:
        shift = gui3.data['shift']
        if len(pfits.shape) == 3:
            pfits[:,:,-1] += shift        # shift all traces
        else:
            pfits[:,-1] += shift          # shift all traces
    
    plt.close()
    return fmask, pfits, widths

def plot_gauss_data_center(datapoints_x, datapoints, label_datapoins, gauss, label_gauss, centroids, label_centroids, filename='', title='', poly=None, label_poly=None):
    """
    Plot a 1D Gaussian to data
    """
    adjust=[0.05,0.95,0.95,0.05, 1.0,1.01]
    size = datapoints.shape
    for step in range(int((size[0]+9)/10)):          # Only plot 10 graphs in one image
        datarange = range(step*10,min((step+1)*10,size[0]))
        fig, frame = plt.subplots(1, 1)
        fig.set_size_inches(8.3, 5.8)   #A5 landscape
        axes = plt.gca()
        plot_img_spec.plot_points(datapoints_x[datarange,:], datapoints[datarange,:], np.array(label_datapoins)[datarange], '', show=False, adjust=adjust, title='', return_frame=True, frame=frame, x_title='Pixel', y_title='Flux')
        ymin1, ymax1 = axes.get_ylim()
        plt.gca().set_prop_cycle(None)          #Reset the color cycle
        gauss_x, gauss_y, poly_x, poly_y = [], [], [], []
        for i in datarange:
            x = np.linspace(min(datapoints_x[i,:]), max(datapoints_x[i,:]), len(datapoints_x[i,:])*20 )
            y = oneD_gauss(x,gauss[i,:])
            gauss_x.append(x)
            gauss_y.append(y)
            if poly is not None:
                mask = ( ( x >= poly[i,0] ) & ( x <= poly[i,1]  ) )                        # Within +- 5* width
                poly_x.append( x[mask] )
                poly_y.append( np.polyval(poly[i,2:],x[mask]) )
        plot_img_spec.plot_points(gauss_x, gauss_y, np.array(label_gauss)[datarange], '', show=False, adjust=adjust, title='', return_frame=True, frame=frame, x_title='Pixel', y_title='Flux', linestyle="-", marker="")
        ymin2, ymax2 = axes.get_ylim()
        xmin, xmax = axes.get_xlim()
        if poly is not None:
            plt.gca().set_prop_cycle(None)          #Reset the color cycle
            label = np.array(label_poly)[datarange]
            label[1:] = ''
            plot_img_spec.plot_points(poly_x, poly_y, label, '', show=False, adjust=adjust, title='', return_frame=True, frame=frame, x_title='Pixel', y_title='Flux', linestyle=":", marker="")
        plt.gca().set_prop_cycle(None)          #Reset the color cycle
        ymin, ymax = min(ymin1,ymin2), max(ymax1,ymax2)
        centr_x = np.repeat([centroids[datarange]],2,axis=0).T
        centr_y = np.repeat([[ymin, ymax]], len(datarange), axis=0)
        label = np.array(label_centroids)[datarange]
        label[1:] = ''
        plot_img_spec.plot_points(centr_x, centr_y, label, '', show=False, adjust=adjust, title='', return_frame=True, frame=frame, x_title='Pixel', y_title='Flux', linestyle="--", marker="")
        axes.set_ylim(ymin,ymax)
        axes.set_xlim(xmin,xmax)
        frame.set_title(title, fontsize=16)
        if filename == '':
            plt.show()
        else:
            plt.savefig(filename.replace('.png','_{0}-{1}.png'.format('%3.3i'%datarange[0],'%3.3i'%datarange[-1])), bbox_inches='tight')
            plt.close()

def add_specinfo_head(spectra, s_spec, n_spec, w_spec, im_head):
    """
    Add information from the spectra to the header
    :param spectra:
    :param s_spec: 2d-array of floats, spectrum with the signal
    :param n_spec: 2d-array of floats, spectrum with the noise
    """
    spectra = np.array(spectra)
    signal = np.nansum(spectra, axis=1)         # Sum of flux in each order
    snr = np.array(s_spec)/np.array(n_spec)     # snr for each pixel
    cenpx = spectra.shape[1]/2.
    px_range = list(range(int(round(0.8*cenpx)),int(round(1.2*cenpx))+1))
    snr1 = np.nanmedian( snr[:,px_range], axis=1)                   # median SNR in the centre of each order: +-10% around 
    medflux = np.nanmedian( s_spec[:,px_range], axis=1)             # median flux in the centre of each order: +-10% around 
    medflux[medflux < 0] = 0
    snr2 = np.sqrt( medflux )                                       # sqrt of flux is also a SNR
    snr = np.nanmin( np.vstack(( snr1, snr2 )), axis=0 )            # SNR can't be better than sqrt of flux
    snr[np.isnan(snr)] = -1
    if w_spec is not None:
        width_med = np.nanmedian(w_spec, axis=1)
        width_std = np.nanstd(w_spec, ddof=1, axis=1)
        width_min = np.nanmin(w_spec, axis=1)
        width_max = np.nanmax(w_spec, axis=1)
        width_sum = np.nansum(w_spec, axis=1)
    im_head['HIERARCH HiFLEx fmin'] = (round_sig(np.nanmin(spectra),3), 'Minimum Flux per pixel')
    im_head['HIERARCH HiFLEx fmax'] = (round_sig(np.nanmax(spectra),5), 'Maximum Flux per pixel')
    im_head['HIERARCH HiFLEx fsum_all'] = (round_sig(np.nansum(spectra, axis=None),5), 'Total flux')
    if w_spec is not None:
        im_head['HIERARCH HiFLEx EXTR_PX'] = ( round(np.nansum(w_spec),1), 'Total Number of pixel extracted' )
    for order in range(spectra.shape[0]):
        im_head['HIERARCH HiFLEx fsum_order{0}'.format('%2.2i'%order)] = ( round_sig(signal[order],5), 'Flux (counts) in order {0}'.format('%2.2i'%order) )
        im_head['HIERARCH HiFLEx SN_order{0}'.format('%2.2i'%order)] = ( round_sig(snr[order],4), 'median SNR (+-10% of central px)' )
        if w_spec is not None:
            im_head['HIERARCH HiFLEx WIDTH_median{0}'.format('%2.2i'%order)] = ( round_sig(width_med[order],3), 'median width, stdev = {0}'.format(round(width_std[order],1)) )
            im_head['HIERARCH HiFLEx WIDTH_minmax{0}'.format('%2.2i'%order)] = ( '{0}, {1}'.format(round_sig(width_min[order],3), round_sig(width_max[order],3)), 'min and max width' )
            im_head['HIERARCH HiFLEx EXTR_PX{0}'.format('%2.2i'%order)] = ( round(width_sum[order],1), 'Total Number of pixel extracted' )
    
    return im_head

def normalise_continuum(spec, wavelengths, nc=8, ll=2., lu=4., frac=0.3, semi_window=10, nc_noise=15):      # without modal noise nc=4,ll=1.,lu=5. might work
    """
    Normalises the spectrum using the area with continuum and excluding the area with lines
    Adapted from the CERES pipeline, workes for HARPS data and modal noise corrected (blaze corrected) EXOhSPEC data
    ! Possible improvement: use many extracted spectra to define the areas with lines better: See 
    :param frac: fraction of the original spectrum (per order) which needs to remain to fit the polynomial
    :return contspec: 2d array of floats, continuum corrected spectrum
    :return sn_cont: 2d array of floats, Noise in the continuum spectrum
    """
    specs = spec.shape
    ccoefs, sn_coefs, cen_waves, offsets = [], [], [], []
    wavelengths = copy.deepcopy(wavelengths)                            # otherwise wavelengths will be replaced globaly by the next line
    wavelengths[ ( np.isnan(spec) ) ] = np.nan                          # wavelength solution can cover area, in which no good flat extracted spectra exist
    data = np.empty((7,specs[1]), dtype=float)                          # Array with the data: 0: index, 1: wavelength, 2: spec, 3: medfilt spec, 4: medfilt+offset, 5: residuals, 6: fit of the residuals
    data[0,:] = np.arange(specs[1])
    for order in range(specs[0]):
        #print "order", order
        data[1:,:] = np.nan                                             # Initialise
        # Array with the subselections of the data: 0: nonnan, 1: nonnan after medfilt, 2: exclude wavelengths, 3, use for continuum flux, 
        #                                           4: no lines for medfilt, 5: no lines for residuals, 6: zeros at this position of the SN-fit removed
        sub = np.zeros((7,specs[1]), dtype=bool)                        
        data[1,:] = wavelengths[order,:]
        data[2,:] = spec[order,:] 
        data[3,:] = scipy.signal.medfilt(spec[order,:], 2*semi_window+1)            # This function replaces nans by  non-nans
        data[5,:] = data[0,:]*np.nan
        sub[0,:] = ( (~np.isnan(data[1,:])) & (~np.isnan(data[2,:])) )              # nonnan values for spectrum
        sub[1,:] = ( (~np.isnan(data[1,:])) & (~np.isnan(data[3,:])) )              # nonnan values for medfiltered spectrum
        sub[2,sub[0,:]] = ( (data[1,sub[0,:]] < 6561.8) | (data[1,sub[0,:]] > 6563.8) )                  # exclude Halpha line
        sub[3,:] = (sub[1,:] & sub[2,:])                                            # data used for fitting the continuum
        sub[4,:] = sub[3,:] * 1                                                     # used to exclude to high residuals
        ori_len = np.sum(sub[3,:])
        if np.sum(ori_len) < 5:                                                     # Not enough useful data, but fill the arrays
            ccoefs.append( np.zeros(nc+1)*np.nan )
            cen_waves.append( np.nan )
            sn_coefs.append( np.zeros(nc+1)*np.nan )
            offsets.append( 0 )
            continue
        offset = -1 * min(0, np.percentile(data[2,sub[0,:]], 1) )                   # correct the data before fitting to avoid negative values for the fit
        data[4,:] = data[3,:]+offset
        x_cen = (np.nanmin(data[1,sub[1,:]]) + np.nanmax(data[1,sub[1,:]]))/2.        # central wavelength
        if np.std(spec[order,:]) < 1E-10:
            ccoefs.append( np.array( [0]*nc + [np.median(spec[order,:])] ) )
            cen_waves.append( x_cen )
            sn_coefs.append( np.zeros(nc+1) )
            offsets.append( 0 )
            continue
        old_good_values = sum(sub[4,:])
        indexes_search_range = np.zeros((specs[1], specs[1] ), dtype=bool)            # Array to store the indexes which are in a given wavelength difference for each wavelength
        for i in data[0,sub[0,:]].astype(int):                                      # Only use values, for which the wavelength is available
            indexes_search_range[i,sub[0,:]] = ( np.abs(data[1,sub[0,:]] - data[1,i]) <= 0.9)    #+- in Angstrom
            #print i, sum(indexes_search_range[i,:]), data[1,sub[0,:]] - data[1,i]
        for step in range(1,int(old_good_values/10)):
            #print order, step, len(data[1,sub[4,:]]), len(data[4,sub[4,:]]), sum(np.isnan(data[1,sub[4,:]])), sum(np.isnan(data[4,sub[4,:]])), data[1,sub[4,:]]-x_cen, data[4,sub[4,:]], min(step,nc)
            # New, puttting sigmaclipping into own thing
            sub[4,:], p = sigmaclip(data[4,:], x=data[1,:], nc=min(step,nc), sub=sub[4,:], ll=ll, lu=lu, repeats=1, exclude_absorption=True, x_cen=x_cen)
            sub[4,:] *= sub[0,:]            # Remove the nan values again
            # data[5] is not needed
            
            # Old, having sigmaclipping here
            """p = polyfit_adjust_order(data[1,sub[4,:]]-x_cen, data[4,sub[4,:]], min(step,nc) )        # Fit the data
            data[5,:] = data[4,:] - np.polyval(p,data[1,:]-x_cen)                   # Residuals for the whole array, >0 for emission lines
            IU = np.where(data[5,sub[4,:]] > 0)[0]                                  # Positive residuals, use only the good data
            dev = np.mean(data[5,sub[4,:]][IU])                                     # deviation using the positive residuals, use only the good data -> without absorption lines
            sub[4,:] = False
            sub[4,sub[3,:]] = ( (data[5,sub[3,:]] < lu*dev) & (data[5,sub[3,:]] > -ll*dev) )  # Sigmacliping, but allow already clipped data back in"""
            for i in data[0, ~sub[4,:]*sub[1,:]].astype(int):                       # Points that will be removed
                #print i
                if np.sum( sub[4,indexes_search_range[i,:]] ) <= 0:                 # Not enough datapoints remain in a given area
                    #print 'changed',i
                    sub[4,i] = True                                                 # Don't create too big gaps as the high order polynomial nc>4 will not behave well there
            if ( sum(sub[4,:]) < frac*ori_len or sum(sub[4,:]) == old_good_values ) and step >= nc:         # stop if too much data is removed or no changes, but not if only linear equation is fitted
                #if order == 29:            # To check results
                #    plot_img_spec.plot_points([data[1,:],data[1,:],data[1,:],data[1,sub[4,:]]], [data[2,:], data[3,:], np.polyval(p,data[1,:]-x_cen), data[3,sub[4,:]]], \
                #                      ['data ori', 'data medfilt', 'fit {0} orders'.format(len(p)-1), 'remaining'], '', show=True, return_frame=False, \
                #                      x_title='wavelength [Angstrom]', y_title='flux [ADU]', title='order {1}, step={0}'.format(step, order))
                break
            old_good_values = sum(sub[4,:])
        ccoefs.append(p)
        cen_waves.append(x_cen)
        #x_full = wavelengths[order,:]                           #only for testing
        #y_full = scipy.signal.medfilt(spec[order,:], window)    #only for testing
        #plot_img_spec.plot_points([x_full,x_full,x_range], [y_full,np.polyval(p,x_full-x_cen),y_range], ['full data', 'fit', 'data for fit'], '', show=True, return_frame=False, x_title='wavelength [Angstrom]', y_title='flux [ADU]')
        
        # Fit the residuals of the individual data points, excluding lines
        sub[5,:] = sub[4,:] * 1
        data[5,:] = data[2,:]+offset - np.polyval(p, data[1,:]-x_cen)          # Residual
        for step in range(5):
            IU = np.where(data[5,sub[5,:]] > 0)[0]                              # Positive residuals
            dev = np.mean(data[5,sub[5,:]][IU])                                 # deviation using the positive residuals, use only the good data -> without absorption lines
            sub[5,:] = False
            sub[5,sub[3,:]] = ( (data[5,sub[3,:]] < lu*dev) & (data[5,sub[3,:]] > -ll*dev) )  # Sigmacliping, but allow already clipped data back in
            for i in data[0, ~sub[5,:]*sub[1,:]].astype(int):                   # Points that will be removed
                if np.sum( sub[5,indexes_search_range[i,:]] ) <= 0:             # Not enough datapoints remain in a given area
                    sub[5,i] = True                                             # Don't create too big gaps as the high order polynomial nc>4 will not behave well there

        p_sn = polyfit_adjust_order(data[1,sub[5,:]]-x_cen, np.abs(data[5,sub[5,:]]), nc_noise)          # Fit the noise using the the residuals
        # The fit of the sn should never create values equal/below 0, the steps below are used to reach this
        data[6,:] = np.polyval(p_sn, data[1,:]-x_cen)
        for nc_step in range(1,nc_noise)[::-1]:                                 # Use a lower order, if fit remains stupid
            #print "nc_step", nc_step
            for dummy in range(100):
                with np.errstate(invalid='ignore'):
                    if sum(~sub[5,:]*sub[1,:]) > 0:                                     # There are values that can be replaced
                        index = data[0,sub[1,:]][np.abs(data[6,sub[1,:]] - np.nanmin(data[6,~sub[5,:]]) ) < 1E-4].astype(int)   # Find the position where the fit has the smallest value and the residual there isn't used for the fit
                    else:
                        if sub[6,:].all():      # ~sub[6,:] are all False -> things will fail
                            continue
                        index = data[0,~sub[6,:]][np.abs(data[6,~sub[6,:]] - np.nanmin(data[6,~sub[6,:]]) ) < 1E-4].astype(int) # Find the smalles value of the fit, ignore already changed values
                #print order,index, np.max(data[5,sub[5,:]]), data[5,index], sub[5,index]
                data[5,index] = np.max(np.abs(data[5,sub[5,:]]))                                            # Replace with the maximum residual as worst case
                sub[5:7,index] = True                                                                       # use the value for the fit
                p_sn = polyfit_adjust_order(data[1,sub[5,:]]-x_cen, np.abs(data[5,sub[5,:]]), nc_step)        # Fit the noise using the the residuals
                data[6,:] = np.polyval(p_sn, data[1,:]-x_cen)
                if np.nanmin( data[6,sub[1,:]] ) > 0:                                                       # End the loop
                    break
            if np.nanmin( data[6,sub[1,:]] ) > 0:                                                       # End the loop
                break
        sn_coefs.append(p_sn)
        offsets.append(offset)
        if order == 129:
            plot_img_spec.plot_points([data[1,:], data[1,:], data[1,sub[5,:]], data[1,:], data[1,:], data[1,sub[4,:]] ], [data[2,:], np.abs(data[5,:]), np.abs(data[5,sub[5,:]]), \
                                        np.polyval(ccoefs[-1],data[1,:]-cen_waves[-1]), data[6,:], data[3,sub[4,:]] ], \
                                  ['data', 'abs res', 'abs res used', 'fit flux', 'fit_sn', 'medfilt used'], '', show=True, return_frame=False, \
                                  x_title='wavelength [Angstrom]', y_title='flux [ADU]', title='order {0}'.format(order))
    # Transform the coefficents into data array
    contspec, sn_cont = [], []
    showorder = 129
    for order in range(specs[0]):
        fit_flux = np.polyval(ccoefs[order],wavelengths[order,:]-cen_waves[order])
        fit_sn = np.polyval(sn_coefs[order],wavelengths[order,:]-cen_waves[order])
        for i in range(len(fit_flux)-1):
            if fit_flux[i]*fit_flux[i+1] < 0:                           # change of sign
                fit_flux[max(0,i-semi_window):min(len(fit_flux),i+semi_window+1)] = np.nan     # ignore this data as artificial peaks are introduced due to the division close to 0 in the flat spectrum
        with np.errstate(invalid='ignore'):
            bad_value = ( (fit_sn <= 0) | (fit_flux <= 0) )                                  # Errors smaller than 0, or Flux smaller than 0 (will make absorption lines into emission lines)
        if np.sum(bad_value) > 5:
            logger(('Warn: Normalisation of order {0} had a problem. {1} pixel of the fitted flux are below zero and {2} pixel of the fitted signal-to-noise is below zero. ' +\
                   'Values below 0 should not hapen, use a polynomial with lower order by adjusting measure_noise_orders').format(order, 
                                                    np.sum(fit_flux[~np.isnan(fit_flux)] <= 0), np.sum(fit_sn[~np.isnan(fit_sn)] <= 0) ))
            showorder = order
        fit_sn[bad_value] = np.nan
        fit_flux[bad_value] = np.nan
        contspec.append( (spec[order]+offsets[order]+0.0) / fit_flux)
        sn_cont.append(fit_sn)
        #if order == showorder:
        #    plot_img_spec.plot_points([wavelengths[order,:],wavelengths[order,:],wavelengths[order,:],wavelengths[order,:]], [spec[order,:], fit_flux, fit_sn, contspec[-1]], \
        #                          ['data', 'fit flux', 'fit_sn', 'continuum'], '', show=True, return_frame=False, \
        #                          x_title='wavelength [Angstrom]', y_title='flux [ADU]', title='order {0}, offset {1}'.format(order, offsets[order]))
    
    """ doesn't improve things in the overlapping area
    spec = np.array(contspec)
    ccoefs = get_cont_ceres(wavelengths ,spec)      # doesn't really work for solar spectra without a correct flat correction    
    contspec = []
    for order in range(specs[0]):
        contspec.append(spec[order] / np.polyval(ccoefs[order],wavelengths[order,:]))
    """
    return np.array(contspec), np.array(sn_cont)

def linearise_wavelength_spec(params, wavelengths, spectra, method='sum', weight=[]):    
    """
    :param spectra: list of 2d array of floats, contains the spectrum to be linearised
    """
    specs = spectra.shape
    dwave = params['wavelength_scale_resolution']
    px_range = np.arange(specs[1])
    data = np.array([]).reshape(0,3)
    for order in tqdm(range(specs[0]), desc='Info: Linearise the spectrum'):
        wave_range = wavelengths[order,:]
        wave_diff = np.nanmax(np.abs(wave_range[1:] - wave_range[:-1]))        # Difference in wavelength between 2 data points, abs needed for dodgy wavelength solutions
        start = int(min(wave_range)/dwave) * dwave          # Minumum wavelength should a multiple of the wavelength_scale_resolution
        lwave_range = np.arange(start, max(wave_range)+dwave, dwave )
        """ # This solution takes long and has the disadvantage described in the comments
        spec = spectra[order,:]
        nonnan = ~np.isnan(spec)
        # spec can be replaced by Spline fitted to it
        # https://stackoverflow.com/questions/17913330/fitting-data-using-univariatespline-in-scipy-python
        # https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.interpolate.UnivariateSpline.integral.html
        for s in [0.1]:
            start = time.time()
            # !!! Spline could create wrong values: "A theoretically impossible result was found during the iteration"
            # s=0.1: relative difference between individual points is less than 1e-6 for just-extracted data.
            # For continuum corrected data it's 0.5%!!! -> s should be relative to the data range (excluding outliers)
            f_spec = inter.UnivariateSpline (wave_range[nonnan], spec[nonnan], s=s)       
            #print time.time()-start, np.nansum(np.abs(f_spec(wave_range)-spec))
        lspec = f_spec(lwave_range)
        for wave in wave_range[~nonnan]:        # Ignore values which were not defined
            diff = np.abs(lwave_range - wave)   # Wavelength difference to the non defined values
            lspec[diff < wave_diff] = np.nan
        #plot_img_spec.plot_points([wave_range, wave_range, lwave_range], [spec, f_spec(wave_range), lspec], ['orig', 'spline', 'lin'], '', show=True, return_frame=False, x_title='wavelength [Angstrom]', y_title='flux [ADU]')
        lweight = np.ones(lwave_range.shape)
        if len(weight) > 0 :
            if weight.shape != specs:
                logger('Error: The weights used for overlapping areas when linearising the wavelength solution have a different shape than the spectrum. This is most likely a coding error')
            for i,wave in enumerate(lwave_range):
                index = np.argmin(np.abs(wave_range - wave))    # find the position of the closest original wavelength 
                lweight[i] = weight[order,index]
        """
        # Short and simple solution, but the minimum and maximum points of curves get lost, if the original data point ends up between 2 linear wavelengths
        # Using a interpolation can create much bigger maxima or minima
        data_o = np.ones(( lwave_range.shape[0], 3 ))
        data_o[:,0] = lwave_range
        data_o[:,1] = np.nan
        for i,wave in enumerate(lwave_range):
            diff = np.abs(wave_range - wave)                # Wavelength difference to the original wavelengths
            indexes = np.where( diff < wave_diff )[0]
            #print( order,i, diff.shape, indexes.shape, diff, np.min(diff), wave_diff, wave, spectra[order,indexes].shape, diff[indexes] )
            if indexes.shape[0] == 0:
                continue
            weights = 1./(diff[indexes]+0.01*dwave)
            data_o[i,1] = np.average( spectra[order,indexes], weights=weights )
            data_o[i,2] = np.average( weight[order,indexes], weights=weights )
        
        data_o = data_o[ ~np.isnan(data_o[:,1]),:]
        data = np.vstack([data, data_o ])       # wavelengths in first column, flux in second, weights in third
    # Combine the overlapping orders
    clr_dat = np.round(data[:,0]/dwave).astype(int)                             # make into integer, as both 7391.58 and 7391.58000001 can be in the wavelength column
    uniq_wave, number = np.unique(clr_dat,return_counts=True)
    todel = np.array([])
    nolog = False
    #print 'data.shape, uniq_wave.shape', data.shape, uniq_wave.shape
    for wave in uniq_wave[number > 1]:      # check the wavelengths with a count > 1
        indexes = np.where(clr_dat == wave)[0]
        todel = np.append(todel, indexes[1:], axis=None)
        if method == 'sum':
            data[indexes[0],1] = np.nansum(data[indexes,1])
        elif method == 'mean':
            data[indexes[0],1] = np.nanmean(data[indexes,1])
        elif method == 'weight':
            notnan = ~np.isnan(data[indexes,1])
            if notnan.any() > 0:
                #if wave*dwave>6343 and wave*dwave<6348:
                #    print( data[indexes,:], data[indexes,1][notnan], data[indexes,2][notnan], sum(data[indexes,2][notnan]) )
                data[indexes[0],1] = np.average(data[indexes,1][notnan], weights=data[indexes,2][notnan])
        else:
            if not nolog:
                logger('Warn: method {0} is not known, using sum'.format(method))
                nolog = True
            data[indexes[0],1] = np.nansum(data[indexes,1])
    data = np.delete(data, todel.astype(int), axis=0)
    #print 'data.shape, len(todel)',data.shape, len(todel)
    # Add wavelengths in case there are gaps between orders
    lwave_range = np.arange(np.nanmin(data[:,0]), np.nanmax(data[:,0])+dwave/2., dwave )
    #missing_wave = np.setxor1d(lwave_range, data[:,0])         # Doesn't work because of 7391.58 and 7391.58000001 -> used the the following 3 lines
    clr_dat = np.round(data[:,0]/dwave).astype(int)             # make into integer, as both 7391.58 and 7391.58000001 can be in the wavelength column
    lwave_range = np.round(lwave_range/dwave).astype(int)       # make into integer, as both 7391.58 and 7391.58000001 can be in the wavelength column
    missing_wave = np.setxor1d(lwave_range, clr_dat)*dwave
    temp = np.empty([2, len(missing_wave)])
    temp.fill(np.nan)
    missing_wave = np.vstack([ missing_wave, temp ]).T
    data = np.vstack([data, missing_wave])
    data = data[data[:,0].argsort()]            # sort by wavelength
    #print 'data.shape, lwave_range.shape, missing_wave.shape', data.shape, lwave_range.shape, missing_wave.shape, np.nanmin(data, axis=0), np.nanmax(data, axis=0)
    #plot_img_spec.plot_points([data[:,0],data[:,0]], [data[:,1], data[:,2]], ['flux', 'weight'], '', show=True, return_frame=False, x_title='wavelength [Angstrom]', y_title='flux [ADU]')
    
    return data[:,0], data[:,1]

"""def mjd_fromheader(params, head):                   # Not necessary anymore
    ""#"
    :return mjd: modified Julian date from header, in UT of the mid of the exposure
    :return mjd0: 2400000.5
    ""#"
    secinday = 24*3600.0
    
    obsdate, obsdate_float, exposure_time, obsdate_begin, exposure_fraction, jd_midexp = get_obsdate(params, head)               # obsdate, obsdate_float, in UTC, mid of the exposure
    # Mid exppsure
    mjd0,mjd = iau_cal2jd(obsdate.year, obsdate.month, obsdate.day)
    ut       = obsdate.hour*3600. + obsdate.minute*60. + obsdate.second
    mjd += ut / secinday                    # Mid-exposure MJD
    # Start exposure
    mjd0,mjd_begin  = iau_cal2jd(obsdate_begin.year, obsdate_begin.month, obsdate_begin.day)
    ut              = obsdate_begin.hour*3600. + obsdate_begin.minute*60. + obsdate_begin.second
    mjd_begin       += ut / secinday        # Start-exposure MJD
    
    return mjd, mjd0, mjd_begin

def iau_cal2jd(IY,IM,ID):               # from CERES, modified                   # Not necessary anymore
    IYMIN = -4799.
    MTAB = np.array([ 31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])
    if IY < IYMIN:
        logger('Error: The year of the observation is not defined: {0}'.format(IY))
    else:
        if IM>=1 and IM <= 12:
            if IY%4 == 0:
                MTAB[1] = 29.
            else:
                MTAB[1] = 28.
            if IY%100 == 0 and IY%400!=0:
                MTAB[1] = 28.
            if ID < 1 or ID > MTAB[IM-1]:
                J = -3
            a = int( ( 14 - IM ) / 12 )
            y = IY + 4800 - a
            m = IM + 12*a -3
            DJM0 = 2400000.5
            DJM = ID + int((153*m + 2)/5) + 365*y + int(y/4) - int(y/100) + int(y/400) - 32045 - 2400001.
        else:
            logger('Error: The month of the observation is not defined: {0}'.format(IM))
    return DJM0, DJM"""

def getcoords_from_file(obnames, mjd, filen='coords.txt', warn_notfound=True, ignore_values=False):               # from CERES, heavily modified
    """
    1-  name of the target as specified in the image header.
    2-  right ascension of the target (J2000) with format hh:mm:ss or as float in degrees.
    3-  declination of the target (J2000) with format dd:mm:ss or as float in degrees.
    4-  proper motion in RA [mas/yr].
    5-  proper motion in DEC [mas/yr].
    6-  integer (0 or 1). If 0, the code uses the cordinates given in the image header.
        If 1, the code uses the coordinates given in this file.
    7-  mask that will be used to compute the CCF. Allowed entries are G2, K5 and M2.
    8-  velocity width in km/s that is used to broaden the lines of the binary mask.
        It should be similar to the standard deviation of the Gaussian that is fitted to the CCF.
    9-  epoch (J2000)
    
    :params ignore_values: Just checks for the object name, but not any of the other fields, just creates a list of strings
    """
    RA0, DEC0, PMRA, PMDEC, epoch = 0., 0., 0., 0., 2000.
    
    if not os.path.isfile(filen) and warn_notfound:
        logger('Warn: Reference coordinates file {0} does not exist.'.format(filen))
        return RA0, DEC0, epoch, PMRA, PMDEC, obnames
    # !!! Use read files with split already available
    lines_txt = read_text_file(filen, no_empty_lines=True, warn_missing_file=False)
    if ignore_values:
        lines = convert_readfile(lines_txt, [str, str, str, str,   str,   str,   str,   str, str], delimiter=',', expand_input=True)    #replaces=[['\t',',']]
    else:                       # Standard
        lines = convert_readfile(lines_txt, [str, str, str, float, float, float, str, str, float], delimiter=',', expand_input=True, ignore_badlines=True)  #replaces=[['\t',',']]
        if len(lines) < len(lines_txt):
            logger('Warn: {1} line(s) could not be read in the reference coordinates file: {0}. Please check that columns 4 to 6 and 9 (starting counting with 1) are numbers'.format(filen, len(lines_txt)-len(lines) ))
    found = False
    for cos in lines:
        if not ignore_values:
            if abs(cos[5]) < 0.9:         # disabled
                continue
        for obname in obnames:
            if cos[0].lower() == obname.lower() or cos[0].lower().replace('_','') == obname.lower() or cos[0].lower().replace('-','') == obname.lower() or cos[0].lower().replace(' ','') == obname.lower():
                if ignore_values:
                    RA0, DEC0, PMRA, PMDEC, epoch = cos[1], cos[2], cos[3], cos[4], cos[8]
                else:
                    if cos[1].find(':') > 0 or cos[1].find(' ') > 0:
                        cos[1] = cos[1].replace(' ',':')
                        cos1 = cos[1].split(':')
                        try:
                            RA0 = (float(cos1[0]) + float(cos1[1])/60. + float(cos1[2])/3600.) * 360. / 24.
                        except:
                            logger('Warn: Problem with the right ascension of entry {1} in the reference coordinates file {0}.'.format(filen,cos))
                            break
                    else:           # Already stored as float (asuming degree)
                        try:
                            RA0 = float(cos[1])
                        except:
                            logger('Warn: Problem with the right ascension of entry {1} in the reference coordinates file {0}.'.format(filen,cos))
                            break
                    if cos[2].find(':') > 0 or cos[2].find(' ') > 0:
                        cos[2] = cos[2].replace(' ',':')
                        cos2 = cos[2].split(':')
                        try:
                            DEC0 = np.absolute(float(cos2[0])) + float(cos2[1])/60. + float(cos2[2])/3600.
                        except:
                            logger('Warn: Problem with the declination of entry {1} in the reference coordinates file {0}.'.format(filen,cos))
                            break
                        if cos2[0][0] == '-':       # negative Declination
                            DEC0 = -DEC0
                    else:
                        try:
                            DEC0 = float(cos[2])
                        except:
                            logger('Warn: Problem with the declination of entry {1} in the reference coordinates file {0}.'.format(filen,cos))
                            break
                    PMRA = float(cos[3])    # mas/yr
                    PMDEC = float(cos[4])   # mas/yr
                    try:
                        epoch = float(cos[8])
                    except:
                        logger('Warn: Problem with the epoch of entry {1} in the reference coordinates file {0}.'.format(filen,cos))
                    """ # Steps to adjust the coordinates to the correct position, not necessary anymore because not using CERES routines to get BCVel and the barycorrpy package takes care of it
                    mjdepoch = 2451545.0 - constants.MJD0 + (float(cos[5]) - 2000.)
        
                    RA  = RA0 + (PMRA/1000./3600.)*(mjd-mjdepoch)/365.
                    DEC = DEC0 + (PMDEC/1000./3600.)*(mjd-mjdepoch)/365."""
                obnames = [obname]
                found = True
                break
        if found:
            break       # otherwise cos will change
    if not found:
        cos = []
        if warn_notfound:
            logger('Warn: Object was not found in the reference coordinates file {0}.'.format(filen))
    return RA0, DEC0, epoch, PMRA, PMDEC, obnames, cos

def JPLR0(lat, altitude):               # from CERES, not needed anymore
    "the corrections due to earth rotation"
    "this function returns the velocity in m/s, the projected distance of the observatory in the equator axis and in the earth spin axis and \
    also returns the distance from the rotation axis, and the distance from the equator plane to the observatory position"
    "The arguments are latitude,altitude and hour angle at observatory , dec: the declination of the star\
    the variables are: \
    Raxis: the distance from the observatory to the earth axis, and Psid: the period of sidereal day"
    lat = Constants.degtorad*lat
    e2 = Constants.f*(2.0-Constants.f)
    c1 = 1.0 - e2*(2.0 - e2)*np.sin(lat)**2
    c2 = 1.0 - e2*np.sin(lat)**2

    #radius at 0 elevation    
    R0 = Constants.Req*np.sqrt(c1/c2) 
    
    #the geocentric latitude 
    c1 = e2*np.sin(2.0*lat)
    c2 = 2.0*c2
    geolat = lat - np.arctan(c1/c2)     
    #Calculate geocentric radius at altitude of the observatory 
    GeoR = R0*np.cos(geolat) + altitude*np.cos(lat)
    
    # the R0 vector is now the distance from the observatory to the declination 0 deg plane
    R0 = R0*np.sin(abs(geolat))+altitude*np.sin(lat)
    return GeoR,R0

"""def obspos(longitude,obsradius,R0):               # from CERES, not needed anymore
    "#""
        Set the observatory position respect to geocenter in the coordinates(x,y,z)required by jplepem, 
        x to equator/greenwich intersection, 
        y 90 degrees east, 
        z positive to north
        "#""
    obpos = []    
    x = obsradius*np.cos( (np.pi / 180.0) * longitude )
    obpos.append(x)    
    y = obsradius*np.sin( (np.pi / 180.0) * longitude )
    obpos.append(y)
    z = R0
    obpos.append(z)
    return obpos"""

def get_object_site_from_header(params, im_head, obnames, obsdate):
    # Define the header keywords, if available
    site_keys       = ['TELESCOP']
    altitude_keys   = ['HIERARCH ESO TEL GEOELEV', 'ESO TEL GEOELEV']       # HIERARCH will be removed from header keywords in python
    latitude_keys   = ['HIERARCH ESO TEL GEOLAT', 'ESO TEL GEOLAT']
    longitude_keys  = ['HIERARCH ESO TEL GEOLON', 'ESO TEL GEOLON']
    ra_keys         = ['RA', 'RA-DEG']
    dec_keys        = ['DEC', 'DEC-DEG']
    epoch_keys      = ['HIERARCH ESO TEL TARG EQUINOX', 'ESO TEL TARG EQUINOX']
    pmra_keys       = ['HIERARCH ESO TEL TARG PMA', 'ESO TEL TARG PMA']
    pmdec_keys      = ['HIERARCH ESO TEL TARG PMD', 'ESO TEL TARG PMD']
    params['epoch'] = 2000.0
    params['pmra']  = 0.
    params['pmdec'] = 0.
    if 'altitude' not in params.keys():
        params['altitude'] = 0.
    source_obs = 'The site coordinates from the configuration file were used.'
    settings = []
    settings.append( [0, site_keys, 'site'] )
    settings.append( [0, altitude_keys, 'altitude'] )
    settings.append( [0, latitude_keys, 'latitude'] )
    settings.append( [0, longitude_keys, 'longitude'] )     # Order is important, at this stage 'site', 'altitude', and 'latitude' need to be defined, but 'ra' and 'dec' will be overwritten later
    settings.append( [0, ra_keys,    'ra'] )
    settings.append( [0, dec_keys,   'dec'] )
    settings.append( [0, epoch_keys, 'epoch'] )
    settings.append( [0, pmra_keys,  'pmra'] )
    settings.append( [0, pmdec_keys, 'pmdec'] )
    if type(obnames).__name__ not in ['list', 'ndarray']:
        obnames = [obnames]
    for [i, header_key_words, parentr] in settings:
        if parentr == 'longitude':                  # Enough information to calculate the ephemerides of sun and moon.
            gobs = ephem.Observer()  
            gobs.name = copy.copy(params['site'])
            gobs.lat  = math.radians(params['latitude'])     # lat/long in decimal degrees  
            gobs.long = math.radians(params['longitude'])
            gobs.elevation = params['altitude']
            gobs.date = obsdate.strftime('%Y-%m-%d %H:%M:%S')
            mephem    = ephem.Moon()
            mephem.compute(gobs)
            sephem    = ephem.Sun()
            sephem.compute(gobs)
            jephem    = ephem.Jupiter()
            jephem.compute(gobs)
            params['ra'] = -999                                             # If not changed this means the object coordinates were made up
            for obname in obnames:
                if obname.lower().find('sun') == 0:     # check with file_assignment.py: calibration_parameters_coordinates_UI(), creating the widgets
                    params['ra']          = sephem.ra
                    params['dec']         = sephem.dec
                    source_radec = 'The object coordinates are derived from the calculated solar ephermeris.'
                    obnames.insert(0,'Sun')
                    break
                elif obname.lower().find('moon') == 0:
                    params['ra']          = mephem.ra
                    params['dec']         = mephem.dec
                    source_radec = 'The object coordinates are derived from the calculated lunar ephermeris.'
                    obnames.insert(0,'Moon')
                    break
                elif obname.lower().find('jupiter') == 0:
                    params['ra']          = jephem.ra
                    params['dec']         = jephem.dec
                    source_radec = 'The object coordinates are derived from the calculated Jupiters ephermeris.'
                    obnames.insert(0,'Jupiter')
                    break
            if params['ra'] == -999:
                params['ra']          = mephem.ra       # To fill in the parameter, might be overwritten later by [0, ra_keys, 'ra']
                params['dec']         = mephem.dec      # To fill in the parameter, might be overwritten later by [0, dec_keys, 'dec']
                source_radec = 'Warn: The object coordinates were made up!'
        # Overwrite hard coded values with the information from the header
        for entry in header_key_words:                  # HIERARCH will be removed from header keywords in python
            if entry.replace('HIERARCH ','') not in header_key_words:
                header_key_words.append(entry.replace('HIERARCH ',''))          # header keyword without HIERARCH
        for entry in header_key_words:
            if entry in im_head.keys():
                if parentr not in params.keys():   # only, if the header_key is not empty (or if the the key is missing in params)
                    params[parentr] = im_head[entry]    # Get the information from the header, overrides the manual values
                elif type(im_head[entry]) == str:           # Can't replace a float
                    if im_head[entry].replace(' ','') != "":
                        params[parentr] = im_head[entry]    # Get the information from the header, overrides the manual values
                else:
                    params[parentr] = im_head[entry]    # Get the information from the header, overrides the manual values
                if parentr == 'ra':                 # Assume that dec is coming from the same source
                    source_radec = 'The object coordinates are derived from the image header.'
                if parentr == 'latitude':           # Assume that latitute is coming from the same source
                    source_obs = 'The site coordinates are derived from the image header.'
    return params, source_radec, source_obs, mephem

def get_barycent_cor(params, im_head, obnames, ra2, dec2, epoch, pmra, pmdec, obsdate, leap_update=False):
    """
    Calculates the barycentric correction.
    To do this, position of the telescope and pointing of the telescope need to be known.
    Header vaulues are read, if they are not available, then using 
    Set leap_update to true, once the servers are available again, see https://github.com/shbhuk/barycorrpy/issues/27
    """
    
    params, source_radec, source_obs, mephem = get_object_site_from_header(params, im_head, obnames, obsdate)
    # mjd,mjd0, mjd_begin = mjd_fromheader(params, im_head)
    mjd = im_head['HIERARCH HiFLEx MJD']
    jd  = im_head['HIERARCH HiFLEx JD']
    
    #ra2, dec2, epoch, pmra, pmdec, obnames, dummy = getcoords_from_file(obnames, mjd, filen=reffile, warn_notfound=False)     # obnames will be a list with only one entry: the matching entry 
    if ra2 !=0 and dec2 != 0:
        params['ra']  = ra2
        params['dec'] = dec2
        params['epoch'] = epoch
        params['pmra'] = pmra
        params['pmdec'] = pmdec
        source_radec = 'The object coordinates are derived from the reference file.'
    if type(params['ra']) == ephem.Angle:
        params['ra'] = math.degrees(params['ra'])           # convert radians(float) into a float
    if type(params['dec']) == ephem.Angle:
        params['dec'] = math.degrees(params['dec'])         # convert into a float
    ra, dec = params['ra'], params['dec']

    bcvel_baryc, bjd = 0.0, 0.0
    if source_radec != 'Warn: The object coordinates were made up!':
        # Using barycorrpy (https://github.com/shbhuk/barycorrpy), pip install barycorrpy
        site = ''
        if params['site'] in barycorrpy.EarthLocation.get_site_names():
            site = params['site']
        ephemeris2 = 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp'
        ephemeris2 = 'de430'        # See https://ssd.jpl.nasa.gov/?planet_eph_export
        # Calculate the barycentric corrections for 60s intervals
        jd_start = im_head['HIERARCH HiFLEx JD_START']
        exposure = im_head['HIERARCH HiFLEx EXPOSURE']
        jd_range = [jd] + list(np.arange(jd_start, jd_start+exposure/(3600.*24), 60/(3600.*24))) + [jd_start+exposure/(3600.*24)]   # Every minute
        success = False
        for ii in range(5):
            if sys.version_info[0] < 3:     # Python 2
                try:
                    bcvel_baryc_range = barycorrpy.get_BC_vel(JDUTC=jd_range,ra=ra,dec=dec,obsname=site,lat=params['latitude'],longi=params['longitude'],alt=params['altitude'],
                                                      pmra=params['pmra'],pmdec=params['pmdec'],px=0,rv=0.0,zmeas=0.0,epoch=params['epoch'],
                                                      ephemeris=ephemeris2,leap_dir=params['logging_path'], leap_update=leap_update)
                    success = True
                    break
                except (urllib2.URLError, ValueError, astropy.utils.iers.iers.IERSRangeError) as e:
                    try:
                        logger('Warn: Problem downloading file for barycentric correction. Will try {0} more times. Error: {1}, Reason: {2}'.format(4-ii, e, e.reason))
                    except:
                        logger('Warn: Problem downloading file for barycentric correction. Will try {0} more times. Error: {1}'.format(4-ii, e))
            else:
                try:
                    bcvel_baryc_range = barycorrpy.get_BC_vel(JDUTC=jd_range,ra=ra,dec=dec,obsname=site,lat=params['latitude'],longi=params['longitude'],alt=params['altitude'],
                                                      pmra=params['pmra'],pmdec=params['pmdec'],px=0,rv=0.0,zmeas=0.0,epoch=params['epoch'],
                                                      ephemeris=ephemeris2,leap_dir=params['logging_path'], leap_update=leap_update)
                    success = True
                    break
                except (urllib.error.URLError, ValueError, astropy.utils.iers.iers.IERSRangeError) as e:
                    logger('Warn: Problem downloading file for barycentric correction. Will try {0} more times. Error: {1}, Reason: {2}'.format(4-ii, e, e.reason))

        if not success:
            logger('Error: Barycentric velocities could not be calculated.')
        bcvel_baryc_range = bcvel_baryc_range[0] / 1E3       # in km/s
        bcvel_baryc = bcvel_baryc_range[0]
        
        im_head['HIERARCH HiFLEx RA'] = (round(ra,6),           'RA in degrees, used to calculated BCV, BJD')
        im_head['HIERARCH HiFLEx DEC'] = (round(dec,6),         'DEC in degrees, used to calculated BCV, BJD')
        im_head['HIERARCH HiFLEx PMRA'] = (round(pmra,3),       'proper motion for RA in mas/yr, for BCV, BJD')
        im_head['HIERARCH HiFLEx PMDEC'] = (round(pmdec,3),     'proper motion for DEC in mas/yr, for BCV, BJD')
        im_head['HIERARCH HiFLEx BCV'] = (round(bcvel_baryc,4), 'Barycentric velocity in km/s')
        im_head['HIERARCH HiFLEx BCV MAX'] = (round(max(bcvel_baryc_range),4), 'Maximum BCV')
        im_head['HIERARCH HiFLEx BCV MIN'] = (round(min(bcvel_baryc_range),4), 'Minimum BCV')
        
        """# Using jplephem
        bjd = jd_corr(mjd, ra, dec, params['epoch'], params['latitude'], params['longitude'], jd_type='bjd')
        bjd = bjd[0]"""
        
        # Using barycorrpy (https://github.com/shbhuk/barycorrpy), pip install barycorrpy
        bjdtdb = barycorrpy.utc_tdb.JDUTC_to_BJDTDB(JDUTC=jd,ra=ra,dec=dec,obsname=site,lat=params['latitude'],longi=params['longitude'],alt=params['altitude'],pmra=params['pmra'],
                                        pmdec=params['pmdec'],px=0,rv=0.0,epoch=params['epoch'],ephemeris=ephemeris2,leap_update=leap_update)           # only precise to 0.2s 
        bjd = bjdtdb[0][0]

        im_head['HIERARCH HiFLEx BJDTDB'] = (round(bjd,6), 'Baryc. cor. JD (incl leap seconds)')     # without leap seconds: remove 32.184+N leap seconds after 1961'

    logger(('Info: Using the following data for object name(s) {10}, Observatory site {9}, mid exposure MJD {11}: '+\
                    'altitude = {0}, latitude = {1}, longitude = {2}, ra = {3}, dec = {4}, epoch = {5}, pmra = {13}, pmdec = {14}. {8} {6} '+\
                    'This leads to a barycentric velocity of {7} km/s and a mid-exposure BJD-TDB of {12}').format(params['altitude'], params['latitude'], params['longitude'], 
                         round(ra,6), round(dec,6), params['epoch'], source_radec, round(bcvel_baryc,4), source_obs, 
                         params['site'], obnames, mjd, round(bjd,5), params['pmra'], params['pmdec'] ))
       
    return params, bcvel_baryc, mephem, im_head

def find_shift_between_wavelength_solutions(params, wave_sol_1, wave_sol_lines_1, wave_sol_2, wave_sol_lines_2, spectra, names=['first','second']):
    """
    Finds the conversion polynomial between the two solutions
    :param wave_sol_1, wave_sol_2: 2d array of floats, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    :param wave_sol_lines_1, wave_sol_lines_2: 2d array of floats with one line for each order. Each order contains the wavelengths of the identified reference lines, 
                                        sorted by the brightest to faintest (in spectrum from which the solution is derived).
                                        To create the same number of lines for each order, the array is filled with 0
    :param spectra: 2d array of floats, spectrum, only needed for the shape of the wavelengths array
    """
    # make it into wavelength, find the wavelengths of the first arclines in both, find the wavelengths of the second arclines in both, find the overlapping lines -> directly to lambda-shift -> RV
    # Make the solutions into wavelengths
    wavelengths1, dummy = create_wavelengths_from_solution(params, wave_sol_1, spectra)
    wavelengths2, dummy = create_wavelengths_from_solution(params, wave_sol_2, spectra)
    # find the closest pixel of each catalog line
    result = []
    for order in range(spectra.shape[0]):
        catalog_lines1 = wave_sol_lines_1[order,:]
        catalog_lines1 = catalog_lines1[catalog_lines1 > 100]          # ignore the entries necessary to fill the array
        catalog_lines2 = wave_sol_lines_2[order,:]
        catalog_lines2 = catalog_lines2[catalog_lines2 > 100]          # ignore the entries necessary to fill the array
        catalog_lines3 = list( set(catalog_lines1) & set(catalog_lines2) )      # overlapping entries
        for [i, catalog_lines] in [ [0, catalog_lines3], [1, catalog_lines1], [2, catalog_lines2] ]:    # use the different catalogs
            for catalog_line in catalog_lines:
                posi = [0, 0]
                for [j, wavelengths] in [ [0, wavelengths1], [1, wavelengths2] ]:                       # use both wavelength solutions
                    diff = catalog_line - wavelengths[order,:]
                    pos = np.argmin(np.abs(diff))                                                       # closest full pixel to the catalog line
                    temp = wavelengths[ order, max(0,pos-1):min(pos+1,wavelengths.shape[1])+1 ]
                    res = (temp[-1] - temp[0]) / (len(temp) - 1.)                                        # resolution around the pixel
                    posi[j] = pos + diff[pos]/res                                                       # fraction of the pixel
                    #print res, pos, posi[j], i, j, len(posi), len(diff)
                    """ compares dlamba/lambda -> RV -> big scatter
                    if (i == 1 and j == 1) or (i == 2 and j == 0):                                      # don't apply catalog lines from the first solution to the second solution
                        continue
                    diff = np.abs(catalog_line - wavelengths[order,:])
                    pos = np.argmin(diff)                                 # closest pixel
                    rel_diff = 2. * (wavelengths1[order,pos] - wavelengths2[order,pos]) / (wavelengths1[order,pos] + wavelengths2[order,pos])
                    result.append([ rel_diff,i,j,order,catalog_line,pos ])
                    """
                    #print rel_diff,i,j,order,catalog_line
                result.append([ posi[1]-posi[0], i, 0, order, catalog_line, np.mean(posi) ])
    result = np.array(result)
    """for i in range(10):
        if i==0:    good_values = range(len(result))        # all data
        if i==1:    good_values = result[:,1] == 0          # only same lines in both solutions
        if i==2:    good_values = result[:,1] == 1          # lines of solution1
        if i==3:    good_values = result[:,1] == 2          # lines of solution2
        if i==4:    good_values = (result[:,1] == 0) & (result[:,5] < 500)          # only same lines in both solutions and left side 
        if i==5:    good_values = (result[:,1] == 0) & (result[:,5] > 1500)          # only same lines in both solutions and right side 
        if i==6:    good_values = (result[:,1] == 0) & (result[:,5] > 700) & (result[:,5] < 1300)          # only same lines in both solutions and middle
        if i==7:    good_values = (result[:,1] == 0) & (result[:,5] > 700) & (result[:,5] < 1300) & (result[:,3] < 20)          # only same lines in both solutions and middle and red orders
        if i==8:    good_values = (result[:,1] == 0) & (result[:,5] > 700) & (result[:,5] < 1300) & (result[:,3] > 60)          # only same lines in both solutions and middle and blue orders
        if i==9:    good_values = (result[:,1] == 0) & (result[:,5] > 700) & (result[:,5] < 1300) & (result[:,3] > 30) & (result[:,3] < 50)         # only same lines in both solutions and middle and middle orders
        print 'i, size, median, average, std:', i, result[good_values,0].shape, np.median(result[good_values,0]), np.mean(result[good_values,0]), np.std(result[good_values,0], ddof=1),
        for pctl in [1,5,10,20]:
            temp = percentile_list(result[good_values,0], pctl/100.)
            print '\tpctl, size, median, mean, std:', pctl, len(temp), np.median(temp), np.mean(temp), np.std(temp, ddof=1),
        print "" """
    if np.sum(result[:,1] == 0) >= 0.1* result.shape[0] and np.sum(result[:,1] == 0) >= 300:        # enough data
        good_values = result[:,1] == 0          # only same lines in both solutions
        text = 'The {0} overlapping identified lines between both solutions have been used.'.format(np.sum(good_values) )
    else:
        good_values = (result[:,1] != -1)         # all data, means the overlapping lines have twice the weight
        text = 'All available lines ({0}) been used, the overlapping lines weighted double.'.format(np.sum(good_values) )
    shift = np.median(result[good_values,0])
    shift_err = np.std(result[good_values,0], ddof=1)
    logger('Info: The shift between the two wavelength solutions ({4} - {3}) is {0} +/- {1} pixel. {2}'.format( round(shift,4), round(shift_err,4), text, names[0], names[1] ))
    
    return shift, shift_err

def remove_orders_low_flux(params, flat_spec_norm):
    number_no_data = np.sum(np.isnan(flat_spec_norm[2,:,:]), axis=1)                                # What orders don't contain enough data
    limit_no_data = (1-params['extraction_min_ident_part_of_trace_percent']/100.)*flat_spec_norm.shape[2]
    keep_orders  = np.where(number_no_data <= limit_no_data)[0]         # remove all orders where more than x% of the pixel are not available
    remove_orders = np.where(number_no_data > limit_no_data)[0]
    if len(remove_orders) > 0:
        logger('Warn: The blaze function for {0} orders contains not enough data points (limit: {2}). The following orders have been removed: {1}'.format(len(remove_orders), remove_orders, limit_no_data))
        #print wavelength_solution.shape, wavelength_solution_arclines.shape, sci_tr_poly.shape, xlows.shape, xhighs.shape, widths.shape, cal_tr_poly.shape, axlows.shape, axhighs.shape, awidths.shape, flat_spec_norm.shape
    return remove_orders, keep_orders

def header_results_to_texfile(params, header_keywords=[]):
    if len(header_keywords) == 0:
        header_keywords.append(['HIERARCH HiFLEx OBJNAME',      'Object name',                  ''  ])
        header_keywords.append(['HIERARCH HiFLEx EXPOSURE',     'Exposure time',                '[s]'])
        header_keywords.append(['HIERARCH HiFLEx DATE-OBS',     'UTC, Begin of expose',         ''  ])
        header_keywords.append(['HIERARCH HiFLEx DATE-MID',     'Middle of exposure',           ''  ])
        header_keywords.append(['HIERARCH HiFLEx JD',           'JD at middle of exposure',     '[d]'])
        header_keywords.append(['HIERARCH HiFLEx fsum_all',     'Total extracted flux',         '[ADU]'])
        header_keywords.append(['HIERARCH HiFLEx BCKNOISE',     'Background noise',             '[ADU]'])
        header_keywords.append(['HIERARCH HiFLEx CD_SHIFT',     'Offset in Cross-Dispersion',   '[px]'])
        header_keywords.append(['HIERARCH HiFLEx CD_S_WDTH',    'Width of the offset in CD',    '[px]'])
        header_keywords.append(['HIERARCH HiFLEx D_SHIFT',      'Offset in Dispersion (incldung DT)',           '[px]'])
        header_keywords.append(['HIERARCH HiFLEx D_SHIFT_ERR',  'Uncertainty of the offset in D',               '[px]'])
        header_keywords.append(['HIERARCH HiFLEx D_SHIFT_NUMBER_LINES', 'Number of lines used to calculate the shift', ''])
        header_keywords.append(['HIERARCH HiFLEx D_WIDTH',      'Gaussian width of the calibration lines',      '[px]'])
        header_keywords.append(['HIERARCH HiFLEx D_SHIFT_KMS',  'Offset in Dispersion',                         '[km/s]'])
        header_keywords.append(['HIERARCH HiFLEx DT_SHIFT',     'Offset between wavelength solutions',          '[px]'])
        header_keywords.append(['HIERARCH HiFLEx DT_SHIFT_KMS', 'Offset between wavelength solutions',          '[km/s]'])
        header_keywords.append(['HIERARCH HiFLEx BJDTDB',       'Barycentric correct JD (incl. leap seconds)',  '[d]'])
        header_keywords.append(['HIERARCH HiFLEx CERES RV',     'CERES RV (not corrected for BVC)', '[km/s]'])
        header_keywords.append(['HIERARCH HiFLEx CERES RV_BARY','CERES RV (corrected for BVC)', '[km/s]'])
        header_keywords.append(['HIERARCH HiFLEx CERES RV_ERR', 'CERES RV error',               '[km/s]'])
        header_keywords.append(['HIERARCH HiFLEx BCV',          'Barycentric Velocity',         '[km/s]'])
        header_keywords.append(['HIERARCH HiFLEx CERES BS',     'CERES Bisector',               ''  ])
        header_keywords.append(['HIERARCH HiFLEx CERES BS_ERR', 'CERES Bisector error',         ''  ])
        header_keywords.append(['HIERARCH HiFLEx TERRA RV',     'TERRA RV relative to template',    'm/s'  ])
        header_keywords.append(['HIERARCH HiFLEx TERRA RV_ERR', 'TERRA RV error',                   'm/s'  ])
        header_keywords.append(['HIERARCH HiFLEx SERVAL RV',    'SERVAL RV relative to template',   'm/s'  ])
        header_keywords.append(['HIERARCH HiFLEx SERVAL RV_ERR', 'SERVAL RV error',                 'm/s'  ])
        
    results = []
    # Create the table header
    for ii, start in enumerate(['Header keyword:', 'Filename', '']):
        result = [start]
        for header_keyword in header_keywords:
            result.append(header_keyword[ii])
        results.append(result)
    # Read the files
    files = sorted(os.listdir(params['path_extraction']))
    files = [os.path.join(params['path_extraction'], f) for f in files] # add path to each file
    #files.sort(key=os.path.getmtime)       # not useful with the RV added in the end
    for file in files:
        if not file.endswith(".fits"):
            continue
        if file.find('wavelength_solution') != -1 or file.find('_wavesol_') != -1:
            continue
        result = [ file.replace(params['path_extraction'],'') ]         # file name
        im_head = fits.getheader(file)
        for header_keyword in header_keywords:
            result.append(str(im_head.get(header_keyword[0], '')))
        results.append(result)
         
    if len(results) > 2:
        with open('measurement_table.csv','w') as file:
            for entry in results:
                file.write("\t".join(entry)+'\n')
        logger('Info: Created {0}'.format('measurement_table.csv'))
     
def rv_results_to_hiflex(params):
    """
    Transforming the TERRA and SERVAL results into spectra compatible with HiFLEx
    """
    obj_names = []
    for root, dirs, files in os.walk(params.get('path_rv_terra', 'dummypath'), followlinks=True):                       # Find all the objects again, as won't be added to obj_names when re-run
            for file in files:
                if file == 'template_order_00':                    # TERRA template
                    [root1, obj_name, dummy] = (os.sep+root).rsplit(os.sep,2)
                    if obj_name not in obj_names:   obj_names.append(obj_name)
                    break
                if file == 'synthetic.rv':                  # TERRA resuts
                    [root1, obj_name, dummy] = (os.sep+root).rsplit(os.sep,2)
                    if obj_name not in obj_names:   obj_names.append(obj_name)
                    break
    #for obname in obj_names:
    #    convert_terra_master_hiflex(params, obname)
       
    #obj_names = []
    for root, dirs, files in os.walk(params.get('path_rv_serval', 'dummypath'), followlinks=True):                       # Find all the objects again, as won't be added to obj_names when re-run
            for file in files:
                if file == 'template.fits':                 # SERVAL template
                    [root1, obj_name] = (os.sep+root).rsplit(os.sep,1)
                    if obj_name not in obj_names:   obj_names.append(obj_name)
                    break
                if file.endswith('.rvc.dat'):               # SERVAL results
                    [root1, obj_name] = (os.sep+root).rsplit(os.sep,1)
                    if obj_name not in obj_names:   obj_names.append(obj_name)
                    break 
    terra_rvs, serval_rvs = [], []
    for obname in obj_names:
        convert_serval_master_hiflex(params, obname)
        convert_terra_master_hiflex(params, obname)
        terra_rvs  += get_terra_results(params, obname)     # Object name, index in file, mid-exposure JD, RV, error RV
        serval_rvs += get_serval_results(params, obname)
        
    add_terra_serval_results_to_header(params, terra_rvs, serval_rvs)
        
def add_terra_serval_results_to_header(params, terra_rvs, serval_rvs):
    terra_jd = np.array( [terra_results[2] for terra_results in terra_rvs] )
    terra_obj = np.array( [terra_results[0] for terra_results in terra_rvs] )
    serval_bjd = np.array( [serval_results[2] for serval_results in serval_rvs] )
    serval_obj = np.array( [serval_results[0] for serval_results in serval_rvs] )
    for file_RV in sorted(os.listdir(params['path_extraction'])):
        if not file_RV.endswith(".fits"):           continue
        im = fits.open(params['path_extraction']+file_RV)
        im_head = im[0].header
        if not 'HiFLEx OBJNAME' in im_head.keys():  continue
        if not 'HiFLEx JD' in im_head.keys():       continue
        jd = im_head['HiFLEx JD']
        bjd = im_head.get('HiFLEx BJDTDB', jd)
        if len(terra_obj) > 0:
            good_data_terra = np.where( terra_obj == im_head['HiFLEx OBJNAME'].lower() )[0]
        else:       good_data_terra = np.zeros(0)
        if len(serval_obj) > 0:
            good_data_serval = np.where( serval_obj == im_head['HiFLEx OBJNAME'].lower() )[0]
        else:       good_data_serval = np.zeros(0)
        changed_header = False
        if good_data_terra.shape[0] > 0:
            diff = np.abs( terra_jd[good_data_terra] - jd )
            if np.min(diff) < 1E-6:         # 1E-4: 9s
                index = good_data_terra[diff.argmin()]
                data = terra_rvs[index]
                im_head['HIERARCH HiFLEx TERRA RV'] = (round(data[3],1), 'RV in m/s (rel. to template: {0})'.format(data[0]))
                im_head['HIERARCH HiFLEx TERRA RV_ERR'] = (round(data[4],1), 'Uncertainty RV in m/s')
                changed_header = True
                terra_rvs[index][1] = -99
            elif False:
                print('not found in terra results', im_head['HiFLEx OBJNAME'], jd, file_RV, diff)#, terra_jd[good_data_terra])
        if good_data_serval.shape[0] > 0:
            diff = np.abs( serval_bjd[good_data_serval] - bjd )
            if np.min(diff) < 5E-6:         # 1E-4: 9s
                index = good_data_serval[diff.argmin()]
                data = serval_rvs[index]
                im_head['HIERARCH HiFLEx SERVAL RV'] = (round(data[3],1), 'RV in m/s (rel. to template: {0})'.format(data[0]))
                im_head['HIERARCH HiFLEx SERVAL RV_ERR'] = (round(data[4],1), 'Uncertainty RV in m/s')
                changed_header = True
                serval_rvs[index][1] = -99
            elif False:
                print('not found in serval results', im_head['HiFLEx OBJNAME'], bjd, file_RV, diff)#, serval_bjd[good_data_serval] )
        if changed_header:
            save_multispec(im[0].data, params['path_extraction']+file_RV, im_head)
    for entry in terra_rvs:
        if entry[1] >= 0:
            print("problem for Ronny: not assigned the TERRA result:", entry)
    for entry in serval_rvs:
        if entry[1] >= 0:
            print("problem for Ronny: not assigned the SERVAL result:", entry)

def get_terra_results(params, obname):
    """
    Reads the result file from TERRA and puts it together with the corrected JD into a list
    :return dataf: list with following entries: obname, index in file, JD, RV, RV error
    """
    filename = params['path_rv_terra']+obname+'/results/synthetic.rv'
    if not os.path.isfile(filename):
        logger('Warn: TERRA result file {0} does not exist'.format(filename))
        return []
    dataf = read_text_file(filename, no_empty_lines=True)
    dataf = convert_readfile(dataf, [float, float, float], delimiter=' ', replaces=[['  ',' ']]*20, ignorelines=['NaN'])     # MJD is unfortunatelly not the same as JD-2400000.5 and the offsets varies between nights or within one data set
    midexpJD = []
    for file in sorted(os.listdir(params['path_rv_terra']+obname+'/data/')):
        if file.endswith(".csv"):
            date = datetime.datetime.strptime(file.replace('.csv',''), '%Y-%m-%d%H%M%S')
            midexpJD.append(get_julian_datetime(date))
    midexpJD = sorted(midexpJD)
    diff = len(midexpJD) - len(dataf)
    if diff == 0:                 # all files were used -> wrong MJD is not a problem
        for ii in range(len(dataf)):
            #print midexpJD[ii]-dataf[ii][0]        # See the offset between wrong MJD and JD
            dataf[ii] = [obname, ii, midexpJD[ii]] + dataf[ii][1:3]
    else:                                           # TERRA excluded fles -> Now it really does get difficult because of the wrong MJD
        midexpJD = np.array(midexpJD)
        possible_results = [0.0, 0.6666666, 0.8333333]
        prev_offset = 0
        for ii in range(len(dataf)):
            diffs = midexpJD[ii:min(ii+diff+1,midexpJD.shape)] - dataf[ii][0]
            found = []
            for jj, entry in enumerate(possible_results):
                found += list(np.where(np.abs( np.mod(diffs, 1) - entry ) < 1E-5)[0])       # find the data which could be right
            if len(found) == 1:
                found = found[0]
                dataf[ii] = [obname, ii, midexpJD[ii+found]] + dataf[ii][1:3]
                prev_offset = found
            else:
                dataf[ii] = [obname, ii, midexpJD[ii+prev_offset]] + dataf[ii][1:3]
                print( "Struggled to assign the MJDs from the TERRA. This is for Ronny to check things:",ii, found, diffs )

    return dataf

def get_serval_results(params, obname):
    """
    Reads the SERVAL results and provides the results in a list
    :return dataf: list with following entries: obname, index in file, BJD, RV, RV error
    """
    filename = params['path_rv_serval']+obname+'/'+obname+'.rvc.dat'
    if not os.path.isfile(filename):
        logger('Warn: SERVAL result file {0} does not exist'.format(filename))
        return []
    dataf = read_text_file(filename, no_empty_lines=True)
    dataf = convert_readfile(dataf, [float, float, float], delimiter=' ', shorten_input=True)
    for ii in range(len(dataf)):
        dataf[ii] = [obname, ii] + dataf[ii][:3]
        
    return dataf   

def convert_terra_master_hiflex(params, obname):
    """
    Reads the TERRA template and safes it in the HiFLEx format
    """
    spec = []
    wave = []
    terra_template = params['path_rv_terra']+obname+'/template/template_order_'
    for ii in range(1000):
        if os.path.isfile(terra_template+'%2.2i'%ii):
            dataf = read_text_file(terra_template+'%2.2i'%ii, no_empty_lines=True)
            dataf = convert_readfile(dataf, [float, float], delimiter=' ', replaces=[['  ',' ']]*10, shorten_input=True)
            dataf = np.array(dataf)
            spec.append(dataf[:,1])
            wave.append(dataf[:,0])
        else:
            break
    if ii > 0:
        wave, spec = np.array(wave), np.array(spec)
        hdu = fits.PrimaryHDU()
        spec = convert_rv_templates_hilfex_normalise(params, spec, wave)        # Creates the list with all the files
        im_name = 'terra_template_hiflexformat_'+obname
        im_head = hdu.header
        
        save_multispec(np.array(spec), params['path_rv_terra']+im_name, im_head, bitpix=params['extracted_bitpix'])
        
        spectra = spec[1]
        im_head_bluefirst = copy.copy(im_head)
        im_head = add_specinfo_head(spectra, spectra, spec[6], None, im_head)
        im_head_bluefirst = add_specinfo_head(spectra[::-1,:], spectra[::-1,:], spec[6][::-1,:], None, im_head_bluefirst)
        doo = dict(spec=True, blaze=False, norm=True, blue=False, harps=True)
        params['path_extraction_single'] = params['path_rv_terra']
        im_name = 'terra_template_single_'+obname
        save_single_files(params, [obname], im_name, im_head, im_head_bluefirst, spectra, spec[3], spec[5], [], spec[0], calimages['wave_sol_dict_sci']['wavesol'], doo)
        
        logger('Info: Created {0} from {1}*'.format(params['path_rv_terra']+'terra_template_hiflexformat_'+obname+'.fits', terra_template))
    else:
        logger('Warn: TERRA template {0} does not exist'.format(terra_template))

def convert_serval_master_hiflex(params, obname):
    """
    Reads the SERVAL template and safes it in the HiFLEx format
    """
    serval_template = params['path_rv_serval']+obname+'/template.fits'
    if not os.path.isfile(serval_template):
        logger('Warn: SERVAL template {0} does not exist'.format(serval_template))
        return
    im = fits.open(serval_template)
    im_head = im[0].header
    spec = np.array(im[1].data, dtype=np.float64)
    wave = np.array(im[2].data, dtype=np.float64)
    if np.max(wave, axis=None) < 15:                            # If np.log(wave)
        wave = np.exp(wave)
    wave = wavelength_vacuum_to_air(wave)                       # Use Air Wavelengths
    #wave *= (1+im_head.get('HiFLEx BCV', 0)/(Constants.c/1e3))  # Use barycentric velocity, but it's already in BCV
    spec = convert_rv_templates_hilfex_normalise(params, spec, wave)
    im_name = 'serval_template_hiflexformat_'+obname
    save_multispec(np.array(spec), params['path_rv_serval']+im_name, im_head, bitpix=params['extracted_bitpix'])
    
    spectra = spec[1]
    im_head_bluefirst = copy.copy(im_head)
    im_head = add_specinfo_head(spectra, spectra, spec[6], None, im_head)
    im_head_bluefirst = add_specinfo_head(spectra[::-1,:], spectra[::-1,:], spec[6][::-1,:], None, im_head_bluefirst)
    doo = dict(spec=True, blaze=False, norm=True, blue=False, harps=True)
    params['path_extraction_single'] = params['path_rv_serval']
    im_name = 'serval_template_single_'+obname
    save_single_files(params, [obname], im_name, im_head, im_head_bluefirst, spectra, spec[3], spec[5], [], spec[0], calimages['wave_sol_dict_sci']['wavesol'], doo)
        
    logger('Info: Created {0} from {1}'.format(params['path_rv_serval']+'serval_template_hiflexformat_'+obname+'.fits', serval_template))
    
def convert_rv_templates_hilfex_normalise(params, spec, wave, limit_for_flat=10):
    """
    Creates an array with al the important information to keep compatible with HiFLEx output, e.g.: Wave, Spec, Error Spec, Flatcorrected, Error Flat, Normalised, noise
    !!! flat correction doesn't work:   TERRA, SERVAL: different amount of x-pixel
    """
    global calimages
    if not'flat_spec_norm' in calimages.keys():
        flat_spec_norm = np.array(fits.getdata(params['master_blaze_spec_norm_filename']))
        remove_orders, keep_orders = remove_orders_low_flux(params, flat_spec_norm)
        flat_spec_norm = flat_spec_norm[:,keep_orders,:] 
        calimages['flat_spec_norm'] = flat_spec_norm
    
    if np.nanmedian(spec, axis=None) > limit_for_flat and spec.shape[0] == calimages['flat_spec_norm'].shape[1] and spec.shape[1] == calimages['flat_spec_norm'].shape[2]:     # same number of orders, same number of pixel
        fspec = spec/(calimages['flat_spec_norm'][2]+0.0)        # 1: extracted flat, 2: low flux removed
    else:
        fspec = spec
    measure_noise_orders = 16
    measure_noise_semiwindow = 10                   # in pixel
    efspec = measure_noise(fspec, p_order=measure_noise_orders, semi_window=measure_noise_semiwindow)             # Noise will be high at areas wih absorption lines
    cspec, ecspec = normalise_continuum(fspec, wave, nc=6, semi_window=measure_noise_semiwindow, nc_noise=measure_noise_orders) 
    #if np.nanmedian(spec, axis=None) <= limit_for_flat:
    #    cspec = spec
    
    return [wave, spec, spec*0, fspec, efspec, cspec, ecspec ]

def load_ceres_modules(params):
        global GLOBALutils
        global correlation
        global lowess
         
        base = params['path_ceres']
        # Import the routines from CERES:
        sys.path.append(base+"utils/Correlation")
        sys.path.append(base+"utils/GLOBALutils")
        sys.path.append(base+"utils/OptExtract")    # for Marsh, at least
        sys.path.append(base+"utils/CCF")           # needed by GLOBALutils.py
        with np.errstate(invalid='ignore'):
            import GLOBALutils
        import correlation
        # Import other stuff
        import statsmodels.api as sm
        lowess = sm.nonparametric.lowess

def stellar_params_ceres_multicore(kwargs):
    
    [ii, spec2, Rx] = kwargs
    
    good_values = ~np.isnan(spec2[5,ii,:])
    if sum(good_values) < 10:                               # If not enough values, then ignore the order
        good_orders = False
    else:
        good_orders = True
        spec2[5,ii,good_values] = GLOBALutils.convolve(spec2[0,ii,good_values],spec2[5,ii,good_values],Rx)
        
    return good_orders, spec2[5,ii,:]

def stellar_params_ceres(params, spec, fitsfile, wavelength_solution):
    
    base = params['path_ceres']
    specs = spec.shape
    models_path = base + 'data/COELHO_MODELS/R_40000b/'
    fsim = fitsfile
    #RESI = 120000.
    RESI = np.median(wavelength_solution[:,-1]/wavelength_solution[:,-2])     # lambda/d_lambda for the central pixel, for MRES that doesn't seem to make difference
    npools = max(1, params['use_cores'])

    force_stellar_pars = False
   
    """ Data description of the file
    0: wavelength for each order and pixel'
    1: extracted spectrum'
    2: measure of error (photon noise, read noise)'
    3: flat corrected spectrum'
    4: error of the flat corrected spectrum (residuals to a {0} order polynomial)'.format(measure_noise_orders)
    5: continuum normalised spectrum'
    6: error in continuum (fit to residuals of {0} order polynomial)'.format(measure_noise_orders)
    7: Mask with good areas of the spectrum: 0.1=saturated_px, 0.2=badpx'
    8: spectrum of the emission line lamp'
    """
    
    for order in range(spec.shape[1]):
        L  = np.where( (spec[1,order,:] != 0) & (~np.isnan(spec[1,order,:])) )              # good values
        with np.errstate(invalid='ignore'):
            LL = np.where(spec[5,order,:] > 1 + 10. / scipy.signal.medfilt(spec[2,order,:],21))[0]          # remove emission lines and cosmics
        spec[5,order,LL] = 1.
    #plot_img_spec.plot_spectra_UI(np.array([spec]))
    T_eff, logg, Z, vsini, vel0 = 5777, 4.4374, 0.0134, 2, 0
    hardcoded = False
    pars_file = params['path_rv_ceres'] + fsim+'_stellar_pars.txt'                    # calculate stellar parameters for each spectrum
    if ( not os.path.isfile(pars_file) or force_stellar_pars ):
        req1 = (np.nanmax(spec[0,:,:],axis=1) < 5150).any()         # At least one order with values below 5150
        req2 = (np.nanmax(spec[0,:,:],axis=1) > 5600).any()
        if req1 and req2:
            good_orders = []
            Rx = np.around(1./np.sqrt(abs(1./40000.**2 - 1./RESI**2)))
            spec2 = copy.copy(spec)
            kwargs = [[ii, spec2, Rx] for ii in range(specs[1])]
            results = []
            if params['use_cores'] > 1:
                p = multiprocessing.Pool(params['use_cores'])
                results = p.map(stellar_params_ceres_multicore, kwargs)
                p.terminate()
            else:
                for kwarg in kwargs:
                    results.append( stellar_params_ceres_multicore(kwarg) )
            for ii,result in enumerate(results):
                good_orders.append(result[0])
                spec2[5,ii,:] = result[1]
            spec2 = spec2[:,good_orders,:]
            
            T_eff, logg, Z, vsini, vel0, ccf = correlation.CCF(spec2,model_path=models_path,npools=npools, base=base+'utils/Correlation/')    # uses scipy.integrate.simps which can create negative values which makes math.sqrt(<0) fail (https://stackoverflow.com/questions/36803745/python-simpsons-rule-negative-answer-for-positive-area-under-the-curve); but don't use try - except, as correlation.CCF is called with Pool and crashes. When it crashes the pool isn't closed/terminated, hence processes are building up, try except doesn't solve this
            line = "%6d %4.1f %4.1f %8.1f %8.1f\n" % (T_eff,logg, Z, vsini, vel0)
            with open(pars_file, 'w') as f:
                f.write(line)
            text = 'Info: Using the following atmosperic parameters for {0}'.format(fsim)
        else:
            text = 'Warn: could not determine the stelar parameters as the wavelength range below 5150 is not available. Using the hard coded values'
            hardcoded = True
    elif os.path.isfile(pars_file):
        T_eff, logg, Z, vsini, vel0 = np.loadtxt(pars_file,unpack=True)
        text = 'Info: Atmospheric parameters loaded from file {0}'.format(pars_file)
    else:
        hardcoded = True
    logger('{0}: T_eff, logg, Z, vsini, vel0: {1}, {2}, {3}, {4} {5}.'.format( text, T_eff, logg, Z, vsini, round(vel0,4) ))
    
    return T_eff, logg, Z, vsini, vel0, hardcoded
    
def rv_analysis_ceres(params, spec, fitsfile, obname, reffile, mephem, vsini):
    
    def get_spec(spectra, waves, cen_wave, range_px):
        """
        Used to get the order in which the cen_wave is most central
        :param spectra: 1d or 2d array of floats with the spectral data
        :param waves: same format as spectra, corresponding wavelengths
        :param cen_wave: wavelength to get the spectra around
        :param range_px: how many pixel around that wavelength
        :return spec, waves: returns the spectra and the wavelengths of the area specified
        """
        if len(spectra.shape) < 2:                      # make into a 2d spectrum, if necessary
            spectra = spectra.reshape(1, spectra.shape[0])
            waves = waves.reshape(1, waves.shape[0])
        wavediff = np.nanmax( np.abs( waves[:,1:] - waves[:,:-1] ))     # maximum wavelengths difference between consecutive pixel
        diff = waves - cen_wave
        pos = np.where( np.abs(diff) <= wavediff )      # gives a 2d array for the positions (order, px)
        diffp = np.abs(pos[1] - spectra.shape[1]/2)     # what is the most central px
        for dummy in np.unique(pos[0]):                 # unique ignores if wavelength is not covered in the spectrum
            pos2 = int(np.argmin(diffp))                # what is the most central pixel
            order, px = pos[0][pos2], pos[1][pos2]
            spec = spectra[order, max(0,px-range_px) : min(spectra.shape[1], px+range_px+1) ]
            if len( spec[~np.isnan(spec)] ) == 0:           # no spectrum for this data
                pos[1][pos[0] == order] = 1E9               # -> check the others orders
                diffp = np.abs(pos[1] - spectra.shape[1]/2)
            else:
                #print spec, waves[order, max(0,px-range_px) : min(spectra.shape[1], px+range_px+1) ]
                return spec, waves[order, max(0,px-range_px) : min(spectra.shape[1], px+range_px+1) ]
    
        return [], []                                      # everything went wrong

    base = params['path_ceres']
    specs = spec.shape
    spec = np.vstack(( spec, np.zeros([11-specs[0], specs[1], specs[2]]) ))
    spec[7:11,:,:] = np.nan
    lbary_ltopo = 1
    refvel = 0
    know_moon = False
    here_moon = False
    fsim = fitsfile
    
    #ron1,gain1 = h[1].header['HIERARCH ESO DET OUT1 RON'],h[1].header['HIERARCH ESO DET OUT1 GAIN']
    #ron2,gain2 = h[2].header['HIERARCH ESO DET OUT1 RON'],h[2].header['HIERARCH ESO DET OUT1 GAIN']
    #halfcounts = h[0].header['HIERARCH ESO INS DET1 TMMEAN']

    """ Data description of the file
    0: wavelength for each order and pixel'
    1: extracted spectrum'
    2: measure of error (photon noise, read noise)'
    3: flat corrected spectrum'
    4: error of the flat corrected spectrum (residuals to a {0} order polynomial)'.format(measure_noise_orders)
    5: continuum normalised spectrum'
    6: error in continuum (fit to residuals of {0} order polynomial)'.format(measure_noise_orders)
    7: Mask with good areas of the spectrum: 0.1=saturated_px, 0.2=badpx'
    8: spectrum of the emission line lamp'
    """
    #spec[6,np.isnan(spec[6,:])] = 0     # doesn't make a difference
    #spec[2,np.isnan(spec[2,:])] = 0     # doesn't make a difference
    for order in range(spec.shape[1]):
        L  = np.where( (spec[1,order,:] != 0) & (~np.isnan(spec[1,order,:])) )              # good values
        #ratio              = np.polyval(ccoefs[order],spec[0,order,:][L])*Rnorms[order]
        with np.errstate(invalid='ignore'):
            ratio = spec[1,order,:][L] / spec[5,order,:][L]                                     # ratio between extracted spectrum and continuum normalised spectrum -> blaze function, cancels absorption lines
        spec[7,order,:][L] = ratio
        #spec[8,order,:][L] = spec[6,order,:][L]                                             # error continuum (first guess), but not good. #sn_order=8
        spec[8,order,:][L] = spec[2,order,:][L]                                             # error of the extracted data, depending on what is used, the RV changes by few 100 km/s -> > several \AA
        #spec[8,order,:][L] = ratio * R_flat_ob_n[order,1,:][L] / np.sqrt( ratio * R_flat_ob_n[order,1,:][L] / gain2 + (ron2/gain2)**2 )        # something like S/N -> used as this by XCor
        #spl           = scipy.interpolate.splrep(np.arange(WavSol.shape[0]), WavSol,k=3)
        #dlambda_dx    = scipy.interpolate.splev(np.arange(WavSol.shape[0]), spl, der=1)
        #NN            = np.average(dlambda_dx)
        #dlambda_dx    /= NN
        with np.errstate(invalid='ignore'):
            LL = np.where(spec[5,order,:] > 1 + 10. / scipy.signal.medfilt(spec[8,order,:],21))[0]          # remove emission lines and cosmics
        spec[5,order,LL] = 1.
        spec[9,order,:][L] = spec[5,order,:][L]# * (dlambda_dx[L] ** 1)         # used for the analysis in XCor (spec_order=9, iv_order=10)
        spec[10,order,:][L] = spec[2,order,:][L]# / (dlambda_dx[L] ** 2)        # used for the analysis in XCor (spec_order=9, iv_order=10)
    #plot_img_spec.plot_spectra_UI(np.array([spec]))
    
    
    # assign mask, obname is the name of the object
    sp_type, mask = GLOBALutils.get_mask_reffile(obname,reffile=reffile,base=base+'data/xc_masks/')     # !!! Warn: upper and lower case matters for obname
    logger('Info: Will use {0} mask for CCF.'.format(sp_type))
    
    spec[5,np.isnan(spec[5,:])] = 0     # 5 -> 9, GLOBALutils.XCor expects 0, not NaN -> otherwise higher RV scatter; but CCF expects nan
    
    # Read in mask
    ml, mh, weight = np.loadtxt(mask,unpack=True)
    ml_v = GLOBALutils.ToVacuum( ml )
    mh_v = GLOBALutils.ToVacuum( mh )
    av_m = 0.5*( ml_v + mh_v )
    mask_hw_kms = (GLOBALutils.Constants.c/1e3) * 0.5*(mh_v - ml_v) / av_m

    disp = GLOBALutils.get_disp(obname, reffile=reffile)        # !!! Warn: upper and lower case matters, disp is "velocity width in km/s that is used to broaden the lines of the binary mask. It should be similar to the standard deviation of the Gaussian that is fitted to the CCF."
    if disp == 0:
        known_sigma = False
        if vsini != -999 and vsini != 0.:
            disp = vsini
        else:
            disp = 3.
    else:
        known_sigma = True

    mask_hw_wide = av_m * disp / (GLOBALutils.Constants.c/1.0e3)
    ml_v = av_m - mask_hw_wide
    mh_v = av_m + mask_hw_wide 

    logger('Step: Computing the CCF... for the RV measurement')
    cond = True

    if sp_type == 'M5':
        moon_sig = 4.5
    elif sp_type == 'K5':
        moon_sig = 4.2
    else:
        moon_sig = 4.0

    while (cond):
        # first rough correlation to find the minimum
        #   spec: spectrum in the form [data type, orders, pixel]
        #   ml_v, mh_v: masks 
        #   lbary_ltopo = 1.0 + res['frac'][0]
        # vel_width=600 instead of vel_width=300 changes the results by up to a few 100m/s in MRES
        # max_vel_rough=600 instead of max_vel_rough=300 changes the results by up to a few 100m/s in MRES
        vels, xc_full, sn, nlines_ccf, W_ccf = \
                GLOBALutils.XCor(spec, ml_v, mh_v, weight,\
                0, lbary_ltopo, vel_width=300, vel_step=3,\
                spec_order=9, iv_order=10, sn_order=8, max_vel_rough=300)
        #print(vels, xc_full, sn, nlines_ccf, W_ccf)
        # W_ccf is a weight
        
        xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=3.0, Simple=True, W=W_ccf)
        #print 'xc_av, vels, xc_full, sn, nlines_ccf, W_ccf',xc_av, vels, xc_full, sn, nlines_ccf, W_ccf
        
        # Normalize the continuum of the CCF robustly with lowess     
        yy = scipy.signal.medfilt(xc_av,11)
        pred = lowess(yy, vels,frac=0.4,it=10,return_sorted=False)
        tck1 = scipy.interpolate.splrep(vels,pred,k=1)
        xc_av_orig = xc_av.copy()
        xc_av /= pred
        vel0_xc = vels[ np.argmin( xc_av ) ] 
            
        rvels, rxc_av, rpred, rxc_av_orig, rvel0_xc = \
                vels.copy(), xc_av.copy(), pred.copy(),\
                xc_av_orig.copy(), vel0_xc

        xc_av_rough = xc_av
        vels_rough  = vels
            
        vel_width = np.maximum( 20.0, 6*disp )                      # Adjusted in order to avoid crashes because of unphysical disp
        #print vel_width, disp, vsini       # problem with vel_width, due to disp, due to p1gau below
        vels, xc_full, sn, nlines_ccf, W_ccf =\
                GLOBALutils.XCor(spec, ml_v, mh_v, weight,\
                vel0_xc, lbary_ltopo, vel_width=vel_width,\
                vel_step=0.1, spec_order=9, iv_order=10, sn_order=8,max_vel_rough=300)

        xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=3.0, Simple=True, W=W_ccf)
        pred = scipy.interpolate.splev(vels,tck1)
        xc_av /= pred
        #print 'min(vels),max(vels),len(vels), min(xc_av),max(xc_av),len(xc_av), refvel, moon_sig', min(vels),max(vels),len(vels), min(xc_av),max(xc_av),len(xc_av), refvel, moon_sig
        if max(np.abs(xc_av)) > 10:
            logger('Warn: stopped RV fit because of too big absolute value in xc_av: min(xc_av) = {0} , max(xc_av) = {1} , len(xc_av) = {2}'.format(min(xc_av),max(xc_av),len(xc_av)))
            return -999.0, 999.0, -999.0, 999.0
        # vels are in 100 m/s steps, see "vel_step=0.1" above, changing it to 0.01 changed the RVs only by <1.4m/s, average 0.14m/s (HARPS blue on CERES)
        p1,XCmodel,p1gau,XCmodelgau,Ls2 = \
                GLOBALutils.XC_Final_Fit( vels, xc_av, sigma_res = 4,\
                 horder=8, moonv=refvel, moons=moon_sig, moon=False)
        #print 'plgau 345', p1gau

        moonmatters = False
        if (know_moon and here_moon):
            moonmatters = True
            ismoon = True
            confused = False
            p1_m,XCmodel_m,p1gau_m,XCmodelgau_m,Ls2_m = GLOBALutils.XC_Final_Fit( vels, xc_av, \
            sigma_res = 4, horder=8, moonv = refvel, moons = moon_sig, moon = True)
            moon_flag = 1
        else:
            confused = False
            ismoon = False
            p1_m,XCmodel_m,p1gau_m,XCmodelgau_m,Ls2_m = p1,XCmodel,p1gau,XCmodelgau,Ls2
            moon_flag = 0

        bspan = GLOBALutils.calc_bss(vels,xc_av)
        SP = bspan[0]
        
        if (not known_sigma):
            disp = np.floor(p1gau[2])
            if (disp < 3.0): 
                disp = 3.0
            mask_hw_wide = av_m * disp / (GLOBALutils.Constants.c/1.0e3)
            ml_v = av_m - mask_hw_wide
            mh_v = av_m + mask_hw_wide            
            known_sigma = True
        else:
            cond = False
            
        if p1gau[2] > 1E3:
            cond = False
            
    xc_dict = {'vels':vels,'xc_av':xc_av,'XCmodelgau':XCmodelgau,'Ls2':Ls2,'refvel':refvel,\
           'rvels':rvels,'rxc_av':rxc_av,'rpred':rpred,'rxc_av_orig':rxc_av_orig,\
           'rvel0_xc':rvel0_xc,'xc_full':xc_full, 'p1':p1, 'sn':sn, 'p1gau':p1gau,\
           'p1_m':p1_m,'XCmodel_m':XCmodel_m,'p1gau_m':p1gau_m,'Ls2_m':Ls2_m,\
           'XCmodelgau_m':XCmodelgau_m}

    #moon_dict = {'moonmatters':moonmatters,'moon_state':moon_state,'moonsep':moonsep,\
    #     'lunation':lunation,'mephem':mephem,'texp':im_head['EXPTIME']}
    moon_dict = {'moonmatters':moonmatters,'moon_state':'dummy','moonsep':0,\
         'lunation':0,'mephem':mephem,'texp':0}

    #pkl_xc = params['path_rv_ceres'] + fsim.split(os.sep)[-1][:-4]+obname+'_XC_'+sp_type+'.pkl'
    pkl_xc = params['path_rv_ceres'] + fsim+'_XC_'+sp_type+'.pkl'      # without filename ending
    pickle.dump( xc_dict, open( pkl_xc, 'w' ) )

    #ccf_pdf = params['logging_path'] + fsim.split(os.sep)[-1][:-4] + obname + '_XCs_' + sp_type + '.pdf'       # dirout + 'logging/'
    #ccf_pdf = params['logging_path'] + fsim + '_XCs_' + sp_type + '.pdf'       # without filename ending
    ccf_pdf = params['path_rv_ceres'] + fsim + '_XCs_' + sp_type + '.pdf'             # without filename ending
    GLOBALutils.plot_CCF(xc_dict,moon_dict,path=ccf_pdf)

    #SNR_5130 = np.median(spec[8,28,1900:2101] ) This wavelength is not covered in exohspec
    #SNR_5130 = np.nanmedian(spec[8,0,1900:2101] ) # Set to order 20 artificially
    SNR_5130 = np.nan
    for cen_wave in [5130,6200,4000,7300,3000,8400,9500]:       # try different wavelengths in case one is not covered
        subspec, wave = get_spec(spec[8,:,:], spec[0,:,:], cen_wave, 100)
        if len(subspec) > 0:
            SNR_5130 = np.nanmedian(subspec)       # 5130 shows significant absorption lines
            break
        logger('Warn: The spectra around the wavelength of {0} Angstrom is not covered'.format(cen_wave))
    
    B,A = -0.00257864,0.07765779            # from Monte Carlo Simulation, different for each instrument (from CERES for HARPS)
    #B,A = 0.005, 0.2                       # Tested Ronny before 2/12/2019
    #print SNR_5130
    RVerr  =  B + ( 1.6 + 0.2 * p1gau[2] ) * A / np.round(SNR_5130)
    depth_fact = 1. + p1gau[0]/(p1gau[2]*np.sqrt(2*np.pi))
    if depth_fact < 0.6:
        depth_fact = 0.6
    depth_fact = (1 - 0.6) / (1 - depth_fact)
    RVerr *= depth_fact
    #print RVerr, depth_fact, p1gau, SNR_5130
    if RVerr < 0.002:
        RVerr = .002

    B,A = -0.00348879, 0.10220848           # from Monte Carlo Simulation, different for each instrument (from CERES for HARPS)
    BSerr = B + ( 1.6 + 0.2 * p1gau[2] ) * A / np.round(SNR_5130)
    if BSerr<0.002:
        BSerr = .002

    RV     = np.around(p1gau_m[1],4)  
    BS     = np.around(SP,4)   
    RVerr2 = np.around(RVerr,4)
    BSerr  = np.around(BSerr,4)
    
    return RV, RVerr2, BS, BSerr

def prepare_for_rv_packages(params):
    """
    Prepare files for TERRA and SERVAL and CERES
    
    """
    global calimages
    if np.max(calimages['wave_sol_dict_sci']['wavesol'][:,-1]) <= 100:
        return
    run_RV = True
    files_RV = []
    headers = dict()
    for file_RV in sorted(os.listdir(params['path_extraction'])):
        if not file_RV.endswith(".fits"):
            continue
        do_RV = True
        for no_RV_name in params['no_RV_names']:
            if file_RV.lower().find(no_RV_name) in [0,1,2,3,4,5]:
                do_RV = False
                break
        if not do_RV:
            continue
        im = fits.open(params['path_extraction']+file_RV)
        im_head = im[0].header
        headers[file_RV] = im_head      # for CERES
        files_RV.append(file_RV)        # for CERES
        spec = im[0].data
        if 'HiFLEx BCV' not in im_head.keys():          # BCV can be only calculated if coordinates are available -> only do RV on these stars
            continue
        obj_name = im_head['HiFLEx OBJNAME'].lower()
        obsdate_midexp = datetime.datetime.strptime(im_head['HiFLEx DATE-MID'],"%Y-%m-%dT%H:%M:%S.%f")
        if 'path_rv_terra' in params.keys():
            # CSV file for TERRA
            fname = params['path_rv_terra']+obj_name+'/data/'+obsdate_midexp.strftime('%Y-%m-%d%H%M%S')
            save_spec_csv(spec[params['dataset_rv_analysis'][0]], spec[0], spec[7], fname)        # spec[1]: Flux, spec[5]: Continuum corrected
        if 'path_rv_serval' in params.keys():
            # Store in a text file for serval
            numbers_levels = params['path_rv_serval'].count(os.sep, 2)  # start at 2 to not count './'
            add_text_to_file('../'*numbers_levels+params['path_extraction']+file_RV, 
                         params['path_rv_serval']+'filelist_{0}.txt'.format(obj_name), warn_missing_file=False )
    
    return files_RV, headers
                         
def run_terra_multicore(kwargs):
    [obj_name, params] = kwargs
    do_RV = True
    for no_RV_name in params['no_RV_names']:
        if obj_name.lower().find(no_RV_name) in [0,1,2,3,4,5]:
            do_RV = False
            break
    if not do_RV:
        return
    os.system('rm -f {0}{1}/results/synthetic.rv'.format(params['path_rv_terra'],obj_name) )     # Delete the old solution, as won't be created otherwise
    newfile = 'echo "0998     synthetic         LAB                LAB                    0.0          0.0       0.0       {0}/" > astrocatalog{0}.example'.format(obj_name)
    logger('For TERRA: creating a new astrocatalog.example with: '+newfile)
    os.system('rm -f astrocatalog{0}.example; '+newfile)
    cmd = 'java -jar {1} -ASTROCATALOG astrocatalog{2}.example 998 -INSTRUMENT CSV {0}'.format(calimages['wave_sol_dict_sci']['wavesol'].shape[0],params['terra_jar_file'], obj_name )
    log = 'logTERRA_{0}'.format(obj_name)
    logger('For TERRA: running TERRA: '+cmd)
    logger('Info: TERRA output and errors can be watched in {0}'.format(log))
    if True:
        with open(log, 'a') as logf:
            p = subprocess.Popen(cmd, stdin=subprocess.PIPE, shell=True, stdout=logf, stderr=subprocess.STDOUT)            # This doesn't wait until the procress has been finished
            p.communicate()                                                     # This makes it wait until it finished running
            log_returncode(p.returncode, 'Please check the logfiles in {0}. Problem occured for object {1}'.format(os.getcwd(), obj_name))
    else:
        logger('Warn: TERRA commented out')
    resultfile = '{2} {0}/results/synthetic.rv'.format(obj_name, params['path_rv_terra'], params['editor'])
    logger('For TERRA: results can be opened with: '+resultfile)
        
def run_terra_rvs(params):
    # Do the TERRA RVs
    if 'path_rv_terra' not in params.keys() or 'terra_jar_file' not in params.keys():
        logger('Warn: parameters "path_rv_terra" or "terra_jar_file" are not set.')
        return
    if not( os.path.isfile(params['terra_jar_file'])  and os.path.exists(params['path_rv_terra']) ):
        logger('Warn: TERRA is not installed or the path to the terra_jar_file wrong (currently: {0}), or the path for TERRA csv files does not exist ({1})'.format(params['terra_jar_file'], params['path_rv_terra']))
        return
        
    global calimages
    obj_names = []   
    logger('Info: Preparing for the TERRA analysis.')
    for root, dirs, files in os.walk(params['path_rv_terra'], followlinks=True):                       # Find all the objects again, as won't be added to obj_names when re-run
        for file in files:
            if file.endswith('.csv'):                       # has the file the correct ending?
                filename = os.path.join(root, file).replace(params['path_rv_terra'],'')                # Only relative folder and filename
                obj_name = filename.split(os.sep)[0]
                if obj_name not in obj_names:
                    obj_names.append(obj_name)
    logger('For TERRA: changing directory to '+params['path_rv_terra']+' . The steps to run TERRA are given in the logfile in that folder.')
    os.chdir(params['path_rv_terra'])
    logger('Info: All data logged in this file is relative to '+params['path_rv_terra'])
    
    kwargs = [[obj_name, params] for obj_name in obj_names]
    if params['use_cores'] > 1:
        logger('Info using multiprocessing on {0} cores'.format(params['use_cores']))
        p = multiprocessing.Pool(params['use_cores'])
        p.map(run_terra_multicore, kwargs)
        p.terminate()
    else:
        for kwarg in kwargs:
            run_terra_multicore(kwarg)
    os.chdir(params['path_run'])        # Go back to previous folder
    print('')
    logger('Info: Some errors reported by TERRA are expected (reading DRS ephemeris). The results are stored in {0}<object name>/results/synthetic.rv'.format(params['path_rv_terra']))

def run_serval_multicore(kwargs):           # create a procedure to run on multicore
    [obj_name, params] = kwargs
    do_RV = True
    for no_RV_name in params['no_RV_names']:
        if obj_name.lower().find(no_RV_name) in [0,1,2,3,4,5]:
            do_RV = False
            break
    if not do_RV or not os.path.isfile('filelist_{0}.txt'.format(obj_name)):
        return
        #continue
    cmd = 'ulimit -n 4096 ; {4}serval/src/serval.py {3} filelist_{3}.txt -inst HIFLEX -targrv 0 -pmin {0} -pmax {1} -oset {2} -safemode 2'.format(params['pmin'], 
                                    params['pmax'], params['oset'], obj_name, params['path_serval'])
    log = 'logSERVAL_{0}'.format(obj_name)
    logger('For SERVAL: running SERVAL: '+cmd)
    logger('Info: SERVAL output and errors can be watched in {0}'.format(log))
    if True:
        with open(log, 'a') as logf:
            p = subprocess.Popen(cmd, stdin=subprocess.PIPE, shell=True, stdout=logf, stderr=subprocess.STDOUT)    # shell=True is necessary
            p.communicate(input=os.linesep.encode())                                       # This waits until the process needs the enter
            log_returncode(p.returncode, 'Please check the logfiles in {0}. Problem occured for object {1}'.format(os.getcwd(), obj_name))
    else:
        logger('Warn: SERVAL commented out')
    resultfile = '{2} {0}/{0}.rvc.dat'.format(obj_name, params['path_rv_serval'], params['editor'])  
    logger('For SERVAL: results can be opened with: '+resultfile)

def run_serval_rvs(params,):
    # Do the SERVAL RVs
    if 'path_serval' not in params.keys() or 'path_rv_serval' not in params.keys():
        logger('Warn: parameters "path_serval" or "path_rv_serval" are not set.')
        return
    if not( os.path.exists(params['path_serval']+'serval') and os.path.exists(params['path_serval']+'python') and os.path.exists(params['path_rv_serval']) ):
        logger('Warn: SERVAL is not installed or path_serval is wrong (currently: {0}), or the folder for the "filelist_* are missing in {1}'.format(params['path_serval'], params['path_rv_serval']))
        return
    if sys.version_info[0] > 2:
        logger('Warn: This is a python3 environment, however, SERVAL requires a python2 environment')
        return
    global calimages
    obj_names = []
    logger('Info: Preparing for the SERVAL analysis.')
    for root, dirs, files in os.walk(params['path_rv_serval'], followlinks=True):                       # Find all the objects again, as won't be added to obj_names when re-run
        for file in files:
            if file.endswith('.txt') and file.find('filelist_') == 0:                       # has the filename the correct format?
                obj_name = file.split('filelist_')[-1].split('.txt')[0]
                if obj_name not in obj_names:
                    obj_names.append(obj_name)                   
    hiflex_file = params['path_rv_serval']+'conf_hiflex_serval.txt'
    servalparams = dict()
    if os.path.isfile(hiflex_file):
        try:
            keys, values = np.genfromtxt(hiflex_file, dtype=str, comments='#', delimiter='=', filling_values='', autostrip=True, unpack=True)
            servalparams = dict(zip(keys, values))
        except ValueError as error:
            print('Problems when reading {0}. Will create a new file with standard settings'.format(hiflex_file))
    else:
        print('Warn: No {0} found.'.format(hiflex_file))
    for key in servalparams.keys():
        try:
            servalparams[key] = int(servalparams[key])
        except:
            del servalparams[key]
    servalparams['orders'] = servalparams.get('orders', calimages['wave_sol_dict_sci']['wavesol'].shape[0])
    servalparams['data_dataset'] = servalparams.get('data_dataset', params['dataset_rv_analysis'][0])
    servalparams['error_dataset'] = servalparams.get('error_dataset', params['dataset_rv_analysis'][1])
    servalparams['wave_dataset'] = servalparams.get('wave_dataset', 9)
    servalparams['mask_dataset'] = servalparams.get('mask_dataset', 7)
    servalparams['order_snr'] = servalparams.get('order_snr', int(calimages['wave_sol_dict_sci']['wavesol'].shape[0]/2) )
    with open(hiflex_file, 'w') as file:
        file.write('# This file is used to control which data is used in SERVAL. It is read when inst_HIFLEX is used by SERVAL.'+os.linesep)
        file.write('orders = {1}{0}data_dataset = {2}{0}error_dataset = {3}{0}wave_dataset = {4}{0}mask_dataset = {5}{0}order_snr = {6}{0}'.format(os.linesep,
                            servalparams['orders'], servalparams['data_dataset'], servalparams['error_dataset'], servalparams['wave_dataset'], 
                            servalparams['mask_dataset'], servalparams['order_snr'] ))

    # Set a path if necessary
    pypath = ''
    if 'PYTHONPATH' in os.environ.keys():
        pypath = os.environ["PYTHONPATH"] + os.pathsep      # Probably should exclude all python3 paths
    pypath += params['path_serval']+'python'
    logger('For SERVAL: set variable: bash: export PYTHONPATH={0} ; csh: setenv PYTHONPATH {0}'.format(pypath))
    os.environ["PYTHONPATH"] = pypath
    xlows, xhighs = calimages['sci_trace'][1:3]
    xran = np.max(xhighs) - np.min(xlows)
    exclude_x = 0.1         # 0.1: 10% of outermost pixel on each side are excluded
    params['pmin'] = int(np.min(xlows) + 0.1*xran)
    params['pmax'] = int(np.max(xhighs) - 0.1*xran)
    params['oset'] = '{0}:{1}'.format(0,calimages['wave_sol_dict_sci']['wavesol'].shape[0])     # if shape[0] is 5, then oset will lead to 0,1,2,3,4 being used
    logger('For SERVAL: changing directory to '+params['path_rv_serval']+' . The steps to run SERVAL are given in the logfile in that folder.')
    os.chdir(params['path_rv_serval'])
    logger('Info: All data logged in this file is relative to '+params['path_rv_serval'])
    kwargs = [[obj_name, params] for obj_name in obj_names]
    #with multiprocessing.Pool(params['use_cores']) as p:       # only possible in python3
    #    p.map(run_serval, obj_names)
    if params['use_cores'] > 1:
        logger('Info using multiprocessing on {0} cores'.format(params['use_cores']))
        p = multiprocessing.Pool(params['use_cores'])
        p.map(run_serval_multicore, kwargs)
        p.terminate()
    else:
        for kwarg in kwargs:
            run_serval_multicore(kwarg)
    os.chdir(params['path_run'])        # Go back to previous folder
    print('')
    logger(('Info: Finished the SERVAL analysis. Some errors reported by SERVAL are expected.'+\
          ' The results are stored in {0}<object name>/<object name>.rvc.dat.'+\
          ' If serval failed (the result file is missing), run it again using less orders by setting oset to a smaller range (especially orders with low SN).'+\
          ' You can also modify the parameters in {0}conf_hiflex_serval.txt , e.g. to select a different dataset or a different order to measure the SNR.'+\
          ' The command history can be found in {0}cmdhistory.txt. Before running serval: cd {0}').format(params['path_rv_serval']))

def run_ceres_multicore(kwargs):
    [file_RV, params, headers, force_rvs, stellar_par] = kwargs
    #params = kwargs['params']
    #stellar_par = kwargs['stellar_par']
    #spec = kwargs['spec']
    # = kwargs['']
    # = kwargs['']
    # = kwargs['']
    if file_RV not in headers.keys():
        return
    im_head = headers[file_RV]
    if ( 'HiFLEx CERES RV_BARY' in im_head.keys() or 'HiFLEx CERES RV' in im_head.keys() ) and not force_rvs:
        return
    if 'HiFLEx OBJNAME' not in im_head.keys():
        return
    obname = im_head['HiFLEx OBJNAME']
    im = fits.open(params['path_extraction']+file_RV)
    spec = np.array(im[0].data, dtype=np.float64)
    spec[0,:,:] = wavelength_air_to_vacuum(spec[9,:,:])     # Vaccuum wavelength from the wavelengths not corrected by barycentric velocity
    spec = clip_noise(spec)
    if obname not in stellar_par.keys():            # hardcoded values
        T_eff, logg, Z, vsini, vel0 = 5777, 4.4374, 0.0134, 2, 0
    else:
        [T_eff, logg, Z, vsini, vel0] = list(stellar_par[obname])
    gobs = ephem.Observer()  
    gobs.name = copy.copy(params['site'])
    gobs.lat  = math.radians(params['latitude'])     # lat/long in decimal degrees  
    gobs.long = math.radians(params['longitude'])
    gobs.elevation = params['altitude']
    gobs.date = im_head['HiFLEx DATE-MID'].replace('T',' ')  # to make into ('%Y-%m-%d %H:%M:%S')
    mephem    = ephem.Moon()
    mephem.compute(gobs)
    
    # Change the path to the object_file to the result_path, if necessary
    object_file = params['object_file']         # needed later
    object_files, object_files_full = find_file_in_allfolders(object_file, [params['result_path']] + params['raw_data_paths'])      # Check also the result and raw data paths for object names
    ra2, dec2, pmra, pmdec, epoch = 0., 0., 0., 0., 2000.
    for object_file in object_files:        # Will be run in any case, even if empty
        ra2, dec2, epoch, pmra, pmdec, dummy, dummy = getcoords_from_file([obname], 0, filen=object_file, warn_notfound=False)        # mjd=0 because because not using CERES to calculated BCV
        if ra2 !=0 or dec2 != 0:                                           # Found the objects -> obnames becomes the single entry which matched entry in params['object_file']
            break
    #if ra2 ==0 and dec2 == 0:
    #    'Was not found, but no problem'
    RV, RVerr2, BS, BSerr = rv_analysis_ceres(params, spec, file_RV.replace('.fits',''), obname, object_file, mephem, vsini)
    bjdtdb = im_head.get('HiFLEx BJDTDB', np.nan)
    bcvel_baryc = im_head.get('HiFLEx BCV', np.nan)
    if np.isnan(bcvel_baryc):
        RV_BARY = round(RV,4)
        im_head['HIERARCH HiFLEx CERES RV'] = (RV_BARY, 'RV in km/s (without BCV)')
    else:
        RV_BARY = round(RV+bcvel_baryc,4)
        im_head['HIERARCH HiFLEx CERES RV'] = (round(RV,4), 'RV in km/s (without BCV)')
        im_head['HIERARCH HiFLEx CERES RV_BARY'] = (RV_BARY, 'barycentric RV in km/s (measured RV+BCV)')
    # Air to Vacuum wavelength difference only causes < 5 m/s variation: https://www.as.utexas.edu/~hebe/apogee/docs/air_vacuum.pdf (no, causes 83.15 km/s shift, < 5m/s is between the different models)
    logger(('Info: The radial velocity (including barycentric correction) for {0} at [ {6} , JD= {8} , BJDTDB= {9} (center)] gives: '+\
            'RV = {1} +- {2} km/s, Barycentric velocity = {5} km/s, and BS = {3} +- {4} km/s. '+\
            'The total extracted flux is {7} counts.').format(file_RV, RV_BARY, round(RVerr2,4), round(BS,4), 
                    round(BSerr,4), bcvel_baryc, im_head['HiFLEx DATE-MID'], np.nansum(spec[1,:,:]), 
                    im_head['HiFLEx JD'], bjdtdb ))
    im_head['HIERARCH HiFLEx CERES RV_ERR']  = (round(RVerr2,4), 'Uncertainty RV in km/s')
    im_head['HIERARCH HiFLEx CERES BS']      = (round(BS,4), 'Bisector')
    im_head['HIERARCH HiFLEx CERES BS_ERR']  = (round(BSerr,4), 'Uncertainty BS')
    
    save_multispec(im[0].data, params['path_extraction']+file_RV, im_head)

def run_ceres_rvs(params, files_RV, headers):
    # Do the CERES RVs
    force_stellar_pars = False
    force_rvs = False
    if not( os.path.exists(params['path_ceres']+'utils/Correlation') and os.path.exists(params['path_ceres']+'utils/GLOBALutils') \
            and os.path.exists(params['path_ceres']+'utils/OptExtract') and os.path.exists(params['path_ceres']+'utils/CCF') ):         # if necessary files exist
        logger('Warn: CERES RV analysis did not run. If this a mistake, please check that it is installed under {0} and includes the folders: utils/Correlation, utils/GLOBALutils, utils/OptExtract, utils/CCF'.format(params['path_ceres']))
        return
    if sys.version_info[0] > 2:
        logger('Warn: This is a python3 environment, however, CERES requires a python2 environment')
        return    
    base = params['path_ceres']
    models_path = base + 'data/COELHO_MODELS/R_40000b/'
    load_ceres_modules(params)
    
    stellar_par = dict()
    for file_RV in files_RV:
        im_head = headers[file_RV]
        if 'HiFLEx OBJNAME' not in im_head.keys():
            continue
        T_eff =  im_head.get('HiFLEx CERES Teff', np.nan)
        logg  =  im_head.get('HiFLEx CERES logg', np.nan)
        Z     =  im_head.get('HiFLEx CERES Z', np.nan)
        vsini =  im_head.get('HiFLEx CERES vsini', np.nan)
        vel0  =  im_head.get('HiFLEx CERES vel0', np.nan)
        hardcoded =  np.isnan([T_eff, logg, Z, vsini, vel0]).any()      # If one is np.nan -> measure again
        obname = im_head['HiFLEx OBJNAME']
        if (hardcoded or force_stellar_pars) and os.path.exists(models_path):       # Don't use the header keywords
            im = fits.open(params['path_extraction']+file_RV)
            spec = np.array(im[0].data, dtype=np.float64)
            spec[0,:,:] = wavelength_air_to_vacuum(spec[9,:,:])     # Vaccuum wavelength
            spec = clip_noise(spec)
            T_eff, logg, Z, vsini, vel0, hardcoded = stellar_params_ceres(params, spec, file_RV.replace('.fits',''), calimages['wave_sol_dict_sci']['wavesol'])
            if not hardcoded:
                im_head['HIERARCH HiFLEx CERES Teff']  = (T_eff, 'Teff measured by CERES CCF')
                im_head['HIERARCH HiFLEx CERES logg']  = (logg,  'logg measured by CERES CCF')
                im_head['HIERARCH HiFLEx CERES Z']     = (Z,     'Z measured by CERES CCF')
                im_head['HIERARCH HiFLEx CERES vsini'] = (vsini, 'vsini measured by CERES CCF')
                im_head['HIERARCH HiFLEx CERES vel0']  = (vel0,  'vel0 measured by CERES CCF')
                headers[file_RV] = im_head
        else:
            if obname in stellar_par.keys():
                stellar_par[obname] = np.vstack(( stellar_par[obname], [T_eff, logg, Z, vsini, vel0] ))
            else:
                stellar_par[obname] = np.array([T_eff, logg, Z, vsini, vel0])
                stellar_par[obname].shape = (1,5)       # Otherwise the median later will be a problem

    for obname in stellar_par.keys():
        nspec = stellar_par[obname].shape[0]
        ori = copy.deepcopy(stellar_par[obname])
        stellar_par[obname] = np.round(np.median(stellar_par[obname], axis=0), 2)
        if nspec >= 2:
            stdev = np.round(np.std(ori, axis=0, ddof=1), 2)
            text = 'Teff = {0} +- {5} , logg = {1} +- {6} , Z = {2} +- {7} , vsini = {3} +- {8} , vel0 = {4} +- {9}'.format(stellar_par[obname][0], 
                     stellar_par[obname][1], stellar_par[obname][2], stellar_par[obname][3], stellar_par[obname][4], stdev[0], stdev[1], stdev[2], stdev[3], stdev[4])
        else:
            text = 'Teff = {0} , logg = {1} , Z = {2} , vsini = {3} , vel0 = {4}'.format(stellar_par[obname][0],
                     stellar_par[obname][1], stellar_par[obname][2], stellar_par[obname][3], stellar_par[obname][4])
        logger('Info: With the CERES routine the following parameters have been determined for {0} using {1} spectra: {2}'.format(obname, nspec, text))

    kwargs = [[file_RV, params, headers, force_rvs, stellar_par] for file_RV in files_RV]
    if params['use_cores'] > 1 and not params['started_from_p3']:
        logger('Info using multiprocessing on {0} cores'.format(params['use_cores']))
        p = multiprocessing.Pool(params['use_cores'])
        p.map(run_ceres_multicore, kwargs)      # this can get stuck in "n_Edlen" when the python2 environment is started by the pipeline in python3
        p.terminate()
        p.join()
    else:
        for kwarg in kwargs:
            run_ceres_multicore(kwarg)

     






