#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
import astropy.time as asttime
import astropy.coordinates as astcoords
import astropy.units as astunits
import astropy.constants as astconst
import matplotlib    # To avoid crashing when ssh into Narit using putty
if not 'DISPLAY' in os.environ:
    matplotlib.use('agg')    # To avoid crashing when ssh into Narit using putty, however this means plots are not shown (test when working in front of uhppc30)
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import time
import datetime
import operator
import copy
import random
import warnings
from scipy.optimize import curve_fit
import scipy.interpolate as inter
import scipy
import scipy.signal
from tqdm import tqdm
import tkcanvas as tkc
import json
# detect python version
if sys.version_info[0] < 3:
    import Tkinter as Tk
    from collections import OrderedDict as dict
else:
    import tkinter as Tk
import plot_img_spec
import psutil
import ephem
from math import radians as rad
import barycorrpy
import glob

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
oldorder = 12

class Constants:                    # from Ceres
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
    file = open(logfile, 'a')
    file.write('{0} - {1} - {2}\n'.format( time.strftime("%Y%m%d%H%M%S", time.localtime()), os.getpid(), message ))
    if printarrayformat != [] and printarray != []:
        for line in printarray:
            text = ''
            for i,printformat in enumerate(printarrayformat):
                #print printformat,line[i]
                text += printformat%line[i] + '\t'
            file.write(text[:-1]+'\n')
            if show:
                print(text)
    file.close()
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
    list_of_files = glob.glob('{0}/*.py'.format(os.path.realpath(__file__).rsplit('/',1)[0]))
    list_of_files = sorted(list_of_files, key=os.path.getmtime, reverse=True)
    for fname in list_of_files:
        filedata1 = read_text_file(fname, no_empty_lines=False)
        filedata2 = read_text_file(fname, no_empty_lines=True)
        text += '\n    {1}  {2}  {3} {4}  {0}'.format( fname, 
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
    # Extra steps to check that the user didn't make mistakes in the file:
    data = read_text_file(textfile, no_empty_lines=True)
    data = convert_readfile(data, [str,str], delimiter='=', replaces=[' '], ignorelines=['#'])
    for line in data:
        if len(line) != 2:
            logger(('Error: The line in file {0} containing the entries {1} has the wrong format. Expected was "parameter = value(s)" .'+\
                    'Please check if the "=" sign is there or if a comment "#" is missing.').format(textfile, line))
    data = np.genfromtxt(textfile, dtype=str, comments='#', delimiter='=')
    if data.shape[0] == 0 or len(data.shape) < 2:
        logger('Error: No values found when reading the parameterfile {0}.'.format(textfile))
    keys, values = data[:, 0], data[:, 1]
    for k in range(len(keys)):
        keys[k] = keys[k].replace(' ', '')
        values[k] = values[k].replace(' ', '')
    textparams = dict(zip(keys, values))
    
    return textparams

def textfileargs(params, textfile=None):
    """
    Imports parameters from text file located at textfile, parameters given in the command line and the parameters given with the programm code
    Parameters in the programm code will be overwritten by parameters in the text file and parameters in the text file will be overwritten by parameters from the command line
    :param params: dictionary, default parameter dictionary (the default values)
    :param textfile: string, location of config file
    :return params: dictionary with added parameters and updated default parameters 
    
    Please note, the parameters are only lists or single entries. No numpy arrays allowed
    """
    
    # Set up standard parameters
    params['in_shift'] = -0
    
    emsg = 'Error in config file: '
    
    textparams = read_parameterfile(textfile)
    params.update(textparams)
    
    if 'configfile_fitsfiles' in params.keys():
        if os.path.exists(params['configfile_fitsfiles']):
            textparams = read_parameterfile(params['configfile_fitsfiles'])
            params.update(textparams)
    
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
        elif arg in ['nocheck', 'prepare']:
            continue        # Do nothing, just prevent the warning below
        else:
            logger('Warn: I dont know how to handle command line argument: {0}'.format(arg))
    params.update(cmdparams)
    list_txt = ['use_catalog_lines', 'raw_data_file_endings', 'raw_data_mid_exposure_keys']
    list_int = ['subframe', 'arcshift_range', 'order_offset', 'px_offset', 'px_offset_order', 'polynom_order_traces', 'polynom_order_intertraces',
             'bin_search_apertures', 'bin_adjust_apertures', 'polynom_bck']
    list_float = ['opt_px_range', 'background_width_multiplier']
    list_abs = ['arcshift_range']
    ints = ['polynom_order_apertures', 'rotate_frame']
    floats = ['max_good_value', 'catalog_file_wavelength_muliplier', 'extraction_width_multiplier', 'arcextraction_width_multiplier',
              'resolution_offset_pct', 'diff_pxs', 'maxshift', 'wavelength_scale_resolution', 'width_percentile', 'raw_data_timezone_cor',
              'altitude', 'latitude', 'longitude', 'in_shift']
    bools = ['flip_frame', 'update_widths', 'GUI']
    sides = ['arcshift_side']
    results = ['path_extraction', 'path_extraction_single', 'logging_path', 'path_reduced', 'path_rv_ceres', 'path_csv_terra', 
               'path_harpsformat', 'object_file']     #, 'configfile_fitsfiles' (excluded, as should be handled as conf.txt
    full_filenames = ['badpx_mask_filename', 'original_master_traces_filename', 'original_master_wavelensolution_filename', 'reference_catalog',
                'configfile_fitsfiles', 'raw_data_file_list', 'terra_jar_file']    # deal with full filenames -> nothing to do
    texts = ['editor', 'extracted_bitpix', 'site']      # -> nothing to do
    list_raw, paths, loggings = [], [], []
    for entry in params.keys():                         # make the live easier by adding entries automatically to the lists above
        if entry not in list_txt and (entry.find('_rawfiles') > 0 or entry.find('calibs_') > -1):           # add to the list
            list_txt.append(entry)
        if entry not in list_raw and (entry.find('_rawfiles') > 0):
            list_raw.append(entry)
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
    all_parameters = list(np.unique( list_txt + list_int + list_float + list_abs + ints + floats + bools + sides + results + full_filenames + list_raw + paths + loggings + texts ))
    trues = ['yes', 'true', '1']
    falses = ['no', 'false', '0']
    lefts = ['left','l']
    rights = ['right','r']
    centers = ['center','c']
    undeclared_params = ''
    # Important settings first, as other things depend on them:
    for entry in paths:
        if entry in params.keys():
            params[entry] = (params[entry]+'/').replace('//', '/')      # Add a / at the end in case the user didn't
    for entry in results:
        if entry in params.keys():                                            # deal with result filenames/folders -> add result_path
            params[entry] = params['result_path'] + params[entry]
    # All parameters
    for entry in params.keys():
        if entry in list_txt+list_int+list_float:
            temp = params[entry]
            for i in ['[', ' ', ']']:
                temp = temp.replace(i,'')
            if len(temp) == 0:
                temp = []
            else:
                temp = temp .split(',')
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
            if params[entry].lower() in trues:
                params[entry] = True
            elif params[entry].lower() in falses:
                params[entry] = False
            else:
                logger(emsg + 'Parameter "{0}" (value of "{1}") must be any of the following: {2}'.format(entry, params[entry], trues+falses))
        if entry in sides:
            if params[entry].lower() in lefts:
                params[entry] = -1
            elif params[entry].lower() in rights:
                params[entry] = 1
            elif params[entry].lower() in centers:
                params[entry] = 0
            else:
                logger(emsg + 'Parameter "{0}" (value of "{1}") must be one of the following: {2}'.format(entr, params[entry], lefts+rights+centers))
        if entry in paths:
            if not os.path.exists(params[entry]) and entry not in ['raw_data_path', 'path_ceres', 'terra_jar_file'] and params[entry].lower() not in ['na/', params['result_path']+'na/', params['result_path']+'/na/']:
                try:                                                    # Create folders, if necessary
                    os.makedirs(params[entry])
                except:
                    logger('Warn: Cannot create directory {0}'.format(params[entry]))
        if entry in list_raw:                                           # deal with lists of raw data filenames -> add path
            for i in range(len(params[entry])):
                params[entry][i] = params['raw_data_path'] + params[entry][i]
        if entry in loggings:                                           # deal with logging filenames/folders -> add logging_path
            params[entry] = params['logging_path'] + params[entry]
        if entry not in all_parameters:
            undeclared_params += '{0}, '.format(entry)

    overdeclared_params = ''
    for entry in all_parameters:
        if entry not in params.keys():
            overdeclared_params += '{0}, '.format(entry)
    if len(overdeclared_params) > 0:
        logger('Warn: The following parameters are expected, but do not appear in the configuration files. : {0}\n\t!!! This is likely to cause problems later !!!'.format(overdeclared_params[:-2]))
    if len(undeclared_params) > 0:
        logger('Warn: The following parameters appear in the configuration files, but the programm did not expect them: {0}'.format(undeclared_params[:-2]))
        
    return params

def update_calibration_memory(key,value):
    """
    Add new information to the global variable calimages, which is not accessable from other python files
    :param key: string, key for the dictionary
    :param value: string, number, array, or anything: value for the dictionary
    """
    global calimages
    calimages[key] = copy.deepcopy(value)

def read_text_file(filename, no_empty_lines=False):
    """
    Read a textfile and put it into a list, one entry per line
    Empty lines means they consist only of tabs and spaces
    """
    text = []
    if os.path.isfile(filename) == True:
        text1 = open(filename,'r').readlines()
        for line in text1:
            line = line.replace('\n', '')
            linetemp = line.replace('\t', '')
            linetemp = linetemp.replace(' ', '')
            if ( line == '' or linetemp == '') and no_empty_lines:
                continue
            text.append(line)
    else:
        logger('Warn: File {0} does not exist, assuming empty file'.format(filename))    
    return text

def add_text_file(text, filename):
    """
    Adds a line or lines of text to a file
    """
    file = open(filename,'a')
    file.write(text+'\n')
    file.close()

def add_text_to_file(text, filename):
    """
    If the text is not yet in the file, the test is added
    """
    oldtext = read_text_file(filename)
    exists = False
    for line in oldtext:
        if line.find(text) != -1:
            exists = True
            break
    if exists == False:
        add_text_file(text, filename)
        
def convert_readfile(input_list, textformats, delimiter='\t', replaces=[], ignorelines=[]):
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
    :param replaces: 1d list or array of strings, contains the strings which should be replaced by ''
    :param ignorelines: List of strings and/or lists. A list within ignorelines needs to consist of a string and the maximum position this string can be found.
                        If the string is found before the position, the entry of the input list is ignored
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
        if len(entry) < len(textformats):
            continue
        for index in range(len(textformats)):
            for textformat in textformats[index]:
                if type(textformat) == type:
                    if type(entry[index]) == datetime.datetime:
                        epoch = datetime.datetime.utcfromtimestamp(0)
                        entry[index] = (entry[index] - epoch).total_seconds()         # (obsdate - epoch) is a timedelta
                    entry[index] = textformat(entry[index])
                elif type(textformat) == str:
                    try:
                        entry[index] = datetime.datetime.strptime(entry[index], textformat)
                    except:
                        logger('Error: Cannot convert {0} into a datetime object using the format {1}'.format(entry[index], textformat))
                else:
                    logger('Error: Programming error: textformats which where hand over to convert_readfile are wrong. It has to be a type or a string')
        result_list.append(entry)
    return result_list

def read_badpx_mask(params): 
    """
    Reads the bad-pixel-mask and applies the corrections to it
    :param params: Dictionary with all the parameters: 
                badpx_mask_filename is used in order to get the bad px mask. If it doesn't exist, than all pixels are fine
                calibs is used in order to if a subframe of the image is used. In this case
                subframe is used to find the area used
    :return badpx_mask: 2d numpy array of the bad pixel mask
    """
    filename = params['badpx_mask_filename']
    subframe = params['subframe']
    if os.path.isfile(filename) == True:
        badpx_mask = np.array(fits.getdata(filename), dtype=np.float64)
        badpx_mask = rotate_flip_frame(badpx_mask, params )
        if subframe != [] and badpx_mask.shape != (subframe[0],subframe[1]):
            badpx_mask = badpx_mask[subframe[2]: subframe[0]+subframe[2], subframe[3]: subframe[1]+subframe[3]]
        logger('Info: badpixel mask loaded: {0}'.format(filename))
    else:
        logger('Warn: No bad pixelmask found ({0}). Assuming no bad pixels'.format(filename))
        badpx_mask = np.ones((subframe[0],subframe[1]))
        
    return badpx_mask

def read_background(params, filename): 
    """
    Reads the background map and applies the corrections to it
    :param params: Dictionary with all the parameters
    :param filename: path and name of the background file
    :return background: 2d numpy array of the background, normalised to 1s exposure time
    """
    imtype = 'background'
    params['calibs'] = params['calibs_read']
    if '{0}_calibs_read'.format(imtype) in params.keys():
        params['calibs'] = params['{0}_calibs_read'.format(imtype)]
    im, im_head = read_file_calibration(params, filename)
    exptime = im_head[params['raw_data_exptim_keyword']]
    background = im/exptime
    
    return background

def warn_images_not_same(ims, names):
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

def read_file_calibration(params, filename, level=0):
    """
    Reads the filename and applies the calibrations as given in params['calibs']
    This works also if the files are stored with the following header: BITPIX = 16, BZERO = 32768
    :param params: Dictionary with all the parameters
    :param filename: path and name of the file
    :return: 2D numpy array of the file, and the header of the file read
    """
    global calimages
    if os.path.isfile(filename) == False:
        logger('Error: File {0} is missing.'.format(filename))
    im = fits.open(filename)
    #print im.info()
    im_head = im[0].header
    im = np.array(im[0].data, dtype=np.float64)
    ims = im.shape
    while len(ims) > 2:
        if ims[0] == 1:                 # e.g. MRES from NARIT
            im = im[0,:,:]
        elif ims[-1] == 1:
            im = im[:,:,0]
        else:
            logger(('Error: The file is stored in a multi-demensional array, which I do not know how to handle. The size of the image is {0}. '+\
                    'This requires a small adjustment to the code in procedure read_file_calibration.').format(ims))
        ims = im.shape
    exptime = im_head[params['raw_data_exptim_keyword']]        #saved as float
    logger('Info: image loaded: {0}'.format(filename))
    im = rotate_flip_frame(im, params )
    for entry in params['calibs']:
        if entry == '':
            continue
        logtxt, headtxt = [], []
        if entry.lower().find('subframe') > -1:
            subframe = params[entry]
            if im.shape != (subframe[0],subframe[1]):                   # only apply subframe if the file doesn't have the size already
                im = im[subframe[2]: subframe[0]+subframe[2], subframe[3]: subframe[1]+subframe[3]]
            logger('Info: {1}: subframe applied: {0} ({2})'.format(entry, level, subframe))
            im_head['HIERARCH EXO_PIPE redu{0}a'.format(level)] = 'Subframe: {0}'.format(entry)
        elif entry.lower().find('bias') > -1:
            if entry not in calimages:
                 create_image_general(params, entry, level=level+1)
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
            warn_images_not_same([im, calimages[entry]], [filename,entry])
            if np.percentile(calimages[entry], 90) > 2000 or np.percentile(calimages[entry], 10) < -100:
                logger('Warn: The dark ({0}) has unphysical values: 10%-percentile = {1}, 90%-percentile = {2}'.format(entry, np.percentile(calimages[entry], 10), np.percentile(calimages[entry], 90)))
            im = im - calimages[entry]
            logtxt, headtxt = ['dark correction applied'], ['redu{0}c'.format(level), 'Dark']
        elif entry.lower().find('rflat') > -1:
            if entry not in calimages:
                 create_image_general(params, entry, level=level+1)
            warn_images_not_same([im, calimages[entry]], [filename,entry])
            im = im / (calimages[entry] / np.median(calimages[entry]) )
            logtxt, headtxt = ['flat correction applied with normalised flat (rflat)'], ['redu{0}d'.format(level), 'Flat']
        elif entry.lower().find('background') > -1 and not entry.lower().find('localbackground') > -1:
            if entry not in calimages:
                 calimages[entry] = read_background(params, params[entry+'filename'])
            warn_images_not_same([im, calimages[entry]], [filename,params[entry+'filename']])
            im = im - calimages[entry]*exptime
            logger('Info: {1}: background correction applied: {0}'.format(params[entry], level))
            im_head['HIERARCH EXO_PIPE redu{0}e'.format(level)] = 'Background: {0}'.format(params[entry])
        elif entry.lower().find('badpx_mask') > -1:
            if entry not in calimages:
                calimages[entry] = read_badpx_mask(params)
            warn_images_not_same([im, calimages[entry]], [filename,entry])
            im = im * calimages[entry]
            ims =  im.shape
            nonzeroind = np.nonzero(1-calimages[entry])       #find the indexes where the badpx mask is zero
            for i in range(len(nonzeroind[0])):         #find a suitable amound of surrounding pixels
                for j in range(1,5):
                    section = im[max(0,nonzeroind[0][i]-j) : min(ims[0],nonzeroind[0][i]+j+1) , max(0,nonzeroind[1][i]-j) : min(ims[1],nonzeroind[1][i]+j+1)]
                    section = section[section!=0]           #only use areas <> 0
                    if len(section) >= 2:
                        break
                if len(section) == 0:
                    logger('Warn: cannot replace bad pixel ({0}, {1}) with surrounding area in {2}'.format(nonzeroind[0][i],nonzeroind[1][i],filename))
                else:
                    im[nonzeroind[0][i],nonzeroind[1][i]] = np.median(section)  #replace bad px with the median of each surrounding area
            logger('Info: {1}: badpx correction applied: {0}'.format(entry, level))
            im_head['HIERARCH EXO_PIPE redu{0}f'.format(level)] = 'Bad-pixel-mask: {0}'.format(entry)
        elif entry.lower().find('localbackground') > -1:
            if 'sci_trace' in calimages.keys() and 'cal_trace' in calimages.keys():
                logger('Step: Performing the background fit')
                sci_tr_poly, xlows, xhighs, widths = calimages['sci_trace']
                cal_tr_poly, axlows, axhighs, awidths = calimages['cal_trace']
                bck_px_sci = find_bck_px(im, sci_tr_poly, xlows, xhighs, widths, params['background_width_multiplier'][0])
                bck_px_cal = find_bck_px(im, cal_tr_poly, axlows, axhighs, awidths, params['background_width_multiplier'][1])
                bck_px = bck_px_sci * bck_px_cal        # mask with 0 for all the traces and 1 for background data
                bad_values = ( im*bck_px > np.percentile(im[bck_px==1],95) )        # exclude brightest data in the background data (emission lines or borders of the traces
                bck_px[bad_values] = 0
                
                # Some deviation at the red side with not many lines, computational heavy
                bck_im = fit_2d_image(im, params['polynom_bck'][1], params['polynom_bck'][0], w=bck_px)
                #plot_img_spec.plot_image(im, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'orig')
                #plot_img_spec.plot_image(bck_px, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'bck_px')
                #plot_img_spec.plot_image(bck_im, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'bck_im')
                #plot_img_spec.plot_image(im-bck_im, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'diff')
                plot_img_spec.plot_image((im - bck_im)*bck_px, \
                                        [params['logging_path']+'background_subtracted-'+filename.rsplit('/',1)[-1]+'.png'],\
                                         1, False, [0.05,0.95,0.95,0.05], 'difference between image and background fit')
                im = im - bck_im
                bck_noise_std, bck_noise_var = measure_background_noise(im * bck_px)
                if np.isnan(bck_noise_var):
                    bck_noise_var = -1
                logger('Info: {1}: background correction applied: {0}'.format(entry, level))
                im_head['HIERARCH EXO_PIPE redu{0}f'.format(level)] = 'Background: {0}'.format(entry)
                im_head['HIERARCH EXO_PIPE BCKNOISE'] = round(bck_noise_std,8)
                im_head['HIERARCH EXO_PIPE BNOISVAR'] = (round(bck_noise_var,8), 'Variation of noise through image')             # Background noise variation can be very high, because some light of the traces remains
            else:
                logger('Warn: Could not apply the calibration step {0} because the science and/or calibration traces are not yet known.'.format(entry))
        elif entry.lower().find('combine_sum') > -1 or entry.lower().find('combine_mean') > -1 or entry.lower().find('normalise') > -1:
            'nothing to do, as for a different step'
        else:
            logger('Warn: do not know what to do with this correction: {0}'.format(entry))
        if logtxt != [] and headtxt != []:
            im_median, im_std = int(round(np.median(calimages[entry]))), int(round(np.std(calimages[entry], ddof=1)))
            logger('Info: {4}: {3}: {0} (median={1}, std={2})'.format(entry, im_median, im_std, logtxt[0], level))
            im_head['HIERARCH EXO_PIPE '+headtxt[0]] = '{3}: {0}, median={1}, std={2}'.format(entry, im_median, im_std, headtxt[1])
    #logger('Info: image loaded and processed: {0}'.format(filename))
    if os.path.exists(params['path_reduced']) and params['path_reduced'].lower() != 'na/':       # Save the reduced image
        fname = filename.rsplit('/',1)
        save_im_fits(params, im, im_head,  params['path_reduced']+fname[-1])
    return im, im_head

def create_image_general(params, imtype, level=0):
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
            logger('Info: Using existing {0}: {1}'.format(imtype,params['master_{0}_filename'.format(imtype)]))
            params['calibs'] = params['calibs_read']
            if '{0}_calibs_read'.format(imtype) in params.keys():
                params['calibs'] = params['{0}_calibs_read'.format(imtype)]
            im, im_head = read_file_calibration(params, params['master_{0}_filename'.format(imtype)], level=level)
            loaded = True
    if loaded == False:
        if '{0}_calibs_create'.format(imtype) not in params.keys():
            if 'standard_calibs_create' not in params.keys():
                logger('Error: Missing entry in the configuration file. Neigther "{0}_calibs_create" nor "standard_calibs_create" is given. Please update the configuration file(s).'.format(imtype))
            params['{0}_calibs_create'.format(imtype)] = params['standard_calibs_create']
        for i in range(len(params['{0}_calibs_create'.format(imtype)])):                                                                    # make it safe from different user input
            params['{0}_calibs_create'.format(imtype)][i] = params['{0}_calibs_create'.format(imtype)][i].lower()
            params['{0}_calibs_create'.format(imtype)][i] = params['{0}_calibs_create'.format(imtype)][i].replace('normaliz', 'normalis')
            params['{0}_calibs_create'.format(imtype)][i] = params['{0}_calibs_create'.format(imtype)][i].replace('normalisation', 'normalise')
        im, med_fluxes, std_fluxes = [], [], []
        if '{0}_rawfiles'.format(imtype) not in params.keys():
            logger('Error: The list of raw files for image type {0} is not defined in the configuration. Please check the configuration files.'.format(imtype))
        if len(params['{0}_rawfiles'.format(imtype)]) == 0:
            logger('Error: The list of raw files for image type {0} is empty. Please check the configuration files.'.format(imtype))
        num_imgs = len(params['{0}_rawfiles'.format(imtype)])                # how many images are expected
        header_updates = np.zeros((num_imgs,2))
        for im_index, imf in enumerate(params['{0}_rawfiles'.format(imtype)]):                   # Only works for maximum 40 images on neils machine
            params['calibs'] = params['{0}_calibs_create'.format(imtype)]       # get's overwritten when other files are being read
            img, im_head = read_file_calibration(params, imf, level=level)
            obsdate, obsdate_float, exposure_time, obsdate_begin, exposure_fraction = get_obsdate(params, im_head)    # unix_timestamp of mid exposure time
            header_updates[im_index,:] = [exposure_time, obsdate_float]
            med_flux = np.median(img, axis=None)
            med_fluxes.append(med_flux)
            std_fluxes.append(np.std(img, axis=None, ddof=1))
            if 'normalise' in params['{0}_calibs_create'.format(imtype)]:
                img = img/(med_flux+0.0)
            if im == []:                                                # Initiate the array with correct precission to avoid swapping
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
        for i in range(len(med_fluxes)):
            im_head['HIERARCH EXO_PIPE NORM_{0}'.format(i)] = med_fluxes[i]
        for i in range(len(std_fluxes)):
            im_head['HIERARCH EXO_PIPE STDV_{0}'.format(i)] = std_fluxes[i]
        if 'combine_mean' in params['{0}_calibs_create'.format(imtype)]:
            im = combine_sum(im)/(len(im)+0.0)
            im_head['HIERARCH EXO_PIPE redu07'] = 'Average of {0} images'.format(len(params['{0}_rawfiles'.format(imtype)]))
            exposure_time = np.mean(header_updates[:,0])                     # Average of the exposure times
        elif 'combine_sum' in params['{0}_calibs_create'.format(imtype)]:
            im = combine_sum(im)
            im_head['HIERARCH EXO_PIPE redu07'] = 'Sum of {0} images'.format(len(params['{0}_rawfiles'.format(imtype)]))
            exposure_time = np.sum(header_updates[:,0])                      # Sum of the exposure times
        else:           # Median combine
            im = combine_median(im)
            im_head['HIERARCH EXO_PIPE redu07'] = 'Median of {0} images'.format(len(params['{0}_rawfiles'.format(imtype)]))
            exposure_time = np.median(header_updates[:,0])                   # Median of the exposure times
        if 'normalise' in params['{0}_calibs_create'.format(imtype)]:
            norm_factor = np.median(med_fluxes)
            im = im * norm_factor
            im_head['HIERARCH EXO_PIPE NORM_MED'] = norm_factor
        im_head['HIERARCH EXO_PIPE '+params['raw_data_dateobs_keyword']] = datetime.datetime.utcfromtimestamp(np.median(header_updates[:,1])).strftime('%Y-%m-%dT%H:%M:%S.%f')
        im_head['HIERARCH EXO_PIPE '+params['raw_data_exptim_keyword']] = exposure_time
        first, last = np.argmin(header_updates[:,1]), np.argmax(header_updates[:,1])
        first = header_updates[first,1]-header_updates[first,0]/2.
        last  = header_updates[last, 1]+header_updates[last, 0]/2.
        im_head['HIERARCH EXO_PIPE BEGIN FIRST'] = datetime.datetime.utcfromtimestamp(first).strftime('%Y-%m-%dT%H:%M:%S.%f')
        im_head['HIERARCH EXO_PIPE END LAST']    = datetime.datetime.utcfromtimestamp(last ).strftime('%Y-%m-%dT%H:%M:%S.%f')
        im_head['HIERARCH EXO_PIPE EXP_RANGE']   = (last - first, 'sec, from BEGIN to END')
        if 'master_{0}_filename'.format(imtype) in params.keys():
            if params['master_{0}_filename'.format(imtype)] != '':
                save_im_fits(params, im, im_head,  params['master_{0}_filename'.format(imtype)])
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
        if np.sum(arr - arr.astype(np.float16), axis=None) == 0:
            arr.dtype = np.float16
        elif np.sum(arr - arr.astype(np.float32), axis=None) == 0:
            arr.dtype = np.float32
    return arr

def save_im_fits(params, im, im_head, filename):
    """
    Saves an image to a fits file                       # This procedure can possibly be combined with save_multispec
    :param params: Dictionary with all the parameters
    :param im: 2d array to be written into fits file
    :param im_head: header of the 2d array
    :param filename: filename in which to save the fits file
    """
    if len(filename.rsplit('/',1)) == 1:     # no path is in the filename
        logger('Warn: no folder to save {0} was given, using the current folder ({1}).'.format( filename, os.getcwd() ))
    elif not os.path.exists(filename.rsplit('/',1)[0]):
        logger('Error: Folder to save {0} does not exists.'.format(filename))
        
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
    
    if 'cen_poly0' not in atable.colnames:                                # Old way of storing the data, can be removed at a few months after July 2018
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
        return np.array(pfits), np.array(xlows), np.array(xhighs), np.array(widths)
    
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
    if not os.path.exists(filename.rsplit('/',1)[0]):
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

"""def save_fits(pfits, xlows, xhighs, filename):                      # not needed anymore
    ""
    save the fits to file, for use in extracting future spectra
    :param pfits: list, length same as number of orders, polynomial values
                      (output of np.polyval)
                      i.e. p where:
                      p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
    :param xlows: list, length same as number of orders, the lowest x pixel
                   (wavelength direction) used in each order
    :param xhighs: list, length same as number of orders, the highest x pixel
                   (wavelength direction) used in each order
    :param filename: string, location and file name of the file containing the polynomial fits
    ""
    data = []
    porder = len(pfits[0])                                                  #added by ronny
    cols = ['Order'] + list(range(porder + 1))[::-1]                        #added by ronny
    cols += ['low_x', 'high_x']

    # Loop around each order and plot fits
    for pp, pf in enumerate(pfits):
        row = [pp] + list(np.zeros(porder + 1 - len(pf)))               #added by ronny
        row += list(pf)
        row += [xlows[pp], xhighs[pp]]
        data.append(row)
    data = np.array(data)

    # convert to astropy table
    atable = Table()
    for c, col in enumerate(cols):
        atable[str(col)] = data[:, c]
    atable.write(filename, overwrite=True)                                     #added by ronny"""

def save_wavelength_solution_to_fits(wavelength_solution, wavelength_solution_arclines, filename):
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
    if not os.path.exists(filename.rsplit('/',1)[0]):
        logger('Error: Folder to save {0} does not exists.'.format(filename))
    
    data = []
    
    # Loop around each order and plot fits
    for pp, pf in enumerate(wavelength_solution):
        row = list(pf)
        row += list(wavelength_solution_arclines[pp])
        data.append(row)
    data = np.array(data)
    
    cols = ['real_order', 'central_px']
    cols += list(range(len(wavelength_solution[0]) -2))[::-1]
    cols += list([ 'arclin' + str(i) for i in range(len(wavelength_solution_arclines[0])) ])
    
    # convert to astropy table
    atable = Table()
    for c, col in enumerate(cols):
        atable[str(col)] = data[:, c]
    atable.write(filename, overwrite=True)
    logger('Info: Master file for arc solution written: {0}'.format(filename))

def read_wavelength_solution_from_fits(filename):
    """
    Reads the wavelength solution from a file
    :param filename: string, location and file name of the file
    fits file should look like the following:
    real_order central_px               2 ...             0 arclin1 ... arclin345
       float64    float64         float64 ...       float64  loat64 ...   float64
    ---------- ---------- --------------- ... ------------- ------- ... ---------
        -101.0     2150.0 -7.82221887e-07 ... 5.6899272e+03 5678.91 ...  6543.987
        -100.0     2150.0 -7.39532842e-07 ... 5.7469279e+03 6542.10 ...   6677.32
    ...
    where n, ..., 3, 2, 1, 0 are the polynomial powers in p
        i.e. p where:
            p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
    :return wavelength_solution: 2d array of floats, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    :return wavelength_solution_arclines: 2d array of floats, same length as number of orders, each line contains the wavelength of the used reference lines. 
                                          Each order must have the same number of entries, if less reference lines are available in the order, the array is filled with 0.
    """
    if os.path.isfile(filename) == False:
        logger('Error: The file {0} with the wavelength solution does not exist'.format(filename)) 
    # convert to astropy table
    atable = Table.read(filename)
    orders = np.array(atable['real_order'], dtype=int)

    wavelength_solution, wavelength_solution_arclines = [], []
    for order in range(len(orders)):
        entries, arclines = [], []
        for p in atable.colnames:
            if p.find('arclin') == 0:
                #if atable[p][order] != 0:      # zeros shouldn't be a problem for later usage
                    arclines.append(atable[p][order])
            else:
                entries.append(atable[p][order])        
        wavelength_solution.append(entries)
        wavelength_solution_arclines.append(arclines)
    logger('Info: Master file for arc solution read: {0}'.format(filename))
    return np.array(wavelength_solution), np.array(wavelength_solution_arclines)

def save_multispec(data, fname, im_head, bitpix='-32'):
    """
    Saves a n-d array into a fits file
    :param data: n-d array of data or list of arrays
    :param fname: filename in which the data should be written. The filename will end in .fits
    :param im_head: header, which should be written to the fits file
    :param bitpix: Precission in which the data should be stored
    """
    if len(fname.rsplit('/',1)) == 1:
        logger('Warn: No folder is given for the file {0}. File will be stored in the current working directory.'.format(fname))
    elif not os.path.exists(fname.rsplit('/',1)[0]):
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

def bin_im(im, binxy):
    """
    :param im: 2d numpy array
    :param binxy: list of integers (entries), contains the number of pixel which should be binned into one in x and y direction
    :return: 2d numpy arrays of the image, the number of elements which are not NaN, and of the standard deviation of each combined pixel
    """
    [binx, biny] = binxy
    ims = im.shape
    if binx <1 or biny <1 or (binx == 1 and biny == 1):
        logger('Warn: no binning possible: {0},{1}'.format(binx,biny))
        return(im)
    nim, gim, sim = [], [], []
    if binx > 1 and biny > 1:
        for i in range(int((ims[0]+binx-1)/binx)):
            nline, gline, sline = [], [], []
            iline = im[i*binx:min(ims[0],(i+1)*binx),:]
            for j in range(int((ims[1]+biny-1)/biny)):
                temdata = iline[:,j*biny:min(ims[1],(j+1)*biny)]
                if np.sum( np.isnan(temdata) ) == 0:
                    nline.append(np.median(temdata))    # median binning
                else:
                    nline.append(np.nan) 
                gtemp = np.sum( ~np.isnan(temdata) )
                gline.append(gtemp)    # number of the elements not nan
                if gtemp > 1:
                    sline.append(np.nanstd(temdata, ddof=1))    # standard deviation
                elif gtemp == 1:
                    sline.append(0)
                else:
                    sline.append(np.nan)   
                #nline.append(sum(percentile_list(list((im[i*binx:(i+1)*binx,j*biny:(j+1)*biny]).flat),.2)))   #sum after 80% clipping -> more flux -> more lines find with gauss of height 50
            nim.append(nline)
            gim.append(gline)
            sim.append(sline)
    elif binx > 1:
        for i in range(int((ims[0]+binx-1)/binx)):
            temdata = im[i*binx:min(ims[0],(i+1)*binx),:]
            nim.append(np.median(temdata,0))
            gim.append(np.sum( ~np.isnan(temdata), axis=0))  # number of the elements not nan
            sim.append(np.nanstd(temdata, ddof=min(1,temdata.shape[0]-1), axis=0))      # standard deviation, if only one column of pixels left, then std with ddof=0
    elif biny > 1:
        for i in range(int((ims[1]+biny-1)/biny)):
            temdata = im[:,i*biny:min(ims[1],(i+1)*biny)]
            nim.append(np.median(temdata,1))
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
        if np.min(im) <= 1:
            add_value = 1-np.min(im)
            im_plot = np.log10(im+add_value) - np.log10(add_value)
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

def twoD_Gaussian(x, y, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    """
    Calculates the Gauss in 2 dimensions
    :param (x,y): lists of the x- and y- values for which the Gauss should be calculated
    :param amplitude: Amplitude of the Gauss function
    :param x0: center of the Gauss in x direction
    :param y0: center of the Gauss in y direction
    :param sigma_x: Gaussian width on x direction
    :param sigma_y: Gaussian width on y direction
    :param theta: rotation of the Gauss compared to the x-axis
    :param offset: zero level of the Gauss
    """
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2.*sigma_x**2) + (np.sin(theta)**2)/(2.*sigma_y**2)
    b = -(np.sin(2*theta))/(4.*sigma_x**2) + (np.sin(2*theta))/(4.*sigma_y**2)
    c = (np.sin(theta)**2)/(2.*sigma_x**2) + (np.cos(theta)**2)/(2.*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

def unify_pord(pfits):                          # Only required because of the use of polynomialfit
    """
    Different orders are tested for each polynom. In order to create a np.array from of them zeros need to be added for the missing parameters
    :param pfits: 2d list of floats, the parameters for the polynoms
    :return: 2d array of floats, the parameters for the polynoms
    """
    # What is the maximum numbers of orders used for the polynomial fit
    porder = 0
    for pfit in pfits:
        porder = max(porder,len(pfit)-1)
    polyfits = []
    for pfit in pfits:
        polyfit = pfit
        while len(polyfit)-1 < porder:
            polyfit.insert(0,0)
        polyfits.append(polyfit)
    return np.array(polyfits)

"""def polynomialfit_notusedanymore(x, y, nmin=1, nmax=3):        # Not really necessary, the result is always that nmax fits best
    ""
    Select the minimum chi squared polynomial fit between nmin and nmax
    :param x: x axis array
    :param y: y axis array
    :param nmin: minimum order to fit
    :param nmax: maximum order to fit
    return ps: array of the parameters from the fit
    return nrange: number of parameters for the best fit, e.g. len(ps)
    ""
    chis, ps = [], []
    nrange = range(nmin, nmax+1)
    for n in nrange:
        p = np.polyfit(x, y, n)
        ps.append(p)
        chis.append(np.sum((y-np.polyval(p, x))**2))
    argmin = np.argmin(chis)
    return ps[argmin], nrange[argmin]"""

def polyfit_adjust_order(xarr, yarr, p_orders, w=None):
    """
    Finds the polynomial with the highest posible order to fit the data
    """
    xarr = np.array(xarr)
    yarr = np.array(yarr)
    poly = np.array([np.mean(yarr)])            # if the order = 0 is failing, then report at least the average (which is a polynom of order 0)
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
            Azz = np.hstack((A, zz.reshape((zz.shape[0],1)) ))
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

def fit_2d_image(im, xord, yord, w=[]):
    """
    Returns the fitted image
    :params im: 2d numpy array, i.e. the image
    :params xord, yord: number of orders for the polynomial to be used in each direction, yord is along dispersion axis
    :params w: weights of the individual points.
    """
    ims = im.shape
    cen = [ ims[0]/2. , ims[1]/2. ]
    x, y = np.meshgrid(np.arange(ims[1])-cen[1], np.arange(ims[0])-cen[0] )       # reverse order, otherwise the flatten doesn't create the right values for z
    xx, yy, zz = x.flatten(), y.flatten(), im.flatten()
    if len(w) > 0:
        w = w.flatten()

    poly2d_params = polynomial_fit_2d_norm(xx, yy, zz, xord, yord, w=w)       #x,y,z,order_x,order_y, w=weights
    im_fit = polynomial_value_2d(xx, yy, xord, yord, poly2d_params)
    #print np.array([xx,yy,zz, im_fit]).T
    im_fit.shape = ims
    #plot_img_spec.plot_image(im, [''], 1, True, [0.05,0.95,0.95,0.05], 'orig')
    #plot_img_spec.plot_image(im_fit, [''], 1, True, [0.05,0.95,0.95,0.05], 'fit')
    #plot_img_spec.plot_image(im-im_fit, [''], 1, True, [0.05,0.95,0.95,0.05], 'residuals')
    
    return im_fit

def sigma_clip(xarr, yarr, p_orders, sigma_l, sigma_h, repeats = 1):
    """
    Performs a sigmaclipping after fitting a polynomial against data
    :param xarr: 1d array with the x-values of the data
    :param yarr: 1d array with the data
    :p_orders: orders of the polynomial
    :sigma_l: Data off by this sigma are rejected on the lower side of the fit
    :sigma_h: Data off by this sigma are rejected on the higher side of the fit
    :repeats: redo the fit with the (cleaned) data how many times? 
    :return goodvalues: 1d array with True/False. The values which are inside the limits are True
    :return p: parameters of the last polynomial fit
    """
    xarr, yarr = np.array(xarr), np.array(yarr)
    if xarr.shape[0] == 0 or yarr.shape[0] == 0:
        logger('Warn: empty array for sigma clipping')
        return np.array([]), np.repeat([0], p_orders)
    elif xarr.shape[0] != yarr.shape[0]:
        logger('Warn: got different sized arrays for sigma clipping. This is a programming error and should not happen')
        return [], np.repeat([0], p_orders)
    goodvalues = (yarr*0 == 0)
    old_values = [0,0]
    for i in range(repeats):
        poly = polyfit_adjust_order(xarr[goodvalues], yarr[goodvalues], p_orders)
        stddiff = np.std(yarr[goodvalues] - np.polyval(poly, xarr[goodvalues]), ddof=p_orders+1)
        diff = (yarr - np.polyval(poly, xarr))/(stddiff+0.0)
        goodvalues = ( (diff >= -sigma_l) & (diff <= sigma_h) )    #average should be 0
        if stddiff == old_values[0] and poly[0] == old_values[1]:
            break
        old_values = [stddiff,poly[0]]
    #plot_img_spec.plot_spectra(np.array([xarr,xarr]),np.array([yarr,np.polyval(poly, xarr)]),['data','fit'], ['savepaths'], True, [0.06,0.92,0.95,0.06, 1.0,1.01], 'line {0}'.format(i))
    return goodvalues, poly

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

def centroid_order(x, y, center, width, significance=3, blended_gauss=False, bugfix=False):
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
    bounds=((0.2*rangey/5., center-max(1.2,0.3*width), 0.1*width, min(y)-0.2*rangey), (5.*rangey, center+max(1.2,0.3*width), 2.*width, min(y)+0.2*rangey))
    if blended_gauss:
        p0 = p0 + p0
        bounds = (bounds[0] + bounds[0], bounds[1] + bounds[1])
    stdlin = np.std(y, ddof=1)
    significant = False
    for border_sub in [ [0,0], [0,1], [1,0], [1,1], [2,1], [1,2], [2,2] ]:    
        range_data = list(range(border_sub[0],len(x)-border_sub[1]))
        if len(range_data) <= 6:
            break
        if not blended_gauss:
            try:
                popt,pcov = curve_fit(oneD_gauss,x[range_data],y[range_data],p0=p0, bounds=bounds)            #a, x0, sigma, b: a*np.exp(-(x-x0)**2/(2*sigma**2))+b
            except:
                # print 'curve fit failed'
                continue
            stdfit = np.std(oneD_gauss(x[range_data],popt)-y[range_data], ddof=len(popt))           # if the average is 0, than this gives the stray of the residuals
            if stdlin/(stdfit+0.0) >= significance or popt[0]/(stdfit+0.0) >= significance:
                significant = True
                break
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
                break
            elif bugfix:
                print('Gauss not significant', p0,bounds, stdfit, border_sub)
                plot_img_spec.plot_spectra(np.array([x,x]),np.array([y,oneD_blended_gauss( x, popts)]),['data','fit'], ['savepaths'], True, [0.06,0.92,0.95,0.06, 1.0,1.01], 
                                           'No significant fit, stdlin={0}, stdGauss={1}, height={2}, needed significance={3}'.format(stdlin,stdfit,popt[0],significance))
            
    #print 'stdfit,stdlin', stdfit,stdlin,stdlin/stdfit < significance, popt[0]/stdfit < significance, popt
    if not significant:   # fit is not significant
        return  np.array([0,0,0,0])
    return popt

def find_border_pctl(data, border_pctl=50):
    """
    Calculates the position where the signal falls below a certain threashold. If enough data points are available, then a linear fit is done in the transition area
    :param data: 1d list or array of floats, contains a Gauss-like signal
    :param border_pctl: integer or float, defines the percentage of the signal at which the cut should be done
    :return result[0]: float, position on the left of the maximum where the signal in data falls below border_pctl percent
    :return result[1]: float, position on the right of the maximum where the signal in data falls below
 border_pctl percent
    """
    min_dat, max_dat = min(data), max(data)
    ran_dat = max_dat - min_dat
    higher = np.where( data > ran_dat*border_pctl/100. + min_dat)[0]         # Find the indexes above the percentile value
    higher = [ max(1, min(higher)) , min(len(data)-2, max(higher)) ]    # To allow a 2px fit, exclude the outermost pixels
    lower = np.where( data < ran_dat*border_pctl/100. + min_dat)[0]
    lower1 = lower[ lower < higher[0] ]                                 # The left side lower than the percentile value
    if len(lower1) < 1:                                                 # happens is the line is cut before it reaches the lower level
        lower1 = [0]
    lower2 = lower[ lower > higher[1] ]                                 # The righ side lower than the percentile value
    if len(lower2) < 1:                                                 # happens is the line is cut before it reaches the lower level
        lower2 = [len(data)-1]
    ran1 = list(range( max(lower1), higher[0]+1 ))                      # This should be a 2px (or more) transition around the transition
    ran2 = list(range( higher[1], min(lower2)+1 ))                      # This should be a 2px (or more) transition around the transition
    result = []
    for ran in [ran1, ran2]:
        if np.std(data[ran], ddof=1) == 0.0:                            # same flux values in all pixel -> fit is not possible
            result.append(np.median(ran))
        else:
            pfit = np.polyfit( data[ran] - min_dat, ran, 1)             # linear fit px over data (corrected by minimum)
            result.append( np.polyval(pfit, ran_dat*border_pctl/100.) )
        
    return result[0], result[1]
    
def find_center(imslice, oldcenter, x, maxFWHM, border_pctl=0, significance=3.0, bugfix=False):
    """
    Determines the center of a gauss
    :param imslice: 1d array to search for the center
    :param oldcenter: integer, initial guess of the center (e.g. maximum)
    :param x: integer, index of the 1d array in the overal image, only used for error investigation
    :param maxFWHM: integer, estimated width of the order
    :param significance: float, Significance of the Gaussian fit compared to a simple average
    :return centerfit: float, center from the centroid fit
    :return width: float, half width of the order at the base, set to 0 in case the fit wasn't good
    :return leftmin: integer, left border of the line
    :return rightmin: integer, right border of the line
    """
    maxFWHM = max(3, maxFWHM)           # To have enough pixel (in order to fit the Gauss)
    oldcenter = max(0, oldcenter)       # The oldcenter needs to be in the data
    oldcenter = min(len(imslice)-1,oldcenter)
    width, center, leftmin, rightmin = width_order(imslice,oldcenter,maxFWHM)
    if width > maxFWHM or width == 0:                                        
        width = 0
        centerfit = center
    else:
        while rightmin-leftmin < max(5,2*maxFWHM):             # To have enough pixel to fit the Gauss, 2 times the maximum FWHM
            leftmin = max(0, leftmin-1)
            rightmin = min(len(imslice)-1, rightmin+1)
        sub_slice = imslice[leftmin:rightmin+1]
        if len(range(leftmin,rightmin+1)) != len(sub_slice):
            print('Error: not same length', x, oldcenter, len(imslice), leftmin, rightmin)
        for border in [0,1]:
            popt = centroid_order(list(range(leftmin+border,rightmin+1-border)), sub_slice[border:len(sub_slice)-border], center, width, significance=significance, bugfix=bugfix)
            if popt[2] > 0:                 # run the second fit only when the first fails
                break 
        centerfit = popt[1]
        width = popt[2]
        if border_pctl > 0:
            left, right = find_border_pctl(sub_slice, border_pctl=border_pctl)          # left, right in terms of sub_slices pixels
            rightmin = right + leftmin
            leftmin = left + leftmin

    return centerfit, width, leftmin,rightmin

def measure_noise(spectra, p_order=12, semi_range=10):
    """
    Measures the noise in a spectrum by fitting a high order polynomial against the data and deriving the residuals
    :param spectra: 2d array of floats, spectrum
    :param p_order: integer, order of the polynomial
    :param semi_range: integer, the noise is calculated using a sample of +/- semi_range pixels around the current value
    :return noise: 2d array of floats, same dimension as spectra
    """
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
        
    return np.array(noise)

def estimate_width(im):
    """
    Finding a guess of the Gaussian width of the traces in an image by using a set of random positions
    :param im: 2d array of float, image with the orders
    :return width: float, median of the measured Gaussian widths
    """
    ims = im.shape
    widths = []
    logger('Info: Estimate the width of the traces')
    maxtests = max( 200, int(ims[0]*ims[1]/50) )
    for i in range(maxtests):     # Use at least 200 positions in order to be sure to find at least 10
        x = random.randint(int(0.05*ims[0]), int(0.95*ims[0]))
        y = random.randint(int(min(50,0.05*ims[1])), int(max(ims[1]-50,0.95*ims[1])))
        yr = list(range(max(0,y-30), min(ims[1],y+31)))
        data = im[x,yr]
        pos_max = data.argmax() + max(0,y-30)       # Maximum in total y (from im)
        widths1 = []
        for w1 in range(2,10,1):
            for w2 in range(1,4,1):
                yy = im[x, max(0,pos_max-w1*w2):min(ims[1],pos_max+w1*w2+1)]
                if len(yy) < 5:
                    continue
                xx = list(range(len(yy)))
                popt = centroid_order(xx, yy, int(len(xx)/2), w1)    # x,y,pos,width ; result: a,x0,sigma,b in a*np.exp(-(x-x0)**2/(2*sigma**2))+b
                
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

def find_trace_orders(params, im):
    """
    :param params: Dictionary with all the parameters. 'binx', 'biny', 'subframe', and 'polynom_order_apertures' are required
    :param im: 2d numpy array
    return: array of the parameters for the polynomial, orders ordered from left
    return: minimum and maximum position of the trace in dispersion direction
    """
    binx = params['binx']
    biny = params['biny']
    maxFWHM = max(1, int(round(estimate_width(im)*2.35482*1.5)))   # The average Gaussian width, transformed into a FWHM and extendet, as this is the maximum FWHM
    maxshift = 0.15*binx                    # Only 0.15px shift per pixel along one order
    ims = im.shape
    cen_px = int(ims[0]*binx/2)             # central pixel to base the fit on
    breakvalue = np.percentile(im,25)       # The brightes pixel of one order needs to be brighter than this value, otherwise it won't be identified (lower values -> darker orders are found)
    im_orig = copy.deepcopy(im)
    im_traces = np.zeros(ims)               # Masks the orders in order to avoid overlapping traces
    traces = []
    trace_pos = [list(range(ims[0]))]     # fill array with index of slice, later for each order found the y values calculated from the polyfit are added -- this entry is not needed
    for dummy in tqdm(range(max( 2600, int(ims[0]*ims[1] * maxFWHM*2 / 1000) )), desc='Searching for traces'):      # 2600 is necessary for the tight orders using Clark's lens and a low diffraction prism
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
        for i in range(pos_max[0],-1,-1):               # check positions upwards
            center, width, leftmin,rightmin = find_center(im_orig[i,:], int(round(oldcenter)), i, maxFWHM, significance=3.5)       # significance=4.0 tested as useful for HARPS, EXOhSPEC
            if pos_max[1] >= 2000 and pos_max[1] <= 260:            # bugfixing
                if not ( width != 0 and ( abs(center-oldcenter)<maxshift or (positions.shape[0]==0 and abs(center-oldcenter)<maxshift*2) ) ):
                    find_center(im_orig[i,:], int(round(oldcenter)), i, maxFWHM, significance=3.5, bugfix=True)
                #print(pos_max, i, center, width, leftmin,rightmin, oldcenter, abs(center-oldcenter), maxshift, len(positions))
            if width != 0 and ( abs(center-oldcenter)<maxshift or (positions.shape[0]==0 and abs(center-oldcenter)<maxshift*2) ):        #first entry can be further off
                if im_traces[i,min(ims[1]-1,int(center))] != 0 or im_traces[i,min(ims[1]-1,int(center+1))] != 0:        # this order shouldn't cross another order
                    order_overlap = True
                    break
                positions = np.vstack([positions, [i, center, 0]])
                widths.append([center-leftmin,rightmin-center, width])
                oldcenter, lastadd, last_trustworth_position, no_center  = center, i, len(positions)-1, 0
            else:#if expected_positions != []:
                no_center += 1
                if abs(oldcenter - 0.) < maxshift/2. or abs(oldcenter - ims[1]) < maxshift/2.:  # if the trace leaves shortly the CCD
                    no_center -= 0.6
                if no_center >= max(3, int(40/binx), 1*maxFWHM):         # stop searching if too many fits are unseccessful, as otherwise the fit might drift off
                    break
                center1, width, leftmin,rightmin = find_center(im_orig[i,:], int(round(oldcenter)), i, maxFWHM, significance=3.0, bugfix=False)
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
            center, width, leftmin,rightmin = find_center(im_orig[i,:], int(round(oldcenter)), i, maxFWHM, significance=3.5)       # significance=4.0 tested as useful for HARPS, EXOhSPEC
            if pos_max[1] >= 2000 and pos_max[1] <= 260:            # bugfixing
                if not ( width != 0 and ( abs(center-oldcenter)<maxshift or (positions.shape[0]==0and abs(center-oldcenter)<maxshift*2) ) ):
                    find_center(im_orig[i,:], int(round(oldcenter)), i, maxFWHM, significance=3.5, bugfix=True)
                print(pos_max, i, center, width, leftmin,rightmin, oldcenter, abs(center-oldcenter), maxshift, len(positions), len(expected_positions), last_trustworth_position)
            if width != 0 and ( abs(center-oldcenter)<maxshift or (positions.shape[0]==0and abs(center-oldcenter)<maxshift*2) ):
                if im_traces[i,min(ims[1]-1,int(center))] != 0 or im_traces[i,min(ims[1]-1,int(center+1))] != 0 or order_overlap == True:
                    order_overlap = True
                    break
                positions = np.vstack([positions, [i, center, 0]])
                widths.append([center-leftmin,rightmin-center, width])
                oldcenter, lastadd, last_trustworth_position, no_center  = center, i, len(positions)-1, 0
            else:#if expected_positions != []:
                no_center += 1
                if abs(oldcenter - 0) < maxshift/2. or abs(oldcenter - ims[1]) < maxshift/2.:  # if the trace leaves shortly the CCD
                    no_center -= 0.6
                if no_center >= max(3, int(40/binx), 1*maxFWHM):         # stop searching if too many fits are unseccessful, as otherwise the fit might drift off
                    break
                center1, width, leftmin,rightmin = find_center(im_orig[i,:], int(round(oldcenter)), i, maxFWHM, significance=3.0, bugfix=False)
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
        width = [np.mean(percentile_list(np.array(widths)[:,0],0.1)), np.mean(percentile_list(np.array(widths)[:,1],0.1)), np.mean(percentile_list(np.array(widths)[:,2],0.1))]      #average left, right, and gaussian width
        #print positions
        # Fit to the original image
        #pfs, ns = polynomialfit(positions[:,0]*binx, positions[:,1]*biny, 1, params['polynom_order_apertures'])     # (positions[:,1]-ims[0]/2)*biny for fit in the center of the image, but this needs to be taken into account for the next step
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
        traces.append([[cen_px]+list(pfs), min(positions[:,0])*binx, min(params['subframe'][0],(max(positions[:,0])+1)*binx-1), np.polyval(pfs, cen_px)])
        logger('Step: order found at central pixel: {0} / {3} (frame/trace). The trace was identified between Pixels {1} and {2}'.format(round(np.polyval(pfs, cen_px),1),
                                                round(traces[-1][1]), round(traces[-1][2]), round(np.polyval(pfs, np.mean(traces[-1][1:3]) ),1) ))
    # sort
    traces = np.array(traces)
    cen_px2 = np.mean([ np.max(traces[:,1]), np.min(traces[:,2]) ])      # Center of the area covered by all orders
    for order in range(traces.shape[0]):
        pfs = traces[order, 0]
        if cen_px2 >= traces[order, 1] and cen_px2 <= traces[order, 2]:
            traces[order, 3] = np.polyval(pfs[1:], cen_px2 - pfs[0])     # Position in Cross-dispersion direction
        else:
            logger('Warn: Please check that all orders have been identified correctly')
            if abs(cen_px2 - traces[order, 1]) < abs(cen_px2 - traces[order, 2]):     # Use the closest value for which the trace is still defined
                traces[order, 3] = np.polyval(pfs[1:], traces[order, 1] - pfs[0])
            else:
                traces[order, 3] = np.polyval(pfs[1:], traces[order, 2] - pfs[0])
    traces = np.array(sorted(traces, key=operator.itemgetter(3)))
    #plot_img_spec.plot_image(im+im_traces*np.max(im), ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'cleared')
    if traces.shape[0] < 2:
        logger('Error: Only found {0} traces. Please check that the binned image ({1}) looks as expected.'.format(traces.shape[0], params['logging_trace1_binned'])+\
               'Please check that the right image was used to search for the traces (e.g. {0}) '.format(params['logging_traces_im_binned'])+\
               'After you selected different file(s) for parameter "trace1_rawfiles", please run the following command before restarting the pipeline'+\
               '\nrm {0}'.format( params['master_trace1_filename'] ))
    logger('Info: {0} orders found and traced'.format(traces.shape[0]))
    return traces[:,0], traces[:,1].astype(int), traces[:,2].astype(int)            # polyfits, xlows, xhighs

def adjust_trace_orders(params, im, im_unbinned, pfits, xlows, xhighs):
    """
    Re-traces the position of the orders
    :param params: Dictionary with all the parameters. 'binx', 'biny', 'subframe', and 'polynom_order_apertures' are required
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
    cen_px = int(ims[0]*binx/2)
    # maxFWHM = 25/biny for old lens
    # maxFWHM = 15/biny+2            # !! Determine in values of space between orders
    maxFWHM = max(1, int(round(estimate_width(im)*2.35482*1.5)))   # The average Gaussian width, transformed into a FWHM and extendet, as this is the maximum FWHM
    if len(pfits) > 1:                              # move at maximum by 20% of the minimum space between the orders, but at old binning in cross-dispersion direction is allowed
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
    logger('Step: Adjusting orders', show=False)
    for order in tqdm(range(trace_pos.shape[0]-1), desc='Adjusting orders'):
        lastadd = trace_pos[0,0]
        positions, widths_o, shifts = [], [], []
        for j,i in enumerate(trace_pos[0,:].astype(int)):           # each pixel, is is the position in the array and i the real position on the CCD
            oldcenter = trace_pos[order+1,j]                        # +1 because first entry are pixels
            if oldcenter > -maxFWHM:                                # center is expected to be in the image
                center, width, leftmin,rightmin = find_center(im[i,:], int(round(oldcenter)), i, maxFWHM, border_pctl=params['width_percentile'], significance=3.0)
            else:
                center, width, lastadd = 0,0, i
            #if width == 0 and oldcenter > -maxFWHM:
            #if width != 0 and abs(center-oldcenter)>=3:
            #    print center,oldcenter,i, leftmin,rightmin
            if width != 0 and abs(center-oldcenter)<maxshift:
                #print 'maxFWHM,order,j,i,leftmin,center, rightmin, width',maxFWHM,order,j,i,leftmin,center, rightmin, width
                positions.append([i, center,0])
                if leftmin < center-width*2.35482*3 or leftmin > center:        # unreasonable small position
                    leftmin = np.nan
                if rightmin > center+width*2.35482*3 or rightmin < center:      # unreasonable large position
                    rightmin = np.nan
                #leftmin = max(leftmin, center-width*2.35482*3)          # anything bigger than 
                #rightmin = min(rightmin, center+width*2.35482*3)
                widths_o.append([i, leftmin, rightmin, width])         # leftmin and rightmin are positions, they are transformed into a width below
                lastadd = i
                shifts.append(center-oldcenter)
            elif oldcenter > -maxFWHM and i-lastadd == 10:          # If inside the image add the center after every 10 data points
                positions.append([i, oldcenter,1])
                lastadd = i
            #else:
            #    print center, oldcenter, maxshift, maxFWHM
                
        #print len(positions), 
        if widths_o == [] or len(shifts) < ims[0]/4.:                # Not enough data for that order
            #print 'len(positions),len(shifts)', len(positions),len(shifts)
            continue
        positions = np.array(positions)
        widths_o = np.array(widths_o)
        # Fit to the original image
        """
        pfs, ns = polynomialfit(positions[:,0]*binx, positions[:,1]*biny, 1, params['polynom_order_apertures'])
        # Remove values which were added only to follow the curve
        positions = positions[positions[:,2]==0,0:2]
        traces.append([ list(pfs) ,ns, max(0,min(positions[:,0])*binx-5), min(params['subframe'][0],(max(positions[:,0])+1)*binx+5), list(width) ])
        #print len(positions), min(positions[:,0]),(max(positions[:,0])+1)*binx-1
        logger('Step: order {0} adjusted: avg shift = {1}, min shift = {2}, max shift = {3}, avg width: gauss = {6}, left = {4}, right = {5}'.format(order, round(np.mean(shifts),2),round(min(shifts),2),round(max(shifts),2), round(width[0],2),round(width[1],2),round(width[2],2) ))"""
        # Fit the center
        pfs = np.polyfit(positions[:,0]*binx - cen_px, positions[:,1]*biny, params['polynom_order_apertures'])
        #centerfit.append(pfs)
        # Fit the borders
        good_vaules = ~np.isnan(widths_o[:,1])
        if sum(good_vaules) > len(good_vaules)/20.:
            leftfit =  np.polyfit(widths_o[good_vaules,0]*binx - cen_px, widths_o[good_vaules,1]*biny, params['polynom_order_apertures'])
            widthleft = np.nanmean(percentile_list(widths_o[:,1],0.1))
        else:                                                   # left side at FWHM of line
            widthleft = np.nanmean(percentile_list(widths_o[:,3],0.1)) * 2.35482
            leftfit = copy.deepcopy(pfs)
            leftfit[-1] -= widthleft
        good_vaules = ~np.isnan(widths_o[:,2])
        if sum(good_vaules) > len(good_vaules)/20.:
            rightfit = np.polyfit(widths_o[good_vaules,0]*binx - cen_px, widths_o[good_vaules,2]*biny, params['polynom_order_apertures'])
            widthright = np.nanmean(percentile_list(widths_o[:,2],0.1))
        else:                                                   # left side at FWHM of line
            widthright = np.nanmean(percentile_list(widths_o[:,3],0.1)) * 2.35482
            rightfit = copy.deepcopy(pfs)
            rightfit[-1] += widthright
        centerfit.append([[cen_px] + list(pfs), [cen_px] + list(leftfit), [cen_px] + list(rightfit)])
        center = np.polyval(pfs, widths_o[:,0]*binx - cen_px)
        widths_o[:,1] = center - widths_o[:,1]*biny
        widths_o[:,2] = widths_o[:,2]*biny - center
        width = [widthleft, widthright, np.nanmean(percentile_list(widths_o[:,3],0.1))]      #average left, right, and gaussian width
        xlows.append(max(0,min(positions[:,0])*binx-5))
        xhighs.append(min(im_unbinned.shape[0],(max(positions[:,0])+1)*binx+5))
        widths.append(width)
        avg_shifts.append(np.mean(shifts))
    """traces = np.array(traces)"""
    if len(widths) == 0:
        logger('Error: When adjusting the traces, all traces were rejected due to poor fit. '+\
                'Please check the binned image {0}: Is the orientation and the binning right? Are the orders covering at least half of the CCD (in dispersion correction)'.format(params['logging_trace1_binned'])+\
               'Please check that the right image was used to search for the traces (e.g. {0}) '.format(params['logging_traces_im_binned'])+\
               'After you selected different file(s) for parameter "trace1_rawfiles", please run the following command before restarting the pipeline'+\
               '\nrm {0} {1}'.format( params['logging_traces_binned'], params['master_trace1_filename'] ))
    logger('Info: traces of the {0} apertures adjusted. The average shift of the individual apertures was between {1} and {2} pixel between the searching of the traces and this solution. The maximum allowed shift was {3} pixel.'.format(len(centerfit), np.round(np.min(avg_shifts),1), np.round(np.max(avg_shifts),1), np.round(maxshift,1) ))

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

def asign_bad_px_value(params, input_value, badpx_mask, im, data_range):
    """
    :param input_value: already asigned value for the bad pixel mask
    :param badpx_mask: n-d array (normally 2d array) of int or float, with 0 to mark the bad pixels
    :param im: n-d array (normally 2d array) of int or float, image from which the extraction will be done
    :param data_range: number, list, or array of integers, to mark the extracted area
    """
    if np.nanmin(badpx_mask[data_range]) == 0:
        input_value = 0.2
    elif np.nanmax(im[data_range]) >= params['max_good_value']:
        input_value = 0.1
    
    return input_value

def no_tqdm(input, desc=''):
    """
    In order to switch between tqdm on and off, this is needed
    """
    return input

def extract_orders(params, image, pfits, xlows, xhighs, widths, w_mult, offset = 0, var='fast', plot_tqdm=True):
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
    :param var: can be 'fast' or 'prec'.
        'fast': extraction of the fraction of the border pixel
        'prec': fit a polynomial around the border of the extrion width to extract with higher precission
    :return spec: list,  length same as number of orders, each containing
                  a 1D numpy array with the spectrum of each order in

    """
    if plot_tqdm:
        plot_tqdm = tqdm
    else:
        plot_tqdm = no_tqdm
    logger('Step: extracting spectra', show=False)
    badpx_mask_name = 'badpx_mask'
    if badpx_mask_name not in calimages:
        calimages[badpx_mask_name] = read_badpx_mask(params)
    badpx_mask = calimages[badpx_mask_name]
    spec, good_px_mask = [], []
    xarr = range(0, image.shape[0])
    maxpx = image.shape[1]-1                    # End of the image in cross-dispersion direction
    for pp in plot_tqdm(range(pfits.shape[0]), desc='Extract Spectrum'):        # pp is order
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
        lowers, uppers = adjust_width_orders(yarr, lowers, uppers, [w_mult, w_mult])          # Adjust the width of the orders
        #print pp,np.nanmean(uppers-lowers), np.nanmedian(uppers-lowers), np.nansum(uppers-lowers)
        if w_mult == 0.0:                                                       # only the central full pixel
            order_in_frame = ((yarr[xarr1] <= maxpx) & (yarr[xarr1] >= 0) )     # Due to order shifting the order might be shifted to the ouside of the images, ignore these values
            for xrow in xarr1[order_in_frame]:
                ssum = image[xrow, int(round(yarr[xrow]))]                      # Verified 20190122
                good_px = 1
                good_px = asign_bad_px_value(params, good_px, badpx_mask, image, (xrow, int(round(yarr[xrow]))) )
                ospecs[xrow] = ssum
                ogood_px_mask[xrow] = good_px  
        elif var == 'prec':         # new solution, better than the one below, tested with many values in excel
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
                    good_px = asign_bad_px_value(params, good_px, badpx_mask, image, (xrow, range(lowersf[xrow],uppersf[xrow]+1)) )
                # Get fractions of pixel from a polynom, if they are not outside the frame borders (spline actes weired between points, when there is an outlier)
                #if w[2]*w_mult < 0.5:             # There will be no full pixel, which messes up everything
                if uppers[xrow] - lowers[xrow] < 1.0:       # There will be no full pixel, which messes up everything
                    fracpixparam = [[yarr[xrow],'b']]
                else:
                    fracpixparam = [[lowers[xrow], 'l'] , [uppers[xrow],'u']]
                for [borderpos, pos] in fracpixparam:
                    x = np.arange(max(0,np.floor(borderpos-1)), min(maxpx,np.ceil(borderpos+1))+1, dtype=int)
                    good_px = asign_bad_px_value(params, good_px, badpx_mask, image, (xrow, x) )
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
                        print('Programming error around line 2100')
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
                    good_px = asign_bad_px_value(params, good_px, badpx_mask, image, (xrow, range(lowersf[xrow],uppersf[xrow]+1)) )
                    # Tested on 15/3/19 for individual entries and ranges of entries: image[xrow, lowersf[xrow]:uppersf[xrow]+1] is the same as data_range = (xrow, range(lowersf[xrow],uppersf[xrow]+1)) ; image[data_range]
                    """if min(badpx_mask[xrow, lowersf[xrow]:uppersf[xrow]+1]) == 0:                       # Bad pixel in extracted data
                        good_px = 0.2
                    elif max(image[xrow, lowersf[xrow]:uppersf[xrow]+1]) >= params['max_good_value']:   # Saturated pixel in extracted data
                        good_px = 0.1"""
                #print 'xrow, ssum, ...', xrow, ssum, lowers[xrow], uppers[xrow], lowersf[xrow], uppersf[xrow]+1, 
                # Get fractions of pixel, if they are not outside the frame borders
                if lowersr[xrow] >= 0:
                    ssum += image[xrow, lowersr[xrow]]*( (1.5 - lowers[xrow]%1)%1 )     # lowers[xrow]%1=[0.5,0.0,0.6,0.4] -> [1.0,1.5,0.9,1.1] -> [.0,.5,.9,.1]
                    good_px = asign_bad_px_value( params, good_px, badpx_mask, image, (xrow, lowersr[xrow]) )
                #print 'l',ssum, lowersr[xrow],(1.5 - lowers[xrow]%1)%1,
                if uppersr[xrow] <= maxpx:      # maxpx is shape-1
                    ssum += image[xrow, uppersr[xrow]]*( (0.5 + uppers[xrow]%1)%1 )     # uppers[xrow]%1)=[0.5,0.0,0.6,0.4] -> [1.0,0.5,1.1,0.9] -> [.0,.5,.1,.9]
                    good_px = asign_bad_px_value( params, good_px, badpx_mask, image, (xrow, uppersr[xrow]) )
                #print 'u', ssum, uppersr[xrow],(0.5 + uppers[xrow]%1)%1
                if lowersr[xrow] == uppersr[xrow] and lowersr[xrow] >= 0 and uppersr[xrow] <= maxpx:
                    ssum += image[xrow, lowersr[xrow]]*( - 1 )     # remove the one added too much when using the fractions on the same pixel
                    #print 'f', ssum, uppersr[xrow],-1
                ospecs[xrow] = ssum
                ogood_px_mask[xrow] = good_px
                #if pp == 57 and xrow == 1350:
                #    print pp,xrow, image[xrow, lowersf[xrow]:uppersf[xrow]+1], image[xrow, lowersr[xrow]], (0.5-lowers[xrow]%1), image[xrow, uppersr[xrow]], (uppers[xrow]%1-.5), ssum

        """
        # new solution, not helpful solution:
        yarr = np.polyval(pf, xarr)+offset
        xarr1 = np.arange(xlows[pp], xhighs[pp])
        width = widths[pp][2]*w_mult
        # Define the borders of the extraction area
        lowers = yarr - width
        uppers = yarr + width
        lowers[lowers<0] = 0
        uppers[uppers>maxpx] = maxpx
        # Define the borders of the fitting area
        lowerw = np.floor(yarr - max(3, width*2))   #mindestens 3 px in each direction
        upperw = np.ceil(yarr + max(3, width*2))    #mindestens 3 px in each direction
        lowerw[lowerw<0] = 0
        upperw[upperw>maxpx] = maxpx
        order_in_frame = ((lowerw[xarr1] <= maxpx) & (upperw[xarr1] >= 0) & (upperw[xarr1] - lowerw[xarr1] > 2*width))     # Due to order shifting the order might be shifted to the ouside of the images, ignore these values
        ospecs = np.repeat([np.NaN], image.shape[0])
        for xrow in xarr1[order_in_frame]:
            x = np.arange(lowerw[xrow], upperw[xrow]+1, dype=int)
            y = image[xrow, x]
            # Better than polynom is spline, and better than fitting everything is fitting only the ends: 
            # https://stackoverflow.com/questions/17913330/fitting-data-using-univariatespline-in-scipy-python
            # https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.interpolate.UnivariateSpline.integral.html
            #p = np.polyfit(x, y, len(x)/2)
            #diff = y - np.polyval(p, x)
            #if xrow in range(500,501):
            #    print y
            #    print diff/y
            #    print np.mean(diff), np.std(diff, ddof=len(x)/2), np.sum(diff),
            #    diff = np.abs(diff)
            #    print np.mean(diff), np.std(diff, ddof=len(x)/2), np.sum(diff)
            #poly = np.poly1d(p)
            #polyint = poly.integ()
            #ospecs[xrow] = polyint(uppers[xrow]) - polyint(lowers[xrow])
        """    
        spec.append(ospecs)
        good_px_mask.append(ogood_px_mask)
    return np.array(spec), np.array(good_px_mask)

def shift_orders(im, params, sci_tr_poly, xlows, xhighs, oldwidths, in_shift = 0):
    """
    Determines the shift in spacial direction between the current flat {im} and an available solution. This is done by fitting the center again
    :param im: 2d numpy array with the flat lines
    :param params: Dictionary with all the parameters. 'maxshift' is required
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
    shifts, twidths, positions = [], [], []
    problem_order = []
    shift_map = np.zeros((im.shape[0],sci_tr_poly.shape[0]))
    logger('Step: Searching for shifts of the orders', show=False) 
    for i in tqdm(range(sci_tr_poly.shape[0]), desc='Searching for shifts of the orders'):       # For each order
        xarr = range(int(xlows[i]),int(xhighs[i]),steps)
        yarr = np.polyval(sci_tr_poly[i, 0, 1:], xarr-sci_tr_poly[i, 0, 0])+in_shift          #apply input shift
        widths = []
        for j,oldcenter in enumerate(yarr):        #each pixel
            center, width, leftmin,rightmin = find_center(im[xarr[j],:], int(round(oldcenter)), xarr[j], 2*params['maxshift'], border_pctl=params['width_percentile'], significance=3.0)
            #if width != 0:
            #    print 'i,j,center, width, leftmin,rightmin',i,j,center, width, leftmin,rightmin
            if width != 0 and abs(center-oldcenter) < params['maxshift']:
                widths.append([center-leftmin,rightmin-center, width])
                shifts.append(center-oldcenter)
                shift_map[max(0,int(xarr[j]-steps/2)):min(im.shape[0],int(xarr[j]+steps/2))+1,i] = center-oldcenter
        if widths == []:
            twidths.append(oldwidths[i])
            problem_order.append(i)
            continue
        widths = np.array(widths)
        width = [np.mean(percentile_list(widths[:,0],0.1)), np.mean(percentile_list(widths[:,1],0.1)), np.mean(percentile_list(widths[:,2],0.1))]      #average left, right, and gaussian width
        twidths.append(width)
    if shifts != []:
        shift = np.mean(percentile_list(np.array(shifts),0.1))
        shift_error = np.std(percentile_list(np.array(shifts),0.1), ddof=1)
    else:
        shift = 0
    twidths = np.array(twidths)
    printarrayformat = ['%1.1i', '%4.2f', '%4.2f']
    printarray = np.array([ range(sci_tr_poly.shape[0]), np.array(oldwidths)[:,2], twidths[:,2] ]).T
    problem_text = ''
    if problem_order != []:
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
    :param im: 2d array of float, the traces are already removed
    :return median(sim): float, median of the standard deviations from parts of the images
    :return std(sim): float, variation of the standard deviations from parts of the images
    """
    im[im == 0.0] = np.nan                                                              # Replace zeros with NaNs
    im[im == 0] = np.nan                                                              # Replace zeros with NaNs
    #plot_img_spec.plot_image(im, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'orig')
    dummy, gim, sim = bin_im(im, [10, 10])          # gim: number of eleements, sim: standard deviation
    #plot_img_spec.plot_image(dummy, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'bined_im')
    #plot_img_spec.plot_image(gim, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'bined_datapoints')
    #plot_img_spec.plot_image(sim, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'bined_std')
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

def find_shift_images(params, im, im_ref, sci_tr_poly, xlows, xhighs, widths, w_mult, cal_tr_poly, extract=True):
    """
    Finds the shift between two images by cross corellating both. The images are substracted and the total flux [sum(abs())] is calculated. At the best position, a minimum will be reached
    :param im: 2d array with the image, for which the shift should be calculated
    :param im_ref: 2d array with the reference image
    :return shift: shift of the image (float)
    """
    # Find the maximum shift by using the space between science and calibration orders
    shifts = []
    if sci_tr_poly.shape[0] == 1:   # only one order, e.g. when using live_extraction
        shift = 10 * max(widths[:,2])
    else:                           # normal case
        shifts.append( np.abs(sci_tr_poly[ :  ,0,-1] - cal_tr_poly[ :  ,0,-1]) )
        shifts.append( np.abs(sci_tr_poly[1:  ,0,-1] - cal_tr_poly[ :-1,0,-1]) )        # calibration fiber could be left or right of science fiber
        shifts.append( np.abs(sci_tr_poly[ :-1,0,-1] - cal_tr_poly[1:  ,0,-1]) )        # calibration fiber could be left or right of science fiber
        if sci_tr_poly.shape[0] >= 3:
            start = 0
            if sum(shifts[0]) != 0.0:                                   # don't fit if all values are 0
                start = 1
            for i in range(start,len(shifts)):
                poly = np.polyfit(range(len(shifts[i])), shifts[i], 2)
                shifts[i] = np.polyval(poly, range(len(shifts[i])) )
        if sum(shifts[0]) != 0.0:
            shift = min( min(shifts[0]), min(shifts[1]), min(shifts[2]) )
        else:                                                           # when science and calibration at the same position
            shift = min( min(shifts[1]), min(shifts[2]) )
    shift = min(shift, 20 * max(widths[:,2]))
    range_shifts = [-int(shift/2),int(shift/2)]
    ims1 = im.shape[1]
    #shifts = range(min(range_shifts),max(range_shifts)+1)
    fluxdiff, shifts = [], []
    oldcen, olderr = np.nan, np.nan
    if extract:                                                     # takes a bit longer
        logger('Step: checking if the current image is shifted compared to the reference frame in which the traces were searched (in cross-dispersion direction)')
    for shift in range(max(np.abs(range_shifts))+1):
        for pm in [-1, +1]:
            if extract:         # extract and find the maximum flux
                spec, good_px_mask = extract_orders(params, im, sci_tr_poly, xlows, xhighs, widths, w_mult, pm*shift, plot_tqdm=False)       
                fluxes = np.nansum(spec, axis=1)
                fluxes[np.isnan(fluxes)] = np.nanmin(fluxes)          # replace nans with the minimum flux to avoid problems caused at the borders at the image
            else:               # Use the difference of the 2 images
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

    if popt[2] > 0.3 * (max(range_shifts) - min(range_shifts)) or popt[1] < min(range_shifts) or popt[1] > max(range_shifts):
        comment = ' The calculated values (shift = {0} px, width = {1} px) were not very precise or outside of the allowed range ([{2},{3}]).'.format(round(popt[1],2), round(popt[2],2), min(range_shifts), max(range_shifts) )
        shift, width = 0.0, 0.0
    elif not extract:                                                                               # keep everything if not extracting
        shift, width, comment = round(popt[1],2), round(popt[2],2), ''
    else:                                                                                           # Find the place of the maximum flux (and redo the gaussian fit)
            logger('Step: Do the finetuning of the center and find the maximum', show=False)
            #nshifts, fluxes = [], []
            nshifts, fluxes = shifts, fluxdiff
            fshifts = np.linspace(popt[1] - 0.1*popt[2], popt[1] + 0.1*popt[2], 7)                  # Search in the central area again
            for dummy in range(5):                                                                  # repeat several times in order to get closer and closer to the maximum
                for fshift in fshifts:
                    if np.min( np.abs( np.array(nshifts) - fshift ) ) < 0.005:                  # Don't do the fit again, if it was done before
                        continue
                    spec, good_px_mask = extract_orders(params, im, sci_tr_poly, xlows, xhighs, widths, w_mult, fshift, plot_tqdm=False)
                    fluxes.append( np.mean( np.nansum(spec, axis=1) ) )
                    nshifts.append(fshift)
                sort_arr = np.argsort(fluxes)
                best_shifts = np.array(nshifts)[sort_arr][-3:]               # best 3 values
                #print shifts, fluxes, best_shifts
                fshifts = np.linspace( min(best_shifts), max(best_shifts), 7)
            # Do the gaussian fit again, although changes are less than 1%
            popt = centroid_order(shifts,fluxdiff,shifts[np.argmax(fluxdiff)],max(shifts)-min(shifts))    #a,x0,sigma,b in a*np.exp(-(x-x0)**2/(2*sigma**2))+b
            shift, width = round(popt[1],2), round(popt[2],2)
            mshift = round( best_shifts[-1], 2)                 # find the maximum
            comment = ' The center of the gauss was shifted by {0} px and the maximum flux found at a shift of {1} px.'.format(shift, mshift)
            #shift = np.mean([shift,mshift]     # if not symetrical traces then, a big shift
            shift = mshift                                      # use the maximum flux as shift
    logger('Info: The shift between this frame and the reference frame is {0} px. The gaussian width is {1} px. A shift between {2} and {3} pixel was tested.{4}'.format(shift, width, min(shifts), max(shifts), comment ))
    return shift

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
    w_mult = 2                                              # Fix it to 2 in order to make the shift independent of the user values
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
    fluxes = []
    logger('Step: Search for the shift of the calibration traces, compared to the science traces')
    # plot_img_spec.plot_image(im, ['savename'], pctile=0, show=True, adjust=[0.05,0.95,0.95,0.05], title='Title', autotranspose=False)
    for shift in tqdm(arcshifts, desc='Search for the shift of the calibration traces, compared to the science traces'):
        # Shift the solution -> give parameter offset; and extract the arc_spectrum with the offset
        arc_spec, good_px_mask = extract_orders(params, im, pfits, xlows, xhighs, widths, w_mult/2., shift)          # The smaller w_mult -> the better the gauss (hopefully)
        flux = np.nansum(arc_spec, axis=1)
        flux[np.isnan(flux)] = np.nanmin(flux)          # replace nans with the minimum flux
        fluxes.append(flux)
    fluxes = np.array(fluxes)
    # Find the center in each order: where is the flux the highest
    gauss, goodvalues = [], []
    label_datapoins, label_gauss, label_centroids = [],[],[]
    for order in orders:
        if cen_pos_diff != []:                          # Automatic determination of the search area
            for i in range(5):                          # If the cen_pos_diff is NaN for this order then check neightboring orders
                orderi = min(len(cen_pos_diff)-1,order+i)
                if not np.isnan(cen_pos_diff[orderi]):
                    break
                orderi = max(0,order-i)
                if not np.isnan(cen_pos_diff[orderi]):
                    break
            goodpos = [0]
            w_mult_test = w_mult*2./0.9
            while len(arcshifts[goodpos]) < widths[order,2]*w_mult*2:           # Enough data to have at least one extra order between the 2 orders
                w_mult_test *= 0.9
                goodpos = ( (np.abs(arcshifts) > widths[order,2]*w_mult_test) & (np.abs(arcshifts) < cen_pos_diff[orderi]-widths[order,2]*w_mult_test) )   # factor 2 to avoid the science fiber
                #print arcshifts[goodpos], widths[order,2]*w_mult_test, cen_pos_diff[orderi]-widths[order,2]*w_mult_test
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
            gauss.append(popts[0][1])      #shift=popt[1]
        else:
            goodvalues.append(False)
            gauss.append([0,0,0,0])
        label_datapoins.append('extracted flux in order {0}'.format('%2.2i'%order))
        label_gauss.append('fitted gauss in order {0}'.format('%2.2i'%order))
        label_centroids.append('final center in order {0}'.format('%2.2i'%order))
    
    gauss = np.array(gauss)
    if (orders[goodvalues]).shape[0] < 3:
        shift = round(np.mean(params['arcshift_side']*np.array(params['arcshift_range'])),2)
        logger('Warn: A shift to the calibration orders could not be measured. The average value of the input parameters "arcshift_range" are used instead: {0}'.format(shift))
        shifts = np.repeat([shift], orders.shape[0])
        poly, diff, min_gauss, max_gauss, used_orders = [0], [0], 0, 0, 0
    else:
        # Fit a polynomial to the center
        goodval, poly = sigma_clip(orders[goodvalues], gauss[goodvalues,1], 2, sigma, sigma, repeats=5)
        shifts = np.round(np.polyval(poly, orders),2)
        diff = shifts[goodvalues][goodval] - gauss[goodvalues,1][goodval]
        min_gauss, max_gauss = round(min(gauss[goodvalues,2][goodval]),2), round(max(gauss[goodvalues,2][goodval]),2)
        used_orders = np.sum(goodval)
    # Log the results
    arcshifts = np.repeat([arcshifts],len(orders),axis=0)
    title = 'Determining the shift of the arc orders'
    plot_gauss_data_center(arcshifts, fluxes.T, label_datapoins, gauss, label_gauss, shifts, label_centroids, params['logging_find_arc_traces'], title=title)
    logger('Info: The shift between science traces and calibration orders is between {0} and {1} px. The standard difference of the measured values to the 2nd order polynomial is {4} px, {5} orders were used. The gaussian width of the arc lines in spacial direction is between {2} and {3} px'.format(min(shifts), max(shifts), min_gauss, max_gauss, round(np.std(diff, ddof=len(poly)),2), used_orders ))
    logger('Info: The parameters of the polynomial are: {0}'.format(poly), show=False)
    return shifts
    
def identify_lines(params, im, im_short=None, im_badpx=None, im_short_badpx=None):
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
    FWHM = 4
    for order in tqdm(range(ims[0]), desc='Identify lines in the arc'):
        xarr = np.arange(ims[1])
        yarr = im[order,:]
        notnans = ~np.isnan(yarr)
        yarr1 = yarr[notnans]
        xarr1 = xarr[notnans]
        if len(yarr1) < 10:         # At the border of the chip the arc traces might be outside of the CCD
            continue
        notline_pos, pfit = sigma_clip(xarr1, yarr1, 12, 10, 2.2, repeats=20)  #orders, sigma low, sigma high ; tested that 2.5 is too high for some lines (UNe), If lines are not identified correctly, the problem isn't here, probably
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
            range_arr = range(max(0,pos_real-FWHM*2),min(ims[1],pos_real+FWHM*2+1))
            # check if line is saturated and if so use the short exposure time
            y_data = yarr[range_arr]
            if im_short is not None and im_badpx is not None and not (im == im_short).all():
                #print order, pos_max, np.sum(im_badpx[order,range_arr]==0), np.sum(im_badpx[order,range_arr]==0.1), np.sum(im_badpx[order,range_arr]==0.2), np.sum(im_badpx[order,range_arr]==1), np.sum(im_badpx[order,range_arr]!=-1)
                if 0.1 in im_badpx[order,range_arr]:
                    if im_short_badpx is not None:
                        if 0.1 in im_short_badpx[order,range_arr]:          # ignore saturated lines in the short exposure
                            continue
                    y_data = im_short[order,range_arr]
                    if np.sum(already_found[range_arr]) == 0:           # Only print once in the area
                        print('Used the short exposure time to find the centroid of the line in order {1} @ px {0}'.format(pos_real, order))
            # Fit a double Gauss 
            #popts = oneD_blended_gauss
            # fit the gauss, pos_real is the expected center
            popt = centroid_order(xarr[range_arr],y_data,pos_real,FWHM*2, significance=2, blended_gauss=False)    #a,x0,sigma,b in a*np.exp(-(x-x0)**2/(2*sigma**2))+b
            #if order == 0 and pos_real>100 and pos_real<13150:
            #    print xarr[range_arr],y_data,pos_real, popt
            if popt[1] == 0 or popt[2] > FWHM*1.5 or int(round(popt[1])) >= ims[1] or int(round(popt[1])) < 0:
                #print 'not used',order, pos_real,popt
                already_found[max(0,pos_real-1):min(len(already_found),pos_real+1)+1] = 1
                continue
            #print int(round(popt[1])),round(popt[1]), popt[1]
            if already_found[int(round(popt[1]))] == 1:
                #print 'already found',order, pos_real,popt
                continue
            #print 'used',order, pos_real,popt
            #print order,popt, pos_real #, yarr[pos_real-5:pos_real+6].astype(int)
            #yarr -= oneD_gauss(xarr,popt)  #line is not really a gauss
            already_found[max(0,int(popt[1]-popt[2]*3)):min(len(yarr),int(popt[1]+popt[2]*3))+1] = 1
            lines.append([order, popt[1], popt[2], popt[0] ])

    lines = np.array(lines)
    # Remove the lines which are too wide/narrow
    good_values, pfit = sigma_clip(lines[:,2]*0, lines[:,2], 0, 3, 3)  #orders, sigma low, sigma high
    lines = lines[good_values, :]
    
    return lines

def read_reference_catalog(filename, wavelength_muliplier, arc_lines):
    """
    Reads the reference catalogue from a file and extracts the lines which should be used
    :param filename: text, filename of the catalogue file. The file needs to consist of 3 columns: wavelength of the line, 
                            strength/heigt of the line [can be empty or followed by text (as from the NIST database)], and name of the line
    :param wavelength_muliplier: float, if the resolution is not given in Angstrom, than it can be converted into Angstrom with this factor
    :param arc_lines: list of text, subsection of the arc lines, which should be extracted
    :return reference_catalog: 2d array with one entry for each line. Each entry contains the wavelength, the intensity of the line, and the index in the catalogue
    :return reference_names: list with same length as reference_catalog, name of each line
    """
    for i in range(len(arc_lines)):
        arc_lines[i] = arc_lines[i].replace(' ','')
    reference_catalog = []
    reference_names = []
    if os.path.isfile(filename) == False:
        logger('Error: file {0} does not exist'.format(filename) )
    file = open(filename, 'r')
    for line in file:
        line = line[:-1].split('\t')
        if len(line) < 3:
            continue
        if line[2].replace(' ','') not in arc_lines:        # if Ar IV
            continue
        try:
            line[0] = float(line[0])*wavelength_muliplier   # wavelength in Angstrom
        except:
            logger('Warn: Wavelength {0} of line {1} in file {2} cannot be transformed to float'.format(line[0], line, filename))
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
            logger('Warn: Line intensity cannot be smaller than 0.0, check line {0} in {1}'.format(line, filename))
            line[1] = 1.
        if line[1] == 0:
            line[1] = 0.1
        reference_names.append(line[2])
        line = line[0:2]
        if reference_names[-1].find('Ar I') == 0:
            line[1] *= 50
        line.append(len(reference_catalog))
        reference_catalog.append(line)   # wavelength in Angstrom, relative intensity of the line, index of line in reference_names
    file.close()
    if reference_catalog == []:
        logger('Error: no reference lines found in {0} for the requested lines {1}'.format(filename, arc_lines))
    reference_catalog = np.array(reference_catalog)
    arcs = reference_catalog.shape
    # Remove the faintest lines, if too many lines are in the catalogue
    if arcs[0] > 100000:
        breakvalue = np.percentile(reference_catalog[:,1],min(90.,arcs[0]/120.))
        keep = np.logical_or(reference_catalog[:,1] >= breakvalue , reference_catalog[:,1] == 1)
        for i in range(arcs[0])[::-1]:
            if keep[i] == False:
                del reference_names[i]                      # remove the names
        reference_catalog = reference_catalog[keep,:]       # remove the wavelengths, intensities
        logger('Info The faintest {0} of {1} entries in the arc reference file {2} will not be used '.format(arcs[0]-reference_catalog.shape[0], arcs[0], filename ))
    return reference_catalog, reference_names

def shift_wavelength_solution(params, aspectra, wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, xlows, xhighs, obsdate_float, sci_tr_poly, cal_tr_poly, objname, maxshift=2.5, in_shift=0, fib='cal'):
    """
    Determines the pixelshift between the current arc lines and the wavelength solution
    Two ways for a wavelength shift can happen:
        1: The CCD moves -> all ThAr lines are moved by the same pixel distance -> this was implemented from the beginning -> This implemented
        2: The light hits the grating at a different position and passes through different are of the lens -> Pixelshift depends on wavelength, pixel
            maybe replace by a cross correlation between arc spectra: x1, y1 from wavelength solution, x2, y2 from the current file -> y2'(x) = y1(x+dx)*a+b so that y2 and y2' match best
        3: Re-fit the wavelength solution each time -> replace the lines that have been identified in this aspectra
            if not enough lines available -> linear shift (method 1)
    :param aspectra: spectrum of the reference orders
    :param wavelength_solution: 2d array of floats, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    :param wavelength_solution_arclines: 2d array of floats with one line for each order. Each order contains the wavelengths of the identified reference lines, 
                                        sorted by the brightest to faintest (in spectrum from which the solution is derived) and 0 to make into an array
    :param reference_catalog: 2d array of floats with one entry for each line. Each entry contains the wavelength and the intensity of the line
    :param reference_names: list of strings with same length as reference_catalog, name of each line
    :param xlows: list of floats, length same as number of orders, the lowest x pixel (wavelength direction) used in each order
    :param xhighs: list of floats, length same as number of orders, the highest x pixel (wavelength direction) used in each order
    :param in_shift: integer, gives an offset of where to expect the lines (positive -> means line in aspectra is to the right, compared to wavelength solution)
    :param fib: for which fiber is wavelength solution
    :return wavelength_solution: 2d array of floats, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    """
    #logger('Step: Finding the wavelength shift for this exposure')
    if np.max(wavelength_solution[:,-1]) < 100:                 # pseudo solution, wavelength of the central pixel is < 100 Angstrom
        return copy.deepcopy(wavelength_solution), 0.               # No shift is necessary
    if np.nansum(np.abs(sci_tr_poly - cal_tr_poly)) == 0.0 and np.nansum(aspectra) == 0:         # science and calibration traces are at the same position and it's not the calibration spectrum
        wavelength_solution_shift, shift = shift_wavelength_solution_times(params, wavelength_solution, obsdate_float, objname)
        return wavelength_solution_shift, shift
    
    FWHM = 3.5
    ratio_lines_identified = 0.15       # if more than ratio_lines_identified of the checked_arc_lines has been identified, then a sigma_clipping will be applied. If less than this number of lines remain after sigma clipping, then the assumption is, that the calibration fiber wasn't used and therefore no wavelength shift is applied
    # In each order get the approx pixel of each reference line, fit a gaussian against the position, calculate the wavelength for the gaussian center, compare the wavelength of the line center with the reference line wavelength
    if in_shift == 0 and params['wavelength_solution_type'] == 'sci-fiber' and fib == 'cal':
        in_shift = params['master_shift']              # if the wavelength solution is from the science fiber, then all the lines calibration fiber will be shifted
    if in_shift == 0 and params['wavelength_solution_type'] == 'cal-fiber' and fib == 'sci':
        in_shift = params['master_shift']                   # if the wavelength solution is from the calibration fiber and emission spectrum from the science fiber, then the shift is necessary
    in_shift_int = int(round(in_shift))
    #print 'input_shift', in_shift, in_shift_int
    
    aspectra = np.array(aspectra)
    ass = aspectra.shape
    shifts = []
    checked_arc_lines = 0
    for order_index in range(wavelength_solution.shape[0]):
        xarr = np.arange(ass[1])
        warr = np.polyval(wavelength_solution[order_index,2:], xarr-wavelength_solution[order_index,1])                 # Wavelength of each pixel in the order
        for arcline in wavelength_solution_arclines[order_index]:
            ref_line_index = np.argmin(np.abs( reference_catalog[:,0] - arcline ))
            if abs( reference_catalog[ref_line_index,0] - arcline ) > 0.0001:           # check that it is in the reference catalog, allowing for uncertainty; catches also zeros used to fill the array
                continue                                                                # it's not in the reference catalog
            diff = np.abs(warr-reference_catalog[ref_line_index,0])                     # Diff in wavelengths
            if min(diff) <= wavelength_solution[order_index,-2]*maxshift:               # Distance should be less than 2.5 px
                checked_arc_lines += 1
                pos = np.argmin(diff)                                               # Position of the line in the array
                range_arr = range( max(0,pos+in_shift_int-int(FWHM*3)), min(ass[1],pos+in_shift_int+int(FWHM*3)+1) )     # Range to search for the arc line
                #print range_arr, xarr, aspectra.shape, order_index
                popt = centroid_order(xarr[range_arr],aspectra[order_index,range_arr], pos+in_shift_int, FWHM*3, significance=3.5)    #a,x0,sigma,b in a*np.exp(-(x-x0)**2/(2*sigma**2))+b
                if popt[1] == 0 or popt[2] > FWHM*1.5:
                    #if order_index == 50:
                    #    print 'not used',order_index, pos,popt, reference_catalog[ref_line_index,0]
                    #    plot_img_spec.plot_points([xarr[range_arr]],[aspectra[order_index,range_arr]],[str(pos)],'path',show=True, x_title='Pixel', y_title='Flux')
                    continue
                #if order_index == 50:
                #    print 'used',order_index, pos,popt, reference_catalog[ref_line_index,0]
                #    x = xarr[range_arr]
                #    plot_img_spec.plot_spectra(np.array([x, x]),np.array([aspectra[order_index,range_arr], oneD_gauss(x,popt)]),['data','fit'], ['savepaths'], True, [0.06,0.92,0.95,0.06, 1.0,1.01])
                xarr_fine = np.arange(xarr[pos]-2, xarr[pos]+2.01, 0.01)                                                        # fine check of where the reference line is lockated
                warr_fine = np.polyval(wavelength_solution[order_index,2:], xarr_fine-wavelength_solution[order_index,1])       # Wavelength of the fine array
                diff = np.abs(warr_fine-reference_catalog[ref_line_index,0])
                pos_fine = xarr_fine[np.argmin(diff)]
                #print 'used',order_index, pos, reference_catalog[ref_line_index,0], popt[1], pos_fine, popt[1] - (pos_fine+in_shift), maxshift
                if np.abs(popt[1] - (pos_fine+in_shift)) <= maxshift:                                   # Miximal 2.5 px shift
                    # Index of order, Index of reference line, wavelength of reference line, px of arc line, \n wavelength of arc line, \n height of arc line, sigma/width of gauss, original pixel position of reference line
                    shifts.append([ order_index, ref_line_index, reference_catalog[ref_line_index,0], popt[1], \
                                    np.polyval(wavelength_solution[order_index,2:], popt[1]-wavelength_solution[order_index,1]), \
                                    popt[0], popt[2], pos_fine ])
    shifts = np.array(shifts)
    if len(shifts) >= ratio_lines_identified * checked_arc_lines and len(shifts) != 0:    
        # Using the wavelenght shift, but this is wrong
        #good_values, pfit = sigma_clip(shifts[:,3]*0, shifts[:,4]-shifts[:,2], 0, 3, 3, repeats=20)  #orders, sigma low, sigma high
        # Using the pixel shift
        good_values, pfit = sigma_clip(shifts[:,3]*0, shifts[:,3]-shifts[:,7], 0, 3, 3, repeats=20)  #orders, sigma low, sigma high
        shifts = shifts[good_values,:]
    shift_avg, shift_std, width_avg, width_std = in_shift,0,0,0
    if len(shifts) >= ratio_lines_identified * checked_arc_lines and len(shifts) != 0:                          # Only if enough lines have been detected
        #shift_avg, shift_std = np.mean( (shifts[:,4]-shifts[:,2]) / shifts[:,2] ), np.std( (shifts[:,4]-shifts[:,2]) / shifts[:,2], ddof=1)    # Will vary along the CCD
        shift_avg, shift_std = np.mean( (shifts[:,3]-shifts[:,7]) ), np.std( (shifts[:,3]-shifts[:,7]), ddof=1)
        width_avg, width_std = np.mean(shifts[:,6]), np.std(shifts[:,6], ddof=1)
        
        # Save the shift for later use
        if params['extract_wavecal']:            # This values are not correct, as later the shift has to be applied
            add_text_to_file('{0}\t{1}\t{2}\t{3}'.format(obsdate_float, shift_avg-in_shift, shift_std, fib), params['master_wavelengths_shift_filename'] )
        
    logger('Info: The shift between the lines used in the wavelength solution and the current calibration spectrum of file {9} is {0} +- {1} px ({8} km/s). {2} reference lines have been used, {7} reference lines have been tested. The arc lines have a Gaussian width of {3} +- {4} px, which corresponds to a FWHM of {5} +- {6} px'\
                .format(round(shift_avg,4), round(shift_std,4), shifts.shape[0], round(width_avg,3), \
                        round(width_std,3), round(width_avg*2.35482,3), round(width_std*2.35482,3), checked_arc_lines,
                        round(shift_avg*np.median(wavelength_solution[:,-2]/wavelength_solution[:,-1])*Constants.c/1000.,4), objname ))
    if len(shifts) >= ratio_lines_identified * checked_arc_lines and len(shifts) != 0:                          # Statistics only if enough lines were detected
        statistics_arc_reference_lines(shifts, [0,1,6,2], reference_names, wavelength_solution, xlows, xhighs, show=False)
    # correction in the other side of the shift

    wavelength_solution_new = copy.deepcopy(wavelength_solution)
    #wavelength_solution_new[:,1] -= shift_avg                       # shift the central pixel, - sign is right, tested before 19/9/2018
    wavelength_solution_new[:,1] += shift_avg                       # shift the central pixel, + sign is right, tested on 19/9/2018
    
    # In case of pixel shift available -> linear interpolation of pixel shift
    if not params['extract_wavecal']:         # science and calibration traces are at the same position and it's not the calibration spectrum
        wavelength_solution_new, shift_stored = shift_wavelength_solution_times(params, wavelength_solution_new, obsdate_float, objname)
        shift_avg += shift_stored
    
    if False:
        shift_avg -= in_shift
        logger('Info: Corrected for input shift. The shift of the currect spectrum is {0}'.format(round(shift_avg,4) ))
    #print 'return_shift', shift_avg
    return wavelength_solution_new, shift_avg
    
    """ old and wrong
    wavelength_solution_new = []
    for wls in wavelength_solution:
        #wls[-1] -= shift_avg        #This cause the wavelength_solution to change globally, even if copy.copy(wavelength_solution) is used in the for loop or in the call of the procedure -> copy.deepcopy might solve this
        wls_neu = wls[0:-1]
        wls_neu = np.append(wls_neu, wls[-1] - shift_avg * wls[-1])     # Alternatively shift the central pixel wls[1] ?
        wavelength_solution_new.append(wls_neu)
    
    return np.array(wavelength_solution_new)"""

def shift_wavelength_solution_times(params, wavelength_solution, obsdate_float, objname):
    """
    In case no calibration spectrum was taken at the same time as the science spectra use the stored information to aply a shift in the lines
    !!! shift_avg = np.average( shifts[:,1], weights=weight )  -> Tested, will work fine for obsdate in range(all_shifts), but will fail for extrapolation. Maybe it's necessary to replace this by a linear fit?
    :return wavelength_solution: 2d array of floats, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    """
    all_shifts = read_text_file(params['master_wavelengths_shift_filename'], no_empty_lines=True)
    all_shifts = convert_readfile(all_shifts, [float, float, float, str], delimiter='\t', replaces=['\n'])       # obsdate_float, shift_avg, shift_std, can contain duplicate obsdate_float (last one is the reliable one)
    if len(all_shifts) == 0:
        logger('Warn: No pixel-shift for the wavelength solution is available. Please re-run prepare_file_list.py and asign "w" to the emission line spectra.')
        return copy.deepcopy(wavelength_solution), 0.               # No shift is possible
    all_shifts_str = np.array(all_shifts)                           # Becomes an array of strings
    # Remove the double entries, e.g. if better shifts were added later
    all_shifts_str_n = []
    for i in range(all_shifts_str.shape[0])[::-1]:
        if np.sum( ( all_shifts_str[i:,0] == all_shifts_str[i,0] ) & ( all_shifts_str[i:,3] == all_shifts_str[i,3] ) ) == 1:    # Find only exactly this entry with same obsdate and same fiber
            all_shifts_str_n.append(all_shifts_str[i,:])
    all_shifts_str = np.array(all_shifts_str_n)                     # Only contains the latest entries, newest entry first
    # select the right data, depending on params['two_solutions']
    shifts = []
    if not params['two_solutions']:                                 # single soution: use all lines
        all_shifts = all_shifts_str[:,0:3].astype(float)
    else:                                                           # two wavelength solutions -> compare the two
        all_shifts_cal = all_shifts_str[all_shifts_str[:,3]=='cal',0:3].astype(float)
        all_shifts_sci = all_shifts_str[all_shifts_str[:,3]=='sci',0:3].astype(float)
        if len(all_shifts_cal)*len(all_shifts_sci) == 0:
            logger('Warn: Offset between the fibers cannot be determined as the offset is only available for one fiber. '+\
                    'Please re-run prepare_file_list.py and asign "w" and "w2" to the emission line spectra of the calibration and science fiber, respectively.')
            return copy.deepcopy(wavelength_solution), 0.               # No shift is possible
        all_shifts = []
        for entry in all_shifts_sci:
            diff = np.abs(entry[0] - all_shifts_cal[:,0])           # difference in obsdate
            if np.min(diff) > 1 * 3600:                             # difference less than one hour
                continue
            index = np.argmin(diff)                                    # the closest entry
            all_shifts.append([ np.mean([entry[0], all_shifts_cal[index,0]]), entry[1]-all_shifts_cal[index,1], np.sqrt(entry[1]**2+all_shifts_cal[index,1]**2) ])
            # mean of obsdate, diff between two fibers (sci-cal), error of the shift
        all_shifts = np.array(all_shifts)
        if params['wavelength_solution_type'] == 'sci-fiber':       # Test on 20190307 indicates that it should be sci-fiber
            all_shifts[:,1] = -1 * all_shifts[:,1]
    # Find the closest entries in time
    for diff_array in [obsdate_float - all_shifts[:,0], all_shifts[:,0] - obsdate_float]:           # difference in two different directions
        diff_array[diff_array < 0] = diff_array.max()                                               # Only values before/after obsdate_float
        index = np.where( diff_array == diff_array.min() )[0][-1]                                   # Get the index in shifts for the last minimum
        shifts.append(all_shifts[index,:])
    shifts = np.array(shifts)
    weight = np.abs(shifts[:,0] - obsdate_float)
    if np.nansum(weight) == 0:                                                                      # When using the border position
        weight += 1.
    weight /= (np.nansum(weight)+0.0)
    weight = 1 - weight
    shift_avg = np.average( shifts[:,1], weights=weight )           # Tested, will work fine for obsdate in range(all_shifts), but will fail for extrapolation. Maybe it's necessary to replace this by a linear fit?
    logger('Info: The shift between the wavelength solution and the current file {7} (center of exposure is {1}) is {0} px ({6} km/s). The stored shifts {2} ({3}) and {4} ({5}) from file {8} were used.'.format(\
                        round(shift_avg,4), datetime.datetime.utcfromtimestamp(obsdate_float).strftime('%Y-%m-%d %H:%M:%S'), 
                        round(shifts[0,1],4), datetime.datetime.utcfromtimestamp(shifts[0,0]).strftime('%Y-%m-%d %H:%M:%S'),
                        round(shifts[-1,1],4), datetime.datetime.utcfromtimestamp(shifts[-1,0]).strftime('%Y-%m-%d %H:%M:%S'),
                        round(shift_avg*np.median(wavelength_solution[:,-2]/wavelength_solution[:,-1])*Constants.c/1000.,4), objname, params['master_wavelengths_shift_filename'] ) )

    wavelength_solution_new = copy.deepcopy(wavelength_solution)
    #wavelength_solution_new[:,1] -= shift_avg                       # shift the central pixel, - sign is right, tested before 19/9/2018
    wavelength_solution_new[:,1] += shift_avg                       # shift the central pixel, + sign is right, tested on 19/9/2018
    
    return wavelength_solution_new, shift_avg

def create_wavelengths_from_solution(wavelength_solution, spectra):
    """
    Converts the wavelength solution into a 2d array with the wavelength for each pixel and order
    :param wavelength_solution: 2d array of floats, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    :param spectra: 2d array of floats, spectrum, only needed for the shape of the wavelengths array
    :return wavelengths: 2d array of floats, same dimensions as spectra, contains the wavelength for each pixel
    """
    wavelengths = []
    xarr = np.arange(np.array(spectra).shape[1])
    for wls in wavelength_solution:
        wavelengths.append(np.polyval(wls[2:], xarr - wls[1]))
    return np.array(wavelengths)

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

def trace_orders(params, im_sflat, im_sflat_head):
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
        params['binx'], params['biny'] = params['bin_search_apertures']
        sim_sflat, dummy, dummy = bin_im(im_sflat, params['bin_search_apertures'] )
        save_im_fits(params, sim_sflat, im_sflat_head, params['logging_trace1_binned'])
     
        # Search for orders in the small image
        polyfits, xlows, xhighs = find_trace_orders(params, sim_sflat)
        polyfits = unify_pord(polyfits)
        # save parameters of the polynoms into a fitsfile (from Neil)
        save_fits_width(polyfits, xlows, xhighs, [], params['logging_traces_binned'])
        plot_traces_over_image(im_sflat, params['logging_traces_im_binned'], polyfits, xlows, xhighs)
        
    # retrace orders in the original image to finetune the orders
    params['binx'], params['biny'] = params['bin_adjust_apertures']
    sim_sflat, dummy, dummy = bin_im(im_sflat, params['bin_adjust_apertures'])        # Not saved
    polyfits, xlows, xhighs, widths = adjust_trace_orders(params, sim_sflat, im_sflat, polyfits, xlows, xhighs)
    if params['GUI']:
        logger('Step: Allowing user to remove orders')
        fmask = run_remove_orders_UI(np.log10(im_sflat), polyfits, xlows, xhighs, userinput=params['GUI'])
        polyfits, xlows, xhighs, widths = polyfits[fmask], xlows[fmask], xhighs[fmask], widths[fmask]
    # save parameters of the polynoms into a fitsfile (from Neil, width added)
    """
        save_fits_width(polyfits, xlows, xhighs, widths, params['master_trace_sci_filename'])
        plot_traces_over_image(im_sflat, params['logging_traces_im'], polyfits, xlows, xhighs, widths)
        data = np.insert(widths, 0, range(len(polyfits)), axis=1)
        positio = []
        for order, pfit in enumerate(polyfits[:,0,:]):                      # For the central data
            positio.append(np.polyval(pfit[1:], im_sflat.shape[0]/2 - pfit[0]))
        data = np.append(data, np.array(positio)[:,None],axis=1)
        data = np.append(data, np.array(xlows)[:,None],axis=1)
        data = np.append(data, np.array(xhighs)[:,None],axis=1)
        printarrayformat = ['%1.1i','%3.1f', '%3.1f', '%4.2f\t','%1.1i','%1.1i','%1.1i']
        logger('\t\torder\tleft\tright\tgausswidth\tpositio\tmin_tr\tmax_tr\t(positio at the center of the image)',printarrayformat=printarrayformat, printarray=data)
    """        
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
    bck_px[bad_values] = 0
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
        espectra.append( np.sqrt( pos_spectra[order,:] + bck_noise[order]**2 ) )
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
    #np.warnings.filterwarnings('ignore')
    #low_flux = ( flat_spec < minflux )          # only use pixels with enough flux, (e.g. flat_spec_norm[1] < 0.1 means flux needs to be at least 10% of median flux)
    #np.warnings.resetwarnings()
    # Replace single values which might have a bit more flux than the minflux, but are surrounded by pixels with not enough flux
    for order in range(fss[0]):                                                 # each order
        # Find the low flux in each order
        np.warnings.filterwarnings('ignore')
        low_flux[order,:] = ( flat_spec[order,:] < minflux*np.nanmax(flat_spec[order,:]) )          # only use pixels with enough flux, (e.g. flat_spec_norm[1] < 0.1 means flux needs to be at least 10% of median flux)
        np.warnings.resetwarnings()
        for i in range(1, int(fss[1]/10)):
            np.warnings.filterwarnings('ignore')
            nonnan = np.sum(~np.isnan(flat_spec[order,:i+1]))
            np.warnings.resetwarnings()
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
    np.warnings.filterwarnings('ignore')
    high_noise = (spectra[noisedataset] > maxnoise * np.nanmedian(spectra[noisedataset]) )
    np.warnings.resetwarnings()
    # Smooth the results to have fixed borders
    for i in correctdatasets:
        spectra[i, high_noise] = np.nan
    #print np.nanmedian(spectra[noisedataset]), np.sum(high_noise), np.sum(high_noise, axis=1).shape, np.sum(high_noise, axis=1)
    return spectra

def get_obsdate(params, im_head):
    """
    Get the observation date and time using the header and exposure time
    
    :return obsdate_midexp: datetime.datetime object, containing the the center of the observation (in UTC)
    :return obsdate_float: float, unix timestamp of the mid observation time (in UTC)
    :return exposure_time: float, exposure time
    """
    obsformats = ['%Y-%m-%dT%H:%M:%S.%f','%Y-%m-%dT%H:%M:%S']           # Put into the parameters?
    if 'raw_data_mid_exposure_keys' in params.keys():                   # To stay backwards compatible, can be removed a few versions after v0.4.1
        exp_fraction_keys = params['raw_data_mid_exposure_keys']
    else:
        exp_fraction_keys = ['HIERARCH ESO INS DET1 TMMEAN', 'ESO INS DET1 TMMEAN']     # HIERARCH will be removed when python reads the header
    # Get the obsdate
    obsdate = -1
    for header_key in [ 'EXO_PIPE '+params['raw_data_dateobs_keyword'], params['raw_data_dateobs_keyword'] ]:
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
    for header_key in [ 'EXO_PIPE '+params['raw_data_exptim_keyword'], params['raw_data_exptim_keyword'] ]:
        if header_key in im_head.keys():
            exposure_time = im_head[header_key]
    if exposure_time == -1:
        logger('Warn: Cannot find the raw_data_exptim_keyword = {0} in the header. Assuming 0 seconds.'.format(params['raw_data_exptim_keyword'] ))
        exposure_time = 0
        
    obsdate_midexp = obsdate + datetime.timedelta(0, fraction*exposure_time)      # days, seconds, then other fields.
    
    epoch = datetime.datetime.utcfromtimestamp(0)
    obsdate_float = (obsdate_midexp - epoch).total_seconds()                                   # (obsdate - epoch) is a timedelta
    return obsdate_midexp, obsdate_float, exposure_time, obsdate, fraction
    
def extraction_wavelengthcal(params, im, im_name, im_head, sci_tr_poly, xlows, xhighs, widths, cal_tr_poly, axlows, axhighs, awidths, wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, flat_spec_norm, im_trace, objname):
    """
    Extracts the wavelength calibration
    
    """
    shift = 0
    obsdate, obsdate_float, exposure_time, obsdate_begin, exposure_fraction = get_obsdate(params, im_head)               # in UTC, mid of the exposure
    
    #shift = find_shift_images(params, im, im_trace, sci_tr_poly, xlows, xhighs, widths, 1, cal_tr_poly)     # w_mult=1 so that the same area is covered as for the find traces
    if im_name[-8:] == '_wavecal':
        aspectra, agood_px_mask = extract_orders(params, im, cal_tr_poly, axlows, axhighs, awidths, params['arcextraction_width_multiplier'], offset=shift, var='fast', plot_tqdm=False)
        fib = 'cal'
    elif im_name[-8:] == '_wavesci':   
        aspectra, agood_px_mask = extract_orders(params, im, sci_tr_poly, xlows, xhighs, widths, params['extraction_width_multiplier'], offset=shift, var='fast', plot_tqdm=False)
        fib = 'sci'
    else:
        logger('Error: The filename does not end as expected: {0} . It should end with _wavecal or _wavesci. This is probably a programming error.'.format(im_name))
    wavelength_solution_shift, shift = shift_wavelength_solution(params, aspectra, wavelength_solution, wavelength_solution_arclines, reference_catalog, 
                                                              reference_names, xlows, xhighs, obsdate_float, sci_tr_poly, cal_tr_poly, objname, fib=fib)   # This is only for a shift of the pixel, but not for the shift of RV
    wavelengths = create_wavelengths_from_solution(wavelength_solution_shift, aspectra)
    im_head['Comment'] = 'File contains a 3d array with the following data in the form [data type, order, pixel]:'
    im_head['Comment'] = ' 0: wavelength for each order and pixel in barycentric coordinates'
    im_head['Comment'] = ' 1: spectrum of the emission line lamp'
    im_head['Comment'] = ' 2: Mask with good areas of the spectrum: 0.1=saturated_px, 0.2=badpx'
    ceres_spec = np.array([wavelengths, aspectra, agood_px_mask])
    save_multispec(ceres_spec, params['path_extraction']+im_name, im_head, bitpix=params['extracted_bitpix'])

def get_possible_object_names(filename, replacements=['_arc','arc', '_thar','thar', '_une','une']):
    """
    Analyses the filename in order to find the possible Name of the Object (for exampled stored in parameter object_file
    The filename is subsequently stripped from the "_" or "-" separated entings
    :param filename: string, filename with removed path, file ending, and \n
    :param replacements: list of strings with entries to be removed from the filename
    return obnames: list of strings with possible object names
    """
    first_entry = filename.replace('-','_').split('_')      # most likely object is without any _ and -
    obnames = [first_entry[0]]
    obname = filename + '-'                   # Have at least one run
    while obname.find('_') != -1 or obname.find('-') != -1:
        for splitter in ['-', '_']:
            obname = obname.rsplit(splitter, 1)
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
                        break
                if i == 0 and len(obnametemp) < 5:
                    continue
                if obnametemp not in obnames:
                    obnames.append(obnametemp)
                if i == 4:                            # Use the stripped filename as new filename
                    obname = obnametemp
                    
    return obnames

def extraction_steps(params, im, im_name, im_head, sci_tr_poly, xlows, xhighs, widths, cal_tr_poly, axlows, axhighs, awidths, wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, flat_spec_norm, im_trace):
    """
    Extracts the spectra and stores it in a fits file
    
    """
    if 'EXO_PIPE BCKNOISE' not in im_head.keys():        # if not already done because localbackground is in parameters
        bck_noise_std, bck_noise_var = prepare_measure_background_noise(params, im, sci_tr_poly, xlows, xhighs, widths, cal_tr_poly, axlows, axhighs, awidths)
        if np.isnan(bck_noise_var):
            bck_noise_var = -1
        im_head['HIERARCH EXO_PIPE BCKNOISE'] = round(bck_noise_std,8)
        im_head['HIERARCH EXO_PIPE BNOISVAR'] = (round(bck_noise_var,8), 'Variation of noise through image')             # Background noise variation can be very high, because some light of the traces remains
    if im_head['EXO_PIPE BCKNOISE'] <= 0 or np.isnan(im_head['EXO_PIPE BCKNOISE']):
        logger('Warn: Measured an unphysical background noise in the data: {0}. Set the noise to 1'.format(im_head['EXO_PIPE BCKNOISE']))
        im_head['HIERARCH EXO_PIPE BNOISVAR'] = (1., '1, because of unphysical measurement')
    obsdate, obsdate_float, exposure_time, obsdate_begin, exposure_fraction = get_obsdate(params, im_head)               # in UTC, mid of the exposure
    im_head['HIERARCH EXO_PIPE DATE-OBS'] = (obsdate_begin.strftime("%Y-%m-%dT%H:%M:%S.%f"), 'UTC, Begin of expose')
    im_head['HIERARCH EXO_PIPE EXP_FRAC'] = (exposure_fraction, 'Normalised mean exposure time')
    
    im_name = im_name.replace('.fits','').replace('.fit','')                # to be sure the file ending was removed
    # Object name needs to be split by '_', while numbering or exposure time needs to be split with '-'
    obname = im_name.split('/')    # get rid of the path
    im_head['HIERARCH EXO_PIPE NAME'] = (obname[-1], im_name)       # TO know later what was the original filename
    obname = obname[-1].replace('\n','')
    #not necessary anymore obname = obname.split('-')  # remove the numbering and exposure time from the filename
    #not necessary anymore obname = obname[0]              # contains only the name, e.g. ArturArc, SunArc
    obnames = get_possible_object_names(obname)
    
    shift = find_shift_images(params, im, im_trace, sci_tr_poly, xlows, xhighs, widths, 1, cal_tr_poly)     # w_mult=1 so that the same area is covered as for the find traces
    spectra, good_px_mask = extract_orders(params, im, sci_tr_poly, xlows, xhighs, widths, params['extraction_width_multiplier'], offset=shift)#, var='prec')
    orders = range(spectra.shape[0])
    if np.nansum(np.abs(sci_tr_poly - cal_tr_poly)) == 0.0:                                                 # science and calibration traces are at the same position
        aspectra, agood_px_mask = spectra*0, copy.copy(good_px_mask)
    else:
        aspectra, agood_px_mask = extract_orders(params, im, cal_tr_poly, axlows, axhighs, awidths, params['arcextraction_width_multiplier'], offset=shift, var='fast', plot_tqdm=False)
    wavelength_solution_shift, shift = shift_wavelength_solution(params, aspectra, wavelength_solution, wavelength_solution_arclines, reference_catalog, 
                                                              reference_names, xlows, xhighs, obsdate_float, sci_tr_poly, cal_tr_poly, im_name)   # This is only for a shift of the pixel, but not for the shift of RV
    wavelengths = create_wavelengths_from_solution(wavelength_solution_shift, spectra)
    wavelengths_vac = wavelength_air_to_vacuum(wavelengths)                             # change into vacuum wavelengths
    
    espectra = combine_photonnoise_readnoise(spectra, im_head['EXO_PIPE BCKNOISE'] * np.sqrt(widths[:,2]) )
    fspectra = spectra/(flat_spec_norm[2]+0.0)        # 1: extracted flat, 2: low flux removed
    # Doing a wavelength shift for the flat_spec_norm is probably not necessay, as it's only few pixel
    measure_noise_orders = 12
    measure_noise_semiwindow = 10                   # in pixel
    efspectra = measure_noise(fspectra, p_order=measure_noise_orders, semi_range=measure_noise_semiwindow)             # Noise will be high at areas wih absorption lines
    cspectra, noise_cont = normalise_continuum(fspectra, wavelengths, nc=6, semi_window=measure_noise_semiwindow, nc_noise=measure_noise_orders)      
    #logger('Warn: !!!! no blaze correction before continuum correction')
    #cspectra, noise_cont = normalise_continuum(spectra, wavelengths, nc=4, semi_window=measure_noise_semiwindow, nc_noise=measure_noise_orders)      
    # normalise_continuum measures the noise different than measure_noise
    ceres_spec = np.array([wavelengths_vac, spectra, espectra, fspectra, efspectra, cspectra, noise_cont, good_px_mask, aspectra])
    ceres_spec = clip_noise(ceres_spec)
    
    # Change the path to the object_file to the result_path, if necessary
    if not os.path.isfile(params['object_file']) and os.path.isfile( params['object_file'].replace(params['result_path'], params['raw_data_path']) ):
        #shutil.copy2( params['object_file'].replace(params['result_path'], params['raw_data_path']), params['result_path'])         # Copies the file from the raw_data_path to the result_path, keeping the metadata
        params['object_file'] = params['object_file'].replace(params['result_path'], params['raw_data_path'])
    
    # Get the baycentric velocity
    params, bcvel_baryc, mephem, obnames, im_head = get_barycent_cor(params, im_head, obnames, params['object_file'])       # obnames becomes the single entry which matched entry in params['object_file']
    im_head['HIERARCH EXO_PIPE OBJNAME'] = (obnames[0], 'Used object name')
    
    # Do the Ceres-pipeline Radial velocity analysis 
    no_RV_names = ['flat', 'tung', 'whili', 'thar', 'th_ar', 'th-ar']
    do_RV = True
    for no_RV_name in no_RV_names:
        if im_name.lower().find(no_RV_name) in [0,1,2,3,4,5]:
            do_RV = False
            break
    if np.max(wavelength_solution[:,-1]) > 100 and do_RV and os.path.exists(params['path_ceres']) == True \
            and os.path.exists(params['path_ceres']+'utils/Correlation') and os.path.exists(params['path_ceres']+'utils/GLOBALutils') \
            and os.path.exists(params['path_ceres']+'utils/OptExtract') and os.path.exists(params['path_ceres']+'utils/CCF'):           # if not pseudo-solution and not flat and necessary files exist
        """
        Ceres mask in /home/ronny/software/ceres-master/data/xc_masks/ contain only limited wavelength range (+ additional gaps):
        G2: 3751 to 6798, 2625 data points/lines
        K5: 3782 to 6798, 5129 data points/lines
        M2: 4401 to 6862, 9493 data points/lines
        Given in the mask is the lower and the higher end of the line, and the weight
        """
        RV, RVerr2, BS, BSerr = rv_analysis(params, ceres_spec, im_head, im_name, obnames[0], params['object_file'], mephem)
        # Air to Vacuum wavelength difference only causes < 5 m/s variation: https://www.as.utexas.edu/~hebe/apogee/docs/air_vacuum.pdf (no, causes 83.15 km/s shift, < 5m/s is between the different models)
        logger('Info: The radial velocity (including barycentric correction) for {0} gives: RV = {1} +- {2} km/s, Barycentric velocity = {5} km/s, and BS = {3} +- {4} km/s'.format(\
                                    im_name, round(RV+bcvel_baryc,4), round(RVerr2,4), round(BS,4), round(BSerr,4), round(bcvel_baryc,4) ))
        im_head['HIERARCH EXO_PIPE RV_BARY'] = (round(RV+bcvel_baryc),'RV including BCV in km/s (measured RV+BCV)')
        im_head['HIERARCH EXO_PIPE RV_ERR'] = round(RVerr2,4)
    
    # Correct wavelength by barycentric velocity
    wavelengths_bary = wavelengths * (1 + bcvel_baryc/(Constants.c/1000.) )
    
    im_head_bluefirst = copy.copy(im_head)
    im_head = add_specinfo_head(spectra, spectra, noise_cont, im_head)
    im_head_bluefirst = add_specinfo_head(spectra[::-1,:], spectra[::-1,:], noise_cont[::-1,:], im_head_bluefirst)
    im_head_wave, im_head_weight = copy.copy(im_head), copy.copy(im_head)
    im_head_harps_format = copy.copy(im_head_bluefirst)
    im_head_iraf_format = copy.copy(im_head)
    im_head_iraf_format_bluefirst = copy.copy(im_head_bluefirst)
    
    ## Save in a easier way
    # Single files
    im_head_wave['Comment'] = 'Contains the wavelength per order and exctracted pixel for file {0}'.format(im_name+'_extr/_blaze')
    im_head_weight['Comment'] = 'Contains the weights per order and exctracted pixel for file {0}'.format(im_name+'_extr/_blaze')
    im_head_iraf_format = wavelength_solution_iraf(params, im_head_iraf_format, wavelengths, wavelength_solution_shift, norder=params['polynom_order_traces'][-1]+2)
    im_head_iraf_format_bluefirst = wavelength_solution_iraf(params, im_head_iraf_format_bluefirst, wavelengths[::-1,:], wavelength_solution_shift[::-1,:], norder=params['polynom_order_traces'][-1]+2)
    save_multispec(spectra,                         params['path_extraction_single']+im_name+'_extr', im_head_iraf_format, bitpix=params['extracted_bitpix'])
    save_multispec(fspectra,                        params['path_extraction_single']+im_name+'_blaze', im_head_iraf_format, bitpix=params['extracted_bitpix'])
    save_multispec(wavelengths,                     params['path_extraction_single']+im_name+'_wave', im_head_wave, bitpix=params['extracted_bitpix'])
    save_multispec(good_px_mask* flat_spec_norm[1], params['path_extraction_single']+im_name+'_weight', im_head_weight, bitpix=params['extracted_bitpix'])  # Weight shouldn't be the espectra, the flat provides a smoother function
    save_multispec(spectra[::-1,:],                 params['path_extraction_single']+im_name+'_extr_bluefirst',  im_head_iraf_format_bluefirst, bitpix=params['extracted_bitpix'])
    save_multispec(fspectra[::-1,:],                params['path_extraction_single']+im_name+'_blaze_bluefirst', im_head_iraf_format_bluefirst, bitpix=params['extracted_bitpix'])
    
    # CSV file for terra
    fname = params['path_csv_terra']+obnames[0]+'/data/'+obsdate.strftime('%Y-%m-%d%H%M%S')
    os.system('rm -f {0}{1}/results/synthetic.rv'.format(params['path_csv_terra'],obnames[0]) )     # Delete the old solution, as won't be created otherwise
    save_spec_csv(cspectra, wavelengths_bary, good_px_mask, fname)
    
    # Harps format
    im_head_harps_format = wavelength_solution_harps(params, im_head_harps_format, wavelengths[::-1,:])        # [::-1,:] -> Blue orders first      # 20190509: wavelengths instead of wavelengths_bary
    serval_keys = []
    #serval_keys.append(['INSTRUME', 'HARPS',                                                    'added for Serval'])
    serval_keys.append(['INSTRUME', 'EXOHSPEC',                                                    'added for Serval'])
    serval_keys.append(['EXPTIME',  im_head_harps_format[params['raw_data_exptim_keyword']],    'Exposure time, for Serval'])
    serval_keys.append(['DATE-OBS', im_head_harps_format['EXO_PIPE DATE-OBS'],                  'UT start, for Serval'])
    serval_keys.append(['MJD-OBS',  im_head_harps_format['EXO_PIPE MJD-START'],                 'MJD start ({0})'.format(im_head_harps_format['EXO_PIPE DATE-OBS']) ])
    if 'EXO_PIPE RA' in im_head_harps_format.keys():
        serval_keys.append(['RA',       im_head_harps_format['EXO_PIPE RA'],                        'RA start, for Serval'])
    if 'EXO_PIPE DEC' in im_head_harps_format.keys():
        serval_keys.append(['DEC',      im_head_harps_format['EXO_PIPE DEC'],                       'DEC start, for Serval'])
    if 'EXO_PIPE RV_BARY' in im_head_harps_format.keys():
        serval_keys.append(['HIERARCH ESO DRS BERV',     im_head_harps_format['EXO_PIPE BCV'],      'Barycentric Earth Radial Velocity'])
    serval_keys.append(['HIERARCH ESO DPR TECH',         'ECHELLE ',        'Observation technique'])
    serval_keys.append(['HIERARCH ESO INS MODE',         'HARPS',           'Instrument mode used.'])
    serval_keys.append(['HIERARCH ESO DRS CAL LOC NBO',  spectra.shape[0],  'nb orders localised'])
    serval_keys.append(['HIERARCH ESO OBS TARG NAME',    obnames[0],        'OB target name'])
    serval_keys.append(['OBJECT',                        obnames[0],        'OB target name'])
    serval_keys.append(['HIERARCH ESO INS DET1 TMMEAN',  im_head_harps_format['EXO_PIPE EXP_FRAC'],     'Normalised mean exposure time'])
    serval_keys.append(['HIERARCH ESO INS DET2 TMMEAN',  im_head_harps_format['EXO_PIPE EXP_FRAC'],     'Normalised mean exposure time'])
    #serval_keys.append(['HIERARCH ESO DRS BLAZE FILE',   'HARPS.2007-04-03T20:57:37.400_blaze_A.fits',  'Bla'])        # not necessary
    #serval_keys.append(['HIERARCH ESO DRS DRIFT RV USED',0.,                'Used RV Drift [m/s]'])                    # not necessary
    #serval_keys.append(['',         '',        ''])
    for order in orders:
        serval_keys.append([ 'HIERARCH ESO DRS SPE EXT SN{0}'.format(order), im_head['EXO_PIPE SN_order{0}'.format('%2.2i'%order)], 'S_N order center{0}'.format(order) ])
    for [newkey, value, comment] in serval_keys:
        if newkey not in im_head_harps_format.keys():
            im_head_harps_format[newkey] = (value, comment)
    if 'COMMENT' in im_head_harps_format.keys():
        del im_head_harps_format['COMMENT']                 # Serval can't read comments
    for entry in im_head_harps_format.keys():
        #print "key, value, comment",(entry, im_head_harps_format[entry], im_head_harps_format.comments[entry])
        if im_head_harps_format.comments[entry] == '':
            im_head_harps_format.comments[entry] = '/'      # Serval can't read header keywords without comment, for NAXISj this needs to be done in save_multispec
    fname = params['path_harpsformat']+obnames[0]+'/'+'HARPS.{0}_e2ds_A.fits'.format(im_head_harps_format['EXO_PIPE DATE-OBS'][:-3])
    if not os.path.exists(fname.rsplit('/',1)[0]):
        try:
            os.makedirs(fname.rsplit('/',1)[0])
        except:
            logger('Error: Folder to save {0} does not exists and cannot be created.'.format(fname))
    #for entry in im_head_harps_format.keys():
    #    if entry.find('AXIS') != -1:
    #        print "key, value, comment",(entry, im_head_harps_format[entry], im_head_harps_format.comments[entry]) 
    save_multispec(spectra[::-1,:], fname, im_head_harps_format, bitpix=params['extracted_bitpix'])                 # [::-1,:] -> Blue orders first
    
    # Create a linearised solution for the input spectrum and the continuum corrected spectrum
    logger('Step: Linearising the spectrum (commented out)')
    #wavelenghts_lin, spectrum_lin = linearise_wavelength_spec(params, wavelength_solution_shift, spectra, method='sum', weight=espectra)
    #save_multispec([wavelenghts_lin,spectrum_lin], params['path_extraction']+im_name+'_lin', im_head)
    #wavelenghts_lin, spectrum_lin = linearise_wavelength_spec(params, wavelength_solution_shift, cspectra, method='weight', weight=espectra)
    #save_multispec([wavelenghts_lin,spectrum_lin], params['path_extraction']+im_name+'_lin_cont', im_head)
    
    # For easier plotting
    add_text_to_file(params['path_extraction']+im_name+'.fits', 'plot_files.lst')

    ceres_spec = np.array([wavelengths_bary, spectra, espectra, fspectra, efspectra, cspectra, noise_cont, good_px_mask, aspectra])        
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
    save_multispec(ceres_spec, params['path_extraction']+im_name, im_head, bitpix=params['extracted_bitpix'])

    return obnames[0]
    
def save_spec_csv(spec, wavelengths, good_px_mask, fname):
    """
    Save the spectra in a csv file to be compatible with Terra. The files need to have the form:
    ONUM,WAVELENGTH (in Angstrongs!),FLUX (arbitrary units, normalization to order each mean recommended).
    Note: all orders must of of the SAME LENGTH
    [email from Guillem Anglada, 18/10/2018 12:25
    """
    if not os.path.exists(fname.rsplit('/',1)[0]):
        try:
            os.makedirs(fname.rsplit('/',1)[0])
        except:
            logger('Error: Folder to save {0} does not exists and cannot be created.'.format(fname))
        
    specs = spec.shape
    spec_cor = spec * good_px_mask
    spec_cor[np.isnan(spec_cor)] = 0
    wave = copy.copy(wavelengths)
    wave[np.isnan(wave)] = 0
    fname = fname.replace('.csv','') + '.csv'
    file = open(fname, 'w')
    for order in range(specs[0]):
        for px in range(specs[1]):
            file.write('{0},{1},{2}\n'.format(order, wave[order,px], spec_cor[order,px]) )
    file.close()
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
                            order+1, order, int(wavelength_solution[order,0]), round(wavelengths[order,0],6), round(avg_dwave,6), ws[1], 
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
        if (np.max(np.absolute(wave_new - wave_prev), axis=None) < 1e-10):
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
    n = 1 + 8.34213e-5 + 2.406030e-2/(130-s2) + 1.5997e-4/(38.9-s2)    # applied from https://www.as.utexas.edu/~hebe/apogee/docs/air_vacuum.pdf same formular as in Ceres
    return n

def plot_traces_over_image(im, fname, pfits, xlows, xhighs, widths=[], w_mult=1, mask=[], frame=None, return_frame=False):
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
    title = 'Plot the traces in the image (log10 of image).'
    if np.mean(widths, axis=None) == 0 and np.std(widths, axis=None, ddof=1) == 0:
        title += ' The marked width (dashed lines) are shown for an extraction width multiplier of {0}.'.format(w_mult)
    plot_img_spec.plot_image(im, [], pctile=1, show=False, adjust=[0.05,0.95,0.95,0.05], title=title, return_frame=True, frame=frame, autotranspose=False, colorbar=colorbar)
    for pp, pf in enumerate(pfits):
            if mask[pp] == False:
                continue
            xarr = np.arange(xlows[pp], xhighs[pp], 1)
            yarr = np.polyval(pf[1:], xarr-pf[0])
            if len(leftfit) == 0:                        # old 2d array, use gaussian width (single value) as width
                yarrl = yarr - widths[pp][2]*w_mult
                yarrr = yarr + widths[pp][2]*w_mult
            else:
                yarrl = np.polyval(leftfit[pp, 1:], xarr - leftfit[pp, 0]) - widths[pp][0]*(w_mult-1)
                yarrr = np.polyval(rightfit[pp, 1:], xarr- rightfit[pp, 0]) + widths[pp][1]*(w_mult-1)
                yarrl, yarrr = adjust_width_orders(yarr, yarrl, yarrr, [w_mult, w_mult])              # Adjust width
            # Take into account the boundaries of the image
            yarr[yarr < 0] = 0
            yarr[yarr > ims[1]] = ims[1]
            yarrl[yarrl < 0] = 0
            yarrl[yarrl > ims[1]] = ims[1]
            yarrr[yarrr < 0] = 0
            yarrr[yarrr > ims[1]] = ims[1]

            ymid = np.median(yarr)
            xmid = np.median(xarr)
            #if pp == 24:
            #    print 'xarr[0],xarr[-1],yarr[0],yarr[-1],yarrl[0],yarrl[-1],yarrr[0],yarrr[-1],pf',xarr[0],xarr[-1],yarr[0],yarr[-1],yarrl[0],yarrl[-1],yarrr[0],yarrr[-1],pf

            frame.plot(yarr, xarr, color='r', linewidth=1)
            frame.plot(yarrl, xarr, color='r', linewidth=1, linestyle='dashed')
            frame.plot(yarrr, xarr, color='r', linewidth=1, linestyle='dashed')
            frame.text(ymid, xmid, 'Order{0}'.format(pp),
                       horizontalalignment='center', verticalalignment='center',
                       rotation=90, color='g', zorder=5)
    if return_frame:
        return frame
    fig.set_size_inches(42, 42)        
    plt.savefig(fname, bbox_inches='tight')
    plt.close()

def plot_wavelength_solution_form(fname, xlows, xhighs, wavelength_solution):
    """
    Creates a map of the wavelength solution
    :param fname: Filename to which the image is saved
    :param xlows: list, length same as number of orders, the lowest x pixel (wavelength direction) used in each order
    :param xhighs: list, length same as number of orders, the highest x pixel (wavelength direction) used in each order
    :param wavelength_solution: 2d array of floats, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    """
    step = 25
    im = np.zeros([max(xhighs)/step, len(xlows)])+np.nan
    text_order = ''
    for order in range(len(wavelength_solution)):
        xarr = np.arange(xlows[order], xhighs[order], step)
        yarr = np.polyval(wavelength_solution[order,2:], xarr-wavelength_solution[order,1])
        diff_wave = yarr[1:] - yarr[:-1]
        im[xarr[:-1]/step, order ] = diff_wave		# im[{floats} will work
        diff_wave2 = diff_wave[1:] - diff_wave[:-1]
        if np.sum(diff_wave2>1E-3) > 10 and np.sum(diff_wave2<-1E-3) > 10:     # There should'nt be any stationary points (10 points are ok, a change of less than 1E-3 as well)
            text_order += '{0},'.format(order)
    if text_order != '':
        logger('Warn: The wavelength solution for oder(s) {0} contain(s) at least one stationary point. Please check {1}'.format(text_order[:-1], fname))
    colorbar = True  
    title = 'Plot the wavelength difference between every {0} pixel (Angstrom/{0}px)'.format(step)
    axis_name = ['Order', 'Position in Dispersion direction [{0} px]'.format(step), 'wavelength difference [Angstrom/{0}px]'.format(step)]
    plot_img_spec.plot_image(im, [fname], pctile=0, show=False, adjust=[0.05,0.95,0.95,0.05], title=title, autotranspose=False, colorbar=colorbar, axis_name=axis_name)

def plot_wavelength_solution_width_emmission_lines(fname, specs, arc_lines_wavelength):
    """
    Creates a map of the Gaussian width of the emmission lines
    :param fname: Filename to which the image is saved
    :param specs: 1d list with two integers: shape of the extracted spectrum, first entry gives the number of orders and the second the number of pixel
    :param wavelength_solution: 2d array of floats, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    """
    step = 20
    gauss_ima = np.empty((specs[0], int(specs[1]/step)+1)) * np.nan
    for order in range(specs[0]):
        for i,px in enumerate(range(0, specs[1], step)):
            good_values = (arc_lines_wavelength[:,0] == order) & (arc_lines_wavelength[:,1] >= px) & (arc_lines_wavelength[:,1] < px+step)
            if np.sum(good_values) > 0:
                gauss_ima[order,i] = np.median(arc_lines_wavelength[good_values,9])
    gauss_ima[np.isnan(gauss_ima)] = np.nanmax(gauss_ima)+1
    title = 'Plot of the Gaussian width every {0} pixel (white shows areas with no data available)'.format(step)
    axis_name = ['Position in Dispersion direction [{0} px]'.format(step), 'Order', 'Gaussian width of the emission lines [px] (white shows areas with no data available)'.format(step)]
    plot_img_spec.plot_image(gauss_ima, [fname], pctile=0, show=False, adjust=[0.05,0.95,0.95,0.05], title=title, autotranspose=False, colorbar=True, axis_name=axis_name)

def plot_wavelength_solution_spectrum(spec1, spec2, fname, wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, plot_log=False):
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
    
    spec1, spec2 = copy.copy(spec1), copy.copy(spec2)
    if plot_log:
        y_title = 'log10 ( {0} )'.format(y_title)
        if np.nanmin(spec1) < 1:
            spec1 += 1 - np.nanmin(spec1)       # minus to cancel negative minimum
        if np.nanmin(spec2) < 1:
            spec2 += 1 - np.nanmin(spec2)       # minus to cancel negative minimum
        spec1 = np.log10(spec1)
        spec2 = np.log10(spec2)
    reference_catalog = np.array(sorted(reference_catalog, key=operator.itemgetter(1), reverse=True ))          # Sort by intensity in reverse order
    wavelengths = create_wavelengths_from_solution(wavelength_solution, spec1)
    
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
            x_data = [wavelengths[order,:], wavelengths[order,:]]
            y_data = [spec1[order,:], spec2[order,:]]
            color = ['tab:blue','tab:orange']
            
            label = copy.deepcopy(labels)
            if np.min(wavelengths) > 1000:      # real wavelength solution, not pseudo solution
                if order < wavelength_solution.shape[0]-1:
                    x_data.insert(0, wavelengths[order+1,:])
                    y_data.insert(0, spec1[order+1,:])
                    label.insert(0, 'next. order')
                    color.insert(0,'lightgrey')
                if order > 0:
                    x_data.insert(0, wavelengths[order-1,:])
                    y_data.insert(0, spec1[order-1,:])
                    label.insert(0, 'prev. order')
                    color.insert(0,'silver')
            plot_img_spec.plot_points(x_data, y_data, label, [], show=False, adjust=[0.12,0.88,0.85,0.12, 1.0,1.01], title=title, 
                                      return_frame=True, frame=frame, x_title=x_title, y_title=y_title, linestyle="-", marker="",
                                      color=color)
            axes = plt.gca()
            #ymin, ymax = axes.get_ylim()
            #y_range = ymax - ymin
            #xmin, xmax = axes.get_xlim()        # wavelengths contains also values for nans
            # Redifine the axis using only the current order
            good_data = ~np.isnan(spec1[order,:])
            xmin, xmax = np.nanmin(wavelengths[order,good_data]), np.nanmax(wavelengths[order,good_data])
            ymin1, ymax1 = np.nanmin(spec1[order,good_data]), np.nanmax(spec1[order,good_data])
            ymin2, ymax2 = np.nanmin(spec2[order,good_data]), np.nanmax(spec2[order,good_data])
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
        d['Author'] = 'EXOhSPEC pipeline'
        d['Subject'] = ''
        d['Keywords'] = 'Spectra atlas EXOhSPEC pipeline'
        d['CreationDate'] = datetime.datetime.today()
        d['ModDate'] = datetime.datetime.today()
    
def plot_wavelength_solution_image(im, fname, pfits, xlows, xhighs, wavelength_solution, wavelength_solution_arclines, reference_catalog_full):
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


            
""" not used anymore def find_bck_fit(im, im_bck_px, p_orders, GUI=True):         # Needs to be rewritten using 2d image fit.
    ims = im.shape
    im_bck = []
    xarr = np.array(range(ims[1]))
    for i in range(ims[0]):           #crossing the orders
        bck = im_bck_px[i,:] == 1
        xarr1 = xarr[bck]
        yarr = im[i,bck]
        sigmaclip = (yarr*0 == 0)
        for j in range(3):
            pbck = polyfit_adjust_order(xarr1[sigmaclip], yarr[sigmaclip], p_orders[1])
            difffit = yarr - np.polyval(pbck, xarr1)
            stddiff = np.std(difffit, ddof=p_orders[1]+1)
            sigmaclip = np.abs(difffit) <= 3*stddiff
        im_bck.append(np.polyval(pbck, xarr))
        #if i/100 == i/100.:
        #if i >1510 and i< 1600 and i/10 == i/10.:
        #    plot_img_spec.plot_spectra(np.array([xarr,xarr]),np.array([(im*im_bck_px)[i,:],np.polyval(pbck, xarr)]),['data','fit'], ['savepaths'], True, [0.06,0.92,0.95,0.06, 1.0,1.01], 'line {0}'.format(i))
    im_bck = np.array(im_bck)
    if GUI == True:
        plot_img_spec.plot_image(im_bck, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'fit perpendicular to traces')
    im_bck1 = []
    yarr = np.array(range(ims[0]))
    for i in range(ims[1]):         #crossing the fits
        pbck = polyfit_adjust_order(yarr, im_bck[:,i], p_orders[0])
        im_bck1.append(np.polyval(pbck, yarr))
    im_bck1 = np.array(im_bck1).T
    if GUI == True:
        plot_img_spec.plot_image(im_bck1, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'Background map')
    return im_bck1
"""

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
                        title='Determine area for the background', widgets=widgets,
                        widgetprops=wprops)
    
    gui3.master.mainloop()
    
    params['width_multiplier'] = float(gui3.data['width_multiplier'])

    im_bck_px = bck_px(im_orig, pfits, xlows, xhighs, widths, params['width_multiplier'])
    #gui3.destroy 
    plt.close()
    return im_bck_px, params

def shift_orders_UI(im, in_shift, pfits, xlows, xhighs, widths):
    w_mult = 1
    fig, frame = plt.subplots(1, 1)
    def plot(frame, im, pfits, xlows, xhighs, widths, w_mult):
        frame.clear()
        try:
            w_mult = float(gui3.data['width_multiplier'])
        except:
            w_mult = w_mult
        try:
            shift = float(gui3.data['shift'])
        except:
            shift = in_shift
        title = ('Determine shifts between science and arc traces')
        plot_img_spec.plot_image(scale_image_plot(im,'log10'), ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], title, return_frame=True, frame=frame, autotranspose=False, colorbar=False)
        if widths is None:
            widths = np.repeat([0], len(pfits))
        for pp, pf in enumerate(pfits):
            xarr = np.arange(xlows[pp], xhighs[pp], 1)
            yarr = np.polyval(pf, xarr) + shift
            yarrl = yarr - widths[pp][2]*w_mult
            yarrr = yarr + widths[pp][2]*w_mult

            yarr[yarr < 0] = 0
            yarr[yarr > len(yarr)] = len(yarr)
            yarrl[yarrl < 0] = 0
            yarrl[yarrl > len(yarrl)] = len(yarrl)
            yarrr[yarrr < 0] = 0
            yarrr[yarrr > len(yarrr)] = len(yarrr)

            ymid = np.median(yarr)
            xmid = np.median(xarr)

            frame.plot(yarr, xarr, color='r', linewidth=1)
            frame.plot(yarrl, xarr, color='r', linewidth=1, linestyle='dashed')
            frame.plot(yarrr, xarr, color='r', linewidth=1, linestyle='dashed')
            frame.text(ymid, xmid, 'Order{0}'.format(pp),
                       horizontalalignment='center', verticalalignment='center',
                       rotation=90, color='g', zorder=5)
        #frame.set_title(title)
    # get kwargs
    pkwargs = dict(frame=frame, im = im, pfits = pfits, xlows = xlows, xhighs = xhighs, widths = widths, w_mult = w_mult)
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
            value = float(xs)
            return True, value
        except:
            return False, ('Error, input must consist of a float')
    # define widgets
    widgets = dict()
    widgets['shift'] = dict(label='Shift science -> arc',
                                comment='Shift between science traces and arc traces' ,
                                kind='TextEntry', minval=0, maxval=None,
                                fmt=str, start=str(in_shift), valid_function=vfunc,
                                width=10)
    widgets['width_multiplier'] = dict(label='Multiplier to the measured width',
                                comment='The width for extractingthe arc traces,\ncompared to the width used for\nextraction of the science traces' ,
                                kind='TextEntry', minval=0, maxval=None,
                                fmt=str, start=str(w_mult), valid_function=vfunc,
                                width=10)

    widgets['accept'] = dict(label='Accept Map', kind='ExitButton', position=Tk.BOTTOM)
    widgets['update'] = dict(label='Update', kind='UpdatePlot', position=Tk.BOTTOM)

    wprops = dict(orientation='v', position=Tk.RIGHT)

    gui3 = tkc.TkCanvas(figure=fig, ax=frame, func=plot, kwargs=pkwargs,
                        title='Determine shifts between science and arc traces', widgets=widgets,
                        widgetprops=wprops)
    
    gui3.master.mainloop()
    
    shift = float(gui3.data['shift'])
    w_mult = float(gui3.data['width_multiplier'])
    plt.close()
    return shift, w_mult

def correlate_UI(im, order, arc_settings, reference_catalog, reference_names, adjust=[0.07,0.93,0.94,0.06, 1.0,1.01]):
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
            return False, ('Error, input must be integer')
            
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
                        title='Plot spectral orders', widgets=widgets,
                        widgetprops=wprops)
    
    gui3.master.mainloop()
    
    #print 'arc_settings',arc_settings
    return arc_settings

def correlate_px_wave_result_UI(im, arc_lines_wavelength, reference_catalog, arc_lines_px, reference_names, wavelength_solution, adjust=[0.07,0.93,0.94,0.06, 1.0,1.01]):
    global oldorder
    order = 10
    oldorder = order
    ims = im.shape
            
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
            high_limit = 100
        try:
            arc_stretch = gui3.data['arc_stretch']
        except:
            arc_stretch = 1
        if order < 0 or order > im.shape[0]-1:
            print('Warn: order outside the allowed area: 0...{0}'.format(im.shape[0]-1))
            order = oldorder
        #print 'order, arc_setting', order, arc_setting
        title = 'Order {0}: identified lines: b [px], catalog lines: r [0.1px, name and wavelength at the botom],\nidentified lines: g [px, wavelength at top]'.format(order)
        
        # Plot the extracted arc spectrum
        xarr = np.arange(len(im[order,:]))
        yarr = im[order,:]
        notnans = ~np.isnan(yarr)
        yarr = yarr[notnans]
        xarr = xarr[notnans]
        yarr_max = max(0,np.percentile(yarr,high_limit))
        yarr[yarr > yarr_max] = yarr_max
        frame.plot(xarr, yarr)
        
        # Plot the identified lines in the extracted spectrum
        for arc_line in arc_lines_px[(arc_lines_px[:,0] == order),:]:
            frame.plot( [arc_line[1],arc_line[1]],[min(yarr),max(yarr)],color='b' )
        #for arc_line_order in arc_lines_px: #arc_lines_px[(arc_lines_px[:,0] == order),:]: doesn't work
        #    if arc_line_order[0] == order:
        #        for arc_line in arc_line_order[1]:
        #            frame.plot( [arc_line[0],arc_line[0]],[min(yarr),max(yarr)],color='b' )
        
        # Plot the catalog lines using the fit
        xarr_catalog = np.arange(-100,ims[1]+100, 0.1)
        if len(wavelength_solution.shape) == 1:
            xarr_wave = polynomial_value_2d(xarr_catalog-wavelength_solution[2], 1.0/(order+wavelength_solution[3]), int(wavelength_solution[0]), int(wavelength_solution[1]), wavelength_solution[4:])
        else:
            xarr_wave = np.polyval(wavelength_solution[order,2:], xarr_catalog-wavelength_solution[order,1])
        for line_index in range(reference_catalog.shape[0]):
                diff = abs(reference_catalog[line_index,0] - xarr_wave)
                if min(diff) < 1:
                    xarc = xarr_catalog[np.argmin(diff)]
                    frame.plot( [xarc,xarc],[min(yarr),min(yarr)+reference_catalog[line_index,1]*arc_stretch],color='r' )
                    frame.text(xarc, min(yarr)+reference_catalog[line_index,1]*arc_stretch, '{0} {1}'.format(reference_names[line_index],reference_catalog[line_index,0]),
                            horizontalalignment='center', verticalalignment='center', rotation=90, color='k', zorder=5)
        
                       
        # Plot the identified catalog lines
        for arc_line in arc_lines_wavelength[(arc_lines_wavelength[:,0] == order),:]:
            xarc = arc_line[1]
            frame.plot( [xarc,xarc],[min(yarr),max(yarr)],color='g' )
            frame.text(xarc, max(yarr), '{0} +/- {1}'.format(arc_line[2], round(arc_line[3],3) ), horizontalalignment='center', verticalalignment='top', rotation=90, color='k', zorder=5)

        frame.set_xlabel('x [px]', fontsize=14)
        frame.set_ylabel('flux [ADU]', fontsize=14)
        frame.set_title(title, fontsize=16)
        #frame.legend(loc='upper left', bbox_to_anchor=(adjust[4], adjust[5]))
        oldorder = order

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
            return False, ('Error, input must be integer')
            
    # define widgets
    widgets = dict()
    widgets['order'] = dict(label='Order',
                                comment='which order to show' ,
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=str, start=str(oldorder), valid_function=vfunc_int,
                                width=10)
    widgets['high_limit'] = dict(label='High limit',
                                comment='Reject highest pixels' ,
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=str, start=str(98), valid_function=vfunc_float,
                                width=10)
    widgets['arc_stretch'] = dict(label='Scale of Arc lines',
                                comment='float value' ,
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=str, start=str(1), valid_function=vfunc_float,
                                width=10)
    widgets['accept'] = dict(label='Save and Close', kind='ExitButton', position=Tk.BOTTOM)
    widgets['update'] = dict(label='Update', kind='UpdatePlot', position=Tk.BOTTOM)

    wprops = dict(orientation='v', position=Tk.RIGHT)

    gui3 = tkc.TkCanvas(figure=fig, ax=frame, func=plot, kwargs=pkwargs,
                        title='Plot spectral orders', widgets=widgets,
                        widgetprops=wprops)
    
    gui3.master.mainloop()

"""
def fit_poly_between_orders(polyfit_pl, polyfit_o, polynom_order_intertrace, sigma):       # only used by correlate_px_wavelength, which is only used by ident_arc.py
    polyfitparams = []
    available_orders = polyfit_pl.shape[1]
    good_values_all = np.repeat([True], polyfit_pl.shape[0])
    for i in range(available_orders):
        polynom_order_intertrace_i = int(min(polynom_order_intertrace, polynom_order_intertrace/(available_orders-1.0)*i+1)) # high orders need less curvature, while central wavelength needs a lot
        good_values, pfit = sigma_clip(polyfit_o, polyfit_pl[:,i], polynom_order_intertrace_i, sigma, sigma, repeats=20)  #orders, sigma low, sigma high
        good_values_all = (good_values_all & good_values)
    for i in range(available_orders):
        polynom_order_intertrace_i = int(min(polynom_order_intertrace, polynom_order_intertrace/(available_orders-1.0)*i+1)) # high orders need less curvature, while central wavelength needs a lot
        p = np.polyfit(polyfit_o[good_values_all], polyfit_pl[good_values_all,i], polynom_order_intertrace_i)
        temp=polyfit_pl[good_values_all,i]-np.polyval(p, polyfit_o[good_values_all])
        #print 'i,p, [len(temp), np.std(temp, ddof=polynom_order_intertrace_i+1),min(temp),max(temp)]',i,p, [len(temp), np.std(temp, ddof=polynom_order_intertrace_i+1),min(temp),max(temp)]
        polyfitparams.append(p)
        # Replace the values, which are not in good_values_all with the fit
        polyfit_pl[~good_values_all,i] = np.polyval(p, polyfit_o[~good_values_all])
    return polyfit_pl, polyfitparams

def fit_poly_between_orders_grating_equation(polyfit_pl, polyfit_o, polynom_order_intertrace, sigma):       # only used by correlate_px_wavelength, which is only used by ident_arc.py
    polyfitparams = []
    available_orders = polyfit_pl.shape[1]
    good_values_all = np.repeat([True], polyfit_pl.shape[0])
    for i in range(available_orders):
        good_values, pfit = sigma_clip(polyfit_o, polyfit_pl[:,i], polynom_order_intertrace, sigma, sigma, repeats=20)  #orders, sigma low, sigma high
        good_values_all = (good_values_all & good_values)
    for i in range(available_orders):
        p = np.polyfit(polyfit_o[good_values_all], polyfit_pl[good_values_all,i], polynom_order_intertrace)
        temp=polyfit_pl[good_values_all,i]-np.polyval(p, polyfit_o[good_values_all])
        #print 'i,p, [len(temp), np.std(temp, ddof=polynom_order_intertrace+1),min(temp),max(temp)]',i,p, [len(temp), np.std(temp, ddof=polynom_order_intertrace+1),min(temp),max(temp)]
        polyfitparams.append(p)
        # Replace the values, which are not in good_values_all with the fit
        polyfit_pl[~good_values_all,i] = np.polyval(p, polyfit_o[~good_values_all])
    return polyfit_pl, polyfitparams
    
def remove_bad_identified_lines(arc_lines_wavelength, polyfit_o, polyfit_p, cen_px, sigma_diff_fit):       # only used by correlate_px_wavelength, which is only used by ident_arc.py
    ""
    Sigmaclips the identified reference lines
    :param arc_lines_wavelength
    ""
    arc_lines_wavelength[:, 3] = 0
    for i,order in enumerate(polyfit_o):
        orderpos = (arc_lines_wavelength[:,0] == order)
        xarr = arc_lines_wavelength[orderpos, 1] - cen_px[order]   #to get the wavelength in the center of the order
        yarr = arc_lines_wavelength[orderpos, 2]        # wavelength in Angstrom
        diff_fit = yarr-np.polyval(polyfit_p[i], xarr)
        arc_lines_wavelength[orderpos, 3] = diff_fit
    std_diff_fit = np.std(arc_lines_wavelength[(arc_lines_wavelength[:, 3]!=0), 3], ddof=1) * sigma_diff_fit            # Is ddof=1 right?
    for i in range(arc_lines_wavelength.shape[0])[::-1]:
        if abs(arc_lines_wavelength[i,3]) > std_diff_fit:
            #print 'deleted', arc_lines_wavelength[i], arc_lines_wavelength.shape
            arc_lines_wavelength = np.delete(arc_lines_wavelength, i, axis=0)
    return arc_lines_wavelength

def correlate_px_wavelength(im, orig_arc_lines_wavelength, reference_catalog, reference_names, arc_lines_px, polynom_order_trace, polynom_order_intertrace, max_diff_px, sigma_diff_fit):       # only used by ident_arc.py
    iter_break = 10
    cen_px_curv_range = np.arange(-5,-2.9,0.2)
    cen_px_offset_range = range(-5000,5001,200)
    cen_px_curv_range, cen_px_offset_range = [0], [0]
    ims = im.shape
    orders = np.arange(ims[0])
    result=[]
    for cen_px_curv in cen_px_curv_range:
        for cen_px_offset in cen_px_offset_range:
            cen_px = int(ims[1]/2) + cen_px_curv*orders + cen_px_offset
            arc_lines_wavelength = copy.copy(orig_arc_lines_wavelength)
            old_arc_lines_wavelength = arc_lines_wavelength
            iteration = 0
            while True and iteration < iter_break:
                iteration += 1
                # Fit a polynomial to each order with enough data
                polyfit_o, polyfit_p = [], []
                arc_lines_wavelength[:, 3] = 0
                for order in range(int(np.max(arc_lines_wavelength[:,0]))+1):    #only until the highest order with data
                    orderpos = (arc_lines_wavelength[:,0] == order)
                    xarr = arc_lines_wavelength[orderpos, 1]-cen_px[order]   #to get the wavelength in the center of the order
                    yarr = arc_lines_wavelength[orderpos, 2]        # wavelength in Angstrom
                    if xarr.shape[0] >= polynom_order_trace+2:
                        p = np.polyfit(xarr, yarr, polynom_order_trace)      #lambda from px
                        diff_fit = yarr-np.polyval(p, xarr)
                        #print 'order,p, [len(diff_fit), np.std(diff_fit, ddof=polynom_order_trace+1),min(diff_fit),max(diff_fit)]', order,p, [len(diff_fit), np.std(diff_fit, ddof=polynom_order_trace+1),min(diff_fit),max(diff_fit)]
                        #print xarr, yarr, np.polyval(p, xarr)
                        polyfit_o.append(order)
                        polyfit_p.append(p)
                        arc_lines_wavelength[orderpos, 3] = diff_fit 
                    #else:
                    #    print 'not enough lines, order,xarr.shape', order,xarr.shape
                if len(polyfit_o) < 3:
                    logger('Error: only {2} orders contain enough lines, please provide at least {0} lines in at least {1} orders'.format(polynom_order*2, 4, len(polyfit_o)))
                polyfit_ol = np.array(polyfit_o)
                polyfit_o = np.array(polyfit_o)
                polyfit_pl = np.array(polyfit_p)
                polyfit_p = np.array(polyfit_p)
                #print 'polyfit_o, polyfit_ol, polyfit_p', polyfit_o, polyfit_ol, polyfit_p
                
                # How do the parameters of the fit to the traces change over the detector
                polyfit_pl, polyfitparams = fit_poly_between_orders(polyfit_pl, polyfit_o, polynom_order_intertrace*2,3)
                
                # Add the solution to neightboring orders
                for order in range(arc_lines_px[-1][0]+1):
                    if order in polyfit_o:
                        continue
                    if (order+1 in polyfit_o and order-1 in polyfit_o) or (order+1 in polyfit_o and order+2 in polyfit_o) or (order-1 in polyfit_o and order-2 in polyfit_o):
                        new_polyfit = []
                        for i in range(polyfit_pl.shape[1]):
                            new_polyfit.append(np.polyval(polyfitparams[i],order))      #what is the parameter for this order
                        polyfit_pl = np.insert(polyfit_pl, 0, new_polyfit, axis=0)
                        polyfit_ol = np.insert(polyfit_ol, 0, order, axis=0)
                        #print 'order, new_polyfit, polyfit_ol, polyfit_pl.shape', order, new_polyfit, polyfit_ol, polyfit_pl.shape
                # Delete lines which are far off from fit
                arc_lines_wavelength = remove_bad_identified_lines(arc_lines_wavelength, polyfit_o, polyfit_p, cen_px, sigma_diff_fit)
                
                # Correlate more lines from the spectrum (arc_lines_px) with the catalog (reference_catalog)
                #arc_lines_wavelength = []
                #print 'polyfit_ol, polyfit_pl', polyfit_ol, polyfit_pl
                for arc_line_px in arc_lines_px:                # for each order with lines in the spectrum
                    order = arc_line_px[0]
                    if order in polyfit_ol:                      # if a solution for the order exists
                        order_pos = (polyfit_ol == order)
                        for arc_line in arc_line_px[1]:         # for each line in the order
                            #print 'polyfit_ol, polyfit_pl, arc_line[0], order_pos, polyfit_pl[order_pos,:][0]', polyfit_ol, polyfit_pl, arc_line[0], order_pos, polyfit_pl[order_pos,:][0]
                            arc_line_wave = np.polyval(polyfit_pl[order_pos,:][0], arc_line[0]-cen_px[order])
                            #print 'arc_line_wave.shape, arc_line_wave', arc_line_wave.shape, arc_line_wave
                            diff_catalog = abs(reference_catalog[:,0] - arc_line_wave)
                            if min(diff_catalog) <= max_diff_px*polyfit_pl[order_pos,-2]:   #calculate the approximate wavelength difference for max_diff_px
                                catalog_line = reference_catalog[np.argmin(diff_catalog),:]
                                #print arc_line,arc_line_wave,catalog_line, arc_lines_wavelength.shape,
                                exist_in_arc_lines = abs(arc_lines_wavelength[:,0] - order)*1e7 + abs(arc_lines_wavelength[:,1] - arc_line[0])
                                #print exist_in_arc_lines.shape, min(exist_in_arc_lines), exist_in_arc_lines[0]
                                if min(exist_in_arc_lines) < 5:
                                    #print 'replace', arc_lines_wavelength[np.argmin(exist_in_arc_lines)], [order, arc_line[0], catalog_line[0]]
                                    arc_lines_wavelength[np.argmin(exist_in_arc_lines)] = [order, arc_line[0], catalog_line[0], 0]
                                else:
                                    #print 'insert', [order, arc_line[0], catalog_line[0]]
                                    arc_lines_wavelength = np.insert(arc_lines_wavelength, -1, [order, arc_line[0], catalog_line[0], 0], axis=0)
                            #else:
                                #print 'no line, min(diff_catalog), order, arc_line', min(diff_catalog), order, arc_line
                arc_lines_wavelength = arc_lines_wavelength[(arc_lines_wavelength[:,0]*1e9+arc_lines_wavelength[:,1]).argsort()]    #sort first by column 0 and then by column 1
                #print arc_lines_wavelength.shape,old_arc_lines_wavelength.shape
                if arc_lines_wavelength.shape[0] == old_arc_lines_wavelength.shape[0]:
                    #print arc_lines_wavelength[arc_lines_wavelength != old_arc_lines_wavelength], old_arc_lines_wavelength[arc_lines_wavelength != old_arc_lines_wavelength]
                    if (arc_lines_wavelength[:,0:3] == old_arc_lines_wavelength[:,0:3]).all():
                        break
                old_arc_lines_wavelength = arc_lines_wavelength
                #print 'arc_lines_wavelength.shape, arc_lines_wavelength', arc_lines_wavelength.shape, arc_lines_wavelength
        
            # Get the real order numbers using the grating equation: wavelength propotional to 1/(real_order) -> y(order)=(order_offset+order)*central_wavelength should have the smallest slope
            # Tested that the slope is independend of the central pixel (e.g ims[1]/2)
            order_offset = find_order_offset(polyfit_o, polyfit_p[:,-1])
            order_grating_equation = 1      #should be one
            good_values, pfit = sigma_clip(polyfit_p[:,-1], 1.0/(order_offset+polyfit_o), order_grating_equation, 3, 3, repeats=20)  #orders, sigma low, sigma high
            p_real_cent = np.polyfit(1.0/(order_offset+polyfit_o[good_values]), polyfit_p[good_values,-1], order_grating_equation)
            diff_real_cent=polyfit_p[good_values,-1]-np.polyval(p_real_cent, 1.0/(order_offset+polyfit_o[good_values]))
            #print 'p_real_cent, [len(diff_real_cent), np.std(diff_real_cent, ddof=order_grating_equation+1),min(diff_real_cent),max(diff_real_cent)]',p_real_cent, [len(diff_real_cent), np.std(diff_real_cent, ddof=order_grating_equation+1),min(diff_real_cent),max(diff_real_cent)]
            text = ''
            if order_grating_equation > 1:
                text = '(including distortion, order = {0}) '.format(order_grating_equation)
            logger('Info: the offset to real orders is {0}. With this offset, the standard deviation between the central wavelengths and the grating equation {2}is {1} Angstrom'.format(order_offset, round(np.std(diff_real_cent, ddof=order_grating_equation+1),4),text) )
            # -> even with different cen_px_curv cen_px_offset the std doesn't decrease. Setup with much smaller p_real_cent[-1] exist (<0.1), e.g. -3.8,3000, -3.6,-2600
            
            # See the results
            #correlate_px_wave_result_UI(im, arc_lines_wavelength, reference_catalog, arc_lines_px, reference_names, polyfit_pl, polyfit_ol, cen_px, adjust=[0.07,0.93,0.94,0.06, 1.0,1.01])
            
            # Create the solution over all orders
            cen_wavelength = np.polyval(p_real_cent, 1.0/(order_offset+orders))
            #print cen_wavelength
            # Fit a polynomial to each order with enough data
            polyfit_o, polyfit_p = [], []
            for order in orders:    #only until the highest order with data
                orderpos = (arc_lines_wavelength[:,0] == order)
                xarr = arc_lines_wavelength[orderpos, 1]-cen_px[order]   #to get the central wavelength at the center pixel
                yarr = arc_lines_wavelength[orderpos, 2]        # wavelength in nm
                if xarr.shape[0] >= polynom_order_trace+2:
                    # add the same value a few times in order to be sure the central wavelength from the grating equation is used
                    #xarr = np.insert(xarr, np.zeros(10000).astype(int), 0, axis=0)          # np.zeros(10000) if cen_wavelength can be trusted
                    #yarr = np.insert(yarr, np.zeros(10000).astype(int), cen_wavelength[order], axis=0)
                    p = np.polyfit(xarr, yarr, polynom_order_trace)      #lambda from px
                    diff_fit = yarr-np.polyval(p, xarr)
                    #print order,p, [len(diff_fit), np.std(diff_fit, ddof=polynom_order_trace+1),min(diff_fit),max(diff_fit)]
                    polyfit_o.append(order)
                    polyfit_p.append(p)
            polyfit_o = np.array(polyfit_o)
            polyfit_p = np.array(polyfit_p)
            # How do the parameters of the fit to the traces change over the detector
            polyfit_p, polyfitparams = fit_poly_between_orders_grating_equation(polyfit_p, 1.0/(polyfit_o+order_offset), polynom_order_intertrace, 2.5)
            # Add the solution to missing orders
            for order in orders:
                if order in polyfit_o:
                    continue
                new_polyfit = []
                for i in range(polyfit_p.shape[1]):
                    new_polyfit.append(np.polyval(polyfitparams[i],1.0/(order+order_offset)))      #what is the parameter for this order
                polyfit_p = np.insert(polyfit_p, order, new_polyfit, axis=0)
            # Get some statistics:
            arc_lines_wavelength = remove_bad_identified_lines(arc_lines_wavelength, orders, polyfit_p, cen_px, 5)
            good_arc_lines = arc_lines_wavelength[(arc_lines_wavelength[:, 3]!=0), :]
            #print good_arc_lines.shape,arc_lines_wavelength.shape, min(good_arc_lines[:,3]), max(good_arc_lines[:,3]), np.mean(good_arc_lines[:,3]), np.median(good_arc_lines[:,3]), good_arc_lines[abs(good_arc_lines[:,3])>0.2,:]
            std_diff_fit = np.std(good_arc_lines[:,3], ddof=1)
            logger('Info: used {0} lines. The standard deviation of the fit is {1} Angstrom'.format(good_arc_lines.shape[0], round(std_diff_fit,4)) )
            wavelength_solution = np.insert(polyfit_p, 0, order_offset+orders, axis=1)
            wavelength_solution = np.insert(wavelength_solution, 1, cen_px, axis=1)
            #print wavelength_solution[:,-1] - cen_wavelength
            #wavelength_solution[:,-1] = cen_wavelength     # only if cen_wavelength can be trusted
            temp=percentile_list(good_arc_lines[:,3],0.1)
            result.append([cen_px_curv,cen_px_offset,order_offset,p_real_cent, round(np.std(diff_real_cent, ddof=1),2),good_arc_lines.shape[0], round(std_diff_fit,4),len(temp),round(np.std(temp, ddof=1),4) ])
    if len(result)>1:
        for line in result:
            print line
    # See the results
    correlate_px_wave_result_UI(im, arc_lines_wavelength, reference_catalog, arc_lines_px, reference_names, polyfit_p, orders, cen_px, adjust=[0.07,0.93,0.94,0.06, 1.0,1.01])
    
    return arc_lines_wavelength, wavelength_solution
"""

def create_pseudo_wavelength_solution(number_orders):
    """
    Creates a pseudo wavelength solution with one Angstrom per pixel
    :param number_orders: interger, number of orders for which the solution should be created
    
    """
    wavelength_solution, wavelength_solution_arclines = [], []
    for order in range(number_orders):
        wavelength_solution.append([order, 0., 1., 0. ])
        wavelength_solution_arclines.append([0])     #The wavelength of the reference lines
    return np.array(wavelength_solution), np.array(wavelength_solution_arclines)

def adjust_wavelength_solution(params, spectrum, arc_lines_px, wavelength_solution_ori, wavelength_solution_arclines_ori, reference_catalog_full, reference_names, xlows, xhighs, show_res=False):
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
    
    specs = spectrum.shape
    orders = np.arange(specs[0])
    if len(arc_lines_px) <= 10:
        logger('Warn: no arc lines available -> creating a pseudo solution (1 step per px)')
        wavelength_solution, wavelength_solution_arclines = create_pseudo_wavelength_solution(specs[0])
        return wavelength_solution, wavelength_solution_arclines
    reference_catalog = reference_catalog_full[ reference_catalog_full[:,1] >= np.percentile(reference_catalog_full[:,1], 50) ,:]       # only use the brightest lines at the beginning
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
    pxdiffords = range(min(params['px_offset_order'][0:2]),max(params['px_offset_order'][0:2])+1,params['px_offset_order'][2])
    resolution_offset = params['resolution_offset_pct']/100.
    res_steps = 11
    resdiffs = np.linspace(1-resolution_offset, 1+resolution_offset, res_steps)
    best_matching = []
    #logger('Step: Compare old wavelength solution with current arc lines')
    # Rescaling the original wavelength solution, calculating the wavelength for the identied lines with this solution
    for resdiff in tqdm(resdiffs, desc='Compare old wavelength solution with current arc lines'):
        wavelength_solution_1 = copy.deepcopy(wavelength_solution_ori)
        wavelength_solution_1[:,-2] *= resdiff                                                      # Scale the resolution
        for pxdiff in pxdiffs:
            for pxdifford in pxdiffords:
                wavelength_solution_2 = copy.deepcopy(wavelength_solution_1)
                wavelength_solution_2[:,1] += pxdiff + pxdifford * np.arange(len(wavelength_solution_2))        # that is equivalent to changing the px position in arc_lines_order_px
                for orderdiff in orderdiffs:
                    matching_lines = np.empty(( len(arc_lines_px)*10, 4 ))
                    index_matching = 0
                    data_available = [0, 0, 0]
                    # Caluculate the wavelength of all arclines with the current settings
                    for order_arcline in np.arange(min(arc_lines_px[:,0]), max(arc_lines_px[:,0])+1, dtype=int):        # The new orders
                        if order_arcline+orderdiff < 0 or order_arcline+orderdiff > len(wavelength_solution_ori)-1 or order_arcline in ignoreorders:    # ignore orders that are not covered by orig solution
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
                        continue
                    #matching_lines = np.array(matching_lines)
                    matching_lines = matching_lines[:index_matching,:]                                  # get rid of unfilled array
                    #best_matching: orderdiff, pxdiff, pxdifford, resdiff, number of identified reference lines, 
                    #               stdev of the wavelength difference between current solution and reference wavelength, sum of the flux of the reference lines, 
                    #               number of covered pixels, number of covered reference lines, sum of covered flux
                    best_matching.append([orderdiff, pxdiff, pxdifford, resdiff, len(matching_lines), np.std(matching_lines[:,2], ddof=1), np.sum(matching_lines[:,3])] + data_available )
                    #print best_matching[-1]
    if best_matching == []:
        if len(wavelength_solution_ori) == specs[0]:
            logger('Warn: No matching configuration of the lines in the arc with the old wavelength solution found. Therefore the old solution will be used')
            return np.array(wavelength_solution_ori), wavelength_solution_arclines_ori
        else:
            logger('Warn: No matching configuration of the lines in the arc with the old wavelength solution found. Additionally the number of orders in the original solution and this setting do not match -> creating a pseudo solution (1 step per px)')
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
    logger('Info: To match the most lines in the arc with the old wavelength solution, a shift of {0} orders, a multiplier to the resolution of {1}, a shift of {2} px, and a shift of {3} px per order needs to be applied. {4} lines were identified. The deviation is {5} Angstrom.'.format(orderdiff, resdiff, pxdiff, pxdifford, int(best_matching[0,4]), round(best_matching[0,5],4) ))
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
    for step in tqdm(range(steps), desc='Finding the new wavelength solution'):
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
                goodvalues, p = sigma_clip(np.arange(arc_lines_wavelength.shape[0]), arc_lines_wavelength[:,3], 0, sigma, sigma, repeats=20)    # poly_order=0 to sigma clip around the average
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
                plot_wavelength_solution_spectrum(spectrum, spectrum, params['logging_arc_line_identification_spectrum'], wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names)
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
    x = arc_lines_wavelength[:,1]-cen_pxs
    y = 1.0/(arc_lines_wavelength[:,0]+order_offset)
    
    weight = arc_lines_wavelength[:,4] * arc_lines_wavelength[:,5] * arc_lines_wavelength[:,6] * abs(arc_lines_wavelength[:,1])**(polynom_order_trace-1)
    median = []
    for order in np.unique(arc_lines_wavelength[:,0]):
        median.append( np.max( weight[arc_lines_wavelength[:,0] == order] ))
    divisor = np.median(median)
    if divisor == 0:
        divisor = 1.
    median = np.array(median)/(divisor+0.0)
    median[median == 0] = 1.
    for i, order in enumerate(np.unique(arc_lines_wavelength[:,0])):
        weight[arc_lines_wavelength[:,0] == order] /= (median[i]+0.0)
    weight = []         # No weight fits the data better
    poly2d_params = polynomial_fit_2d_norm(x, y, arc_lines_wavelength[:,2], polynom_order_trace, polynom_order_intertrace, w=weight)
    
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
    # Remove the most scattered lines (highest 1% and lowest 1%
    todel = np.argsort(arc_lines_wavelength[:,3])
    todel = np.delete(todel, np.arange(int(len(todel)/100), int(len(todel)*99/100)), axis=0)
    arc_lines_wavelength = np.delete(arc_lines_wavelength, todel, axis=0)
    #print np.histogram(arc_lines_wavelength[:,3], 20), arc_lines_wavelength.shape
    avg_diff_fit = np.mean(np.abs(arc_lines_wavelength[:,3]))           # Diff between catalog wavelength and fitted wavelength
    std_diff_fit = np.std(arc_lines_wavelength[:,3], ddof=(polynom_order_trace+polynom_order_intertrace+1))             # real standard deviation, asume values +-1 to convice yourself
    # Resolution from the precission of the fit
    std_R_fit = 1.0/np.std(arc_lines_wavelength[:,3]/(arc_lines_wavelength[:,2]+0.0), ddof=(polynom_order_trace+polynom_order_intertrace+1))    # lambda/d_lambda
    # Resolution using the Gaussian width of the arc lines
    R_gauss     = arc_lines_wavelength[:,2]/(1.*arc_lines_wavelength[:,-1]*arc_lines_wavelength[:,-3])    # lambda/d_lambda ; -1 is in px, -3 is resolution
    good_values, p = sigma_clip(R_gauss, R_gauss, 0, 3, 3)   #x, y, order of polynom (if 0 than average -> x doesn't matter), sigma, sigma
    R_gauss     = R_gauss[good_values]
    avg_R_gauss = np.mean(R_gauss)
    std_R_gauss = np.std(R_gauss, ddof=1)
    # 2px Resolution (using only the identified arc lines
    R_2px     = arc_lines_wavelength[:,2]/(2.*arc_lines_wavelength[:,-3])    # lambda/d_lambda
    avg_R_2px = np.mean(R_2px)
    std_R_2px = np.std(R_2px, ddof=1) 
    text = 'Info: used {0} lines. The standard deviation of the residuals between the lines and the fit is {1} Angstrom '+\
                  '(the average of the abs of the residuals is {7} Angstrom). This converts into a resolution R = {2}. '+\
                  'The Gaussian width of the emission lines results in an R = {3} +- {4}. The 2-pixel resolution (around the identified lines) is R = {5} +- {6}.'
    logger(text.format(arc_lines_wavelength.shape[0], round(std_diff_fit,4), int(std_R_fit), int(avg_R_gauss), int(std_R_gauss), int(avg_R_2px), int(std_R_2px), round(avg_diff_fit,4) ))
    p_cen_px = np.round(p_cen_px,3)
    text = 'Info: A 2d polynom fit with {0} orders in dispersion direction (along the traces) and {1} orders in cross-dispersion direction was used. '+\
                  'With this solution, the offset between aperture and real orders is {2}. To fulfil the grating equation the central pixel of the individual orders needs to be {5} + {6}*order + {7}*order**2.'+\
                  'With this values the standard deviation of the residuals between the central wavelengths and the grating equation is {3} Angstrom. Using the original solution gives an offset of {4}.'
    logger(text.format(polynom_order_trace, polynom_order_intertrace, int(order_offset), round(np.std(diff_real_cent, ddof=len(p_real_cent)),4), order_offset_old, p_cen_px[2], p_cen_px[1], p_cen_px[0] ))
    
    # Create a wavelength solution
    wavelength_solution = np.array( [polynom_order_trace, polynom_order_intertrace, np.mean(cen_px), order_offset] + list(poly2d_params) )
    text = 'Info: Wavelenght solution in 2d (for pixel and order at the same time) is [No of orders1, No of orders2, mean central pixel, '+\
                  'offset to real order, parameters of the 2d polynomial(px^0*m^0, px^1*m^0, px^0*m^1, px^2*m^0, px^1*m^1, px^0*m^2, ....)]: \n'+\
                  '{0}'
    logger(text.format(wavelength_solution), show=False)
    
    printarrayformat = ['%1.1i', '%3.1f', '%9.4f', '%6.4f']
    logger('order\tpixel\twavelength\tdwavel', show=False, printarrayformat=printarrayformat, printarray=arc_lines_wavelength[:,0:4], logfile=params['logging_identified_arc_lines'])
    
    # Create an image of the gaussian widths of the identified lines
    plot_wavelength_solution_width_emmission_lines(params['logging_em_lines_gauss_width_form'], specs, arc_lines_wavelength)
    
    x,y,l=[],[],[]
    max_number_reflines = 0             # will be needed later in order to save the data correctly into a fits file
    for order in orders:
        x.append(arc_lines_wavelength[(arc_lines_wavelength[:,0]==order),1])    # Pixel
        y.append(arc_lines_wavelength[(arc_lines_wavelength[:,0]==order),3])    # Residuals
        l.append('{0}'.format(order))
        max_number_reflines = max(max_number_reflines, len(arc_lines_wavelength[(arc_lines_wavelength[:,0]==order),2]))
    text = 'Residuals of the identificated arc lines: identified {0} lines for a 2d polynom fit with {1} and {2} orders'.format(arc_lines_wavelength.shape[0], polynom_order_trace, polynom_order_intertrace)
    plot_img_spec.plot_points(x,y,l,[params['logging_arc_line_identification_residuals']],show=False, title=text, x_title='Pixel', y_title='Residuals [Angstrom]', marker=['o','s','*','P','^','v','>','<','x'])

    # Transform the wavelength solution into the old wavelength solution
    wavelength_solution, wavelength_solution_arclines = [], []
    xarr = np.arange(0,specs[1], 0.1)
    for i,order in enumerate(orders):
        yarr = polynomial_value_2d(xarr-cen_px[i], 1.0/(order+order_offset), polynom_order_trace, polynom_order_intertrace, poly2d_params)
        polyfit = np.polyfit(xarr-cen_px[i], yarr, polynom_order_trace)      #lambda from px
        solution = [order_offset+order, cen_px[i] ] + list(polyfit)
        wavelength_solution.append(solution)
        reflines = list(arc_lines_wavelength[(arc_lines_wavelength[:,0]==order),2])     # Wavelength
        for j in range(max_number_reflines - len(reflines)):                            # Fill up with zeros in order to create an array
            reflines.append(0)
        wavelength_solution_arclines.append(reflines)     #The wavelength of the reference lines
        
        #print order,polynomial_value_2d(0, 1.0/(order+order_offset), polynom_order_trace, polynom_order_intertrace, poly2d_params), solution
        #diff_fit = yarr-np.polyval(polyfit, xarr-cen_px[0])        #<1e-10
        #print [len(diff_fit), np.std(diff_fit, ddof=polynom_order_trace+1),min(diff_fit),max(diff_fit)]
    wavelength_solution = np.array(wavelength_solution)
    wavelength_solution_arclines = np.array(wavelength_solution_arclines)
    statistics_arc_reference_lines(arc_lines_wavelength, [0,-2,-1,2], reference_names, wavelength_solution, xlows, xhighs)
        
    if std_diff_fit > 2*max(wavelength_solution[:,-2]) or std_diff_fit < 1E-8:        # if residuals are bigger than 1px or unreasonable small
        plot_wavelength_solution_spectrum(spectrum, spectrum, params['logging_arc_line_identification_spectrum'], wavelength_solution, 
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
        correlate_px_wave_result_UI(spectrum, arc_lines_wavelength, reference_catalog, arc_lines_px, reference_names, wavelength_solution, adjust=[0.07,0.93,0.94,0.06, 1.0,1.01])
    
    return wavelength_solution, wavelength_solution_arclines

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

def read_fit_wavelength_solution(params, filename, im):
    """
    Reads the file with pixel - wavelength corellation and fits a 2d array against it to create a rough wavelength solution
    :params filename: string to a textfile with the following data: order (starting at 0), real order, pixel, wavelength. If one of the information is missing, this line will be skipped.
    """
    arc_lines_wavelength = []
    file = open(filename, 'r')
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
    file.close()
    if len(arc_lines_wavelength) == 0:
        logger(('Error: No useful information was found in {0}. Please make sure the entries in the file contain of the tab-separated '+\
                '(exactly one tab between each column) values: order (starting at 0), real order, pixel, wavelength').format(filename))
    arc_lines_wavelength = np.array(arc_lines_wavelength)
    if 0.0 in arc_lines_wavelength[:,1]:
        logger('Error: The second coloumn (real order) in file {1} was not set correctly as it contains a zero. '+\
               'Please use the (estimated) real order (from grating equation).'.format(filename))
    orig_lines = arc_lines_wavelength.shape[0]
    
    order_offset = arc_lines_wavelength[:,1] - arc_lines_wavelength[:,0]
    orders = np.arange(max(arc_lines_wavelength[:,0])+1, dtype=int)
    if np.std(order_offset, ddof=1) != 0:
        logger('Error: There is an inconsistency between coloumn 1 (order/aperture) and coloumn 2 (real/physical order) in the file {0}. '+\
               'Please check that the difference between these two columns is always the same.'.format(filename))
    order_offset = order_offset[0]
    polynom_order_trace = max(2, max(params['polynom_order_traces']) )
    polynom_order_intertrace = max(1, max(params['polynom_order_intertraces']) )
    
    cen_px = np.repeat([im.shape[0]/2.], len(orders))       # cen_px needs to be float, cen_px defines the zeropoint of np.polyfit
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
    logger( ('Info: The standard deviation of the residual of the fit to the manual wavelength solution is {1} Angstrom (average is {0} Angstrom). '+\
             'Only input data, for which the residuals were less than 1 Angstrom have been used. '+\
             '{2} of {5} lines have been used to calculate the solution for a 2d polynom fit with {3} orders in dispersion direction (along the traces) '+\
             'and {4} orders in cross-dispersion direction. To fulfil the grating equation the central pixel of the individual orders needs to be '+\
             '{6} + {7}*order + {8}*order**2.').format( round(np.mean(arc_line_res),4), 
                            round(np.std(arc_line_res, ddof=polynom_order_trace+polynom_order_intertrace+1),4), arc_lines_wavelength.shape[0], 
                            polynom_order_trace, polynom_order_intertrace, orig_lines, p_cen_px[2], p_cen_px[1], p_cen_px[0] ))
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

def remove_orders(pfits, rm_orders):
    """
    Remove orders by masking them
    :param pfits:
    :param rm_orders: list of integers, contains the orders to be removed
    return mask: array of bool with the same length as original orders
    """
    mask = np.repeat([True], len(pfits))
    if type(rm_orders) is not list:
        return mask
    # Remove orders (if rm_orders is populated)
    for r in range(len(pfits)):
        if r in rm_orders:
            mask[r] = False
    return mask

def run_remove_orders_UI(im1, pfits, xlows, xhighs, userinput=True):
    """
    Removes orders from the pfits array in a GUI
    """
    # Plot fitted orders
    if not userinput:
        return pfits, xlows, xhighs, True

    # convert to numpy arrays
    pfits = np.array(pfits)
    xlows, xhighs = np.array(xlows), np.array(xhighs)

    # set up plot
    fig, frame = plt.subplots(1, 1)

    rm_orders = []
    # get kwargs
    pkwargs = dict(frame=frame, im1=im1, pfits=pfits, xlows=xlows,
                   xhighs=xhighs, rm_orders=rm_orders)

    # define update plot function
    def plot(frame, im1, pfits, xlows, xhighs, rm_orders):
        frame.clear()
        title = ('Removing bad orders  \n(Largest order '
                 'number = {0})'.format(len(pfits)))
        mask = remove_orders(pfits, rm_orders)
        frame = plot_traces_over_image(im1, 'dummy_filename', pfits, xlows, xhighs, mask=mask, frame=frame, return_frame=True)
        frame.set_title(title)

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

    # define widgets
    widgets = dict()
    widgets['rm_orders'] = dict(label='Select orders to remove',
                                comment='Enter all order numbers to remove \n'
                                        'separated by a whitespace or comma \n'
                                        'to undo just delete the entered '
                                        'number',
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=str, start=" ", valid_function=vfunc,
                                width=60)
    widgets['accept'] = dict(label='Accept Orders', kind='ExitButton',
                             position=Tk.BOTTOM)
    widgets['update'] = dict(label='Update', kind='UpdatePlot',
                             position=Tk.BOTTOM)
                             
    wprops = dict(orientation='v', position=Tk.RIGHT)

    gui3 = tkc.TkCanvas(figure=fig, ax=frame, func=plot, kwargs=pkwargs,
                        title='Locating orders', widgets=widgets,
                        widgetprops=wprops)
    gui3.master.mainloop()

    if 'rm_orders' in gui3.data:
        if gui3.data['rm_orders'] is not None:
            if type(gui3.data['rm_orders']) == list:
                rm_orders = gui3.data['rm_orders']

    fmask = remove_orders(pfits, rm_orders)
    plt.close()
    return fmask

def plot_gauss_data_center(datapoints_x, datapoints, label_datapoins, gauss, label_gauss, centroids, label_centroids, filename='', title=''):
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
        gauss_x, gauss_y = [], []
        for i in datarange:
            x = np.linspace(min(datapoints_x[i,:]), max(datapoints_x[i,:]), len(datapoints_x[i,:])*20 )
            y = oneD_gauss(x,gauss[i,:])
            gauss_x.append(x)
            gauss_y.append(y)
        plot_img_spec.plot_points(gauss_x, gauss_y, np.array(label_gauss)[datarange], '', show=False, adjust=adjust, title='', return_frame=True, frame=frame, x_title='Pixel', y_title='Flux', linestyle="-", marker="")
        ymin2, ymax2 = axes.get_ylim()
        xmin, xmax = axes.get_xlim()
        plt.gca().set_prop_cycle(None)          #Reset the color cycle
        ymin, ymax = min(ymin1,ymin2), max(ymax1,ymax2)
        centr_x = np.repeat([centroids[datarange]],2,axis=0).T
        centr_y = np.repeat([[ymin, ymax]], len(datarange), axis=0)
        plot_img_spec.plot_points(centr_x, centr_y, np.array(label_centroids)[datarange], '', show=False, adjust=adjust, title='', return_frame=True, frame=frame, x_title='Pixel', y_title='Flux', linestyle="--", marker="")
        axes.set_ylim(ymin,ymax)
        axes.set_xlim(xmin,xmax)
        frame.set_title(title, fontsize=16)
        if filename == '':
            plt.show()
        else:
            plt.savefig(filename.replace('.png','_{0}-{1}.png'.format('%3.3i'%datarange[0],'%3.3i'%datarange[-1])), bbox_inches='tight')
            plt.close()

def add_specinfo_head(spectra, s_spec, n_spec, im_head):
    """
    Add information from the spectra to the header
    :param spectra:
    :param s_spec: 2d-array of floats, spectrum with the signal
    :param n_spec: 2d-array of floats, spectrum with the noise
    """
    spectra = np.array(spectra)
    signal = np.nansum(spectra, axis=1)         # Sum of flux in each order
    snr = np.array(s_spec)/np.array(n_spec)
    cenpx = spectra.shape[1]/2.
    px_range = list(range(int(round(0.8*cenpx)),int(round(1.2*cenpx))+1))
    snr1 = np.nanmedian( snr[:,px_range], axis=1)                   # median SNR in the centre of each order: +-10% around 
    snr2 = np.sqrt( np.nanmedian( s_spec[:,px_range], axis=1) )     # sqrt of flux
    snr = np.nanmin( np.vstack(( snr1, snr2 )), axis=0 )            # SNR can't be better than sqrt of flux
    snr[np.isnan(snr)] = -1
    im_head['HIERARCH EXO_PIPE fmin'] = (np.nanmin(spectra), 'Minimum Flux per pixel')
    im_head['HIERARCH EXO_PIPE fmax'] = (np.nanmax(spectra), 'Maximum Flux per pixel')
    im_head['HIERARCH EXO_PIPE fsum_all'] = (np.nansum(spectra, axis=None), 'Total flux')
    for order in range(spectra.shape[0]):
        im_head['HIERARCH EXO_PIPE fsum_order{0}'.format('%2.2i'%order)] = ( round(signal[order],2), 'Flux (counts) in order {0}'.format('%2.2i'%order) )
        im_head['HIERARCH EXO_PIPE SN_order{0}'.format('%2.2i'%order)] = ( round(snr[order],1), 'median SNR (+-10% of centre px)' )
        
    return im_head


def normalise_continuum(spec, wavelengths, nc=8, ll=2., lu=4., frac=0.3, semi_window=10, nc_noise=15):      # without modal noise nc=4,ll=1.,lu=5. might work
    """
    Normalises the spectrum using the area with continuum and excluding the area with lines
    Adapted from the ceres pipeline, workes for HARPS data and modal noise corrected (blaze corrected) EXOhSPEC data
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
        old_good_values = sum(sub[4,:])
        indexes_search_range = np.zeros((specs[1], specs[1] ), dtype=bool)            # Array to store the indexes which are in a given wavelength difference for each wavelength
        for i in data[0,sub[0,:]].astype(int):                                      # Only use values, for which the wavelength is available
            indexes_search_range[i,sub[0,:]] = ( np.abs(data[1,sub[0,:]] - data[1,i]) <= 0.9)    #+- in Angstrom
            #print i, sum(indexes_search_range[i,:]), data[1,sub[0,:]] - data[1,i]
        for step in range(1,int(old_good_values/10)):
            #print order, step, len(data[1,sub[4,:]]), len(data[4,sub[4,:]]), sum(np.isnan(data[1,sub[4,:]])), sum(np.isnan(data[4,sub[4,:]])), data[1,sub[4,:]]-x_cen, data[4,sub[4,:]], min(step,nc)
            p = polyfit_adjust_order(data[1,sub[4,:]]-x_cen, data[4,sub[4,:]], min(step,nc) )        # Fit the data, !!! In Ceres Pipeline without x_cen -> polynomial fit isn't as good as here
            data[5,:] = data[4,:] - np.polyval(p,data[1,:]-x_cen)                   # Residuals for the whole array, >0 for emission lines
            IU = np.where(data[5,sub[4,:]] > 0)[0]                                  # Positive residuals, use only the good data
            dev = np.mean(data[5,sub[4,:]][IU])                                     # deviation using the positive residuals, use only the good data -> without absorption lines
            sub[4,:] = False
            sub[4,sub[3,:]] = ( (data[5,sub[3,:]] < lu*dev) & (data[5,sub[3,:]] > -ll*dev) )  # Sigmacliping, but allow already clipped data back in
            for i in data[0, ~sub[4,:]*sub[1,:]].astype(int):                       # Points that will be removed
                #print i
                if np.sum( sub[4,indexes_search_range[i,:]] ) <= 0:                 # Not enough datapoints remain in a given area
                    #print 'changed',i
                    sub[4,i] = True                                                 # Don't create too big gaps as the high order polynomial nc>4 will not behave well there
            if ( sum(sub[4,:]) < frac*ori_len or sum(sub[4,:]) == old_good_values ) and step >= nc:         # stop if too much data is removed or no changes, but not if only linear equation is fitted
                if order == 129:
                    plot_img_spec.plot_points([data[1,:],data[1,:],data[1,:],data[1,sub[4,:]]], [data[2,:], data[3,:], np.polyval(p,data[1,:]-x_cen), data[3,sub[4,:]]], \
                                      ['data ori', 'data medfilt', 'fit {0} orders'.format(len(p)-1), 'remaining'], '', show=True, return_frame=False, \
                                      x_title='wavelength [Angstrom]', y_title='flux [ADU]', title='order {1}, step={0}'.format(step, order))
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
            for dummy in range(100):
                np.warnings.filterwarnings('ignore')
                if sum(~sub[5,:]*sub[1,:]) > 0:                                     # There are values that can be replaced
                    index = data[0,sub[1,:]][np.abs(data[6,sub[1,:]] - np.nanmin(data[6,~sub[5,:]]) ) < 1E-4].astype(int)   # Find the position where the fit has the smallest value and the residual there isn't used for the fit
                else:
                    index = data[0,~sub[6,:]][np.abs(data[6,~sub[6,:]] - np.nanmin(data[6,~sub[6,:]]) ) < 1E-4].astype(int) # Find the smalles value of the fit, ignore already changed values
                np.warnings.resetwarnings()
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
    #ccoefs = np.array(ccoefs)
    contspec, sn_cont = [], []
    showorder = 129
    for order in range(specs[0]):
        fit_flux = np.polyval(ccoefs[order],wavelengths[order,:]-cen_waves[order])
        fit_sn = np.polyval(sn_coefs[order],wavelengths[order,:]-cen_waves[order])
        for i in range(len(fit_flux)-1):
            if fit_flux[i]*fit_flux[i+1] < 0:                           # change of sign
                fit_flux[max(0,i-semi_window):min(len(fit_flux),i+semi_window+1)] = np.nan     # ignore this data as artificial peaks are introduced due to the division close to 0 in the flat spectrum
        np.warnings.filterwarnings('ignore')
        bad_value = ( (fit_sn <= 0) | (fit_flux <= 0) )                                  # Errors smaller than 0, or Flux smaller than 0 (will make absorption lines into emission lines)
        np.warnings.resetwarnings()
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

"""def get_cont_ceres(W,F,nc=3,ll=1.,lu = 5.,frac=0.05,window=21):     # Based on ceres pipeline (copied, added comments, maybe changed code), not used anymore
    blns = [[6755,6769],[6530,6600],[4840,4880],[4320,4360],[4085,4120],[3950,3990],[3880,3910],[3825,3850]]        #These wavelength ranges will be ignored
    flx = F[0]
    wav = W[0]
    I = np.where( (flx!=0) & (~np.isnan(flx)) )[0]
    wav,flx = wav[I],flx[I]      #only order 0, as preparation for next step

    for i in range(F.shape[0]-1):           # For each-1 order
        f = F[i+1]
        w = W[i+1]
        I = np.where( (f!=0) & (~np.isnan(f)) )[0]
        w,f = w[I],f[I]                     # wav, flux for order 1 and larger
        I = np.where(w>=wav[0])[0]          # only wavelengths bigger than the smallest in the previous order   (all in the wrong setup, overlapping in the harps)
        
        if len(I)>5:                        # check for the overlap with previous order
            J = np.where(wav<w[-1])[0]      # in the previous order get the wavelengths which are smaller than the biggest one in this order    (all in the wrong setup, overlapping in the harps)
            tck = scipy.interpolate.splrep(w[I],f[I],k=3)       # Get the B-sline of the overlapping flux in the current order
            nf = scipy.interpolate.splev(wav[J],tck)            # Apply the B-spline to the previous order
            flx[J] = .5*(flx[J]+nf)         # Adjust the overlapping flux in the previous order to the average of the flux and the flux of the bspline
        I = np.where(w<wav[0])[0]           # only wavelengths smaller than the smallest in the previous order  (none in the wrong setup, non-overlapping the the harps
        wav = np.hstack((w[I],wav))         # Add the missing wavelengths
        flx = np.hstack((f[I],flx))         # ... and the according flux
    # The steps until here are done in order to allow smooth overlap between orders (I think)

    wavo,flxo = wav.copy(),flx.copy()
    for lns in blns:                                    # Get rid off flux in problematic wavelength ranges
        I = np.where((wav>lns[0])&(wav<lns[1]))[0]
        wav,flx = np.delete(wav,I),np.delete(flx,I)
    for i in range(F.shape[0]):             # for each order get some section of the wavelengths (+-2 orders if possible), full wavelength range
        if i == 0:
            wi = W[i+2,0]                   # (harps: smallest wavelength)
            wf = W[i,-1]                    # (harps: biggest wavelength)
        elif i == 1:
            wi = W[i+2,0]
            wf = W[i-1,-1]                  # changed line to remove bug
        elif i == W.shape[0] - 1:
            wi = W[i,0]
            wf = W[i-2,-1]
        elif i == W.shape[0] - 2:
            wi = W[i+1,0]
            wf = W[i-2,-1]
        else:
            wi = W[i+2,0]
            wf = W[i-2,-1]

        IM = np.where((wav>wi) & (wav<wf))[0]
        tw,tf = wav[IM],scipy.signal.medfilt(flx[IM],window)    # wavelength and median filtered flux (with 21px window) for the flux in +-2 orders     # changed line to remove bug
        #JJJ = np.where((wav>W[i,0])&(wav<W[i,-1]))[0]           # only wavelength in this order     # not used
        
        # Fit the polynomial only through the selected set of data +-2 orders
        ori = len(tw)
        coef = np.polyfit(tw,tf,nc)     # it says "RankWarning: Polyfit may be poorly conditioned", even though that more than 6000 entries are in tw, tf. Does happen for nc=5 or nc=10, but not nc=3
        #plot_img_spec.plot_points([tw,tw], [tf, np.polyval(coef,tw)], ['data', 'fit'], '', show=True, return_frame=False, x_title='wavelength [Angstrom]', y_title='flux [ADU]')
        res = tf - np.polyval(coef,tw)
        IU = np.where(res>0)[0]
        dev = np.mean(res[IU])
        I = np.where((res<lu*dev)&(res>-ll*dev))[0]
        cond = True

        while cond:
            tw,tf = tw[I],tf[I]
            coef = np.polyfit(tw,tf,nc)
            #plot_img_spec.plot_points([tw,tw], [tf, np.polyval(coef,tw)], ['data', 'fit'], '', show=True, return_frame=False, x_title='wavelength [Angstrom]', y_title='flux [ADU]')
            res = tf - np.polyval(coef,tw)
            IU = np.where(res>0)[0]
            dev = np.mean(res[IU])
            #print dev
            #plot_img_spec.plot_points([tw], [res], ['res'], '', show=True, return_frame=False, x_title='wavelength [Angstrom]', y_title='flux [ADU]')
            I = np.where((res<lu*dev)&(res>-ll*dev))[0]
            J1 = np.where(res>=lu*dev)[0]
            J2 = np.where(res<=-ll*dev)[0]
        
            if (len(J1)==0 and len(J2)==0) or len(tw)<frac*ori:
                cond = False
        #plot_img_spec.plot_points([tw,tw], [tf, np.polyval(coef,tw)], ['data', 'fit'], '', show=True, return_frame=False, x_title='wavelength [Angstrom]', y_title='flux [ADU]')
        if i == 0:
            coefs = coef
        else:
            coefs = np.vstack((coefs,coef))

    return coefs"""

def linearise_wavelength_spec(params, wavelength_solution, spectra, method='sum', weight=[]):    
    """
    :param spectra: list of 2d array of floats, contains the spectrum to be linearised
    """
    specs = spectra.shape
    dwave = params['wavelength_scale_resolution']
    px_range = np.arange(specs[1])
    data = []
    for order in range(specs[0]):
        wave_range = np.polyval(wavelength_solution[order,2:], px_range-wavelength_solution[order,1])
        wave_diff = np.nanmax(wave_range[1:] - wave_range[:-1])        # Difference in wavelength between 2 data points
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
        lspec = np.empty(lwave_range.shape)
        lspec.fill(np.nan)
        lweight = np.ones(lwave_range.shape)
        for i,wave in enumerate(lwave_range):
            diff = np.abs(wave_range - wave)                # Wavelength difference to the original wavelengths
            indexes = np.where( diff < wave_diff )[0]
            #print i, wave, spectra[order,indexes].shape, diff[indexes].shape
            lspec[i] = np.average( spectra[order,indexes], weights=1./diff[indexes] )
            lweight[i] = np.average( weight[order,indexes], weights=1./diff[indexes] )
        #print np.nanmean(spectra[order]), np.nanmean(lspec)
        #plot_img_spec.plot_points([wave_range, lwave_range], [spectra[order], lspec], ['orig', 'lin'], '', show=True, return_frame=False, x_title='wavelength [Angstrom]', y_title='flux [ADU]')
        
        
        if len(data) == 0:
            data = np.vstack([lwave_range,lspec, lweight])
        else:
            data = np.hstack([data, np.vstack([lwave_range,lspec, lweight]) ])
    data = data.T                # wavelengths in first column, flux in second
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
            if len(data[indexes,1][notnan]) > 0:
                #print data[indexes,:], data[indexes,1][notnan], data[indexes,2][notnan], sum(data[indexes,2][notnan])
                data[indexes[0],1] = np.average(data[indexes,1][notnan], weights=data[indexes,2][notnan])
        else:
            if not nolog:
                logger('Warn: method {0} is not known, using sum'.format(method))
                nolog = True
            data[indexes[0],1] = np.nansum(data[indexes,1])
    data = np.delete(data, todel, axis=0)
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

def mjd_fromheader(params, head):
    """
    :return mjd: modified Julian date from header, in UT of the mid of the exposure
    :return mjd0: 2400000.5
    """
    secinday = 24*3600.0
    
    obsdate, obsdate_float, exposure_time, obsdate_begin, exposure_fraction = get_obsdate(params, head)               # obsdate, obsdate_float, in UTC, mid of the exposure
    # Mid exppsure
    mjd0,mjd = iau_cal2jd(obsdate.year, obsdate.month, obsdate.day)
    ut       = obsdate.hour*3600. + obsdate.minute*60. + obsdate.second
    mjd += ut / secinday                    # Mid-exposure MJD
    # Start exposure
    mjd0,mjd_begin  = iau_cal2jd(obsdate_begin.year, obsdate_begin.month, obsdate_begin.day)
    ut              = obsdate_begin.hour*3600. + obsdate_begin.minute*60. + obsdate_begin.second
    mjd_begin       += ut / secinday        # Start-exposure MJD
    
    return mjd, mjd0, mjd_begin

def iau_cal2jd(IY,IM,ID):               # from CERES, modified
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
    return DJM0, DJM

def getcoords_from_file(obnames,mjd,filen='coords.txt'):               # from CERES, heavily modified
    """
    1-  name of the target as specified in the image header.
    2-  right ascension of the target (J2000) with format hh:mm:ss or as float in degrees.
    3-  declination of the target (J2000) with format dd:mm:ss or as float in degrees.
    4-  proper motion in RA [mas/yr].
    5-  proper motion in DEC [mas/yr].
    6-  epoch (J2000)
    7-  integer (0 or 1). If 0, the code uses the cordinates given in the image header.
        If 1, the code uses the coordinates given in this file.
    old7-  mask that will be used to compute the CCF. Allowed entries are G2, K5 and M2.
    old8-  velocity width in km/s that is used to broaden the lines of the binary mask.
        It should be similar to the standard deviation of the Gaussian that is fitted to the CCF.
    """
    RA0, DEC0, PMRA, PMDEC, epoch = 0., 0., 0., 0., 2000.
    
    if not os.path.isfile(filen):
        logger('Warn: Reference coordinates files {0} does not exist.'.format(filen))
        return RA0, DEC0, epoch, PMRA, PMDEC, obnames
    # !!! Use read files with split already available
    lines_txt = read_text_file(filen, no_empty_lines=True)
    lines = convert_readfile(lines_txt, [str, str, str, float, float, float, float], delimiter=',', replaces=[['\t',',']])
    if len(lines) < len(lines_txt):
        logger('Warn: {1} line(s) could not be read in the reference coordinates file: {0}. Please check that columns 4 to 7 (starting counting with 1) are numbers'.format(filen, len(lines_txt)-len(lines) ))
    found = False
    for cos in lines:
        if abs(cos[6]) < 1E-5:         # disabled
            continue
        for obname in obnames:
            if cos[0].lower() == obname.lower() or cos[0].lower().replace('_','') == obname.lower() or cos[0].lower().replace('-','') == obname.lower() or cos[0].lower().replace(' ','') == obname.lower():
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
                
                """ # Steps to adjust the coordinates to the correct position, not necessary anymore because not using CERES routines to get BCVel and the barycorrpy package takes care of it
                mjdepoch = 2451545.0 - 2400000.5 + (float(cos[5]) - 2000.)
    
                RA  = RA0 + (PMRA/1000./3600.)*(mjd-mjdepoch)/365.
                DEC = DEC0 + (PMDEC/1000./3600.)*(mjd-mjdepoch)/365."""
                obnames = [obname]
                found = True
                break
                
    if not found:
        logger('Warn: Object was not found in the reference coordinates file {0}.'.format(filen))
    return RA0, DEC0, epoch, PMRA, PMDEC, obnames

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

def obspos(longitude,obsradius,R0):               # from CERES, not needed anymore
    """
        Set the observatory position respect to geocenter in the coordinates(x,y,z)required by jplepem, 
        x to equator/greenwich intersection, 
        y 90 degrees east, 
        z positive to north
        """
    obpos = []    
    x = obsradius*np.cos( (np.pi / 180.0) * longitude )
    obpos.append(x)    
    y = obsradius*np.sin( (np.pi / 180.0) * longitude )
    obpos.append(y)
    z = R0
    obpos.append(z)
    return obpos

def get_barycent_cor(params, im_head, obnames, reffile):
    """
    Calculates the barycentric correction.
    To do this, position of the telescope and pointing of the telescope need to be known.
    Header vaulues are read, if they are not available, then using 
    
    """
    # Define the header keywords, if available
    site_keys       = ['TELESCOP']
    altitude_keys   = ['HIERARCH ESO TEL GEOELEV', 'ESO TEL GEOELEV']       # HIERARCH will be removed from header keywords in python
    latitude_keys   = ['HIERARCH ESO TEL GEOLAT', 'ESO TEL GEOLAT']
    longitude_keys  = ['HIERARCH ESO TEL GEOLON', 'ESO TEL GEOLON']
    ra_keys         = ['RA']
    dec_keys        = ['DEC']
    epoch_keys      = ['HIERARCH ESO TEL TARG EQUINOX', 'ESO TEL TARG EQUINOX']
    pmra_keys       = ['HIERARCH ESO TEL TARG PMA', 'ESO TEL TARG PMA']
    pmdec_keys      = ['HIERARCH ESO TEL TARG PMD', 'ESO TEL TARG PMD']
    params['epoch'] = 2000.0
    params['pmra']  = 0.
    params['pmdec'] = 0.
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
    for [i, header_key_words, parentr] in settings:
        if parentr == 'longitude':                  # Enough information to calculate the ephemerides of sun and moon.
            gobs = ephem.Observer()  
            gobs.name = copy.copy(params['site'])
            gobs.lat  = rad(params['latitude'])     # lat/long in decimal degrees  
            gobs.long = rad(params['longitude'])
            obsdate, obsdate_float, exposure_time, obsdate_begin, exposure_fraction = get_obsdate(params, im_head)               # in UTC, mid of the exposure
            gobs.date = obsdate.strftime('%Y-%m-%d %H:%M:%S')
            mephem    = ephem.Moon()
            mephem.compute(gobs)
            sephem    = ephem.Sun()
            sephem.compute(gobs)
            jephem    = ephem.Jupiter()
            jephem.compute(gobs)
            params['ra'] = -999                                             # If not changed this means the object coordinates were made up
            for obname in obnames:
                if obname.lower().find('sun') == 0:
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

    mjd,mjd0, mjd_begin = mjd_fromheader(params, im_head)
    
    ra2, dec2, epoch, pmra, pmdec, obnames = getcoords_from_file(obnames, mjd, filen=reffile)     # obnames will be a list with only one entry: the matching entry 
    if ra2 !=0 and dec2 != 0:
        params['ra']  = ra2
        params['dec'] = dec2
        params['epoch'] = epoch
        params['pmra'] = pmra
        params['pmdec'] = pmdec
        source_radec = 'The object coordinates are derived from the reference file: {0}.'.format(reffile)
    
    ra, dec = params['ra'], params['dec']
    
    bcvel_baryc, bjd = 0.0, 0.0
    if source_radec != 'Warn: The object coordinates were made up!':
        """# Using jplephem from CERES pipeline
        iers          = JPLiers( baryc_dir, mjd-999.0, mjd+999.0 )      # updates the iers.tab file so it contains only 2000 entries
        obsradius, R0 = JPLR0( params['latitude'], params['altitude'])
        obpos         = obspos( params['longitude'], obsradius, R0 )
        man_jplephem.set_ephemeris_dir( baryc_dir , ephemeris )
        man_jplephem.set_observer_coordinates( obpos[0], obpos[1], obpos[2] )
        res         = man_jplephem.doppler_fraction(ra/15.0, dec, int(mjd), mjd%1, 1, 0.0)
        lbary_ltopo = 1.0 + res['frac'][0]
        bcvel_baryc = ( lbary_ltopo - 1.0 ) * (Constants.c/1000.)       # c in m/s"""
        
        # Using barycorrpy (https://github.com/shbhuk/barycorrpy), pip install barycorrpy
        site = ''
        if params['site'] in barycorrpy.EarthLocation.get_site_names():
            site = params['site']
        ephemeris2='https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp'
        bcvel_baryc = barycorrpy.get_BC_vel(JDUTC=mjd+mjd0,ra=ra,dec=dec,obsname=site,lat=params['latitude'],longi=params['longitude'],alt=params['altitude'],pmra=params['pmra'],
                           pmdec=params['pmdec'],px=0,rv=0.0,zmeas=0.0,epoch=params['epoch'],ephemeris=ephemeris2,leap_update=True)
        bcvel_baryc = bcvel_baryc[0][0] / 1E3       # in km/s
        #print "result3", bcvel_baryc
        #print "result3 RV", bcvel_baryc[0][0]
        
        im_head['HIERARCH EXO_PIPE RA'] = (round(ra,6),           'RA in degrees, used to calculated BCV, BJD')
        im_head['HIERARCH EXO_PIPE DEC'] = (round(dec,6),         'DEC in degrees, used to calculated BCV, BJD')
        im_head['HIERARCH EXO_PIPE BCV'] = (round(bcvel_baryc,4), 'barycentric velocity in km/s')
        im_head['HIERARCH EXO_PIPE PMRA'] = (round(pmra,3),       'proper motion for RA in mas/yr, for BCV, BJD')
        im_head['HIERARCH EXO_PIPE PMDEC'] = (round(pmdec,3),     'proper motion for DEC in mas/yr, for BCV, BJD')
        
        """# Using jplephem
        bjd = jd_corr(mjd, ra, dec, params['epoch'], params['latitude'], params['longitude'], jd_type='bjd')
        bjd = bjd[0]"""
        
        # Using barycorrpy (https://github.com/shbhuk/barycorrpy), pip install barycorrpy
        bjdtdb = barycorrpy.utc_tdb.JDUTC_to_BJDTDB(JDUTC=mjd+mjd0,ra=ra,dec=dec,obsname=site,lat=params['latitude'],longi=params['longitude'],alt=params['altitude'],pmra=params['pmra'],
                                        pmdec=params['pmdec'],px=0,rv=0.0,epoch=params['epoch'],ephemeris=ephemeris2,leap_update=True)           # only precise to 0.2s 
        bjd = bjdtdb[0][0]
        #print "bjdtdb", bjdtdb
        #print "bjdtdb0", bjdtdb[0][0]
        
        # end test
        im_head['HIERARCH EXO_PIPE BJDTDB'] = (round(bjd,5), 'Baryc. cor. JD (incl leap seconds)')     # without leap seconds: remove 32.184+N leap seconds after 1961'

    im_head['HIERARCH EXO_PIPE JD']         = (round(mjd+mjd0,5), 'mid-exp Julian date')             # MJD = JD - 2400000.5 from http://www.csgnetwork.com/julianmodifdateconv.html
    im_head['HIERARCH EXO_PIPE MJD']        = (round(mjd,5), 'mid-exp modified JD')        # round 5 -> precision is 1 second, timing is not more precise
    im_head['HIERARCH EXO_PIPE MJD-START']  = (round(mjd_begin,5), 'mid-exp modified JD')

    logger(('Info: Using the following data for object name(s) {10}, Observatory site {9}, mid exposure MJD {11}: '+\
                    'altitude = {0}, latitude = {1}, longitude = {2}, ra = {3}, dec = {4}, epoch = {5}, pmra = {13}, pmdec = {14}. {8} {6} '+\
                    'This leads to a barycentric velocity of {7} km/s and a mid-exposure BJD-TDB of {12}').format(params['altitude'], params['latitude'], params['longitude'], 
                         round(ra,6), round(dec,6), params['epoch'], source_radec, round(bcvel_baryc,4), source_obs, 
                         params['site'], str(obnames).replace("'","").replace('"',''), mjd, round(bjd,5), params['pmra'], params['pmdec'] ))
       
    return params, bcvel_baryc, mephem, obnames, im_head

def JPLiers(path, mjdini, mjdend):                  # from CERES, updates the iers.tab file so it contains only 2000 entries, # not needed anymore
    output    = open(path+'iers.tab','w')
    filename  = path+'finals2000A.data'
    finaldata = open(filename,'r')

    for line in finaldata:
        mj = line[7:15]
        if float(mj) >= float(mjdini) and float(mj) <= float(mjdend) and len(line.split()) > 5:    
            c1 = line[18:27]
            c2 = line[37:46]
            c3 = line[58:68]
            l  = ' '+mj+' '+c1+' '+c2+' '+c3+' '+'\n'
            output.write(l)
            if mj == mjdini+999: print("estoy en el dia D")
        if float(mj) > float(mjdend):
            break
    finaldata.close()
    output.close()
    pass

def jd_corr(mjd, ra, dec, epoch, lat, lon, jd_type='bjd'):          # not used anymore
    """
    Return BJD or HJD for input MJD(UTC).
    Adapted from https://mail.python.org/pipermail/astropy/2014-April/002843.html
    Ignores the position of the observatory at the moment, but this will change it by less than 0.1s
    Comparing to http://astroutils.astronomy.ohio-state.edu/time/utc2bjd.html gives a difference of 69s, which is probably due to using BJD_TDB there and BJD_UTC here (32.184+37 leap seconds (Feb 2019))
    return new_jd.jd: 1d array of floats, same length as mjd
    """
    # Initialise ephemeris from jplephem
    eph = jplephem.Ephemeris(de423)

    # Source unit-vector
    ## Set distance to unit (kilometers)
    src_vec = astcoords.SkyCoord(ra=ra*astunits.degree, dec=dec*astunits.degree, frame=astcoords.FK5(equinox='J{0}'.format(epoch)), distance=1*astunits.km)
    # Convert epochs to astropy.time.Time
    ## Assume MJD(UTC)
    t = asttime.Time(mjd, scale='utc', format='mjd')#, lat=lat, lon=lon)    # Lat, lon not supported anymore?

    # Get Earth-Moon barycenter position
    ## NB: jplephem uses Barycentric Dynamical Time, e.g. JD(TDB)
    ## and gives positions relative to solar system barycenter
    barycenter_earthmoon = eph.position('earthmoon', t.tdb.jd)

    # Get Moon position vectors
    moonvector = eph.position('moon', t.tdb.jd)

    # Compute Earth position vectors
    pos_earth = (barycenter_earthmoon - moonvector * eph.earth_share)*astunits.km

    if jd_type == 'bjd':
        # Compute BJD correction
        ## Assume source vectors parallel at Earth and Solar System Barycenter
        ## i.e. source is at infinity
        corr = np.dot(pos_earth.T, src_vec.cartesian.xyz)/astconst.c            # light travel time in seconds
    elif jd_type == 'hjd':
        # Compute HJD correction via Sun ephemeris
        pos_sun = eph.position('sun', t.tdb.jd)*astunits.km
        sun_earth_vec = pos_earth - pos_sun
        corr = np.dot(sun_earth_vec.T, src_vec.cartesian.xyz)/astconst.c        # light travel time in seconds

    # TDB is the appropriate time scale for these ephemerides
    dt = asttime.TimeDelta(corr, scale='tdb', format='jd')                      # light travel time in days
    t.format='jd'
    # Compute and return HJD/BJD as astropy.time.Time
    new_jd = t + dt

    return new_jd.jd            # Return as float-array

def find_shift_between_wavelength_solutions(wave_sol_1, wave_sol_lines_1, wave_sol_2, wave_sol_lines_2, spectra, names=['first','second']):     # Doesn't do the job as expected, it's a pixel shift, not an RV shift
    """
    :param wave_sol_1, wave_sol_2: 2d array of floats, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    :param wave_sol_lines_1, wave_sol_lines_2: 2d array of floats with one line for each order. Each order contains the wavelengths of the identified reference lines, 
                                        sorted by the brightest to faintest (in spectrum from which the solution is derived).
                                        To create the same number of lines for each order, the array is filled with 0
    :param spectra: 2d array of floats, spectrum, only needed for the shape of the wavelengths array
    """
    # make it into wavelength, find the wavelengths of the first arclines in both, find the wavelengths of the second arclines in both, find the overlapping lines -> directly to lambda-shift -> RV
    # Make the solutions into wavelengths
    wavelengths1 = create_wavelengths_from_solution(wave_sol_1, spectra)
    wavelengths2 = create_wavelengths_from_solution(wave_sol_2, spectra)
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
        if i==0:
            good_values = range(len(result))        # all data
        if i==1:
            good_values = result[:,1] == 0          # only same lines in both solutions
        if i==2:
            good_values = result[:,1] == 1          # lines of solution1
        if i==3:
            good_values = result[:,1] == 2          # lines of solution2
        if i==4:
            good_values = (result[:,1] == 0) & (result[:,5] < 500)          # only same lines in both solutions and left side 
        if i==5:
            good_values = (result[:,1] == 0) & (result[:,5] > 1500)          # only same lines in both solutions and right side 
        if i==6:
            good_values = (result[:,1] == 0) & (result[:,5] > 700) & (result[:,5] < 1300)          # only same lines in both solutions and middle
        if i==7:
            good_values = (result[:,1] == 0) & (result[:,5] > 700) & (result[:,5] < 1300) & (result[:,3] < 20)          # only same lines in both solutions and middle and red orders
        if i==8:
            good_values = (result[:,1] == 0) & (result[:,5] > 700) & (result[:,5] < 1300) & (result[:,3] > 60)          # only same lines in both solutions and middle and blue orders
        if i==9:
            good_values = (result[:,1] == 0) & (result[:,5] > 700) & (result[:,5] < 1300) & (result[:,3] > 30) & (result[:,3] < 50)         # only same lines in both solutions and middle and middle orders
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


def rv_analysis(params, spec, im_head, fitsfile, obname, reffile, mephem):
    base = params['path_ceres']
    # Import the routines from Ceres:
    sys.path.append(base+"utils/Correlation")
    sys.path.append(base+"utils/GLOBALutils")
    sys.path.append(base+"utils/OptExtract")    # for Marsh, at least
    sys.path.append(base+"utils/CCF")           # needed by GLOBALutils.py
    import GLOBALutils
    import correlation
    # Import other stuff
    import statsmodels.api as sm
    lowess = sm.nonparametric.lowess
    import pickle

    def get_spec(spectra, waves, cen_wave, range_px):
        """
        :param spectra: 1d or 2d array of floats with the spectral data
        :param waves: same format as spectra, corresponding wavelengths
        :param cen_wave: wave_length to get the spectra around
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

    specs = spec.shape
    spec = np.vstack(( spec, np.zeros([11-specs[0], specs[1], specs[2]]) ))
    spec[7:11,:,:] = np.nan
    lbary_ltopo = 1
    npools = 1
    refvel = 0
    know_moon = False
    here_moon = False
    models_path = base + 'data/COELHO_MODELS/R_40000b/' # "/home/ronny/software/ceres-master/data/COELHO_MODELS/R_40000b/"
    #dirout = './'          # replaced by params['path_rv_ceres']
    fsim = fitsfile
    RESI = 120000.
    force_stellar_pars = False
    
    
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
    for order in range(spec.shape[1]):
        L  = np.where( (spec[1,order,:] != 0) & (~np.isnan(spec[1,order,:])) )              # good values
        #ratio              = np.polyval(ccoefs[order],spec[0,order,:][L])*Rnorms[order]
        np.warnings.filterwarnings('ignore')
        ratio = spec[1,order,:][L] / spec[5,order,:][L]                                     # ratio between extracted spectrum and continuum normalised spectrum -> blaze function, cancels absorption lines
        np.warnings.resetwarnings()
        spec[7,order,:][L] = ratio
        #spec[8,order,:][L] = spec[6,order,:][L]                                             # error continuum (first guess), but not good. #sn_order=8
        spec[8,order,:][L] = spec[2,order,:][L]                                             # error of the extracted data, depending on what is used, the RV changes by few 100 km/s -> > several \AA
        #spec[8,order,:][L] = ratio * R_flat_ob_n[order,1,:][L] / np.sqrt( ratio * R_flat_ob_n[order,1,:][L] / gain2 + (ron2/gain2)**2 )        # something like S/N -> used as this by XCor
        #spl           = scipy.interpolate.splrep(np.arange(WavSol.shape[0]), WavSol,k=3)
        #dlambda_dx    = scipy.interpolate.splev(np.arange(WavSol.shape[0]), spl, der=1)
        #NN            = np.average(dlambda_dx)
        #dlambda_dx    /= NN
        np.warnings.filterwarnings('ignore')
        LL = np.where(spec[5,order,:] > 1 + 10. / scipy.signal.medfilt(spec[8,order,:],21))[0]          # remove emission lines and cosmics
        np.warnings.resetwarnings()
        spec[5,order,LL] = 1.
        spec[9,order,:][L] = spec[5,order,:][L]# * (dlambda_dx[L] ** 1)         # used for the analysis in XCor (spec_order=9, iv_order=10)
        spec[10,order,:][L] = spec[2,order,:][L]# / (dlambda_dx[L] ** 2)        # used for the analysis in XCor (spec_order=9, iv_order=10)
    #plot_img_spec.plot_spectra_UI(np.array([spec]))
    T_eff, logg, Z, vsini, vel0 = 5777, 4.4374, 0.0134, 2, 0
    if True:  
        if np.nanmax(spec[0,:,:], axis=None) > 5500:       # only run for wavelengths bigger than 5500A, as otherwise problems in correlation.CCF
            #pars_file = params['path_rv_ceres'] + fsim.split('/')[-1][:-4]+'_stellar_pars.txt'
            pars_file = params['path_rv_ceres'] + obname+'_stellar_pars.txt'                  # same stellar parameters for same object
            pars_file = params['path_rv_ceres'] + fsim+'_stellar_pars.txt'                    # calculate stellar parameters for each spectrum

            if os.access(pars_file,os.F_OK) == False or force_stellar_pars:
                Rx = np.around(1./np.sqrt(1./40000.**2 - 1./RESI**2))
                spec2 = spec.copy()
                for i in tqdm(range(spec.shape[1]), desc="Step: Estimating atmospheric parameters"):
                    IJ = np.where(spec[5,i]!=0.)[0]
                    spec2[5,i,IJ] = GLOBALutils.convolve(spec[0,i,IJ],spec[5,i,IJ],Rx)
                try:
                    T_eff, logg, Z, vsini, vel0, ccf = correlation.CCF(spec2,model_path=models_path,npools=npools, base=base+'utils/Correlation/')     # Fails, because our spectrum doesn't cover the hard coded wavelength
                    line = "%6d %4.1f %4.1f %8.1f %8.1f\n" % (T_eff,logg, Z, vsini, vel0)
                    f = open(pars_file,'w')
                    f.write(line)
                    f.close()
                    loadtxt = ''
                except:
                    loadtxt = 'Warn: could not determine the stelar parameters. This is probably caused because of the availble wavelength range. Do not worry about it'
                    logger(loadtxt)
            else:
                T_eff, logg, Z, vsini, vel0 = np.loadtxt(pars_file,unpack=True)
                loadtxt = ' (Atmospheric parameters loaded from file {0})'.format(pars_file)
            if loadtxt.find('Warn') != 0:
                logger('Info: Using the following atmosperic parameters for T_eff, logg, Z, vsini, vel0: {1}, {2}, {3}, {4}{0}'.format(loadtxt, T_eff, logg, Z, vsini, vel0))
        
        logger('Step: Radial Velocity analysis:')
        # assign mask
        #   obname is the name of the object
        #   reffile is a reference file: reffile = dirin+'reffile.txt' -> this file doesn't exist
        sp_type, mask = GLOBALutils.get_mask_reffile(obname,reffile=reffile,base=base+'data/xc_masks/')     # !!! Warn: upper and lower case matters
        logger("Will use {0} mask for CCF.".format(sp_type))

        # Read in mask
        ml, mh, weight = np.loadtxt(mask,unpack=True)
        ml_v = GLOBALutils.ToVacuum( ml )
        mh_v = GLOBALutils.ToVacuum( mh )
        av_m = 0.5*( ml_v + mh_v )
        mask_hw_kms = (GLOBALutils.Constants.c/1e3) * 0.5*(mh_v - ml_v) / av_m

        disp = GLOBALutils.get_disp(obname, reffile=reffile)        # !!! Warn: upper and lower case matters disp is "velocity width in km/s that is used to broaden the lines of the binary mask. It should be similar to the standard deviation of the Gaussian that is fitted to the CCF."
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

        logger('Computing the CCF...')
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
            vels, xc_full, sn, nlines_ccf, W_ccf = \
                    GLOBALutils.XCor(spec, ml_v, mh_v, weight,\
                    0, lbary_ltopo, vel_width=300, vel_step=3,\
                    spec_order=9, iv_order=10, sn_order=8, max_vel_rough=300)
            
            # W_ccf is a weigth
            
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

        #pkl_xc = params['path_rv_ceres'] + fsim.split('/')[-1][:-4]+obname+'_XC_'+sp_type+'.pkl'
        pkl_xc = params['path_rv_ceres'] + fsim+'_XC_'+sp_type+'.pkl'      # without filename ending
        pickle.dump( xc_dict, open( pkl_xc, 'w' ) )

        #ccf_pdf = params['logging_path'] + fsim.split('/')[-1][:-4] + obname + '_XCs_' + sp_type + '.pdf'       # dirout + 'logging/'
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
        
        B,A = -0.00257864,0.07765779            # from Monte Carlo Simulation, different for each instrument
        B,A = 0.005, 0.2
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

        B,A = -0.00348879, 0.10220848
        BSerr = B + ( 1.6 + 0.2 * p1gau[2] ) * A / np.round(SNR_5130)
        if BSerr<0.002:
            BSerr = .002

        RV     = np.around(p1gau_m[1],4)  
        BS     = np.around(SP,4)   
        RVerr2 = np.around(RVerr,4)
        BSerr  = np.around(BSerr,4)

        return RV, RVerr2, BS, BSerr













