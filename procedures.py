#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
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

calimages = dict()  # dictionary for all calibration images used by create_image_general and read_file_calibration
oldorder = 12

def logger(message, show=True, printarrayformat=[], printarray=[], logfile='logfile'):
    """
    Saves the status information to a logfile
    :param message: Text to log
    :param show: if False than it will be logged only to the logfile but not printed on the screen
    :param printarrayformat: Format of the columns in printarray, e.g. '%3.1f'
    :param printarray: Array with values to log
    :param logfile: filename in which the information will be written
    """
    if show:
        print message
    file = open(logfile, 'a')
    file.write(time.strftime("%Y%m%d%H%M%S - ", time.localtime()) + message + '\n')
    if printarrayformat <> [] and printarray <> []:
        for line in printarray:
            text = '\t'
            for i,printformat in enumerate(printarrayformat):
                #print printformat,line[i]
                text += '\t'+printformat%line[i]
            file.write(text+'\n')
            if show:
                print text
    file.close()
    if message.find('Error') == 0:
        print '\t-> exiting'
        exit(1)

def log_params(params):
    """
    formats the dictionary to be saved in the logfile
    """
    logger('params: '+json.dumps(params,sort_keys = False, indent = 4), show=False, logfile='logfile_params')

def read_parameterfile(textfile):
    # load text file (remove all white spaces)
    if not os.path.exists(textfile):
        logger('Error: The parameterfile {0} does not exist -> exit'.format(textfile))
    data = np.genfromtxt(textfile, dtype=str, comments='#', delimiter='=')
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
        else:
            logger('Warn: I dont know how to handle command line argument {0}'.format(arg))
    
    params.update(cmdparams)
    
    # deal with paths, create folders
    for entry in params.keys():
         if entry.find('path') <> -1:
            params[entry] = (params[entry]+'/').replace('//', '/')
            if not os.path.exists(params[entry]):
                try:
                    os.makedirs(params[entry])
                except:
                    logger('Warn: Cant create directory {0}'.format(params[entry]))
    
    # deal with lists
    lists = ['subframe', 'arcshift_range', 'order_offset', 'px_offset', 'px_offset_order', 'use_catalog_lines', 'polynom_order_traces', 'polynom_order_intertraces', 'opt_px_range', 'bin_search_apertures', 'bin_adjust_apertures', 'raw_data_file_endings', 'background_width_multiplier', 'polynom_bck']
    for entry in params.keys():
        if entry.find('_rawfiles') > 0 or entry.find('calibs_') > -1:
            lists.append(entry)
    for entry in lists:
        temp = params[entry]
        for i in ['[', ' ', ']']:
            temp = temp.replace(i,'')
        if len(temp) == 0:
            temp = []
        else:
            temp = temp .split(',')
        params[entry] = temp
    
    # deal with lists of integers
    lists = ['subframe', 'arcshift_range', 'order_offset', 'px_offset', 'px_offset_order', 'polynom_order_traces', 'polynom_order_intertraces', 'bin_search_apertures', 'bin_adjust_apertures', 'polynom_bck']
    for entry in lists:
        if params[entry] <> ['']:           # only make useful data into integers
            for i in range(len(params[entry])):
                try:
                    params[entry][i] = int(params[entry][i])
                except:
                    emsg2 = 'Parameter "{0}" (value of "{1}")'.format(entry, params[entry])
                    logger(emsg + emsg2 + ' must be a list of integers')
            
    # deal with lists of floats
    lists = ['opt_px_range', 'background_width_multiplier']
    for entry in lists:
        if params[entry] <> ['']:           # only make useful data into floats
            for i in range(len(params[entry])):
                try:
                    params[entry][i] = float(params[entry][i])
                except:
                    emsg2 = 'Parameter "{0}" (value of "{1}")'.format(entry, params[entry])
                    logger(emsg + emsg2 + ' must be a list of floats')
    
    # deal with ints
    floats = ['polynom_order_apertures', 'rotate_frame']
    for entry in floats:
        try:
            params[entry] = int(params[entry])
        except:
            emsg2 = 'Parameter "{0}" (value of "{1}")'.format(entry, params[entry])
            logger(emsg + emsg2 + ' must be an integer')
            
    # deal with floats
    floats = ['max_good_value', 'catalog_file_wavelength_muliplier', 'extraction_width_multiplier', 'arcextraction_width_multiplier', 'resolution_offset_pct', 'diff_pxs', 'maxshift', 'wavelength_scale_resolution']
    for entry in floats:
        if entry not in params.keys():
            emsg2 = 'Parameter "{0}" '.format(entry)
            logger(emsg + emsg2 + ' missing')
        try:
            params[entry] = float(params[entry])
        except:
            emsg2 = 'Parameter "{0}" (value of "{1}")'.format(entry, params[entry])
            logger(emsg + emsg2 + ' must be a float')
    
    # deal with True/False:
    bools = ['flip_frame', 'update_widths', 'GUI']
    trues = ['yes', 'true', '1']
    for entry in bools:
        if params[entry].lower() in trues:
            params[entry] = True
        else:
            params[entry] = False
    
    """# deal with order direction
    svar = 'order_direction'
    sargs = ['h', 'horizontal','v', 'vertical']
    sargs2 = 'H, h, Horizontal, V, v, or Vertical'
    params['transpose'] = False
    if svar in params:
        if params[svar].lower() not in sargs:
            emsg2 = '"{0}" (value of "{1}")'.format(svar, params[svar])
            raise Exception(emsg + emsg2 + ' must contain: ' + sargs2)
        elif params[svar].lower() in sargs[2:4]:
            params['transpose'] = True"""
    
    # deal with arc orders
    params['arcshift_range'] = list(np.abs(params['arcshift_range']))
    svar = 'arcshift_side'
    if params[svar].lower() in ['left','l']:
        params[svar] = -1
    elif params[svar].lower() in ['right','r']:
        params[svar] = 1
    elif params[svar].lower() in ['center','c']:
        params[svar] = 0
    else:
        emsg2 = 'Parameter "{0}" (value of "{1}")'.format(svar, params[svar])
        logger(emsg + emsg2 + ' must be one of the following: left, l, right, r, center, c')
    
    # deal with lists of raw data filenames -> add path
    filenlists = []
    for entry in params.keys():
        if entry.find('_rawfiles') > 0:
            filenlists.append(entry)
    for entry in filenlists:
        for i in range(len(params[entry])):
            params[entry][i] = params['raw_data_path'] + params[entry][i]
    
    # deal with result filenames/folders -> add path
    filenames = ['path_extraction', 'logging_path']     #, 'configfile_fitsfiles' (excluded, as should be handled as conf.txt
    for entry in params.keys():
        if (entry.find('master_') == 0 or entry.find('background_') == 0) and entry.find('_filename') > 0:
            filenames.append(entry)
    for entry in filenames:
        params[entry] = params['result_path'] + params[entry]
        
    # deal with logging filenames/folders -> add path
    filenames = []
    for entry in params.keys():
        if entry.find('logging_') == 0 and entry.find('logging_path') == -1:
            filenames.append(entry)
    for entry in filenames:
        params[entry] = params['logging_path'] + params[entry]
    
    # deal with full filenames -> nothing to do
    filenames = ['badpx_mask_filename', 'original_master_order_filename']
    
    return params

def update_calibration_memory(key,value):
    global calimages
    calimages[key] = value

def convert_readfile(input_list, textformat, delimiter='\t', replaces=[]):
    result_list = []
    for entry in input_list:
        for replce in replaces:
            entry = entry.replace(replce,'')
        entry = entry.split(delimiter)
        if len(entry) < len(textformat):
            continue
        for i in range(len(textformat)):
            entry[i] = textformat[i](entry[i])
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
        if subframe <> [] and badpx_mask.shape <> (subframe[0],subframe[1]):
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
    if len(ims) <> len(names):
        print 'len(ims) <> len(names), that seems like a coding error'
    problems = ''
    for i in range(len(ims)-1):
        for j in range(i,len(ims)):
            if ims[i].shape <> ims[j].shape:
                problems += '\tImage: {0} ({2}) and Image {1} ({3})\n'.format(names[i], names[j], ims[i].shape, ims[j].shape)
    if problems <> '':
        logger('Error: The following images have not the same size, but should have. This is most likely caused by a missing "subframe" in one of the parameters "calib*". Please check.\n {0}'.format(problems[:-1]) )
    return

def read_file_calibration(params, filename, level=0):
    """
    Reads the filename and applies the calibrations as given in params['calibs']
    :param params: Dictionary with all the parameters
    :param filename: path and name of the file
    :return: 2D numpy array of the file, and the header of the file read
    """
    global calimages
    if os.path.isfile(filename) == False:
        logger('Error: File {0} is missing'.format(filename))
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
            logger('Error: The file is stored in a multi-demensional array, which I do not know how to handle. The size of the image is {0}. This requires a small adjustment to the code in procedure read_file_calibration'.format(ims))
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
            if im.shape <> (subframe[0],subframe[1]):                   # only apply subframe if the file doesn't have the size already
                im = im[subframe[2]: subframe[0]+subframe[2], subframe[3]: subframe[1]+subframe[3]]
            logger('Info: {1}: subframe applied: {0}'.format(entry, level))
            im_head['redu{0}a'.format(level)] = 'Subframe: {0}'.format(entry)
        elif entry.lower().find('bias') > -1:
            if entry not in calimages:
                 create_image_general(params, entry, level=level+1)
            warn_images_not_same([im, calimages[entry]], [filename,entry])
            im = im - calimages[entry]
            logtxt, headtxt = ['bias correction applied'], ['redu{0}b'.format(level), 'Bias']
        elif entry.lower().find('dark') > -1:
            if entry == 'dark':             #only add exposure time, if just dark is given
                entry = entry+str(exptime)
            if entry not in calimages:
                 create_image_general(params, entry, level=level+1)
            warn_images_not_same([im, calimages[entry]], [filename,entry])
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
            im_head['redu{0}e'.format(level)] = 'Background: {0}'.format(params[entry])
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
                    section = section[section<>0]           #only use areas <> 0
                    if len(section) >= 2:
                        break
                if len(section) == 0:
                    logger('Warn: cannot replace bad pixel ({0}, {1}) with surrounding area in {2}'.format(nonzeroind[0][i],nonzeroind[1][i],filename))
                else:
                    im[nonzeroind[0][i],nonzeroind[1][i]] = np.median(section)  #replace bad px with the median of each surrounding area
            logger('Info: {1}: badpx correction applied: {0}'.format(entry, level))
            im_head['redu{0}f'.format(level)] = 'Bad-pixel-mask: {0}'.format(entry)
        elif entry.lower().find('localbackground') > -1:
            if 'sci_trace' in calimages.keys() and 'cal_trace' in calimages.keys():
                logger('Step: Performing the background fit')
                sci_tr_poly, xlows, xhighs, widths = calimages['sci_trace']
                cal_tr_poly, axlows, axhighs, awidths = calimages['cal_trace']
                bck_px_sci = find_bck_px(im, sci_tr_poly, xlows, xhighs, widths, 3)#params['background_width_multiplier'])
                bck_px_cal = find_bck_px(im, cal_tr_poly, axlows, axhighs, awidths, 1)#params['background_width_multiplier'])
                bck_px = bck_px_sci * bck_px_cal
                bad_values = ( im*bck_px > np.percentile(im[bck_px==1],95) )
                bck_px[bad_values] = 0
                
                # Some deviation at the red side with not many lines, computational heavy
                bck_im = fit_2d_image(im, params['polynom_bck'][1], params['polynom_bck'][0], w=bck_px)
                #plot_img_spec.plot_image(im, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'orig')
                #plot_img_spec.plot_image(bck_px, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'bck_px')
                #plot_img_spec.plot_image(bck_im, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'bck_im')
                #plot_img_spec.plot_image(bck_im-im, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'diff')
                plot_img_spec.plot_image((im - bck_im)*bck_px, [params['logging_path']+filename.replace(params['raw_data_path'],'background_subtracted-')+'.png'], 1, False, [0.05,0.95,0.95,0.05], 'difference between image and background fit')
                im = im - bck_im
                logger('Info: {1}: background correction applied: {0}'.format(entry, level))
                im_head['redu{0}f'.format(level)] = 'Background: {0}'.format(entry)
            else:
                logger('Warn: Could not apply the calibration step {0} because the science and/or calibration traces are not yet known.'.format(entry))
        elif entry.lower().find('combine_sum') > -1 or entry.lower().find('normalise') > -1:
            'nothing to do, as for a different step'
        else:
            logger('Warn: do not know what to do with this correction: {0}'.format(entry))
        if logtxt <> [] and headtxt <> []:
            im_median, im_std = int(round(np.median(calimages[entry]))), int(round(np.std(calimages[entry], ddof=1)))
            logger('Info: {4}: {3}: {0} (median={1}, std={2})'.format(entry, im_median, im_std, logtxt[0], level))
            im_head[headtxt[0]] = '{3}: {0}, median={1}, std={2}'.format(entry, im_median, im_std, headtxt[1])
            
    #logger('Info: image loaded and processed: {0}'.format(filename))
    return im, im_head

def create_image_general(params, imtype, level=0):
    """
    Reads or creates the imtype file. If the key and file for master_<imtype>_filename exists the file is read, otherwise the file is created by combining the <imtype>_rawfiles
    :param params: Dictionary with all the parameters
    :param imtype: type of images, e.g. flat, dark5.0, bias
    :return im: 2d array of the combined file
    :return im_head: header of the last read file
    """
    global calimages
    loaded = False
    if 'master_{0}_filename'.format(imtype) in params.keys():
        if params['master_{0}_filename'.format(imtype)] <> '' and os.path.isfile(params['master_{0}_filename'.format(imtype)]) == True:
            logger('Info: Using exiting {0}: {1}'.format(imtype,params['master_{0}_filename'.format(imtype)]))
            params['calibs'] = params['calibs_read']
            if '{0}_calibs_read'.format(imtype) in params.keys():
                params['calibs'] = params['{0}_calibs_read'.format(imtype)]
            im, im_head = read_file_calibration(params, params['master_{0}_filename'.format(imtype)], level=level)
            loaded = True
    if loaded == False:
        if '{0}_calibs_create'.format(imtype) not in params.keys():
            if 'standard_calibs_create' not in params.keys():
                logger('Error: Missing entry in the configuration file. Neigther "{0}_calibs_create" nor "standard_calibs_create" is given. Please update the configuration file.'.format(imtype))
            params['{0}_calibs_create'.format(imtype)] = params['standard_calibs_create']
        for i in range(len(params['{0}_calibs_create'.format(imtype)])):                                                                    # make it safe from different user input
            params['{0}_calibs_create'.format(imtype)][i] = params['{0}_calibs_create'.format(imtype)][i].lower()
            params['{0}_calibs_create'.format(imtype)][i] = params['{0}_calibs_create'.format(imtype)][i].replace('normaliz', 'normalis')
            params['{0}_calibs_create'.format(imtype)][i] = params['{0}_calibs_create'.format(imtype)][i].replace('normalisation', 'normalise')
        im, med_fluxes, std_fluxes = [], [], []
        head_variation = [[params['raw_data_exptim_keyword']], [params['raw_data_dateobs_keyword']], ['JD-HELIO']]
        if '{0}_rawfiles'.format(imtype) not in params.keys():
            logger('Error: The list of raw files for image type {0} is not defined in the configuration. Please check the configuration files.'.format(imtype))
        if len(params['{0}_rawfiles'.format(imtype)]) > 145:                    # make more general, divide expected memory usage by 0.44, as test have shown that im.nbytes/memory will be about 0.45, when ram is fully used
            logger('Warn: The ammount of pictures will cause swapping')                                              # will cause swapping eventually
            prec = np.float16
        elif len(params['{0}_rawfiles'.format(imtype)]) > 70:
            prec = np.float16                                            # only about 145 files of size 4250,2838 before swapping starts on a 8GB machine
        elif len(params['{0}_rawfiles'.format(imtype)]) > 35:
            prec = np.float32                                            # only about 70 files of size 4250,2838 before swapping starts on a 8GB machine
        else:
            prec = np.float64                                            # only about 35 files of size 4250,2838 before swapping starts on a 8GB machine
        for imf in params['{0}_rawfiles'.format(imtype)]:                   # Only works for maximum 40 images on neils machine
            params['calibs'] = params['{0}_calibs_create'.format(imtype)]   #get's overwritten when other files are being read
            img, im_head = read_file_calibration(params, imf, level=level)
            med_flux = np.median(img, axis=None)
            med_fluxes.append(med_flux)
            std_fluxes.append(np.std(img, axis=None, ddof=1))
            if 'normalise' in params['{0}_calibs_create'.format(imtype)]:
                img = img/(med_flux+0.0)
            if im == []:
                # -> create the full array at once, insertering len(rawfiles) into img.shape at 0
                im = np.array([img], dtype=prec)
            else:
                im = np.append(im, np.array([img], dtype=prec), axis=0)
            for i in range(len(head_variation)):
                value = 0
                if head_variation[i][0] in im_head.keys():
                    value = im_head[head_variation[i][0]]
                head_variation[i].append(value)
            #print im.dtype, im.itemsize, im.nbytes, sys.getsizeof(im), im.nbytes/7979408000.
        for i in range(len(med_fluxes)):
            im_head['NORM_{0}'.format(i)] = med_fluxes[i]
        for i in range(len(std_fluxes)):
            im_head['STDV_{0}'.format(i)] = std_fluxes[i]
        if 'combine_sum' in params['{0}_calibs_create'.format(imtype)]:
            im = combine_sum(im)
            im_head['redu07'] = 'Sum of {0} images'.format(len(params['{0}_rawfiles'.format(imtype)]))
            im_head[params['raw_data_exptim_keyword']] = sum(head_variation[0][1:])
        else:
            im = combine_median(im)
            im_head['redu07'] = 'Median of {0} images'.format(len(params['{0}_rawfiles'.format(imtype)]))
        if 'normalise' in params['{0}_calibs_create'.format(imtype)]:
            norm_factor = np.median(med_fluxes)
            im = im * norm_factor
            im_head['NORM_MED'] = norm_factor
        im_head['JD-HELIO'] = np.mean(head_variation[2][1:])
        im_head[params['raw_data_dateobs_keyword']] = min(head_variation[1][1:])
        if 'master_{0}_filename'.format(imtype) in params.keys():
            if params['master_{0}_filename'.format(imtype)] <> '':
                save_im_fits(params, im, im_head,  params['master_{0}_filename'.format(imtype)])
    calimages[imtype] = im
    calimages['{0}_head'.format(imtype)] = im_head
    return im, im_head

def get_minimum_data_type(arr, allow_unsigned=True):
    arr = np.array(arr)
    if np.sum(arr - arr.astype(int), axis=None) == 0:
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
    im = rotate_flip_frame(im, params, invert=True)
    #im = get_minimum_data_type(im, allow_unsigned=False)   #That doen't make the fits files smaller for uint16, int16, or float32. Additionally, plotting a file with uint16 or int16 with ds9 or gaia doesn't show the data correctly
    fits.writeto(filename, im, im_head, overwrite=True)
    logger('Info: image saved: {0}'.format(filename))

def read_fits_width(filename):
    """
    read in the polynomial fit fits file
    :param filename: string, location and file name of the file containing the polynomial fits
    fits file should look like the following:
    Order          6                  5          ...       0        low_x   high_x  width_left width_right width_gauss
    float64      float64            float64       ...    float64    float64 float64    float64     float64     float64
    ------- ------------------ ------------------ ... ------------- ------- ------- ---------- ----------- -----------
        1.0  7.80380819018e-20 -7.28455089523e-16 ... 1086.73647399     0.0  3028.0       13.5        11.2      2.234
        2.0  2.09138850722e-19 -2.01802656266e-15 ... 1123.38227429     0.0  3082.0       11.5        12.7      2.546
    ...
    where 6, 5, 4, 3, 2, 1, 0 are the polynomial powers in p
        i.e. p where:
            p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
    :return: array of the parameters for the polynomial fit
    :return: minimum and maximum position of the trace in dispersion direction
    :return: array with the width of the trace in cross-dispersion direction, giving the width to the left, to the right, and the Gaussian width
    """
    # convert to astropy table
    if os.path.isfile(filename) == False:
        logger('Error: File {0} is missing'.format(filename))
    atable = Table.read(filename)
    orders = np.array(atable['Order'], dtype=int)

    pfits, xlows, xhighs, widths = [], [], [], []
    for order in orders:
        pfit = []
        for p in atable.colnames:
            if p not in ['Order', 'low_x', 'high_x', 'width_left', 'width_right', 'width_gaus', 'width_gauss']:  
                pfit.append(atable[p][order-1])        
        pfits.append(pfit)
        xlows.append(int(atable['low_x'][order-1]))
        xhighs.append(int(atable['high_x'][order-1]))
        try:
            widths.append([atable['width_left'][order-1], atable['width_right'][order-1], atable['width_gauss'][order-1]])
        except:
            widths = []
    logger('Info: Master file to trace orders read: {0}'.format(filename))
    return np.array(pfits), np.array(xlows), np.array(xhighs), np.array(widths)

def save_fits_width(pfits, xlows, xhighs, widths, filename):
    """
    save the fits to file, for use in extracting future spectra
    :param pfits: list, length same as number of orders, polynomial values
                      (output of np.polyval)
                      i.e. p where:
                      p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
    :param xlows: list, length same as number of orders, the lowest x pixel
                   (wavelength direction) used in each order
    :param xhighs: list, length same as number of orders, the highest x pixel
                   (wavelength direction) used in each order
    :param widths: 2d list, length same as number of orders, each entry contains left border, right border, and Gaussian width of the lines, as estimated in the master flat
    :param filename: string, location and file name of the file containing the polynomial fits
    """
    data = []
    temp = filename.rsplit('.',1)
    if temp[1] not in ['fit', 'fits']:
        logger('Warn: the filename does not end in fits: {0} . The proper ending will be added to save the file without problems, but this might cause issues when loading data. Please check your parameters.'.format(filename))
        filename += '.fits'
    
    # What is the maximum numbers of orders used for the polynomial fit
    porder = len(pfits[0])
    
    # Loop around each order and plot fits
    for pp, pf in enumerate(pfits):
        row = [pp + 1] + list(np.zeros(porder + 1 - len(pf)))
        row += list(pf)
        row += list([xlows[pp], xhighs[pp]])
        row += list(widths[pp])
        data.append(row)
    data = np.array(data)
    
    cols = ['Order'] + list(range(porder + 1))[::-1]
    cols += ['low_x', 'high_x', 'width_left', 'width_right', 'width_gauss']

    # convert to astropy table
    atable = Table()
    for c, col in enumerate(cols):
        atable[str(col)] = data[:, c]
    atable.write(filename, overwrite=True)
    logger('Info: master file to trace orders written: {0}'.format(filename))

def save_fits(pfits, xlows, xhighs, filename):
    """
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
    """
    data = []
    porder = len(pfits[0])                                                  #added by ronny
    cols = ['Order'] + list(range(porder + 1))[::-1]                        #added by ronny
    cols += ['low_x', 'high_x']

    # Loop around each order and plot fits
    for pp, pf in enumerate(pfits):
        row = [pp + 1] + list(np.zeros(porder + 1 - len(pf)))               #added by ronny
        row += list(pf)
        row += [xlows[pp], xhighs[pp]]
        data.append(row)
    data = np.array(data)

    # convert to astropy table
    atable = Table()
    for c, col in enumerate(cols):
        atable[str(col)] = data[:, c]
    atable.write(filename, overwrite=True)                                     #added by ronny

def save_arc_fits(wavelength_solution, wavelength_solution_arclines, filename):
    """
    save the arc solution in a file for later use
    :param wavelength_solution: list, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
                      (output of np.polyval)
                      i.e. p where:
                      p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
    :param wavelength_solution_arclines: list, same length as number of orders, each line contains the wavelength of the used reference lines. Each line must have the same length --> how???
    :param filename: string, location and file name of the file
    """
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

def read_arc_fits(filename):
    """
    reads the arc solution from a file
    :param filename: string, location and file name of the file
    fits file should look like the following:
    real_order central_px               2 ...             0 arclin1 ... arclin345
       float64    float64         float64 ...       float64  loat64 ...   float64
    ---------- ---------- --------------- ... ------------- ------- ... ---------
        -101.0     2150.0 -7.82221887e-07 ... 5.6899272e+03 5678.91 ...  6543.987
        -100.0     2150.0 -7.39532842e-07 ... 5.7469279e+03 6542.10 ...   6677.32
    ...
    where 6, 5, 4, 3, 2, 1, 0 are the polynomial powers in p
        i.e. p where:
            p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
    :return wavelength_solution: list, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    :return wavelength_solution_arclines: list, same length as number of orders, each line contains the wavelength of the used reference lines. Each line must have the same length
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
                if atable[p][order] <> 0:
                    arclines.append(atable[p][order])
            else:
                entries.append(atable[p][order])        
        wavelength_solution.append(entries)
        wavelength_solution_arclines.append(arclines)
    logger('Info: Master file for arc solution read: {0}'.format(filename))
    return np.array(wavelength_solution), wavelength_solution_arclines                # Make this into an array!

def save_multispec(data, fname, head, bitpix='-64'):
    """
    Saves a n-d array into a fits file
    :param data: n-d array of data
    :param fname: filename in which the data should be written. The filename will end in .fits
    :param head: header, which should be written to the fits file
    """
    fname = fname.replace('.npy', '')   
    fname = fname.replace('.fits', '')
    fname = fname.replace('.fit', '')
    if bitpix == '-32':
        # relative difference is less than 1E-6 of the wavelength/flux value compared to float64, only needs half the size
        hdu = fits.PrimaryHDU(np.array(data).astype(np.float32), header=head)
    elif bitpix == '-64':
        hdu = fits.PrimaryHDU(np.array(data).astype(np.float64), header=head)
    else:
        hdu = fits.PrimaryHDU(np.array(data), header=head)
        logger('Warn: The bitpix parameter ({0}) is unknown, using the data suggested one ({1})'.format(bitpix,hdu.header['BITPIX']))
        
    hdu.writeto(fname+'.fits', clobber=True)    # as .fits
    
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

def bin_im(params, im):
    """
    :param params: Dictionary with all the parameters. 'binx' and 'biny' are required
    :param im: 2d numpy array
    :return: 2d numpy array of the image
    """
    binx = params['binx']
    biny = params['biny']
    ims = im.shape
    if binx <1 or biny <1 or (binx == 1 and biny == 1):
        logger('Warning: no binning possible: {0},{1}'.format(binx,biny))
        return(im)
    nim=[]
    if binx > 1 and biny > 1:
        for i in range((ims[0]+binx-1)/binx):
            nline=[]
            iline=im[i*binx:min(ims[0],(i+1)*binx),:]
            for j in range((ims[1]+biny-1)/biny):
                nline.append(np.median(iline[:,j*biny:min(ims[1],(j+1)*biny)]))    #median binning
                #nline.append(sum(percentile_list(list((im[i*binx:(i+1)*binx,j*biny:(j+1)*biny]).flat),.2)))   #sum after 80% clipping -> more flux -> more lines find with gauss of height 50
            nim.append(nline)
    elif binx > 1:
        for i in range((ims[0]+binx-1)/binx):
            nim.append(np.median(im[i*binx:min(ims[0],(i+1)*binx),:],0))
    elif biny > 1:
        for i in range((ims[1]+biny-1)/biny):
            nim.append(np.median(im[:,i*biny:min(ims[1],(i+1)*biny)],1))
        nim = np.transpose(np.array(nim))
    logger('Info: image binned')
    return np.array(nim)

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
    :param prcentl: gives the number of how much data to remove
    :return data: sorted list without the highest and smallest values
    """
    if prcentl < 0 or prcentl > 0.5:
        print 'Warn: Percentile must be in the range [0,0.5]'
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
    if x0==0 and sigma==0 and offset==0:     #Necessary in order for the fit to work
        [amplitude, x0, sigma, offset] = amplitude
    if sigma == 0:
        return np.repeat([np.NaN],len(x))
    return amplitude * np.exp(-(x-x0)**2/(2*sigma**2))+offset

def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
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
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
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

def polynomialfit(x, y, nmin=1, nmax=3):        # Not really necessary, the result is always that nmax fits best
    """
    Select the minimum chi squared polynomial fit between nmin and nmax
    :param x: x axis array
    :param y: y axis array
    :param nmin: minimum order to fit
    :param nmax: maximum order to fit
    return ps: array of the parameters from the fit
    return nrange: number of parameters for the best fit, e.g. len(ps)
    """
    chis, ps = [], []
    nrange = range(nmin, nmax+1)
    for n in nrange:
        p = np.polyfit(x, y, n)
        ps.append(p)
        chis.append(np.sum((y-np.polyval(p, x))**2))
    argmin = np.argmin(chis)
    return ps[argmin], nrange[argmin]

def polyfit_adjust_order(xarr, yarr, p_orders, w=None):
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
    :param xx: array of floats, flattened results of a meshgrid of the x-coordinates
    :param yy: array of floats, same dimension of xx, flattened results of a meshgrid of the y-coordinates
    :params xord, yord: number of orders for the polynomial to be used in each direction, order = 0 means constant value, same conversion as np.polyfit
    :return: Arry with the dimension of xx or yy in one axis and number possible combinations in the other direction (i.e. 16 for xord=3 and yord=3, or 25 for 4,4)
    """
    
    """ # old, slow solution   
    ims = np.insert(xx.shape, 0, 1)
    print ims, xx.shape
    b = np.array([xx*0+1])
    if len(ims) > 1:                # If the input is an array of values
        b.shape = ims
    for k in tqdm(range(1, (xord+1)*(yord+1) )):             # All combinations, starting with the lowest order ones in order
        for j in range(yord+1):
            for i in range(xord+1):
                if i+j == k:
                    data = xx**i * yy**j
                    #if k==3:
                    #    print 'xx**{0} * yy**{1}'.format(i,j), data
                    if len(ims) > 1:                # If the input is an array of values
                        data.shape = ims
                        b = np.append(b, data, axis=0)
                    else:
                        b = np.append(b, data)      # If the input is a single values
                    print b.shape, data.shape
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

def polynomial_fit_2d_norm(xx, yy, zz, xord, yord, w=[]):
    """
    Fits a 2d poynormial to the data.
      Before the fit the axes are normalised in order to weight both axis the same way
    :params xx,yy,zz: all data needs to have the same entries: x-value, y-value, and z-value
    :params xord, yord: number of orders for the polynomial to be used in each direction, order = 0 means constant value, same conversion as np.polyfit (was order=1 before 20171207)
    :params w: weights of the individual points. Weights are created by copying the value the respective times, depending on the weight
    """
    ww = w
    # Get rid of the data points with weight == 0
    if len(ww) > 0:
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
        A[:,i] = A[:,i]/norm
    # Add the weights -> more entries if the weight is higher (crude solution)
    if len(ww) > 0:
        if min(ww)/max(ww) < 0.5:
            # Normalise so that the smallest weight has a value of 1
            #ww = ((ww/min(ww))**2 * 10).astype(int)       # Normalise so that the smallest weight has a value of 10, but this increaese the array so massively
            ww = ((ww/min(ww))).astype(int)
            for i in range(ww.shape[0])[::-1]:
                if ww[i] > 1:
                    A = np.concatenate( (A, A[[i],:].repeat(ww[i]-1,axis=0)) , axis=0)
                    zz = np.concatenate( (zz, zz[[i]].repeat(ww[i]-1,axis=0)) , axis=0)
            
    if np.prod(A.shape) > 1E7:
        logger('Step: Performing the fit, that might take a bit')
    coeff, res, rank, s = np.linalg.lstsq(A, zz, rcond=1E-50)        # for wavelengths fit
    #fit = polynomial_value_2d(xx, yy, xord, yord, coeff)             # calculates the fit
    #print np.array([xx,yy,zz, fit]).T                 to see the residuals between zz and the fit
    
    return coeff/norm_f

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

    poly2d_params = polynomial_fit_2d_norm(xx, yy, zz, xord, yord, w=w)       #x,y,z,order_x,order_y, w=weigths
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
    elif xarr.shape[0] <> yarr.shape[0]:
        logger('Warn: got different sized arrays for sigma clipping. This is a programming error and should not happen')
        return [], np.repeat([0], p_orders)
    goodvalues = (yarr*0 == 0)
    old_values = [0,0]
    for i in range(repeats):
        poly = polyfit_adjust_order(xarr[goodvalues], yarr[goodvalues], p_orders)
        stddiff = np.std(yarr[goodvalues] - np.polyval(poly, xarr[goodvalues]), ddof=p_orders)
        diff = (yarr - np.polyval(poly, xarr))/stddiff
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
    if len(rangediff) <> 0 and len(range2pct) <> 0:
        left = max(min(rangediff),min(range2pct))
        right = min(max(rangediff)+smoothdiff,max(range2pct))
    elif len(rangediff) <> 0:
        left = min(rangediff)
        right = max(rangediff)+smoothdiff
    elif len(range2pct) <> 0:
        left = min(range2pct)
        right = max(range2pct)
    #if right - left < 2:
    #    print 'leftmin,rightmin,linemax,linemin,range30pct,range2pct,rangediff,y[leftmin:rightmin+1].astype(int),left,right',leftmin,rightmin,linemax,linemin,range30pct,range2pct,rangediff,y[leftmin:rightmin+1].astype(int),left,right
    return int(round((right-left+.1)/2)), int(round((right+left+.1)/2)), leftmin, rightmin

def centroid_order(x,y,center,width, significance=3, bugfix=False):
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
    if (x[pos_max] - center) < width/2:
        center = x[pos_max]
    p0=[rangey,center,width,min(y)]
    bounds=((0.2*rangey/5, center-max(1.2,0.3*width), 0.1*width, min(y)-0.2*rangey), (5*rangey, center+max(1.2,0.3*width), 2*width, min(y)+0.2*rangey))
    #print 'p0,bounds',p0,bounds
    stdlin = np.std(y, ddof=1)
    significant = False
    for border_sub in [ [0,0], [0,1], [1,0], [1,1], [2,1], [1,2], [2,2] ]:    
        range_data = range(border_sub[0],len(x)-border_sub[1])
        if len(range_data) <= 6:
            break
        try:
            popt,pcov = curve_fit(oneD_gauss,x[range_data],y[range_data],p0=p0, bounds=bounds)            #a, x0, sigma, b: a*np.exp(-(x-x0)**2/(2*sigma**2))+b
        except:
            # print 'curve fit failed'
            continue
        stdfit = np.std(oneD_gauss(x[range_data],popt)-y[range_data], ddof=len(popt))           # if the average is 0, than this gives the stray of the residuals
        if stdlin/stdfit >= significance or popt[0]/stdfit >= significance:
            significant = True
            break
        elif bugfix:
            print 'Gauss not significant', p0,bounds, stdfit, border_sub
            plot_img_spec.plot_spectra(np.array([x,x]),np.array([y,oneD_gauss(x,popt)]),['data','fit'], ['savepaths'], True, [0.06,0.92,0.95,0.06, 1.0,1.01], 'No significant fit, stdlin={0}, stdGauss={1}, height={2}, needed significance={3}'.format(stdlin,stdfit,popt[0],significance))
    #print 'stdfit,stdlin', stdfit,stdlin,stdlin/stdfit < significance, popt[0]/stdfit < significance, popt
    if not significant:   # fit is not significant
        return  np.array([0,0,0,0])
    return popt

def find_center(imslice, oldcenter, x, maxFWHM, significance=3.0, bugfix=False):
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
    maxFWHM = max(3, maxFWHM)         # To have enough pixel (in order to fit the Gauss)
    oldcenter = max(0, oldcenter)
    oldcenter = min(len(imslice)-1,oldcenter)
    width, center, leftmin, rightmin = width_order(imslice,oldcenter,maxFWHM)
    if width > maxFWHM or width == 0:                                        
        #print 'bad', width, center, leftmin, rightmin
        #print 'x,oldcenter,width,center,leftmin,rightmin,imslice[leftmin:rightmin]', x, oldcenter, width, center, leftmin,rightmin, imslice[leftmin:rightmin+1]
        #print 'bad fit: ',index,value,len(deleted)
        width = 0
        centerfit = center
    else:
        while rightmin-leftmin < max(5,2*maxFWHM):             # To have enough pixel to fit the Gauss, 2 times the maximum FWHM
            leftmin = max(0, leftmin-1)
            rightmin = min(len(imslice)-1, rightmin+1)
        if len(range(leftmin,rightmin+1)) <> len(imslice[leftmin:rightmin+1]):
            print 'Error, not same length', x, oldcenter, len(imslice), leftmin, rightmin
        popt = centroid_order(range(leftmin,rightmin+1), imslice[leftmin:rightmin+1], center, width, significance=significance, bugfix=bugfix)
        #print 'good', width, center, leftmin, rightmin, popt
        centerfit = popt[1]
        #center_index = int(round(centerfit))
        #if popt[2] < maxFWHM/25. and oldcenter < 850:
        #    print 'x,oldcenter,width,center,leftmin,rightmin,popt,imslice[leftmin:rightmin]', x, oldcenter, width, center, leftmin,rightmin,popt
        #    print imslice[leftmin:rightmin+1]
        width = popt[2]
        #if popt[2] < maxFWHM/25./2:  # at least 1px FWHM, -> disabled because won't work for heavy binning
        #    #print 'x,oldcenter,width,center,leftmin,rightmin,popt,imslice[leftmin:rightmin]', x, oldcenter, width, center, leftmin,rightmin,popt, imslice[leftmin:rightmin+1]
        #    width = 0
    #if abs(oldcenter - centerfit) > 2 and oldcenter < 850 and width <> 0:
    #    print 'x,oldcenter,width,center,leftmin,rightmin,popt,imslice[leftmin:rightmin]', x, oldcenter, width, center, leftmin,rightmin,popt
    #    print imslice[leftmin:rightmin+1]
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
    noise = []
    for order in range(specs[0]):
        xarr = np.arange(specs[1])
        spec = spectra[order]
        nonnan = ~np.isnan(spec)
        poly = polyfit_adjust_order(xarr[nonnan], spec[nonnan], p_order)
        res = spec - np.polyval(poly, xarr)
        noise_order = np.repeat([np.NaN], specs[1])
        for index in xarr:
            res_sub = res[ max(0,index-semi_range) : min(specs[1],index+semi_range) + 1]
            res_sub = res_sub[~np.isnan(res_sub)]
            if len(res_sub) >= p_order+2:
                noise_order[index] = np.std(res_sub, ddof=p_order)
        noise.append(noise_order)
        
    return np.array(noise)

def estimate_width(im):
    ims = im.shape
    widths = []
    logger('Info: Estimate the width of the traces')
    for i in range(ims[0]*ims[1]/50):
        x = random.randint(int(0.05*ims[0]), int(0.95*ims[0]))
        y = random.randint(int(min(50,0.05*ims[1])), int(max(ims[1]-50,0.95*ims[1])))
        yr = range(max(0,y-30), min(ims[1],y+31))
        data = im[x,yr]
        pos_max = data.argmax() + max(0,y-30)       # Maximum in total y (from im)
        widths1 = []
        for w1 in range(2,10,1):
            for w2 in range(1,4,1):
                yy = im[x, max(0,pos_max-w1*w2):min(ims[1],pos_max+w1*w2+1)]
                if len(yy) < 5:
                    continue
                xx = range(len(yy))
                popt = centroid_order(xx, yy, len(xx)/2, w1)    # x,y,pos,width ; result: a,x0,sigma,b in a*np.exp(-(x-x0)**2/(2*sigma**2))+b
                
                """title = '{}, {}, {}, {} - {}, {}, {}, {}'.format(x,y,len(yr),pos_max, popt[0],popt[1],popt[2],popt[3])
                label_datapoins = '{}, '.format(popt[2])
                gauss = popt
                label_gauss = 'fit'
                shifts = popt[1]
                label_centroids = 'position'
                plot_gauss_data_center(np.array([xx]), np.array([yy]), [label_datapoins], np.array([gauss]), [label_gauss], np.array([shifts]), [label_centroids], filename='', title=title)"""
                if popt[2] <> 0:
                    widths1.append(popt[2])
        if len(widths1) < 1:
            continue
        widths.append(np.median(widths1))
        #print i, x,y, np.median(widths1), np.std(widths1, ddof=1)
        if len(widths) >= 50:
            break
    if len(widths) < 10:
        return 1
    width = np.median(widths)
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
    cen_px = ims[0]*binx/2
    breakvalue = np.percentile(im,40)       # The brightes pixel of one order needs to be brighter than this value, otherwise it won't be identified (lower values -> darker orders are found)
    im_orig = copy.copy(im)
    im_traces = np.zeros(ims)               # Masks the orders in order to avoid overlapping traces
    traces = []
    trace_pos = [range(ims[0])]     # fill array with index of slice, later for each order found the y values calculated from the polyfit are added -- this entry is not needed
    for dummy in tqdm(range(max( 1000, ims[0]*ims[1]/(1000/binx*2*maxFWHM*2) ))):      # 
        pos_max = np.unravel_index(im.argmax(), ims)
        if im_orig[pos_max] <= breakvalue:
            break
        for i in range(max(0,pos_max[0]-100/binx), min(ims[0],pos_max[0]+100/binx+1)):
            im[i,max(0,pos_max[1]-maxFWHM): min(ims[1],pos_max[1]+maxFWHM+1)] = breakvalue              #don't find this and sourinding values again
        if pos_max[0] < 100/binx or pos_max[0] > ims[0]-100/binx:                                   # The brightest point of an order shouldn't be close to the border
            continue                # the maximum is too close to the border
        #if np.max(im_traces[pos_max[0],max(0,pos_max[1]-maxFWHM/3):min(ims[1],pos_max[1]+maxFWHM/3+1)]) <> 0:        # old, which could result in too wide exclution area
        if np.max(im_traces[pos_max[0],max(0,pos_max[1]-2):min(ims[1],pos_max[1]+3)]) <> 0:
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
        no_center = 0
        last_trustworth_position = 0
        for i in range(pos_max[0],-1,-1):               # check positions upwards
            center, width, leftmin,rightmin = find_center(im_orig[i,:], int(round(oldcenter)), i, maxFWHM, significance=4)       # significance=4.0 tested as useful for HARPS, EXOhSPEC
            #if pos_max[1] > 50 and pos_max[1] < 70:
            #    #find_center(im_orig[i,:], int(round(oldcenter)), i, maxFWHM, significance=3.5, bugfix=True)
            #    print pos_max, i, center, width, leftmin,rightmin, oldcenter, abs(center-oldcenter), maxshift, len(positions)
            if width <> 0 and ( abs(center-oldcenter)<maxshift or (positions.shape[0]==0 and abs(center-oldcenter)<maxshift*2) ):        #first entry can be further off
                if im_traces[i,int(center)] <> 0 or im_traces[i,min(ims[1]-1,int(center+1))] <> 0:        # this order shouldn't cross another order
                    order_overlap = True
                    break
                positions = np.vstack([positions, [i, center, 0]])
                widths.append([center-leftmin,rightmin-center, width])
                oldcenter, lastadd, last_trustworth_position, no_center  = center, i, len(positions)-1, 0
            else:#if expected_positions <> []:
                no_center += 1
                if no_center >= max(3, 150/binx, 10*maxFWHM):         # stop searching if too many fits are unseccessful, as otherwise the fit might drift off
                    break
                center1, width, leftmin,rightmin = find_center(im_orig[i,:], int(round(oldcenter)), i, maxFWHM, significance=2.5, bugfix=False)
                if width <> 0 and abs(center1-oldcenter) < maxshift:
                    positions = np.vstack([positions, [i, center1, 1]])
                if lastadd-i == 5 and expected_positions <> []:        #add entries every 10 empty line in order to avoid the fit of the trace going off
                    positions = np.vstack([positions, [i, oldcenter, 1]])
                    lastadd = i
                if positions.shape[0] > 0:
                    if expected_positions <> []:                        # use the solution of the other traces to keep following this trace
                        oldcenter = positions[last_trustworth_position][1] - expected_positions[int(positions[last_trustworth_position,0])] + expected_positions[max(0,i-1)]
                    else:                                               # use the last values of this trace to keep following it
                        good_values = np.where( (positions[:,0] <= i+30) & (positions[:,2] == 0) )[0]
                        if len(good_values) > 5:
                            poly = np.polyfit(positions[good_values,0], positions[good_values,1], 1)
                            oldcenter = np.polyval(poly, i-1)
                    if abs(oldcenter - positions[last_trustworth_position,1]) > maxFWHM * 1.5:    # Stop, if the shift is in danger of touching the next order
                        #print 'too far off', oldcenter, positions[last_trustworth_position,1], no_center
                        break
        if positions.shape[0] <> 0:
            oldcenter = positions[0][1]
        else:
            oldcenter = pos_max[1]
        lastadd = pos_max[0]
        last_trustworth_position = 0
        no_center = 0
        for i in range(pos_max[0]+1,ims[0]):               # check positions downwards
            center, width, leftmin,rightmin = find_center(im_orig[i,:], int(round(oldcenter)), i, maxFWHM, significance=4)       # significance=4.0 tested as useful for HARPS, EXOhSPEC
            #if pos_max[1] > 50 and pos_max[1] < 70:
            #    #find_center(im_orig[i,:], int(round(oldcenter)), i, maxFWHM, significance=3.5, bugfix=True)
            #    print pos_max, i, center, width, leftmin,rightmin, oldcenter, abs(center-oldcenter), maxshift, len(positions), len(expected_positions), last_trustworth_position
            if width <> 0 and ( abs(center-oldcenter)<maxshift or (positions.shape[0]==0and abs(center-oldcenter)<maxshift*2) ):
                if im_traces[i,int(center)] <> 0 or im_traces[i,min(ims[1]-1,int(center+1))] <> 0 or order_overlap == True:
                    order_overlap = True
                    break
                positions = np.vstack([positions, [i, center, 0]])
                widths.append([center-leftmin,rightmin-center, width])
                oldcenter, lastadd, last_trustworth_position, no_center  = center, i, len(positions)-1, 0
            else:#if expected_positions <> []:
                no_center += 1
                if no_center >= max(3, 150/binx, 10*maxFWHM):         # stop searching if too many fits are unseccessful, as otherwise the fit might drift off
                    break
                center1, width, leftmin,rightmin = find_center(im_orig[i,:], int(round(oldcenter)), i, maxFWHM, significance=2.5, bugfix=False)
                if width <> 0 and abs(center1-oldcenter) < maxshift:
                    positions = np.vstack([positions, [i, center1, 1]])
                if i-lastadd == 5 and expected_positions <> []:
                    positions = np.vstack([positions, [i, oldcenter, 1]])
                    lastadd = i
                if positions.shape[0] > 0:
                    if expected_positions <> []:
                        oldcenter = positions[last_trustworth_position,1] - expected_positions[int(positions[last_trustworth_position,0])] + expected_positions[min(ims[0]-1,i+1)]
                    else:
                        good_values = np.where( (positions[:,0] >= i-30) & (positions[:,2] == 0) )[0]
                        if len(good_values) > 5:
                            poly = np.polyfit(positions[good_values,0], positions[good_values,1], 1)
                            oldcenter = np.polyval(poly, i+1)
                    if abs(oldcenter - positions[last_trustworth_position,1]) > maxFWHM * 1.5:    # Stop, if the shift is in danger of touching the next order
                        break
        if len(positions) < ims[0]/4:                                   # first this check and afterwards the check for order_overlap == True, as otherwise the order_overlap floats the output
            if len(positions) > ims[0]/10 and order_overlap == False:
                logger('Warn: the order around {0}, {1}, was skipped, as only {2} centroids along the order have been found'.format(pos_max[0]*binx,pos_max[1]*biny, len(positions)), show=False)
            continue
        if order_overlap == True:
            logger('Warn: the order around {0}, {1}, was skipped, as it overlapped with another order'.format(pos_max[0]*binx,pos_max[1]*biny))
            continue
        positions = np.array(sorted(positions, key=operator.itemgetter(0)))
        width = [np.mean(percentile_list(np.array(widths)[:,0],0.1)), np.mean(percentile_list(np.array(widths)[:,1],0.1)), np.mean(percentile_list(np.array(widths)[:,2],0.1))]      #average left, right, and gaussian width
        #print positions
        # Fit to the original image
        pfs, ns = polynomialfit(positions[:,0]*binx, positions[:,1]*biny, 1, params['polynom_order_apertures'])
        # Remove values which were added only to follow the curve
        positions = positions[positions[:,2]==0,0:2]
        # Calculate the values from the fitted curve for each slices
        yarr = np.polyval(pfs, positions[:,0]*binx)/biny
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
        traces.append([list(pfs), ns, min(positions[:,0])*binx, min(params['subframe'][0],(max(positions[:,0])+1)*binx-1), np.polyval(pfs, cen_px)])
        logger('Step: order found at central pixel: {0}. The trace was identified between Pixels {1} and {2}'.format(round(np.polyval(pfs, cen_px),1),  round(traces[-1][2]), round(traces[-1][3]) ))
    traces = np.array(sorted(np.array(traces), key=operator.itemgetter(4)))
    #plot_img_spec.plot_image(im+im_traces*np.max(im), ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'cleared')
    if traces.shape[0] < 1:
        logger('Error: Not enough traces found, only found {0} traces. Please check that the binned image looks as expected.'.format(traces.shape[0]))
    logger('Info: {0} orders found and traced'.format(traces.shape[0]))
    #print 'traces[:,0], traces[:,2].astype(int), traces[:,3].astype(int)', traces[:,0], traces[:,2].astype(int), traces[:,3].astype(int)
    return traces[:,0], traces[:,2].astype(int), traces[:,3].astype(int)            # polyfits, xlows, xhighs

def adjust_trace_orders(params, im, pfits, xlows, xhighs):
    """
    Re-traces the position of the orders
    :param params: Dictionary with all the parameters. 'binx', 'biny', 'subframe', and 'polynom_order_apertures' are required
    :param im: 2d numpy array
    :param pfits: list, length same as number of orders, polynomial values
                      (output of np.polyval)
                      i.e. p where:
                      p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]
    :param xlows: list, length same as number of orders, the lowest x pixel
                   (wavelength direction) used in each order
    :param xhighs: list, length same as number of orders, the highest x pixel
                   (wavelength direction) used in each order
    :return polyfits: 2d numpy array of floats, length same as number of orders, each line contains the polynomial values
    :return xlows: 1d array of integers, length same as number of orders, minimum position of the trace in dispersion direction
    :return xhighs: 1d array of integers, length same as number of orders, maximum position of the trace in dispersion direction
    :return widths: 2d array of floats, length same as number of orders, width of the trace in cross-dispersion direction. 
                    Each line contains the giving the width to the left, to the right, and the Gaussian width
    """
    binx = params['binx']
    biny = params['biny']
    ims = im.shape
    # maxFWHM = 25/biny for old lens
    # maxFWHM = 15/biny+2            # !! Determine in values of space between orders
    maxFWHM = max(1, int(round(estimate_width(im)*2.35482*1.5)))   # The average Gaussian width, transformed into a FWHM and extendet, as this is the maximum FWHM
    if len(pfits) > 1:
        maxshift = min(np.abs(pfits[1:,-1] - pfits[:-1,-1]))/5             # move at maximum by 20% of the minimum space between the orders
    else:
        maxshift = 10    #old binning in x
    trace_pos = [range(int(min(xlows)),int(max(xhighs)),binx)]      #x-values as first entry in the array
    for order, pfit in enumerate(pfits):
        trace_pos.append(np.polyval(pfit, trace_pos[0]))
        trace_pos[-1][:(xlows[order]-min(xlows))/binx] = -1e10     # do not try on bad data
        trace_pos[-1][(xhighs[order]-min(xlows)+1)/binx:] = -1e10  # do not try on bad data
    trace_pos = np.array(trace_pos)
    trace_pos[0,:] = trace_pos[0,:]/binx
    trace_pos = trace_pos[:, trace_pos[0,:] < ims[0] ]              # Due to binning, entries outside the binned images can appear in trace_pos[0,:]
    #traces = []
    polyfits, xlows, xhighs, widths, avg_shifts = [], [], [], [], []
    logger('Step: adjusting orders')
    for order in tqdm(range(trace_pos.shape[0]-1)):
        lastadd = trace_pos[0,0]
        positions, widths_o, shifts = [], [], []
        for j,i in enumerate(trace_pos[0,:].astype(int)):        #each pixel, is is the position in the array and i the real position on the CCD
            oldcenter = trace_pos[order+1,j]                        # +1 because first entry are pixels
            if oldcenter > -maxFWHM:
                center, width, leftmin,rightmin = find_center(im[i,:], int(round(oldcenter)), i, maxFWHM, significance=3.0)
            else:
                center, width, lastadd = 0,0, i
            #if width == 0 and oldcenter > -maxFWHM:
            #if width <> 0 and abs(center-oldcenter)>=3:
            #    print center,oldcenter,i, leftmin,rightmin
            if width <> 0 and abs(center-oldcenter)<maxshift:
                positions.append([i, center,0])
                widths_o.append([center-leftmin,rightmin-center, width])
                lastadd = i
                shifts.append(center-oldcenter)
            elif oldcenter > -maxFWHM and i-lastadd == 10:
                positions.append([i, oldcenter,1])
                lastadd = i
        #!! maybe: try to follow the order in areas, which were not defined before
        #print len(positions), 
        if widths_o == [] or len(shifts) < ims[0]/4:
            #print 'len(positions),len(shifts)', len(positions),len(shifts)
            continue
        positions = np.array(positions)
        widths_o = np.array(widths_o)
        width = [np.mean(percentile_list(widths_o[:,0],0.1)), np.mean(percentile_list(widths_o[:,1],0.1)), np.mean(percentile_list(widths_o[:,2],0.1))]      #average left, right, and gaussian width
        # Fit to the original image
        """
        pfs, ns = polynomialfit(positions[:,0]*binx, positions[:,1]*biny, 1, params['polynom_order_apertures'])
        # Remove values which were added only to follow the curve
        positions = positions[positions[:,2]==0,0:2]
        traces.append([ list(pfs) ,ns, max(0,min(positions[:,0])*binx-5), min(params['subframe'][0],(max(positions[:,0])+1)*binx+5), list(width) ])
        #print len(positions), min(positions[:,0]),(max(positions[:,0])+1)*binx-1
        logger('Step: order {0} adjusted: avg shift = {1}, min shift = {2}, max shift = {3}, avg width: gauss = {6}, left = {4}, right = {5}'.format(order, round(np.mean(shifts),2),round(min(shifts),2),round(max(shifts),2), round(width[0],2),round(width[1],2),round(width[2],2) ))"""
        pfs = np.polyfit(positions[:,0]*binx, positions[:,1]*biny, params['polynom_order_apertures'])
        polyfits.append(pfs)
        xlows.append(max(0,min(positions[:,0])*binx-5))
        xhighs.append(min(params['subframe'][0],(max(positions[:,0])+1)*binx+5))
        widths.append(width)
        avg_shifts.append(np.mean(shifts))
    """traces = np.array(traces)"""
    logger('Info: traces of the {0} apertures adjusted. The average shift of the individual apertures was between {1} and {2} pixel between the searching of the traces and this solution'.format(len(polyfits), round(min(avg_shifts),1), round(max(avg_shifts),1) ))
    if len(polyfits) == 0:
        logger('Error: no traces of orders found. Please check the binned image: Is the orientation and the binning right? Are the orders covering at least half of the CCD (in dispersion correction)')
    """polyfits, xlows, xhighs, widths = np.array(traces[:,0]), np.array(traces[:,2].astype(int)), np.array(traces[:,3].astype(int)), np.array(traces[:,4])"""
    return np.array(polyfits), np.array(xlows).astype(int), np.array(xhighs).astype(int), np.array(widths)

def extract_orders(params, image, pfits, xlows, xhighs, widths, w_mult, offset = 0, var='fast'):
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
    :param widths: 2d list, length same as number of orders, each entry contains left border, right border, and Gaussian width of the lines, as estimated in the master flat
    :param w_mult: w_mult * widths (Gaussian) defines the range on either side of the trace to be extracted
    :param offset: allows to shift the spectrum
    :param var: can be 'fast' or 'prec'.
        'fast': extraction of the fraction of the border pixel
        'prec': fit a polynomial around the border of the extrion width to extract with higher precission
    :return spec: list,  length same as number of orders, each containing
                  a 1D numpy array with the spectrum of each order in

    """
    logger('Step: extracting spectra', show=False)
    badpx_mask_name = 'badpx_mask'
    if badpx_mask_name not in calimages:
        calimages[badpx_mask_name] = read_badpx_mask(params)
    badpx_mask = calimages[badpx_mask_name]
    spec, good_px_mask = [], []
    xarr = range(0, image.shape[0])
    maxpx = image.shape[1]-1
    for pp, pf in tqdm(enumerate(pfits)):
        #print('\t Order {0}'.format(pp + 1))
        # find nearest y pixel
        yarr = np.polyval(pf, xarr)+offset
        xarr1 = np.arange(xlows[pp], xhighs[pp])
        ospecs = np.repeat([np.NaN], image.shape[0])
        ogood_px_mask = np.repeat([0.0], image.shape[0])            # Array needs to be float in order to save 0.1, 0.2, ... correctly
        w = widths[pp]
        # getting the lower and upper bounds of the trace, doing everything in np, to speed up things
        lowers = yarr - w[2]*w_mult
        uppers = yarr + w[2]*w_mult
        if var == 'prec':         # new solution, better than the one below, tested with many values in excel
            lowersf = (np.ceil(lowers+.5)).astype(int)       #Full lowest pixel of the order
            uppersf = (np.floor(uppers-.5)).astype(int)      #Full highest pixel of the order
            lowers[lowers<-0.5] = -0.5
            uppers[uppers>maxpx+0.5] = maxpx+0.5
            lowersf[lowersf<0] = 0
            uppersf[uppersf>maxpx] = maxpx
            order_in_frame = ((lowersf[xarr1] <= maxpx) & (uppersf[xarr1] >= 0))     # Due to order shifting the order might be shifted to the ouside of the images, ignore these values
            # Loop around x-px, this is necessary, as the number of full pixels between lowersf and uppersf varies 
            for xrow in xarr1[order_in_frame]:
                ssum = np.sum(image[xrow, lowersf[xrow]:uppersf[xrow]+1])
                good_px = 1
                if min(badpx_mask[xrow, lowersf[xrow]:uppersf[xrow]+1]) == 0:
                    good_px = 0.2
                elif max(image[xrow, lowersf[xrow]:uppersf[xrow]+1]) >= params['max_good_value']:
                    good_px = 0.1
                # Get fractions of pixel from a polynom, if they are not outside the frame borders (spline actes weired between points, when there is an outlier)
                if w[2]*w_mult < 0.5:             # There will be no full pixel, which messes up everything
                    fracpixparam = [[yarr[xrow],'b']]
                else:
                    fracpixparam = [[lowers[xrow], 'l'] , [uppers[xrow],'u']]
                for [borderpos, pos] in fracpixparam:
                    x = np.arange(max(0,np.floor(borderpos-1)), min(maxpx,np.ceil(borderpos+1))+1).astype(int)
                    y = image[xrow, x]
                    if min(badpx_mask[xrow, x]) == 0:
                        good_px = 0.2
                    elif max(y) > params['max_good_value']:
                        good_px = 0.1
                    weight = 1/(np.abs(x-borderpos)+0.1)**2                   # weight the values so that the data around the fraction of the pixel is used most
                    p = polyfit_adjust_order(x, y, max(1,len(x)-3), w=weight)
                    poly = np.poly1d(p)
                    polyint = poly.integ()
                    if pos == 'l':
                        if lowersf[xrow] <> 0:
                            ssum += polyint(lowersf[xrow]-0.5) - polyint(lowers[xrow])
                    elif pos == 'u':
                        if uppersf[xrow] <> maxpx:
                            ssum += polyint(uppers[xrow]) - polyint(uppersf[xrow]+0.5)
                    elif pos == 'b':
                        ssum += polyint(uppers[xrow]) - polyint(lowers[xrow])
                    else:
                        print 'Programming error around line 1480'
                ospecs[xrow] = ssum
                ogood_px_mask[xrow] = good_px
        else:           # old solution
            lowersf = (np.round(lowers+.5)).astype(int)
            uppersf = (np.round(uppers-.5)).astype(int)
            lowersf[lowersf<0] = 0
            uppersf[uppersf>maxpx] = maxpx
            order_in_frame = ((lowersf[xarr1] <= maxpx) & (uppersf[xarr1] >= 0))     # Due to order shifting the order might be shifted to the ouside of the images, ignore these values
            lowersr = (np.round(lowers)).astype(int)
            uppersr = (np.round(uppers)).astype(int)
            # Loop around x-px, this is necessary, as the number of full pixels between lowersf and uppersf varies 
            for xrow in xarr1[order_in_frame]:
                ssum = np.sum(image[xrow, lowersf[xrow]:uppersf[xrow]+1])
                good_px = 1
                if min(badpx_mask[xrow, lowersf[xrow]:uppersf[xrow]+1]) == 0:
                    good_px = 0.2
                elif max(image[xrow, lowersf[xrow]:uppersf[xrow]+1]) >= params['max_good_value']:
                    good_px = 0.1
                # Get fractions of pixel, if they are not outside the frame borders
                if lowersr[xrow] >= 0:
                    ssum += image[xrow, lowersr[xrow]]*(0.5-lowers[xrow]%1)
                    if badpx_mask[xrow, lowersr[xrow]] == 0:
                        good_px = 0.2
                    elif image[xrow, lowersr[xrow]] > params['max_good_value']:
                        good_px = 0.1
                if uppersr[xrow] <= maxpx:
                    ssum += image[xrow, uppersr[xrow]]*(uppers[xrow]%1-.5)
                    if badpx_mask[xrow, uppersr[xrow]] == 0:
                        good_px = 0.2
                    elif  image[xrow, uppersr[xrow]] >= params['max_good_value']:
                        good_px = 0.1
                ospecs[xrow] = ssum
                ogood_px_mask[xrow] = good_px

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
            x = np.arange(lowerw[xrow], upperw[xrow]+1).astype(int)
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
    steps = 10
    shifts, twidths, positions = [], [], []
    problem_order = []
    shift_map = np.zeros((im.shape[0],len(sci_tr_poly)))
    logger('Step: searching for shifts of the orders') 
    for i in tqdm(range(len(sci_tr_poly))):       # For each order
        xarr = range(int(xlows[i]),int(xhighs[i]),steps)
        yarr = np.polyval(sci_tr_poly[i], xarr)+in_shift          #apply input shift
        widths = []
        for j,oldcenter in enumerate(yarr):        #each pixel
            center, width, leftmin,rightmin = find_center(im[xarr[j],:], int(round(oldcenter)), xarr[j], 2*params['maxshift'], significance=3.0)
            #if width <> 0:
            #    print 'i,j,center, width, leftmin,rightmin',i,j,center, width, leftmin,rightmin
            if width <> 0 and abs(center-oldcenter) < params['maxshift']:
                widths.append([center-leftmin,rightmin-center, width])
                shifts.append(center-oldcenter)
                shift_map[max(0,xarr[j]-steps/2):min(im.shape[0],xarr[j]+steps/2)+1,i] = center-oldcenter
        if widths == []:
            twidths.append(oldwidths[i])
            problem_order.append(i)
            continue
        widths = np.array(widths)
        width = [np.mean(percentile_list(widths[:,0],0.1)), np.mean(percentile_list(widths[:,1],0.1)), np.mean(percentile_list(widths[:,2],0.1))]      #average left, right, and gaussian width
        twidths.append(width)
    if shifts <> []:
        shift = np.mean(percentile_list(np.array(shifts),0.1))
        shift_error = np.std(percentile_list(np.array(shifts),0.1), ddof=1)
    else:
        shift = 0
    printarrayformat = ['%1.1i', '%4.2f', '%4.2f']
    printarray = np.array([ range(len(sci_tr_poly)), np.array(oldwidths)[:,2], np.array(twidths)[:,2] ]).T
    problem_text = ''
    if problem_order <> []:
        problem_text = ' The shift in the following orders could not be measured: {0}'.format(problem_order)
        if len(problem_order) > 0.5*len(sci_tr_poly):
            problem_text += '. No shift could be measured for {0} out of {1} orders.'.format(len(problem_order), len(sci_tr_poly))
            shift_error = -1
    logger('Info: The shift of the traces was determinded to {0} +- {2}, including the inputshift of {1}.{3}'.format(round(shift+in_shift,2), in_shift, round(shift_error,2), problem_text ))
    logger('Info: The gaussian width of the orders changed: order\told\tnew', printarrayformat=printarrayformat, printarray=printarray)
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
    for pp,pfit in enumerate(pfits):
        xarr = range(max(0,xlows[pp]-10),min(ims[0],xhighs[pp]+10))
        yarr = np.polyval(pfit, xarr)
        w = widths[pp]
        for i,xrow in enumerate(xarr):
            lower = max(0,int(yarr[i]-w[2]*w_mult))
            upper = min(ims[1],int(yarr[i]+w[2]*w_mult+1))
            if upper < 0:               # Otherwise it goes from e.g. 0 to -19, which will be nearly the whole image
                upper = 0
            im[xrow, lower:upper] = 0
    return im

def find_shift_images(im, im_ref, sci_tr_poly, cal_tr_poly):
    """
    Finds the shift between two images by cross corellating both. The images are substracted and the total flux [sum(abs())] is calculated. At the best position, a minimum will be reached
    :param im: 2d array with the image, for which the shift should be calculated
    :param im_ref: 2d array with the reference image
    :return shift: shift of the image (float)
    """
    # Find the maximum shift by using the space between science and calibration orders
    shifts1 = np.abs(sci_tr_poly[:,-1] - cal_tr_poly[:,-1])
    shifts2 = np.abs(sci_tr_poly[1:,-1] - cal_tr_poly[:-1,-1])      # calibration fiber could be left or right of science fiber
    shifts3 = np.abs(sci_tr_poly[:-1,-1] - cal_tr_poly[1,-1])       # calibration fiber could be left or right of science fiber
    shift = min( min(shifts1), min(shifts2), min(shifts3) )
    
    range_shifts = [-int(shift/2),int(shift/2)]
    ims1 = im.shape[1]
    shifts = range(min(range_shifts),max(range_shifts)+1)
    fluxdiff = []
    for shift in shifts:
        fluxes = []
        diffim = im[:,max(0,shift):ims1-max(0,-shift)] - im_ref[:,max(0,-shift):ims1-max(0,shift)]      # Tested that the minus signs are at the right place
        ims2 = diffim.shape[1]
        temp = max(np.abs(range_shifts).astype(int))-abs(shift)
        for i in range(temp+1):         # For a small shift, the resulting image is bigger than for a big shift -> avoid impact total flux
            #print shift, i, diffim.shape, diffim[:,i:ims2-(temp-i)].shape, np.sum(np.abs(diffim[:,i:ims2-(temp-i)]))
            fluxes.append(-np.sum(np.abs(diffim[:,i:ims2-(temp-i)])))                   # using a negative flux, so the minimum becomes a maximum, which only can be handled by the script in order to fit a gaussian
        fluxdiff.append(np.mean(fluxes))
        #print 'shift, fluxdiff[-1], fluxes', shift, fluxdiff[-1], fluxes
    #print shifts,fluxdiff,shifts[np.argmax(fluxdiff)],max(shifts)-min(shifts)
    popt = centroid_order(shifts,fluxdiff,shifts[np.argmax(fluxdiff)],max(shifts)-min(shifts))    #a,x0,sigma,b in a*np.exp(-(x-x0)**2/(2*sigma**2))+b
    shift = round(popt[1],2)
    logger('Info: The shift between this frame and the reference frame is {0} px. The gaussian width is {1} px.'.format(shift, round(popt[2],2) ))
    if popt[2] > 0.3 * (max(range_shifts) - min(range_shifts)):
        logger('Warn: The center is not very precise or the measured shift is outside of the allowed range ({0} to {1})  -> using a shift of 0.0'.format(min(range_shifts),max(range_shifts)))
        shift = 0.0
    return shift

def arc_shift(params, im_arc, pfits, xlows, xhighs, widths):
    """
    Determines the difference between the science orders and arc orders in cross-dispersion direction
    :param params: Dictionary with all the paramerters. arcshift_range (is made positive when reading the configuration), arcshift_side, and arcextraction_width_multiplier will be used
    :param im_arc: 2d image in which the arc orders are iluminated
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
    sigma = 3
    orders = np.arange(len(pfits))
    cen_pos, cen_pos_diff = [], []
    arcshift_range = params['arcshift_range']
    if len(pfits) >= 10:
        for pf in pfits:
            cen_pos.append(np.polyval(pf, im_arc.shape[0]/2))       # The central pixel of each order
        cen_pos = np.array(cen_pos)
        cen_pos_diff = np.abs(cen_pos[1:] - cen_pos[:-1])           # The difference between the central pixels
        arcshift_range = [int(np.median(widths[:,2])*2), int(np.max([np.median(cen_pos_diff[:5]), np.median(cen_pos_diff[-5:])]) - np.median(widths[:,2])*2 + 1) ]  # Test everything between the 3 times Gauss width and next Order-3*width
        cen_pos_diff = np.insert(cen_pos_diff, 0, cen_pos_diff[0])                      # create the same length array
    arcshifts = np.arange(min(params['arcshift_side']*np.array(arcshift_range)),max(params['arcshift_side']*np.array(arcshift_range))+1)
    fluxes = []
    print 'Searching for the shift of the arc traces, compared to the science traces:'
    for shift in tqdm(arcshifts):
        # Shift the solution -> give parameter offset; and extract the arc_spectrum with the offset
        arc_spec, good_px_mask = extract_orders(params, im_arc, pfits, xlows, xhighs, widths, params['arcextraction_width_multiplier'], shift)
        flux = []
        for order in orders:
            # Search for significant lines
            xarr = np.arange(arc_spec.shape[1])
            yarr = arc_spec[order,:]
            notnans = ~np.isnan(yarr)
            yarr1 = yarr[notnans]
            xarr1 = xarr[notnans]
            if len(xarr1) < 5:
                flux.append(-100)
                continue
            notline_pos, pfit = sigma_clip(xarr1, yarr1, xarr1.shape[0]/300, 10, 2.5, repeats=10)  #orders, sigma low, sigma high
            line_pos = ~notline_pos
            for i in range(1,len(line_pos)-1):
                if (line_pos[i-1:i+2]==[False, True, False]).all():     #exclude lines over only one px
                    line_pos[i] = False
            # Get the flux of the order
            flux.append(sum(yarr1[line_pos]) )
        fluxes.append(flux)
    fluxes = np.array(fluxes)
    # Find the center in each order: where is the flux the highest
    gauss, goodvalues = [], []
    label_datapoins, label_gauss, label_centroids = [],[],[]
    for order in orders:
        if cen_pos_diff <> []:
            goodpos = ( (np.abs(arcshifts) > widths[order,2]*2) & (np.abs(arcshifts) < cen_pos_diff[order]-widths[order,2]*2) )
            x, y = arcshifts[goodpos], fluxes[goodpos,order]
            if len(x) < 5:                      # not enough values, use standard settings
                goodpos = ( (np.abs(arcshifts) > min(arcshift_range)) & (np.abs(arcshifts) < max(arcshift_range)) )
                x, y = arcshifts[goodpos], fluxes[goodpos,order]
        else:
            x, y = arcshifts, fluxes[:,order]
        popts = []
        for centeroffset in range(0,int(len(x)/2-widths[order,2]),3):                # Test for different positions of the line
            for sign in [-1,+1]:
                popt = centroid_order(x, y, np.mean(x)+sign*centeroffset, min(15,max(x)-min(x)) )    #input: x,y,center, width; result: a,x0,sigma,b in a*np.exp(-(x-x0)**2/(2*sigma**2))+b
                if np.all(popt <> [0,0,0,0]):
                    diff = np.sum(np.abs(oneD_gauss(x,popt) - y))
                    popts.append([diff, list(popt)])
                    break
        if len(popts) <> 0:
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
        logger('Warn: A shift of the arc orders could not be measured. The average value of the input parameters "arcshift_range" are used instead')
        shifts = np.repeat([round(np.mean(params['arcshift_side']*np.array(params['arcshift_range'])),2)], orders.shape[0])
        poly, diff, min_gauss, max_gauss = 0, [0], 0, 0
    else:
        # Fit a polynomial to the center
        goodval, poly = sigma_clip(orders[goodvalues], gauss[goodvalues,1], 2, sigma, sigma, repeats=5)
        shifts = np.round(np.polyval(poly, orders),2)
        diff = shifts[goodvalues]-gauss[goodvalues,1]
        min_gauss, max_gauss = round(min(gauss[goodvalues,2]),2), round(max(gauss[goodvalues,2]),2)
    # Log the results
    arcshifts = np.repeat([arcshifts],len(orders),axis=0)
    title = 'Determining the shift of the arc orders'
    plot_gauss_data_center(arcshifts, fluxes.T, label_datapoins, gauss, label_gauss, shifts, label_centroids, params['logging_find_arc_orders'], title=title)
    logger('Info: The shift between flat traces and arc lines is between {0} and {1} px. The standard difference of the measured values to the 2nd order polynomial is {4} px. The parameters of the polynomial are: {5}. The gaussian width of the arc lines in spacial direction is between {2} and {3} px'.format(min(shifts), max(shifts), min_gauss, max_gauss, round(np.std(diff, ddof=orders),2), np.round(poly,5) ))
    return shifts
    
def identify_lines(params, im, im_short=None, im_badpx=None, new_format=False):         # Old format is only needed for ident_arc.py, this routine needs revision
    """
    Identifies the lines in a spectrum by searching for the significant outliers in a polynomial fit to the data and subsequent fitting of Gaussian profiles to this positions
    :param im: 2d array with the extracted spectra
    :param im_short: 2d array with the extracted spectra in which the saturated lines of {im} should be identified
    :param im_badpx: The bad-pixel-mask for the extracted spectrum {im}. Used to identify the lines with saturated pixels
    :param new_format: For compatibility, now always new_format=True should be used
    :return lines: (new format): 2d array with one line for each identified line, sorted by order and amplitude of the line. For each line the following informaiton is given:
                    order, pixel, width of the line, and height of the line
    :return lines: (old format) list with one entry for each order in which lines have been found.
                    Each entry consists of the number of the order and an array which contains all lines. 
                            Each entry in the array consists of the sorted pixel, width of the line, and height of the line
    """
    ims = im.shape
    lines = []
    lines1 = []
    FWHM = 4
    print 'Identify lines in the arc'
    for order in tqdm(range(ims[0])):
        xarr = np.arange(ims[1])
        yarr = im[order,:]
        notnans = ~np.isnan(yarr)
        yarr1 = yarr[notnans]
        xarr1 = xarr[notnans]
        notline_pos, pfit = sigma_clip(xarr1, yarr1, 12, 10, 2.2, repeats=20)  #orders, sigma low, sigma high ; tested that 2.5 is too high for some lines (UNe), If lines are not identified correctly, the problem isn't here, probably
        line_pos = ~notline_pos
        for i in range(1,len(line_pos)-1):
            if (line_pos[i-1:i+2]==[False, True, False]).all():     #exclude lines over only one px
                line_pos[i] = False
        #print 'xarr1[line_pos], pfit',xarr1[line_pos], pfit
        lines_order = []
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
                #print order, pos_max, np.sum(im_badpx[order,range_arr]==0), np.sum(im_badpx[order,range_arr]==0.1), np.sum(im_badpx[order,range_arr]==0.2), np.sum(im_badpx[order,range_arr]==1), np.sum(im_badpx[order,range_arr]<>-1)
                if 0.1 in im_badpx[order,range_arr]:
                        y_data = im_short[order,range_arr]
                        if np.sum(already_found[range_arr]) == 0:           # Only print once in the area
                            print 'Used the short exposure time to find the centroid of the line in order {1} @ px {0}'.format(pos_real, order)
            # fit the gauss, pos_real is the expected center
            popt = centroid_order(xarr[range_arr],y_data,pos_real,FWHM*2, significance=2)    #a,x0,sigma,b in a*np.exp(-(x-x0)**2/(2*sigma**2))+b
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
            lines_order.append([popt[1], popt[2], popt[0] ])
            lines1.append([order, popt[1], popt[2], popt[0] ])
            #centerfit, width, leftmin,rightmin = find_center(yarr, pos_real, order, 10)
            #if width == 0:
            #    continue
            #if already_found[int(round(centerfit))] == 1:
            #    continue
            #already_found[leftmin:rightmin+1] = 1
            #lines_order.append([centerfit, width])
            #print order,centerfit, width, leftmin,rightmin
        if lines_order <> []:
            lines_order = np.array(sorted(lines_order, key=operator.itemgetter(0) ))
            lines.append([order, lines_order])
    lines1 = np.array(lines1)
    # Remove the lines which are too wide/narrow
    good_values, pfit = sigma_clip(lines1[:,2]*0, lines1[:,2], 0, 3, 3)  #orders, sigma low, sigma high
    lines1 = lines1[good_values, :]
    
    logger('Info: Identified {0} lines in the arc spectrum. These lines are stored in file {1}'.format(len(lines1), params['logging_found_arc_lines']))
    printarrayformat = ['%1.1i', '%3.2f', '%3.2f', '%3.1f']
    logger('order\tpixel\twidth\theight of the line', show=False, printarrayformat=printarrayformat, printarray=lines1, logfile=params['logging_found_arc_lines'])
    if new_format:
        return lines1
    else:
        return lines

def read_reference_catalog(filename, wavelength_muliplier, arc_lines):
    """
    Reads the reference catalogue from a file and extracts the lines which should be used
    :param filename: text, filename of the catalogue file. The file needs to consist of 3 columns: wavelength of the line, strength/heigt of the line [can be empty], and name of the line
    :param wavelength_muliplier: float, if the resolution is not given in Angstrom, than it can be converted into Angstrom with this factor
    :param arc_lines: list of text, subsection of the arc lines, which should be extracted
    :return reference_catalog: 2d array with one entry for each line. Each entry contains the wavelength, the intensity of the line, and the index in the catalogue
    :return reference_names: list with same length as reference_catalog, name of each line
    """
    for i in range(len(arc_lines)):
        arc_lines[i] = arc_lines[i].replace(' ','')
    reference_catalog = []
    reference_names = []
    file = open(filename, 'r')
    for line in file:
        line = line[:-1].split('\t')
        if len(line) < 3:
            continue
        try:
            line[0] = float(line[0])*wavelength_muliplier    #wavelength in Angstrom
        except:
            continue
        for i in range(len(line[1])+1)[::-1]:
            try:
                line[1] = float(line[1][:i])
                break
            except:
                continue
        if i == 0:
            line[1] = 1
        if line[2].replace(' ','') not in arc_lines:
            continue
        reference_names.append(line[2])
        line = line[0:2]
        line.append(len(reference_catalog))
        reference_catalog.append(line)   # wavelength in Angstrom, relative intensity of the line, index of line in reference_names
    file.close()
    if reference_catalog == []:
        logger('Error: no reference lines found in {0} for the requested lines {1}'.format(filename, arc_lines))
    reference_catalog = np.array(reference_catalog)
    arcs = reference_catalog.shape
    if arcs[0] > 100000:
        breakvalue = np.percentile(reference_catalog[:,1],min(90,arcs[0]/120))
        keep = np.logical_or(reference_catalog[:,1] >= breakvalue , reference_catalog[:,1] == 1)
        for i in range(arcs[0])[::-1]:
            if keep[i] == False:
                del reference_names[i]
        reference_catalog = reference_catalog[keep,:]
        logger('Info The faintest {0} of {1} entries in the arc reference file {2} will not be used '.format(arcs[0]-reference_catalog.shape[0], arcs[0], filename ))
    return reference_catalog, reference_names

def shift_wavelength_solution(aspectra, wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, xlows, xhighs):
    """
    Determines the pixelshift between the current arc lines and the wavelength solution
    :param aspectra: spectrum of the reference orders
    :param wavelength_solution: list, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    :param wavelength_solution_arclines: list of floats, same length as number of orders, each line contains the wavelength of the used reference lines.
    :param reference_catalog: 2d array of floats with one entry for each line. Each entry contains the wavelength and the intensity of the line
    :param reference_names: list of strings with same length as reference_catalog, name of each line
    :param xlows: list of floats, length same as number of orders, the lowest x pixel (wavelength direction) used in each order
    :param xhighs: list of floats, length same as number of orders, the highest x pixel (wavelength direction) used in each order
    :return wavelength_solution: list, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    """
    #logger('Step: Finding the wavelength shift for this exposure')
    FWHM = 6
    in_shift = 0       #0 
    ratio_lines_identified = 0.15       # if more than ratio_lines_identified of the checked_arc_lines has been identified, then a sigma_clipping will be applied. If less than this number of lines remain after sigma clipping, then the assumption is, that the calibration fiber wasn't used and therefore no wavelength shift is applied
    # In each order get the approx pixel of each reference line, fit a gaussian against the position, calculate the wavelength for the gaussian center, compare the wavelength of the line center with the reference line wavelength
    aspectra = np.array(aspectra)
    ass = aspectra.shape
    shifts = []
    checked_arc_lines = 0
    for order_index in range(wavelength_solution.shape[0]):
        xarr = np.arange(ass[1])
        warr = np.polyval(wavelength_solution[order_index,2:], xarr-wavelength_solution[order_index,1])
        for arcline in wavelength_solution_arclines[order_index]:
            ref_line_index = np.argmin(np.abs( reference_catalog[:,0] - arcline ))
            if abs( reference_catalog[ref_line_index,0] - arcline ) > 0.00001:      # check that it is in the reference catalog
                continue
            diff = abs(warr-reference_catalog[ref_line_index,0])
            if min(diff) <= wavelength_solution[order_index,-2]*2:    #distance should be less than 2 px
                checked_arc_lines += 1
                pos = xarr[np.argmin(diff)]     #position of the line in px
                range_arr = range(max(0,pos-FWHM*2),min(ass[1],pos+FWHM*2+1))
                popt = centroid_order(xarr[range_arr],aspectra[order_index,range_arr],pos+in_shift,FWHM*2)    #a,x0,sigma,b in a*np.exp(-(x-x0)**2/(2*sigma**2))+b
                if popt[1] == 0 or popt[2] > FWHM*1.5:
                    #print 'not used',order_index, pos,popt
                    #plot_img_spec.plot_points([xarr[range_arr]],[aspectra[order_index,range_arr]],[str(pos)],'path',show=True, x_title='Pixel', y_title='Flux')
                    continue
                # Index of order, Index of reference line, wavelength of reference line, px of arc line, wavelength of arc line, height of arc line, sigma of gauss
                shifts.append([ order_index, ref_line_index, reference_catalog[ref_line_index,0], popt[1], np.polyval(wavelength_solution[order_index,2:], popt[1]-wavelength_solution[order_index,1]), popt[0], popt[2] ])
    shifts = np.array(shifts)
    #print shifts[:,4]-shifts[:,2]
    if len(shifts) >= ratio_lines_identified * checked_arc_lines and shifts <> []:    
        good_values, pfit = sigma_clip(shifts[:,3]*0, shifts[:,4]-shifts[:,2], 0, 3, 3, repeats=20)  #orders, sigma low, sigma high
        shifts = shifts[good_values,:]
    #print shifts[:,4]-shifts[:,2]
    shift_avg, shift_std, width_avg, width_std = 0,0,0,0
    if len(shifts) >= ratio_lines_identified * checked_arc_lines and shifts <> []:
        shift_avg, shift_std = np.mean(shifts[:,4]-shifts[:,2]), np.std(shifts[:,4]-shifts[:,2], ddof=1)
        width_avg, width_std = np.mean(shifts[:,6]), np.std(shifts[:,6], ddof=1)
    logger('Info: The wavelength shift is {0} +- {1} Angstrom. {2} reference lines have been used, {7} reference lines have been tested. The arc lines have a Gaussian width of {3} +- {4}, which corresponds to a FWHM of {5} +- {6} px'.format(round(shift_avg,4), round(shift_std,4), shifts.shape[0], round(width_avg,2), round(width_std,2), round(width_avg*2.35482,2), round(width_std*2.35482,2), checked_arc_lines ))
    if len(shifts) >= ratio_lines_identified * checked_arc_lines and shifts <> []:
        statistics_arc_reference_lines(shifts, [0,1,-1,2], reference_names, wavelength_solution, xlows, xhighs)
    # correction in the other side of the shift
    wavelength_solution_new = []
    for wls in wavelength_solution:
        #wls[-1] -= shift_avg        #This cause the wavelength_solution to change globally, even if copy.copy(wavelength_solution) is used in the for loop or in the call of the procedure
        wls_neu = wls[0:-1]
        wls_neu = np.append(wls_neu, wls[-1]-shift_avg)
        wavelength_solution_new.append(wls_neu)
    
    return np.array(wavelength_solution_new)

def create_wavelengths_from_solution(wavelength_solution, spectra):
    """
    Converts the wavelength solution into a 2d array with the wavelength for each pixel and order
    :param wavelength_solution: list, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    :param spectra: 2d array of floats, spectrum
    :return wavelengths: 2d array of floats, same dimensions as spectra, contains the wavelength for each pixel
    """
    wavelengths = []
    for wls in wavelength_solution:
        #wls[-1] -= shift_avg        #This cause the wavelength_solution to change globally, even if copy.copy(wavelength_solution) is used in the for loop or in the call of the procedure
        wls_neu = wls[2:-1]
        wls_neu = np.append(wls_neu, wls[-1])
        xarr = range(np.array(spectra).shape[1])-wls[1]
        wavelengths.append(np.polyval(wls_neu, xarr))
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
    if os.path.isfile(params['master_order_filename']) == True:
        logger('Info: Master order file already exists: {0}'.format(params['master_order_filename']))
        polyfits, xlows, xhighs, widths = read_fits_width(params['master_order_filename'])
    else:
        # Find the traces in a binned file
        if os.path.isfile(params['logging_traces_binned']) == True:
            logger('Info: Using exiting order file: {0}'.format(params['logging_traces_binned']))
            polyfits, xlows, xhighs, dummy = read_fits_width(params['logging_traces_binned'])
        else:
            # make flat smaller: combine px in spatial direction
            params['binx'], params['biny'] = params['bin_search_apertures']
            sim_sflat = bin_im(params, im_sflat)
            save_im_fits(params, sim_sflat, im_sflat_head, params['logging_sflat_binned'])
        
            # Search for orders in the small image
            polyfits, xlows, xhighs = find_trace_orders(params, sim_sflat)
            polyfits = unify_pord(polyfits)
            # save parameters of the polynoms into a fitsfile (from Neil)
            save_fits(polyfits, xlows, xhighs, params['logging_traces_binned'])
            plot_traces_over_image(im_sflat, params['logging_orders_binned'], polyfits, xlows, xhighs)
        
        # retrace orders in the original image to finetune the orders
        params['binx'], params['biny'] = params['bin_adjust_apertures']
        sim_sflat = bin_im(params, im_sflat)        # Not saved
        polyfits, xlows, xhighs, widths = adjust_trace_orders(params, sim_sflat, polyfits, xlows, xhighs)
        if params['GUI']:
            logger('Step: Allowing user to remove orders')
            fmask = run_remove_orders(np.log10(im_sflat), polyfits, xlows, xhighs, userinput=params['GUI'])
            polyfits, xlows, xhighs, widths = polyfits[fmask], xlows[fmask], xhighs[fmask], widths[fmask]
        # save parameters of the polynoms into a fitsfile (from Neil, width added)
        save_fits_width(polyfits, xlows, xhighs, widths, params['master_order_filename'])
        plot_traces_over_image(im_sflat, params['logging_orders'], polyfits, xlows, xhighs, widths)
        data = np.insert(widths, 0, range(len(polyfits)), axis=1)
        positio = []
        for order, pfit in enumerate(polyfits):
            positio.append(np.polyval(pfit, im_sflat.shape[0]/2))
        data = np.append(data, np.array(positio)[:,None],axis=1)
        data = np.append(data, np.array(xlows)[:,None],axis=1)
        data = np.append(data, np.array(xhighs)[:,None],axis=1)
        printarrayformat = ['%1.1i','%3.1f', '%3.1f', '%4.2f\t','%1.1i','%1.1i','%1.1i']
        logger('\t\torder\tleft\tright\tgausswidth\tpositio\tmin_tr\tmax_tr\t(positio at the center of the image)',printarrayformat=printarrayformat, printarray=data)
        
    return polyfits, xlows, xhighs, widths

def extraction_steps(params, im, im_name, im_head, polyfits, xlows, xhighs, widths, apolyfits, axlows, axhighs, awidths, wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, flat_spec_norm, im_flat):
    """
    Extracts the spectra and stores it in a fits file
    
    """    
    shift = find_shift_images(im, im_flat, polyfits, apolyfits)
    spectra, good_px_mask = extract_orders(params, im, polyfits, xlows, xhighs, widths, params['extraction_width_multiplier'], offset=shift)#, var='prec')
    aspectra, agood_px_mask = extract_orders(params, im, apolyfits, axlows, axhighs, awidths, params['arcextraction_width_multiplier'], offset=shift, var='fast')
    espectra = measure_noise(spectra)
    #print wavelength_solution[0]
    wavelength_solution_shift = shift_wavelength_solution(aspectra, wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, xlows, xhighs)
    wavelengths = create_wavelengths_from_solution(wavelength_solution_shift, spectra)
    #print wavelength_solution[0]
    fspectra = spectra/flat_spec_norm[1]    # a routine is necessary to count for shifts in the arc_solution!!!!
    efspectra = measure_noise(fspectra)
    #cspectra, sn_cont = normalise_continuum(fspectra, wavelengths)      # if flat correction works, than this is better, Problem, if flat is in part with nearly no light
    cspectra, sn_cont = normalise_continuum(spectra, wavelengths)        # !!! Testing, as flat correction doesn't work
    logger('Step: Saving extracted spectra')
    im_name=im_name.replace('.fits','')
    im_name=im_name.replace('.fit','')
    im_head = add_specinfo_head(spectra, im_head)
    im_head['Comment'] = 'File contains a 3d array with the following data in the form [data type, order, pixel]: '
    im_head['Comment'] = ' 0: wavelength for each order and pixel'
    im_head['Comment'] = ' 1: extracted spectrum'
    im_head['Comment'] = ' 2: measure of error (missing)'
    im_head['Comment'] = ' 3: flat corrected spectrum'
    im_head['Comment'] = ' 4: error of the flat corrected spectrum (missing)'
    im_head['Comment'] = ' 5: continuum normalised spectrum'
    im_head['Comment'] = ' 6: S/N in the continuum'
    im_head['Comment'] = ' 7: Mask with good areas of the spectrum: 0.1=saturated_px, 0.2=badpx'
    im_head['Comment'] = ' 8: spectrum of the emission line lamp'
    save_multispec([wavelengths, spectra, espectra, fspectra, efspectra, cspectra, sn_cont, good_px_mask, aspectra], params['path_extraction']+im_name, im_head, bitpix=params['extracted_bitpix'])
    """
    # Create a linearised solution
    wavelenght_solution_lin, spectra_lin = linearise_wavelength_spec(params, wavelength_solution_shift, [spectra, fspectra, aspectra])
    wavelengths = create_wavelengths_from_solution(wavelength_solution_lin, spectra_lin[0])
    print spectra_lin.shape
    spectra_lin = np.append(spectra_lin, wavelengths, axis=0)
    print spectra_lin.shape
    save_multispec(spectra_lin, params['path_extraction']+im_name+'_lin', im_head)
    """
    
    # For easier plotting
    add_text_to_file(params['path_extraction']+im_name+'.fits', 'plot_files.lst')
    
def plot_traces_over_image(im, fname, pfits, xlows, xhighs, widths=[], w_mult=1, mask=[], frame=None, return_frame=False):
    """
    Plot the found traces over the CCD image in order to allow error dedection
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
    title = 'Plot the traces in the image (log10 of image)'
    plot_img_spec.plot_image(im, [], pctile=1, show=False, adjust=[0.05,0.95,0.95,0.05], title=title, return_frame=True, frame=frame, autotranspose=False, colorbar=colorbar)
    for pp, pf in enumerate(pfits):
            if mask[pp] == False:
                continue
            xarr = np.arange(xlows[pp], xhighs[pp], 1)
            yarr = np.polyval(pf, xarr)
            yarrl = yarr - widths[pp][2]*w_mult
            yarrr = yarr + widths[pp][2]*w_mult
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
    :param wavelength_solution: list, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    """
    step = 25
    im = np.zeros([max(xhighs)/step, len(xlows)])+np.nan
    for order in range(len(xlows)):
        xarr = np.arange(xlows[order], xhighs[order], step)
        yarr = np.polyval(wavelength_solution[order,2:], xarr-wavelength_solution[order,1])
        diff_wave = yarr[1:] - yarr[:-1]
        im[xarr[:-1]/step, order ] = diff_wave
        diff_wave2 = diff_wave[1:] - diff_wave[:-1]
        if len(diff_wave2[diff_wave2>0]) > 0 and len(diff_wave2[diff_wave2<0]) > 0:     # There should'nt be any stationary points
            logger('Warn: The wavelength solution for oder {0} contains at least one stationary point. Please check {1}'.format(order, fname))
    colorbar = True  
    title = 'Plot the wavelength difference between every {0} pixel (Angstrom/{0}px)'.format(step)
    axis_name = ['Order', 'Position in Dispersion direction [{0} px]'.format(step), 'wavelength difference [Angstrom/{0}px]'.format(step)]
    plot_img_spec.plot_image(im, [fname], pctile=0, show=False, adjust=[0.05,0.95,0.95,0.05], title=title, autotranspose=False, colorbar=colorbar, axis_name=axis_name)

def plot_wavelength_solution_spectrum(spec1, spec2, fname, wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names):
    """
    Creates a pdf with one page for each ThAr spectrum and adds the identified and further catalogue lines to the plots
    :param spec1: 2d array of floats, extracted spectrum of the long exposed arc
    :param spec2: 2d array of floats, extracted spectrum of the short exposed arc
    :param fname: string, Filename to which the image is saved
    :param wavelength_solution: 2d array, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    :param wavelength_solution_arclines: list, same length as number of orders, each line contains the wavelength of the used reference lines
    :param reference_catalog: 2d array of floats with one entry for each line. Each entry contains the wavelength and the intensity of the line (set to 1 if not provided by the catalog file), and the index in the file
    :param reference_names: list of strings, same length as reference_catalog. Names of the reference lines
    """
    max_reflines = 30
    plot_pos = [0.01, 0.035, 0.04]      # begin line, end line, beginn of text; all relative to the the flux at the position and proportional to the range of the flux
    
    reference_catalog = np.array(sorted(reference_catalog, key=operator.itemgetter(1), reverse=True ))          # Sort by intensity in reverse order
    wavelengths = create_wavelengths_from_solution(wavelength_solution, spec1)
    
    # multiple pdf pages from https://matplotlib.org/examples/pylab_examples/multipage_pdf.html
    with PdfPages(fname) as pdf:
        labels = ['long exp','short exp']
        x_title = 'Wavelength [Angstroms]'
        y_title = 'extracted flux [ADU]'
        titel_f = 'Order {0}, real Order {1}\nmarked (proportional to line strength) are the identified reference lines (red, {2} lines) and\na subset of the omitted reference lines (green, showing the brightes {3} out of {4} lines)'
        for order in range(wavelength_solution.shape[0]):
            x_range = wavelengths[ order, ~np.isnan(spec1[order,:]) ]
            good_values = np.where((reference_catalog[:,0] >= min(x_range)) & (reference_catalog[:,0] <= max(x_range)) )[0]     # lines in the right wavelength range
            reference_lines = reference_catalog[good_values,:]                                                                  # possible reference lines in the order
            
            num_ident = len(wavelength_solution_arclines[order])
            if num_ident == 0:                                          # Add a dummy value in case no lines are available for this order, as otherwise it will crash
                wavelength_solution_arclines[order].append(-1E5)
            num_notident = len(reference_lines)
            
            title = titel_f.format(order, int(wavelength_solution[order,0]), num_ident, min(max_reflines,num_notident-num_ident), num_notident-num_ident)
            fig, frame = plt.subplots(1, 1)
            x_data = [wavelengths[order,:], wavelengths[order,:]]
            y_data = [spec1[order,:], spec2[order,:]]
            plot_img_spec.plot_points(x_data, y_data, labels, [], show=False, adjust=[0.12,0.88,0.85,0.12, 1.0,1.01], title=title, 
                                      return_frame=True, frame=frame, x_title=x_title, y_title=y_title, linestyle="-", marker="")
            axes = plt.gca()
            ymin, ymax = axes.get_ylim()
            y_range = ymax - ymin
            #xmin, xmax = axes.get_xlim()        # wavelengths contains also values for nans
            
            if len(reference_lines) > 0:
                y_scale = (plot_pos[1] - plot_pos[0]) / max(reference_lines[:,1])
                num_notident = 1
                for refline in reference_lines:
                    if np.min(np.abs(refline[0] - wavelength_solution_arclines[order])) < 1E-5:
                        color = 'r'
                    else:
                        if num_notident > max_reflines:
                            continue
                        num_notident += 1    
                        color = 'g'
                    x_position = np.argmin(np.abs(refline[0] - wavelengths[order,:]))
                    if np.isnan(spec1[order, x_position]):
                        continue
                    y_position = np.nanmax(spec1[order, max(0,x_position-2):min(len(spec1[order,:]),x_position+2)])
                    ymax = max(ymax, y_position+0.23*y_range)
                    frame.plot( [refline[0],refline[0]], [y_position+plot_pos[0]*y_range,y_position+(plot_pos[0]+y_scale*refline[1])*y_range], color=color )
                    frame.text( refline[0], y_position+plot_pos[2]*y_range, '{0} - {1}'.format(round(refline[0],4),reference_names[int(refline[2])]), 
                                horizontalalignment='center', verticalalignment='bottom', rotation=90, color=color, zorder=5 )
                    
            axes.set_ylim(ymin,ymax)
            #axes.set_xlim(xmin,xmax)
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
    :param wavelength_solution: list, same length as number of orders, each line consists of the real order, central pixel, and the polynomial values of the fit
    :param wavelength_solution_arclines: list, same length as number of orders, each line contains the wavelength of the used reference lines. The lines don't have the same length
    :param reference_catalog: 2d array of floats with one entry for each line. Each entry contains the wavelength and the intensity of the line
    """
    logger('Step: logging results to {0}'.format(fname))
    # Prepare the refernece catalog to show only limit number of lines
    wavelength_solution_arclines_flattened = [item for sublist in wavelength_solution_arclines for item in sublist]
    min_wav = np.polyval(wavelength_solution[-1,2:], np.polyval(pfits[-1], xlows[-1])-wavelength_solution[-1,1])        # minimum wavelength on the chip
    max_wav = np.polyval(wavelength_solution[0,2:], np.polyval(pfits[0], xhighs[0])-wavelength_solution[0,1])           # maximum wavelength on the chip
    d_wav = 0.05*(max_wav - min_wav) 
    good_values = np.where((reference_catalog_full[:,0] >= min_wav -d_wav) & (reference_catalog_full[:,0] <= max_wav+d_wav) )[0]     # lines in the right wavelength regime
    reference_catalog_full = reference_catalog_full[good_values,:]                                                      # lines in the right wavelength regime
    reference_catalog_sorted = np.array(sorted(reference_catalog_full, key=operator.itemgetter(1), reverse=True ))      # Sort the reference catalogue by brightness of lines
    reference_catalog_not_ident = []                                                                                    # Will store only the lines not identified
    for entry in reference_catalog_sorted:
        if entry[0] not in wavelength_solution_arclines_flattened:
            reference_catalog_not_ident.append(entry)                                                                   # non-identified line
    reference_catalog_not_ident = np.array(reference_catalog_not_ident)
    reference_catalog = reference_catalog_not_ident[:min(len(reference_catalog_not_ident),len(pfits)*50) , :]           # only use about 70 lines per order
    
    # Plot the date
    im = scale_image_plot(im,'log10')
    ims = im.shape
    fig, frame = plt.subplots(1, 1)
    colorbar = True  
    title = 'Plot of the identified and used reference lines (red, {2} lines) and a subset of the omitted reference lines (green, showing the brightes {0} out of {1} lines) to the image (log10 of image)'.format(len(reference_catalog), len(reference_catalog_not_ident), len(wavelength_solution_arclines_flattened))
    plot_img_spec.plot_image(im, [], pctile=1, show=False, adjust=[0.05,0.95,0.95,0.05], title=title, return_frame=True, frame=frame, autotranspose=False, colorbar=colorbar)
    axes = plt.gca()
    ymin, ymax = axes.get_ylim()
    xmin, xmax = axes.get_xlim()
    for pp in range(len(pfits)):
            xarr_fine = np.arange(xlows[pp], xhighs[pp], 0.2)
            yarr_fine_p5 = np.polyval(pfits[pp], xarr_fine) + 5
            yarr_fine_m5 = np.polyval(pfits[pp], xarr_fine) - 5
            # wavelength_solution_arclines contains the catalog wavelength -> make them into px
            xarr_wave = np.polyval(wavelength_solution[pp,2:], xarr_fine-wavelength_solution[pp,1])
            min_wave, max_wave = min(xarr_wave), max(xarr_wave)
            for source in [1,2]:
                if source == 1:
                    source_data = reference_catalog[:,0]
                    color = 'g'
                elif source == 2:
                    source_data = wavelength_solution_arclines[pp]
                    color = 'r'
                for arc_line_wave in source_data:
                    if arc_line_wave < min_wave or arc_line_wave > max_wave:          # Dummy lines to fill the array
                        continue
                    min_diff = (np.abs(arc_line_wave - xarr_wave)).argmin()
                    xarc = xarr_fine[min_diff]
                    yarc_p5 = yarr_fine_p5[min_diff]
                    yarc_m5 = yarr_fine_m5[min_diff]
                    #print pp, arc_line_wave, xarc, yarc        # Can be used to get the positions of the lines on the CCD
                    frame.plot( [yarc_m5,yarc_p5],[xarc,xarc],color=color )
                    frame.text(yarc_p5, xarc, '{0}'.format(round(arc_line_wave,2)), horizontalalignment='left', verticalalignment='center', rotation=0, color=color, zorder=5)
    axes.set_ylim(ymin,ymax)
    axes.set_xlim(xmin,xmax)
    fig.set_size_inches(42, 42)
    plt.savefig(fname, bbox_inches='tight')
    plt.close()


            
def find_bck_fit(im, im_bck_px, p_orders, GUI=True):         # Needs to be rewritten using 2d image fit.
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
            stddiff = np.std(difffit, ddof=p_orders[1])
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

def bck_px_UI(params, im_orig, pfits, xlows, xhighs, widths, w_mult, userinput=True):
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
            print 'Warn: order outside the allowed area: 0...{0}'.format(im.shape[0]-1)
            order = oldorder
        if order <> oldorder:
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
        xarr_arc = (reference_catalog[:,0]-arc_setting[1])/arc_setting[0]+arc_setting[2]
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
        #if x_range <> (0.0, 1.0) and y_range <> (0.0, 1.0):
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
            print 'Warn: order outside the allowed area: 0...{0}'.format(im.shape[0]-1)
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
        #print 'i,p, [len(temp), np.std(temp, ddof=polynom_order_intertrace_i),min(temp),max(temp)]',i,p, [len(temp), np.std(temp, ddof=polynom_order_intertrace_i),min(temp),max(temp)]
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
        #print 'i,p, [len(temp), np.std(temp, ddof=polynom_order_intertrace),min(temp),max(temp)]',i,p, [len(temp), np.std(temp, ddof=polynom_order_intertrace),min(temp),max(temp)]
        polyfitparams.append(p)
        # Replace the values, which are not in good_values_all with the fit
        polyfit_pl[~good_values_all,i] = np.polyval(p, polyfit_o[~good_values_all])
    return polyfit_pl, polyfitparams

def remove_bad_identified_lines(arc_lines_wavelength, polyfit_o, polyfit_p, cen_px, sigma_diff_fit):       # only used by correlate_px_wavelength, which is only used by ident_arc.py
    """
    Sigmaclips the identified reference lines
    :param arc_lines_wavelength
    """
    arc_lines_wavelength[:, 3] = 0
    for i,order in enumerate(polyfit_o):
        orderpos = (arc_lines_wavelength[:,0] == order)
        xarr = arc_lines_wavelength[orderpos, 1] - cen_px[order]   #to get the wavelength in the center of the order
        yarr = arc_lines_wavelength[orderpos, 2]        # wavelength in Angstrom
        diff_fit = yarr-np.polyval(polyfit_p[i], xarr)
        arc_lines_wavelength[orderpos, 3] = diff_fit
    std_diff_fit = np.std(arc_lines_wavelength[(arc_lines_wavelength[:, 3]<>0), 3], ddof=1) * sigma_diff_fit            # Is ddof=1 right?
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
            cen_px = ims[1]/2 + cen_px_curv*orders + cen_px_offset
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
                        #print 'order,p, [len(diff_fit), np.std(diff_fit, ddof=polynom_order_trace),min(diff_fit),max(diff_fit)]', order,p, [len(diff_fit), np.std(diff_fit, ddof=polynom_order_trace),min(diff_fit),max(diff_fit)]
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
                    #print arc_lines_wavelength[arc_lines_wavelength <> old_arc_lines_wavelength], old_arc_lines_wavelength[arc_lines_wavelength <> old_arc_lines_wavelength]
                    if (arc_lines_wavelength[:,0:3] == old_arc_lines_wavelength[:,0:3]).all():
                        break
                old_arc_lines_wavelength = arc_lines_wavelength
                #print 'arc_lines_wavelength.shape, arc_lines_wavelength', arc_lines_wavelength.shape, arc_lines_wavelength
        
            # Get the real order numbers using the grating equation: wavelength propotional to 1/(real_order) -> y(order)=(order_offset+order)*central_wavelength should have the smallest slope
            # Tested that the slope is independend of the central pixel (e.g ims[1]/2)
            slope_real_order = []
            for order_offset in range(-200,201):
                y = (order_offset+polyfit_o)*polyfit_p[:,-1]       #(order_offset+order)*central_wavelength
                good_values, pfit = sigma_clip(polyfit_o, y, 1, 3, 3, repeats=20)  #orders, sigma low, sigma high
                p = np.polyfit(polyfit_o[good_values], y[good_values], 1)
                slope_real_order.append([abs(p[0]), order_offset, len(polyfit_o[good_values]) ])
            slope_real_order = sorted(slope_real_order)
            order_offset = slope_real_order[0][1]
            order_grating_equation = 2      #should be one
            good_values, pfit = sigma_clip(polyfit_p[:,-1], 1.0/(order_offset+polyfit_o), order_grating_equation, 3, 3, repeats=20)  #orders, sigma low, sigma high
            p_real_cent = np.polyfit(1.0/(order_offset+polyfit_o[good_values]), polyfit_p[good_values,-1], order_grating_equation)
            diff_real_cent=polyfit_p[good_values,-1]-np.polyval(p_real_cent, 1.0/(order_offset+polyfit_o[good_values]))
            #print 'p_real_cent, [len(diff_real_cent), np.std(diff_real_cent, ddof=order_grating_equation),min(diff_real_cent),max(diff_real_cent)]',p_real_cent, [len(diff_real_cent), np.std(diff_real_cent, ddof=order_grating_equation),min(diff_real_cent),max(diff_real_cent)]
            text = ''
            if order_grating_equation > 1:
                text = '(including distortion, order = {0}) '.format(order_grating_equation)
            logger('Info: the offset to real orders is {0}. With this offset, the standard deviation between the central wavelengths and the grating equation {2}is {1} Angstrom'.format(order_offset, round(np.std(diff_real_cent, ddof=order_grating_equation),4),text) )
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
                    #print order,p, [len(diff_fit), np.std(diff_fit, ddof=polynom_order_trace),min(diff_fit),max(diff_fit)]
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
            good_arc_lines = arc_lines_wavelength[(arc_lines_wavelength[:, 3]<>0), :]
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


def create_pseudo_wavelength_solution(number_orders):
    wavelength_solution, wavelength_solution_arclines = [], []
    for order in range(number_orders):
        wavelength_solution.append([order, 0, 1, 0 ])
        wavelength_solution_arclines.append([])     #The wavelength of the reference lines
    return np.array(wavelength_solution), wavelength_solution_arclines

def adjust_wavelength_solution(params, spectrum, arc_lines_px, wavelength_solution_ori, wavelength_solution_arclines_ori, reference_catalog, reference_names, xlows, xhighs, show_res=False):
    """
    :return wavelength_solution_arclines: list with one entry for each order. Each order contains the wavelengths of the identified reference lines, sorted by the arc_lines_px (brightest to faintest)
    """
    ignoreorders = [] #[0,1,2,3,4,5,6,7,8]
    iter_break = 10     # Normally the script stops already at step 4-8, as no further improvements are possible
    steps = 5           # Standard: 5
    sigma = 3.5         # Standard: 3.5, the smaller the better the solution, but lines at the corners of the images might be rejected
    specs = spectrum.shape
    orders = np.arange(specs[0])
    if len(arc_lines_px) <= 10:
        logger('Warn: no arc lines available -> creating a pseudo solution (1 step per px)')
        wavelength_solution, wavelength_solution_arclines = create_pseudo_wavelength_solution(specs[0])
        return wavelength_solution, wavelength_solution_arclines
    max_diff = wavelength_solution_ori[-1][-2] * params['px_offset'][2]  # Assign lines only to resolution * px_offset
    orderdiffs = range(min(params['order_offset'][0:2]),max(params['order_offset'][0:2])+1)
    pxdiffs = range(min(params['px_offset'][0:2]),max(params['px_offset'][0:2])+1,params['px_offset'][2])
    pxdiffords = range(min(params['px_offset_order'][0:2]),max(params['px_offset_order'][0:2])+1,params['px_offset_order'][2])
    resolution_offset = params['resolution_offset_pct']/100.
    res_steps = 11
    resdiffs = np.linspace(1-resolution_offset, 1+resolution_offset, res_steps)
    best_matching = []
    print 'Compare old wavelength solution with current arc lines'
    for orderdiff in tqdm(orderdiffs):
        for pxdiff in pxdiffs:
            for pxdifford in pxdiffords:
                for resdiff in resdiffs:
                    matching_lines = []
                    # Caluculate the wavelength of all arclines with the current settings
                    arc_lines_wavelength = []
                    for order_arcline in range(int(min(arc_lines_px[:,0])),int(max(arc_lines_px[:,0])+1)):
                        if order_arcline+orderdiff < 0 or order_arcline+orderdiff > len(wavelength_solution_ori)-1 or order_arcline in ignoreorders:
                            # No solution is available for this shifted order
                            continue
                        wls = list(wavelength_solution_ori[order_arcline+orderdiff][:-2])     # These are the higher orders, which could cause problems if the shift is -700 px
                        #wls = wavelength_solution_ori[order_arcline+orderdiff][:2]     # this produces worse results
                        wls.append(wavelength_solution_ori[order_arcline+orderdiff][-2]*resdiff)
                        wls.append(wavelength_solution_ori[order_arcline+orderdiff][-1])
                        arc_lines_order_px = arc_lines_px[arc_lines_px[:,0] == order_arcline,1]    #array with px position for this order
                        arc_lines_order_wl = np.polyval(wls[2:], arc_lines_order_px-wls[1]+pxdiff+pxdifford*order_arcline)
                        # Get the offset to the closest lines
                        for arc_line_wave in arc_lines_order_wl:
                            diff_catalog = reference_catalog[:,0] - arc_line_wave
                            diff_catalog = diff_catalog[(abs(diff_catalog) <= max_diff)]
                            for entry in diff_catalog:
                                matching_lines.append([order_arcline, arc_line_wave, entry])
                    if matching_lines == []:
                        continue
                    matching_lines = np.array(matching_lines)
                    #print matching_lines.shape, orderdiff, pxdiff, resdiff, np.median(matching_lines[:,2]), np.std(matching_lines[:,2], ddof=1)
                    best_matching.append([orderdiff, pxdiff, pxdifford, resdiff, len(matching_lines), np.std(matching_lines[:,2], ddof=1)])
    if best_matching == []:
        if len(wavelength_solution_ori) == specs[0]:
            logger('Warn: No matching configuration of the lines in the arc with the old wavelength solution found. Therefore the old solution will be used')
            return np.array(wavelength_solution_ori), wavelength_solution_arclines_ori
        else:
            logger('Warn: No matching configuration of the lines in the arc with the old wavelength solution found. Additionally the number of orders in the original solution and this setting does not match -> creating a pseudo solution (1 step per px)')
            wavelength_solution, wavelength_solution_arclines = create_pseudo_wavelength_solution(specs[0])
            return wavelength_solution, wavelength_solution_arclines
    best_matching = np.array(best_matching)
    best_matching = best_matching[np.argsort(best_matching[:,4], axis=0),:]     #sort by number of identified arclines -> best at the end
    #print best_matching[-10:,:]
    best_matching = best_matching[best_matching.shape[0]*3/4:,:]     #only 25% of the best values
    best_matching = best_matching[np.argsort(best_matching[:,5], axis=0),:]     #sort by standard deviation of the shift of the arc lines -> best at the front
    #print best_matching[:10]
    # Use this solution to identify lines and fit polynomials in each order
    orderdiff, pxdiff, pxdifford, resdiff = int(best_matching[0,0]), int(best_matching[0,1]), best_matching[0,2], best_matching[0,3]
    logger('Info: To match the most lines in the arc with the old wavelength solution, a shift of {0} orders, a multiplier to the resolution of {1}, a shift of {2} px, and a shift of {3} px per order needs to be applied. {4} lines were identified. The deviation is {5} Angstrom.'.format(orderdiff, resdiff, pxdiff, pxdifford, int(best_matching[0,4]), round(best_matching[0,5],4) ))
    ## assign lines in the arcline with lines from the reference catalog
    arc_lines_wavelength = []       #order, central pixel of the arc line, wavelength of the assigned reference line, diff between both solutions, height of the arc line, intensity of the reference line, 1/width of the arc line
    for order_arcline in range(int(min(arc_lines_px[:,0])),int(max(arc_lines_px[:,0])+1)):
        if order_arcline+orderdiff < 0 or order_arcline+orderdiff > len(wavelength_solution_ori)-1 or order_arcline in ignoreorders:
            # No solution is available for this shifted order
            continue
        wls = list(wavelength_solution_ori[order_arcline+orderdiff][:-2])     # These are the higher orders, which could cause problems if the shift is -700 px
        #wls = wavelength_solution_ori[order_arcline+orderdiff][:2]
        wls.append(wavelength_solution_ori[order_arcline+orderdiff][-2]*resdiff)
        wls.append(wavelength_solution_ori[order_arcline+orderdiff][-1])
        arc_lines_order_px = arc_lines_px[arc_lines_px[:,0] == order_arcline,:]    #[:,1] is px position, [:,2] is width of lines, [:,3] is height of the line
        arc_lines_order_wl = np.polyval(wls[2:], arc_lines_order_px[:,1]-wls[1]+pxdiff+pxdifford*order_arcline)
        for i, arc_line_wave in enumerate(arc_lines_order_wl):
            diff_catalog = reference_catalog[:,0] - arc_line_wave
            good_values = (abs(diff_catalog) <= max_diff)
            for entry in reference_catalog[good_values,:]:
                arc_lines_wavelength.append([order_arcline, arc_lines_order_px[i,1], entry[0],0, arc_lines_order_px[i,3], entry[1], arc_lines_order_px[i,2]])
    ## Fit polynoms in and between orders
    cen_px_curv, cen_px_offset = 0, 0                           # remains to check which central pixel fits the whole solution best
    order_offset = wavelength_solution_ori[0][0]+orderdiff      # checked that it is "+orderdiff"
    opt_px_range = (1 - np.linspace(params['opt_px_range'][0], params['opt_px_range'][1], steps)) * specs[1]/2
    polynom_order_traces = (np.linspace(params['polynom_order_traces'][0], params['polynom_order_traces'][1], steps)).astype(int)
    polynom_order_intertraces = params['polynom_order_intertraces']
    if len(polynom_order_intertraces) == 1:
        polynom_order_intertraces.append(params['polynom_order_intertraces'])
    polynom_order_intertraces = (np.linspace(polynom_order_intertraces[0], polynom_order_intertraces[1], steps)).astype(int)
    max_diff_pxs = np.linspace(params['px_offset'][2], params['diff_pxs'], steps)
    cen_px = specs[1]/2 + cen_px_curv*orders + cen_px_offset    # remains to check which central pixel fits the whole solution best
    arc_lines_wavelength = np.array(arc_lines_wavelength)
    old_arc_lines_wavelength = arc_lines_wavelength
    for step in range(steps):
        max_diff_px = max_diff_pxs[step]
        px_range = opt_px_range[step]
        polynom_order_trace = polynom_order_traces[step]
        polynom_order_intertrace = polynom_order_intertraces[step]
        #print step,max_diff_px,px_range,polynom_order_trace
        for iteration in range(iter_break):
                x = arc_lines_wavelength[:,1]-cen_px[0]
                y = 1.0/(arc_lines_wavelength[:,0]+order_offset)
                # Weight: height of the arc line + log10(intensity of the reference line) + 1/width of the arc line
                medians = [np.median(arc_lines_wavelength[:,4]), np.median(arc_lines_wavelength[:,5]), np.median(arc_lines_wavelength[:,6])]
                for i in range(len(medians)):
                    if medians[i] == 0:
                        medians[i] = np.mean(arc_lines_wavelength[:,i+4])
                    if medians[i] == 0:
                        medians[i] = 1
                weight = arc_lines_wavelength[:, 4]/medians[0]+arc_lines_wavelength[:,5]/medians[1]+arc_lines_wavelength[:,6]/medians[2]
                poly2d_params = polynomial_fit_2d_norm(x, y, arc_lines_wavelength[:,2], polynom_order_trace, polynom_order_intertrace, w=weight)
                x = arc_lines_px[:,1]-cen_px[0]
                y = 1.0/(arc_lines_px[:,0]+order_offset)
                arc_line_wave = polynomial_value_2d(x, y, polynom_order_trace, polynom_order_intertrace, poly2d_params)
                arc_line_res = np.abs(polynomial_value_2d(x+1, y, polynom_order_trace, polynom_order_intertrace, poly2d_params) - polynomial_value_2d(x-1, y, polynom_order_trace, polynom_order_intertrace, poly2d_params))/2
                arc_lines_wavelength = []
                for i in range(arc_line_wave.shape[0]):
                    if arc_lines_px[i,1] < px_range or arc_lines_px[i,1] > specs[1]-px_range or arc_lines_px[i,0] in ignoreorders:
                        # Do not use lines in the outer parts of the spectrum at the beginning
                        continue
                        # Wavelength difference to the catalog lines
                    diff_catalog = reference_catalog[:,0] - arc_line_wave[i]    # Diff between catalog wavelength and fitted wavelength
                    #print min(abs(diff_catalog)),arc_line_res[i]
                    if min(abs(diff_catalog)) <= max_diff_px*arc_line_res[i]:   # calculate the approximate wavelength difference for max_diff_px
                        # Best matching line
                        pos_min = np.argmin(abs(diff_catalog))
                        catalog_line = reference_catalog[pos_min,:]
                        weight1, weight2, weight3 = arc_lines_px[i,3], np.log10(catalog_line[1]), 1/arc_lines_px[i,2]     #height of the arc line, log10(intensity of the reference line), 1/width of the arc line
                        arc_lines_wavelength.append([arc_lines_px[i,0], arc_lines_px[i,1], catalog_line[0], diff_catalog[pos_min], weight1, weight2, weight3, arc_line_res[i], pos_min, arc_lines_px[i,2] ])
                arc_lines_wavelength = np.array(arc_lines_wavelength)   
                # arc_lines_wavelength: order, px in spectrum, wavelength in reference catalog, difference in wavelength between wavelength solution and reference line, 
                #                       3x weights, 1px resolution at that position, index in refence catalog, width of the line
                goodvalues, p = sigma_clip(np.arange(arc_lines_wavelength.shape[0]), arc_lines_wavelength[:,3], 0, sigma, sigma, repeats=20)    # poly_order=0 to sigma clip around the average
                arc_lines_wavelength = arc_lines_wavelength[goodvalues,:]
    
                #print max_diff_px, arc_lines_wavelength.shape[0],old_arc_lines_wavelength.shape[0],arc_lines_wavelength.shape[0] == old_arc_lines_wavelength.shape[0],np.mean(arc_lines_wavelength[:,3]), np.mean(old_arc_lines_wavelength[:]),np.mean(arc_lines_wavelength[:,3]) == np.mean(old_arc_lines_wavelength[:]), np.min(arc_lines_wavelength[:,3]), np.max(arc_lines_wavelength[:,3])
                if arc_lines_wavelength.shape[0] == old_arc_lines_wavelength.shape[0]:
                    if np.mean(arc_lines_wavelength[:,3]) == np.mean(old_arc_lines_wavelength[:]):
                        break
                old_arc_lines_wavelength = copy.copy(arc_lines_wavelength[:,3])
                
                #print arc_lines_wavelength.shape, poly2d_params
                """ # Only for printing the progress
                x,y,l=[],[],[]
                for order in orders:
                    x.append(arc_lines_wavelength[(arc_lines_wavelength[:,0]==order),1])
                    y.append(arc_lines_wavelength[(arc_lines_wavelength[:,0]==order),3])
                    l.append('{0}'.format(order))
                if iteration == 100:
                    plot_img_spec.plot_points(x,y,l,'path',show=True, x_title='Pixel', y_title='Residual')          """
    
    # Get the real order numbers using the grating equation: wavelength propotional to 1/(real_order) -> y(order)=(order_offset_new+order)*central_wavelength should have the smallest slope
    order_offset_old = int(order_offset)
    slope_real_order = []
    cenwave = polynomial_value_2d(orders*0, 1.0/(orders+order_offset), polynom_order_trace, polynom_order_intertrace, poly2d_params)    #old order_offset
    for order_offset in range(-200,201):
        y = (order_offset+orders)*cenwave   #(order_offset+order)*central_wavelength
        p = np.polyfit(orders, y, 1)
        slope_real_order.append([abs(p[0]), order_offset, len(orders) ])
    slope_real_order = sorted(slope_real_order)
    order_offset = slope_real_order[0][1]
    if 0.0 in order_offset+orders:            # Definitively the wrong order offset
        order_offset = order_offset_old
    p_real_cent = np.polyfit(1.0/(order_offset+orders), cenwave, 1)
    diff_real_cent = cenwave-np.polyval(p_real_cent, 1.0/(order_offset+orders))
    # Redo the 2d polynomial fit (it should also re-create the arc_lines_wavelength table, but most times the difference is 0)
    x = arc_lines_wavelength[:,1]-cen_px[0]
    y = 1.0/(arc_lines_wavelength[:,0]+order_offset)
    
    medians = [np.median(arc_lines_wavelength[:,4]), np.median(arc_lines_wavelength[:,5]), np.median(arc_lines_wavelength[:,6])]
    for i in range(len(medians)):
        if medians[i] == 0:
            medians[i] = np.mean(arc_lines_wavelength[:,i+4])
        if medians[i] == 0:
            medians[i] = 1
    weight = arc_lines_wavelength[:, 4]/medians[0]+arc_lines_wavelength[:,5]/medians[1]+arc_lines_wavelength[:,6]/medians[2]
    poly2d_params = polynomial_fit_2d_norm(x, y, arc_lines_wavelength[:,2], polynom_order_trace, polynom_order_intertrace, w=weight)
    
    if len(arc_lines_wavelength) <= polynom_order_trace+polynom_order_intertrace:
        logger('Error: Not enough lines are remaining for the fit. Please check the parameters "order_offset", "px_offset", and "px_offset_order".')
    
    # Log the results
    # print np.vstack([ arc_lines_wavelength[:,2], arc_lines_wavelength[:,3], arc_lines_wavelength[:,-3], arc_lines_wavelength[:,2]/arc_lines_wavelength[:,3] ]).T
    avg_diff_fit = np.mean(np.abs(arc_lines_wavelength[:,3]))           # Diff between catalog wavelength and fitted wavelength
    std_diff_fit = np.std(arc_lines_wavelength[:,3], ddof=(polynom_order_trace+polynom_order_intertrace))
    # Resolution from the precission of the fit
    std_R_fit = 1/np.std(arc_lines_wavelength[:,3]/arc_lines_wavelength[:,2], ddof=(polynom_order_trace+polynom_order_intertrace))    # lambda/d_lambda
    # Resolution using the Gaussian width of the arc lines
    R_gauss     = arc_lines_wavelength[:,2]/(arc_lines_wavelength[:,-1]*arc_lines_wavelength[:,-3])    # lambda/d_lambda ; -1 is in px, -3 is resolution
    good_values, p = sigma_clip(R_gauss, R_gauss, 0, 3, 3)   #x, y, order of polynom (if 0 than average -> x doesn't matter), sigma, sigma
    R_gauss     = R_gauss[good_values]
    avg_R_gauss = np.mean(R_gauss)
    std_R_gauss = np.std(R_gauss, ddof=1)
    # 2px Resolution (using only the identified arc lines
    R_2px     = arc_lines_wavelength[:,2]/(2*arc_lines_wavelength[:,-3])    # lambda/d_lambda
    avg_R_2px = np.mean(R_2px)
    std_R_2px = np.std(R_2px, ddof=1) 
    text = 'Info: used {0} lines. The standard deviation of the residuals between the lines and the fit is {1} Angstrom (the average of the abs of the residuals is {7} Angstrom). This converts into a resoution R = {2}. The Gaussian width of the emssion lines results in an R = {3} +- {4}. The 2-pixel resolution (around the identified lines) is R = {5} +- {6}.'
    logger(text.format(arc_lines_wavelength.shape[0], round(std_diff_fit,4), int(std_R_fit), int(avg_R_gauss), int(std_R_gauss), int(avg_R_2px), int(std_R_2px), round(avg_diff_fit,4) ))
    text = 'Info: A 2d polynom fit with {0} orders in dispersion direction (along the traces) and {1} orders in cross-dispersion direction was used. With this solution, the offset between aperture and real orders is {2}. With this offset, the standard deviation of the residuals between the central wavelengths and the grating equation is {3} Angstrom (abs of average is {5} Angstrom). Using the original solution gives an offset of {4}.'
    logger(text.format(polynom_order_trace, polynom_order_intertrace, order_offset, round(np.std(diff_real_cent, ddof=len(p_real_cent)),4), order_offset_old, round(np.mean(np.abs(diff_real_cent)),4) ))
    
    # Create a wavelength solution
    wavelength_solution = [polynom_order_trace, polynom_order_intertrace, cen_px[0], order_offset]
    wavelength_solution = np.append(wavelength_solution, poly2d_params, axis=0)
    logger('Info: Wavelenght solution in 2d (for pixel and order at the same time) is [No of orders1, No of orders2, central pixel, offset to real order, parameters of the 2d polynomial(px^0*m^0, px^1*m^0, px^0*m^1, px^2*m^0, px^1*m^1, px^0*m^2, ....)]: \n{0}'.format(wavelength_solution), show=False)
    
    x,y,l=[],[],[]
    max_number_reflines = 0             # will be needed later in order to save the data correctly into a fits file
    for order in orders:
        x.append(arc_lines_wavelength[(arc_lines_wavelength[:,0]==order),1])    # Pixel
        y.append(arc_lines_wavelength[(arc_lines_wavelength[:,0]==order),3])    # Residuals
        l.append('{0}'.format(order))
        max_number_reflines = max(max_number_reflines, len(arc_lines_wavelength[(arc_lines_wavelength[:,0]==order),2]))
    text = 'Residuals of the identificated arc lines: identified {0} lines for a 2d polynom fit with {1} and {2} orders'.format(arc_lines_wavelength.shape[0], polynom_order_trace, polynom_order_intertrace)
    plot_img_spec.plot_points(x,y,l,[params['logging_arc_line_identification_residuals']],show=False, title=text, x_title='Pixel', y_title='Residual')

    # Transform the wavelength solution into the old wavelength solution
    wavelength_solution, wavelength_solution_arclines = [], []
    xarr = np.arange(0,specs[1], 0.1)
    for order in orders:
        yarr = polynomial_value_2d(xarr-cen_px[0], 1.0/(order+order_offset), polynom_order_trace, polynom_order_intertrace, poly2d_params)
        polyfit = np.polyfit(xarr-cen_px[0], yarr, polynom_order_trace)      #lambda from px
        solution = [order_offset+order, cen_px[0] ]
        solution = np.append(solution, polyfit, axis=0)
        wavelength_solution.append(solution)
        reflines = list(arc_lines_wavelength[(arc_lines_wavelength[:,0]==order),2])     # Wavelength
        for i in range(max_number_reflines - len(reflines)):                            # Fill up with zeros in order to create an array
            reflines.append(0)
        wavelength_solution_arclines.append(reflines)     #The wavelength of the reference lines
        
        #print order,polynomial_value_2d(0, 1.0/(order+order_offset), polynom_order_trace, polynom_order_intertrace, poly2d_params), solution
        #diff_fit = yarr-np.polyval(polyfit, xarr-cen_px[0])        #<1e-10
        #print [len(diff_fit), np.std(diff_fit, ddof=polynom_order_trace),min(diff_fit),max(diff_fit)]
    wavelength_solution = np.array(wavelength_solution)
    statistics_arc_reference_lines(arc_lines_wavelength, [0,-2,-1,2], reference_names, wavelength_solution, xlows, xhighs)
        
    if std_diff_fit > max(wavelength_solution[:,-2]) or std_diff_fit < 1E-8:
        logger('Error: The wavelength solution seems wrong. Please check the parameters "order_offset", "px_offset", and "px_offset_order". It might be useful to compare the file with the emission lines ({0}) in the current folder and the folder folder with the previous wavelength solution (see parameter "original_master_arc_solution_filename")'.format(params['master_arc_l_filename']))
    
    # See the results
    if show_res:
        correlate_px_wave_result_UI(spectrum, arc_lines_wavelength, reference_catalog, arc_lines_px, reference_names, wavelength_solution, adjust=[0.07,0.93,0.94,0.06, 1.0,1.01])
    
    return wavelength_solution, wavelength_solution_arclines

def statistics_arc_reference_lines(data, positions, reference_names, wavelength_solution, xlows, xhighs):
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
            if index_r <> []:
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
                result.append([ order, cenwave, minwave, maxwave, maxwave-minwave, wavelength_solution[order,-2], refname, len(data_o[index_r,positions[2]]), np.average(data_o[index_r,positions[2]]), std, min_refline, max_refline, range_refline ])
                #print result[-1]

    for refname in uniq_refnames:
        index_r = []
        for i,entry in enumerate(data[:,positions[1]].astype(int)):
            if reference_names[entry] == refname:
                index_r.append(i)
        if index_r <> []:
            min_refline = min(data[index_r,positions[3]])
            max_refline = max(data[index_r,positions[3]])
            result.append([ -1, -1, -1, -1, -1, -1, refname, len(data[index_r,positions[2]]), np.average(data[index_r,positions[2]]), np.std(data[index_r,positions[2]], ddof=1), min_refline, max_refline, max_refline-min_refline ])
    printarrayformat = ['%1.1i', '%3.1f', '%3.1f', '%3.1f', '%3.1f', '%5.3f', '%s', '%1.1i', '%4.2f', '\t%4.2f', '\t%4.1f', '\t%4.1f', '\t%4.1f']
    logger('\t\tapert\tcenwave\tminwave\tmaxwave\tranwave\tAng/px\tname\tnumber\tgausswidth_avg\tgausswidth_std\tmin_refline\tmax_refline\trange_reflines_whole_order',printarrayformat=printarrayformat, printarray=result)

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
        if line[0]<>'' and line[1]<>'' and line[2]<>'' and line[3]<>'':
            # get the order, real order, pixel, and wavelength
            try:
                arc_lines_wavelength.append([ int(line[0]), int(line[1]), float(line[2]), float(line[3]),0 ])
            except:
                print 'Problems to convert to int/float:', line
                a=1
    file.close()
    if len(arc_lines_wavelength) == 0:
        logger('Error: No useful information was found in {0}. Please make sure the entries in the file contain of the tab-separated values: order (starting at 0), real order, pixel, wavelength'.format(filename))
    arc_lines_wavelength = np.array(arc_lines_wavelength)
    if 0.0 in arc_lines_wavelength[:,1]:
        logger('Error: The second coloumn (real order) in file {1} was not set correctly as it contains a zero. Please use the real order (from grating equation)'.format(filename))
    
    #order_offset = 117      #blue chip
    #order_offset = 89      #red chip
    order_offset = arc_lines_wavelength[:,1] - arc_lines_wavelength[:,0]
    if np.std(order_offset, ddof=1) <> 0:
        logger('Error: There seems to be a inconsistency between coloumn 1 (order) and coloumn 2 (real order) in the file {0}. Please check'.format(filename))
    order_offset = order_offset[0]
    polynom_order_trace = max(params['polynom_order_traces'])-1
    polynom_order_intertrace = max(params['polynom_order_intertraces'])-1
    
    #polynom_order_trace = 3
    #polynom_order_intertrace = 2
    cen_px = im.shape[0]/2
    #arc_line_res = arc_lines_wavelength[:,3]*0      # Keep all data at the first step
    #for i in range(5):
    #    arc_lines_wavelength = arc_lines_wavelength[(arc_line_res < 1),:]      # Problem with bad wavelengths -> rejecting nearly everything
    for i in range(arc_lines_wavelength.shape[0] - polynom_order_trace*polynom_order_intertrace):
        x = arc_lines_wavelength[:,2]-cen_px
        y = 1.0/(arc_lines_wavelength[:,1])
        weight = arc_lines_wavelength[:,0]*0+1
        poly2d_params = polynomial_fit_2d_norm(x, y, arc_lines_wavelength[:,3], polynom_order_trace, polynom_order_intertrace, w=weight)
        arc_line_wave_fit = polynomial_value_2d(x, y, polynom_order_trace, polynom_order_intertrace, poly2d_params)
        arc_line_res = np.abs(arc_line_wave_fit - arc_lines_wavelength[:,3])
        if max(arc_line_res) < 1:       # only good data left, probably should be less than 1 Angstrom
            break
        arc_lines_wavelength = np.delete(arc_lines_wavelength, arc_line_res.argmax(), 0)
    arc_line_res = np.abs(arc_line_wave_fit - arc_lines_wavelength[:,3])
    logger('Info: The average residual of the fit to the manual wavelength solution is {0} +/- {1} Angstrom. Only input data, for which the residuals were less than 1 Angstrom have been used. {2} lines have been used to calculate the solution for a 2d polynom fit with {3} orders in dispersion direction (along the traces) and {4} orders in cross-dispersion direction.'.format( round(np.mean(arc_line_res),4), round(np.std(arc_line_res, ddof=polynom_order_trace+polynom_order_intertrace),4), arc_lines_wavelength.shape[1], polynom_order_trace, polynom_order_intertrace ))
    if arc_lines_wavelength.shape[0] <= polynom_order_trace*polynom_order_intertrace:
        logger('Error: The solution seems unphysical, at least {0} lines should be used (degrees of freedom)'.format(polynom_order_trace*polynom_order_intertrace+1))
    # Transform the wavelength solution into the old wavelength solution
    wavelength_solution, wavelength_solution_arclines = [], []
    xarr = np.arange(0,max(x), 0.1)
    for order in range(int(max(arc_lines_wavelength[:,0]))+1):
        yarr = polynomial_value_2d(xarr-cen_px, 1.0/(order+order_offset), polynom_order_trace, polynom_order_intertrace, poly2d_params)
        polyfit = np.polyfit(xarr-cen_px, yarr, polynom_order_trace)      #lambda from px
        solution = [order+order_offset, cen_px ]
        solution = np.append(solution, polyfit, axis=0)
        wavelength_solution.append(solution)
        wavelength_solution_arclines.append([])
        
    wavelength_solution = np.array(wavelength_solution)
    
    return wavelength_solution, wavelength_solution_arclines

def read_text_file(filename, no_empty_lines=False):
    text = []
    if os.path.isfile(filename) == True:
        text1 = open(filename,'r').readlines()
        for line in text1:
            line = line.replace('\n', '')
            linetemp = line.replace('\t', '')
            if ( line == '' or linetemp == '') and no_empty_lines:
                continue
            text.append(line)
    else:
        logger('Warn: File {0} does not exist, assuming empty file'.format(filename))    
    return text

def add_text_file(text, filename):
    file = open(filename,'a')
    file.write(text+'\n')
    file.close()

def add_text_to_file(text, filename):
    oldtext = read_text_file(filename)
    exists = False
    for line in oldtext:
        if line.find(text) <> -1:
            exists = True
            break
    if exists == False:
        add_text_file(text, filename)
        
def remove_orders(pfits, rm_orders):
    mask = np.repeat([True], len(pfits))
    if type(rm_orders) is not list:
        return mask
    # Remove orders (if rm_orders is populated)
    for r in range(len(pfits)):
        if r in rm_orders:
            mask[r] = False
    return mask

def run_remove_orders(im1, pfits, xlows, xhighs, userinput=True):
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
    adjust=[0.05,0.95,0.95,0.05, 1.0,1.01]
    size = datapoints.shape
    for step in range((size[0]+9)/10):
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
            plt.savefig(filename.replace('.png','_{0}-{1}.png'.format('%2.2i'%datarange[0],'%2.2i'%datarange[-1])), bbox_inches='tight')
            plt.close()

def add_specinfo_head(spectra, im_head):
    spectra = np.array(spectra)
    im_head['fmin'] = np.nanmin(spectra)
    im_head['fmax'] = np.nanmax(spectra)
    im_head['fsum_all'] = np.nansum(spectra, axis=None)
    for order in range(spectra.shape[0]):
        im_head['fsum_{0}'.format('%2.2i'%order)] = np.nansum(spectra[order,:], axis=None)
    
    return im_head


def normalise_continuum(spec, wavelengths, nc=5, ll=2., lu=5., frac=0.1, window=21):      # without modal noise nc=4,ll=1.,lu=5. might work
    """
    Normalises the spectrum using the area with continuum and excluding the area with lines
    Copied from the ceres pipeline, workes for HARPS data, but not for EXOhSPEC spectra (before 20171201)
    """
    specs = spec.shape
    ccoefs, sn_coefs, cen_waves = [], [], []
    for order in range(specs[0]):
        x_range = wavelengths[order,:]
        y_range = scipy.signal.medfilt(spec[order,:], window)
        I = np.where( (y_range != 0) & (~np.isnan(y_range)) )[0]
        x_range, y_range = x_range[I], y_range[I]                       # use only good values
        x_cen = (min(x_range) + max(x_range))/2
        ori_len = len(y_range)
        cond = True
        while cond:
            #print len(y_range)
            p = polyfit_adjust_order(x_range-x_cen, y_range, nc)        # Fit the data, !!! In Ceres Pipeline without x_cen -> polynomial fit isn't as good as here
            res = y_range - np.polyval(p,x_range-x_cen)                 # Residuals
            IU = np.where(res>0)[0]                                     # Positive residuals
            dev = np.mean(res[IU])                                      # deviation using the positive residuals -> without absorption lines
            I = np.where((res<lu*dev)&(res>-ll*dev))[0]                 # Sigmacliping
            J1 = np.where(res>=lu*dev)[0]
            J2 = np.where(res<=-ll*dev)[0]
            if (len(J1)==0 and len(J2)==0) or len(y_range)<frac*ori_len:
                cond = False
            else:
                #plot_img_spec.plot_points([x_range,x_range,x_range[I]], [y_range, np.polyval(p,x_range-x_cen), y_range[I]], ['data', 'fit', 'remaining'], '', show=True, return_frame=False, x_title='wavelength [Angstrom]', y_title='flux [ADU]', title=str(order))
                x_range, y_range = x_range[I], y_range[I]               # Sigmacliping
        ccoefs.append(p)
        cen_waves.append(x_cen)
        #x_full = wavelengths[order,:]                           #only for testing
        #y_full = scipy.signal.medfilt(spec[order,:], window)    #only for testing
        #plot_img_spec.plot_points([x_full,x_full,x_range], [y_full,np.polyval(p,x_full-x_cen),y_range], ['full data', 'fit', 'data for fit'], '', show=True, return_frame=False, x_title='wavelength [Angstrom]', y_title='flux [ADU]')
        p = polyfit_adjust_order(x_range-x_cen, np.abs(res), nc)
        sn_coefs.append(p)
    ccoefs = np.array(ccoefs)
    contspec, sn_cont = [], []
    for order in range(specs[0]):
        fit_flux = np.polyval(ccoefs[order],wavelengths[order,:]-cen_waves[order])
        fit_sn = np.polyval(sn_coefs[order],wavelengths[order,:]-cen_waves[order])
        for i in range(len(fit_flux)-1):
            if fit_flux[i]*fit_flux[i+1] < 0:                           # change of sign
                fit_flux[max(0,i-window/2):min(len(fit_flux),i+window/2+1)] = np.nan     # ignore this data as artificial peaks are introduced due to the division close to 0
        bad_value = (fit_sn <= 0)                                       # Errors smaller than 0
        fit_sn[bad_value] = np.nan
        fit_flux[bad_value] = np.nan
        contspec.append(spec[order] / fit_flux)
        sn_cont.append(fit_sn)
    
    """ doesn't improve things in the overlapping area
    spec = np.array(contspec)
    ccoefs = get_cont_ceres(wavelengths ,spec)      # doesn't really work for solar spectra without a correct flat correction    
    contspec = []
    for order in range(specs[0]):
        contspec.append(spec[order] / np.polyval(ccoefs[order],wavelengths[order,:]))
    """
    return np.array(contspec), np.array(sn_cont)

def get_cont_ceres(W,F,nc=3,ll=1.,lu = 5.,frac=0.05,window=21):     # Based on ceres pipeline (copied, added comments, maybe changed code)
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

    return coefs

def linearise_wavelength_spec(params, wavelength_solution, spectra):    
    """
    Is some idea, but unfished.
    Issue: Doesn't create comparability, as the wavelength axis might be different between two spectra
    
    """
    spectra = np.array(spectra)
    if len(spectra.shape) == 2:
        spectra = np.array([spectra])
    specs = spectra.shape
    dwave = params['wavelength_scale_resolution']
    px_range = np.arange(specs[-1])
    for order in tqdm(range(specs[1])):
        wave_range = np.polyval(wavelength_solution[order,2:], px_range-wavelength_solution[order,1])
        lwave_range = np.arange(min(wave_range), max(wave_range)+dwave, dwave )
        for i in range(specs[0]):
            print order,i, wave_range.shape, lwave_range.shape
            spec = spectra[i,order,:]
            nonnan = ~np.isnan(spec)
            spec2 = spec[nonnan]
            wave2 = wave_range[nonnan]
            lwave2 = lwave_range[ (lwave_range >= min(wave2)) & (lwave_range <= max(wave2)) ]
            # spec can be replaced by Spline fitted to it
            # https://stackoverflow.com/questions/17913330/fitting-data-using-univariatespline-in-scipy-python
            # https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.interpolate.UnivariateSpline.integral.html
            spec = inter.UnivariateSpline (wave2, spec2, s=0.1)
            lspec = spec(lwave2)
            print np.mean(spec2), np.mean(lspec)
            plot_img_spec.plot_points([wave_range, wave_range, lwave_range], [spectra[i,order,:], spec(wave_range), lspec], ['orig', 'spline', 'lin'], '', show=True, return_frame=False, x_title='wavelength [Angstrom]', y_title='flux [ADU]')
    
    return wavelenght_solution_lin, spectra_lin




    
