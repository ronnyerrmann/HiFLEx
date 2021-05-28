#!/usr/bin/env python

import numpy as np
import os
from procedures import *
from astroquery.simbad import Simbad

# =============================================================================
# Define variables
# =============================================================================
params = dict()     # default param dictionary
# location of config file
CONFIGFILE = 'conf.txt'

# Start of code
# deal with arguments from a text file
params = textfileargs(params, CONFIGFILE)

def create_parameters(conf_data, warn, param, textparam, parcals, exist_bias, exist_rflat, exp_darks, entry):
    """
    Add information to conf_data dictionary
    :param conf_data: dictionary
    :param warn: list of strings: Collection of warnings
    :param param: string:           used for param+'_rawfiles' and 'master_'+param+'_filename'
    :param textparam: string:       used for 'master_'+textparam+'.fits'
    :param parcals: list of strings: which general calibrations to use: parcal+'_calibs_create_g'
    :param exist_bias: bool: Do Biases exist?
    :param exist_rflat: bool: Do rflats exist?
    :param exp_darks: list of floats
    :param entry: list of [str, str, str, float, float]: entry as in params['raw_data_file_list']
    """
    calibs = []
    for parcal in parcals:                                  # A list of calib_create can be given in case one wasn't defined by the user
        calibs.append(parcal+'_calibs_create_g')
    calibs.append('standard_calibs_create')                 # Standard calibration as backup
    calibs.append('')                                       # No calibration at all as very last option
    if param == '':
        return conf_data, warn
    if param+'_rawfiles' not in conf_data.keys():           # Do all the long steps as first file for this kind
        conf_data[param+'_rawfiles'] = entry[0]             # Filename
        if (param.find('extract') == -1 or param.find('extract_combine') == 0) and param.find('waveoffset') == -1:      # wavecal and normal extraction files don't need a master, only the combined extraction needs one
            conf_data['master_'+param+'_filename'] = 'master_'+textparam+'.fits'
        temp_param = []
        for calib in calibs:
            if calib in params.keys():
                temp_param = params[calib]
                break
        if calib == '':
            temp_param = []
            logger('Warn: Missing entry in the configuration file: {0}. Please update the configuration file.'.format(calibs[:-1]))
        text = ''
        for param_entry in temp_param:
            if (param_entry == 'bias' and exist_bias == False) or (param_entry == 'rflat' and exist_rflat == False) or (param_entry == 'dark' and entry[3] not in exp_darks):          # If the calibration data is not available, then don't do the calibration
                warn_text = 'Warn: The parameter {0}_calibs_create_g in the configuration file requires {1} correction, but this calibration data is not available.'.format(parcal,param_entry)
                if param_entry == 'dark':
                    warn_text = 'Warn: The parameter {0}_calibs_create_g in the configuration file requires {1} correction, but a dark with exposure time of {2}s is not available. \n\t\tPlease change the dark parameter to a fixed exposure time in {3}'.format(parcal, param_entry, entry[3], params['configfile_fitsfiles'] )
                if warn_text not in warn:
                    warn.append(warn_text)
                continue
            text += ',' + param_entry
        if len(text) > 0:
            text = text[1:]
        conf_data[param+'_calibs_create'] = text
    else:                                                   # Just add the file
        conf_data[param+'_rawfiles'] += ', ' + entry[0]

    return conf_data, warn

def add_new_rawfiles_file_list(params, file_list=[]):
    for raw_data_path in params['raw_data_paths']:
        for root, dirs, files in os.walk(raw_data_path, followlinks=True):
            if root.find(' ') != -1:
                logger('Warn: Folder contains a space and will be ignored: {0}'.format(root))
                continue
            for file in files:
                matchin_fileending = False                      # has the file the correct ending?
                for fileending in params['raw_data_file_endings']:
                    if file.endswith(fileending):
                        matchin_fileending = True               # right file ending
                        break
                if not matchin_fileending:
                    continue                                    # wrong file ending
                        
                #filename = os.path.join(root, file).replace(params['raw_data_path'],'')     # Only relative folder and filename
                #if not os.path.exists(params['raw_data_path'] + filename):                  # that should never be a problem
                #    continue
                filename = os.path.join(root, file)                                         # Full folders and filenames
                if file.find(' ') != -1:
                    logger('Warn: File contains a space and will be ignored: {0}'.format(file))
                    continue
                new = True                                      # is it a new file
                for entry in file_list:
                    if entry[0].find(filename) == 0 or entry[0].find(' '+filename) >= 0 or entry[0].find('#'+filename) >= 0:    # spaces are replaced in convert_readfile, so shouldn't be a problem
                        new = False
                        break
                if not new:
                    continue                                    # if not new, then no further analysis is needed
                im_head = fits.getheader(filename)
                # Identify the image type and find out what the fibers show
                fiber1, fiber2 = 'none', 'none'
                fnlow = filename.rsplit(os.sep,1)[-1].lower()   # get rid of the path
                if fnlow.find('bias') >= 0:                     # hardcoded: Bias
                    fiber1, fiber2 = 'bias', 'bias'
                if fnlow.find('dark') >= 0:                     # hardcoded: Dark
                    fiber1, fiber2 = 'dark', 'dark'
                posi  = [fnlow.find('flat') , fnlow.find('whli') , fnlow.find('white') , fnlow.find('tung') ]      # hardcoded: White light (Tungston) spectrum in science fiber
                posi2 = [fnlow.find('flat2'), fnlow.find('whli2'), fnlow.find('white2'), fnlow.find('tung2')]      # hardcoded: White light (Tungston) spectrum in calibration fiber
                for i in range(len(posi)):
                    if posi2[i] >= 0:
                        fiber2 = 'cont'
                        if posi[i] != posi2[i]:                 # both fibers
                            fiber1 = 'cont'
                    elif posi[i] >= 0:
                        fiber1 = 'cont'
                if fnlow.find('rflat') >= 0:                    # hardcoded: real flat
                    fiber1, fiber2 = 'rflat', 'rflat'
                if params['raw_data_imtyp_keyword'] in im_head.keys():
                    if im_head[params['raw_data_imtyp_keyword']].replace(' ','') == params['raw_data_imtyp_flat']:      # replace because when reading the conf file spaces are placed
                        if fnlow.find('rflat') >= 0:
                            fiber1, fiber2 = 'rflat', 'rflat'
                        else:
                            fiber1, fiber2 = 'cont', 'cont'
                posi  = [fnlow.find('arc') , fnlow.find('thar') , fnlow.find('th_ar') , fnlow.find('thorium') , fnlow.find('une') ]        # hardcoded: emission line spectrum in calibration/science fiber
                posi1 = [fnlow.find('arc1'), fnlow.find('thar1'), fnlow.find('th_ar1'), fnlow.find('thorium1'), fnlow.find('une1')]        # hardcoded: emission line spectrum in calibration fiber
                posi2 = [fnlow.find('arc2'), fnlow.find('thar2'), fnlow.find('th_ar2'), fnlow.find('thorium2'), fnlow.find('une2')]        # hardcoded: emission line spectrum in calibration fiber
                for i in range(len(posi)):
                    if posi1[i] >= 0:                           # name contains "ThAr1"
                        fiber1 = 'wave'
                    if posi2[i] >= 0:                           # name contains "ThAr2"
                        fiber2 = 'wave'
                    if posi[i] != posi1[i] and posi[i] != posi2[i]:
                        if posi[i] in [0,1,2]:                  # name starts with ThAr
                            fiber1 = 'wave'
                        elif posi[i] > 2:                       # ThAr comes so late that it's probably in calibration fiber
                            fiber2 = 'wave'
                if params['raw_data_imtyp_keyword'] in im_head.keys():
                    if im_head[params['raw_data_imtyp_keyword']].replace(' ','') == params['raw_data_imtyp_bias']:
                        fiber1, fiber2 = 'bias', 'bias'
                    if im_head[params['raw_data_imtyp_keyword']].replace(' ','') == params['raw_data_imtyp_dark']:
                        fiber1, fiber2 = 'dark', 'dark'
                    if im_head[params['raw_data_imtyp_keyword']].replace(' ','') == params['raw_data_imtyp_trace1']:
                        fiber1, fiber2 = 'cont', 'dark'
                    if im_head[params['raw_data_imtyp_keyword']].replace(' ','') == params['raw_data_imtyp_blaze']:
                        fiber1, fiber2 = 'cont', 'wave'        # for HARPS it is cont, cont
                    if im_head[params['raw_data_imtyp_keyword']].replace(' ','') == params['raw_data_imtyp_trace2']:
                        fiber1, fiber2 = 'wave', 'wave'
                if (fnlow.find('/arc') == -1 and fnlow.find('/thar') == -1  and fnlow.find('/th_ar') == -1  and fnlow.find('/thorium') == -1 and fnlow.find('/une') == -1) and \
                            not (fnlow.find('arc') == 0 or fnlow.find('thar') == 0 or fnlow.find('th_ar') == 0 or fnlow.find('thorium') == 0 or fnlow.find('une') == 0) and \
                            fiber1 not in ['rflat', 'cont', 'dark', 'bias', 'wave']:
                    fiber1 = 'science'
                    if fnlow.find('harps') != -1:
                        fiber2 = 'wave'
                im_head, obsdate_midexp, obsdate_mid_float, jd_midexp = get_obsdate(params, im_head)    # dateobs: unix_timestamp of mid exposure time
                extract = ''
                file_list.append([filename, fiber1, fiber2, im_head['HIERARCH HiFLEx EXPOSURE'], obsdate_mid_float, extract])
    
    return file_list

def check_raw_files_infos(params, file_list):
    """
    Reads the available information from the light in fiber1 and fiber2 to determine what calibration is possible with the data
    :param file_list: list of lists, with each entry containing: filename, fiber1, fiber2, exposure time, middle of observation time, string of extractions
    """
    # Check what data is available
    cal2_l_exp, cal2_s_exp, cal1_l_exp, cal1_s_exp = 0, 1E10, 0, 1E10
    arc_fib1 = ''
    sflat_fib2 = 'wave'
    blazecor_fib2 = 'none'
    exist_bias, exist_rflat, exp_darks = False, False, []
    for entry in file_list:
        if entry[1] == 'cont' and entry[2] == 'none' and blazecor_fib2 != 'dark':              # use entry[2] == 'wave' for cont only if entry[2] == 'none' not available
            sflat_fib2 = 'none'
        if entry[1] == 'cont' and entry[2] == 'dark':              # use entry[2] == 'dark' for cont only if entry[2] == 'none' not available
            sflat_fib2 = 'dark'
        if entry[1] == 'cont' and entry[2] == 'wave':               # use entry[2] == 'none' for arcflat only if entry[2] == 'wave' not available
            blazecor_fib2 = 'wave'
        if entry[1] == 'cont' and entry[2] == 'dark' and blazecor_fib2 != 'wave':   # use entry[2] == 'dark' for blazecor only if entry[2] == 'wave' not available
            blazecor_fib2 = 'dark'
        if (entry[1] == 'none' or entry[1] == 'dark' or entry[1] == 'wave') and entry[2] == 'wave':
            cal2_l_exp = max(entry[3], cal2_l_exp)
            cal2_s_exp = min(entry[3], cal2_s_exp)
        if (entry[2] == 'none' or entry[2] == 'dark' or entry[2] == 'wave') and entry[1] == 'wave':
            cal1_l_exp = max(entry[3], cal1_l_exp)
            cal1_s_exp = min(entry[3], cal1_s_exp)
        if entry[1] == 'bias' and entry[2] == 'bias':
            exist_bias = True
        if entry[1] == 'rflat' and entry[2] == 'rflat':
            exist_rflat = True
        if entry[1] == 'dark' and entry[2] == 'dark':
            if entry[3] not in exp_darks:
                exp_darks.append(entry[3])
    if cal2_l_exp == 0 and cal2_s_exp == 1E10:         # no single arcs taken, use cont and arc
        arc_fib1 = 'cont'
        for entry in file_list:                     # re-run to find exposure times
            if (entry[1] == arc_fib1) and entry[2] == 'wave':
                cal2_l_exp = max(entry[3], cal2_l_exp)
                cal2_s_exp = min(entry[3], cal2_s_exp)
    
    params['cal2_l_exp'] = cal2_l_exp
    params['cal2_s_exp'] = cal2_s_exp
    params['cal1_l_exp'] = cal1_l_exp
    params['cal1_s_exp'] = cal1_s_exp
    params['arc_fib1'] = arc_fib1
    params['sflat_fib2'] = sflat_fib2
    params['blazecor_fib2'] = blazecor_fib2
    params['exist_bias'] = exist_bias
    params['exist_rflat'] = exist_rflat
    params['exp_darks'] = exp_darks
    
    return params
    
def check_assigned_infos(params, file_list):
    """
    Same as check_raw_files_infos, but instead of reading what is in fiber1, fiber2 use the assigned extraction flags
    """
    cal2_l_exp, cal2_s_exp, cal1_l_exp, cal1_s_exp = 0, 1E10, 0, 1E10
    exist_bias, exist_rflat, exp_darks = False, False, []
    
    for entry in file_list:
        flags = entry[5].lower().replace(' ','').split(',')
        if 'w2' in flags:
            cal2_l_exp = max(entry[3], cal2_l_exp)
            cal2_s_exp = min(entry[3], cal2_s_exp)
        if 'w2l' in flags:
            cal2_l_exp = min(entry[3], cal2_l_exp)
        if 'w2s' in flags:
            cal2_s_exp = min(entry[3], cal2_s_exp)
        if 'w1' in flags:
            cal1_l_exp = max(entry[3], cal1_l_exp)
            cal1_s_exp = min(entry[3], cal1_s_exp)
        if 'w1l' in flags:
            cal1_l_exp = min(entry[3], cal1_l_exp)
        if 'w1s' in flags:
            cal1_s_exp = min(entry[3], cal1_s_exp)
        if 'b' in flags:
            exist_bias = True
        if 'a' in flags:
            exist_rflat = True
        if 'd' in flags:
            if entry[3] not in exp_darks:
                exp_darks.append(entry[3])

    params['cal2_l_exp'] = cal2_l_exp
    params['cal2_s_exp'] = cal2_s_exp
    params['cal1_l_exp'] = cal1_l_exp
    params['cal1_s_exp'] = cal1_s_exp
    params['exist_bias'] = exist_bias
    params['exist_rflat'] = exist_rflat
    params['exp_darks'] = exp_darks
    
    return params

def add_extraction_parameters_file_list(params, file_list, start_index):
    warn = []
    easy_assignments = []
    easy_assignments.append(['bias'  , 'bias'                 , 'b' ])       # Bias
    easy_assignments.append(['dark'  , 'dark'                 , 'd' ])       # Dark
    easy_assignments.append(['rflat' , 'rflat'                , 'a' ])       # real Flat
    easy_assignments.append(['cont'  , params['blazecor_fib2'] , 'z' ])       # blaze correction (continuum source)
    easy_assignments.append(['cont'  , params['sflat_fib2']   , 't1' ])      # trace of the science fiber
    for index, entry in enumerate(file_list):                                         # Extraction is set up below
        if index < start_index:                                     # Dont't touch the old files
            continue
        param = ''
        extract = ''
        for easy_assignment in easy_assignments:
            if entry[1] == easy_assignment[0] and entry[2] == easy_assignment[1]:                 # Fiber1 and Fiber2
                extract += easy_assignment[2]+', '
        if (entry[1] == 'none' or entry[1] == 'dark' or entry[1] == 'wave' or entry[1] == params['arc_fib1']) and entry[2] == 'wave':               # Fiber1 and Fiber2
            if params['cal2_l_exp']*0.7 <= entry[3] <= params['cal2_l_exp']*1.3:
                extract += 'w2l, '
                extract += 't2, '
            elif params['cal2_s_exp']*0.7 <= entry[3] <= params['cal2_s_exp']*1.3:                             # In case only one exposure time -> cal2_l is copied in cal2_s automatically
                extract += 'w2s, '
        if (entry[2] == 'none' or entry[2] == 'dark' or entry[2] == 'wave') and entry[1] == 'wave' and params['cal1_l_exp'] >= params['cal1_s_exp']:               # wave in science fiber
            parcal = 'arc'
            if params['cal1_l_exp']*0.7 <= entry[3] <= params['cal1_l_exp']*1.3:
                extract += 'w1l, '
            elif params['cal1_s_exp']*0.7 <= entry[3] <= params['cal1_s_exp']*1.3:                             # In case only one exposure time -> cal1_l is copied in cal1_s
                extract += 'w1s, '
        if entry[1] == 'science':
                extract += 'e, '
        if len(extract) >=2:
            extract = extract[:-2]
        file_list[index][5] = extract
    
    return file_list

def create_configuration_file(params, file_list):
    """
    Translates the file_list into fits_conf
    """
    exist_bias   = params['exist_bias']
    exist_rflat  = params['exist_rflat']
    exp_darks    = params['exp_darks']
    
    conf_data = dict()
    # Check if master files exist:
    if not exist_rflat:
        exist_rflat = os.path.isfile('master_rflat.fits')
        conf_data['master_rflat_filename'] = 'master_rflat.fits'
    if not exist_bias:
        exist_bias = os.path.isfile('master_bias.fits')
        conf_data['master_bias_filename'] = 'master_bias.fits'
    for file in os.listdir(params['result_path']):
        if file.endswith(".fits"):
            filename = os.path.join("", file)
            if filename.find('master_dark') == 0:
                im_head = fits.getheader(filename)
                im_head, obsdate_midexp, obsdate_mid_float, jd_midexp = get_obsdate(params, im_head)    # obsdate_mid_float: unix_timestamp of mid exposure time
                exptime = im_head['HIERARCH HiFLEx EXPOSURE']
                """exptime = filename.replace('master_dark','').replace('s.fits','').replace('p','.')
                try:
                    exptime = float( exptime )
                except:
                    continue"""
                if exptime not in exp_darks:
                    exp_darks.append(exptime)
                    conf_data['master_dark{}_filename'.format(exptime)] = filename
    
    # Create the configuartion file
    warn = []
    easy_assignments = []
    easy_assignments.append(['b'   , 'bias'    , 'bias'   ])       # Bias
    easy_assignments.append(['a'   , 'rflat'   , 'rflat'  ])       # real Flat
    easy_assignments.append(['z'   , 'blazecor' , 'blazecor'])       # blaze correction (continuum source)
    easy_assignments.append(['t1'  , 'trace1'  , 'trace1' ])       # trace of the science fiber
    easy_assignments.append(['t2'  , 'trace2'  , 'trace2' ])       # trace of the calibration fiber
    easy_assignments.append(['w1s' , 'cal1_s'  , 'arc'    ])       # Wavelength calibration (in science fiber)
    easy_assignments.append(['w1l' , 'cal1_l'  , 'arc'    ])       # Wavelength calibration (in science fiber)
    easy_assignments.append(['w2s' , 'cal2_s'  , 'arc'    ])       # Wavelength calibration (in calibration fiber)
    easy_assignments.append(['w2l' , 'cal2_l'  , 'arc'    ])       # Wavelength calibration (in calibration fiber)
    for entry in file_list:                                         # Extraction is set up below
        #if entry[0].find('#') != -1:        # comment # found
        #    continue
        #print entry, conf_data
        param = ''
        extract = entry[5].lower().replace(' ','').split(',')
        for easy_assignment in easy_assignments:
            if easy_assignment[0] in extract:
                conf_data, warn = create_parameters(conf_data, warn, easy_assignment[1], easy_assignment[1], [easy_assignment[2]], exist_bias, exist_rflat, exp_darks, entry)
        print(params['cal1_l_exp'], params['cal2_l_exp'], entry[3], 'w1' in extract, 'w2' in extract, entry )
        if 'w1' in extract:
            if params['cal1_l_exp']*0.7 <= entry[3] <= params['cal1_l_exp']*1.3:
                conf_data, warn = create_parameters(conf_data, warn, 'cal1_l', 'cal1_l', ['arc'], exist_bias, exist_rflat, exp_darks, entry)
            elif params['cal1_s_exp']*0.7 <= entry[3] <= params['cal1_s_exp']*1.3:                             # In case only one exposure time -> cal1_l is copied in cal1_s
                conf_data, warn = create_parameters(conf_data, warn, 'cal1_s', 'cal1_s', ['arc'], exist_bias, exist_rflat, exp_darks, entry)
        if 'w2' in extract:
            if params['cal2_l_exp']*0.7 <= entry[3] <= params['cal2_l_exp']*1.3:
                conf_data, warn = create_parameters(conf_data, warn, 'cal2_l', 'cal2_l', ['arc'], exist_bias, exist_rflat, exp_darks, entry)
            elif params['cal2_s_exp']*0.7 <= entry[3] <= params['cal2_s_exp']*1.3:                             # In case only one exposure time -> cal1_l is copied in cal1_s
                conf_data, warn = create_parameters(conf_data, warn, 'cal2_s', 'cal2_s', ['arc'], exist_bias, exist_rflat, exp_darks, entry)
        if 'd' in extract:               # Fiber1 and Fiber2
            param = 'dark{0}'.format(entry[3])                      # Exposure time
            textparam = param.replace('.','p')+'s'
            parcal = 'dark'
            conf_data, warn = create_parameters(conf_data, warn, param, textparam, [parcal], exist_bias, exist_rflat, exp_darks, entry)
        if entry[5].lower().find('e') == -1 and entry[5].lower().find('w') == -1:          # No extraction and no wavelength calibration
            continue
        for extraction in entry[5].replace(' ','').split(','):
            if extraction.lower() == 'wc':       # use this file for wavelength calibration between spectra
                param = 'waveoffsetcal'
                parcal = 'wavelengthcal'        # !!! Replace this by "waveoffset", but this also needs change in conf.txt: wavelengthcal_calibs_create_g to waveoffset_calibs_create_g
                conf_data, warn = create_parameters(conf_data, warn, param, param, [parcal], exist_bias, exist_rflat, exp_darks, entry)
                continue  
            elif extraction.lower() == 'ws':       # use this file for wavelength calibration between spectra
                param = 'waveoffsetsci'
                parcal = 'wavelengthcal'
                conf_data, warn = create_parameters(conf_data, warn, param, param, [parcal], exist_bias, exist_rflat, exp_darks, entry)
                continue                        # otherwise Warn from below will be triggered
            if extraction.lower().find('e') != 0:        # use find as len(extraction) might be 0
                # warn.append('Warn: I dont know what to do with the extraction parameter {0} (it doesnt begin with "e") for file {1}. This spectrum will therefore not be extracted. Please check {2} and run the script again, if you perform changes.'.format(extraction, entry[0], params['raw_data_file_list'] ))
                continue                                # The lines remaining in the for loop are only to be executed when it's about extraction, otherwise something might readded another time
            if extraction.lower().find('ec_') == 0:            # combine data before extraction
                param = 'extract_combine'+extraction[2:]
            elif extraction.lower().find('els_') == 0:
                param = 'extract_lin_sum'+extraction[3:]
            elif extraction.lower().find('elw_') == 0:
                param = 'extract_lin_weight'+extraction[3:]
            elif extraction.lower().find('elm_') == 0:
                param = 'extract_lin_med'+extraction[3:]
            elif extraction.lower().find('e') == 0:           # just extraction
                param = 'extract'+extraction[1:]
            conf_data, warn = create_parameters(conf_data, warn, param, param, [param, 'extract'], exist_bias, exist_rflat, exp_darks, entry)
    for fiber in [2,1]:             # Copy the long exposures into the short exposures, if not available and vice versa
        if 'cal{0}_l_rawfiles'.format(fiber) in conf_data.keys() and 'cal{0}_s_rawfiles'.format(fiber) not in conf_data.keys():
            conf_data['cal{0}_s_rawfiles'.format(fiber)]        = conf_data['cal{0}_l_rawfiles'.format(fiber)]
            conf_data['master_cal{0}_s_filename'.format(fiber)] = conf_data['master_cal{0}_l_filename'.format(fiber)].replace('cal{0}_l'.format(fiber),'cal{0}_s'.format(fiber))
            conf_data['cal{0}_s_calibs_create'.format(fiber)]   = conf_data['cal{0}_l_calibs_create'.format(fiber)]
        elif 'cal{0}_s_rawfiles'.format(fiber) in conf_data.keys() and 'cal{0}_l_rawfiles'.format(fiber) not in conf_data.keys():
            conf_data['cal{0}_l_rawfiles'.format(fiber)]        = conf_data['cal{0}_s_rawfiles'.format(fiber)]
            conf_data['master_cal{0}_l_filename'.format(fiber)] = conf_data['master_cal{0}_s_filename'.format(fiber)].replace('cal{0}_l'.format(fiber),'cal{0}_s'.format(fiber))
            conf_data['cal{0}_l_calibs_create'.format(fiber)]   = conf_data['cal{0}_s_calibs_create'.format(fiber)]

    for entry in warn:
        logger(entry)
    
    return conf_data

def file_list_UI(file_list):
    # Extract a common path
    maxlen_files = 0
    for entry in file_list:
        maxlen_files = max(maxlen_files, len(entry[0].replace('#','').replace(' ','')) )
    common_path = ''
    if file_list[0][0][2:].find(os.sep) != -1:         # it contains subfolders
        fnames = []
        for entry in file_list:
            fnames.append( entry[0].replace('#','').replace(' ','') )
        common_path = fnames[0].rsplit(os.sep,1)[0]    # just remove the filename
        for dummy in range(len(fnames[0].split(os.sep))-1):
            found_common_path = True
            for entry in fnames:
                if entry.find(common_path) == -1:   # This path is not the same
                    found_common_path = False
                    break
            if found_common_path:
                common_path += os.sep
                break
            else:
                common_path = common_path.rsplit(os.sep,1)     # remove the next subfolder
                if len(common_path) == 2:                   # Another subfolder in the path
                    common_path = common_path[0]
                else:                                       # no more subfolders
                    common_path = ''
                    break
    # Split the common_path to make the GUI shorter (in x), if necessary
    maxlen_files -= len(common_path)
    if len(common_path) == 0:
        common_path_text = '{0}\n(no common path)'.format(' '*int(maxlen_files*1.5))
    else:
        common_path_text = common_path
    if len(common_path) > maxlen_files:
        common_path_temp = common_path[:-1].split(os.sep)          # without the last '/'
        common_path_text = [common_path_temp[0]+os.sep]            # start with the first entry
        ii = 0
        for entry in common_path_temp[1:]:
            entry += os.sep
            if len(common_path_text[ii]+entry) > maxlen_files and len(common_path_text[ii]) > maxlen_files/2.:     # If too long to add
                common_path_text.append(entry)
                ii += 1
            else:
                common_path_text[ii] += entry
        common_path_text = '\n'.join(common_path_text)
    
    # define widgets
    pkwargs = dict()
    widgets = dict()
    widgets['comment'] = dict(label='Comm-\nented\nout ', kind='Label', row=0, column=0, columnspan=1, orientation=Tk.W)
    widgets['mid_exp'] = dict(label='   Observation Time    \n(mid exposure)\nUTC', kind='Label', row=0, column=1, orientation=Tk.W)
    widgets['exp']     = dict(label='Expo-\nsure \n[s]  ', kind='Label', row=0, column=2, orientation=Tk.W)
    widgets['name']    = dict(label='Path and folder\n{0}'.format(common_path_text), kind='Label', row=0, column=3, orientation=Tk.W)
    #widgets['fib1']    = dict(label='Science\nfiber', kind='Label', row=0, column=5)
    #widgets['fib2']    = dict(label='Calibration\nfiber', kind='Label', row=0, column=6)
    widgets['b']       = dict(label='Bias', kind='Label', row=0, column=6)
    widgets['d']       = dict(label='Dark', kind='Label', row=0, column=7)
    widgets['a']       = dict(label='Real\nFlat', kind='Label', row=0, column=8)
    widgets['t1']      = dict(label='Sci.\ntra-\nce', kind='Label', row=0, column=9)
    widgets['t2']      = dict(label='Cal.\ntra-\nce', kind='Label', row=0, column=10)
    widgets['z']       = dict(label='Blaze', kind='Label', row=0, column=11)
    widgets['w1']      = dict(label='Wave\nSci.', kind='Label', row=0, column=12)
    #widgets['w1l']     = dict(label='Wave\nSci.\nlong', kind='Label', row=0, column=12)
    #widgets['w1s']     = dict(label='Wave\nSci.\nshort', kind='Label', row=0, column=13)
    widgets['w2']     = dict(label='Wave\nCal.', kind='Label', row=0, column=14)
    #widgets['w2l']     = dict(label='Wave\nCal.\nlong', kind='Label', row=0, column=14)
    #widgets['w2s']     = dict(label='Wave\nCal.\nshort', kind='Label', row=0, column=15)
    widgets['ws']      = dict(label='Wave\nshft\nSci', kind='Label', row=0, column=16)
    widgets['wc']      = dict(label='Wave\nshft\nCal', kind='Label', row=0, column=17)
    widgets['e']       = dict(label='Ex-\ntract', kind='Label', row=0, column=18)
    widgets['extra']   = dict(label='Further usage\nof files\n(comma separated)', kind='Label', row=0, column=19)
    for ii, entry in enumerate(file_list):
        pkwargs['comment_{0}'.format(ii)] = ( 0 <= entry[0].find('#') < 20 )            # Commented out
        widgets['comment_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['comment_{0}'.format(ii)], row=ii+1, column=0)
        expstr = '%1.2f'%entry[3]       # '%9.2f'%entry[3] has not long enough spaces
        if len(expstr) > 7:
            expstr = '  '+expstr[:-3]        # without fractions at the end
        else:
            #expstr = ' '*{4:7, 5:5, 6:4, 7:2, 8:1, 9:1, 10:1}[len(expstr)] + expstr            # This has the 100s not perfectly alligned
            expstr = ' '*{4:6, 5:4, 6:2, 7:1, 8:1, 9:1, 10:1}[len(expstr)] + expstr             # This has the 100s and 100s not perfectly alligned
        text = '{0}{1}   {2}'.format( datetime.datetime.utcfromtimestamp(entry[4]).strftime('%Y-%m-%d %H:%M:%S'), 
                                     expstr, entry[0].replace('#','').replace(' ','').replace(common_path,'') )
        widgets['name_{0}'.format(ii)] = dict(label=text, kind='Label', row=ii+1, column=1, columnspan=3, orientation=Tk.W)
        #fname = entry[0].replace('#','').replace(' ','').replace(common_path,'')
        #widgets['mid_exp_{0}'.format(ii)] = dict(label=datetime.datetime.utcfromtimestamp(entry[4]).strftime('%Y-%m-%dT%H:%M:%S'), kind='Label', row=ii+1, column=5)
        #widgets['exp_{0}'.format(ii)]  = dict(label=entry[3], kind='Label', row=ii+1, column=4)
        #widgets['name_{0}'.format(ii)] = dict(label=fname, kind='Label', row=ii+1, column=1, columnspan=1, orientation=Tk.W)
        #pkwargs['fib1_{0}'.format(ii)] = entry[1]     # Allows modification of the fiber content
        #widgets['fib1_{0}'.format(ii)] = dict(kind='TextEntry', minval=None, maxval=None, fmt=str, start=pkwargs['fib1_{0}'.format(ii)], width=8, row=ii+1, column=5)     # Allows modification of the fiber content
        #pkwargs['fib2_{0}'.format(ii)] = entry[2]     # Allows modification of the fiber content
        #widgets['fib2_{0}'.format(ii)] = dict(kind='TextEntry', minval=None, maxval=None, fmt=str, start=pkwargs['fib2_{0}'.format(ii)], width=8, row=ii+1, column=6)     # Allows modification of the fiber content
        flags = entry[5].replace(' ','').split(',')
        pkwargs['b_{0}'.format(ii)] = ( 'b' in flags )
        widgets['b_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['b_{0}'.format(ii)], row=ii+1, column=6)
        pkwargs['d_{0}'.format(ii)] = ( 'd' in flags )
        widgets['d_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['d_{0}'.format(ii)], row=ii+1, column=7)
        pkwargs['a_{0}'.format(ii)] = ( 'a' in flags )
        widgets['a_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['a_{0}'.format(ii)], row=ii+1, column=8)
        pkwargs['t1_{0}'.format(ii)] = ( 't1' in flags )
        widgets['t1_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['t1_{0}'.format(ii)], row=ii+1, column=9)
        pkwargs['t2_{0}'.format(ii)] = ( 't2' in flags )
        widgets['t2_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['t2_{0}'.format(ii)], row=ii+1, column=10)
        pkwargs['z_{0}'.format(ii)] = ( 'z' in flags )
        widgets['z_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['z_{0}'.format(ii)], row=ii+1, column=11)
        pkwargs['w1_{0}'.format(ii)] = ( ('w1' in flags) | ('w1l' in flags) | ('w1s' in flags) )
        widgets['w1_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['w1_{0}'.format(ii)], row=ii+1, column=12)
        #pkwargs['w1l_{0}'.format(ii)] = ( 'w1l' in flags )
        #widgets['w1l_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['w1l_{0}'.format(ii)], row=ii+1, column=12)
        #pkwargs['w1s_{0}'.format(ii)] = ( 'w1s' in flags )
        #widgets['w1s_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['w1s_{0}'.format(ii)], row=ii+1, column=13)
        pkwargs['w2_{0}'.format(ii)] = ( ('w2' in flags) | ('w2l' in flags) | ('w2s' in flags) )
        widgets['w2_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['w2_{0}'.format(ii)], row=ii+1, column=14)
        #pkwargs['w2l_{0}'.format(ii)] = ( 'w2l' in flags )
        #widgets['w2l_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['w2l_{0}'.format(ii)], row=ii+1, column=14)
        #pkwargs['w2s_{0}'.format(ii)] = ( 'w2s' in flags )
        #widgets['w2s_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['w2s_{0}'.format(ii)], row=ii+1, column=15)
        pkwargs['ws_{0}'.format(ii)] = ( 'ws' in flags )
        widgets['ws_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['ws_{0}'.format(ii)], row=ii+1, column=16)
        pkwargs['wc_{0}'.format(ii)] = ( 'wc' in flags )
        widgets['wc_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['wc_{0}'.format(ii)], row=ii+1, column=17)
        pkwargs['e_{0}'.format(ii)] = ( 'e' in flags )
        widgets['e_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['e_{0}'.format(ii)], row=ii+1, column=18)
        extra = ''
        for flag in flags:          # Add extra flags, if necessary
            if flag not in ['b', 'd', 'a', 't1', 't2', 'z', 'w2l', 'w2s', 'w2', 'w1l', 'w1s', 'w1', 'ws', 'wc', 'e']:
                extra += ','+flag
        if len(extra) > 0:
            extra = extra[1:]
        pkwargs['extra_{0}'.format(ii)] = extra
        widgets['extra_{0}'.format(ii)] = dict(kind='TextEntry', minval=None, maxval=None, fmt=str, start=pkwargs['extra_{0}'.format(ii)], width=15, row=ii+1, column=19)
    
    explain = 'Explanation of the columns:\n'+\
              '- Tick first column to not use some files at all\n'+\
              '- Mark the files to be used for calibration:\n'+\
              '-- Bias: These file are combined into a master bias\n'+\
              '-- Dark: exposure time will be automatically taken\n   into account\n'+\
              '-- Real Flat: Evenly exposed detector to calibrate\n   pixel-to-pixel sensitivity variation\n'+\
              '-- Science trace: To trace the science orders\n'+\
              '-- Calibration trace: To trace the calibration orders\n'+\
              '-- Blaze: To derive the blaze function\n'+\
              '-- Wavelength solition for science fiber \n   (long and short expsoure time)\n'+\
              '-- Wavelength solution for calibration fiber (*)\n'+\
              '-- Wavelength offset Science fiber (**) to correct for\n   wavelength drit\n'+\
              '-- Wavelength offset between the Science fiber and\n   Calibration fiber (*)\n'+\
              '-- Extract: Extract these files on an individual basis\n'+\
              '-- Further settings (manual): e.g. to combine files\n   before or after extraction\n'+\
              '(*) not for single fiber spectrographs\n'+\
              '(**) important for unstabilised (single) fiber\n     spectrographs\n\n'+\
              'The automatic assignment is based on the parameters\n raw_data_* in {0} (and in procedure\n add_new_rawfiles_file_list). '.format(CONFIGFILE)
              #'- Type of Science and Calibration\n  fibers are derived from header or\n  filename and can be changed here\n  (optional)\n'+\     # Allows modification of the fiber content
    for ii, commentii in enumerate(explain.split('\n')):
        if len(commentii) > 0:
            widgets['explain_{0}'.format(ii)] = dict(label=commentii, kind='Label', row=ii, column=20, rowspan=1, orientation=Tk.W )#, wraplength=100 )      
    widgets['accept'] = dict(label='Accept', kind='ExitButton', row=ii+1, column=20, rowspan=2)
                              
    wprops = dict(fullscreen=False )
    #wprops['width_data'] = 800   # not neccssary, as uses the automatic width
    
    if len(file_list) > 100:
        logger('Info: A GUI with {0} elements will be created, on some machines that can take up to a few minutes.\n'.format(len(widgets))+\
               '\tIf you want to use an editor instead of the GUI, please kill the process (kill {0})\n'.format(os.getpid())+\
               '\tand run the script with parameter "nogui", e.g.\n\t\tpython {0} nogui'.format(sys.argv[0]))

    gui3 = tkc.TkCanvasGrid(title='HiFLEx: Asigning Files to extraction steps (What file contains what data)', 
                            kwargs=pkwargs,widgets=widgets, widgetprops=wprops )
    gui3.master.mainloop()
    
    # Get the information from the GUI
    file_list_new, file_list_full = [], []
    for ii in range(len(file_list)):
        text = ( {True:'b',False:''}[gui3.data['b_{0}'.format(ii)]] +','+
                 {True:'d',False:''}[gui3.data['d_{0}'.format(ii)]] +','+
                 {True:'a',False:''}[gui3.data['a_{0}'.format(ii)]] +','+
                 {True:'t1',False:''}[gui3.data['t1_{0}'.format(ii)]] +','+
                 {True:'t2',False:''}[gui3.data['t2_{0}'.format(ii)]] +','+
                 {True:'z',False:''}[gui3.data['z_{0}'.format(ii)]] +','+
                 {True:'w2',False:''}[gui3.data['w2_{0}'.format(ii)]] +','+
                 {True:'w1',False:''}[gui3.data['w1_{0}'.format(ii)]] +','+
                 {True:'ws',False:''}[gui3.data['ws_{0}'.format(ii)]] +','+
                 {True:'wc',False:''}[gui3.data['wc_{0}'.format(ii)]] +','+
                 {True:'e',False:''}[gui3.data['e_{0}'.format(ii)]] +','+
                 gui3.data['extra_{0}'.format(ii)]
               ).replace(',,',',').replace(',,',',').replace(',,',',').replace(',,',',').replace(',,',',')
        if text[0] == ',':
            text = text[1:]
        if len(text) > 0:
            if text[-1] == ',':
                text = text[:-1]
        #file_list_full.append([ {True:'#',False:''}[gui3.data['comment_{0}'.format(ii)]] + file_list[ii][0].replace('#','').replace(' ',''),
        #                       gui3.data['fib1_{0}'.format(ii)] , gui3.data['fib2_{0}'.format(ii)], file_list[ii][3], file_list[ii][4], 
        #                       text ])     # Allows modification of the fiber content
        file_list_full.append([ {True:'#',False:''}[gui3.data['comment_{0}'.format(ii)]] + file_list[ii][0].replace('#','').replace(' ',''),
                               file_list[ii][1] , file_list[ii][2], file_list[ii][3], file_list[ii][4], text ])
        if not gui3.data['comment_{0}'.format(ii)]:     # Ignore the commented out lines
            file_list_new.append(file_list_full[-1])
            
    return file_list_new, file_list_full

def get_observed_objects(params, conf_data):
    object_information_full = []
    object_information_head = []
    object_files, object_files_full = find_file_in_allfolders(params['object_file'], [params['result_path']] + params['raw_data_paths'])
    Simbad.add_votable_fields("pmra")  # Store proper motion in RA
    Simbad.add_votable_fields("pmdec")  # Store proper motion in Dec.
    for entry in sorted(conf_data):
        if entry.find('extract') != -1 and entry.find('rawfiles') != -1:     # Only check for extract*rawfiles
            list_to_search = conf_data[entry].split(', ')
        elif entry.find('master_extract_combine') != -1 and entry.find('filename') != -1:
            list_to_search = [entry.replace('master_extract_combine_','').replace('master_extract_combine','').replace('_filename','')]
        else:
            continue
        for im_name in list_to_search:
            found = False
            if os.path.isfile(im_name):
                im_head = fits.getheader(im_name)
                obname = im_name.replace('\n','').split(os.sep)    # get rid of the path
                obnames = get_possible_object_names(obname[-1], im_head, params['raw_data_object_name_keys'])
            else:           # not a raw file
                obnames = [im_name]
            # Check if already in the list
            for obname in obnames:
                for index, obname_found in enumerate(object_information_full):        # Same objectname doesn't need to be done again
                    if obname in [obname_found[0], obname_found[0].replace('_',''), obname_found[0].replace('-',''), obname_found[0].replace('_','').replace('-','')]:
                        if im_name not in object_information_full[index][-1]:
                            object_information_full[index][-1].append(im_name)
                        found = True
                        break
                if found:
                    break
            if found:
                continue
            # Search if the object exists in an availible file
            allentries_end = ['', '']
            for object_file in object_files:
                ra2, dec2, epoch, pmra, pmdec, obnames, allentries = getcoords_from_file(obnames, 0, filen=object_file, warn_notfound=False, ignore_values=True)        # mjd=0 because because not using ceres to calculated BCV
                if ra2 !=0 or dec2 != 0:                                           # Found the object
                    found = True
                    if len(allentries) > 0:
                        if len(allentries) >= 8:
                            allentries_end = allentries[6:8]
                    allentries.insert(1,'not searched')
                    break
            if not found:
                # Search the information on Simbad
                for obname in obnames:
                    obnames_space = [obname]
                    for i in range(1,len(obname)-1):    # Simbad requires the space in the object: e.g. Tau Ceti
                        obnames_space.append(obname[:i]+' '+obname[i:])
                    warnings.filterwarnings("ignore", category=UserWarning)
                    for ii in range(3):     # Try 3 times in case the connection is flaaky
                        try:
                            simbad_results = Simbad.query_objects(obnames_space)
                            break
                        except:
                            simbad_results = None
                            logger('Warn: Could not get the coordinates from Simbad for objects {0}. Will try another {1} time(s).'.format(obname, 2-ii))
                    warnings.filterwarnings("default", category=UserWarning)
                    if simbad_results is not None:
                        #print obname, simbad_results
                        obnames = [obname]
                        # Simbad quaries standard in J2000, this might change in future
                        allentries = [ obname, simbad_results[0]['MAIN_ID'], simbad_results[0]['RA'], simbad_results[0]['DEC'], simbad_results[0]['PMRA'], simbad_results[0]['PMDEC'], 1] + allentries_end + [2000]
                        found = True
                        break
            if not found:               # Create an empty list
                allentries = [obnames[0]] + ['']*6 + allentries_end + ['']
            allentries.append([im_name])
            # Get the information from the header
            im_head, obsdate_midexp, obsdate_mid_float, jd_midexp = get_obsdate(params, im_head)
            parms, source_radec, source_obs, mephem, obnames = get_object_site_from_header(params, im_head, obnames, obsdate_midexp)
            head_info = ['']*5
            if source_radec.find('from the image header') != -1:
                head_info = [ params['ra'] , params['dec'], params['pmra'], params['pmdec'], params['epoch'] ]
            #Put the information into a list
            object_information_full.append(allentries)
            object_information_head.append(head_info)
    return object_information_full, object_information_head

def calibration_parameters_coordinates_UI(conf_data, object_information_full, object_information_head):
    max_length = 1024
    def add_widgets(widgets, pkwargs, ii, ftype, doneftype):
        if '{0}_rawfiles'.format(ftype) not in conf_data.keys() or ftype in doneftype:
            return widgets, pkwargs, ii, doneftype
            
        widgets['type_{0}'.format(ii)] = dict(label=ftype, kind='Label', row=ii, column=0, columnspan=1, orientation=Tk.W)
        if 'master_{0}_filename'.format(ftype) in conf_data.keys() and (ftype.find('extract') == -1 or ftype.find('extract_combine') != -1) and ftype.find('wavelenoffset') == -1:
            elname = 'master_{0}_filename'.format(ftype)
            pkwargs[elname] = conf_data[elname]
            widgets[elname] = dict(kind='TextEntry', fmt=str, start=pkwargs[elname], width=30, row=ii, column=1, columnspan=3, orientation=Tk.W)
        elname = '{0}_calibs_create'.format(ftype)
        pkwargs[elname] = conf_data[elname]
        widgets[elname] = dict(kind='TextEntry', fmt=str, start=pkwargs[elname], width=50, row=ii, column=4, columnspan=5, orientation=Tk.W)
        label = conf_data['{0}_rawfiles'.format(ftype)]
        if len(label) > max_length:
            label = label[:max_length-3] + '...'        # Otherwise: X Error of failed request:  BadAlloc (insufficient resources for operation)
        widgets['files_{0}'.format(ii)]  = dict(label=label, kind='Label', row=ii, column=9, columnspan=1000, orientation=Tk.W)
        doneftype.append(ftype)
        ii += 1
        return widgets, pkwargs, ii, doneftype
    
    def vfunc_float(xs):
        try:
            value = float(xs)
            return True, value
        except:
            return False, ('Error, input must be float')
    
    # define widgets
    pkwargs = dict()
    widgets = dict()
    
    # Add widgets for fits_conf.txt
    widgets['type']        = dict(label='Type', kind='Label', row=0, column=0, columnspan=1, orientation=Tk.W)
    widgets['name_master'] = dict(label='Name of Master file', kind='Label', row=0, column=1, columnspan=3, orientation=Tk.W)
    widgets['calibs']      = dict(label='Calibrations to be applied', kind='Label', row=0, column=4, columnspan=5, orientation=Tk.W)
    widgets['files_used']  = dict(label='Files included (see last GUI or {0} to change the assigned files)'.format(params['raw_data_file_list']), kind='Label', row=0, column=9, columnspan=1000, orientation=Tk.W)
    darks = []
    for entry in sorted(params['exp_darks']):
        darks.append('dark{0}'.format(entry))
    doneftype = []
    ii = 1
    for ftype in ['bias'] + darks + ['rflat', 'trace1', 'trace2', 'blazecor', 'cal2_l', 'cal2_s', 'cal1_l', 'cal1_s', 'waveoffsetsci', 'waveoffsetcal']:
            widgets, pkwargs, ii, doneftype = add_widgets(widgets, pkwargs, ii, ftype, doneftype)
    for entry in sorted(conf_data.keys()):
        if entry.find('_rawfiles') >= 0:
            ftype = entry.replace('_rawfiles','')
            widgets, pkwargs, ii, doneftype = add_widgets(widgets, pkwargs, ii, ftype, doneftype)
    
    widgets['accept1'] = dict(label='Accept', kind='ExitButton', row=ii+1, column=6, rowspan=2, columnspan=2)
    ii += 3
    
    # Add widgets for object_list.txt
    widgets['name'] = dict(label='Object\nname', kind='Label', row=ii, column=0, rowspan=2, columnspan=1, orientation=Tk.E)
    widgets['header_info1'] = dict(label='Header information', kind='Label', row=ii, column=1, columnspan=2, orientation=Tk.E)
    widgets['header_ra'] = dict(label='RA', kind='Label', row=ii+1, column=1)
    widgets['header_dec'] = dict(label='DEC', kind='Label', row=ii+1, column=2)
    widgets['header_info2'] = dict(label='Header information', kind='Label', row=ii, column=3, columnspan=3)
    widgets['header_pmra'] = dict(label='PMRA', kind='Label', row=ii+1, column=3)
    widgets['header_pmdec'] = dict(label='PMRA', kind='Label', row=ii+1, column=4)
    widgets['header_epoch'] = dict(label='Epoch', kind='Label', row=ii+1, column=5)
    widgets['use_header'] = dict(label='Use\nhead', kind='Label', row=ii, column=6, rowspan=2, columnspan=1)
    widgets['ra'] = dict(label='RA\n<:>,< >,float', kind='Label', row=ii, column=7, rowspan=2, columnspan=1)
    widgets['dec'] = dict(label='DEC\n<:>,< >,float', kind='Label', row=ii, column=8, rowspan=2, columnspan=1)
    widgets['pmra'] = dict(label='PMRA\n[mas/yr]', kind='Label', row=ii, column=9, rowspan=2, columnspan=1)
    widgets['pmdec'] = dict(label='PMDEC\n[mas/yr]', kind='Label', row=ii, column=10, rowspan=2, columnspan=1)
    widgets['epoch'] = dict(label='Epoch\n(number)', kind='Label', row=ii, column=11, rowspan=2, columnspan=1)
    widgets['simbad_name'] = dict(label='Simbad Name\n(check that correct)', kind='Label', row=ii, column=12, rowspan=2, columnspan=1)
    widgets['mask'] = dict(label='Mask (optional)\nG2,K5,M2', kind='Label', row=ii, column=13, rowspan=2, columnspan=1)
    widgets['rot'] = dict(label='rotation (optional)\n[km/s]', kind='Label', row=ii, column=14, rowspan=2, columnspan=1)
    ii += 2
    jj = 0                              # if no object to extract, jj will not exist, but is needed later
    for jj, entry in enumerate(object_information_full):
        widgets['name_{0}'.format(jj)]  = dict(label=entry[0], kind='Label', row=ii+jj, column=0, columnspan=1, orientation=Tk.E)
        state = None
        if entry[0].lower().find('sun') == 0 or entry[0].lower().find('moon') == 0 or entry[0].lower().find('jupiter') == 0:     # Check with get_object_site_from_header()
            for i in range(1,7):
                entry[i] = ''
            state = Tk.DISABLED
        pkwargs['ra_{0}'.format(jj)]    = entry[2]
        widgets['ra_{0}'.format(jj)]    = dict(kind='TextEntry', fmt=str, start=pkwargs['ra_{0}'.format(jj)], width=13, row=ii+jj, column=7, columnspan=1, state=state)
        pkwargs['dec_{0}'.format(jj)]   = entry[3]
        widgets['dec_{0}'.format(jj)]   = dict(kind='TextEntry', fmt=str, start=pkwargs['dec_{0}'.format(jj)], width=13, row=ii+jj, column=8, columnspan=1, state=state)
        pkwargs['pmra_{0}'.format(jj)]  = entry[4]
        widgets['pmra_{0}'.format(jj)]  = dict(kind='TextEntry', fmt=str, valid_function=vfunc_float, start=pkwargs['pmra_{0}'.format(jj)], width=7, row=ii+jj, column=9, columnspan=1, state=state)
        pkwargs['pmdec_{0}'.format(jj)] = entry[5]
        widgets['pmdec_{0}'.format(jj)] = dict(kind='TextEntry', fmt=str, valid_function=vfunc_float, start=pkwargs['pmdec_{0}'.format(jj)], width=7, row=ii+jj, column=10, columnspan=1, state=state)
        pkwargs['epoch_{0}'.format(jj)] = entry[9]
        widgets['epoch_{0}'.format(jj)] = dict(kind='TextEntry', fmt=str, valid_function=vfunc_float, start=pkwargs['epoch_{0}'.format(jj)], width=6, row=ii+jj, column=11, columnspan=1, state=state)
        # simbad name at end, next to button (button to search again is missing, therefore state=Tk.DISABLED
        pkwargs['simbad_name_{0}'.format(jj)] = entry[1]
        widgets['simbad_name_{0}'.format(jj)] = dict(kind='TextEntry', fmt=str, start=pkwargs['simbad_name_{0}'.format(jj)], width=13, row=ii+jj, column=12, columnspan=1, orientation=Tk.W, state=Tk.DISABLED)

        pkwargs['mask_{0}'.format(jj)]  = entry[7]
        widgets['mask_{0}'.format(jj)]  = dict(kind='TextEntry', fmt=str, start=pkwargs['mask_{0}'.format(jj)], width=5, row=ii+jj, column=13, columnspan=1)
        pkwargs['rot_{0}'.format(jj)]   = entry[8]
        widgets['rot_{0}'.format(jj)]   = dict(kind='TextEntry', fmt=str, valid_function=vfunc_float, start=pkwargs['rot_{0}'.format(jj)], width=5, row=ii+jj, column=14, columnspan=1)
        
        if object_information_head[jj][0] == '':
            pkwargs['use_header_{0}'.format(jj)]   = False
            widgets['use_header_{0}'.format(jj)]   = dict(label=None,  kind='CheckBox', start=pkwargs['use_header_{0}'.format(jj)], row=ii+jj, column=6, state=Tk.DISABLED)
            continue
        entry = object_information_head[jj]
        widgets['header_ra_{0}'.format(jj)]    = dict(label=entry[0], kind='Label', row=ii+jj, column=1)
        widgets['header_dec_{0}'.format(jj)]   = dict(label=entry[1], kind='Label', row=ii+jj, column=2)
        widgets['header_pmra_{0}'.format(jj)]  = dict(label=entry[2], kind='Label', row=ii+jj, column=3)
        widgets['header_pmdec_{0}'.format(jj)] = dict(label=entry[3], kind='Label', row=ii+jj, column=4)
        widgets['header_epoch_{0}'.format(jj)] = dict(label=entry[4], kind='Label', row=ii+jj, column=5)
        pkwargs['use_header_{0}'.format(jj)]   = True
        widgets['use_header_{0}'.format(jj)]   = dict(label=None,  kind='CheckBox', start=pkwargs['use_header_{0}'.format(jj)], row=ii+jj, column=6)
    
    ii += jj+1
    widgets['accept2'] = dict(label='Accept', kind='ExitButton', row=ii, column=6, rowspan=2, columnspan=2)
    ii += 2
    
    explain = 'Explanation of upper half:\n'+\
              '  This assigns the calibration that will be applied to files of the different types before creating the master file or before extracting the spectra. The following comma separated options are possible:\n'+\
              '     subframe, badpx_mask, bias, dark, flat, background, normalise, combine_mean, combine_sum\n'+\
              '     Please check the manual for more information on these options.\n'+\
              '  The assigned calibration steps are read from {0} (*_calibs_create_g)\n'.format(CONFIGFILE)+\
              '  The information from this part of the GUI will be stored in {0}.\n'.format(params['configfile_fitsfiles'])+\
              '\nExplanation of the lower half:\n'+\
              '  For each object to be extracted the coordinates are derived/displayed here. If the header information is available then the information is shown and the user can decide if this information should be used.\n'+\
              '  The editable coordinates are taken from the file {0}, for which is also checked in the result and raw data path. If the object does not exist in the file, then Simbad is searched using the Object name.\n'.format(params['object_file'])+\
              '  The object name for the results from Simbad is shown and should be the same as the object name derived from header/filename\n'+\
              '  The user can modify this information, which should be correct to perform correct barycentric correction. The information is then stored in {0}, overwriting the previous information.\n'.format(params['object_file'])+\
              '     The RA and DEC can be given in hour or degree format. The optional parameters "mask" and "rotation speed" are used if RV analysis with the CERES pipeline is performed.'
    for jj, commentjj in enumerate(explain.split('\n')):
        if len(commentjj) > 0:
            widgets['explain_{0}'.format(jj)] = dict(label=commentjj, kind='Label', row=ii+jj, column=0, columnspan=20, orientation=Tk.W )#, wraplength=100 )      
    
    wprops = dict(fullscreen=False )
    #wprops['width_GUI'] = 800   # not neccssary, as uses the automatic width
    
    gui3 = tkc.TkCanvasGrid(title='HiFLEx: Asigning Calibration and Coordinates', 
                            kwargs=pkwargs,widgets=widgets, widgetprops=wprops )
    gui3.master.mainloop()
    
    # Get the information from the GUI
    for entry in conf_data.keys():      # For the conf data
        if entry in gui3.data.keys():
            conf_data[entry] = gui3.data[entry]
            if entry.find('_calibs_create') != -1:
                if len(conf_data[entry]) > 0:
                    if conf_data[entry].find('[') + conf_data[entry].find(']') == -2:
                        conf_data[entry] = '[{0}]'.format(conf_data[entry])                     # add brackets
    
    object_information = []
    for ii in range(len(object_information_full)):
        enabled = {True:0,False:1}[gui3.data['use_header_{0}'.format(ii)]]
        if len(object_information_full[ii][1]) >= 3:            # The coordinates were read from a file
            try:
                object_information_full[ii][6] = float(object_information_full[ii][6])      # see if we can make it to a number
            except:
                True
            if type(object_information_full[ii][6]).__name__ != 'str':          # it's a number
                if abs(object_information_full[ii][6]) < 0.9:
                    enabled = False                                 # Object was disabled before
        if ( len(gui3.data['mask_{0}'.format(ii)]) < 1 and type(gui3.data['rot_{0}'.format(ii)]).__name__ == 'str' ) and ( len(gui3.data['ra_{0}'.format(ii)]) < 3 or len(gui3.data['ra_{0}'.format(ii)]) < 3 ):
            continue                                            # Both coordinates are not complete and Mask or rotation is not there -> no useful information
        if (len(gui3.data['ra_{0}'.format(ii)]) >= 3 and len(gui3.data['dec_{0}'.format(ii)]) >= 3):     # Fill pm* and epoch if coordinates are full
            for entry in ['pmra', 'pmdec']:
                if type(gui3.data['{1}_{0}'.format(ii, entry)]).__name__ == 'str':
                    gui3.data['{1}_{0}'.format(ii, entry)] = 0
            if type(gui3.data['epoch_{0}'.format(ii)]).__name__ == 'str':
                gui3.data['epoch_{0}'.format(ii)] = 2000
        object_information.append([ object_information_full[ii][0], gui3.data['ra_{0}'.format(ii)], gui3.data['dec_{0}'.format(ii)], 
                                    gui3.data['pmra_{0}'.format(ii)], gui3.data['pmdec_{0}'.format(ii)], gui3.data['epoch_{0}'.format(ii)],
                                    enabled, gui3.data['mask_{0}'.format(ii)], gui3.data['rot_{0}'.format(ii)] ])
    
    return conf_data, object_information

def comment_out_nonexisting_rawfiles_file_list(params, file_list):
    """
    Comments out the file which are in raw_data_file_list, but can't be found anymore
    """
    commented_out = 0
    for ii in range(len(file_list)):
        if file_list[ii][0].find('#') == -1:
            if not os.path.isfile(file_list[ii][0]):
                file_list[ii][0] = '#' + file_list[ii][0]       # files doesn't exist -> comment out
                commented_out += 1
    if commented_out > 0:
        logger('Warn: {0} Files in {1} could not be found and have been commented out.'.format(commented_out, params['raw_data_file_list']))
    
    return file_list

if __name__ == "__main__":
    logger('Info: Preparing a list with the files')
    log_params(params)
    
    # get the available list
    file_list = read_text_file(params['raw_data_file_list'], no_empty_lines=True)
    try:
        #file_list = convert_readfile(file_list, [str, str, str, float, ['%Y-%m-%dT%H:%M:%S', float], str], delimiter='\t', replaces=['\n',' '], ignorelines=[['#',20]])     #new way of storing the data
        file_list = convert_readfile(file_list, [str, str, str, float, ['%Y-%m-%dT%H:%M:%S', float], str], delimiter='\t', replaces=['\n',' '])     #new way of storing the data
    except:
        file_list = convert_readfile(file_list, [str, str, str, float, float, str], delimiter='\t', replaces=['\n',' ']) # old way of reading the data, To stay backwards compatible, can be removed in a few versions after v0.4.1
    file_list = comment_out_nonexisting_rawfiles_file_list(params, file_list)
    number_old_entries = len(file_list)
    
    # get the new files
    file_list = add_new_rawfiles_file_list(params, file_list)           # Check for new entries for file_list
    
    # analyse the global information
    params = check_raw_files_infos(params, file_list)
    
    # add the configuration to extraction parameters
    file_list = add_extraction_parameters_file_list(params, file_list, number_old_entries)
    if len(file_list) == 0:
        logger('Error: No files found. Please check that the folder(s) {0} contain files ending with {1}'.format( params['raw_data_paths'], params['raw_data_file_endings'] ))

    file_list = sorted(file_list, key=operator.itemgetter(4,0))           #itemgetter(1,2,3,4,0)
    
    # Show the results in a GUI
    if 'nogui' in sys.argv or '-nogui' in sys.argv or '--nogui' in sys.argv:
        file_list_commented = file_list
    else:
        file_list, file_list_commented = file_list_UI(file_list)
    
    # Save the list, show to user, so the user can disable files, read the list
    with open(params['raw_data_file_list'], 'w') as file:
        file.write('### This file contains the information for all the fits files in the raw_data_paths: {0} and its/their subfolders.\n'.format(params['raw_data_paths']))
        file.write('### Each line contains the following information, separated by one tab:\n')
        file.write('###   - Fill filename \n')
        file.write('###   - Type of fiber 1 (science fiber)\n')
        file.write('###   - Type of fiber 2 (calibration fiber)\n')
        file.write('###   - Exposure time in seconds (from the header, if the information is not in the header, then from the filename)\n')
        file.write('###   - Mid-exposure observation time in UTC (from header)\n')
        file.write('###   - Flags: Mark the files to be used for calibration the data (comma-separated list):\n')
        file.write('###        "b", Bias.\n')
        file.write('###        "d", Dark.\n')
        file.write('###        "a", real Flat of the detector.\n')
        file.write('###        "z", Spectrum for the blaze correction, e.g. of a continuum source.\n')
        file.write('###        "t1", Spectrum to find the trace [of the science fiber], e.g. a continuum source.\n')
        file.write('###        "t2", Spectrum to find the trace of the calibration fiber.\n')
        #file.write('###        "w1l, w1s", Spectruum to find the wavelength solution of the science fiber (long and short exposure time).\n')
        #file.write('###        "w2l, w2s", Spectruum to find the wavelength solution of the calibration fiber.\n')
        file.write('###        "w1", Spectruum to find the wavelength solution of the science fiber (long and short exposure time).\n')
        file.write('###        "w2", Spectruum to find the wavelength solution of the calibration fiber.\n')
        file.write('###        "e", if the spectra of this file should be extraced. By standard only the science data is extracted.\n')
        file.write('###             Combination of images before the extraction is possible, please refer to the manual for more information\n')
        file.write('###        "ws" or "wc", if the spectrum contains wavelength information to calculate the offset to the wavelength solution and the offset between both fibers of a bifurcated fiber).\n')
        file.write('###             Used to mark the files with calibration spectrum taken before or after the science observation.\n')
        file.write('### To exlude the use of some files please comment the line with a "#" or delete the line. \n\n')
        file.write('### -> When finished with the editing, save the file and close the editor \n\n')
        for entry in file_list_commented:
            file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(entry[0].ljust(50), entry[1], entry[2], entry[3], datetime.datetime.utcfromtimestamp(entry[4]).strftime('%Y-%m-%dT%H:%M:%S'), entry[5] ))
    
    # If necessary show the text file instead of the GUI
    if ('nogui' in sys.argv or '-nogui' in sys.argv or '--nogui' in sys.argv):
        start = time.time()
        rtn = os.system('{1} {0}'.format(params['raw_data_file_list'], params['editor'] ))
        if rtn != 0 or time.time()-start < 10:
            print('Please check that file {0} is correct.'.format(params['raw_data_file_list']))
            raw_input('To continue please press Enter\t\t')
        file_list = read_text_file(params['raw_data_file_list'], no_empty_lines=True)
        file_list = convert_readfile(file_list, [str, str, str, float, ['%Y-%m-%dT%H:%M:%S', float], str], delimiter='\t', replaces=['\n',' '], ignorelines=[['#',20]])
        file_list = sorted(file_list, key=operator.itemgetter(4,0))       # itemgetter(1,2,3,0)
    
    # Reset the list of parameters in case important data was deleted, e.g. all Darks
    del params['cal2_l_exp']
    del params['cal2_s_exp']
    del params['cal1_l_exp']
    del params['cal1_s_exp']
    del params['arc_fib1']
    del params['sflat_fib2']
    del params['blazecor_fib2']
    del params['exist_bias']
    del params['exist_rflat']
    del params['exp_darks']
    params = check_assigned_infos(params, file_list)
    
    # Connect wavelength solutions to the correct exposure time
    # -> done in create_configuration_file
    
    # Create the data for fits_conf.txt
    conf_data = create_configuration_file(params, file_list)
    
    # What objects were observed:
    object_information_full, object_information_head = get_observed_objects(params, conf_data)
    
    # Select the calibration parameters and Object coordinates in a GUI
    conf_data, object_information = calibration_parameters_coordinates_UI(conf_data, object_information_full, object_information_head)

    # Append information to params['object_file']
    lines_txt = read_text_file(params['object_file'], no_empty_lines=True, warn_missing_file=False)
    object_names = convert_readfile(lines_txt, [str], delimiter=',', replaces=[['\t',',']], expand_input=True)
    new_obj_names = []
    if len(object_information) > 0:
        with open(params['object_file'],'w') as file:
            for entry in object_information:
                file.write('{0},{1},{2},{3},{4},{6},{7},{8},{5}\n'.format( entry[0], entry[1], entry[2], entry[3], entry[4], entry[5], entry[6], entry[7], entry[8] ))  # Epoch (5) at end to be compatible with CERES
                new_obj_names.append(entry[0])
            for ii in range(len(object_names)):
                if object_names[ii][0] not in new_obj_names:
                    file.write(lines_txt[ii]+'\n')
                    new_obj_names.append(object_names[ii][0])
    
    # Save the results in a conf_data.txt file
    #print json.dumps(conf_data,sort_keys = False, indent = 4)
    with open(params['configfile_fitsfiles'], 'w') as file:
        file.write('# Description of this file:\n#--------------------------\n') 
        file.write('# This file was created by combining the information in the file {0} and the parameters given in the configuration file {1}\n'.format(params['raw_data_file_list'], CONFIGFILE))
        file.write('# Changes in here will be overwritten the next time prepare_file_list.py is run, therefore we suggest to make changes in {0} and afterwards run prepare_file_list.py again.\n\n'.format(params['raw_data_file_list']))
        file.write('# For each type of calibration filetype (e.g. bias, darks) a <filetype>_calibs_create list all corrections that will be applied to the individual files listed in <filetype>_rawfiles. These files will then combined and stored in master_<filetype>_filename (remove/comment the parameter in order to avoid saving.\n') 
        file.write('#   The following calibrations are possible (case-insensitive): subframe, badpx_mask, bias, dark, flat, <background>, normalisation, combine_sum, localbackground\n') 
        file.write('#   For dark and flat the paramerters can contain the the exposure time in float format (e.g. flat15.0_rawfiles, dark4.0_calibs_create).\n')  
        file.write('#       If a fixed dark should be used, than the parameter needs to contain "dark" and additional text (e.g. "darkfixed", or if a different exposure time should be used "dark5.0")\n') 
        file.write('#   For <background> the calibration needs to contain "background" but can contain more information. The key needs to be defined\n') 
        file.write('#        (e.g. if the "background_image_filename" is used for calibration then the following entry is needed as well here or in the calibration file "background_image_filename = background.fits"\n') 
        file.write('# \n') 
        for paramtype in ['rawfiles', 'calibs_create', 'master']:
            file.write('\n')
            for entry in sorted(conf_data.keys()):
                if entry.find(paramtype) >= 0:
                    file.write('{0} = {1} \n'.format(entry.ljust(24), conf_data[entry]))
    
    logger('Info: The calibration file for handling the raw data has been created. Please check {0} before starting the data reduction (hiflex.py)'.format(params['configfile_fitsfiles']))
    
    

