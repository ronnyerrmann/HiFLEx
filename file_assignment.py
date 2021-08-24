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
    :param entry: list of [str, str, str, float, datetime, str, str]: entry as in params['raw_data_file_list']
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
                    warn_text = 'Warn: The parameter {0}_calibs_create_g in the configuration file requires {1} correction, but a dark with exposure time of {2}s is not available. {4}\t\tPlease change the dark parameter to a fixed exposure time in {3}'.format(parcal, param_entry, entry[3], params['configfile_fitsfiles'], os.linesep )
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
                obname = filename.replace('\n','').replace(os.linesep,'').split(os.sep)    # get rid of the path
                obnames = get_possible_object_names(obname[-1], im_head, params['raw_data_object_name_keys'])
                file_list.append([filename, fiber1, fiber2, im_head['HIERARCH HiFLEx EXPOSURE'], obsdate_mid_float, obnames[0], extract])
    
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
        flags = entry[6].lower().replace(' ','').split(',')
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
        if len(extract) >=2:    # remove the ', ' at the end
            extract = extract[:-2]
        file_list[index][6] = extract
    
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
        extract = entry[6].lower().replace(' ','').split(',')
        for easy_assignment in easy_assignments:
            if easy_assignment[0] in extract:
                conf_data, warn = create_parameters(conf_data, warn, easy_assignment[1], easy_assignment[1], [easy_assignment[2]], exist_bias, exist_rflat, exp_darks, entry)
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
        if entry[6].lower().find('e') == -1 and entry[6].lower().find('w') == -1:          # No extraction and no wavelength calibration
            continue
        for extraction in entry[6].replace(' ','').split(','):
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

def get_observed_objects(params, conf_data, file_list):
    """
    
    """
    file_list =copy.deepcopy(file_list)
    object_information_full = []
    object_information_head = []
    object_files, object_files_full = find_file_in_allfolders(params['object_file'], [params['result_path']] + params['raw_data_paths'])
    Simbad.add_votable_fields("pmra")  # Store proper motion in RA
    Simbad.add_votable_fields("pmdec")  # Store proper motion in Dec.
    for entryc in sorted(conf_data):
        if entryc.find('extract') != -1 and entryc.find('rawfiles') != -1:     # Only check for extract*rawfiles
            list_to_search = conf_data[entryc].split(', ')
        elif entryc.find('master_extract_combine') != -1 and entryc.find('filename') != -1:
            list_to_search = [entryc.replace('master_extract_combine_','').replace('master_extract_combine','').replace('_filename','')]
        else:
            continue
        for im_name in list_to_search:
            if os.path.isfile(im_name):
                im_head = fits.getheader(im_name)
            else:
                im_head = dict()
            found = False
            # Look in the other list for the information about the object - this process will take a while
            for ii, entry in enumerate(file_list):
                if entry[0] == im_name:
                    if entry[5] != '':
                        obnames = entry[5].replace(' ','').split(',')
                        found = True
                    del file_list[ii]
                    break
            # Get the object name from the filename
            if not found:
                if os.path.isfile(im_name):
                    obname = im_name.replace('\n','').replace(os.linesep,'').split(os.sep)    # get rid of the path
                    obnames = get_possible_object_names(obname[-1], im_head, params['raw_data_object_name_keys'])
                else:           # not a raw file
                    obnames = [im_name]
            # Check if already in the list
            found = False
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

def comment_out_nonexisting_rawfiles_file_list(params, file_list):
    """
    Comments out the file which are in raw_data_file_list, but can't be found anymore
    """
    commented_out = 0
    for ii in range(len(file_list)):
        if file_list[ii][0].find('#') == -1:                    # Not already commended out
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
    
    # Only necessary to be backwards compatible to versions 1.5.0 and older
    for ii in range(len(file_list)):
        number_entries = file_list[ii].count('\t')
        if number_entries == 5:
            line_temp = file_list[ii].rsplit('\t',1)        # Split the old extraction
            file_list[ii] = line_temp[0] + '\t\t' + line_temp[-1]   # Add another tab at the place
    
    # Back to normal        
    file_list = convert_readfile(file_list, [str, str, str, float, ['%Y-%m-%dT%H:%M:%S', float], str, str], delimiter='\t', replaces=['\n',' ', os.linesep])     # will fail for lists created before version v0.4.0
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
        file_list, file_list_commented = proc_gui.file_list_UI(file_list, CONFIGFILE)
    
    # Save the list, show to user, so the user can disable files, read the list
    with open(params['raw_data_file_list'], 'w') as file:
        file.write('### This file contains the information for all the fits files in the raw_data_paths: {0} and its/their subfolders.{1}'.format(params['raw_data_paths'],os.linesep))
        file.write('### Each line contains the following information, separated by one tab:'+os.linesep)
        file.write('###   - Fill filename'+os.linesep)
        file.write('###   - Type of fiber 1 (science fiber)'+os.linesep)
        file.write('###   - Type of fiber 2 (calibration fiber)'+os.linesep)
        file.write('###   - Exposure time in seconds (from the header, if the information is not in the header, then from the filename)'+os.linesep)
        file.write('###   - Mid-exposure observation time in UTC (from header)'+os.linesep)
        file.write('###   - Object name (several can be given comma-separated)'+os.linesep)
        file.write('###   - Flags: Mark the files to be used for calibration the data (comma-separated list):'+os.linesep)
        file.write('###        "b", Bias.'+os.linesep)
        file.write('###        "d", Dark.'+os.linesep)
        file.write('###        "a", real Flat of the detector.'+os.linesep)
        file.write('###        "z", Spectrum for the blaze correction, e.g. of a continuum source.'+os.linesep)
        file.write('###        "t1", Spectrum to find the trace [of the science fiber], e.g. a continuum source.'+os.linesep)
        file.write('###        "t2", Spectrum to find the trace of the calibration fiber.'+os.linesep)
        #file.write('###        "w1l, w1s", Spectruum to find the wavelength solution of the science fiber (long and short exposure time).\n')
        #file.write('###        "w2l, w2s", Spectruum to find the wavelength solution of the calibration fiber.\n')
        file.write('###        "w1", Spectruum to find the wavelength solution of the science fiber (long and short exposure time).'+os.linesep)
        file.write('###        "w2", Spectruum to find the wavelength solution of the calibration fiber.'+os.linesep)
        file.write('###        "ws" or "wc", if the spectrum contains wavelength information to calculate the offset to the wavelength solution and the offset between both fibers of a bifurcated fiber).'+os.linesep)
        file.write('###             Used to mark the files with calibration spectrum taken before or after the science observation.'+os.linesep)
        file.write('###        "e", if the spectra of this file should be extraced. By standard only the science data is extracted.'+os.linesep)
        file.write('###        "ec_<obj>", combine files before extration.'+os.linesep)
        file.write('###        "elw_<obj>" or "elm_<obj>" or "els_<obj>", combine the extracted and linearised spectra using weighted average or median or sum.'+os.linesep)
        file.write('###             Please refer to the manual for more information.'+os.linesep)
        file.write('### To exlude the use of some files please comment the line with a "#" or delete the line.'+os.linesep+os.linesep)
        file.write('### -> When finished with the editing, save the file and close the editor'+os.linesep+os.linesep)
        for entry in file_list_commented:
            file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}{7}'.format(entry[0].ljust(50), entry[1], entry[2], entry[3], datetime.datetime.utcfromtimestamp(entry[4]).strftime('%Y-%m-%dT%H:%M:%S'), entry[5], entry[6], os.linesep ))
    
    # If necessary show the text file instead of the GUI
    if ('nogui' in sys.argv or '-nogui' in sys.argv or '--nogui' in sys.argv):
        start = time.time()
        rtn = os.system('{1} {0}'.format(params['raw_data_file_list'], params['editor'] ))
        if rtn != 0 or time.time()-start < 10:
            print('Please check that file {0} is correct.'.format(params['raw_data_file_list']))
            raw_input('To continue please press Enter\t\t')
        file_list = read_text_file(params['raw_data_file_list'], no_empty_lines=True)
        file_list = convert_readfile(file_list, [str, str, str, float, ['%Y-%m-%dT%H:%M:%S', float], str, str], delimiter='\t', replaces=['\n',os.linesep,' '], ignorelines=[['#',20]])
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
    object_information_full, object_information_head = get_observed_objects(params, conf_data, file_list)
    
    # Select the calibration parameters and Object coordinates in a GUI
    conf_data, object_information = proc_gui.calibration_parameters_coordinates_UI(params, conf_data, object_information_full, object_information_head, CONFIGFILE)

    # Append information to params['object_file']
    lines_txt = read_text_file(params['object_file'], no_empty_lines=True, warn_missing_file=False)
    object_names = convert_readfile(lines_txt, [str], delimiter=',', replaces=[['\t',',']], expand_input=True)
    new_obj_names = []
    if len(object_information) > 0:
        with open(params['object_file'],'w') as file:
            for entry in object_information:
                file.write('{0},{1},{2},{3},{4},{6},{7},{8},{5}{9}'.format( entry[0], entry[1], entry[2], entry[3], entry[4], entry[5], entry[6], entry[7], entry[8], os.linesep ))  # Epoch (5) at end to be compatible with CERES
                new_obj_names.append(entry[0])
            for ii in range(len(object_names)):
                if object_names[ii][0] not in new_obj_names:
                    file.write(lines_txt[ii]+os.linesep)
                    new_obj_names.append(object_names[ii][0])
    
    # Save the results in a conf_data.txt file
    #print json.dumps(conf_data,sort_keys = False, indent = 4)
    with open(params['configfile_fitsfiles'], 'w') as file:
        file.write('# Description of this file:'+os.linesep+'#--------------------------'+os.linesep) 
        file.write('# This file was created by combining the information in the file {0} and the parameters given in the configuration file {1}{2}'.format(params['raw_data_file_list'], CONFIGFILE, os.linesep))
        file.write('# Changes in here will be overwritten the next time prepare_file_list.py is run, therefore we suggest to make changes in {0} and afterwards run prepare_file_list.py again.{1}{1}'.format(params['raw_data_file_list'], os.linesep))
        file.write('# For each type of calibration filetype (e.g. bias, darks) a <filetype>_calibs_create list all corrections that will be applied to the individual files listed in <filetype>_rawfiles. These files will then combined and stored in master_<filetype>_filename (remove/comment the parameter in order to avoid saving.'+os.linesep) 
        file.write('#   The following calibrations are possible (case-insensitive): subframe, badpx_mask, bias, dark, flat, <background>, normalisation, combine_sum, localbackground'+os.linesep) 
        file.write('#   For dark and flat the paramerters can contain the the exposure time in float format (e.g. flat15.0_rawfiles, dark4.0_calibs_create).'+os.linesep)  
        file.write('#       If a fixed dark should be used, than the parameter needs to contain "dark" and additional text (e.g. "darkfixed", or if a different exposure time should be used "dark5.0")'+os.linesep) 
        file.write('#   For <background> the calibration needs to contain "background" but can contain more information. The key needs to be defined'+os.linesep) 
        file.write('#        (e.g. if the "background_image_filename" is used for calibration then the following entry is needed as well here or in the calibration file "background_image_filename = background.fits"'+os.linesep) 
        file.write('# '+os.linesep) 
        for paramtype in ['rawfiles', 'calibs_create', 'master']:       # sort: first rawfiles, master file names at end
            file.write(os.linesep)
            for entry in sorted(conf_data.keys()):
                if entry.find(paramtype) >= 0:
                    file.write('{0} = {1} {2}'.format(entry.ljust(24), conf_data[entry], os.linesep))
    
    logger('Info: The calibration file for handling the raw data has been created. Please check {0} before starting the data reduction (hiflex.py)'.format(params['configfile_fitsfiles']))
    
    

