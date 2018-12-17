import numpy as np
import os
from procedures import *


# =============================================================================
# Define variables
# =============================================================================
#global calimages
params = dict()     # default param dictionary
#calimages = dict()  # dictionary for all calibration images
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
        if param.find('extract') == -1 or param.find('extract_combine') == 0 or param.find('wavelengthcal') == -1:      # wavecal and normal extraction files don't need a master, only the combined extraction needs one
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
        conf_data[param+'_rawfiles'] += ',' + entry[0]
    return conf_data, warn

if __name__ == "__main__":
    logger('Info: Preparing a list with the files')
    log_params(params)
    
    file_list = read_text_file(params['raw_data_file_list'], no_empty_lines=True)
    file_list = convert_readfile(file_list, [str, str, str, float, float, str], delimiter='\t', replaces=['\n',' '])
    
    if os.path.isdir(params['raw_data_path']) == False:
        if os.path.isdir(params['raw_data_path'].replace('ncook','ronny')) == True:
            params['raw_data_path'] = params['raw_data_path'].replace('ncook','ronny')

    for root, dirs, files in os.walk(params['raw_data_path'], followlinks=True):
        for file in files:
            for fileending in params['raw_data_file_endings']:
                if file.endswith(fileending):
                    filename = os.path.join(root, file).replace(params['raw_data_path'],'')     # Only relative folder and filename
                    if not os.path.exists(params['raw_data_path'] + filename):
                        continue
                    new = True
                    for entry in file_list:
                        if entry[0].find(filename) == 0 or entry[0].find(' '+filename) >= 0 or entry[0].find('#'+filename) >= 0:    # 3 entries as otherwise problem with SunArc-0... and Arc-0...
                            new = False
                            break
                    if new:
                        im_head = fits.getheader(params['raw_data_path'] + filename)
                        # Identify the image type and find out what the fibers show
                        fiber1, fiber2 = 'none', 'none'
                        if filename.lower().find('bias') >= 0:
                            fiber1, fiber2 = 'bias', 'bias'
                        if filename.lower().find('dark') >= 0:
                            fiber1, fiber2 = 'dark', 'dark'
                        if filename.lower().find('flat') >= 0 or filename.lower().find('whli') >= 0:
                            fiber1= 'cont'
                        if filename.lower().find('rflat') >= 0:
                            fiber1, fiber2 = 'rflat', 'rflat'
                        if params['raw_data_imtyp_keyword'] in im_head.keys():
                            if im_head[params['raw_data_imtyp_keyword']].replace(' ','') == params['raw_data_imtyp_flat']:      # replace because when reading the conf file spaces are placed
                                if filename.lower().find('rflat') >= 0:
                                    fiber1, fiber2 = 'rflat', 'rflat'
                                else:
                                    fiber1, fiber2 = 'cont', 'cont'
                        if filename.lower().find('arc') >= 0 or filename.lower().find('thar') >= 0 or filename.lower().find('une') >= 0:
                            fiber2 = 'wave'
                        if params['raw_data_imtyp_keyword'] in im_head.keys():
                            if im_head[params['raw_data_imtyp_keyword']].replace(' ','') == params['raw_data_imtyp_bias']:
                                fiber1, fiber2 = 'bias', 'bias'
                            if im_head[params['raw_data_imtyp_keyword']].replace(' ','') == params['raw_data_imtyp_dark']:
                                fiber1, fiber2 = 'dark', 'dark'
                            if im_head[params['raw_data_imtyp_keyword']].replace(' ','') == params['raw_data_imtyp_trace1']:
                                fiber1, fiber2 = 'cont', 'dark'
                            if im_head[params['raw_data_imtyp_keyword']].replace(' ','') == params['raw_data_imtyp_flatarc']:
                                fiber1, fiber2 = 'cont', 'wave'        # for HARPS it is cont, cont
                            if im_head[params['raw_data_imtyp_keyword']].replace(' ','') == params['raw_data_imtyp_trace2']:
                                fiber1, fiber2 = 'wave', 'wave'
                        if (filename.lower().find('/arc') == -1 and filename.lower().find('/thar') == -1 and filename.lower().find('/une') == -1) and not (filename.lower().find('arc') == 0 or filename.lower().find('thar') == 0 or filename.lower().find('une') == 0) and fiber1 not in ['rflat', 'cont', 'dark', 'bias', 'wave']:
                            fiber1 = 'science'
                            if filename.lower().find('harps') <> -1:
                                fiber2 = 'wave'
                        # Get the exposure time
                        exptime = -1
                        if params['raw_data_exptim_keyword'] in im_head.keys():
                            exptime = im_head[params['raw_data_exptim_keyword']]
                        else:
                            temp = filename.rsplit('_',1)
                            temp = temp[1].split('s')
                            try:
                                exptime = float(temp[0])
                            except:
                                logger('Warn: Problem with exposure time for file {1}. Header keyword {0} does not exist. Can not transform {2} into a number'.format(params['raw_data_exptim_keyword'], filename, temp[0]))
                        # Get the observation time
                        dateobs = 0
                        if params['raw_data_dateobs_keyword'] in im_head.keys():
                            dateobs = im_head[params['raw_data_dateobs_keyword']]
                            found_time = False
                            for timestring in ['%Y-%m-%dT%H:%M:%S.%f', '%Y-%m-%dT%H:%M:%S']:
                                try:
                                    dateobs = datetime.datetime.strptime(dateobs, timestring)
                                    found_time = True
                                except:
                                    found_time = False
                                if found_time:
                                    break
                            if found_time:        
                                dateobs = time.mktime(dateobs.timetuple())             #Time in seconds
                            else:
                                dateobs = 0
                        
                        extract = ''
                        if fiber1 == 'science':
                            extract = 'e'
                        file_list.append([filename, fiber1, fiber2, exptime, dateobs, extract])
                        #print file_list[-1]
    file_list = sorted(file_list, key=operator.itemgetter(1,2,3,0))
    # Save the list, show to user, so the user can disable files, read the list
    file = open(params['raw_data_file_list'],'w')
    file.write('### This file contains the information for all the fits files in the raw_data_path: {0} and it\'s subfolders.\n'.format(params['raw_data_path']))
    file.write('### Each line contains the following information, separated by one tab:\n')
    file.write('###   - Filename relativ to raw_data_path parameter\n')
    file.write('###   - Type of fiber 1 (science fiber)\n')
    file.write('###   - Type of fiber 2 (calibration fiber)\n')
    file.write('###   - Exposure time in seconds (from the header, if the information is not in the header, then from the filename)\n')
    file.write('###   - Observation time in Unix timestamp (from header)\n')
    file.write('###   - Flags:\n')
    file.write('###        "e", if the spectra of this file should be extraced. By standard only the science data is extracted.\n')
    file.write('###             Combination of images before the extraction is possible, please refer to the manual for more information\n')
    file.write('###        "w", if the spectrum contains wavelength information in case of a spectrograph with single fiber input.\n')
    file.write('###             Used to mark the files with calibration spectrum taken before or after the science observation.\n')
    file.write('### To exlude the use of some files please comment the line with a "#" or delete the line. \n\n')
    file.write('### -> When finished with the editing, save the file and close the editor \n\n')
    for entry in file_list:
        file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(entry[0].ljust(50), entry[1], entry[2], entry[3], entry[4], entry[5] ))
    file.close()
    if 'nocheck' not in sys.argv:
        start = time.time()
        rtn = os.system('{1} {0}'.format(params['raw_data_file_list'], params['editor'] ))
        if rtn <> 0 or time.time()-start < 10:
            print('Please check that file {0} is correct.'.format(params['raw_data_file_list']))
            raw_input('To continue please press Enter\t\t')
    time.sleep(1)
    file_list = read_text_file(params['raw_data_file_list'], no_empty_lines=True)
    file_list = convert_readfile(file_list, [str, str, str, float, float], delimiter='\t', replaces=['\n',' '])
    
    file_list = sorted(file_list, key=operator.itemgetter(1,2,3,0))
    # Check what data is available
    arc_l_exp, arc_s_exp = 0, 1E10
    arc_fib1 = ''
    sflat_fib2 = 'wave'
    flatarc_fib2 = 'none'
    exist_bias, exist_rflat, exp_darks = False, False, []
    for entry in file_list:
        if entry[1] == 'cont' and entry[2] == 'none' and flatarc_fib2 <> 'dark':              # use entry[2] == 'wave' for cont only if entry[2] == 'none' not available
            sflat_fib2 = 'none'
        if entry[1] == 'cont' and entry[2] == 'dark':              # use entry[2] == 'dark' for cont only if entry[2] == 'none' not available
            sflat_fib2 = 'dark'
        if entry[1] == 'cont' and entry[2] == 'wave':               # use entry[2] == 'none' for arcflat only if entry[2] == 'wave' not available
            flatarc_fib2 = 'wave'
        if entry[1] == 'cont' and entry[2] == 'dark' and flatarc_fib2 <> 'wave':   # use entry[2] == 'dark' for flatarc only if entry[2] == 'wave' not available
            flatarc_fib2 = 'dark'
        if (entry[1] == 'none' or entry[1] == 'dark' or entry[1] == 'wave') and entry[2] == 'wave':
            arc_l_exp = max(entry[3], arc_l_exp)
            arc_s_exp = min(entry[3], arc_s_exp)
        if entry[1] == 'bias' and entry[2] == 'bias':
            exist_bias = True
        if entry[1] == 'rflat' and entry[2] == 'rflat':
            exist_rflat = True
        if entry[1] == 'dark' and entry[2] == 'dark':
            if entry[3] not in exp_darks:
                exp_darks.append(entry[3])
    if arc_l_exp == 0 and arc_s_exp == 1E10:         # no single arcs taken, use cont and arc
        arc_fib1 = 'cont'
        for entry in file_list:                     # re-run to find exposure times
            if (entry[1] == arc_fib1) and entry[2] == 'wave':
                arc_l_exp = max(entry[3], arc_l_exp)
                arc_s_exp = min(entry[3], arc_s_exp)
    
    conf_data = dict()
    # Check if master files exist:
    if not exist_rflat:
        exist_rflat = os.path.isfile('master_rflat.fits')
        conf_data['master_rflat_filename'] = 'master_rflat.fits'
    if not exist_bias:
        exist_bias = os.path.isfile('master_bias.fits')
        conf_data['master_bias_filename'] = 'master_bias.fits'
    for file in os.listdir("."):
        if file.endswith(".fits"):
            filename =os.path.join("", file)
            if filename.find('master_dark') == 0:
                exposure = filename.replace('master_dark','').replace('s.fits','').replace('p','.')
                try:
                    exposure = float( exposure )
                except:
                    continue
                if exposure not in exp_darks:
                    exp_darks.append(exposure)
                    conf_data['master_dark{}_filename'.format(exposure)] = filename
    
    # Create the configuartion file
    warn = []
    pos_params = ['bias', 'rflat' ]
    for entry in file_list:                                         # Extraction is set up below
        if entry[0].find('#') <> -1:        # comment # found
            continue
        param = ''
        for par in pos_params:
            if entry[1] == par and entry[2] == par:                 # Fiber1 and Fiber2
                conf_data, warn = create_parameters(conf_data, warn, par, par, [par], exist_bias, exist_rflat, exp_darks, entry)
        if entry[1] == 'cont' and entry[2] == flatarc_fib2:         # Fiber1 and Fiber2    , might be overwritten by cont -> later copy cont into flatarc
            conf_data, warn = create_parameters(conf_data, warn, 'flatarc', 'flatarc', ['flatarc'], exist_bias, exist_rflat, exp_darks, entry)
        if entry[1] == 'cont' and entry[2] == sflat_fib2:          # Fiber1 and Fiber2
            conf_data, warn = create_parameters(conf_data, warn, 'trace1', 'trace1', ['trace1'], exist_bias, exist_rflat, exp_darks, entry)
        if (entry[1] == 'none' or entry[1] == 'dark' or entry[1] == 'wave' or entry[1] == arc_fib1) and entry[2] == 'wave':               # Fiber1 and Fiber2
            parcal = 'arc'
            if entry[3] == arc_l_exp:
                conf_data, warn = create_parameters(conf_data, warn, 'arc_l', 'arc_l', [parcal], exist_bias, exist_rflat, exp_darks, entry)
                conf_data, warn = create_parameters(conf_data, warn, 'trace2', 'trace2', ['trace2'], exist_bias, exist_rflat, exp_darks, entry)
            elif entry[3] == arc_s_exp:                             # In case only one exposure time -> arc_l is copied in arc_s
                conf_data, warn = create_parameters(conf_data, warn, 'arc_s', 'arc_s', [parcal], exist_bias, exist_rflat, exp_darks, entry)
        if entry[1] == 'dark' and entry[2] == 'dark':               # Fiber1 and Fiber2
            param = 'dark{0}'.format(entry[3])                      # Exposure time
            textparam = param.replace('.','p')+'s'
            parcal = 'dark'
            conf_data, warn = create_parameters(conf_data, warn, param, textparam, [parcal], exist_bias, exist_rflat, exp_darks, entry)
        if entry[5].lower().find('e') == -1 and entry[5].lower().find('w') == -1:          # No extraction and no wavelength calibration
            continue
        entry[5] = entry[5].replace(' ','').split(',')
        for extraction in entry[5]:
            if extraction.lower().find('w') == 0:       # use this file for wavelength calibration between spectra
                param = 'wavelengthcal'
                conf_data, warn = create_parameters(conf_data, warn, param, param, [param], exist_bias, exist_rflat, exp_darks, entry)
                continue                        # otherwise Warn from below will be triggered
            if extraction.lower().find('e') <> 0:        # use find as len(extraction) might be 0
                warn.append('Warn: I dont know what to do with the extraction parameter {0} (it doesnt begin with "e") for file {1}. This spectrum will therefore not be extracted. Please check {2} and run the script again, if you perform changes.'.format(extraction, entry[0], params['raw_data_file_list'] ))
                continue
            if extraction.lower().find('c') == 1:            # combine data before extraction
                param = 'extract_combine'+extraction[2:]
            else:
                param = 'extract'+extraction[1:]
            conf_data, warn = create_parameters(conf_data, warn, param, param, [param, 'extract'], exist_bias, exist_rflat, exp_darks, entry)
    if 'arc_l_rawfiles' in conf_data.keys() and 'arc_s_rawfiles' not in conf_data.keys():
        conf_data['arc_s_rawfiles'] = conf_data['arc_l_rawfiles']
        conf_data['master_arc_s_filename'] = conf_data['master_arc_l_filename'].replace('arc_l','arc_s')
        conf_data['arc_s_calibs_create'] = conf_data['arc_l_calibs_create']
    for entry in warn:
        logger(entry)
    # Save the results in a conf_data.txt file
    #print json.dumps(conf_data,sort_keys = False, indent = 4)
    file = open(params['configfile_fitsfiles'],'w')
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
        for entry in sorted(conf_data):
            if entry.find(paramtype) >= 0:
                file.write('{0} = {1} \n'.format(entry.ljust(24), conf_data[entry]))
    file.close()
    
    logger('Info: The calibration file for handling the raw data has been created. Please check {0} before starting the data reduction (reduction_day.py)'.format(params['configfile_fitsfiles']))
    
    

