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

if __name__ == "__main__":
    logger('Info: Preparing a list with the files')
    log_params(params)
    
    file_list = read_text_file(params['raw_data_file_list'], no_empty_lines=True)
    file_list = convert_readfile(file_list, [str, str, str, float, float, str], delimiter='\t', replaces=['\n',' '])
    
    for root, dirs, files in os.walk(params['raw_data_path']):
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
                        if filename.lower().find('flat') >= 0 and not filename.lower().find('sflat') >= 0:
                            fiber1, fiber2 = 'flat', 'flat'
                        if filename.lower().find('sflat') >= 0:
                            fiber1= 'sflat'
                        if params['raw_data_imtyp_keyword'] in im_head.keys():
                            if im_head[params['raw_data_imtyp_keyword']].replace(' ','') == params['raw_data_imtyp_flat'] and not filename.lower().find('sflat') >= 0:      # replace because when reading the conf file spaces are placed
                                fiber1, fiber2 = 'flat', 'flat'
                        if filename.lower().find('arc') >= 0:
                            fiber2 = 'wave'
                        if params['raw_data_imtyp_keyword'] in im_head.keys():
                            if im_head[params['raw_data_imtyp_keyword']].replace(' ','') == params['raw_data_imtyp_bias']:
                                fiber1, fiber2 = 'bias', 'bias'
                            if im_head[params['raw_data_imtyp_keyword']].replace(' ','') == params['raw_data_imtyp_dark']:
                                fiber1, fiber2 = 'dark', 'dark'
                            if im_head[params['raw_data_imtyp_keyword']].replace(' ','') == params['raw_data_imtyp_sflat']:
                                fiber1, fiber2 = 'sflat', 'dark'
                            if im_head[params['raw_data_imtyp_keyword']].replace(' ','') == params['raw_data_imtyp_arc']:
                                fiber1, fiber2 = 'wave', 'wave'
                        if (filename.lower().find('/arc') == -1) and not (filename.lower().find('arc') == 0) and fiber1 not in ['flat', 'sflat', 'dark', 'bias', 'wave']:
                            fiber1 = 'science'
                            if filename.lower().find('hars') <> -1:
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
                                print 'Problem with exptime',filename,temp
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
    file.write('### Each line contains the following information:\n')
    file.write('###   - Filename relativ to raw_data_path parameter\n')
    file.write('###   - Type of fiber 1 (science fiber)\n')
    file.write('###   - Type of fiber 2 (calibration fiber)\n')
    file.write('###   - Exposure time in seconds (from the header, if the information is not in the header, then from the filename)\n')
    file.write('###   - Observation time in Unix timestamp (from header)\n')
    file.write('###   - Flag "e", if the spectra of this file should be extraced. By standard only the science data is extracted.\n###         Combination of images before the extraction is possible, please refer to the manual for more information\n')
    file.write('### To exlude the use of some files please comment the line with a "#" or delete the line. \n\n')
    file.write('### -> When finished with the editing, save the file and close the editor \n\n')
    for entry in file_list:
        file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(entry[0].ljust(50), entry[1], entry[2], entry[3], entry[4], entry[5] ))
    file.close()
    rtn = os.system('gedit {0}'.format(params['raw_data_file_list']))
    if rtn <> 0:
        print('Please check that file {0} is correct.'.format(params['raw_data_file_list']))
        raw_input('To continue please press Enter\t\t')
    time.sleep(1)
    file_list = read_text_file(params['raw_data_file_list'], no_empty_lines=True)
    file_list = convert_readfile(file_list, [str, str, str, float, float], delimiter='\t', replaces=['\n',' '])
    
    file_list = sorted(file_list, key=operator.itemgetter(1,2,3,0))
    # Check what data is available
    arc_l_exp, arc_s_exp = 0, 1E10
    sflat_fib2 = 'wave'
    flatarc_fib2 = 'none'
    exist_bias, exist_flat, exp_darks, exist_flatarc = False, False, [], False
    for entry in file_list:
        if entry[1] == 'sflat' and entry[2] == 'none' and flatarc_fib2 <> 'dark':              # use entry[2] == 'wave' for sflat only if entry[2] == 'none' not available
            sflat_fib2 = 'none'
        if entry[1] == 'sflat' and entry[2] == 'dark':              # use entry[2] == 'dark' for sflat only if entry[2] == 'none' not available
            sflat_fib2 = 'dark'
        if entry[1] == 'sflat' and entry[2] == 'wave':               # use entry[2] == 'none' for arcflat only if entry[2] == 'wave' not available
            flatarc_fib2 = 'wave'
        if entry[1] == 'sflat' and entry[2] == 'dark' and flatarc_fib2 <> 'wave':   # use entry[2] == 'dark' for arcflat only if entry[2] == 'wave' not available
            flatarc_fib2 = 'dark'
        if entry[1] == 'none' and entry[2] == 'wave':
            arc_l_exp = max(entry[3], arc_l_exp)
            arc_s_exp = min(entry[3], arc_s_exp)
        if entry[1] == 'bias' and entry[2] == 'bias':
            exist_bias = True
        if entry[1] == 'flat' and entry[2] == 'flat':
            exist_flat = True
        if entry[1] == 'dark' and entry[2] == 'dark':
            if entry[3] not in exp_darks:
                exp_darks.append(entry[3])
    
    # Create the configuartion file
    conf_data = dict()
    warn = []
    pos_params = ['bias', 'flat' ]
    for entry in file_list:                                         # Extraction is set up below
        if entry[0].find('#') <> -1:
            continue
        param = ''
        for par in pos_params:
            if entry[1] == par and entry[2] == par:                 # Fiber1 and Fiber2
                param = par
                break
        if entry[1] == 'sflat' and entry[2] == flatarc_fib2:         # Fiber1 and Fiber2    , might be overwritten by sflat -> later copy sflat into flatarc
            param = 'flatarc'
        if entry[1] == 'sflat' and entry[2] == sflat_fib2:          # Fiber1 and Fiber2
            param = 'sflat'
        parcal = param
        
        if entry[1] == 'none' and entry[2] == 'wave':               # Fiber1 and Fiber2
            parcal = 'arc'
            if entry[3] == arc_l_exp:
                param = 'arc_l'
            elif entry[3] == arc_s_exp:
                param = 'arc_s'
        textparam = param
        if entry[1] == 'dark' and entry[2] == 'dark':               # Fiber1 and Fiber2
            param = 'dark{0}'.format(entry[3])                      # Exposure time
            textparam = param.replace('.','p')+'s'
            parcal = 'dark'
            
        if param <> '':
            if param+'_rawfiles' not in conf_data.keys():
                conf_data[param+'_rawfiles'] = entry[0]             # Filename
                conf_data['master_'+param+'_filename'] = 'master_'+textparam+'.fits'
                if parcal+'_calibs_create_g' in params.keys():
                    temp_param = params[parcal+'_calibs_create_g']
                elif 'standard_calibs_create' in params.keys():
                    temp_param = params['standard_calibs_create']
                else:
                    temp_param = []
                    logger('Warn: Missing entry in the configuration file. Neigther "{0}_calibs_create_g" nor "standard_calibs_create" is given. Please update the configuration file.'.format(parcal))
                text = ''
                for param_entry in temp_param:
                    if (param_entry == 'bias' and exist_bias == False) or (param_entry == 'flat' and exist_flat == False) or (param_entry == 'dark' and entry[3] not in exp_darks):          # If the calibration data is not available, then don't do the calibration
                        warn_text = 'Warn: The parameter {0}_calibs_create_g in the configuration file requires {1} correction, but this calibration data is not available.'.format(parcal,param_entry)
                        if param_entry == 'dark':
                            warn_text += ' Please assign dark of a different exposure time in {0}.'.format(params['configfile_fitsfiles'])
                        if warn_text not in warn:
                            warn.append(warn_text)
                        continue
                    text += ',' + param_entry
                if len(text) > 0:
                    text = text[1:]
                conf_data[param+'_calibs_create'] = text
            else:
                conf_data[param+'_rawfiles'] += ',' + entry[0]
    if 'sflat_rawfiles' in conf_data.keys() and not  'flatarc_rawfiles' in conf_data.keys():
        conf_data['flatarc_rawfiles'] = conf_data['sflat_rawfiles']
        conf_data['master_flatarc_filename'] = conf_data['master_sflat_filename'].replace('sflat','flatarc')
        conf_data['flatarc_calibs_create'] = conf_data['sflat_calibs_create']

    if 'arc_l_rawfiles' in conf_data.keys():
        conf_data['arc_rawfiles'] = conf_data['arc_l_rawfiles']
        conf_data['master_arc_filename'] = conf_data['master_arc_l_filename'].replace('arc_l','arc')
        conf_data['arc_calibs_create'] = conf_data['arc_l_calibs_create']
        if 'arc_s_rawfiles' not in conf_data.keys():
            conf_data['arc_s_rawfiles'] = conf_data['arc_l_rawfiles']
            conf_data['master_arc_s_filename'] = conf_data['master_arc_l_filename'].replace('arc_l','arc_s')
            conf_data['arc_s_calibs_create'] = conf_data['arc_l_calibs_create']
    # Define extraction
    for entry in file_list:
        if entry[0].find('#') <> -1:
            continue
        if entry[5].find('e') == -1:          # No extraction
            continue
        entry[5] = entry[5].replace(' ','').split(',')
        for extraction in entry[5]:
            if extraction.find('e') <> 0:        # use find as len(extraction) might be 0
                warn.append('Warn: I dont know what to do with the extraction parameter {0} (it doesnt begin with "e") for file {1}. This spectrum will therefore not be extracted. Please check {2} and run the script again, if you perform changes.'.format(extraction, entry[0], params['raw_data_file_list'] ))
                continue
            if extraction.find('c') == 1:            # combine data before extraction
                param = 'extract_combine'+extraction[2:]
            else:
                param = 'extract'+extraction[1:]
            if param+'_rawfiles' not in conf_data.keys():
                conf_data[param+'_rawfiles'] = entry[0]             # Filename
                if param.find('extract_combine') <> -1:                              # Combine only science data, which was requested for
                    conf_data['master_'+param+'_filename'] = 'master_'+param+'.fits'
                calibs = [param+'_calibs_create_g', param.replace('extract','science')+'_calibs_create_g', 'extract_calibs_create_g', 'science_calibs_create_g', 'standard_calibs_create', '']
                for calib in calibs:
                    if calib in params.keys():     # extract1_calibs_create_g
                        temp_param = params[calib]
                        break
                if calib == '':
                    temp_param = []
                    logger('Warn: Missing entry in the configuration file. Neigther "{0}" nor "{1} nor {2}" is given. Please update the configuration file.'.format(calibs[0], calibs[2], calibs[4]))
                text = ''
                for param_entry in temp_param:
                    if (param_entry == 'bias' and exist_bias == False) or (param_entry == 'flat' and exist_flat == False) or (param_entry == 'dark' and entry[3] not in exp_darks):          # If the calibration data is not available, then don't do the calibration
                        warn_text = 'Warn: The parameter {0} in the configuration file requires {1} correction, but this calibration data is not available'.format(calib,param_entry)
                        if param_entry == 'dark':
                            warn_text += ' Please assign dark of a different exposure time in {0}.'.format(params['configfile_fitsfiles'])
                        if warn_text not in warn:
                            warn.append(warn_text)
                        continue
                    text += ',' + param_entry       # Format list into text
                if len(text) > 0:
                    text = text[1:]
                conf_data[param+'_calibs_create'] = text
            else:
                conf_data[param+'_rawfiles'] += ',' + entry[0]

    for entry in warn:
        logger(entry)
    # Save the results in a conf_data.txt file
    #print json.dumps(conf_data,sort_keys = False, indent = 4)
    file = open(params['configfile_fitsfiles'],'w')
    file.write('# Description of this file:\n#--------------------------\n') 
    file.write('# This file was created by combining the information in the file {0} and the parameters given in the configuration file {1}\n'.format(params['raw_data_file_list'], CONFIGFILE))
    file.write('# Changes in here will be overwritten the next time prepare_file_list.py is run, therefore we suggest to make changes in {0} and afterwards run prepare_file_list.py again.\n\n'.format(params['raw_data_file_list']))
    file.write('# For each type of calibration filetype (e.g. bias, darks) a <filetype>_calibs_create list all corrections that will be applied to the individual files listed in <filetype>_rawfiles. These files will then combined and stored in master_<filetype>_filename (remove/comment the parameter in order to avoid saving.\n') 
    file.write('#   The following calibrations are possible (case-insensitive): subframe, badpx_mask, bias, dark, flat, <background>, normalisation\n') 
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
    
    logger('Info: The calibration file for handling the raw data has been created. Please check {0} before starting the calibration (reduction_day.py)'.format(params['configfile_fitsfiles']))
    
    

