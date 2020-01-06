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
        conf_data[param+'_rawfiles'] += ', ' + entry[0]

    return conf_data, warn

def add_new_rawfiles_file_list(params, file_list=[]):
    for raw_data_path in params['raw_data_paths']:
     for root, dirs, files in os.walk(raw_data_path, followlinks=True):
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
            fnlow = filename.lower()
            if fnlow.find('bias') >= 0:                     # hardcoded: Bias
                fiber1, fiber2 = 'bias', 'bias'
            if fnlow.find('dark') >= 0:                     # hardcoded: Dark
                fiber1, fiber2 = 'dark', 'dark'
            posi  = [fnlow.find('flat') , fnlow.find('whli') , fnlow.find('white') , fnlow.find('tung') ]      # hardcoded: White light (Tungston) spectrum in science fiber
            posi2 = [fnlow.find('flat2'), fnlow.find('whli2'), fnlow.find('white2'), fnlow.find('tung2')]      # hardcoded: White light (Tungston) spectrum in calibration fiber
            for i in range(len(posi)):
                if posi2[i] >= 0:
                    fiber2 = 'cont'
                    if posi[i] <> posi2[i]:                 # both fibers
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
            posi  = [fnlow.find('arc') , fnlow.find('thar') , fnlow.find('th_ar') , fnlow.find('thorium') , fnlow.find('une') ]        # hardcoded: emission line spectrum in calibration fiber
            posi2 = [fnlow.find('arc2'), fnlow.find('thar2'), fnlow.find('th_ar2'), fnlow.find('thorium2'), fnlow.find('une2')]        # hardcoded: emission line spectrum in science fiber
            for i in range(len(posi)):
                if posi2[i] >= 0:
                    fiber1 = 'wave'
                    if posi[i] <> posi2[i]:                 # both fibers
                        fiber2 = 'wave'
                elif posi[i] >= 0:
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
                if fnlow.find('harps') <> -1:
                    fiber2 = 'wave'
            im_head, obsdate_midexp, obsdate_mid_float, jd_midexp = get_obsdate(params, im_head)    # dateobs: unix_timestamp of mid exposure time
            extract = ''
            file_list.append([filename, fiber1, fiber2, im_head['HIERARCH HiFLEx EXPOSURE'], obsdate_mid_float, extract])
    
    return file_list

def check_raw_files_infos(params, file_list):
    # Check what data is available
    cal2_l_exp, cal2_s_exp, cal1_l_exp, cal1_s_exp = 0, 1E10, 0, 1E10
    arc_fib1 = ''
    sflat_fib2 = 'wave'
    blazecor_fib2 = 'none'
    exist_bias, exist_rflat, exp_darks = False, False, []
    for entry in file_list:
        if entry[1] == 'cont' and entry[2] == 'none' and blazecor_fib2 <> 'dark':              # use entry[2] == 'wave' for cont only if entry[2] == 'none' not available
            sflat_fib2 = 'none'
        if entry[1] == 'cont' and entry[2] == 'dark':              # use entry[2] == 'dark' for cont only if entry[2] == 'none' not available
            sflat_fib2 = 'dark'
        if entry[1] == 'cont' and entry[2] == 'wave':               # use entry[2] == 'none' for arcflat only if entry[2] == 'wave' not available
            blazecor_fib2 = 'wave'
        if entry[1] == 'cont' and entry[2] == 'dark' and blazecor_fib2 <> 'wave':   # use entry[2] == 'dark' for blazecor only if entry[2] == 'wave' not available
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
            if params['cal2_l_exp']*0.9 <= entry[3] <= params['cal2_l_exp']*1.1:
                extract += 'w2l, '
                extract += 't2, '
            elif params['cal2_s_exp']*0.9 <= entry[3] <= params['cal2_s_exp']*1.1:                             # In case only one exposure time -> cal2_l is copied in cal2_s
                extract += 'w2s, '
        if (entry[2] == 'none' or entry[2] == 'dark' or entry[2] == 'wave') and entry[1] == 'wave' and params['cal1_l_exp'] >= params['cal1_s_exp']:               # wave in science fiber
            parcal = 'arc'
            if params['cal1_l_exp']*0.9 <= entry[3] <= params['cal1_l_exp']*1.1:
                extract += 'w1l, '
            elif params['cal1_s_exp']*0.9 <= entry[3] <= params['cal1_s_exp']*1.1:                             # In case only one exposure time -> cal1_l is copied in cal1_s
                extract += 'w1s, '
        if entry[1] == 'science':
                extract += 'e, '
        if len(extract) >=2:
            extract = extract[:-2]
        file_list[index][5] = extract
    
    return file_list

def create_configuration_file(params, file_list):
    #cal2_l_exp   = params['cal2_l_exp']
    #cal2_s_exp   = params['cal2_s_exp']
    #cal1_l_exp   = params['cal1_l_exp']
    #cal1_s_exp   = params['cal1_s_exp']
    #arc_fib1     = params['arc_fib1']
    #sflat_fib2   = params['sflat_fib2']
    #blazecor_fib2 = params['blazecor_fib2']
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
    for file in os.listdir("."):
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
    easy_assignments.append(['w2s' , 'cal2_s'  , 'arc'    ])       # Wavelength calibration (standard)
    easy_assignments.append(['w2l' , 'cal2_l'  , 'arc'    ])       # Wavelength calibration (standard)
    for entry in file_list:                                         # Extraction is set up below
        #if entry[0].find('#') <> -1:        # comment # found
        #    continue
        #print entry, conf_data
        param = ''
        extract = entry[5].lower().replace(' ','').split(',')
        for easy_assignment in easy_assignments:
            if easy_assignment[0] in extract:
                conf_data, warn = create_parameters(conf_data, warn, easy_assignment[1], easy_assignment[1], [easy_assignment[2]], exist_bias, exist_rflat, exp_darks, entry)
        if 'd' in extract:               # Fiber1 and Fiber2
            param = 'dark{0}'.format(entry[3])                      # Exposure time
            textparam = param.replace('.','p')+'s'
            parcal = 'dark'
            conf_data, warn = create_parameters(conf_data, warn, param, textparam, [parcal], exist_bias, exist_rflat, exp_darks, entry)
        if entry[5].lower().find('e') == -1 and entry[5].lower().find('w') == -1:          # No extraction and no wavelength calibration
            continue
        for extraction in entry[5].replace(' ','').split(','):
            if extraction.lower() == 'w2':       # use this file for wavelength calibration between spectra
                param = 'wavelengthcal2'
                parcal = 'wavelengthcal'
                conf_data, warn = create_parameters(conf_data, warn, param, param, [parcal], exist_bias, exist_rflat, exp_darks, entry)
                continue  
            elif extraction.lower() == 'w':       # use this file for wavelength calibration between spectra
                param = 'wavelengthcal'
                conf_data, warn = create_parameters(conf_data, warn, param, param, [param], exist_bias, exist_rflat, exp_darks, entry)
                continue                        # otherwise Warn from below will be triggered
            if extraction.lower().find('e') <> 0:        # use find as len(extraction) might be 0
                # warn.append('Warn: I dont know what to do with the extraction parameter {0} (it doesnt begin with "e") for file {1}. This spectrum will therefore not be extracted. Please check {2} and run the script again, if you perform changes.'.format(extraction, entry[0], params['raw_data_file_list'] ))
                continue                                # The lines remaining in the for loop are only to be executed when it's about extraction, otherwise something might readded another time
            if extraction.lower().find('ec') == 0:            # combine data before extraction
                param = 'extract_combine'+extraction[2:]
            elif extraction.lower().find('e') == 0:           # just extraction:
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
    common_path = ''
    if file_list[0][0][2:].find('/') != -1:         # it contains subfolders
        fnames = []
        for entry in file_list:
            fnames.append( entry[0].replace('#','').replace(' ','') )
        common_path = fnames[0].rsplit('/',1)[0]    # just remove the filename
        for dummy in range(len(fnames[0].split('/'))-2):
            found_common_path = True
            for entry in fnames:
                if entry.find(common_path) == -1:
                    found_common_path = False
                    break
            if found_common_path:
                common_path += '/'
                break
            else:
                common_path = common_path.rsplit('/',1)[0]      # remove the next subfolder
    # Split the common_path to make the GUI shorter (in x), if necessary
    common_path_text = common_path
    if len(common_path) > 40:
        common_path_temp = common_path[:-1].split('/')          # without the last '/'
        common_path_text = [common_path_temp[0]+'/']            # start with the first entry
        ii = 0
        for entry in common_path_temp[1:]:
            entry += '/'
            if len(common_path_text[ii]+entry) > 40 and len(common_path_text[ii]) > 20:     # If too long to add
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
    widgets['w2l']     = dict(label='Wave\nlong', kind='Label', row=0, column=12)
    widgets['w2s']     = dict(label='Wave\nshort', kind='Label', row=0, column=13)
    widgets['w1l']     = dict(label='Wave\nSci.\nlong', kind='Label', row=0, column=14)
    widgets['w1s']     = dict(label='Wave\nSci.\nshort', kind='Label', row=0, column=15)
    widgets['w']       = dict(label='Wave\nshft\ncal', kind='Label', row=0, column=16)
    widgets['w2']      = dict(label='Wave\nSci\nshft', kind='Label', row=0, column=17)
    widgets['e']       = dict(label='Ex-\ntract', kind='Label', row=0, column=18)
    widgets['extra']   = dict(label='Further usage\nof files\n(comma separated)', kind='Label', row=0, column=19)
    for ii, entry in enumerate(file_list):
        pkwargs['comment_{0}'.format(ii)] = ( 0 <= entry[0].find('#') < 20 )            # Commented out
        widgets['comment_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['comment_{0}'.format(ii)], row=ii+1, column=0)
        text = '{0}{1}   {2}'.format( datetime.datetime.utcfromtimestamp(entry[4]).strftime('%Y-%m-%d %H:%M:%S'), 
                                     '%9.2f'%entry[3], entry[0].replace('#','').replace(' ','').replace(common_path,'') )
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
        pkwargs['w2l_{0}'.format(ii)] = ( 'w2l' in flags )
        widgets['w2l_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['w2l_{0}'.format(ii)], row=ii+1, column=12)
        pkwargs['w2s_{0}'.format(ii)] = ( 'w2s' in flags )
        widgets['w2s_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['w2s_{0}'.format(ii)], row=ii+1, column=13)
        pkwargs['w1l_{0}'.format(ii)] = ( 'w1l' in flags )
        widgets['w1l_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['w1l_{0}'.format(ii)], row=ii+1, column=14)
        pkwargs['w1s_{0}'.format(ii)] = ( 'w1s' in flags )
        widgets['w1s_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['w1s_{0}'.format(ii)], row=ii+1, column=15)
        pkwargs['w_{0}'.format(ii)] = ( 'w' in flags )
        widgets['w_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['w_{0}'.format(ii)], row=ii+1, column=16)
        pkwargs['w2_{0}'.format(ii)] = ( 'w2' in flags )
        widgets['w2_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['w2_{0}'.format(ii)], row=ii+1, column=17)
        pkwargs['e_{0}'.format(ii)] = ( 'e' in flags )
        widgets['e_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['e_{0}'.format(ii)], row=ii+1, column=18)
        extra = ''
        for flag in flags:          # Add extra flags, if necessary
            if flag not in ['b', 'd', 'a', 't1', 't2', 'z', 'w2l', 'w2s', 'w1l', 'w1s', 'w', 'w2', 'e']:
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
              '-- Wavelength solition (long and short expsoure time)\n'+\
              '-- Wavelength solution for science fiber (*)\n   (long or short exposure time)\n'+\
              '-- Wavelength shift calibration (**) to correct for\n   wavelength shift\n'+\
              '-- Wavelength shift between the Science fiber and\n   Calibration fiber (*)\n'+\
              '-- Extract: Extract these files on an individual basis\n'+\
              '-- Further settings (manual): e.g. to combine files\n   before extraction\n'+\
              '(*) not for single fiber spectrographs\n'+\
              '(**) important for unstabilised single fiber\n     spectrographs'
              #'- Type of Science and Calibration\n  fibers are derived from header or\n  filename and can be changed here\n  (optional)\n'+\     # Allows modification of the fiber content
    for ii, commentii in enumerate(explain.split('\n')):
        widgets['explain_{0}'.format(ii)] = dict(label=commentii, kind='Label', row=ii, column=20, rowspan=1, orientation=Tk.W )#, wraplength=100 )      
    widgets['accept'] = dict(label='Accept', kind='ExitButton', row=ii+1, column=20, rowspan=2)
                              
    wprops = dict(fullscreen=False )
    #wprops['width_data'] = 800   # not neccssary, as uses the automatic width
    
    if len(file_list) > 100:
        logger('Info: A GUI with {0} elements will be created, on some machines that can take up to a few minutes.'.format(len(widgets)))

    gui3 = tkc.TkCanvasGrid(title='HiFlEx: Asigning Files to extraction steps (What file contains what data)', 
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
                 {True:'w2l',False:''}[gui3.data['w2l_{0}'.format(ii)]] +','+
                 {True:'w2s',False:''}[gui3.data['w2s_{0}'.format(ii)]] +','+
                 {True:'w1l',False:''}[gui3.data['w1l_{0}'.format(ii)]] +','+
                 {True:'w1s',False:''}[gui3.data['w1s_{0}'.format(ii)]] +','+
                 {True:'w',False:''}[gui3.data['w_{0}'.format(ii)]] +','+
                 {True:'w2',False:''}[gui3.data['w2_{0}'.format(ii)]] +','+
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

if __name__ == "__main__":
    logger('Info: Preparing a list with the files')
    log_params(params)
    
    # get the available list
    file_list = read_text_file(params['raw_data_file_list'], no_empty_lines=True)
    try:
        file_list = convert_readfile(file_list, [str, str, str, float, float, str], delimiter='\t', replaces=['\n',' ']) # old way of reading the data, To stay backwards compatible, can be removed in a few versions after v0.4.1
    except:
        #file_list = convert_readfile(file_list, [str, str, str, float, ['%Y-%m-%dT%H:%M:%S', float], str], delimiter='\t', replaces=['\n',' '], ignorelines=[['#',20]])     #new way of storing the data
        file_list = convert_readfile(file_list, [str, str, str, float, ['%Y-%m-%dT%H:%M:%S', float], str], delimiter='\t', replaces=['\n',' '])     #new way of storing the data
    number_old_entries = len(file_list)
    
    # get the new files
    file_list = add_new_rawfiles_file_list(params, file_list)           # Check for new entries for file_list
    
    # analyse the global information
    params = check_raw_files_infos(params, file_list)
    
    # add the configuration to extraction parameters
    file_list = add_extraction_parameters_file_list(params, file_list, number_old_entries)

    file_list = sorted(file_list, key=operator.itemgetter(4,0))           #itemgetter(1,2,3,4,0)
    
    # Show the results in a GUI
    file_list, file_list_commented = file_list_UI(file_list)
    
    # Save the list, show to user, so the user can disable files, read the list
    file = open(params['raw_data_file_list'],'w')
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
    file.write('###        "w2l, w2s", Spectruum to find the wavelength solution (long and short exposure time) [of the calibration fiber].\n')
    file.write('###        "w1l, w1s", Spectruum to find the wavelength solution of the science fiber.\n')
    file.write('###        "e", if the spectra of this file should be extraced. By standard only the science data is extracted.\n')
    file.write('###             Combination of images before the extraction is possible, please refer to the manual for more information\n')
    file.write('###        "w" or "w2", if the spectrum contains wavelength information (in case of a spectrograph with single fiber input or\n')
    file.write('###                     to find the time dependent offset between both fibers of a bifurcated fiber).\n')
    file.write('###             Used to mark the files with calibration spectrum taken before or after the science observation.\n')
    file.write('### To exlude the use of some files please comment the line with a "#" or delete the line. \n\n')
    file.write('### -> When finished with the editing, save the file and close the editor \n\n')
    for entry in file_list_commented:
        file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(entry[0].ljust(50), entry[1], entry[2], entry[3], datetime.datetime.utcfromtimestamp(entry[4]).strftime('%Y-%m-%dT%H:%M:%S'), entry[5] ))
    file.close()
    
    """
    if not ('nocheck' in sys.argv or '-nocheck' in sys.argv or '--nocheck' in sys.argv):
        start = time.time()
        rtn = os.system('{1} {0}'.format(params['raw_data_file_list'], params['editor'] ))
        if rtn <> 0 or time.time()-start < 10:
            print('Please check that file {0} is correct.'.format(params['raw_data_file_list']))
            raw_input('To continue please press Enter\t\t')
     
    time.sleep(0.3)
    file_list = read_text_file(params['raw_data_file_list'], no_empty_lines=True)
    file_list = convert_readfile(file_list, [str, str, str, float, ['%Y-%m-%dT%H:%M:%S', float], str], delimiter='\t', replaces=['\n',' '], ignorelines=[['#',20]])
    
    file_list = sorted(file_list, key=operator.itemgetter(4,0))       # itemgetter(1,2,3,0)
    """
    
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
    params = check_raw_files_infos(params, file_list)
    
    conf_data = create_configuration_file(params, file_list)
    
    # Create a GUI that let's the user modify the configuration data
    
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
    
    

