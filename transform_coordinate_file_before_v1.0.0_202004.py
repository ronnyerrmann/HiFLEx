#!/usr/bin/env python

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

object_files = [ params['object_file'] ]
object_file_no_path = params['object_file'].split('/')[-1]
for entry in [ params['result_path'] ] + params['raw_data_paths']:      # Check also the result and raw data paths for object names
        object_files.append( entry + object_file_no_path )
for filen in object_files:
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
        
    :params ignore_values: Just checks for the object name, but not any of the other fields, just creates a list of strings
    """
    RA0, DEC0, PMRA, PMDEC, epoch = 0., 0., 0., 0., 2000.
    
    if not os.path.isfile(filen):
        #logger('Warn: Reference coordinates file {0} does not exist.'.format(filen))
        continue
    # !!! Use read files with split already available
    lines_txt = read_text_file(filen, no_empty_lines=True, warn_missing_file=False)
    lines = convert_readfile(lines_txt, [str, str, str, str, str, str, str, str, str], delimiter=',', expand_input=True, replaces=[['\t',',']])    #
    for line in lines:
        lin5 = 0
        lin8 = 0
        try:
            lin5 = float(line[5])
        except:
            dummy = 0
        try:
            lin8 = float(line[8])
        except:
            dummy = 0
        convert = False
        if 1900 <= lin5 <= 2100 and not (1900 <= lin8 <= 2100):
            convert = True
        elif not (1900 <= lin5 <= 2100) and 1900 <= lin8 <= 2100:
            logger('It seems that file {0} was already converted.'.format(filen))
        elif 1900 <= lin5 <= 2100:
            logger('It seems that file {0} was needs conversion, but also column 10 (starting at 1) contains an epoch: {1}\n\t-> Exiting.'.format(filen,line[8]))
        else:
            logger('It seems that file {0} was already converted, as column 6 (starting at 1) does not contain a useful epoch: {1}\n\tBut also column 9 misses a useful epoch: {2}\n\t-> Exiting.'.format(filen,line[5],line[8]))
        if not convert:
            break
    if not convert:
        continue
        
    with open(filen,'w') as file:
        for entry in lines:
            file.write('{0},{1},{2},{3},{4},{6},{7},{8},{5}\n'.format( entry[0], entry[1], entry[2], entry[3], entry[4], entry[5], entry[6], entry[7], entry[8] ))  # Epoch (5) at end to be compatible with CERES
    logger('Converted file {0}'.format(filen))
            
