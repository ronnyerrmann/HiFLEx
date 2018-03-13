import numpy as np
import os
from procedures import *

# =============================================================================
# Define variables
# =============================================================================
params = dict()     # default param dictionary
calimages = dict()  # dictionary for all calibration images
# location of config file
CONFIGFILE = 'conf.txt'

# Start of code
# deal with arguments from a text file
params = textfileargs(params, CONFIGFILE)
params['folder_original_traces'] = 'original_traces'

if __name__ == "__main__":
    logger('Info: Starting removing traces')
    log_params(params)
    
    if not os.path.exists(params['folder_original_traces']):
        try:
            os.makedirs(params['folder_original_traces'])
        except:
            logger('Error: Cant create directory {0}'.format(params['folder_original_traces']))
            exit(1)
    
    im_flat, im_head = create_image_general(params, 'sflat')
    
    # Remove traces in the GUI
    
    polyfits, xlows, xhighs, widths = read_fits_width(params['master_traces_filename'])
    fmask = run_remove_orders_UI(np.log10(im_flat), polyfits, xlows, xhighs, userinput=True)
    if len(polyfits) == len(polyfits[fmask]):
        logger('Info: No change made, therefore finished removing traces')
        exit(0)
    
    polyfits, xlows, xhighs, widths = polyfits[fmask], xlows[fmask], xhighs[fmask], widths[fmask]
    
    # Move original data into a sub-folder
    os.system('mv {0} {1}/'.format(params['master_traces_filename'], params['folder_original_traces']))
    
    for entry in ['background_filename', 'master_tracesarc_filename', 'master_wavelensolution_filename', 'master_flat_spec_norm_filename', 'logging_path', 'path_extraction']:
        if os.path.isfile(params[entry]) == True:
            os.system('mv {0} {1}/'.format(params[entry], params['folder_original_traces']))
    
    # Save the new file and restart the reduction
    logger('Info: New file with the traces was created. All data depending on this will be moved to the folder {0}, the traces will be written and reduction_day.py will be run'.format(params['folder_original_traces']))
    save_fits_width(polyfits, xlows, xhighs, widths, params['master_traces_filename'])
    plot_traces_over_image(im_flat, params['logging_traces_im'], polyfits, xlows, xhighs, widths)
    
    os.system('python {0}/reduction_day.py'.format(os.path.dirname(sys.argv[0])) )
     
    log_params(params)
    logger('Info: Finished removing traces and updating all the files for the day. \n')
    
    
    
    
    
    
    
    
