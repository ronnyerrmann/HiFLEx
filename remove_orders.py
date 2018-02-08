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
params['folder_original_orders'] = 'original_orders'

if __name__ == "__main__":
    logger('Info: Starting removing orders')
    log_params(params)
    
    if not os.path.exists(params['folder_original_orders']):
        try:
            os.makedirs(params['folder_original_orders'])
        except:
            logger('Error: Cant create directory {0}'.format(params['folder_original_orders']))
            exit(1)
    
    im_flat, im_head = create_image_general(params, 'sflat')
    
    # Remove orders in the GUI
    
    polyfits, xlows, xhighs, widths = read_fits_width(params['master_order_filename'])
    fmask = run_remove_orders(np.log10(im_flat), polyfits, xlows, xhighs, userinput=True)
    if len(polyfits) == len(polyfits[fmask]):
        logger('Info: No change made, therefore finished removing orders')
        exit(0)
    
    polyfits, xlows, xhighs, widths = polyfits[fmask], xlows[fmask], xhighs[fmask], widths[fmask]
    # Move original data into a sub-folder and save the data
    os.system('mv {0} {1}/'.format(params['master_order_filename'], params['folder_original_orders']))
    os.system('mv {0} {1}/'.format(params['logging_orders'], params['folder_original_orders']))
    save_fits_width(polyfits, xlows, xhighs, widths, params['master_order_filename'])
    plot_traces_over_image(im_flat, params['logging_orders'], polyfits, xlows, xhighs, widths)
    
    logger('Info: New file with the orders has been written. Now all data depending on this will be moved to the folder {0} and reduction_day.py will be run'.format(params['folder_original_orders']))
    for entry in ['background_image_filename', 'master_orderarc_filename', 'master_arc_solution_filename', 'master_flat_spec_norm_filename', 'logging_path', 'path_extraction']:
        #'logging_background', 'logging_arcorders', 'logging_arc_line_identification_positions', 'logging_arc_line_identification_residuals']:
        if os.path.isfile(params[entry]) == True:
            os.system('mv {0} {1}/'.format(params[entry], params['folder_original_orders']))
    
    os.system('python {0}/reduction_day.py'.format(os.path.dirname(sys.argv[0])) )
     
    log_params(params)
    logger('Info: Finished removing orders and updating all the files for the day. \n')
    
    
    
    
    
    
    
    
