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
    logger('\nInfo: Starting removing traces')
    log_params(params)
    
    if not os.path.exists(params['folder_original_traces']):
        try:
            os.makedirs(params['folder_original_traces'])
        except:
            logger('Error: Cant create directory {0}'.format(params['folder_original_traces']))
            exit(1)
    
    im_flat, im_head = create_image_general(params, 'trace1')
    
    # Remove traces in the GUI
    
    polyfits, xlows, xhighs, widths = read_fits_width(params['master_trace_sci_filename'])
    fmask, dummy, dummy = remove_adjust_orders_UI( scale_image_plot(im_flat, 'log10'), polyfits, xlows, xhighs, userinput=True, do_rm=True)
    if len(polyfits) == len(polyfits[fmask]):
        logger('Info: No change made, therefore finished removing traces')
        exit(0)
    logger('Info: Removed the following orders: {0}'.format( np.where(fmask == False)[0] ))
    polyfits, xlows, xhighs, widths = polyfits[fmask], xlows[fmask], xhighs[fmask], widths[fmask]
    
    # Move original data into a sub-folder
    #os.system('mv {0} {1}/'.format(params['master_trace_sci_filename'], params['folder_original_traces']))
    params['master_wavelensolution_sci_filename'] = params['master_wavelensolution_filename'].replace('.fits','_sci.fits')
    for entry in ['background_filename', 'master_trace_sci_filename', 'master_trace_cal_filename', 'master_wavelensolution_filename', 'master_blaze_spec_norm_filename', 'logging_path', 'path_extraction']:
        if os.path.exists(params[entry]):           # Covers folders, links, files
            os.system('mv {0} {1}/'.format(params[entry], params['folder_original_traces']))
            if entry.find('path_') <> -1 or entry.find('_path') <> -1:
                if not os.path.exists(params[entry]) and entry not in ['raw_data_paths', 'path_ceres', 'terra_jar_file'] and params[entry].lower() not in ['na/', params['result_path']+'na/', params['result_path']+'/na/']:
                    try:                                                    # Create folders, if necessary
                        os.makedirs(params[entry])
                    except:
                        logger('Warn: Cannot create directory {0}'.format(params[entry]))
    
    # Save the new file and restart the reduction
    logger('Info: New file with the traces was created. All data depending on this will be moved to the folder {0}, the traces will be written and reduction_day.py will be run'.format(params['folder_original_traces']))
    save_fits_width(polyfits, xlows, xhighs, widths, params['master_trace_sci_filename'])
    plot_traces_over_image(im_flat, params['logging_traces_im'], polyfits, xlows, xhighs, widths)
    
    os.system('python {0}/reduction_day.py'.format(os.path.dirname(sys.argv[0])) )
     
    log_params(params)
    logger('Info: Finished removing traces and updating all the files for the day. \n')
    
    
    
    
    
    
    
    
