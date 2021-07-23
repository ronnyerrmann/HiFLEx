import numpy as np
import os
from procedures import *

# =============================================================================
# Define variables
# =============================================================================
params = dict()     # default param dictionary
# location of config file
CONFIGFILE = 'conf.txt'

# Start of code
# deal with arguments from a text file
params = textfileargs(params, CONFIGFILE)
params['folder_original_traces'] = 'original_traces'

if __name__ == "__main__":
    logger('\nInfo: Starting modify traces')
    log_params(params)
    
    if os.path.exists(params['folder_original_traces']):
        os.system('mv {0} old_{0}'.format(params['folder_original_traces']))
    try:
            os.makedirs(params['folder_original_traces'])
    except:
            logger('Error: Cant create directory {0}'.format(params['folder_original_traces']))
            exit(1)
    
    im_flat, im_head = create_image_general(params, 'trace1')
    
    # Remove traces in the GUI
    
    polyfits_ori, xlows, xhighs, widths = read_fits_width(params['master_trace_sci_filename'])
    fmask, polyfits, widths, xlows, xhighs = remove_adjust_orders_UI( scale_image_plot(im_flat, 'log10'), polyfits_ori, xlows, xhighs, widths=widths, userinput=True, do_rm=True, do_add=True)
    if polyfits_ori.shape[0] < polyfits.shape[0]:
        logger('Info: Added {0} new orders'.format(polyfits.shape[0]-polyfits_ori.shape[0]))
    if polyfits_ori.shape[0] == polyfits[fmask].shape[0]:
        if np.all(polyfits_ori ==polyfits[fmask]):      # can't be tested if the test before fails
            logger('Info: No change made, therefore finished removing traces')
            exit(0)
    logger('Info: Removed the following orders: {0}'.format( np.where(fmask == False)[0] ))
    polyfits, xlows, xhighs, widths = polyfits[fmask], xlows[fmask], xhighs[fmask], widths[fmask]
    
    # Move original data into a sub-folder
    #os.system('mv {0} {1}/'.format(params['master_trace_sci_filename'], params['folder_original_traces']))
    params['master_wavelensolution_sci_filename'] = params['master_wavelensolution_filename'].replace('.fits','_sci.fits')
    params['master_wavelensolution_cal_filename'] = params['master_wavelensolution_filename'].replace('.fits','_cal.fits')
    for entry in ['master_trace_sci_filename', 'master_trace_cal_filename', 'master_wavelensolution_filename', 'master_wavelensolution_sci_filename', 
                  'master_wavelensolution_cal_filename', 'master_blaze_spec_norm_filename', 'master_wavelengths_shift_filename', 'logging_path', 'path_extraction']:
        if entry in params.keys(): 
            if os.path.exists(params[entry]):           # Covers folders, links, files
                os.system('mv {0} {1}/'.format(params[entry], params['folder_original_traces']))
                if entry.find('path_') != -1 or entry.find('_path') != -1:
                    if not os.path.exists(params[entry]) and entry not in ['raw_data_paths', 'path_ceres', 'terra_jar_file'] and params[entry].lower() not in ['na/', params['result_path']+'na/', params['result_path']+'/na/']:
                        try:                                                    # Create folders, if necessary
                            os.makedirs(params[entry])
                        except:
                            logger('Warn: Cannot create directory {0}'.format(params[entry]))
    
    # Save the new file and restart the reduction
    logger('Info: New file with the traces was created. All data depending on this will be moved to the folder {0}, the traces will be written.'.format(params['folder_original_traces']))
    save_fits_width(polyfits, xlows, xhighs, widths, params['master_trace_sci_filename'])
    plot_traces_over_image(scale_image_plot(im_flat, 'log10'), params['logging_traces_im'].replace('*', 'master_trace1'), polyfits, xlows, xhighs, widths)
    
    cmd = 'python {0}/hiflex.py {1}'.format(os.path.dirname(sys.argv[0]), ' '.join(sys.argv[1:]) )
    logger('Info: Finished modifying traces. Now, please run the pipeline again:\n{0}'.format(cmd))
    
    log_params(params)
    
    """
    logger('Warn: On some systems there might the process might hang when doing multicore extraction. If this is the case, kill the current process with CTRL+C and run hiflex manually with:\n{0}\n'.format(cmd))
    os.system(cmd)
     
    log_params(params)
    logger('Info: Finished removing traces and updating all the files for the day. \n')
    """
    
    
    
    
    
    
    
