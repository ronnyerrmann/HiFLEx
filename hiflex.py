#!/usr/bin/env python

""" This was thought to be a way to easily improve calculation speed. However, only 2s in 5min are saved. It seems that using numpy for all big calculations is enough
try:
    #import cython
    import pyximport
    pyximport.install()
    from cy_procedures import *
except:
    #print('cython not available')
    from procedures import *
"""
from procedures import *

# =============================================================================
# Define variables
# =============================================================================
params = dict()     # default param dictionary
global calimages    # dictionary for all calibration images
# location of config file
CONFIGFILE = 'conf.txt'

# Start of code
# deal with arguments from a text file
params = textfileargs(params, CONFIGFILE)
if __name__ == "__main__":
    logger('Info: Starting routines for a new night of data, including: finding or shifting orders, find calibration orders, create wavelength solution, and create normalised blaze function', params=params)
    params['path_run'] = os.getcwd()
    params['extract_wavecal'] = False
    params['no_RV_names'] = ['flat', 'tung', 'whili', 'thar', 'th_ar', 'th-ar', 'laser']
    if platform.system() == 'Windows':
        params['use_cores'] = 1
        logger('Info: Running on Windows, disabled multicore.')
    log_params(params, start=True)
    
    # create the median combined files
    im_trace1, im_trace1_head = create_image_general(params, 'trace1')
    if params['arcshift_side'] != 0 or 'trace2_rawfiles' in params.keys():      # calibration spectrum at the same time
        im_trace2, im_trace2_head = create_image_general(params, 'trace2')
    else:                                                                       # no calibration spectrum at the same time
        im_trace2, im_trace2_head = im_trace1, im_trace1_head
    #blazecor, cal2_l, cal2_s: at a later step to know the orders for localbackground -> cont1cal2
    
    # Load the reference catalogue. Note that the Ar I lines are rescaled!
    reference_lines_dict = read_reference_catalog(params)
    
    # Create or read the file with the orders for this night
    if os.path.isfile(params['master_trace_sci_filename']) :
        logger('Info: Using exiting trace solution: {0}'.format(params['master_trace_sci_filename']), params=params)
        sci_tr_poly, xlows, xhighs, widths = read_fits_width(params['master_trace_sci_filename'])
    else:
        printresults = False
        # load the original solution
        if os.path.isfile(params['original_master_traces_filename']) :
            sci_tr_poly, xlows, xhighs, widths = read_fits_width(params['original_master_traces_filename'])
            # find the shift between the original solution and the current flat
            xhighs[ xhighs > im_trace1.shape[0] ] = im_trace1.shape[0]      # otherwise problems if the image is cut a bit
            shift, widths_new, shift_map, shift_error = shift_orders(im_trace1, params, sci_tr_poly, xlows, xhighs, widths, params['in_shift'])
            # save the map of the shifts
            save_im_fits(params, shift_map, im_trace1_head, params['logging_map_shift_orders'])
        else:
            shift_error = -1
        if shift_error > 1 or shift_error == -1 or abs(shift) > params['maxshift_orders']:
            logger('Warn: The deviation of the shift of the orders seems too big or no previous solution was available, therefore searching for the position of the orders from scratch:', params=params)
            sci_tr_poly, xlows, xhighs, widths = find_adjust_trace_orders(params, im_trace1, im_trace1_head)
            printresults = True
        else:
            sci_tr_poly[:,:,-1] += shift                # update the sci_tr_poly parameters
            if params['update_width_orders'] :
                sci_tr_poly, widths = update_tr_poly_width_multiplicate(sci_tr_poly, widths, [widths_new[:,0]/widths[:,0], widths_new[:,1]/widths[:,1]], xlows, xhighs)
                widths = widths_new
                logger('Info: widths of the traces have been updated')
            dummy, sci_tr_poly, widths, xlows, xhighs = proc_gui.remove_adjust_orders_UI( scale_image_plot(im_trace1, 'log10'), sci_tr_poly, xlows, xhighs, widths, userinput=params['GUI'], do_adj=True)
            printresults = True
        if printresults:
            # save parameters of the polynoms into a fitsfile (from Neil)
            save_fits_width(sci_tr_poly, xlows, xhighs, widths, params['master_trace_sci_filename'])
            # Produce some useful statistics
            plot_traces_over_image(scale_image_plot(im_trace1, 'log10'), params['logging_traces_im'].replace('*', 'master_trace1'), sci_tr_poly, xlows, xhighs, widths)
            data = np.insert(widths, 0, list(range(len(sci_tr_poly))), axis=1)       # order, left, right, gausswidth
            positio, pctlwidth = [], []
            for order in range(sci_tr_poly.shape[0]):                      # For the central data
                xarr = list(range(xlows[order],xhighs[order]))
                positio.append(np.polyval(sci_tr_poly[order, 0, 1:], int(im_trace1.shape[0]/2) - sci_tr_poly[order, 0, 0]))
                pctlwidth.append( np.median ( np.polyval(sci_tr_poly[order, 2, 1:], xarr - sci_tr_poly[order, 2, 0]) - np.polyval(sci_tr_poly[order, 1, 1:], xarr - sci_tr_poly[order, 1, 0]) ) )   # median pctlwidth
            data = np.append(data, np.array(pctlwidth)[:,None],axis=1)
            data = np.append(data, np.array(positio)[:,None],axis=1)
            data = np.append(data, np.array(xlows)[:,None],axis=1)
            data = np.append(data, np.array(xhighs)[:,None],axis=1)
            printarrayformat = ['%1.1i','%3.1f', '%3.1f', '%4.2f\t','%4.2f\t', '%1.1i','%1.1i','%1.1i']
            logger('\t\torder\tleft\tright\tgausswidth\tpctlwidth\tpositio\tmin_tr\tmax_tr\t(positio: position of the trace at the center of the image, pctlwidth: median full width of the trace at {0}% of maximum)'\
                      .format(params['width_percentile']),printarrayformat=printarrayformat, printarray=data, params=params)
            # Do a bisector analysis: plot position as fuction of flux
            bisector_measurements_orders(im_trace1,  params.get('logging_traces_bisector', params['logging_path']+'bisector_traces1.png'), sci_tr_poly, xlows, xhighs, widths)
                
    # Create the file for the calibration orders, if it doesn't exist
    if os.path.isfile(params['master_trace_cal_filename']):
        logger('Info: Trace for calibration orders already exists: {0}'.format(params['master_trace_cal_filename']), params=params)
        cal_tr_poly, axlows, axhighs, awidths = read_fits_width(params['master_trace_cal_filename'])
    else:
        # use im_trace2 for automatic solution
        shifts = arc_shift(params, im_trace2, sci_tr_poly, xlows, xhighs, widths)
        # update the sci_tr_poly parameters and create the cal_tr_poly
        cal_tr_poly = []
        for order in range(sci_tr_poly.shape[0]):
            new_pfit = []
            for dataset in range(sci_tr_poly.shape[1]):
                new_pfit.append(list(sci_tr_poly[order, dataset, :-1]) + [sci_tr_poly[order, dataset, -1]+shifts[order]] )
            cal_tr_poly.append(new_pfit)
        cal_tr_poly, awidths, axlows, axhighs = np.array(cal_tr_poly), copy.deepcopy(widths), copy.deepcopy(xlows), copy.deepcopy(xhighs)
        
        # check the shift between the original solution and arc using a GUI
        if params['arcshift_side'] != 0:
            dummy, cal_tr_poly, awidths, axlows, axhighs = proc_gui.remove_adjust_orders_UI(scale_image_plot(im_trace2, 'log10'), cal_tr_poly, axlows, axhighs, awidths, shift=0, userinput=params['GUI'], do_adj=True, do_shft=True)
        
        # save parameters of the polynoms into a fitsfile (from Neil)
        save_fits_width(cal_tr_poly, axlows, axhighs, awidths, params['master_trace_cal_filename'])
        plot_traces_over_image(scale_image_plot(im_trace2, 'log10'), params['logging_arctraces_im'], cal_tr_poly, axlows, axhighs, awidths)
    
    # Catch the problem, when the script re-runs with different settings and therefore the number of orders changes.
    if cal_tr_poly.shape[0] != sci_tr_poly.shape[0]:
        logger('Error: The number of traces for the science fiber and for the calibration fiber do not match. '+\
               'Please remove eighter {0} or {1} and re-run the script in order to solve.'.format(params['master_trace_cal_filename'], params['master_trace_sci_filename']), params=params)
    
    calimages['sci_trace'] = copy.deepcopy( [sci_tr_poly, xlows, xhighs, widths] )      # Apertures might be shifted before extraction -> this would also affect the localbackground
    calimages['cal_trace'] = copy.deepcopy( [cal_tr_poly, axlows, axhighs, awidths] )
    
    # Create the wavelength solutions
    calibs = [ ['_sci','cal1', 'sci', ''], ['_cal','cal2','cal', 'arc'] ]
    if params['arcshift_side'] == 0:        # no second solution necessary
        calibs = [ calibs[0] ]
    elif 'cal2_l_rawfiles' not in params.keys():    # no wavelength solution for calibration fiber available
        del calibs[1]
        logger('Warn: No data available to create the wavelength solution in the calibration fiber', params=params)
    if 'cal1_l_rawfiles' not in params.keys():    # no wavelength solution for science fiber available, also check for single fiber spectrographs
        del calibs[0]
        logger('Warn: No data available to create the wavelength solution in the science fiber. This means your spectra probably might not be well wavelength defined (if there is an offset between science and calibration fiber).', params=params)
    if len(calibs) == 0 and params['original_master_wavelensolution_filename'].lower() != 'pseudo':
        params['original_master_wavelensolution_filename'] = 'pseudo'
        logger('Warn: No data available to create a wavelength solution, therefore will create a pseudo solution', params=params)
    if params['original_master_wavelensolution_filename'].lower() == 'pseudo':              # Pseudo solution
        logger('Warning: Using a pseudo solution for the wavelength (1 step per px)', params=params)
        calimages['wave_sol_dict_sci'] = create_pseudo_wavelength_solution(sci_tr_poly.shape[0])
        calimages['wave_sol_dict_cal'] = calimages['wave_sol_dict_sci']
        calimages['wave_sols_sci'] = None
        calibs = []
    prev = ''
    done_sci = False
    for calib in calibs:     # cal2 is the normal wavelength solution in the calibration fiber
        # first entry is the standard wavelength solution (science fiber)
        if 'master_wavelensolution'+calib[0]+'_filename' not in params.keys():
            params['master_wavelensolution'+calib[0]+'_filename'] = params['master_wavelensolution_filename'].replace('.fit',calib[0]+'.fit')
        for pngparam in ['logging_wavelength_solution_form', 'logging_em_lines_gauss_width_form', 'logging_arc_line_identification_residuals',
                         'logging_arc_line_identification_positions', 'logging_arc_line_identification_residuals_hist', 
                         'logging_em_lines_bisector', 'logging_resolution_form']:
            params[pngparam] = params[pngparam].replace(prev+'.png','')+calib[0]+'.png'
        for pdfparam in ['logging_arc_line_identification_spectrum']:
            params[pdfparam] = params[pdfparam].replace(prev+'.pdf','')+calib[0]+'.pdf'
        for txtparam in ['logging_identified_arc_lines']:
            params[txtparam] = params[txtparam].replace(prev+'.txt','')+calib[0]+'.txt'
        
        # Create files once required, but first try to load them
        im_cal_l, im_arclhead = create_image_general(params, calib[1]+'_l')         # This file is necessary later
        if os.path.isfile(params['master_wavelensolution'+calib[0]+'_filename']) :
            logger('Info: wavelength solution already exists: {0}'.format(params['master_wavelensolution'+calib[0]+'_filename']), params=params)
            wave_sol_dict = read_wavelength_solution_from_fits(params['master_wavelensolution'+calib[0]+'_filename'])
        else:           # Need to create the wavelength solution
            # Identify the Emission lines
            """ needs to replace: cal_l_spec, cal_l_gpm, cal_s_spec, cal_s_gpm, im_cal_l, im_arclhead
            for exptime in ['l','s']:
                fname_em_spec = params['path_extraction'] + params.get('master_{0}_{1}_filename'.format(calib[1], exptime), 'master_emission_{0}_{1}.fits'.format(calib[1], exptime)).rsplit(os.sep,1)[-1]
                if os.path.isfile(fname_em_spec) and exptime == 's':      # 'l' file is necessary later  
                    calimages['im_arc_{0}_{1}'.format(calib[2], exptime)], calimages['im_head_arc_{0}_{1}'.format(calib[2], exptime)] = [], []
                elif os.path.isfile(fname_em_spec) and exptime == 'l':    # 'l' file is necessary later 
                    calimages['im_arc_{0}_{1}'.format(calib[2], exptime)], calimages['im_head_arc_{0}_{1}'.format(calib[2], exptime)] = create_image_general(params, calib[1]+'_'+exptime)
                else:                               # Needs to be extracted again
                    calimages['im_arc_{0}_{1}'.format(calib[2], exptime)], calimages['im_head_arc_{0}_{1}'.format(calib[2], exptime)] = create_image_general(params, calib[1]+'_'+exptime)
                    if 'shift_{0}_l'.format(calib[2]) not in calimages.keys():
                        calimages['shift_{0}_l'.format(calib[2])], calimages['im_head_arc_{0}_l'.format(calib[2])] = find_shift_images(params, 
                                                    calimages['im_arc_{0}_{1}'.format(calib[2], exptime)], im_trace1, calimages['sci_trace'], 0, 
                                                    cal_tr_poly, extract=True, im_head=calimages['im_head_arc_{0}_l'.format(calib[2])])  # that needs to be science traces
                [calimages['em_spec_{0}_{1}'.format(calib[2], exptime)], calimages['em_gpm_{0}_{1}'.format(calib[2], exptime)]] = read_create_spec(params, fname_em_spec, im_cal_s, im_head, calimages['{0}_trace'.format(calib[2])], 
                                                           params['{0}extraction_width_multiplier'.format(calib[3])], calimages.get('shift_{0}_l'.format(calib[2]),0))"""
            
            fname_em_spec = params['path_extraction'] + params.get('master_{0}_l_filename'.format(calib[1]), 'master_emission_{0}_long.fits'.format(calib[1])).rsplit(os.sep,1)[-1]
            if 'shift_{0}_l'.format(calib[2]) not in calimages.keys() and not os.path.isfile(fname_em_spec):
                calimages['shift_{0}_l'.format(calib[2])], im_arclhead = find_shift_images(params, im_cal_l, im_trace1, calimages['sci_trace'], 0, 
                                                                            cal_tr_poly, extract=True, im_head=im_arclhead, max_allowed=[-1,1])        # that needs to be science traces, We don't expect more than +-1 pixel
                 
            [cal_l_spec, cal_l_gpm] = read_create_spec(params, fname_em_spec, im_cal_l, im_arclhead, calimages['{0}_trace'.format(calib[2])], 
                                                       params['{0}extraction_width_multiplier'.format(calib[3])], calimages.get('shift_{0}_l'.format(calib[2]),0) )
            # gpm: good-pixel-mask: 1 is good, below 1 not so good
            fname_em_spec = params['path_extraction'] + params.get('master_{0}_s_filename'.format(calib[1]), 'master_emission_{0}_short.fits'.format(calib[1])).rsplit(os.sep,1)[-1]
            if os.path.isfile(fname_em_spec):
                im_cal_s, im_head = [], []
            else:                               # Needs to be extracted again
                im_cal_s, im_head = create_image_general(params, calib[1]+'_s')
                if 'shift_{0}_l'.format(calib[2]) not in calimages.keys():
                    calimages['shift_{0}_l'.format(calib[2])], im_arclhead = find_shift_images(params, im_cal_l, im_trace1, calimages['sci_trace'], 0,
                                                                                               cal_tr_poly, extract=True, im_head=im_arclhead)  # that needs to be science traces
            [cal_s_spec, cal_s_gpm] = read_create_spec(params, fname_em_spec, im_cal_s, im_head, calimages['{0}_trace'.format(calib[2])], 
                                                       params['{0}extraction_width_multiplier'.format(calib[3])], calimages.get('shift_{0}_l'.format(calib[2]),0))
            
            fname_arclines = params['logging_found_arc_lines'].replace('.txt','')+calib[0]+'.txt'
            if os.path.isfile(fname_arclines):
                logger('Info: List of the identified emission lines already exists. Using the information from file: {0}'.format(fname_arclines ), params=params)
                arc_lines_px_txt = read_text_file(fname_arclines, no_empty_lines=True)              # list of strings, first entry is header 
                arc_lines_px = np.array( convert_readfile(arc_lines_px_txt[1:], [int, float, float, float], delimiter='\t', replaces=[os.linesep,' ']) )
            else:       # Need to identify the lines
                
                arc_lines_px = identify_emission_lines(params, cal_l_spec, cal_s_spec, cal_l_gpm, cal_s_gpm)
                logger('Info: Identified {0} lines in the emission line spectrum. These lines are stored in file {1}'.format(len(arc_lines_px), fname_arclines ), params=params)
                printarrayformat = ['%1.1i', '%7.5f', '%4.2f', '%3.1f']     # needs to be 7.5f for position ([1]) as before 20210413 the lines are only precise to 10milli-pixel -> wavelengths in reflines are only that precise -> wavelenth shift can't be more precise than that 
                logger('order\tpixel\twidth\theight of the line', show=False, printarrayformat=printarrayformat, printarray=arc_lines_px, logfile=fname_arclines, params=params)
            
            # After the preparation start the creation of the wavelength solution
            if calib[1] == 'cal1' or ( calib[1] == 'cal2' and not done_sci) :          # for the science fiber, runs first
                limit_search_space = False
                if os.path.isfile(params['original_master_wavelensolution_filename']) == False:                                             # Create a new solution
                    params['px_to_wavelength_file'] = params.get('px_to_wavelength_file', 'pixel_to_wavelength.txt')                            # Backwards compatible
                    wave_sol_dict = proc_gui.create_new_wavelength_UI(params, cal_l_spec, cal_s_spec, arc_lines_px, reference_lines_dict)
                    save_wavelength_solution_to_fits(wave_sol_dict, params['original_master_wavelensolution_filename'])               # For a new wavelength solution
                    plot_wavelength_solution_form(params['logging_wavelength_solution_form'].replace('.png','')+'_manual.png', calimages['{0}_trace'.format(calib[2])], wave_sol_dict['wavesol'])
                    plot_wavelength_solution_spectrum(params, cal_l_spec, cal_s_spec,params['logging_arc_line_identification_spectrum'].replace('.pdf','')+'_manual.pdf',
                                                      wave_sol_dict['wavesol'], wave_sol_dict['reflines'],
                                                      reference_lines_dict['reference_catalog'][0], reference_lines_dict['reference_names'][0], plot_log=True)
                    limit_search_space = True
                elif 0 < params['original_master_wavelensolution_filename'].find('master_wavelength_manual') < 5:
                    limit_search_space = True
                if limit_search_space:      # Limit the search space if solution from GUI or manual solution as indicated by filename
                    params['order_offset'] = [0,0]
                    params['px_offset'] = [-10,10,2]
                    params['px_offset_order'] = [-0.2,0.2,0.1]
                wave_sol_ori_dict = read_wavelength_solution_from_fits(params['original_master_wavelensolution_filename'])
                 # Find the new wavelength solution
                wave_sol_dict = adjust_wavelength_solution(params, cal_l_spec, cal_s_spec, arc_lines_px, wave_sol_ori_dict['wavesol'], 
                                                           wave_sol_ori_dict['reflines'], params['original_master_wavelensolution_filename'],
                                                           reference_lines_dict, calimages['{0}_trace'.format(calib[2])], show_res=params['GUI'], search_order_offset=True)
                done_sci = True
            else:                           # Calibration fiber, runs second
                params['order_offset'] = [0,0]
                params['px_offset'] = [-60,60,10]
                params['px_offset_order'] = [-0.2,0.2,0.1]
                wave_sol_dict = adjust_wavelength_solution(params, cal_l_spec, cal_s_spec, arc_lines_px, wave_sol_dict['wavesol'], 
                                                           wave_sol_dict['reflines'], 'Science fiber wavelength solution',
                                                           reference_lines_dict, calimages['{0}_trace'.format(calib[2])], show_res=params['GUI'], search_order_offset=False)
            save_wavelength_solution_to_fits(wave_sol_dict, params['master_wavelensolution'+calib[0]+'_filename'])
            plot_wavelength_solution_form(params['logging_wavelength_solution_form'], calimages['{0}_trace'.format(calib[2])], wave_sol_dict['wavesol'])
            plot_wavelength_solution_spectrum(params, cal_l_spec, cal_s_spec, params['logging_arc_line_identification_spectrum'], wave_sol_dict['wavesol'], wave_sol_dict['reflines'], 
                                              reference_lines_dict['reference_catalog'][0], reference_lines_dict['reference_names'][0], plot_log=True)
            plot_wavelength_solution_image(im_cal_l, params['logging_arc_line_identification_positions'], calimages['{0}_trace'.format(calib[2])], 
                                           wave_sol_dict['wavesol'], wave_sol_dict['reflines'], reference_lines_dict['reference_catalog'][0])
        im_head, obsdate_midexp, obsdate_mid_float, jd_midexp = get_obsdate(params, im_arclhead)    # in UTC, mid of the exposure
        name = params['master_wavelensolution'+calib[0]+'_filename'].rsplit(os.sep,1)[-1].replace('.fits','')
        calimages['wave_sols_'+calib[2]] = [ [copy.deepcopy( wave_sol_dict['wavesol'] ), copy.deepcopy( wave_sol_dict['reflines'] ), jd_midexp, name] ]          # Create a list of all wavelength solutions
        calimages['wave_sol_dict_'+calib[2]] = copy.deepcopy( wave_sol_dict )
        add_text_to_file('{0}\t{1}\t{2}\t{3}'.format(jd_midexp, 0, 0, calib[2]), params['master_wavelengths_shift_filename'], warn_missing_file=False )
        if params['arcshift_side'] == 0:
            break
        prev = calib[0]
    # calimages['wave_sol_dict_*'] is a dictonary with the keys: wavesol, wavesol2d, reflines, linestat; e.g. check adjust_wavelength_solution
    # calimages['wave_sols_*'] is a list containing lists with the entries: wave_sol_dict['wavesol'], wave_sol_dict['reflines'], jd_midexp, name
    if 'wave_sol_dict_sci' in calimages.keys() and 'wave_sol_dict_cal' in calimages.keys():
        params['two_solutions'] = True                  # Needed for create_wavelengths_from_solution, as shift needs to be applied
        # Catch the problem, when the script re-runs with different settings and therefore the number of orders changes.
        if calimages['wave_sol_dict_sci']['wavesol'].shape[0] != calimages['wave_sol_dict_cal']['wavesol'].shape[0]:
            logger('Error: The number of traces of the wavelength solutions for science ({0}) and calibration ({1}) fibers do not match. Please remove one of the solutions and re-run the script in order to solve.'\
                    .format(calimages['wave_sol_dict_sci']['wavesol'].shape[0], calimages['wave_sol_dict_cal']['wavesol'].shape[0]), params=params)
    else:
        params['two_solutions'] = False
    if 'wave_sol_dict_sci' not in calimages.keys():                             # Use the calibration fiber solution for the science solution
        calimages['wave_sol_dict_sci'] = copy.deepcopy(calimages['wave_sol_dict_cal'])
        calimages['wave_sols_sci'] = copy.deepcopy(calimages['wave_sols_cal'])
        if len(calimages['wave_sols_sci']) != 1:
            print("253 hiflex.py: Ronny did not expect len(calimages['wave_sols_sci']) to be {0}".format(len(calimages['wave_sols_sci']) ))
        calimages['wave_sol_dict_sci']['wavesol'][:,1] += params['pxshift_between_wavesolutions']       # + checked 20210128
        calimages['wave_sols_sci'][0][0][:,1] += params['pxshift_between_wavesolutions']                # + checked 20210128
        logger('Info: Applied a shift of {0} px from parameter pxshift_between_wavesolutions to convert calibration fiber wavelength solution to science fiber wavelength solution.'.format(params['pxshift_between_wavesolutions']), params=params)
    # Catch the problem, when the script re-runs with different settings and therefore the number of orders changes.
    if calimages['wave_sol_dict_sci']['wavesol'].shape[0] != calimages['sci_trace'][0].shape[0]:
        #print('im_trace1.shape', im_trace1.shape)
        logger('Error: The number of traces for extraction and for the wavelength calibration do not match. Please remove eighter {0} ({2}) or {1} ({3}) and re-run the script in order to solve.'\
                    .format(params['master_trace_sci_filename'], params['master_wavelensolution_filename'], calimages['sci_trace'][0].shape[0], calimages['wave_sol_dict_sci']['wavesol'].shape[0]), params=params)
    
    """# Find the conversion between the two solutions (px-shift = f(wavelength) (but not really as seen on 20200813 229? images)
    if params['two_solutions']:
        params['extract_wavecal'] = True
        if calimages['wave_sol_dict_cal']['wavesol'].shape[0] != calimages['wave_sol_dict_sci']['wavesol'].shape[0]:
            logger('Error: The number of traces for the science ({0}) and calibration ({1}) wavelength solution differ. Please delete the wrong, old file ({2} or {3})'.format(\
                        calimages['wave_sol_sci'].shape[0], calimages['wave_sol_cal'].shape[0], params['master_wavelensolution_sci_filename'], params['master_wavelensolution_filename'] ), params=params)
        shift, shift_err = find_shift_between_wavelength_solutions(params, calimages['wave_sol_dict_cal']['wavesol'], calimages['wave_sol_dict_cal']['reflines'], 
                                                                   calimages['wave_sol_dict_sci']['wavesol'], calimages['wave_sol_dict_sci']['reflines'], 
                                                                   np.zeros((calimages['wave_sol_dict_cal']['wavesol'].shape[0],im_trace1.shape[0])), ['calibration fiber','science fiber'] )
        add_text_to_file('{0}\t{1}\t{2}\t{3}'.format(0, round(shift,4), round(shift_err,4), 'sci-cal'), params['master_wavelengths_shift_filename'], warn_missing_file=False )                       
        # shift is positive if lines in science are right of lines in calibration
        "#""if params['wavelength_solution_type'] == 'sci-fiber':
            steps = ['cal-fiber']           # first do it for the not wavelength solution
        else:
            steps = ['sci-fiber']           # first do it for the not wavelength solution
        steps.append(params['wavelength_solution_type'])
        "#""
        if params['wavelength_solution_type'] == 'cal-fiber':           # Calibration fiber
            im_name = 'long_exposure_emision_lamp_in_science_fiber'
            wttype = 'sci'  # opposite
        else:                                                           # Science fiber
            im_name = 'long_exposure_emision_lamp_in_calibration_fiber'
            shift = -shift                          # reverse, as opposite to find_shift_between_wavelength_solutions 
            wttype = 'cal'  # opposite
        aspectra = calimages['arc_l_spec_'+wttype]
        im_arclhead = calimages['arc_l_head_'+wttype]
        im_head, obsdate_midexp, obsdate_mid_float, jd_midexp = get_obsdate(params, im_arclhead)               # in UTC, mid of the exposure

        dummy_shift_wavesoln, master_shift, im_head = shift_wavelength_solution(params, aspectra, wave_sol_dict, reference_lines_dict, 
                                                              xlows, xhighs, obsdate_mid_float, jd_midexp, sci_tr_poly, cal_tr_poly, im_name, maxshift=max(3,2*shift_err), in_shift=shift, im_head=im_head )
        # master shift gives the shift from the wavelength_solution to the aspectra
        params['pxshift_between_wavesolutions'] = master_shift
        # print("0:params['pxshift_between_wavesolutions']", params['pxshift_between_wavesolutions'])
        "#"" # Comparing the shifted solution with the original wavelength solution -> shifts in either direction, depending on polynomial, maximum shift: -> 0.025 \AA shift at 4880 -> 1.5 km/s
        wsci = create_wavelengths_from_solution(calimages['wave_sol_sci'], calimages['arc_l_spec_cal'])
        wcal = create_wavelengths_from_solution(calimages['wave_sol_cal'], calimages['arc_l_spec_cal'])
        wshift = create_wavelengths_from_solution(dummy_shift_wavesoln, calimages['arc_l_spec_cal'])
        save_multispec([wsci, calimages['arc_l_spec_sci']], 'fib1_wsci.fits', im_arclhead, bitpix=params['extracted_bitpix'])
        save_multispec([wsci, calimages['arc_l_spec_cal']], 'fib2_wsci.fits', im_arclhead, bitpix=params['extracted_bitpix'])
        save_multispec([wcal, calimages['arc_l_spec_sci']], 'fib1_wcal.fits', im_arclhead, bitpix=params['extracted_bitpix'])
        save_multispec([wcal, calimages['arc_l_spec_cal']], 'fib2_wcal.fits', im_arclhead, bitpix=params['extracted_bitpix'])
        save_multispec([wshift, calimages['arc_l_spec_sci']], 'fib1_wshift.fits', im_arclhead, bitpix=params['extracted_bitpix'])
        save_multispec([wshift, calimages['arc_l_spec_cal']], 'fib2_wshift.fits', im_arclhead, bitpix=params['extracted_bitpix'])
        #"#""
    else:                           # read the shift to the science fiber
        all_shifts = read_text_file(params['master_wavelengths_shift_filename'], no_empty_lines=True)
        all_shifts = convert_readfile(all_shifts, [float, float, float, str], delimiter='\t', replaces=['\n']+[['  ',' ']]*20)       # jd_midexp, shift_avg, shift_std, can contain duplicate jd_midexp (last one is the reliable one)
        for entry in all_shifts[::-1]:
            if entry[3] == 'sci-cal': 
                params['pxshift_between_wavesolutions'] = - entry[1]
                break"""
    
    params['extract_wavecal'] = False
    
    # Extract the flat spectrum and normalise it
    if os.path.isfile(params['master_blaze_spec_norm_filename']) :
        logger('Info: Normalised blaze already exists: {0}'.format(params['master_blaze_spec_norm_filename']), params=params)
        # The file is read later on purpose
    else:
        create_blaze_norm(params, im_trace1, calimages['sci_trace'], calimages['cal_trace'], calimages['wave_sol_dict_sci'], reference_lines_dict)
    
    flat_spec_norm = np.array(fits.getdata(params['master_blaze_spec_norm_filename']))              # read it again, as the file is different than the data above
    # Catch the problem, when the script re-runs with different settings and therefore the number of orders changes.
    if flat_spec_norm.shape[1] != calimages['wave_sol_dict_sci']['wavesol'].shape[0]:
        #print('im_trace1.shape', im_trace1.shape)
        logger('Error: The number of traces in the blaze of the blaze correction and for the wavelength calibration do not match. Please remove {0} ({1} instead of expected {2}) and re-run the script in order to solve.'\
                        .format(params['master_blaze_spec_norm_filename'], flat_spec_norm.shape[1], calimages['wave_sol_dict_sci']['wavesol'].shape[0]), params=params)
        
    # Use only for extraction -> move down to after the wavelenghs solution have been loaded
    remove_orders, keep_orders = remove_orders_low_flux(params, flat_spec_norm)
    calimages['flat_spec_norm'] = copy.deepcopy( flat_spec_norm )
    if len(remove_orders) > 0:
        remove_orders_from_calimages(params, calimages, keep_orders)            # remove the bad orders
    """calimages['flat_spec_norm'] = copy.deepcopy( flat_spec_norm )
    calimages['wave_sol_dict'] = copy.deepcopy( wave_sol_dict )
    calimages['wavelength_solution'] = copy.deepcopy( wave_sol_dict['wavesol'] )    # Needed for templates later"""
    
    if not params['started_from_p3']:       # if started from python3 environment, this won't be necessary
        logger('Info: Finished routines for a new night of data. Now science data can be extracted. Please check before the output in the logging directory {1}: Are all orders identified correctly for science and calibration fiber, are the correct emission lines identified for the wavelength solution?{0}'.format(os.linesep, params['logging_path']), params=params)
        
        obj_names = []
        extractions = []
        wavelengthcals_cal, wavelengthcals_sci = [], []
        for entry in params.keys():
            if entry.find('extract') >= 0 and entry.find('_rawfiles') >= 0:
                extractions.append(entry.replace('_rawfiles',''))
            if entry.find('waveoffsetsci_') >= 0 and entry.find('_rawfiles') >= 0:          # ThAr in science fiber
                wavelengthcals_sci.append(entry.replace('_rawfiles',''))
            elif entry.find('waveoffsetcal_') >= 0 and entry.find('_rawfiles') >= 0:         # ThAr in calibration fiber
                wavelengthcals_cal.append(entry.replace('_rawfiles',''))
        
        if ( params['arcshift_side'] == 0 or params['two_solutions'] or True ) and len(wavelengthcals_cal)+len(wavelengthcals_sci) > 0:         # or True: added 20201016, do it in any case, as there might be ThAr only sometimes in cal fiber
            def wavecal_multicore(parameter):
                wavelengthcal, fib, im_name_full, im_name_wc, im_name, realrun = parameter
                params['calibs'] = params[wavelengthcal+'_calibs_create']
                im, im_head = read_file_calibration(params, im_name_full, realrun=realrun)
                if realrun:
                    extraction_wavelengthcal(params, im, im_name_wc, im_head, calimages['sci_trace'], calimages['cal_trace'],
                                                        calimages['wave_sol_dict_'+fib], reference_lines_dict, im_trace1, im_name_wc)
            params['extract_wavecal'] = True                                                                                # necessary for shift_wavelength_solution so the shift is stored in a file
            all_wavelengthcals = []
            for [wavelengthcals, fib] in [ [wavelengthcals_cal,'cal'], [wavelengthcals_sci,'sci'] ]:
                if 'wave_sol_dict_'+fib not in calimages.keys():                    # Ignore cal if it's a single fiber spectrograph
                    continue
                for wavelengthcal in wavelengthcals:
                    # !!! Posible improvement: combine a few files if they are taken close to each other
                    for im_name_full in params[wavelengthcal+'_rawfiles']:
                        im_name = im_name_full.rsplit(os.sep)
                        im_name = im_name[-1].rsplit('.',1)         # remove the file ending
                        im_name = im_name[0]
                        im_name_wc = im_name+'_wave'+fib
                        if os.path.isfile(params['path_extraction']+im_name_wc+'.fits'):
                            logger('Info: File {0} was already processed for the calibration of the wavelength solution. If you want to extract again, please delete {1}{0}.fits'.format(im_name_wc, params['path_extraction']), params=params)
                        else:
                            wavecal_multicore([ wavelengthcal, fib, im_name_full, im_name_wc, im_name, False ])      # Try run to get the calibration data
                            all_wavelengthcals.append([ wavelengthcal, fib, im_name_full, im_name_wc, im_name, True ])
            if params['use_cores'] > 1  and len(all_wavelengthcals) > 1:
                logger('Info: Starting to extract wavelength calibrations using multiprocessing on {0} cores, hence output will be for several files in parallel.'.format(params['use_cores']), params=params)
                sort_wavelengthcals = sort_for_multiproc_map(all_wavelengthcals, params['use_cores'])
                p = multiprocessing.Pool(params['use_cores'])
                p.map(wavecal_multicore, sort_wavelengthcals)
                p.terminate()
                p.join()
            elif len(all_wavelengthcals) > 0:
                logger('Info: Starting to extract wavelength calibrations', params=params)
                for all_wavelengthcal in all_wavelengthcals:
                    wavecal_multicore(all_wavelengthcal)
        
        # Load all the wavelength solutions
        last_shape = None
        for wavesol in sorted(os.listdir(params['path_extraction'])):
            if not wavesol.endswith(".fits"):      continue
            if wavesol.find('_wavesol_') == -1:    continue
            fib = wavesol.rsplit('_',1)[-1].replace('.fits','')
            jd_midexp = float( wavesol.split('_wavesol_')[-1].replace('_{0}.fits'.format(fib),'') )
            name = wavesol.split('_wavesol_')[0]
            wavelength_solution_shift_dict = read_wavelength_solution_from_fits(params['path_extraction']+wavesol)
            if wavelength_solution_shift_dict['wavesol'].shape[1] != last_shape and last_shape is not None:
                logger('Warn: This wavelength solution has a different number of freedoms {0} than the previous solution {1}.'.format(wavelength_solution_shift_dict['wavesol'].shape[1], last_shape), params=params)
            last_shape = wavelength_solution_shift_dict['wavesol'].shape[1]
            calimages['wave_sols_'+fib].append( [wavelength_solution_shift_dict['wavesol'], wavelength_solution_shift_dict['reflines'], jd_midexp, name])
        
        # Read the list of files
        file_list = read_text_file(params['raw_data_file_list'], no_empty_lines=True)
        file_list = convert_readfile(file_list, [str, str, str, float, ['%Y-%m-%dT%H:%M:%S', float], str, str], delimiter='\t', replaces=['\n',' ', os.linesep])
            
        params['extract_wavecal'] = False
        if len(extractions) == 0:                                               # no extractions to do
            logger('Warn: Nothing to extract.', params=params)
            #header_results_to_texfile(params)           # Save the results from the header in a logfile
            #exit(0)
        
        def extraction_multicore(all_extractions):
            [extraction, im_name_full, im_name, realrun, obnames] = all_extractions
            if os.path.isfile(params['path_extraction']+im_name+'.fits'):
                logger('Info: File {0} was already processed. If you want to extract again, please delete {1}{0}.fits'.format(im_name, params['path_extraction']), params=params)
                return ''
            if  extraction.find('extract_combine') == -1:     # Single file extraction
                params['calibs'] = params[extraction+'_calibs_create']
                im, im_head = read_file_calibration(params, im_name_full, realrun=realrun)   
            else:                                       # Combine files before extraction
                im, im_head = create_image_general(params, extraction, realrun=realrun)
            if realrun:
                obj_name = extraction_steps(params, im, im_name, obnames, im_head, calimages['sci_trace'], calimages['cal_trace'], 
                                        calimages.get('wave_sol_dict_cal',calimages['wave_sol_dict_sci']), reference_lines_dict, calimages['flat_spec_norm'], im_trace1)
                return obj_name.lower()

        
        all_extractions, list_im_name = [], []
        for extraction in extractions:
            if extraction.find('extract_combine') == -1:     # Single file extraction
                for im_name_full in params[extraction+'_rawfiles']:
                    im_name = im_name_full.rsplit(os.sep)
                    im_name = im_name[-1].rsplit('.',1)         # remove the file ending
                    im_name = im_name[0]
                    if im_name_full in list_im_name:
                        continue
                    if os.path.isfile(params['path_extraction']+im_name+'.fits'):
                        logger('Info: File {0} was already processed. If you want to extract again, please delete {1}{0}.fits'.format(im_name, params['path_extraction']), params=params)
                    else:
                        # Look in the other list for the information about the object - this process will take a while
                        obnames = ['']
                        for ii, entry in enumerate(file_list):
                            if entry[0] == im_name_full:
                                if entry[5] != '':
                                    obnames = entry[5].replace(' ','').split(',')
                                break
                        extraction_multicore([extraction, im_name_full, im_name, False, obnames])        # Dry run to create darks and such
                        all_extractions.append([extraction, im_name_full, im_name, True, obnames])       # Real run
                    list_im_name.append(im_name_full)
            else:                                           # Combine the files before extraction
                if extraction in list_im_name:
                    continue
                if os.path.isfile(params['path_extraction']+extraction+'.fits'):
                    logger('Info: File {0} was already processed. If you want to extract again, please delete {1}{0}.fits'.format(extraction, params['path_extraction']), params=params)
                else:
                    
                    extraction_multicore([extraction, extraction, extraction, False, [''] ])           # Dry run to create darks and such
                    all_extractions.append([extraction, extraction, extraction, True, [''] ])          # Real run
                list_im_name.append(extraction)
        if params['use_cores'] > 1 and len(all_extractions) > 1:
            logger('Info: Starting to extract spectra using multiprocessing on {0} cores, hence output will be for several files in parallel'.format(params['use_cores']), params=params)
            p = multiprocessing.Pool(params['use_cores'])
            #sort_extractions = sort_for_multiproc_map(all_extractions, params['use_cores'])
            obj_names += p.map(extraction_multicore, all_extractions)
            p.terminate()
            p.join()
        elif len(all_extractions) > 0:
            logger('Info: Starting to extract spectra', params=params)
            for all_extraction in all_extractions:
                obj_names.append( extraction_multicore(all_extraction) )
                
        #obj_names = np.unique(obj_names)
        #obj_names = list(obj_names[obj_names != ''])
        
        # Combine files after extraction
        if params['wavelength_scale_resolution'] > 0:
            dwave = params['wavelength_scale_resolution']
            for extraction in extractions:
                if extraction.find('extract_lin_') == 0:
                    if os.path.isfile(params['path_extraction_single']+extraction+'_lin_cont.fits'):
                        logger('Info: File {0} was already created. If you want to extract again, please delete {1}'.format(extraction, params['path_extraction_single']+extraction+'_lin_cont.fits'))
                        continue
                    logger('Info: Combining {0}'.format(extraction), params=params)
                    data, info, text_info = [], [], []
                    for im_name_full in params[extraction+'_rawfiles']:
                        im_name = im_name_full.rsplit(os.sep)
                        im_name = im_name[-1].rsplit('.',1)         # remove the file ending
                        im_name = im_name[0]
                        im_name_full = params['path_extraction_single']+im_name+'_lin_cont.fits'
                        if os.path.isfile(im_name_full):
                            fdata, im_head = read_fits_file(im_name_full)
                            data.append(fdata)
                            info.append([ im_head['HIERARCH HiFLEx JD'], im_head['HIERARCH HiFLEx EXPOSURE'], im_head['HIERARCH HiFLEx fsum_all'] ])
                            text_info.append([ im_name_full, im_name ])
                            logger('Info: Loaded {0}'.format(im_name_full), params=params)
                        else:
                            logger('Warn: Can not find {0}'.format(im_name_full), params=params)
                    info = np.array(info)
                    if info.shape[0] == 0:
                        logger('Warn: Could not find any file for extraction {0}.'.format(extraction), params=params)
                        continue
                    waves = np.unique( np.concatenate( list(entry[0] for entry in data) ).flat, axis=None )        # might be an issue due to float uncertainties
                    waves = np.arange( int(np.nanmin(waves)/dwave) * dwave, np.nanmax(waves)+dwave, dwave )
                    fluxes = np.zeros(( info.shape[0], waves.shape[0] ))*np.nan
                    for ii, entry in enumerate(data):
                        index_s1, index_e1 = np.nanargmin(entry[0,:]), np.nanargmax(entry[0,:])
                        index_s2 = np.argmin(np.abs(waves-entry[0,index_s1]))    # use the first nonnan in the wavelengths of the linear spectrum
                        index_e2 = np.argmin(np.abs(waves-entry[0,index_e1]))
                        if index_e2-index_s2 != index_e1-index_s1:
                            #print('should not happen', index_s1, index_e1, index_s2, index_e2, entry[index_s1:index_s1+2,0], entry[index_e1-2:index_e1,0], waves[index_s2:index_s2+2,0], waves[index_e2-2:index_e2,0]
                            logger('Error: file {0} has not the expected number of data points between wavelengths {1} and {2}. Expected {3}, got {4} data points. Did you change parameter wavelength_scale_resolution during the extraction process. Please delete all extracted files for {5}.'.format(text_info[ii][0], waves[index_s2,0], waves[index_e2,0], index_e2-index_s2, index_e1-index_s1, text_info[ii][1]), params=params)
                        fluxes[ii, index_s2:index_e2+1] = entry[1,index_s1:index_e1+1]
                        #print(ii, index_s2,index_e2,index_s1,index_e1,fluxes[ii, index_s2:index_e2+1] )
                    exptime = np.sum(info[:,1])
                    nans = np.isnan(fluxes)
                    notallnans = ~nans.all(axis=0)
                    flux_comb = np.zeros( waves.shape[0] )*np.nan
                    if extraction.find('extract_lin_sum') == 0:
                        flux_comb[notallnans] = np.nansum(fluxes[:,notallnans], axis=0)
                        jd = np.mean(info[:,0])
                    elif extraction.find('extract_lin_med') == 0:
                        flux_comb[notallnans] = np.nanmedian(fluxes[:,notallnans], axis=0)
                        jd = np.mean(info[:,0])
                    else:
                        weights = info[:,2]**2
                        if np.min(weights) < 0:
                            weights += np.min(weights)
                        jd = np.average(info[:,0], weights=weights)
                        weights = np.tile(weights, (waves.shape[0],1) ).T
                        fluxes[nans] = 0        # replace nans with 0
                        weights[nans] = 0       # replace nans with 0
                        #weights[0,allnans] = 1  # if all weigths are 0 then change one into 1  # not needed anymore due to notallnans
                        flux_comb[notallnans] = np.average(fluxes[:,notallnans], axis=0, weights=weights[:,notallnans])
                    
                    im_head['HIERARCH HiFLEx JD'] = jd
                    im_head['HIERARCH HiFLEx EXPOSURE'] = exptime
                    save_multispec([waves,flux_comb], params['path_extraction_single']+extraction+'_lin_cont', im_head)
                    save_linspec_csv(params['path_extraction_single']+extraction+'_lin_cont.csv', waves, flux_comb)
                    logger('Info: Saved the result in {0}'.format(params['path_extraction_single']+extraction+'_lin_cont.*'), params=params)
        else:
            for extraction in extractions:
                if extraction.find('extract_lin_') == 0:
                    logger('Warn: found settings to combine the linearised extracted files, however parameter wavelength_scale_resolution was set to 0. Please modify wavelength_scale_resolution and delete the files in {0} that you want to be combined. Then run hiflex.py again.'.format(params['path_extraction']), params=params)
        
        logger('')      # To have an empty line
        logger('Info: Finished extraction of the science frames. The extracted {0}*.fits file contains different data in a 3d array in the form: data type, order, and pixel. First data type is the wavelength (barycentric corrected), second is the extracted spectrum, followed by a measure of error. Forth and fifth are the flat corrected spectra and its error. Sixth and sevens are the the continuum normalised spectrum and the S/N in the continuum. Eight is the bad pixel mask, marking data, which is saturated or from bad pixel. The ninth entry is the spectrum of the calibration fiber. The last entry is the wavelength without barycentric correction'.format(params['path_extraction']), params=params)
        header_results_to_texfile(params)           # Save the results from the header in a logfile
        logger(('Info: Will try to do the RV analysis on all files for which the barycentric velocity could be calculated'+\
                ' and for which the filename does NOT start with {0} in the first few letters in a moment').format(params['no_RV_names']), params=params) 
        #time.sleep(2)
        runRV = True
        for arg in sys.argv[1:]:
            if 0 <= arg.lower().find('norv') <= 2:
                runRV = False
                break
        if np.max(calimages['wave_sol_dict_sci']['wavesol'][:,-1]) < 100:
            logger('Info: Using a pseudo wavelength solution -> no RV analysis', params=params)   
        elif not runRV:
            logger('Info: User request to skip RVs')
        else:
            # Run RV analysis that can be run
            files_RV, headers = prepare_for_rv_packages(params)
            kwargs = run_terra_rvs(params)
            kwargs += run_serval_rvs(params)
            run_terra_seral_rvs(params, kwargs)
            # Prepare RV analysis for CERES
            if sys.version_info[0] > 2:             # Try to load a python 2 environment
                anaconda = sys.executable.split('bin{0}python'.format(os.sep))[0]
                if anaconda.lower().find('conda') != -1:
                    log_python2 = 'logfile_python2'
                    if anaconda.lower().find('conda') < anaconda.find('{0}envs{0}'.format(os.sep)):
                        anaconda = anaconda.split('envs{0}'.format(os.sep))[0]
                    cmd  = '__conda_setup="$('
                    cmd += "'{1}bin{0}conda' 'shell.bash' 'hook' 2> /dev/null)".format(os.sep, anaconda)
                    cmd += '" ; if [ $? -eq 0 ]; then eval "$__conda_setup" ; else '
                    cmd += 'if [ -f "{1}etc{0}profile.d{0}conda.sh" ]; then . "{1}etc{0}profile.d{0}conda.sh" ; else export PATH="{1}bin{2}$PATH" ;'.format(os.sep, anaconda, os.pathsep)
                    cmd += ' fi ; fi ; unset __conda_setup ; '
                    cmd += 'conda activate hiflex_p2 ; '
                    cmd += 'python {0}/hiflex.py {1} started_from_p3=True'.format(os.path.dirname(sys.argv[0]), ' '.join(sys.argv[1:]) )
                    logger('Info: Loading a python 2 environment, as CERES requires python 2 and this is a python 3 environment.'+
                           ' The progress of the process can be watched in logfile or {0}'.format(log_python2)+
                           ' Running the following commands to start the other environment:{1}{0}'.format(cmd, os.linesep), params=params)
                    if os.path.isfile(log_python2):
                        os.system('rm {0}'.format(log_python2))
                    with open(log_python2, 'a') as logf:
                        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, shell=True, executable='/bin/bash', stdout=logf, stderr=subprocess.STDOUT)    # shell=True is necessary
                        p.communicate(input=os.linesep.encode())                                       # This waits until the process needs the enter
                        code = log_returncode(p.returncode, 'Please check {0} for the error message. You might want to run the pipeline again in the Python 2 environment.'.format(log_python2))
                    if code > 0:    # It didn't run through:
                        logger('Warn: Loading a python 2 environment has failed or something crashed within the python 2 environment, this will prevent CERES from running. '+\
                               'Please check {1} and that the following line works in a clear terminal (e.g. after running "env -i bash --norc --noprofile"):{2}{2}{0}{2}'.format(cmd, log_python2, os.linesep), params=params)
                else:
                    logger('Warn: CERES might not work. This is a python3 environment, however, it was not started by conda and the script does not know how to get to a python2 environment', params=params)
            else:
                # Run CERES, or if it's a python 2 environment
                run_ceres_rvs(params, files_RV, headers)
            
            # Collect all RV results
            rv_results_to_hiflex(params)
     
            header_results_to_texfile(params)           # Save the results from the header in a logfile 
    else:       # Started from p3
        files_RV, headers = prepare_for_rv_packages(params)     # Won't recreate the folders for TERRA/SERVAL if started_from_p3=True
        run_ceres_rvs(params, files_RV, headers)
        
    logger('Finished', params=params)
    log_params(params)
    

        
    

