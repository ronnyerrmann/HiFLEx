#!/usr/bin/env python

#import numpy as np
#import os
from procedures import *
import subprocess
#import multiprocessing

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
    logger('Info: Starting routines for a new night of data, including: finding or shifting orders, find calibration orders, create wavelength solution, and create normalised blaze function')
    params['path_run'] = os. getcwd()
    params['extract_wavecal'] = False
    params['no_RV_names'] = ['flat', 'tung', 'whili', 'thar', 'th_ar', 'th-ar']
    log_params(params)
    
    # create the median combined files
    im_trace1, im_trace1_head = create_image_general(params, 'trace1')
    if params['arcshift_side'] != 0 or 'trace2_rawfiles' in params.keys():      # calibration spectrum at the same time
        im_trace2, im_trace2_head = create_image_general(params, 'trace2')
    else:                                                                       # no calibration spectrum at the same time
        im_trace2, im_trace2_head = im_trace1, im_trace1_head
    #blazecor, cal2_l, cal2_s: at a later step to know the orders for localbackground -> cont1cal2
    
    # Load the reference catalogue. Note that the Ar I lines are rescaled!
    reference_catalog, reference_names = read_reference_catalog(params)
    
    # Create or read the file with the orders for this night
    if os.path.isfile(params['master_trace_sci_filename']) :
        logger('Info: Using exiting trace solution: {0}'.format(params['master_trace_sci_filename']))
        sci_tr_poly, xlows, xhighs, widths = read_fits_width(params['master_trace_sci_filename'])
    else:
        printresults = False
        # load the original solution
        if os.path.isfile(params['original_master_traces_filename']) :
            sci_tr_poly, xlows, xhighs, widths = read_fits_width(params['original_master_traces_filename'])
            # find the shift between the original solution and the current flat
            shift, widths_new, shift_map, shift_error = shift_orders(im_trace1, params, sci_tr_poly, xlows, xhighs, widths, params['in_shift'])
            # save the map of the shifts
            save_im_fits(params, shift_map, im_trace1_head, params['logging_map_shift_orders'])
        else:
            shift_error = -1
        if shift_error > 1 or shift_error == -1 or abs(shift) > params['maxshift_orders']:
            logger('Warn: The deviation of the shift of the orders seems too big or no previous solution was available, therefore searching for the position of the orders from scratch:')
            sci_tr_poly, xlows, xhighs, widths = find_adjust_trace_orders(params, im_trace1, im_trace1_head)
            printresults = True
        else:
            sci_tr_poly[:,:,-1] += shift                # update the sci_tr_poly parameters
            if params['update_width_orders'] :
                sci_tr_poly, widths = update_tr_poly_width_multiplicate(sci_tr_poly, widths, [widths_new[:,0]/widths[:,0], widths_new[:,1]/widths[:,1]], xlows, xhighs)
                widths = widths_new
                logger('Info: widths of the traces have been updated')
            dummy, sci_tr_poly, widths = remove_adjust_orders_UI( scale_image_plot(im_trace1, 'log10'), sci_tr_poly, xlows, xhighs, widths, userinput=params['GUI'], do_adj=True)
            printresults = True
        if printresults:
            # save parameters of the polynoms into a fitsfile (from Neil)
            save_fits_width(sci_tr_poly, xlows, xhighs, widths, params['master_trace_sci_filename'])
            # Produce some useful statistics
            plot_traces_over_image(im_trace1, params['logging_traces_im'], sci_tr_poly, xlows, xhighs, widths)
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
                      .format(params['width_percentile']),printarrayformat=printarrayformat, printarray=data)
            # Do a bisector analysis: plot position as fuction of flux
            bisector_measurements_orders(im_trace1,  params.get('logging_traces_bisector', params['logging_path']+'bisector_traces1.png'), sci_tr_poly, xlows, xhighs, widths)
        
    """Not really useful
    # Create the background map, if it doesn't exist
    if os.path.isfile(params['background_filename']) :
        logger('Info: Background map already exists: {0}'.format(params['result_path']+params['background_filename']))
    else:
        # create the background map
        if params['GUI'] :
            bck_px, params = bck_px_UI(params, im_trace1, sci_tr_poly, xlows, xhighs, widths, params['background_width_multiplier'][0], params['GUI'])
        else:
            bck_px = find_bck_px(im_trace1, sci_tr_poly, xlows, xhighs, widths, params['background_width_multiplier'][0])
        bad_values = ( im_trace1*bck_px > np.percentile(im_trace1[bck_px==1],95) )
        bck_px[bad_values] = 0
        save_im_fits(params, bck_px, im_trace1_head, params['background_px_filename'])
        save_im_fits(params, bck_px*im_trace1, im_trace1_head, params['logging_orig_for_background'])
        # Create the fitted background map
        #im_bck = find_bck_fit(im_trace1, im_bck_px, params['polynom_bck'], params['GUI'])       #Old
        bck_im = fit_2d_image(im_trace1, params['polynom_bck'][1], params['polynom_bck'][0], w=bck_px)
        save_im_fits(params, bck_im, im_trace1_head, params['background_filename'])"""
        
    # Create the file for the calibration orders, if it doesn't exist
    if os.path.isfile(params['master_trace_cal_filename']) :
        logger('Info: Arc trace solution already exists: {0}'.format(params['master_trace_cal_filename']))
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
        dummy, cal_tr_poly, awidths = remove_adjust_orders_UI((im_trace2), cal_tr_poly, axlows, axhighs, awidths, shift=0, userinput=params['GUI'], do_adj=True, do_shft=True)
        
        # save parameters of the polynoms into a fitsfile (from Neil)
        save_fits_width(cal_tr_poly, axlows, axhighs, awidths, params['master_trace_cal_filename'])
        plot_traces_over_image(im_trace2, params['logging_arctraces_im'], cal_tr_poly, axlows, axhighs, awidths)
    
    # Catch the problem, when the script re-runs with different settings and therefore the number of orders changes.
    if cal_tr_poly.shape[0] != sci_tr_poly.shape[0]:
        logger('Error: The number of traces for the science fiber and for the calibration fiber do not match. Please remove eighter {0} or {1} and re-run the script in order to solve.'.format(params['master_trace_cal_filename'], params['master_trace_sci_filename']))
    
    calimages['sci_trace'] = copy.deepcopy( [sci_tr_poly, xlows, xhighs, widths] )      # Apertures might be shifted before extraction -> this would also affect the localbackground
    calimages['cal_trace'] = copy.deepcopy( [cal_tr_poly, axlows, axhighs, awidths] )
    
    # Do the wavelength solution stuff: Find for calibration fiber (and for science fiber if necessary, in this case some comparison between the solutions is necessay
    params['two_solutions'] = False
    for calib in [ ['','cal2','cal'], ['_sci','cal1', 'sci'] ]:     # cal2 is the normal wavelength solution in the calibration fiber
        # first entry is the standard wavelength solution
        # second entry is for finding the wavelength solution in a bifurcated fiber
        if calib[1] == 'cal1' and ('cal1_l_rawfiles' not in params.keys() or  params['original_master_wavelensolution_filename'].lower() == 'pseudo' or params['arcshift_side'] == 0):
            break                                           # Not set up for bifurcated fiber use or pseudo or only one fiber -> stop after the first step
        elif calib[1] == 'cal1':                            # Update and rename a few parameters
            if 'master_wavelensolution'+calib[0]+'_filename' not in params.keys():
                params['master_wavelensolution'+calib[0]+'_filename'] = params['master_wavelensolution_filename'].replace('.fit',calib[0]+'.fit')
            for pngparam in ['logging_wavelength_solution_form', 'logging_em_lines_gauss_width_form', 'logging_arc_line_identification_residuals', 'logging_arc_line_identification_positions']:
                params[pngparam] = params[pngparam].replace('.png','')+calib[0]+'.png'
            for pdfparam in ['logging_arc_line_identification_spectrum']:
                params[pdfparam] = params[pdfparam].replace('.pdf','')+calib[0]+'.pdf'
        # Create the wavelength solution for the night
        if params['original_master_wavelensolution_filename'].lower() != 'pseudo':                  # Create the master files
            im_cal_l, im_arclhead = create_image_general(params, calib[1]+'_l')
            shift_cal, im_arclhead = find_shift_images(params, im_cal_l, im_trace1, sci_tr_poly, xlows, xhighs, widths, 0, cal_tr_poly, extract=True, im_head=im_arclhead)
            if calib[1] == 'cal2':                          # Calibration fibre
                shift_cal2 = copy.copy(shift_cal)
                cal_l_spec, good_px_mask_l, extr_width = extract_orders(params, im_cal_l, cal_tr_poly, axlows, axhighs, awidths, params['arcextraction_width_multiplier'], offset=shift_cal2, header=im_arclhead)
            else:                                           # Science fibre
                shift_cal1 = copy.copy(shift_cal)
                cal_l_spec, good_px_mask_l, extr_width = extract_orders(params, im_cal_l, sci_tr_poly, xlows, xhighs, widths, params['extraction_width_multiplier'], offset=shift_cal1, header=im_arclhead)
        if os.path.isfile(params['master_wavelensolution'+calib[0]+'_filename']) :
            logger('Info: wavelength solution already exists: {0}'.format(params['master_wavelensolution'+calib[0]+'_filename']))
            wavelength_solution, wavelength_solution_arclines = read_wavelength_solution_from_fits(params['master_wavelensolution'+calib[0]+'_filename'])
        elif params['original_master_wavelensolution_filename'].lower() == 'pseudo':
            logger('Warning: Using a pseudo solution for the wavelength (1 step per px)')
            wavelength_solution, wavelength_solution_arclines = create_pseudo_wavelength_solution(sci_tr_poly.shape[0])
        else:
            im_cal_s, im_arcshead = create_image_general(params, calib[1]+'_s')
            if calib[1] == 'cal2':                              # Calibration fibre
                cal_s_spec, good_px_mask_s, extr_width = extract_orders(params, im_cal_s, cal_tr_poly, axlows, axhighs, awidths, params['arcextraction_width_multiplier'], offset=shift_cal2, header=im_arcshead)
            else:                                               # Science fibre
                cal_s_spec, good_px_mask_s, extr_width = extract_orders(params, im_cal_s, sci_tr_poly, xlows, xhighs, widths, params['extraction_width_multiplier'], offset=shift_cal1, header=im_arcshead)
            # Begin: This bit is not necessary once the procedure has been written
            im_name = 'master_'+calib[1]+'_long'
            if 'master_'+calib[1]+'_l_filename' in params.keys():
                im_name = params['master_'+calib[1]+'_l_filename'].replace('.fits','')
            im_name = im_name.replace('.fit','')
            save_multispec([cal_l_spec, cal_l_spec, cal_l_spec, cal_l_spec], params['path_extraction']+im_name, im_arclhead)                    # This needs updating!!!
            im_name = 'master_'+calib[1]+'_short'
            if 'master_'+calib[1]+'_s_filename' in params.keys():
                im_name = params['master_'+calib[1]+'_s_filename'].replace('.fits','')
            im_name = im_name.replace('.fit','')
            save_multispec([cal_s_spec, cal_s_spec, cal_s_spec, cal_s_spec], params['path_extraction']+im_name, im_arcshead)
            
            # Identify the Emission lines
            fname = params['logging_found_arc_lines'].replace('.txt','')+calib[0]+'.txt'
            if os.path.isfile(fname) :
                logger('Info: List of the identified emission lines already exists. Using the information from file: {0}'.format(fname ))
                arc_lines_px_txt = read_text_file(fname, no_empty_lines=True)              # list of strings, first entry is header 
                arc_lines_px = np.array( convert_readfile(arc_lines_px_txt[1:], [int, float, float, float], delimiter='\t', replaces=['\n',' ']) )
            else:
                arc_lines_px = identify_lines(params, cal_l_spec, cal_s_spec, good_px_mask_l, good_px_mask_s)
                logger('Info: Identified {0} lines in the arc spectrum. These lines are stored in file {1}'.format(len(arc_lines_px), fname ))
                printarrayformat = ['%1.1i', '%3.2f', '%3.2f', '%3.1f']
                logger('order\tpixel\twidth\theight of the line', show=False, printarrayformat=printarrayformat, printarray=arc_lines_px, logfile=fname)
    
            if calib[1] == 'cal2':          # for the calibration fiber
                if os.path.isfile(params['original_master_wavelensolution_filename']) == False:                                                         # Create a new solution
                    if 'px_to_wavelength_file' not in params:   params['px_to_wavelength_file'] = 'pixel_to_wavelength.txt'
                    wavelength_solution, wavelength_solution_arclines = create_new_wavelength_UI(params, cal_l_spec, cal_s_spec, arc_lines_px, reference_catalog, reference_names)
                    """ original, manual solution
                    if os.path.isfile('pixel_to_wavelength.txt') == False:                                                                             # No pixel_to_wavelength.txt available
                        wavelength_solution, wavelength_solution_arclines = create_pseudo_wavelength_solution(cal_l_spec.shape[0])                      # Create a pseudo solution
                        plot_wavelength_solution_spectrum(cal_l_spec, cal_s_spec, params['logging_arc_line_identification_spectrum'].replace('.pdf','')+'_manual.pdf', 
                                                      wavelength_solution, wavelength_solution_arclines, np.array([0,1,0]).reshape(1,3), ['dummy'], plot_log=True)     # Plot the spectrum
                        logger('Error: Files for creating the wavelength solution do not exist: {0}, {1}. Please check parameter {2} or create {1}.'.format(\
                                                params['original_master_wavelensolution_filename'], 'pixel_to_wavelength.txt', 'original_master_wavelensolution_filename'))
                    wavelength_solution, wavelength_solution_arclines = read_fit_wavelength_solution(params, 'pixel_to_wavelength.txt', cal_l_spec)         # For a new wavelength solution
                    """
                    save_wavelength_solution_to_fits(wavelength_solution, wavelength_solution_arclines, params['original_master_wavelensolution_filename'])                   # For a new wavelength solution
                    plot_wavelength_solution_form(params['logging_wavelength_solution_form'].replace('.png','')+'_manual.png', axlows, axhighs, wavelength_solution)
                    plot_wavelength_solution_spectrum(cal_l_spec, cal_s_spec, params['logging_arc_line_identification_spectrum'].replace('.pdf','')+'_manual.pdf', 
                                                      wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, plot_log=True)
                    params['order_offset'] = [0,0]
                    params['px_offset'] = [-10,10,2]
                    params['px_offset_order'] = [-1,1,1]
                wavelength_solution_ori, wavelength_solution_arclines_ori = read_wavelength_solution_from_fits(params['original_master_wavelensolution_filename'])
                # Find the new wavelength solution
                wavelength_solution, wavelength_solution_arclines = adjust_wavelength_solution(params, np.array(cal_l_spec), arc_lines_px, wavelength_solution_ori, 
                                                                                               wavelength_solution_arclines_ori, reference_catalog, reference_names, xlows, xhighs, params['GUI'])
            else:                           # Science fiber
                params['order_offset'] = [0,0]
                params['px_offset'] = [-60,60,6]
                params['px_offset_order'] = [-1,1,1]
                wavelength_solution, wavelength_solution_arclines = adjust_wavelength_solution(params, np.array(cal_l_spec), arc_lines_px, wavelength_solution, 
                                                                                               wavelength_solution_arclines, reference_catalog, reference_names, xlows, xhighs, params['GUI'])
            save_wavelength_solution_to_fits(wavelength_solution, wavelength_solution_arclines, params['master_wavelensolution'+calib[0]+'_filename'])
            plot_wavelength_solution_form(params['logging_wavelength_solution_form'], axlows, axhighs, wavelength_solution)
            plot_wavelength_solution_spectrum(cal_l_spec, cal_s_spec, params['logging_arc_line_identification_spectrum'], wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, plot_log=True)
            plot_wavelength_solution_image(im_cal_l, params['logging_arc_line_identification_positions'], cal_tr_poly, axlows, axhighs, wavelength_solution, wavelength_solution_arclines, reference_catalog)
        calimages['wave_sol_'+calib[2]] = copy.deepcopy( wavelength_solution )
        calimages['wave_sol_lines_'+calib[2]] = copy.deepcopy( wavelength_solution_arclines )   # Store the information for later
        calimages['arc_l_spec_'+calib[2]] = copy.deepcopy( cal_l_spec )
        calimages['arc_l_head_'+calib[2]] = copy.deepcopy( im_arclhead )
        if params['original_master_wavelensolution_filename'].lower() != 'pseudo': 
            im_head, obsdate_midexp, obsdate_mid_float, jd_midexp = get_obsdate(params, im_arclhead)               # in UTC, mid of the exposure
            params['wavelength_solution_type'] = calib[2]+'-fiber'
            add_text_to_file('{0}\t{1}\t{2}\t{3}'.format(jd_midexp, 0, 0, calib[2]), params['master_wavelengths_shift_filename'], warn_missing_file=False )
            if calib[1] == 'cal1':          # Science fiber
                params['two_solutions'] = True

    # Use the better wavelength solution: should be the one of the science fiber (wavelength_solution_sci), but if calibration fiber is better, then use calibration fiber
    if params['two_solutions']:       # solutions for both fibers
        if np.sum(calimages['wave_sol_lines_cal'] > 100) > 2* np.sum(calimages['wave_sol_lines_sci'] > 100):
            # if twice as many lines were identified (should also check that the residuals are not worse, but this is difficult as solution might have been loaded from file
            params['wavelength_solution_type'] = 'cal-fiber'                # Put the calibration fiber wavelength solution back as main solution
            logger('Info: Using the calibration fiber wavelength solution (first solution) as master wavelength solution, as this solution seems to be better')
        else:
            logger('Info: Using the science fiber wavelength solution (second solution) as master wavelength solution')

    #params['wavelength_solution_type'] = 'cal-fiber'        # just a test
    wtype = params['wavelength_solution_type'][:3]
    wavelength_solution, wavelength_solution_arclines = calimages['wave_sol_'+wtype], calimages['wave_sol_lines_'+wtype]
    # Catch the problem, when the script re-runs with different settings and therefore the number of orders changes.
    if wavelength_solution.shape[0] != sci_tr_poly.shape[0]:
        #print('im_trace1.shape', im_trace1.shape)
        logger('Error: The number of traces for extraction and for the wavelength calibration do not match. Please remove eighter {0} ({2}) or {1} ({3}) and re-run the script in order to solve.'\
                    .format(params['master_trace_sci_filename'], params['master_wavelensolution_filename'], sci_tr_poly.shape[0], wavelength_solution.shape[0]))
    
    # Find the pixel-shift between the two solutions
    if params['two_solutions']:
        params['extract_wavecal'] = True
        if calimages['wave_sol_cal'].shape[0] != calimages['wave_sol_sci'].shape[0]:
            logger('Error: The number of traces for the science ({0}) and calibration ({1}) wavelength solution differ. Please delete the wrong, old file ({2} or {3})'.format(\
                        calimages['wave_sol_sci'].shape[0], calimages['wave_sol_cal'].shape[0], params['master_wavelensolution_sci_filename'], params['master_wavelensolution_filename'] ))
        shift, shift_err = find_shift_between_wavelength_solutions(calimages['wave_sol_cal'], calimages['wave_sol_lines_cal'], calimages['wave_sol_sci'], calimages['wave_sol_lines_sci'], 
                                                                            np.zeros((calimages['wave_sol_cal'].shape[0],im_trace1.shape[0])), ['calibration fiber','science fiber'] )
        add_text_to_file('{0}\t{1}\t{2}\t{3}'.format(0, round(shift,4), round(shift_err,4), 'sci-cal'), params['master_wavelengths_shift_filename'], warn_missing_file=False )                       
        # shift is positive if lines in science are right of lines in calibration
        """if params['wavelength_solution_type'] == 'sci-fiber':
            steps = ['cal-fiber']           # first do it for the not wavelength solution
        else:
            steps = ['sci-fiber']           # first do it for the not wavelength solution
        steps.append(params['wavelength_solution_type'])
        """
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

        dummy_shift_wavesoln, master_shift, im_head = shift_wavelength_solution(params, aspectra, wavelength_solution, wavelength_solution_arclines, reference_catalog, 
                                                              reference_names, xlows, xhighs, obsdate_mid_float, jd_midexp, sci_tr_poly, cal_tr_poly, im_name, maxshift=max(3,2*shift_err), in_shift=shift, im_head=im_head )
        # master shift gives the shift from the wavelength_solution to the aspectra
        params['pxshift_between_wavesolutions'] = master_shift
        # print("0:params['pxshift_between_wavesolutions']", params['pxshift_between_wavesolutions'])
        """ # Comparing the shifted solution with the original wavelength solution -> shifts in either direction, depending on polynomial, maximum shift: -> 0.025 \AA shift at 4880 -> 1.5 km/s
        wsci = create_wavelengths_from_solution(calimages['wave_sol_sci'], calimages['arc_l_spec_cal'])
        wcal = create_wavelengths_from_solution(calimages['wave_sol_cal'], calimages['arc_l_spec_cal'])
        wshift = create_wavelengths_from_solution(dummy_shift_wavesoln, calimages['arc_l_spec_cal'])
        save_multispec([wsci, calimages['arc_l_spec_sci']], 'fib1_wsci.fits', im_arclhead, bitpix=params['extracted_bitpix'])
        save_multispec([wsci, calimages['arc_l_spec_cal']], 'fib2_wsci.fits', im_arclhead, bitpix=params['extracted_bitpix'])
        save_multispec([wcal, calimages['arc_l_spec_sci']], 'fib1_wcal.fits', im_arclhead, bitpix=params['extracted_bitpix'])
        save_multispec([wcal, calimages['arc_l_spec_cal']], 'fib2_wcal.fits', im_arclhead, bitpix=params['extracted_bitpix'])
        save_multispec([wshift, calimages['arc_l_spec_sci']], 'fib1_wshift.fits', im_arclhead, bitpix=params['extracted_bitpix'])
        save_multispec([wshift, calimages['arc_l_spec_cal']], 'fib2_wshift.fits', im_arclhead, bitpix=params['extracted_bitpix'])
        #"""
    else:                           # read the shift to the science fiber
        all_shifts = read_text_file(params['master_wavelengths_shift_filename'], no_empty_lines=True)
        all_shifts = convert_readfile(all_shifts, [float, float, float, str], delimiter='\t', replaces=['\n']+[['  ',' ']]*20)       # jd_midexp, shift_avg, shift_std, can contain duplicate jd_midexp (last one is the reliable one)
        for entry in all_shifts[::-1]:
            if entry[3] == 'sci-cal': 
                params['pxshift_between_wavesolutions'] = - entry[1]
                break
    
    params['extract_wavecal'] = False
    
    im_blazecor, im_blazecor_head = create_image_general(params, 'blazecor')    # -> cont1cal2
    
    # Extract the flat spectrum and normalise it
    if os.path.isfile(params['master_blaze_spec_norm_filename']) :
        logger('Info: Normalised flat already exists: {0}'.format(params['master_blaze_spec_norm_filename']))
        # The file is read later on purpose
    else:
        logger('Step: Create the normalised flat for the night')
        im_head, obsdate_midexp, obsdate_mid_float, jd_midexp = get_obsdate(params, im_blazecor_head)
        shift, im_head = find_shift_images(params, im_blazecor, im_trace1, sci_tr_poly, xlows, xhighs, widths, 0, cal_tr_poly, extract=True, im_head=im_head)
        flat_spec, good_px_mask, extr_width = extract_orders(params, im_blazecor, sci_tr_poly, xlows, xhighs, widths, params['extraction_width_multiplier'], var='linfit', offset=shift)
        flat_spec_norm = flat_spec/np.nanmedian(flat_spec)
        flat_spec_norm_cor = correct_blaze(flat_spec_norm, minflux=0.001)         # Ignore all areas where the flux is 0.5% of median flux
        if np.nansum(np.abs(sci_tr_poly - cal_tr_poly)) == 0.0:                                                 # science and calibration traces are at the same position
            blazecor_spec, agood_px_mask = flat_spec*0, copy.copy(good_px_mask)
        else:
            blazecor_spec, agood_px_mask, extr_width = extract_orders(params, im_blazecor, cal_tr_poly, axlows, axhighs, awidths, params['arcextraction_width_multiplier'], offset=shift)
        wavelength_solution_shift, shift, im_head = shift_wavelength_solution(params, blazecor_spec, wavelength_solution, wavelength_solution_arclines, 
                                            reference_catalog, reference_names, xlows, xhighs, obsdate_mid_float, jd_midexp, sci_tr_poly, cal_tr_poly, params['master_blaze_spec_norm_filename'], im_head=im_head)
        wavelengths = create_wavelengths_from_solution(wavelength_solution_shift, blazecor_spec)
        save_multispec([wavelengths, flat_spec_norm, flat_spec_norm_cor, blazecor_spec], params['master_blaze_spec_norm_filename'], im_blazecor_head)
        #save_im_fits(params, flat_spec_norm, im_sflat_head, params['master_flat_spec_norm_filename'])
    
    logger('Info: Finished routines for a new night of data. Now science data can be extracted. Please check before the output in the loging directory {0}: Are all orders identified correctly for science and calibration fiber, are the correct emission lines identified for the wavelength solution?\n'.format(params['logging_path']))
    
    obj_names = []
    extractions = []
    wavelengthcals_cal, wavelengthcals_sci = [], []
    for entry in params.keys():
        if entry.find('extract') >= 0 and entry.find('_rawfiles') >= 0:
            extractions.append(entry.replace('_rawfiles',''))
        if entry.find('wavelengthcal2') >= 0 and entry.find('_rawfiles') >= 0:          # ThAr in science fiber
            wavelengthcals_sci.append(entry.replace('_rawfiles',''))
        elif entry.find('wavelengthcal') >= 0 and entry.find('_rawfiles') >= 0:         # ThAr in calibration fiber
            wavelengthcals_cal.append(entry.replace('_rawfiles',''))
            
    flat_spec_norm = np.array(fits.getdata(params['master_blaze_spec_norm_filename']))              # read it again, as the file is different than the data above
    # Catch the problem, when the script re-runs with different settings and therefore the number of orders changes.
    if flat_spec_norm.shape[1] != wavelength_solution.shape[0]:
        #print('im_trace1.shape', im_trace1.shape)
        logger('Error: The number of traces in the blaze of the blaze correction and for the wavelength calibration do not match. Please remove {0} ({1} instead of expected {2}) and re-run the script in order to solve.'\
                    .format(params['master_blaze_spec_norm_filename'], flat_spec_norm.shape[1], wavelength_solution.shape[0]))
    
    remove_orders, keep_orders = remove_orders_low_flux(params, flat_spec_norm)
    if len(remove_orders) > 0:
        #print wavelength_solution.shape, wavelength_solution_arclines.shape, sci_tr_poly.shape, xlows.shape, xhighs.shape, widths.shape, cal_tr_poly.shape, axlows.shape, axhighs.shape, awidths.shape, flat_spec_norm.shape
        wavelength_solution, wavelength_solution_arclines, sci_tr_poly, xlows, xhighs, widths, cal_tr_poly, axlows, axhighs, awidths, flat_spec_norm = \
                    wavelength_solution[keep_orders,:], wavelength_solution_arclines[keep_orders,:], \
                    sci_tr_poly[keep_orders,:,:], xlows[keep_orders], xhighs[keep_orders], widths[keep_orders,:], \
                    cal_tr_poly[keep_orders,:,:], axlows[keep_orders], axhighs[keep_orders], awidths[keep_orders,:], \
                    flat_spec_norm[:,keep_orders,:]                                                 # remove the bad orders
    
    calimages['flat_spec_norm'] = copy.deepcopy( flat_spec_norm )
    calimages['sci_trace'] = copy.deepcopy( [sci_tr_poly, xlows, xhighs, widths] )
    calimages['cal_trace'] = copy.deepcopy( [cal_tr_poly, axlows, axhighs, awidths] )
    
    if ( params['arcshift_side'] == 0 or params['two_solutions'] ) and len(wavelengthcals_cal)+len(wavelengthcals_sci) > 0:         # no calibration spectrum at the same time
        def wavecal_multicore(parameter):
                    wavelengthcal, fib, im_name_full = parameter
                    # !!! Posible improvement: combine a few files if they are taken close to each other
                    im_name = im_name_full.rsplit(os.sep)
                    im_name = im_name[-1].rsplit('.',1)         # remove the file ending
                    im_name = im_name[0]
                    im_name_wc = im_name+'_wave'+fib
                    if os.path.isfile(params['path_extraction']+im_name_wc+'.fits'):
                        logger('Info: File {0} was already processed for the calibration of the wavelength solution. If you want to extract again, please delete {1}{0}.fits'.format(im_name_wc, params['path_extraction']))
                        return
                    params['calibs'] = params[wavelengthcal+'_calibs_create']
                    im, im_head = read_file_calibration(params, im_name_full)
                    extraction_wavelengthcal(params, im, im_name_wc, im_head, sci_tr_poly, xlows, xhighs, widths, cal_tr_poly, axlows, axhighs, awidths, \
                                                    wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, flat_spec_norm, im_trace1, im_name)
        
        
        logger('Info: Starting to extract wavelength calibrations')
        if params['use_cores'] > 1:
            print('Note: Will use multiple cores, hence output will be for several files in parallel')
        params['extract_wavecal'] = True                                                                                # necessary for shift_wavelength_solution so the shift is stored in a file
        all_wavelengthcals = []
        for [wavelengthcals, fib] in [ [wavelengthcals_cal,'cal'], [wavelengthcals_sci,'sci'] ]:
            for wavelengthcal in wavelengthcals:
                for im_name_full in params[wavelengthcal+'_rawfiles']:
                    all_wavelengthcals.append([ wavelengthcal, fib, im_name_full ])
        if params['use_cores'] > 1:
            logger('Info using multiprocessing on {0} cores'.format(params['use_cores']))
            p = multiprocessing.Pool(params['use_cores'])
            p.map(wavecal_multicore, all_wavelengthcals)
            p.terminate()
        else:
            for all_wavelengthcal in all_wavelengthcals:
                wavecal_multicore(all_wavelengthcal)
        """for [wavelengthcals, fib] in [ [wavelengthcals_cal,'cal'], [wavelengthcals_sci,'sci'] ]:
            for wavelengthcal in wavelengthcals:
                for im_name_full in params[wavelengthcal+'_rawfiles']:
                    # !!! Posible improvement: combine a few files if they are taken close to each other
                    im_name = im_name_full.rsplit(os.sep)
                    im_name = im_name[-1].rsplit('.',1)         # remove the file ending
                    im_name = im_name[0]
                    im_name_wc = im_name+'_wave'+fib
                    if os.path.isfile(params['path_extraction']+im_name_wc+'.fits'):
                        logger('Info: File {0} was already processed for the calibration of the wavelength solution. If you want to extract again, please delete {1}{0}.fits'.format(im_name_wc, params['path_extraction']))
                        continue
                    params['calibs'] = params[wavelengthcal+'_calibs_create']
                    im, im_head = read_file_calibration(params, im_name_full)
                    extraction_wavelengthcal(params, im, im_name_wc, im_head, sci_tr_poly, xlows, xhighs, widths, cal_tr_poly, axlows, axhighs, awidths, \
                                                    wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, flat_spec_norm, im_trace1, im_name)"""
    params['extract_wavecal'] = False
    if len(extractions) == 0:                                               # no extractions to do
        logger('Warn: Nothing to extract. -> Exiting')
        header_results_to_texfile(params)           # Save the results from the header in a logfile
        exit(0)
    logger('Info: Starting to extract spectra')
    if params['use_cores'] > 1:
        print('Note: Will use multiple cores, hence output will be for several files in parallel')
    def extraction_multicore(all_extractions):
        [extraction, im_name_full] = all_extractions
        if  extraction.find('extract_combine') == -1:     # Single file extraction
            #for im_name_full in params[extraction+'_rawfiles']:
            if True:
                im_name = im_name_full.rsplit(os.sep)
                im_name = im_name[-1].rsplit('.',1)         # remove the file ending
                im_name = im_name[0]
                if os.path.isfile(params['path_extraction']+im_name+'.fits'):
                    logger('Info: File {0} was already processed. If you want to extract again, please delete {1}{0}.fits'.format(im_name, params['path_extraction']))
                    return ''
                #print extraction, im_name_full, im_name
                params['calibs'] = params[extraction+'_calibs_create']
                im, im_head = read_file_calibration(params, im_name_full)
                
        else:                                       # Combine files before extraction
            im_name = extraction
            if os.path.isfile(params['path_extraction']+im_name+'.fits'):
                logger('Info: File {0} was already processed. If you want to extract again, please delete {1}{0}.fits'.format(im_name, params['path_extraction']))
                return ''
            im, im_head = create_image_general(params, extraction)
        obj_name = extraction_steps(params, im, im_name, im_head, sci_tr_poly, xlows, xhighs, widths, cal_tr_poly, axlows, axhighs, awidths, 
                                    wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, flat_spec_norm, im_trace1)
        return obj_name.lower()
    
    all_extractions = []
    for extraction in extractions:
        if  extraction.find('extract_combine') == -1:     # Single file extraction
            for im_name_full in params[extraction+'_rawfiles']:
                all_extractions.append([extraction, im_name_full])
        else:
            all_extractions.append([extraction, extraction])
    if params['use_cores'] > 1:
        logger('Info using multiprocessing on {0} cores'.format(params['use_cores']))
        p = multiprocessing.Pool(params['use_cores'])
        obj_names = p.map(extraction_multicore, all_extractions)
        p.terminate()
    else:
        for all_extraction in all_extractions:
            obj_names.append( extraction_multicore(all_extraction) )
            
    obj_names = np.unique(obj_names)
    obj_names = list(obj_names[obj_names != ''])
    
    """for extraction in extractions:
        if  extraction.find('extract_combine') == -1:     # Single file extraction
            for im_name_full in params[extraction+'_rawfiles']:
                im_name = im_name_full.rsplit(os.sep)
                im_name = im_name[-1].rsplit('.',1)         # remove the file ending
                im_name = im_name[0]
                if os.path.isfile(params['path_extraction']+im_name+'.fits'):
                    logger('Info: File {0} was already processed. If you want to extract again, please delete {1}{0}.fits'.format(im_name, params['path_extraction']))
                    continue
                #print extraction, im_name_full, im_name
                params['calibs'] = params[extraction+'_calibs_create']
                im, im_head = read_file_calibration(params, im_name_full)
                obj_name = extraction_steps(params, im, im_name, im_head, sci_tr_poly, xlows, xhighs, widths, cal_tr_poly, axlows, axhighs, awidths, \
                                                    wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, flat_spec_norm, im_trace1)
                if obj_name not in obj_names:
                    obj_names.append(obj_name)
        else:                                       # Combine files before extraction
            im_comb, im_comb_head = create_image_general(params, extraction)
            im_name = extraction
            if os.path.isfile(params['path_extraction']+im_name+'.fits'):
                logger('Info: File {0} was already processed. If you want to extract again, please delete {1}{0}.fits'.format(im_name, params['path_extraction']))
                continue
            obj_name = extraction_steps(params, im_comb, im_name, im_comb_head, sci_tr_poly, xlows, xhighs, widths, cal_tr_poly, axlows, axhighs, awidths, \
                                        wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, flat_spec_norm, im_trace1)
                                        # puts the files in obj_name.lower()
            if obj_name.lower() not in obj_names:
                obj_names.append(obj_name.lower())"""
    logger('')      # To have an empty line
    logger('Info: Finished extraction of the science frames. The extracted {0}*.fits file contains different data in a 3d array in the form: data type, order, and pixel. First data type is the wavelength (barycentric corrected), second is the extracted spectrum, followed by a measure of error. Forth and fith are the flat corrected spectra and its error. Sixth and sevens are the the continium normalised spectrum and the S/N in the continuum. Eight is the bad pixel mask, marking data, which is saturated or from bad pixel. The nineth entry is the spectrum of the calibration fiber. The last entry is the wavelength without barycentric correction'.format(params['path_extraction']))
    logger('Info: Will try to do the RV analysis in a moment') 
    header_results_to_texfile(params)           # Save the results from the header in a logfile
    #time.sleep(2)
    
    # Prepare files for TERRA and SERVAL and CERES
    if np.max(wavelength_solution[:,-1]) > 100:
     run_RV = True
     files_RV = []
     headers = dict()
     for file_RV in sorted(os.listdir(params['path_extraction'])):
        if not file_RV.endswith(".fits"):
            continue
        do_RV = True
        for no_RV_name in params['no_RV_names']:
            if file_RV.lower().find(no_RV_name) in [0,1,2,3,4,5]:
                do_RV = False
                break
        if not do_RV:
            continue
        im = fits.open(params['path_extraction']+file_RV)
        im_head = im[0].header
        headers[file_RV] = im_head      # for CERES
        files_RV.append(file_RV)        # for CERES
        spec = im[0].data
        if 'HiFLEx OBJNAME' not in im_head.keys():
            continue
        obj_name = im_head['HiFLEx OBJNAME'].lower()
        obsdate_midexp = datetime.datetime.strptime(im_head['HiFLEx DATE-MID'],"%Y-%m-%dT%H:%M:%S.%f")
        # CSV file for TERRA
        fname = params['path_rv_terra']+obj_name+'/data/'+obsdate_midexp.strftime('%Y-%m-%d%H%M%S')
        save_spec_csv(spec[params['dataset_rv_analysis'][0]], spec[0], spec[7], fname)        # spec[1]: Flux, spec[5]: Continuum corrected
    
        # Store in a text file for serval
        numbers_levels = params['path_rv_serval'].count(os.sep, 2)  # start at 2 to not count './'
        add_text_to_file('../'*numbers_levels+params['path_extraction']+file_RV, 
                         params['path_rv_serval']+'filelist_{0}.txt'.format(obj_name), warn_missing_file=False )
        
     # Do the TERRA RVs
     if os.path.isfile(params['terra_jar_file'])  and os.path.exists(params['path_rv_terra']):
     
        def run_terra(obj_name):
            do_RV = True
            for no_RV_name in params['no_RV_names']:
                if obj_name.lower().find(no_RV_name) in [0,1,2,3,4,5]:
                    do_RV = False
                    break
            if not do_RV:
                return
            os.system('rm -f {0}{1}/results/synthetic.rv'.format(params['path_rv_terra'],obj_name) )     # Delete the old solution, as won't be created otherwise
            newfile = 'echo "0998     synthetic         LAB                LAB                    0.0          0.0       0.0       {0}/" > astrocatalog{0}.example'.format(obj_name)
            logger('For TERRA: creating a new astrocatalog.example with: '+newfile)
            os.system('rm -f astrocatalog{0}.example; '+newfile)
            cmd = 'java -jar {1} -ASTROCATALOG astrocatalog{2}.example 998 -INSTRUMENT CSV {0}'.format(wavelength_solution.shape[0],params['terra_jar_file'], obj_name )
            log = 'logTERRA_{0}'.format(obj_name)
            logger('For TERRA: running TERRA: '+cmd)
            logger('Info: TERRA output and errors can be watched in {0}'.format(log))
            if run_RV:
                with open(log, 'a') as logf:
                    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, shell=True, stdout=logf, stderr=subprocess.STDOUT)            # This doesn't wait until the procress has been finished
                    p.communicate()                                                     # This makes it wait until it finished running
            else:
                logger('Warn: TERRA commented out')
            resultfile = '{2} {0}/results/synthetic.rv'.format(obj_name, params['path_rv_terra'], params['editor'])
            logger('For TERRA: results can be opened with: '+resultfile)
            
        logger('Info: Preparing for the TERRA analysis.')
        for root, dirs, files in os.walk(params['path_rv_terra'], followlinks=True):                       # Find all the objects again, as won't be added to obj_names when re-run
            for file in files:
                if file.endswith('.csv'):                       # has the file the correct ending?
                    filename = os.path.join(root, file).replace(params['path_rv_terra'],'')                # Only relative folder and filename
                    obj_name = filename.split(os.sep)[0]
                    if obj_name not in obj_names:
                        obj_names.append(obj_name)
        logger('For TERRA: changing directory to '+params['path_rv_terra']+' . The steps to run TERRA are given in the logfile in that folder.')
        os.chdir(params['path_rv_terra'])
        logger('Info: All data logged in this file is relative to '+params['path_rv_terra'])
        if params['use_cores'] > 1:
            logger('Info using multiprocessing on {0} cores'.format(params['use_cores']))
            p = multiprocessing.Pool(params['use_cores'])
            p.map(run_terra, obj_names)
            p.terminate()
        else:
            for obj_name in obj_names:
                run_terra(obj_name)
        os.chdir(params['path_run'])        # Go back to previous folder
        print('')
        logger('Info: Some errors reported by TERRA are expected (reading DRS ephemeris). The results are stored in {0}<object name>/results/synthetic.rv'.format(params['path_rv_terra']))
     else:
        logger('TERRA is not installed or the path to the terra_jar_file wrong (currently: {0}), or the path for TERRA csv files does not exist ({1})'.format(params['terra_jar_file'], params['path_rv_terra']))
    
     # Do the SERVAL RVs
     if os.path.exists(params['path_serval']+'serval') and os.path.exists(params['path_serval']+'python') and os.path.exists(params['path_rv_serval']):
    
        def run_serval(obj_name):           # create a procedure to run on multicore
            do_RV = True
            for no_RV_name in params['no_RV_names']:
                if obj_name.lower().find(no_RV_name) in [0,1,2,3,4,5]:
                    do_RV = False
                    break
            if not do_RV or not os.path.isfile('filelist_{0}.txt'.format(obj_name)):
                return
                #continue
            cmd = 'ulimit -n 4096 ; {4}serval/src/serval.py {3} filelist_{3}.txt -inst HIFLEX -targrv 0 -pmin {0} -pmax {1} -oset {2} -safemode 2'.format(params['pmin'], 
                                            params['pmax'], params['oset'], obj_name, params['path_serval'])
            log = 'logSERVAL_{0}'.format(obj_name)
            logger('For SERVAL: running SERVAL: '+cmd)
            logger('Info: SERVAL output and errors can be watched in {0}'.format(log))
            if run_RV:
                with open(log, 'a') as logf:
                    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, shell=True, stdout=logf, stderr=subprocess.STDOUT)    # shell=True is necessary
                    p.communicate(input="\n")                                       # This waits until the process needs the enter
            else:
                logger('Warn: SERVAL commented out')
            resultfile = '{2} {0}/{0}.rvc.dat'.format(obj_name, params['path_rv_serval'], params['editor'])  
            logger('For SERVAL: results can be opened with: '+resultfile)
        
        logger('Info: Preparing for the SERVAL analysis.')
        for root, dirs, files in os.walk(params['path_rv_serval'], followlinks=True):                       # Find all the objects again, as won't be added to obj_names when re-run
            for file in files:
                if file.endswith('.txt') and file.find('filelist_') == 0:                       # has the filename the correct format?
                    obj_name = file.split('filelist_')[-1].split('.txt')[0]
                    if obj_name not in obj_names:
                        obj_names.append(obj_name)                   
        hiflex_file = params['path_rv_serval']+'conf_hiflex_serval.txt'
        with open(hiflex_file, 'w') as file:
            file.write('# This file is used to control which data is used in SERVAL. It is read by inst_HIFLEX\n')
            file.write('orders = {0}\ndata_dataset = {1}\nerror_dataset = {2}\nwave_dataset = {3}\nmask_dataset = {4}\n'.format(
                                wavelength_solution.shape[0], params['dataset_rv_analysis'][0], params['dataset_rv_analysis'][1], 9, 7))

        pypath = ''
        if 'PYTHONPATH' in os.environ.keys():
            pypath = os.environ["PYTHONPATH"] + os.pathsep
        pypath += params['path_serval']+'python'
        logger('For SERVAL: set variable: bash: export PYTHONPATH={0} ; csh: setenv PYTHONPATH {0}'.format(pypath))
        os.environ["PYTHONPATH"] = pypath
        xran = np.max(xhighs) - np.min(xlows)
        exclude_x = 0.1         # 0.1: 10% of outermost pixel on each side are excluded
        params['pmin'] = int(np.min(xlows) + 0.1*xran)
        params['pmax'] = int(np.max(xhighs) - 0.1*xran)
        params['oset'] = '{0}:{1}'.format(0,wavelength_solution.shape[0])     # if shape is 9, then oset will lead to 0,1,2,3,4,5,6,7,8 being used
        logger('For SERVAL: changing directory to '+params['path_rv_serval']+' . The steps to run SERVAL are given in the logfile in that folder.')
        os.chdir(params['path_rv_serval'])
        logger('Info: All data logged in this file is relative to '+params['path_rv_serval'])
        #with multiprocessing.Pool(params['use_cores']) as p:       # only possible in python3
        #    p.map(run_serval, obj_names)
        if params['use_cores'] > 1:
            logger('Info using multiprocessing on {0} cores'.format(params['use_cores']))
            p = multiprocessing.Pool(params['use_cores'])
            p.map(run_serval, obj_names)
            p.terminate()
        else:
            for obj_name in obj_names:
                run_serval(obj_name)
        os.chdir(params['path_run'])        # Go back to previous folder
        print('')
        logger(('Info: Finished the SERVAL analysis. Some errors reported by SERVAL are expected.'+\
              ' The results are stored in {0}<object name>/<object name>.rvc.dat.'+\
              ' If serval failed (the result file is missing), run it again using less orders by setting oset to a smaller range (especially orders with low SN).'+\
              ' The command history can be found in {0}cmdhistory.txt. Before running serval: cd {0}').format(params['path_rv_serval']))
     else:
        logger('SERVAL is not installed or path_serval is wrong (currently: {0}), or the folder for the "filelist_* are missing in {1}'.format(params['path_serval'], params['path_rv_serval']))
     
     # Do the CERES RVs
     force_stellar_pars = False
     force_rvs = False
     if os.path.exists(params['path_ceres']+'utils/Correlation') and os.path.exists(params['path_ceres']+'utils/GLOBALutils') \
            and os.path.exists(params['path_ceres']+'utils/OptExtract') and os.path.exists(params['path_ceres']+'utils/CCF'):         # if not pseudo-solution and necessary files exist
        
        base = params['path_ceres']
        models_path = base + 'data/COELHO_MODELS/R_40000b/'
        load_ceres_modules(params)
        
        stellar_par = dict()
        for file_RV in files_RV:
            im_head = headers[file_RV]
            if 'HiFLEx OBJNAME' not in im_head.keys():
                continue
            T_eff =  im_head.get('HiFLEx CERES Teff', np.nan)
            logg  =  im_head.get('HiFLEx CERES logg', np.nan)
            Z     =  im_head.get('HiFLEx CERES Z', np.nan)
            vsini =  im_head.get('HiFLEx CERES vsini', np.nan)
            vel0  =  im_head.get('HiFLEx CERES vel0', np.nan)
            hardcoded =  np.isnan([T_eff, logg, Z, vsini, vel0]).any()      # If one is np.nan -> measure again
            obname = im_head['HiFLEx OBJNAME']
            if (hardcoded or force_stellar_pars) and os.path.exists(models_path):       # Don't use the header keywords
                im = fits.open(params['path_extraction']+file_RV)
                spec = np.array(im[0].data, dtype=np.float64)
                spec[0,:,:] = wavelength_air_to_vacuum(spec[9,:,:])     # Vaccuum wavelength
                spec = clip_noise(spec)
                T_eff, logg, Z, vsini, vel0, hardcoded = stellar_params_ceres(params, spec, file_RV.replace('.fits',''), wavelength_solution)
                if not hardcoded:
                    im_head['HIERARCH HiFLEx CERES Teff']  = (T_eff, 'Teff measured by CERES CCF')
                    im_head['HIERARCH HiFLEx CERES logg']  = (logg,  'logg measured by CERES CCF')
                    im_head['HIERARCH HiFLEx CERES Z']     = (Z,     'Z measured by CERES CCF')
                    im_head['HIERARCH HiFLEx CERES vsini'] = (vsini, 'vsini measured by CERES CCF')
                    im_head['HIERARCH HiFLEx CERES vel0']  = (vel0,  'vel0 measured by CERES CCF')
                    headers[file_RV] = im_head
            if not hardcoded:
                    if obname in stellar_par.keys():
                        stellar_par[obname] = np.vstack(( stellar_par[obname], [T_eff, logg, Z, vsini, vel0] ))
                    else:
                        stellar_par[obname] = np.array([T_eff, logg, Z, vsini, vel0])
                        stellar_par[obname].shape = (1,5)       # Otherwise the median later will be a problem

        for obname in stellar_par.keys():
            nspec = stellar_par[obname].shape[0]
            ori = copy.deepcopy(stellar_par[obname])
            stellar_par[obname] = np.round(np.median(stellar_par[obname], axis=0), 2)
            if nspec >= 2:
                stdev = np.round(np.std(ori, axis=0, ddof=1), 2)
                text = 'Teff = {0} +- {5} , logg = {1} +- {6} , Z = {2} +- {7} , vsini = {3} +- {8} , vel0 = {4} +- {9}'.format(stellar_par[obname][0], 
                         stellar_par[obname][1], stellar_par[obname][2], stellar_par[obname][3], stellar_par[obname][4], stdev[0], stdev[1], stdev[2], stdev[3], stdev[4])
            else:
                text = 'Teff = {0} , logg = {1} , Z = {2} , vsini = {3} , vel0 = {4}'.format(stellar_par[obname][0],
                         stellar_par[obname][1], stellar_par[obname][2], stellar_par[obname][3], stellar_par[obname][4])
            logger('Info: With the CERES routine the following parameters have been determined for {0} using {1} spectra: {2}'.format(obname, nspec, text))
        
        def prepare_ceres_rv(file_RV):
            
            #params = kwargs['params']
            #stellar_par = kwargs['stellar_par']
            #spec = kwargs['spec']
            # = kwargs['']
            # = kwargs['']
            # = kwargs['']
            if file_RV not in headers.keys():
                return
            im_head = headers[file_RV]
            if ( 'HiFLEx CERES RV_BARY' in im_head.keys() or 'HiFLEx CERES RV' in im_head.keys() ) and not force_rvs:
                return
            if 'HiFLEx OBJNAME' not in im_head.keys():
                return
            obname = im_head['HiFLEx OBJNAME']
            im = fits.open(params['path_extraction']+file_RV)
            spec = np.array(im[0].data, dtype=np.float64)
            spec[0,:,:] = wavelength_air_to_vacuum(spec[9,:,:])     # Vaccuum wavelength from the wavelengths not corrected by barycentric velocity
            spec = clip_noise(spec)
            
            if obname not in stellar_par.keys():            # hardcoded values
                T_eff, logg, Z, vsini, vel0 = 5777, 4.4374, 0.0134, 2, 0
            else:
                [T_eff, logg, Z, vsini, vel0] = list(stellar_par[obname])
            gobs = ephem.Observer()  
            gobs.name = copy.copy(params['site'])
            gobs.lat  = math.radians(params['latitude'])     # lat/long in decimal degrees  
            gobs.long = math.radians(params['longitude'])
            gobs.elevation = params['altitude']
            gobs.date = im_head['HiFLEx DATE-MID'].replace('T',' ')  # to make into ('%Y-%m-%d %H:%M:%S')
            mephem    = ephem.Moon()
            mephem.compute(gobs)
            
            # Change the path to the object_file to the result_path, if necessary
            object_file = params['object_file']         # needed later
            object_files, object_files_full = find_file_in_allfolders(object_file, [params['result_path']] + params['raw_data_paths'])      # Check also the result and raw data paths for object names
            ra2, dec2, pmra, pmdec, epoch = 0., 0., 0., 0., 2000.
            for object_file in object_files:        # Will be run in any case, even if empty
                ra2, dec2, epoch, pmra, pmdec, dummy, dummy = getcoords_from_file([obname], 0, filen=object_file, warn_notfound=False)        # mjd=0 because because not using CERES to calculated BCV
                if ra2 !=0 or dec2 != 0:                                           # Found the objects -> obnames becomes the single entry which matched entry in params['object_file']
                    break
            #if ra2 ==0 and dec2 == 0:
            #    'Was not found, but no problem'
            RV, RVerr2, BS, BSerr = rv_analysis_ceres(params, spec, file_RV.replace('.fits',''), obname, object_file, mephem, vsini)
            bjdtdb = im_head.get('HiFLEx BJDTDB', np.nan)
            bcvel_baryc = im_head.get('HiFLEx BCV', np.nan)
            if np.isnan(bcvel_baryc):
                RV_BARY = round(RV,4)
                im_head['HIERARCH HiFLEx CERES RV'] = (RV_BARY, 'RV in km/s (without BCV)')
            else:
                RV_BARY = round(RV+bcvel_baryc,4)
                im_head['HIERARCH HiFLEx CERES RV'] = (round(RV,4), 'RV in km/s (without BCV)')
                im_head['HIERARCH HiFLEx CERES RV_BARY'] = (RV_BARY, 'barycentric RV in km/s (measured RV+BCV)')
            # Air to Vacuum wavelength difference only causes < 5 m/s variation: https://www.as.utexas.edu/~hebe/apogee/docs/air_vacuum.pdf (no, causes 83.15 km/s shift, < 5m/s is between the different models)
            logger(('Info: The radial velocity (including barycentric correction) for {0} at [ {6} , JD= {8} , BJDTDB= {9} (center)] gives: '+\
                    'RV = {1} +- {2} km/s, Barycentric velocity = {5} km/s, and BS = {3} +- {4} km/s. '+\
                    'The total extracted flux is {7} counts.').format(file_RV, RV_BARY, round(RVerr2,4), round(BS,4), 
                            round(BSerr,4), bcvel_baryc, im_head['HiFLEx DATE-MID'], np.nansum(spec[1,:,:]), 
                            im_head['HiFLEx JD'], bjdtdb ))
            im_head['HIERARCH HiFLEx CERES RV_ERR']  = (round(RVerr2,4), 'Uncertainty RV in km/s')
            im_head['HIERARCH HiFLEx CERES BS']      = (round(BS,4), 'Bisector')
            im_head['HIERARCH HiFLEx CERES BS_ERR']  = (round(BSerr,4), 'Uncertainty BS')
        
            save_multispec(im[0].data, params['path_extraction']+file_RV, im_head)
        
        if params['use_cores'] > 1:
            logger('Info using multiprocessing on {0} cores'.format(params['use_cores']))
            p = multiprocessing.Pool(params['use_cores'])
            p.map(prepare_ceres_rv, files_RV)
            p.terminate()
        else:
            for file_RV in files_RV:
                prepare_ceres_rv(file_RV)
     else:
        logger('Warn: CERES RV analysis did not run. If this a mistake, please check that it is installed under {0} and includes the folders: utils/Correlation, utils/GLOBALutils, utils/OptExtract, utils/CCF'.format(params['path_ceres']))
            
     rv_results_to_hiflex(params)        # also as new script
    else:
     logger('Info: Using a pseudo wavelength solution -> no RV analysis')
    
    header_results_to_texfile(params)           # Save the results from the header in a logfile
    log_params(params)
    

        
    

