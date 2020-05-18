#!/usr/bin/env python

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
if __name__ == "__main__":
    logger('Info: Starting routines for a new night of data, including: finding or shifting orders, find calibration orders, create wavelength solution, and create normalised blaze function')
    log_params(params)
    params['path_run'] = os. getcwd()
    params['extract_wavecal'] = False
    
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
    if os.path.isfile(params['master_trace_sci_filename']) == True:
        logger('Info: Using exiting trace solution: {0}'.format(params['master_trace_sci_filename']))
        sci_tr_poly, xlows, xhighs, widths = read_fits_width(params['master_trace_sci_filename'])
    else:
        printresults = False
        # load the original solution
        if os.path.isfile(params['original_master_traces_filename']) == True:
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
            if params['update_width_orders'] == True:
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
            bisector_measurements_orders(im_trace1,  params['logging_traces_bisector'], sci_tr_poly, xlows, xhighs, widths)
        
    """Not really useful
    # Create the background map, if it doesn't exist
    if os.path.isfile(params['background_filename']) == True:
        logger('Info: Background map already exists: {0}'.format(params['result_path']+params['background_filename']))
    else:
        # create the background map
        if params['GUI'] == True:
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
    if os.path.isfile(params['master_trace_cal_filename']) == True:
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
    
    update_calibration_memory('sci_trace',[sci_tr_poly, xlows, xhighs, widths])         # Apertures might be shifted before extraction -> this would also affect the localbackground
    update_calibration_memory('cal_trace',[cal_tr_poly, axlows, axhighs, awidths])
    
    # Do the wavelength solution stuff: Find for calibration fiber (and for science fiber if necessary, in this case some comparison between the solutions is necessay
    params['two_solutions'] = False
    for calib in [ ['','cal2'], ['_sci','cal1'] ]:
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
            if calib[1] == 'cal2':
                shift_cal2, im_arclhead = find_shift_images(params, im_cal_l, im_trace1, sci_tr_poly, xlows, xhighs, widths, 0, cal_tr_poly, extract=True, im_head=im_arclhead)
                cal_l_spec, good_px_mask_l, extr_width = extract_orders(params, im_cal_l, cal_tr_poly, axlows, axhighs, awidths, params['arcextraction_width_multiplier'], offset=shift_cal2, header=im_arclhead)
            else:
                shift_cal1, im_arclhead = find_shift_images(params, im_cal_l, im_trace1, sci_tr_poly, xlows, xhighs, widths, 0, cal_tr_poly, extract=True, im_head=im_arclhead)
                cal_l_spec, good_px_mask_l, extr_width = extract_orders(params, im_cal_l, sci_tr_poly, xlows, xhighs, widths, params['extraction_width_multiplier'], offset=shift_cal1, header=im_arclhead)
        if os.path.isfile(params['master_wavelensolution'+calib[0]+'_filename']) == True:
            logger('Info: wavelength solution already exists: {0}'.format(params['master_wavelensolution'+calib[0]+'_filename']))
            wavelength_solution, wavelength_solution_arclines = read_wavelength_solution_from_fits(params['master_wavelensolution'+calib[0]+'_filename'])
        elif params['original_master_wavelensolution_filename'].lower() == 'pseudo':
            logger('Warning: Using a pseudo solution for the wavelength (1 step per px)')
            wavelength_solution, wavelength_solution_arclines = create_pseudo_wavelength_solution(sci_tr_poly.shape[0])
        else:
            im_cal_s, im_arcshead = create_image_general(params, calib[1]+'_s')
            if calib[1] == 'cal2':
                cal_s_spec, good_px_mask_s, extr_width = extract_orders(params, im_cal_s, cal_tr_poly, axlows, axhighs, awidths, params['arcextraction_width_multiplier'], offset=shift_cal2, header=im_arcshead)
            else:
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
            if os.path.isfile(fname) == True:
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
                # Store the information for later
                update_calibration_memory('wave_sol'+calib[0],[wavelength_solution, wavelength_solution_arclines])
            else:
                params['order_offset'] = [0,0]
                params['px_offset'] = [-60,60,6]
                params['px_offset_order'] = [-1,1,1]
                wavelength_solution, wavelength_solution_arclines = adjust_wavelength_solution(params, np.array(cal_l_spec), arc_lines_px, wavelength_solution, 
                                                                                               wavelength_solution_arclines, reference_catalog, reference_names, xlows, xhighs, params['GUI'])
            save_wavelength_solution_to_fits(wavelength_solution, wavelength_solution_arclines, params['master_wavelensolution'+calib[0]+'_filename'])
            plot_wavelength_solution_form(params['logging_wavelength_solution_form'], axlows, axhighs, wavelength_solution)
            plot_wavelength_solution_spectrum(cal_l_spec, cal_s_spec, params['logging_arc_line_identification_spectrum'], wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, plot_log=True)
            plot_wavelength_solution_image(im_cal_l, params['logging_arc_line_identification_positions'], cal_tr_poly, axlows, axhighs, wavelength_solution, wavelength_solution_arclines, reference_catalog)
        if params['original_master_wavelensolution_filename'].lower() != 'pseudo': 
            im_head, obsdate_midexp, obsdate_mid_float, jd_midexp = get_obsdate(params, im_arclhead)               # in UTC, mid of the exposure
            if calib[1] == 'cal2':          # for the calibration fiber
                wavelength_solution_cal, wavelength_solution_arclines_cal = copy.deepcopy(wavelength_solution), copy.deepcopy(wavelength_solution_arclines)
                params['wavelength_solution_type'] = 'cal-fiber'
                cal2_l_spec = copy.deepcopy(cal_l_spec)
                im_arc2lhead = copy.copy(im_arclhead)
                add_text_to_file('{0}\t{1}\t{2}\t{3}'.format(obsdate_mid_float, 0, 0, 'cal'), params['master_wavelengths_shift_filename'] )
            else:
                wavelength_solution_sci, wavelength_solution_arclines_sci = copy.deepcopy(wavelength_solution), copy.deepcopy(wavelength_solution_arclines)
                params['wavelength_solution_type'] = 'sci-fiber'
                cal1_l_spec = copy.deepcopy(cal_l_spec)
                im_arc1lhead = copy.copy(im_arclhead)
                add_text_to_file('{0}\t{1}\t{2}\t{3}'.format(obsdate_mid_float, 0, 0, 'sci'), params['master_wavelengths_shift_filename'] )
                params['two_solutions'] = True

    # Use the better wavelength solution: should be the one of the science fiber (wavelength_solution_sci), but if calibration fiber is better, then use calibration fiber
    if params['two_solutions']:       # solutions for both fibers
        if len(wavelength_solution_arclines_cal[wavelength_solution_arclines_cal > 100]) > 2* len(wavelength_solution_arclines_sci[wavelength_solution_arclines_sci > 100]):
            # if twice as many lines were identified (should also check that the residuals are not worse, but this is difficult as solution might have been loaded from file
            wavelength_solution, wavelength_solution_arclines = wavelength_solution_cal, wavelength_solution_arclines_cal           # Put the calibration fiber wavelength solution back as main solution
            params['wavelength_solution_type'] = 'cal-fiber'
            logger('Info: Using the calibration fiber wavelength solution (first solution) as master wavelength solution, as this solution seems to be better')
        else:
            logger('Info: Using the science fiber wavelength solution (second solution) as master wavelength solution')
    
    # Find the pixel-shift between the two solutions
    params['master_shift'] = 0
    if params['two_solutions']:
        params['extract_wavecal'] = True
        if wavelength_solution_cal.shape[0] != wavelength_solution_sci.shape[0]:
            logger('Error: The number of traces for the science ({0}) and calibration ({1}) wavelength solution differ. Please delete the wrong, old file ({2} or {3}'.format(\
                        wavelength_solution_sci.shape[0], wavelength_solution_cal.shape[0], params['master_wavelensolution_sci_filename'], params['master_wavelensolution_filename'] ))
        shift, shift_err = find_shift_between_wavelength_solutions(wavelength_solution_cal, wavelength_solution_arclines_cal, wavelength_solution_sci, wavelength_solution_arclines_sci, 
                                                                            np.zeros((wavelength_solution_cal.shape[0],im_trace1.shape[0])), ['calibration fiber','science fiber'] )
        # shift is positive if lines in science are right of lines in calibration
        """if params['wavelength_solution_type'] == 'sci-fiber':
            steps = ['cal-fiber']           # first do it for the not wavelength solution
        else:
            steps = ['sci-fiber']           # first do it for the not wavelength solution
        steps.append(params['wavelength_solution_type'])
        """
        if params['wavelength_solution_type'] == 'sci-fiber':
            aspectra = cal2_l_spec
            im_arclhead = im_arc2lhead
            im_name = 'long_exposure_emision_lamp_calibration_fiber'
            shift = -shift                          # reverse, as opposite to find_shift_between_wavelength_solutions
        else:
            aspectra = cal1_l_spec
            im_arclhead = im_arc1lhead
            im_name = 'long_exposure_emision_lamp_science_fiber'
        im_head, obsdate_midexp, obsdate_mid_float, jd_midexp = get_obsdate(params, im_arclhead)               # in UTC, mid of the exposure

        dummy_shift_wavesoln, master_shift, im_head = shift_wavelength_solution(params, aspectra, wavelength_solution, wavelength_solution_arclines, reference_catalog, 
                                                              reference_names, xlows, xhighs, obsdate_mid_float, jd_midexp, sci_tr_poly, cal_tr_poly, im_name, maxshift=max(2,2*shift_err), in_shift=shift, im_head=im_head )
        # master shift gives the shift from the wavelength_solution to the aspectra
        params['master_shift'] = master_shift
        #print("0:params['master_shift']", params['master_shift'])
        """ # Comparing the shifted solution with the original wavelength solution -> shifts in either direction, depending on polynomial, maximum shift: -> 0.025 \AA shift at 4880 -> 1.5 km/s
        wsci = create_wavelengths_from_solution(wavelength_solution_sci, cal2_l_spec)
        wcal = create_wavelengths_from_solution(wavelength_solution_cal, cal2_l_spec)
        wshift = create_wavelengths_from_solution(dummy_shift_wavesoln, cal2_l_spec)
        save_multispec([wsci, cal1_l_spec], 'fib1_wsci.fits', im_arclhead, bitpix=params['extracted_bitpix'])
        save_multispec([wsci, cal2_l_spec], 'fib2_wsci.fits', im_arclhead, bitpix=params['extracted_bitpix'])
        save_multispec([wcal, cal1_l_spec], 'fib1_wcal.fits', im_arclhead, bitpix=params['extracted_bitpix'])
        save_multispec([wcal, cal2_l_spec], 'fib2_wcal.fits', im_arclhead, bitpix=params['extracted_bitpix'])
        save_multispec([wshift, cal1_l_spec], 'fib1_wshift.fits', im_arclhead, bitpix=params['extracted_bitpix'])
        save_multispec([wshift, cal2_l_spec], 'fib2_wshift.fits', im_arclhead, bitpix=params['extracted_bitpix'])
        #"""
    params['extract_wavecal'] = False
    # Catch the problem, when the script re-runs with different settings and therefore the number of orders changes.
    if wavelength_solution.shape[0] != sci_tr_poly.shape[0]:
        #print('im_trace1.shape', im_trace1.shape)
        logger('Error: The number of traces for extraction and for the wavelength calibration do not match. Please remove eighter {0} ({2}) or {1} ({3}) and re-run the script in order to solve.'\
                    .format(params['master_trace_sci_filename'], params['master_wavelensolution_filename'], sci_tr_poly.shape[0], wavelength_solution.shape[0], im_trace1.shape[0]))
    
    im_blazecor, im_blazecor_head = create_image_general(params, 'blazecor')    # -> cont1cal2
    
    # Extract the flat spectrum and normalise it
    if os.path.isfile(params['master_blaze_spec_norm_filename']) == True:
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
        if entry.find('wavelengthcal2') >= 0 and entry.find('_rawfiles') >= 0:          # ThAr in calibration fiber
            wavelengthcals_cal.append(entry.replace('_rawfiles',''))
        elif entry.find('wavelengthcal') >= 0 and entry.find('_rawfiles') >= 0:         # ThAr in science fiber
            wavelengthcals_sci.append(entry.replace('_rawfiles',''))
            
    flat_spec_norm = np.array(fits.getdata(params['master_blaze_spec_norm_filename']))              # read it again, as the file is different than the data above
    number_no_data = np.sum(np.isnan(flat_spec_norm[2,:,:]), axis=1)                                # What orders don't contain enough data
    limit_no_data = (1-params['extraction_min_ident_part_of_trace_percent']/100.)*flat_spec_norm.shape[2]
    keep_orders  = np.where(number_no_data <= limit_no_data)[0]         # remove all orders where more than x% of the pixel are not available
    remove_orders = np.where(number_no_data > limit_no_data)[0]
    if len(remove_orders) > 0:
        logger('Warn: The blaze function for {0} orders contains not data points (limit: {2}). The following orders have been removed: {1}'.format(len(remove_orders), remove_orders, limit_no_data))
        #print wavelength_solution.shape, wavelength_solution_arclines.shape, sci_tr_poly.shape, xlows.shape, xhighs.shape, widths.shape, cal_tr_poly.shape, axlows.shape, axhighs.shape, awidths.shape, flat_spec_norm.shape
        wavelength_solution, wavelength_solution_arclines, sci_tr_poly, xlows, xhighs, widths, cal_tr_poly, axlows, axhighs, awidths, flat_spec_norm = \
                    wavelength_solution[keep_orders,:], wavelength_solution_arclines[keep_orders,:], \
                    sci_tr_poly[keep_orders,:,:], xlows[keep_orders], xhighs[keep_orders], widths[keep_orders,:], \
                    cal_tr_poly[keep_orders,:,:], axlows[keep_orders], axhighs[keep_orders], awidths[keep_orders,:], \
                    flat_spec_norm[:,keep_orders,:]                                                 # remove the bad orders
    
    if ( params['arcshift_side'] == 0 or params['two_solutions'] ) and len(wavelengthcals_cal)+len(wavelengthcals_sci) > 0:         # no calibration spectrum at the same time
        logger('Info: Starting to extract wavelength calibrations')
        params['extract_wavecal'] = True                                                                                # necessary for shift_wavelength_solution so the shift is stored in a file
        for [wavelengthcals, fib] in [ [wavelengthcals_cal,'cal'], [wavelengthcals_sci,'sci'] ]:
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
                                                    wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, flat_spec_norm, im_trace1, im_name)
    params['extract_wavecal'] = False
    if len(extractions) == 0:                                               # no extractions to do
        logger('Warn: Nothing to extract. -> Exiting')
        header_results_to_texfile(params)           # Save the results from the header in a logfile
        exit(0)
    logger('Info: Starting to extract spectra')
    for extraction in extractions:
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
                """print '!!!!!!!!! WARN: image is modified'
                print im.shape
                imtemp = copy.copy(im)
                imtemp[:-1,:] = im[1:,:]
                imtemp[ -1,:] = im[0 ,:]
                im = imtemp"""
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
                obj_names.append(obj_name.lower())
    
    logger('Info: Finished extraction of the science frames. The extracted {0}/*.fits file contains different data in a 3d array in the form: data type, order, and pixel. First data type is the wavelength (barycentric corrected), second is the extracted spectrum, followed by a measure of error. Forth and fith are the flat corrected spectra and its error. Sixth and sevens are the the continium normalised spectrum and the S/N in the continuum. Eight is the bad pixel mask, marking data, which is saturated or from bad pixel. Th last entry is the spectrum of the calibration fiber.'.format(params['path_extraction']))
    header_results_to_texfile(params)           # Save the results from the header in a logfile
    log_params(params)
        
    no_RV_names = ['flat', 'tung', 'whili', 'thar', 'th_ar', 'th-ar']
    # Do the Terra RVs
    if os.path.isfile(params['terra_jar_file']) == True and os.path.exists(params['path_rv_terra']):
        print('\nInfo: Preparing for the Terra analysis.')
        for root, dirs, files in os.walk(params['path_rv_terra'], followlinks=True):                       # Find all the objects again, as won't be added to obj_names when re-run
            for file in files:
                if file.endswith('.csv'):                       # has the file the correct ending?
                    filename = os.path.join(root, file).replace(params['path_rv_terra'],'')                # Only relative folder and filename
                    obj_name = filename.split(os.sep)[0]
                    if obj_name not in obj_names:
                        obj_names.append(obj_name)
        logger('For TERRA: changing directory to '+params['path_rv_terra'], show=False)
        os.chdir(params['path_rv_terra'])
        for obj_name in obj_names:
            do_RV = True
            for no_RV_name in no_RV_names:
                if obj_name.lower().find(no_RV_name) in [0,1,2,3,4,5]:
                    do_RV = False
                    break
            if not do_RV:
                continue
            newfile = 'echo "0998     synthetic         LAB                LAB                    0.0          0.0       0.0       {0}/" > astrocatalog.example'.format(obj_name)
            logger('For TERRA: creating a new astrocatalog.example with: '+newfile, show=False)
            os.system('rm -f astrocatalog.example; '+newfile)
            cmd = 'java -jar {1} -ASTROCATALOG astrocatalog.example 998 -INSTRUMENT CSV {0}'.format(wavelength_solution.shape[0],params['terra_jar_file'] )
            logger('For TERRA: running TERRA: '+cmd, show=False)
            if True:
                os.system(cmd)
            else:
                logger('Warn: TERRA commented out')
            resultfile = '{2} {1}{0}/results/synthetic.rv'.format(obj_name, params['path_rv_terra'], params['editor'])
            logger('For TERRA: results can be opened with: '+resultfile, show=False)
            
        os.chdir(params['path_run'])        # Go back to previous folder
        print('\nInfo: Finished the TERRA analysis. Some errors reported by TERRA are expected. The results are stored in {0}<object name>/results/synthetic.rv'.format(params['path_rv_terra']))
    
    # Do the Serval RVs
    if os.path.exists(params['path_serval']) and os.path.exists(params['path_serval']+'serval') and os.path.exists(params['path_serval']+'python'):
        pypath = ''
        if 'PYTHONPATH' in os.environ.keys():
            pypath = os.environ["PYTHONPATH"] + os.pathsep
        pypath += params['path_serval']+'python'
        logger('For SERVAL: set variable PYTHONPATH to {0}'.format(pypath))
        os.environ["PYTHONPATH"] = pypath
        xran = np.max(xhighs) - np.min(xlows)
        exclude_x = 0.1         # 0.1: 10% of outermost pixel on each side are excluded
        pmin = int(np.min(xlows) + 0.1*xran)
        pmax = int(np.max(xhighs) - 0.1*xran)
        oset = '{0}:{1}'.format(0,wavelength_solution.shape[0])
        logger('For SERVAL: changing directory to '+params['path_rv_serval'], show=False)
        os.chdir(params['path_rv_serval'])
        for obj_name in obj_names:
            do_RV = True
            for no_RV_name in no_RV_names:
                if obj_name.lower().find(no_RV_name) in [0,1,2,3,4,5]:
                    do_RV = False
                    break
            if not do_RV or not os.path.isfile('filelist_{0}.txt'.format(obj_name)):
                continue

            cmd = params['path_serval']+'serval/src/serval.py {3} filelist_{3}.txt -inst HIFLEX -targrv 0 -pmin {0} -pmax {1} -oset {2} -safemode'.format(pmin, pmax, oset, obj_name)
            logger('For SERVAL: running SERVAL: '+cmd, show=False)
            if True:
                os.system(cmd+'\n\n\n')
            else:
                logger('Warn: SERVAL commented out')
            resultfile = '{2} {1}{0}/{0}.rvc.dat'.format(obj_name, params['path_rv_serval'], params['editor'])  
            logger('For SERVAL: results can be opened with: '+resultfile, show=False)
        os.chdir(params['path_run'])        # Go back to previous folder
        print('\n\nInfo: Finished the SERVAL analysis. Some errors reported by SERVAL are expected. The results are stored in {0}<object name>/<object name>.rvc.dat'.format(params['path_rv_serval']))
    
""" ######### Analysis with terra
rm -f astrocatalog.example; echo "0998     synthetic         LAB                LAB                    0.0          0.0       0.0       HD82885/" > astrocatalog.example
rm -f astrocatalog.example; echo "0998     synthetic         LAB                LAB                    0.0          0.0       0.0       HD42807/" > astrocatalog.example
rm -f astrocatalog.example; echo "0998     synthetic         LAB                LAB                    0.0          0.0       0.0       HD140538/" > astrocatalog.example
java -jar /home/ronny/software/terra/PRV.jar -ASTROCATALOG astrocatalog.example 998 -INSTRUMENT CSV 55
cat HD*/results/synthetic.rv

############# Gnuplot analysis of part of the RV data:
f(x) = a2*x**2 + a1*x + a0; a2=0; a1=1; a0=0; fit f(x) 'data_1208' us 1:2 via a2,a1,a0
plot 'data_1208' us 1:2 , f(x)



"""
