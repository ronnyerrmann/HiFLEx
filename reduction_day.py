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

# Best solution:
# Arcs: long exposure and short exposures to have more lines to identify lines
# Flat: no saturation for tracing and background solution
# FlatArc: no saturation for the normalised flat

# Start of code
# deal with arguments from a text file
params = textfileargs(params, CONFIGFILE)
if __name__ == "__main__":
    logger('Info: Starting routines for a new night of data, including: shift orders, find background, find arc orders, create wavelength solution, and create normalised extracted flat')
    log_params(params)
    
    # create the median combined files
    im_trace1, im_trace1_head = create_image_general(params, 'trace1')
    if params['arcshift_side'] <> 0 or 'trace2_rawfiles' in params.keys():      # calibration spectrum at the same time
        im_trace2, im_trace2_head = create_image_general(params, 'trace2')
    else:                                                                       # no calibration spectrum at the same time
        im_trace2, im_trace2_head = im_trace1, im_trace1_head
    #flatarc, cal2_l, cal2_s: at a later step to know the orders for localbackground -> cont1cal2
    
    # Load the reference catalogue. Note that the Ar I lines are rescaled!
    reference_catalog, reference_names = read_reference_catalog(params['reference_catalog'], params['catalog_file_wavelength_muliplier'], params['use_catalog_lines'])
    
    # Create or read the file with the orders for this night
    printresults = False
    if os.path.isfile(params['master_trace_sci_filename']) == True:
        logger('Info: Using exiting trace solution: {0}'.format(params['master_trace_sci_filename']))
        sci_tr_poly, xlows, xhighs, widths = read_fits_width(params['master_trace_sci_filename'])
    else:
        # load the original solution
        if os.path.isfile(params['original_master_traces_filename']) == True:
            sci_tr_poly, xlows, xhighs, widths = read_fits_width(params['original_master_traces_filename'])
            # find the shift between the original solution and the current flat
            shift, widths_new, shift_map, shift_error = shift_orders(im_trace1, params, sci_tr_poly, xlows, xhighs, widths, params['in_shift'])
            # save the map of the shifts
            save_im_fits(params, shift_map, im_trace1_head, params['logging_map_shift_orders'])
        else:
            shift_error = -1
        if shift_error > 1 or shift_error == -1 or abs(shift) > params['maxshift']:
            logger('Warn: The deviation of the shift of the orders seems too big or no previous solution was available, therefore searching for the position of the orders from scratch:')
            sci_tr_poly, xlows, xhighs, widths = trace_orders(params, im_trace1, im_trace1_head)
            printresults = True
        else:
            # update the sci_tr_poly parameters
            sci_tr_poly[:,:,-1] += shift
            if params['update_widths'] == True:
                for order in range(sci_tr_poly.shape[0]):
                    xarr = range(xlows[order], xhighs[order])
                    center = np.polyval(sci_tr_poly[order, 0, 1:], xarr-sci_tr_poly[order, 0, 0])
                    left   = np.polyval(sci_tr_poly[order, 1, 1:], xarr-sci_tr_poly[order, 1, 0])
                    right  = np.polyval(sci_tr_poly[order, 2, 1:], xarr-sci_tr_poly[order, 2, 0])
                    print widths_new[order,0]/widths[order,0], widths_new[order,1]/widths[order,1], widths_new[order,0],widths[order,0], widths_new[order,1],widths[order,1]
                    left, right = adjust_width_orders(center, left, right, [widths_new[order,0]/widths[order,0], widths_new[order,1]/widths[order,1] ] )
                    sci_tr_poly[order, 2, 1:] = np.polyfit(xarr - sci_tr_poly[order, 2, 0], right, len(sci_tr_poly[order, 2, 1:])-1 )
                widths = widths_new
                logger('Info: widths of the traces have been updated')
            printresults = True
        if printresults:
            # save parameters of the polynoms into a fitsfile (from Neil)
            save_fits_width(sci_tr_poly, xlows, xhighs, widths, params['master_trace_sci_filename'])
            # Produce some useful statistics
            plot_traces_over_image(im_trace1, params['logging_traces_im'], sci_tr_poly, xlows, xhighs, widths)
            data = np.insert(widths, 0, range(len(sci_tr_poly)), axis=1)       # order, left, right, gausswidth
            positio, pctlwidth = [], []
            for order in range(sci_tr_poly.shape[0]):                      # For the central data
                xarr = range(xlows[order],xhighs[order])
                positio.append(np.polyval(sci_tr_poly[order, 0, 1:], im_trace1.shape[0]/2 - sci_tr_poly[order, 0, 0]))
                pctlwidth.append( np.median ( np.polyval(sci_tr_poly[order, 2, 1:], xarr - sci_tr_poly[order, 2, 0]) - np.polyval(sci_tr_poly[order, 1, 1:], xarr - sci_tr_poly[order, 1, 0]) ) )   # median pctlwidth
            data = np.append(data, np.array(pctlwidth)[:,None],axis=1)
            data = np.append(data, np.array(positio)[:,None],axis=1)
            data = np.append(data, np.array(xlows)[:,None],axis=1)
            data = np.append(data, np.array(xhighs)[:,None],axis=1)
            printarrayformat = ['%1.1i','%3.1f', '%3.1f', '%4.2f\t','%4.2f\t', '%1.1i','%1.1i','%1.1i']
            logger('\t\torder\tleft\tright\tgausswidth\tpctlwidth\tpositio\tmin_tr\tmax_tr\t(positio: position of the trace at the center of the image, pctlwidth: median full width of the trace at {0}% of maximum)'\
                      .format(params['width_percentile']),printarrayformat=printarrayformat, printarray=data)
        
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
        width_multiplier = 1
        # check the shift between the original solution and arc using a GUI -> This is outdated as it allows only one value for the shift
        shift = round(np.median(shifts),2)
        if params['GUI'] == True:
            shift_gui, width_multiplier = shift_orders_UI(im_trace2, shift, sci_tr_poly, xlows, xhighs, widths)
            shifts += shift_gui - shift
        # update the sci_tr_poly parameters
        cal_tr_poly = []
        for order in range(sci_tr_poly.shape[0]):
            new_pfit = []
            for dataset in range(sci_tr_poly.shape[1]):
                new_pfit.append(list(sci_tr_poly[order, dataset, :-1]) + [sci_tr_poly[order, dataset, -1]+shifts[order]] )
            cal_tr_poly.append(new_pfit)
        cal_tr_poly, awidths, axlows, axhighs = np.array(cal_tr_poly), copy.deepcopy(widths), copy.deepcopy(xlows), copy.deepcopy(xhighs)
        # save parameters of the polynoms into a fitsfile (from Neil)
        save_fits_width(cal_tr_poly, axlows, axhighs, awidths, params['master_trace_cal_filename'])
        plot_traces_over_image(im_trace2, params['logging_arctraces_im'], cal_tr_poly, axlows, axhighs, awidths)
    
    # Catch the problem, when the script re-runs with different settings and therefore the number of orders changes.
    if cal_tr_poly.shape[0] <> sci_tr_poly.shape[0]:
        logger('Error: The number of traces for the science fiber and for the calibration fiber do not match. Please remove eighter {0} or {1} and re-run the script in order to solve.'.format(params['master_trace_cal_filename'], params['master_trace_sci_filename']))
    
    update_calibration_memory('sci_trace',[sci_tr_poly, xlows, xhighs, widths])         # Apertures might be shifted before extraction -> this would also affect the localbackground
    update_calibration_memory('cal_trace',[cal_tr_poly, axlows, axhighs, awidths])
    
    for calib in [ ['','cal2'], ['_sci','cal1'] ]:
        # first entry is the standard wavelength solution
        # second entry is for finding the wavelength solution in a bifurcated fiber
        if calib[1] == 'cal1' and ('cal1_l_rawfiles' not in params.keys() ):
            break                                           # Not set up for bifurcated fiber use -> stop after the first step
        elif calib[1] == 'cal1':                            # Update and rename a few parameters
            if 'master_wavelensolution'+calib[0]+'_filename' not in params.keys():
                params['master_wavelensolution'+calib[0]+'_filename'] = params['master_wavelensolution_filename'].replace('.fit',calib[0]+'.fit')
            for pngparam in ['logging_wavelength_solution_form', 'logging_em_lines_gauss_width_form', 'logging_arc_line_identification_residuals', 'logging_arc_line_identification_positions']:
                params[pngparam] = params[pngparam].replace('.png','')+calib[0]+'.png'
            for pdfparam in ['logging_arc_line_identification_spectrum']:
                params[pdfparam] = params[pdfparam].replace('.pdf','')+calib[0]+'.pdf'
        # Create the wavelength solution for the night
        if params['original_master_wavelensolution_filename'].lower() <> 'pseudo':                  # Create the master files
            im_cal_l, im_arclhead = create_image_general(params, calib[1]+'_l')           # arc_l -> cal2_l
            im_cal_s, im_arcshead = create_image_general(params, calib[1]+'_s')           # arc_s -> cal2_s
        #wavelength_solution, wavelength_solution_arclines = read_fit_wavelength_solution(params, 'arc_lines_wavelength.txt', im_arc_l)         # Testing
        if os.path.isfile(params['master_wavelensolution'+calib[0]+'_filename']) == True:
            logger('Info: wavelength solution already exists: {0}'.format(params['master_wavelensolution'+calib[0]+'_filename']))
            wavelength_solution, wavelength_solution_arclines = read_wavelength_solution_from_fits(params['master_wavelensolution_filename'])
        elif params['original_master_wavelensolution_filename'].lower() == 'pseudo':
            logger('Warning: Using a pseudo solution for the wavelength (1 step per px)')
            wavelength_solution, wavelength_solution_arclines = create_pseudo_wavelength_solution(sci_tr_poly.shape[0])
        else:
            cal_l_spec, good_px_mask_l = extract_orders(params, im_cal_l, cal_tr_poly, axlows, axhighs, awidths, params['arcextraction_width_multiplier'])
            cal_s_spec, good_px_mask_s = extract_orders(params, im_cal_s, cal_tr_poly, axlows, axhighs, awidths, params['arcextraction_width_multiplier'])
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
            if os.path.isfile(params['logging_found_arc_lines']+calib[0]) == True:
                logger('Info: List of the identified emission lines already exists. Using the information from file: {0}'.format(params['logging_found_arc_lines']+calib[0] ))
                arc_lines_px_txt = read_text_file(params['logging_found_arc_lines']+calib[0], no_empty_lines=True)              # list of strings, first entry is header 
                arc_lines_px = np.array( convert_readfile(arc_lines_px_txt[1:], [int, float, float, float], delimiter='\t', replaces=['\n',' ']) )
            else:
                arc_lines_px = identify_lines(params, cal_l_spec, cal_s_spec, good_px_mask_l, good_px_mask_s)
                logger('Info: Identified {0} lines in the arc spectrum. These lines are stored in file {1}'.format(len(arc_lines_px), params['logging_found_arc_lines']+calib[0] ))
                printarrayformat = ['%1.1i', '%3.2f', '%3.2f', '%3.1f']
                logger('order\tpixel\twidth\theight of the line', show=False, printarrayformat=printarrayformat, printarray=arc_lines_px, logfile=params['logging_found_arc_lines']+calib[0])
    
            if calib[1] == 'cal2':          # for the calibration fiber
                if os.path.isfile(params['original_master_wavelensolution_filename']) == False:                                                         # Create a new solution
                    if os.path.isfile('arc_lines_wavelength.txt') == False:                                                                             # No arc_lines_wavelength.txt available
                        wavelength_solution, wavelength_solution_arclines = create_pseudo_wavelength_solution(cal_l_spec.shape[0])                      # Create a pseudo solution
                        plot_wavelength_solution_spectrum(cal_l_spec, cal_s_spec, params['logging_arc_line_identification_spectrum'].replace('.pdf','')+'_manual.pdf', 
                                                      wavelength_solution, wavelength_solution_arclines, np.array([0,1,0]).reshape(1,3), ['dummy'], plot_log=True)     # Plot the spectrum
                        logger('Error: Files for creating the wavelength solution do not exist: {0}, {1}. Please check parameter {2} or create {1}.'.format(\
                                                params['original_master_wavelensolution_filename'], 'arc_lines_wavelength.txt', 'original_master_wavelensolution_filename'))
                    wavelength_solution, wavelength_solution_arclines = read_fit_wavelength_solution(params, 'arc_lines_wavelength.txt', im_cal_l)         # For a new wavelength solution
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
            else:
                params['order_offset'] = [0,0]
                params['px_offset'] = [-10,10,2]
                params['px_offset_order'] = [-1,1,1]
                wavelength_solution, wavelength_solution_arclines = adjust_wavelength_solution(params, np.array(cal_l_spec), arc_lines_px, wavelength_solution, 
                                                                                               wavelength_solution_arclines, reference_catalog, reference_names, xlows, xhighs, params['GUI'])
            save_wavelength_solution_to_fits(wavelength_solution, wavelength_solution_arclines, params['master_wavelensolution'+calib[0]+'_filename'])
            plot_wavelength_solution_form(params['logging_wavelength_solution_form'], axlows, axhighs, wavelength_solution)
            plot_wavelength_solution_spectrum(cal_l_spec, cal_s_spec, params['logging_arc_line_identification_spectrum'], wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, plot_log=True)
            plot_wavelength_solution_image(im_cal_l, params['logging_arc_line_identification_positions'], cal_tr_poly, axlows, axhighs, wavelength_solution, wavelength_solution_arclines, reference_catalog)
    
    # Create the wavelength solution for the science fiber, if using a bifurcated fiber
    
            
    # Catch the problem, when the script re-runs with different settings and therefore the number of orders changes.
    if wavelength_solution.shape[0] <> sci_tr_poly.shape[0]:
        logger('Error: The number of traces for extraction and for the wavelength calibration do not match. Please remove eighter {0} ({2}) or {1} ({3}) and re-run the script in order to solve.'\
                    .format(params['master_trace_sci_filename'], params['master_wavelensolution_filename'], sci_tr_poly.shape[0], wavelength_solution.shape[0] ))
        
    update_calibration_memory('wave_sol',[wavelength_solution, wavelength_solution_arclines])
    
    im_flatarc, im_flatarc_head = create_image_general(params, 'flatarc')    # -> cont1cal2
    
    # Extract the flat spectrum and normalise it
    if os.path.isfile(params['master_flat_spec_norm_filename']) == True:
        logger('Info: Normalised flat already exists: {0}'.format(params['master_flat_spec_norm_filename']))
        # The file is read later on purpose
    else:
        logger('Step: Create the normalised flat for the night')
        obsdate, obsdate_float = get_obsdate(params, im_flatarc_head)
        shift = find_shift_images(params, im_flatarc, im_trace1, sci_tr_poly, xlows, xhighs, widths, 1, cal_tr_poly, extract=True)
        flat_spec, good_px_mask = extract_orders(params, im_flatarc, sci_tr_poly, xlows, xhighs, widths, params['extraction_width_multiplier'], var='prec', offset=shift)
        flat_spec_norm = flat_spec/np.nanmedian(flat_spec)
        flat_spec_norm_cor = correct_blaze(flat_spec_norm, minflux=0.1)
        if np.nansum(np.abs(sci_tr_poly - cal_tr_poly)) == 0.0:                                                 # science and calibration traces are at the same position
            flatarc_spec, agood_px_mask = flat_spec*0, copy.copy(good_px_mask)
        else:
            flatarc_spec, agood_px_mask = extract_orders(params, im_flatarc, cal_tr_poly, axlows, axhighs, awidths, params['arcextraction_width_multiplier'], offset=shift)
        wavelength_solution_shift = shift_wavelength_solution(params, flatarc_spec, wavelength_solution, wavelength_solution_arclines, 
                                            reference_catalog, reference_names, xlows, xhighs, obsdate_float, sci_tr_poly, cal_tr_poly, params['master_flat_spec_norm_filename'])
        wavelengths = create_wavelengths_from_solution(wavelength_solution_shift, flatarc_spec)
        save_multispec([wavelengths, flat_spec_norm, flat_spec_norm_cor, flatarc_spec], params['master_flat_spec_norm_filename'], im_flatarc_head)
        #save_im_fits(params, flat_spec_norm, im_sflat_head, params['master_flat_spec_norm_filename'])
    
    log_params(params)
    logger('Info: Finished routines for a new night of data. Now science data can be extracted. Please check before the output in the loging directory {0}: Are all orders identified correctly for science and calibration fiber, are the correct emission lines identified for the wavelength solution?\n'.format(params['logging_path']))
    
    obj_names = []
    extractions = []
    wavelengthcals = []
    for entry in params.keys():
        if entry.find('extract') >= 0 and entry.find('_rawfiles') >= 0:
            extractions.append(entry.replace('_rawfiles',''))
        if entry.find('wavelengthcal') >= 0 and entry.find('_rawfiles') >= 0:
            wavelengthcals.append(entry.replace('_rawfiles',''))
    if len(extractions) == 0:                                               # no extractions to do
        exit(0)
    flat_spec_norm = np.array(fits.getdata(params['master_flat_spec_norm_filename']))           # read it again, as the file is different than the data above
    if params['arcshift_side'] == 0:                                            # no calibration spectrum at the same time
        logger('Info: Starting to extract wavelength calibrations')
        for wavelengthcal in wavelengthcals:
            for im_name_full in params[wavelengthcal+'_rawfiles']:
                im_name = im_name_full.rsplit('/')
                im_name = im_name[-1].rsplit('.',1)         # remove the file ending
                im_name = im_name[0]
                im_name_wc = im_name+'_wavecal'
                if os.path.isfile(params['path_extraction']+im_name_wc+'.fits'):
                    logger('Info: File {0} was already processed for the calibration of the wavelength solution. If you want to extract again, please delete {1}{0}.fits'.format(im_name_wc, params['path_extraction']))
                    continue
                params['calibs'] = params[wavelengthcal+'_calibs_create']
                im, im_head = read_file_calibration(params, im_name_full)
                extraction_wavelengthcal(params, im, im_name_wc, im_head, sci_tr_poly, xlows, xhighs, widths, cal_tr_poly, axlows, axhighs, awidths, \
                                                    wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, flat_spec_norm, im_trace1, im_name)
    logger('Info: Starting to extract spectra')
    for extraction in extractions:
        if  extraction.find('extract_combine') == -1:     # Single file extraction
            for im_name_full in params[extraction+'_rawfiles']:
                im_name = im_name_full.rsplit('/')
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
            obj_name = extraction_steps(params, im_comb, im_name, im_comb_head, sci_tr_poly, xlows, xhighs, widths, cal_tr_poly, axlows, axhighs, awidths, \
                                        wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, flat_spec_norm, im_trace1)
            if obj_name not in obj_names:
                obj_names.append(obj_name)
    
    logger('Info: Finished extraction of the science frames. The extracted {0}/*.fits file contains different data in a 3d array in the form: data type, order, and pixel. First data type is the wavelength (barycentric corrected), second is the extracted spectrum, followed by a measure of error. Forth and fith are the flat corrected spectra and its error. Sixth and sevens are the the continium normalised spectrum and the S/N in the continuum. Eight is the bad pixel mask, marking data, which is saturated or from bad pixel. Th last entry is the spectrum of the calibration fiber.'.format(params['path_extraction']))
    
    # Do the Terra RVs
    os.system('echo "0998     synthetic         LAB                LAB                    0.0          0.0       0.0       Object1/" > astrocatalog.example')
    if os.path.isfile('/home/ronny/software/terra/PRV.jar') == True:
        os.system('java -jar ~/software/terra/PRV.jar -ASTROCATALOG astrocatalog.example 998 -INSTRUMENT CSV {0}'.format(wavelength_solution.shape[0]) )
        os.system('gedit Object1/results/synthetic.rv &')
        
        
""" Gnuplot analysis of part of the RV data:
f(x) = a2*x**2 + a1*x + a0; a2=0; a1=1; a0=0; fit f(x) 'data_1208' us 1:2 via a2,a1,a0
plot 'data_1208' us 1:2 , f(x)



"""
