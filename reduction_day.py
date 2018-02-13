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
    im_sflat, im_sflat_head = create_image_general(params, 'sflat')
    """zx """
    im_arc, im_archead = create_image_general(params, 'arc')
    im_arc_l, im_arclhead = create_image_general(params, 'arc_l')
    im_arc_s, im_arcshead = create_image_general(params, 'arc_s')
    im_flatarc, im_flatarchead = create_image_general(params, 'flatarc')
    
    reference_catalog, reference_names = read_reference_catalog(params['arc_lines_catalog'], params['catalog_file_wavelength_muliplier'], params['use_catalog_lines'])
    
    # Create or read the file with the orders for this night
    if os.path.isfile(params['master_order_filename']) == True:
        logger('Info: Using exiting trace solution: {0}'.format(params['master_order_filename']))
        sci_tr_poly, xlows, xhighs, widths = read_fits_width(params['master_order_filename'])        
        # plot_traces_over_image(im_sflat, params['logging_orders'], sci_tr_poly, xlows, xhighs, widths)        # Should already exist
    else:
        # load the original solution
        if os.path.isfile(params['original_master_order_filename']) == True:
            sci_tr_poly, xlows, xhighs, widths = read_fits_width(params['original_master_order_filename'])
            # find the shift between the original solution and the current flat
            shift, widths_new, shift_map, shift_error = shift_orders(im_sflat, params, sci_tr_poly, xlows, xhighs, widths, params['in_shift'])
            # save the map of the shifts
            save_im_fits(params, shift_map, im_sflat_head, params['logging_map_shift_orders'])
        else:
            shift_error = -1
        if shift_error > 1 or shift_error == -1 or abs(shift) > params['maxshift']:
            logger('Warn: The deviation of the shift of the orders seems too big or no previous solution was available, therefore running order_trace.py to find the real position of the orders')
            #os.sytem('mv {0} {1}'.format())#master_flat_filename already set in conf
            #ret = os.system('python {1}/order_trace.py original_master_order_filename={0}'.format(params['master_order_filename'], os.path.dirname(sys.argv[0])))
            #if ret > 0:
            #    exit(1)
            sci_tr_poly, xlows, xhighs, widths = trace_orders(params, im_sflat, im_sflat_head)
            #sci_tr_poly, xlows, xhighs, widths = read_fits_width(params['master_order_filename'])
        else:
            if params['update_widths'] == True:
                widths = widths_new
                logger('Info: widths of the traces have been updated')
            # update the sci_tr_poly parameters
            for pfit in sci_tr_poly:
               pfit[-1] += shift
            # save parameters of the polynoms into a fitsfile (from Neil)
            save_fits_width(sci_tr_poly, xlows, xhighs, widths, params['master_order_filename'])
            plot_traces_over_image(im_sflat, params['logging_orders'], sci_tr_poly, xlows, xhighs, widths)
            
    # Create the background map, if it doesn't exist
    if os.path.isfile(params['background_filename']) == True:
        logger('Info: Background map already exists: {0}'.format(params['result_path']+params['background_filename']))
    else:
        # create the background map
        if params['GUI'] == True:
            im_bck_px, params = bck_px_UI(params, im_sflat, sci_tr_poly, xlows, xhighs, widths, params['background_width_multiplier'], params['GUI'])
        else:
            im_bck_px = bck_px(im_sflat, sci_tr_poly, xlows, xhighs, widths, params['background_width_multiplier'])
        save_im_fits(params, im_bck_px, im_sflat_head, params['background_px_filename'])
        save_im_fits(params, im_bck_px*im_sflat, im_sflat_head, params['logging_background'])
        # Create the fitted background map
        im_bck = bck_fit(im_sflat, im_bck_px, params['polynom_bck'], params['GUI'])
        save_im_fits(params, im_bck, im_sflat_head, params['background_filename'])
        
    # Create the file for the arc orders, if it doesn't exist
    if os.path.isfile(params['master_orderarc_filename']) == True:
        logger('Info: Arc trace solution already exists: {0}'.format(params['master_orderarc_filename']))
        cal_tr_poly, axlows, axhighs, awidths = read_fits_width(params['master_orderarc_filename'])
    else:
        # use im_arc for automatic solution
        shifts = arc_shift(params, im_arc, sci_tr_poly, xlows, xhighs, widths)
        width_multiplier = 1
        # check the shift between the original solution and arc using a GUI -> This is outdated as it allows only one value for the shift
        shift = round(np.median(shifts),2)
        if params['GUI'] == True:
            shift_gui, width_multiplier = shift_orders_UI(im_arc, shift, sci_tr_poly, xlows, xhighs, widths)
            shifts += shift_gui - shift
        # update the sci_tr_poly parameters
        cal_tr_poly = []
        for order,pfit in enumerate(sci_tr_poly):
           new_pfit = list(pfit[:-1])
           new_pfit.append(pfit[-1]+shifts[order])
           cal_tr_poly.append(new_pfit)
        awidths = []
        for width in widths:
            new_width = []
            for w in width:
                new_width.append(w*width_multiplier)
            awidths.append(new_width)
        cal_tr_poly, awidths = np.array(cal_tr_poly), np.array(awidths)
        axlows = xlows
        axhighs = xhighs
        # save parameters of the polynoms into a fitsfile (from Neil)
        save_fits_width(cal_tr_poly, axlows, axhighs, awidths, params['master_orderarc_filename'])
        plot_traces_over_image(im_arc, params['logging_arcorders'], cal_tr_poly, axlows, axhighs, awidths)
    
    # Create the wavelength solution for the night
    if os.path.isfile(params['master_arc_solution_filename']) == True:
        logger('Info: wavelength solution already exists: {0}'.format(params['master_arc_solution_filename']))
        wavelength_solution, wavelength_solution_arclines = read_arc_fits(params['master_arc_solution_filename'])
        wavelength_solution = np.array(wavelength_solution)
    else:
        arc_l_spec, good_px_mask_l = extract_orders(params, im_arc_l, cal_tr_poly, axlows, axhighs, awidths, params['arcextraction_width_multiplier'])
        arc_s_spec, good_px_mask_s = extract_orders(params, im_arc_s, cal_tr_poly, axlows, axhighs, awidths, params['arcextraction_width_multiplier'])
        # Begin: This bit is not necessary once the procedure has been written
        im_name = 'master_arc_long'
        if 'master_arc_l_filename' in params.keys():
            im_name = params['master_arc_l_filename'].replace('.fits','')
        im_name = im_name.replace('.fit','')
        save_multispec([arc_l_spec, arc_l_spec, arc_l_spec, arc_l_spec], params['path_extraction']+im_name, im_arclhead)                    # This needs updating!!!
        im_name = 'master_arc_short'
        if 'master_arc_s_filename' in params.keys():
            im_name = params['master_arc_s_filename'].replace('.fits','')
        im_name = im_name.replace('.fit','')
        save_multispec([arc_s_spec, arc_s_spec, arc_s_spec, arc_s_spec], params['path_extraction']+im_name, im_arcshead)
        arc_lines_px = identify_lines(params, arc_l_spec, arc_s_spec, good_px_mask_l, new_format=True)
        if os.path.isfile(params['original_master_arc_solution_filename']) == False:                                                        # Create a new solution
            if os.path.isfile('arc_lines_wavelength.txt') == False:
                logger('Error: Files for creating the wavelength solution do not exist: {0}, {1}. Please check parameter {2} or create {1}.'.format(params['original_master_arc_solution_filename'], 'arc_lines_wavelength.txt', 'original_master_arc_solution_filename'))
            wavelength_solution, wavelength_solution_arclines = read_fit_wavelength_solution(params, 'arc_lines_wavelength.txt', im_arc_l)                        # For HARPS or a new wavelength solution
            save_arc_fits(wavelength_solution, wavelength_solution_arclines, params['original_master_arc_solution_filename'])                   # For HARPS or a new wavelength solution
            params['order_offset'] = [0,0]
            params['px_offset'] = [-20,20,10]
            params['px_offset_order'] = [-1,1,1]
        wavelength_solution_ori, wavelength_solution_arclines_ori = read_arc_fits(params['original_master_arc_solution_filename'])
        #wavelength_solution_ori = np.array(wavelength_solution_ori[::-1])      # if blue and red orders are swapped
        #wavelength_solution_ori[:,0] = np.abs(wavelength_solution_ori[:,0])    # if blue and red orders are swapped
        # Find the new wavelength solution
        wavelength_solution, wavelength_solution_arclines = adjust_wavelength_solution(params, np.array(arc_l_spec), arc_lines_px, wavelength_solution_ori, wavelength_solution_arclines_ori, reference_catalog, reference_names, xlows, xhighs, params['GUI'])
        save_arc_fits(wavelength_solution, wavelength_solution_arclines, params['master_arc_solution_filename'])
        plot_wavelength_solution_form(params['logging_wavelength_solution_form'], axlows, axhighs, wavelength_solution)
        plot_wavelength_solution_image(im_arc_l, params['logging_arc_line_identification_positions'], cal_tr_poly, axlows, axhighs, wavelength_solution, wavelength_solution_arclines, reference_catalog)
    
    # Extract the flat spectrum and normalise it
    if os.path.isfile(params['master_flat_spec_norm_filename']) == True:
        logger('Info: Normalised flat already exists: {0}, {1}'.format(params['master_flat_spec_norm_filename'], params['master_flat_spec_norm_filename']))
    else:
        logger('Step: Create the normalised flat for the night')
        shift = find_shift_images(im_flatarc, im_sflat, sci_tr_poly, cal_tr_poly)
        flat_spec, good_px_mask = extract_orders(params, im_flatarc, sci_tr_poly, xlows, xhighs, widths, params['extraction_width_multiplier'], var='prec', offset=shift)
        flatarc_spec, agood_px_mask = extract_orders(params, im_flatarc, cal_tr_poly, axlows, axhighs, awidths, params['arcextraction_width_multiplier'], offset=shift)
        flat_spec_norm = flat_spec/np.median(flat_spec[~np.isnan(flat_spec)])
        wavelength_solution_shift = shift_wavelength_solution(flatarc_spec, wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, xlows, xhighs)
        wavelengths = create_wavelengths_from_solution(wavelength_solution_shift, flatarc_spec)
        save_multispec([wavelengths, flat_spec_norm, flatarc_spec], params['master_flat_spec_norm_filename'], im_flatarchead)
        #save_im_fits(params, flat_spec_norm, im_sflat_head, params['master_flat_spec_norm_filename'])
    
    log_params(params)
    logger('Info: Finished routines for a new night of data. Now science data can be extracted. Please check before the output in the loging directory {0}: Are all orders identified correctly for science and calibration fiber, are the correct emission lines identified for the wavelength solution?\n'.format(params['logging_path']))
    
    extractions = []
    for entry in params.keys():
        if entry.find('extract') >= 0 and entry.find('_rawfiles') >= 0:
            extractions.append(entry.replace('_rawfiles',''))
    if len(extractions) > 0:
        logger('Info: Starting to extract spectra')
        flat_spec_norm = np.array(fits.getdata(params['master_flat_spec_norm_filename']))           # read it again, as the file is different than the data above
        for extraction in extractions:
            if  extraction.find('extract_combine') == -1:     # Single file extraction
                for im_name_full in params[extraction+'_rawfiles']:
                    im_name = im_name_full.rsplit('/')
                    im_name = im_name[-1].rsplit('.',1)         # remove the file ending
                    im_name = im_name[0]
                    if os.path.isfile(params['path_extraction']+im_name+'.fits'):
                        logger('Info: File {0} was already processed. If you want to extract again, please delete {1}{0}'.format(im_name, params['path_extraction']))
                    else:
                        print extraction, im_name_full, im_name
                        params['calibs'] = params[extraction+'_calibs_create']
                        im, im_head = read_file_calibration(params, im_name_full)
                        extraction_steps(params, im, im_name, im_head, sci_tr_poly, xlows, xhighs, widths, cal_tr_poly, axlows, axhighs, awidths, wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, flat_spec_norm, im_sflat)
            else:                                       # Combine files before extraction
                im_comb, im_comb_head = create_image_general(params, extraction)
                im_name = extraction
                extraction_steps(params, im_comb, im_name, im_comb_head, sci_tr_poly, xlows, xhighs, widths, cal_tr_poly, axlows, axhighs, awidths, wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, flat_spec_norm, im_sflat)
    
        logger('Info: Finished extraction of the science frames. The extracted {0}/*.fits file contains different data in a 3d array in the form: data type, order, and pixel. First data type is the wavelength, second is the extracted spectrum, followed by a measure of error (missing). Forth and fith are the flat corrected spectra and its error. Sixth and sevens are the the continium normalised spectrum and the S/N in the continuum. Eight is the arc spectrum.'.format(params['path_extraction']))
    
    
    
    
