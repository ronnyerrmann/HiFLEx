#!/usr/bin/env python

""" Checks which emission lines were used for different wavelength solutions - reads all master wavelength solutions given as parameter
"""
from procedures import *

# =============================================================================
# Define variables
# =============================================================================
#params = dict()     # default param dictionary
global calimages    # dictionary for all calibration images
# location of config file
#CONFIGFILE = 'conf.txt'

# Start of code
# deal with arguments from a text file
#params = textfileargs(params, CONFIGFILE)





if __name__ == "__main__":
    if len(sys.argv) < 3:
        print('Please start with parameter <percentage of appearances of the emission lines to be used in the new line list> and a list to all files of wavelength solutions to be used.')
        print('e.g. python {0} 40 ..{1}<folder*>{1}master_wavelength_sci.fits ..{1}<folder*>{1}master_wavelength_cal.fits'.format(sys.argv[0], os.sep))
        exit(1)
    try:
        percentage = float(sys.argv[1])
    except:
        print('Error: {0} cannot be transformed into a float. Please fix parameter <percentage of appearances of the emission lines to be used in the new line list>.'.format(sys.argv[1]))
        exit(1)
    filelist = sys.argv[2:]
    logger('Info: Using the lines that appear in {0}% of the following wavelength solutions: {1}'.format(percentage, filelist))
    # For now the orders need to cover the same wavelength, problem without that limit is that one wavelength can appear several times and at the boundaries it might not be good (data_o) - but then it would be nice to fit the same wavelengths fitted everywhere in the spectrum (data_g)
    
    for ii, filename in enumerate(filelist):
        wave_sol_dict = read_wavelength_solution_from_fits(filename)
        wave_px = wave_sol_dict['wavesol'][:,-2]        # Wavelength per pixel
        #print(wave_px/wave_sol_dict['wavesol'][:,-1])
        wave_px_rel = np.median(wave_px/wave_sol_dict['wavesol'][:,-1])    # "Resolution" (dl/dx / lambda)
        
        if ii == 0:         # First entry: get all lines
            real_orders = wave_sol_dict['wavesol'][:,0]       # Real order
            
            orders = wave_sol_dict['reflines'].shape[0]      # Number of orders
            lines = wave_sol_dict['reflines'].shape[1]       # Number of lines 
            data_o = np.zeros((len(filelist), orders, lines), dtype=float)      # first data set: split by order
            data_o[0,:,:] = wave_sol_dict['reflines']
            
            # Needs to be done in each step, in case of same_number_orders = False
            allwave = wave_sol_dict['reflines'][ wave_sol_dict['reflines']>100 ]    # 1D array of all real wavelengths (without fillings of 0)
            allwave = np.sort(allwave)                                  # Sorted, so no distingtion between orders
            """allwave = np.array([ 4027.00712631, 4027.0086481, 4027.0086481, 4027.34344568, 4027.34344568, 4027.54400778, 4027.651185, 5329.37326428, 5329.37761875, 5330.07784637, 5330.08059739, 5330.08059739, 5331.03806054, 5331.03806054, 6700.74893158, 6700.74942946, 6700.74942946, 6700.74942946, 6711.25435287, 6713.9702042, 6715.19138517, 6715.19138517, 6716, 6717, 6718, 6719 ])
            #allwave = np.array([ 4027.00712631, 4027.34344568, 4027.54400778, 4027.651185, 5329.37326428, 5330.07784637, 5330.08059739, 5331.03806054, 6700.74893158, 6711.25435287, 6713.9702042, 6715.19138517, 6716, 6717, 6718, 6719 ])
            allwave = np.array([ 5329.37761875, 4027.0086481, 6700.74942946, 4027.0086481, 4027.34344568, 6711.25435287, 4027.34344568, 4027.54400778, 4027.651185, 5329.37326428, 5330.07784637, 5330.08059739, 5330.08059739, 5331.03806054, 5331.03806054, 6700.74893158, 6700.74942946, 6700.74942946, 6713.9702042, 4027.00712631, 6715.19138517, 6715.19138517, 6716, 6717, 6718, 6719 ])"""
            
            replace_list, double_list = remove_double_entries_from_list(allwave, wave_px_rel, relative=True)
            
            for entry in replace_list:
                allwave[entry[0]] = np.median(allwave[entry])
            allwave[double_list] = 0
            allwave = allwave[allwave > 100]
            data_g = np.zeros((len(filelist), allwave.shape[0]), dtype=float)   # second data set: all lines independent of order
            data_g[0,:] = allwave
            
        else:           # further entries: add lines
            order_range = np.arrange( max(wave_sol_dict['wavesol'][0,0],real_orders[0]), min(wave_sol_dict['wavesol'][-1,0],real_orders[-1]) )
            if True:
                for real_order in order_range:
                    order1 = np.where(real_orders ==real_order)[0]
                    order2 = np.where(wave_sol_dict['wavesol'][:,0] ==real_order)[0]
                    for testwave in wave_sol_dict['reflines'][order2, wave_sol_dict['reflines'][order1,:]>0 ]:    # real wavelengths per order
                        diff = np.abs(testwave - data_o[:ii,order1,:])
                        pos_mins = np.argmin(diff, axis=1)
                        
                        print(pos_mins.shape, diff.shape, data_o.shape)
            
            # If oders are not covered: add them
            if (wave_sol_dict['wavesol'][0,0] < real_orders[0] or wave_sol_dict['wavesol'][-1,0] > real_orders[-1] ) and lines < wave_sol_dict['reflines'].shape[1]:
                data_o_add = np.zeros((len(filelist), orders, wave_sol_dict['reflines'].shape[1]-lines), dtype=float)
                data_o = np.concatenate((data_o, data_o_new), axis=2)
                lines = wave_sol_dict['reflines'].shape[1]       # Number of lines
            if wave_sol_dict['wavesol'][0,0] < real_orders[0]:       # Add redder orders
                data_o_new = np.zeros((len(filelist), real_orders[0]-wave_sol_dict['wavesol'][0,0], lines), dtype=float)      # first data set: split by order
                data_o_new[ii, :, :] = wave_sol_dict['reflines'][:data_o_new.shape[1]+1 ,:]
                data_o = np.concatenate((data_o_new, data_o), axis=1)
                real_orders = np.concatenate((np.arange(wave_sol_dict['wavesol'][0,0], real_orders[0]), real_orders), axis=0)
            if wave_sol_dict['wavesol'][-1,0] > real_orders[-1]:     # Add bluer orders
                data_o_new = np.zeros((len(filelist), wave_sol_dict['wavesol'][-1,0]-real_orders[-1], lines), dtype=float)      # first data set: split by order
                data_o_new[ii, :, :] = wave_sol_dict['reflines'][:data_o_new.shape[1]+1 ,:]
                data_o = np.concatenate((data_o_new, data_o), axis=1)
                real_orders = np.concatenate((real_orders,np.arange(real_orders[-1],wave_sol_dict['wavesol'][-1,0])), axis=0)
                
            # update lines, real_order













