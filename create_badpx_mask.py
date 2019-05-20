import numpy as np
import os
#import sys
from procedures import *
from tqdm import tqdm
import copy
import random

# =============================================================================
# Define variables
# =============================================================================
params = dict()     # default param dictionary
calimages = dict()  # dictionary for all calibration images
# location of config file
CONFIGFILE = 'conf.txt'

#params['exptimes'] = [0.1,0.2,0.3,0.5,0.7, 1,1.5,2,3,4,5,6,8,10,12]
params['exptimes'] = [0.1,0.3,0.5,0.7,0.9, 1.1,1.3,1.5,1.7,1.9, 2.1,2.3,2.5,2.7,2.9]                        #21/7
params['exptimes'] = [0.1,0.4,0.7, 1,1.2,1.4,1.6,1.8, 2,2.2,2.4,2.6,2.8, 3,3.2,3.4,3.6,3.8, 4,4.2,4.5]      #5/12
#params['exptimes'] = [1,2,3,3.4,3.8]      #5/12
#params['exptimes'] = [3.8]      #5/12
#params['exptimes'] = [0.01,0.04,0.07,0.1,0.4,0.5, 1,1.3,1.7, 2,2.4,2.8,3.2,3.6, 4,4.5,5]      #5/1
#params['exptimes'] = [0.1,1,3.2]      #5/1
params['exptimes'] = [0.1,0.15,0.2,0.3,0.5,0.7, 1,1.5,2,3,5,7, 10,15,16,17,18,19,20]     # 20180924
params['exptimes'] = [0.001,0.002,0.003,0.005,0.007, 0.01,0.02,0.03,0.05,0.07, 0.1,0.2,0.3,0.5,0.7, 1,2,3,5,7, 10]     # 20190110
params['exptimes'] = [0.01,0.02,0.03,0.05,0.07, 0.1,0.2,0.3,0.5,0.7, 1,2,3,5,7, 10]     # 20190110
params['exptimes'] = [0.001,0.01,0.02,0.03,0.05,0.07, 0.1,0.2,0.3,0.5,0.7, 1,1.5,2,2.5,3,3.5,4,5,6,7]                  # 20190116
#params['exptimes'] = [0.5,0.7, 1,1.5,2,2.5,3,3.5,4,5,6]                  # 20190116
params['exptimes'] = [1,10,20,40,60]                  # 20190116
params['exptimes'] = [0.01, 0.02,0.03,0.04,0.06,0.08, 0.1,0.12,0.14,0.16,0.18, 0.2,0.22,0.24,0.26,0.28, 0.3,0.4,0.5,0.6,0.7,0.8,0.9, 1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8, 2,2.2,2.4,2.6,2.8, 3,3.2,3.4,3.6,3.8, 4,4.2,4.4,4.6,4.8, 5,5.2,5.4,5.6,5.8, 6,6.2,6.4,6.6,6.8, 7,7.2,7.4,7.6]                  # 20190116
#params['exptimes'] = []
range_bias = range(1,52)
range_dark = range(1,8)
range_flat = range(1,6)
#range_flat = range(0,0)
show_stats = ['bias', 'dark', 'flat']
show_stats = ['flat']
show_stats = ['dark']
show_stats = []
params['gains'] = 'gains.fits'
params['zerop'] = 'zerop.fits'
params['xord'], params['yord'] = 5, 5       # for qsi camera yord is along dispersion axis, 3,4 needs about 3GB memory, 4,4 about 3.5 GB; 1,0 is for bias/dark; 4,4 for flats, gain
params['xord'], params['yord'] = 2, 2       # for biases/darks

#get_statistics = ['Flat-0001_0p001s.fit', 'Flat-0002_0p001s.fit', 'Flat-0003_0p001s.fit', 'Flat-0004_0p001s.fit', 'Flat-0005_0p001s.fit', 'Flat-0001b_0p001s.fit', 'Flat-0002b_0p001s.fit', 'Flat-0003b_0p001s.fit', 'Flat-0004b_0p001s.fit', 'Flat-0005b_0p001s.fit', 'Flat-0001c_0p001s.fit', 'Flat-0002c_0p001s.fit', 'Flat-0003c_0p001s.fit', 'Flat-0004c_0p001s.fit', 'Flat-0005c_0p001s.fit']
get_statistics = []

# Start of code
# deal with arguments from a text file
params = textfileargs(params, CONFIGFILE)

def find_stats(im_stats):
    if len(im_stats) == 0:
        return
    stat_single = []
    for im in tqdm(im_stats):
        im_fit = fit_2d_image(im, params['xord'], params['yord'])
        im_diff = im - im_fit
        #plot_img_spec.plot_image(im_diff, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'residuals between fit of the gain and the gain')
        list_diff = im_diff.flatten()
        bins = int( ( max(list_diff) - min(list_diff) ) / 25 + 1)        # 25 ADU wide bins
        hist, bin_edges = np.histogram(list_diff, bins=bins)
        #print hist,bin_edges, hist.shape,bin_edges.shape
        #plot_img_spec.plot_points([bin_edges], [hist], ['histogram'], '', show=True, return_frame=False, x_title='flux [ADU]', y_title='number')
        index = np.argmax(hist)
        indexl = -1
        indexr = len(bin_edges)-1       # -1 because otherwise it would try to read the index after the last one
        for i in range(index-1)[::-1]:
            if hist[i] == 0:
                indexl = i
                break
        for j in range(index+1,len(hist)):
            if hist[j] == 0:
                indexr = j
                break
        #print indexl,indexr, len(list_diff), len(hist), len(bin_edges), bin_edges[indexl+1], bin_edges[indexr]
        list_diff = list_diff[list_diff >= bin_edges[indexl+1] ]
        list_diff = list_diff[list_diff <= bin_edges[indexr] ]
        #print len(list_diff)
        stat_single.append([len(stat_single), np.mean(im), np.median(im), np.std(im, ddof=1), np.median(im_diff), np.std(im_diff, ddof=1), np.std(list_diff, ddof=1)])
    printarrayformat = ['%1.1i', '%3.1f', '%3.1f', '%4.2f', '%3.1f', '%4.2f', '%4.2f']
    logger('The individual files have the following properties (values given in ADU, the first set without a fit of the data, the second is on the residuals of the data, the third is the residuals, cleared by outlier):\n\t\tindex\taverage\tmedian\tstdev\tmedian\tstdev\tstdev',printarrayformat=printarrayformat, printarray=stat_single)
    stat_diff = []
    for i in range(len(im_stats)-1):
        for j in range(i+1,len(im_stats)):
            im = im_stats[i]-im_stats[j]
            stat_diff.append([i,j, np.mean(im), np.median(im), np.std(im, ddof=1)])
    printarrayformat = ['%1.1i', '%1.1i', '%4.2f', '%4.2f', '%5.3f']
    logger('The difference between 2 files have the following properties (values given in ADU, stdev defines the readout noise):\n\t\tindex1\tindex2\taverage\tmedian\tstdev',printarrayformat=printarrayformat, printarray=stat_diff)

def plot_linearity_UI(exp_times, data, title='', adjust=[0.07,0.90,0.94,0.06, 1.0,1.01]):
    
    fig, frame = plt.subplots(1, 1)
    plt.subplots_adjust(left=adjust[0], right=adjust[1], top=adjust[2], bottom=adjust[3])
    
    datac = np.median(data, axis=(1,2))
        
    def plot(frame, x, data, datac):
        frame.clear()
        datas = data.shape
        try:
            option = gui3.data['option_plot']
        except:
            option = 0
        try:
            number_graphs = gui3.data['number_graphs']
        except:
            number_graphs = 5           # not yet initialised
        for number in range(number_graphs):
            i, j = np.random.randint(0,datas[1]+1), np.random.randint(0,datas[2]+1)     # get random pixel
            if option == 0:
                frame.plot((data[1:,i,j]-data[0,i,j]), (data[1:,i,j]-data[0,i,j])/(x[1:]-x[0]), label='{0},{1}'.format(i,j) )
            elif option == 1:
                frame.plot(x, data[:,i,j], label='{0},{1}'.format(i,j) )
            elif option == 2:
                frame.plot(x, data[:,i,j]/x, label='{0},{1}'.format(i,j) )
            elif option == 3:
                frame.plot(x[1:], (data[1:,i,j]-data[0,i,j])/(x[1:]-x[0]), label='{0},{1}'.format(i,j) )
        if option == 0:
            frame.plot((datac[1:]-datac[0]), (datac[1:]-datac[0])/(x[1:]-x[0]), label='all', linewidth=3, color='black' )
        elif option == 1:
            frame.plot(x, datac[:], label='all', linewidth=3, color='black' )
        elif option == 2:
            frame.plot(x, datac[:]/x, label='all', linewidth=3, color='black' )
        elif option == 3:
            frame.plot(x[1:], (datac[1:]-datac[0])/(x[1:]-x[0]), label='all', linewidth=3, color='black' )
        #xlabel_text = 'exposure time [s]'
        xlabel_text = ['zero-corrected flux [ADU]'                           , 'exposure time [s]', 'exposure time [s]'                     , 'exposure time [s]'                                    ]
        ylabel_text = ['zero-corrected flux devided by exposure time [ADU/s]', 'flux [ADU]'       , 'flux divided by exposure time [ADU/s]' , 'zero-corrected flux divided by exposure time [ADU/s]' ]
        frame.set_xlabel(xlabel_text[option], fontsize=14)
        frame.set_ylabel(ylabel_text[option], fontsize=14)
        frame.set_title(title, fontsize=16)
        frame.legend(loc='upper left', bbox_to_anchor=(adjust[4], adjust[5]))
    # get kwargs
    pkwargs = dict(frame=frame, x=exp_times, data=data, datac=datac)
    # run initial update plot function
    plot(**pkwargs)
    
    # define valid_function
    # input is one variable (the string input)
    # return is either:
    #   True and values
    # or
    #   False and error message
    def vfunc_int(xs):
        try:
            value = int(xs)
            return True, value
        except:
            return False, ('Error, input must be integer')
    # define widgets
    widgets = dict()
    starta = 5
    startb = 0
    widgets['number_graphs'] = dict(label='Number of\nGraphs?', comment=None, #'integer',
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=str, start=starta, valid_function=vfunc_int,
                                width=10)
    widgets['option_plot'] = dict(label='Option', comment='0, 1, 2, or 3',
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=str, start=startb, valid_function=vfunc_int,
                                width=10)
    widgets['accept'] = dict(label='Close', kind='ExitButton', position=Tk.BOTTOM)
    widgets['update'] = dict(label='Update', kind='UpdatePlot', position=Tk.BOTTOM)

    wprops = dict(orientation='v', position=Tk.RIGHT)

    gui3 = tkc.TkCanvas(figure=fig, ax=frame, func=plot, kwargs=pkwargs,
                        title='Plot flux as funktion of exposure time', widgets=widgets,
                        widgetprops=wprops)
    
    gui3.master.mainloop()
    
if __name__ == "__main__":
    logger('Info: Starting routine to create a bad pixel map')
    log_params(params)
    
    if True:        # run everything with and without a badpx mask
        calimages['badpx_mask'] = np.ones((params['subframe'][0],params['subframe'][1]))
    
    biases, im_stats = [], []
    params['calibs'] = ['subframe', 'badpx_mask']
    for i in range_bias:
        biases.append(params['raw_data_path']+'Bias-{0}.fit'.format('%4.4i'%i))        # Change the name of the Bias files, if necessary
        if 'bias' in show_stats:
            im, im_head = read_file_calibration(params, biases[-1])
            im_stats.append(im)
    params['bias_rawfiles'] = biases
    params['bias_calibs_create'] = ['subframe', 'badpx_mask']
    params['master_bias_filename'] = 'master_bias.fits'
    im_bias, im_head_bias = create_image_general(params, 'bias')
    find_stats(im_stats)

    if len(get_statistics) > 0:
        #params['subframe'] = [48,48,1000,1000]      # measure the gain in a subwindow
        params['calibs'] = ['subframe', 'badpx_mask']
        images = []
        for image in get_statistics:
            images.append(params['raw_data_path']+image)
            im, im_head = read_file_calibration(params, images[-1])
            im_stats.append(im)
        find_stats(im_stats)
        

    im_flats = []
    for exptime in params['exptimes']:
        expname = str(exptime)
        expname = expname.replace('.','p')
        flats, darks, im_stats = [], [], []
        for i in range_dark:
            darks.append(params['raw_data_path']+'Dark-{0}_{1}s.fit'.format('%4.4i'%i, expname))        # Change the name of the Dark files, if necessary
            if 'dark' in show_stats:
                params['calibs'] = ['subframe', 'badpx_mask']
                im, im_head = read_file_calibration(params, darks[-1])
                im_stats.append(im)
        params['dark{0}_rawfiles'.format(exptime)] = darks
        params['master_dark{0}_filename'.format(exptime)] = 'master_dark_{0}s.fits'.format(expname)        # Change the name of the Flat files, if necessary
        #im_dark, im_head_dark = create_image_general(params, 'dark{0}'.format(exptime))                # if you need the master dark files
        find_stats(im_stats)
        
        if len(range_flat) == 0:    
            continue
        for i in range_flat:
            flats.append(params['raw_data_path']+'Flat-{0}_{1}s.fit'.format('%4.4i'%i, expname))
            if 'flat' in show_stats:
                params['calibs'] = ['subframe', 'badpx_mask']                                   # is overwritten, if a dark is loaded
                im, im_head = read_file_calibration(params, flats[-1])                          # disable, if only darks should be checked and flats with this exposure time don't exist
                im_stats.append(im)
        params['flatexp_rawfiles'] = flats                                                                  # don't use flat_rawfiles, as this will overwrite the standard flat
        params['flatexp_calibs_create'] = ['subframe', 'badpx_mask', 'normalise']
        params['master_flatexp_filename'] = 'master_flat_{0}s.fits'.format(expname)
        im_flat, im_head_flat = create_image_general(params, 'flatexp')                         # disable this and the following 7 lines, if only darks should be checked and flats with this exposure time don't exist
        ims = im_flat.shape
        ims = np.insert(ims, 0, 1)      #Add one dimension, to append the files
        im_flat.shape = ims
        if len(im_flats) == 0:
            im_flats = im_flat
        else:
            im_flats = np.append(im_flats, im_flat, axis=0)
        find_stats(im_stats)
        
    #exit(100)
    
    exptimes = np.array(params['exptimes'],dtype=float)
    badpx_mask = calimages['badpx_mask']
    try:
        im_head_bias
    except:
        im_head_bias = im_head
    
    # Check the gain and zeropoint for different max_good_values to check the linearity
    print('Check for linearity')
    # Create the data as Mugrauer 2010
    datac = np.median(im_flats, axis=(1,2))
    
    # Plot the data
    plot_linearity_UI(exptimes, im_flats)
    for max_value in [1200,2500,10000,13000,20000,30000,40000,50000,55000,60000,61000,62000,630000,640000]:
        gains = copy.copy(badpx_mask)*0
        zerop = copy.copy(badpx_mask)*0
        for i in tqdm(range(params['subframe'][0]), desc='determine the gain for each pixel for up to {0} ADU: maximum flux, gain, zeropoint'.format(max_value)):
            for j in range(params['subframe'][1]):
                exp_range = (im_flats[:,i,j] < max_value)# & (im_flats[:,i,j] > 100)
                fit = np.polyfit(exptimes[exp_range], im_flats[exp_range,i,j],1)
                gains[i,j] = fit[0]
                zerop[i,j] = fit[1]
        print(max_value, np.median(gains), np.median(zerop))
    
    
    if os.path.isfile(params['result_path']+params['gains']) == True and os.path.isfile(params['result_path']+params['zerop']) == True:
        params['calibs'] = ['subframe', 'badpx_mask']
        #gains, im_head = read_file(params, calimages, params['result_path']+params['gains'])
        #zerop, im_head = read_file(params, calimages, params['result_path']+params['zerop'])
        gains, im_head = read_file_calibration(params, params['result_path']+params['gains'])
        zerop, im_head = read_file_calibration(params, params['result_path']+params['zerop'])
    else:
        gains = copy.copy(badpx_mask)*0
        zerop = copy.copy(badpx_mask)*0
        for i in tqdm(range(params['subframe'][0]), desc='determine the gain for each pixel'):
            for j in range(params['subframe'][1]):
                exp_range = (im_flats[:,i,j] < params['max_good_value'])# & (im_flats[:,i,j] > 100)
                fit = np.polyfit(exptimes[exp_range], im_flats[exp_range,i,j],1)            # linear fit
                gains[i,j] = fit[0]
                zerop[i,j] = fit[1]
        save_im_fits(params, gains, im_head_bias, params['result_path']+params['gains'])
        save_im_fits(params, zerop, im_head_bias, params['result_path']+params['zerop'])

    gain_fit = fit_2d_image(gains, params['xord'], params['yord'])
    gain_diff = gains - gain_fit
    plot_img_spec.plot_image(gain_diff, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'residuals between fit of the gain and the gain')
    plot_img_spec.plot_image(gains, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'gain')
    plot_img_spec.plot_image(gain_fit, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'fit of the gain')
    plot_img_spec.plot_image(gain_diff, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], 'residuals between fit of the gain and the gain')
    
    """ # old solution without fit
    gain80pctl = percentile_list(sum(gains.tolist(),[]),0.1)
    gain, gain_std = np.mean(gain80pctl), np.std(gain80pctl)
    logger('Info: The gain is {0} +- {1} ADU/s'.format(round(gain,2), round(gain_std,2)))
    for sigm in [1,2,3,3.5,4,4.5,5,5.5,6,7,8]:
        print 'With {0} Sigma, {1} pixel would be marked as bad'.format(sigm, np.sum(abs(gains - gain) > sigm * gain_std))
    """
    gain95pctl = percentile_list(gain_diff.flatten(), 0.025)
    average, gain_std = np.mean(gain95pctl), np.std(gain95pctl)
    logger('Info: The average offset between fit and data is {0} ADU/s. The noise of the sensitivity is {1} ADU/s'.format(round(average,2), round(gain_std,2)))
    sigm, badpx = 1., 1
    while badpx > 0:
        badpx = np.sum(abs(gain_diff - average) > sigm * gain_std)
        print('With {0} Sigma, {1} pixel would be marked as bad'.format(sigm, badpx))
        sigm +=.5
    
    sigma = float(raw_input('What Sigma to use?\n>> '))
    
    badpx_mask[abs(gain_diff - average) > sigma * gain_std] = 0
    
    for fname in ['investigation_badpx.cvs', 'investigation_goodpx.cvs']:
        exptimes = params['exptimes']
        text = ['x+1','y+1','gain','zerop']
        for i in exptimes:
            text.append(str(i))
        coords = []
        if fname == 'investigation_badpx.cvs':
            coord = np.where(badpx_mask == 0)
            for i in range(coord[0].shape[0]):
                coords.append([coord[0][i],coord[1][i]])
        else:
            while len(coords)<50:
                i,j = random.randint(0,params['subframe'][0]), random.randint(0,params['subframe'][1])
                if badpx_mask[i,j] != 0:
                    coords.append([i,j])
        for [i,j] in coords:
            text[0] += '\t%1.1i'%(i+1)
            text[1] += '\t%1.1i'%(j+1)
            text[2] += '\t%1.1f'%gains[i,j]
            text[3] += '\t%1.1f'%zerop[i,j]
            for k in range(len(exptimes)):
                text[k+4] += '\t'
                if im_flats[k,i,j] < params['max_good_value']:
                    text[k+4] += '%1.1f'%im_flats[k,i,j]
                
        file = open(fname,'w')
        for line in text:
            file.write(line+'\n')
        file.close()
    
    save_im_fits(params, badpx_mask, im_head_bias, params['badpx_mask_filename'].rsplit('/',1)[1])
    log_params(params)
    logger('Info: Finished creating the bad pixel mask')
    
   
   
