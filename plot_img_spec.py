#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
from procedures import *
from astropy.io import fits
#import matplotlib       # To avoid crashing when ssh into Narit using putty
#if not 'DISPLAY' in os.environ:
#    matplotlib.use('agg')   # To avoid crashing when ssh into Narit using putty, however this means plots are not shown (test when working in front of uhppc30)
import matplotlib.pyplot as plt
import matplotlib.cm as cmaps
import matplotlib.colors as colors
import sys
from astropy.table import Table
from gatspy.periodic import LombScargleFast
import tkcanvas as tkc
if sys.version_info[0] < 3:
    import Tkinter as Tk
    from collections import OrderedDict as dict
else:
    import tkinter as Tk
import copy
import pickle

def plot_image(image, spaths, pctile=0, show=False, adjust=[0.05,0.95,0.95,0.05], title='', return_frame=False, frame=None, autotranspose=True, colorbar=True, axis_name=['x [px]','y [px]','flux [ADU]'], size=32):
    """
    Plots the raw CCD image to files $spaths$
    :param image: numpy 2D array, containing the CCD image
    :param spath: path to save the plot to
    :param pctile: Percentile of the picture to use in the plot. 0 means min and max values are used
    :param show: if True then plot is shown instead of saved
    :param adjust: to avoid big white areas in the figure
    :param title: Titel for the plot
    :return:
    """
    ims = image.shape
    im_good = image[~np.isnan(image)]
    if frame is None:
        fig, frame = plt.subplots(1, 1)
        if type(size).__name__ not in ('list', 'ndarray'):
            size = [size]
        fig.set_size_inches(size[0], size[-1])
    if len(ims) == 2:
        if ims[0] > ims[1] and show == True and autotranspose == True:		#transpose the image
            image = np.transpose(image)
            title = title + ' - transposed'
    plt.subplots_adjust(left=adjust[0], right=adjust[1], top=adjust[2], bottom=adjust[3])
    if type(pctile).__name__ == 'list':
        vmin = pctile[0]
        vmax = pctile[1]
    else:
        vmin = np.percentile(im_good,pctile)
        vmax = np.percentile(im_good,100-pctile)
    cframe = frame.imshow(image, cmap='gray',vmin=vmin, vmax=vmax)
    frame.set_xlabel(axis_name[0], fontsize=14)
    frame.set_ylabel(axis_name[1], fontsize=14)
    frame.set_title(title, fontsize=16)
    if colorbar:
        cbar = plt.colorbar(cframe,fraction=0.024, pad=0.02)
        cbar.set_label(axis_name[2], fontsize=14)
        cbar.ax.invert_yaxis()
    if return_frame:
        return frame
    if show:
        print("\n Please close the graph (Figure 1) to continue (Alt F4).\n")
        plt.show()
    else:
        if type(spaths).__name__ not in ('list', 'ndarray'):
            spaths = [spaths]
        for spath in spaths:
            plt.savefig(spath, bbox_inches='tight')
        plt.close()

def plot_spectra(spectra_x, spectra_y, labels, spaths, show=False, adjust=[0.05,0.95,0.95,0.05, 1.0,1.01], title='', return_frame=False, frame=None, size=[16.2, 10]):
    """
    Plots the Spectra to files $spaths$
    :param spectra: ???
    :param spath: path to save the plot to
    :param show: if True then plot is shown instead of saved
    :param adjust: to avoid big white areas in the figure: left, right, top, bottom, legend left, legend top
    :param title: Titel for the plot
    :return:
    """
    #spectra_x = np.array(spectra_x)
    if frame is None:
        fig, frame = plt.subplots(1, 1)
        if type(size).__name__ not in ('list', 'ndarray'):
            size = [size]
        fig.set_size_inches(size[0], size[-1])
    plt.subplots_adjust(left=adjust[0], right=adjust[1], top=adjust[2], bottom=adjust[3])
    dx = (np.max(spectra_x) - np.min(spectra_x))*0.01
    dy = (np.max(spectra_y) - np.min(spectra_y))*0.01
    plt.axis([np.min(spectra_x)-dx,np.max(spectra_x)+dx, np.min(spectra_y)-dy,np.max(spectra_y)+dy])
    #frame.plot(spectra_x,spectra_y)#, label=labels)
    if len(spectra_x.shape) == 1:
        frame.plot(spectra_x, spectra_y, label=labels)
    else:
        for i in range(len(spectra_x)):
            frame.plot(spectra_x[i], spectra_y[i], label=labels[i])
    frame.set_xlabel('x [px]', fontsize=14)
    frame.set_ylabel('flux [ADU]', fontsize=14)
    frame.set_title(title, fontsize=16)
    frame.legend(loc='upper left', bbox_to_anchor=(adjust[4], adjust[5]))
    if return_frame:
        return frame
    elif show:
        print("\n Please close the graph (Figure 1) to continue (Alt F4).\n")
        plt.show()
    else:
        if type(spaths).__name__ not in ('list', 'ndarray'):
            spaths = [spaths]
        for spath in spaths:
            plt.savefig(spath, bbox_inches='tight')
        plt.close()

def create_plot_linestyle_marker_markersize_color(length, linestyle=None, marker=None, markersize=None, color=None):
    """
    :param length: integer, number of individual plots in the graph
    :param linestyle, marker, color: List, Array, String, or None: use the given properties for the plot. String will be transformed to a list, and the list will multiplied by the number of entries
                If None is given for color then the matplotlib standards are used
    :return cycler_plt: cycler.Cycler, contains the information to make a easy looking graph
    """
    cycler_ori = plt.rcParams["axes.prop_cycle"]            # Get the current properties of the cycler
    #if linestyle == None:                       # replace None by empty string
    #    linestyle = ['']
    #if marker == None:                          # replace None by empty string
    #    marker = ['']
    if type(linestyle).__name__ == 'str':       # Make string to list
        linestyle = [linestyle]
    if type(marker).__name__ == 'str':          # Make string to list
        marker = [marker]
    if type(color).__name__ == 'str':           # Make string to list
        color = [color]
    if type(markersize).__name__ == 'int' or type(markersize).__name__ == 'float':           # Make number to list
        markersize = [markersize]
    if linestyle != None:
        for i in range(len(linestyle)):         # replace None by empty string
            if linestyle[i] == None:
                linestyle[i] = ''
    if marker != None:
        for i in range(len(marker)):            # replace None by empty string
            if marker[i] == None:
                marker[i] = ''
    if color != None:                           # user color
        cycler_plt = plt.cycler(color=color*length)
    else:                                       # No color given
        if 'color' in cycler_ori.by_key():      # user original colors
            cycler_plt = plt.cycler(color=cycler_ori.by_key()['color']*length)
        else:                                   # use dummy colors
            cycler_plt = plt.cycler(color=[u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']*length)
    if type(marker) == list or type(marker) == np.ndarray:
        cycler_plt += plt.cycler(marker=marker*length)
    elif 'marker' in cycler_ori.by_key():
        cycler_plt += plt.cycler(marker=cycler_ori.by_key()['marker']*length)
    if type(linestyle) == list or type(linestyle) == np.ndarray:
        cycler_plt += plt.cycler(linestyle=linestyle*length)
    elif 'linestyle' in cycler_ori.by_key():
        cycler_plt += plt.cycler(linestyle=cycler_ori.by_key()['linestyle']*length)
    if type(markersize) == list or type(markersize) == np.ndarray:
        cycler_plt += plt.cycler(markersize=markersize*length)
    
    return cycler_plt

def plot_points(data_x, data_y, labels, spaths, show=False, adjust=[0.05,0.95,0.95,0.05, 1.0,1.01], title='', return_frame=False, frame=None, x_title='x', y_title='y', linestyle="", marker="o", markersize=None, color=None, size=[16.2, 10]):
    """
    Plots the Spectra to files $spaths$
    :param spectra: ???
    :param spath: path to save the plot to
    :param show: if True then plot is shown instead of saved
    :param adjust: to avoid big white areas in the figure: left, right, top, bottom, legend left, legend top
    :param title: Titel for the plot
    :param size: size of the plot. For display on screen: 16.2, 10 ; for publication (single column): 6.7, 4.1; for publication, 2 columns: 3.3, 2.1
    :return:
    """
    if frame is None:
        fig, frame = plt.subplots(1, 1)
        if type(size).__name__ not in ('list', 'ndarray'):
            size = [size]
        fig.set_size_inches(size[0], size[-1])
    plt.subplots_adjust(left=adjust[0], right=adjust[1], top=adjust[2], bottom=adjust[3])
    minx, miny, maxx, maxy = 1e10, 1e10, -1e10, -1e10
    for i in range(len(data_x)):
        if type(data_x[i]) in [np.ndarray, list]:
            if sum(~np.isnan(list(data_x[i]))) > 0:
                minx = min(minx, np.nanmin(data_x[i]) )
                maxx = max(maxx, np.nanmax(data_x[i]) )
            data_yi = np.array(data_y[i], dtype=float)           # Necessary, because data_y[i][ ~np.isnan(data_y[i]) ] fails, float necessary, otherwise data_yi.dtype will object, which fails np.isnan
            #print type(data_yi), data_yi.shape, data_yi.dtype, np.isnan(data_yi)
            if sum(~np.isnan(data_yi)) > 0:       # only search min/max, if non-nan data exists
                miny = min(miny, np.nanmin(data_yi) )
                maxy = max(maxy, np.nanmax(data_yi) )
        else:
            minx = min(minx, data_x[i] )
            maxx = max(maxx, data_x[i] )
            miny = min(miny, data_y[i] )
            maxy = max(maxy, data_y[i] )
    #miny,maxy = -0.3,0.3
    #minx,maxx = 0,4250
    dx = max(1,maxx - minx)*0.01
    dy = max(1,maxy - miny)*0.01
    plt.axis([minx-dx,maxx+dx, miny-dy,maxy+dy])

    cycler_plt = create_plot_linestyle_marker_markersize_color(len(data_x), linestyle, marker, markersize, color)
    frame.set_prop_cycle(cycler_plt)
    label = ''
    for i in range(len(data_x)):
        if len(labels) > 0:
            label = labels[i]
            if label == '':
                label = None
        if type(data_x[i]) is not np.ndarray and type(data_x[i]) is not list:
            frame.plot(data_x[i], data_y[i], label=label)#, linestyle=linestyle[i], marker=marker[i])
        elif len(data_x[i]) - len(data_y[i]) == 1:        # Histogram (one data point less in y
            width = (( 1 - i/len(data_x) ) * 0.7 + 0.1) * (data_x[i][1] - data_x[i][0])
            center = (data_x[i][:-1] + data_x[i][1:]) / 2
            plt.bar(center, data_y[i], align='center', width=width, label=label)
        elif len(data_x[i]) > 0:
            frame.plot(data_x[i], data_y[i], label=label)#, linestyle=linestyle[i], marker=marker[i])
            
    frame.set_xlabel(x_title, fontsize=13)
    frame.set_ylabel(y_title, fontsize=13)
    frame.set_title(title, fontsize=15)
    frame.legend(loc='upper left', bbox_to_anchor=(adjust[4], adjust[5]))
    if return_frame:
        return frame
    elif show:
        print("\n Please close the graph (Figure 1) to continue (Alt F4).\n")
        plt.show()
    else:
        if type(spaths).__name__ not in ('list', 'ndarray'):
            spaths = [spaths]
        for spath in spaths:
            plt.savefig(spath, bbox_inches='tight')
        plt.close()

def save_obj(obj, name ):
    try:
        with open(name + '.pkl', 'wb') as f:
            pickle.dump(obj, f, 0)
    except:
        logger('Warn: Cannot save {0}.'.format(name + '.pkl'))

def load_obj(name ):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

def plot_spectra_UI(im, title=''):

    ims = im.shape
    iml = len(ims)

    settings = dict()
    if os.path.isfile('plot_settings.pkl') == True:
        settings = load_obj('plot_settings')
    
    fig, frame = plt.subplots(1, 1, figsize=(8,6))
    #plt.subplots_adjust(left=adjust[0], right=adjust[1], top=adjust[2], bottom=adjust[3])

    def plot(**param):
        frame, im = param['frame'], param['im']
        x_range, y_range, old_xlabel_text = copy.copy(frame.get_xlim()), copy.copy(frame.get_ylim()), copy.copy(frame.get_xlabel())
        frame.clear()
        size_inch = copy.deepcopy(fig.get_size_inches())
        adjust=[0.04,0.97,0.97,0.03, 1.0,1.01]
        adjust[0] = max(0.04, 0.2 - 0.04682*(size_inch[0] - 3.8)** 0.4416)          # left
        adjust[3] = max(0.03, 0.0 + 0.47866*(size_inch[1] - 0.0)**-0.9368)          # bottom
        if param['draw_filenames']:
            frame.set_title(title, fontsize=14)
            adjust[2] = 0.75
        if param['draw_legend']:
            adjust[1] = min(0.95, 0.7 + 0.07518*(size_inch[0] - 3.8)** 0.4289)          # right
        plt.subplots_adjust(left=adjust[0], right=adjust[1], top=adjust[2], bottom=adjust[3])
        
        param['name_data_-1'] = ['data']
        param['data'] = im
        for j in range(0,iml-1):
            param['name_data_{0}'.format(j)]=[]
            for data_entry in param['name_data_{0}'.format(j-1)]:
                data = param[data_entry]
                for i in range(ims[j]):
                    if ( param['exclude_{0}'.format(j)] and i not in param['plot_sub_{0}'.format(j)] ) or (not param['exclude_{0}'.format(j)] and i in param['plot_sub_{0}'.format(j)]):
                        if data_entry+'_{0}'.format(i) not in param.keys():     # Only load if necessary
                            param[data_entry+'_{0}'.format(i)] = data[i]
                        param['name_data_{0}'.format(j)].append(data_entry+'_{0}'.format(i))
        #print 'param.keys()',param.keys()
        #print "param['name_data_{0}'.format(iml-1-1)]",param['name_data_{0}'.format(iml-1-1)]
        settings = dict()
        for entry in param.keys():
            if entry.find('plot_sub_') == 0 or entry.find('exclude_') == 0:
                settings[entry] = param[entry]
        save_obj(settings, 'plot_settings')
        plot_ranges = [1E10,0,1E10,0]
        xlabel_text = 'Nothing to plot -> modify parameters or remove plot_settings.pkl'
        if 'name_data_{0}'.format(iml-1-1) in param.keys():
            for data_entry in tqdm(param['name_data_{0}'.format(iml-1-1)]):
                data_label = data_entry.replace('data_','')
                data = np.array(param[data_entry], float)
                if param['wavelength_plot'].find('w') != -1:
                    w_pos = 0
                    searchpos = param['wavelength_plot'].find('w')+1
                    for searchlength in range(len(param['wavelength_plot'])-searchpos, 0,-1):
                        try:
                            w_pos = int(param['wavelength_plot'][searchpos:searchpos+searchlength+1])
                            break
                        except:
                            w_pos = 0
                    data_entry_w = data_entry.rsplit('_',iml-2)                     # iml-2 should be 2 for nomal images
                    if iml == 4:                                      # for normal images
                        #data_entry_w = data_entry_w[0]+'_0_'+data_entry_w[2]        # Select the wavelength in the data file
                        data_entry_w = '{0}_{1}_{2}'.format(data_entry_w[0],w_pos,data_entry_w[2])        # Select the wavelength in the data file
                    elif iml ==3:                                                           # for images with linearised wavelength
                        #data_entry_w = data_entry_w[0]+'_0'     # barycentric
                        data_entry_w = '{0}_{1}'.format(data_entry_w[0],w_pos)    # uncorrected wavelength
                    else:
                        print('not implemented')
                    if data_entry_w not in param.keys():
                        text = data_entry_w.split('_')
                        if iml == 4:
                            param[data_entry_w] = im[int(text[1]), int(text[2]), int(text[3]) ]
                        elif iml ==3:
                            param[data_entry_w] = im[int(text[1]), int(text[2]) ]
                    x_axis = np.array(param[data_entry_w], float)
                    xlabel_text = 'wavelength [Angstrom] (dataset {0})'.format(w_pos)
                else:
                    x_axis = np.arange(len(data), dtype=float)
                    xlabel_text = 'x [px]'
                x_axis = x_axis[~np.isnan(data)]
                data = data[~np.isnan(data)]
                """ Create a test signal
                A = 27. # amplitude
                period = 31 # in x_axis uits
                omega = 2. * np.pi / period # radians per second
                phi = 0.5 * np.pi # radians
                # Create the full signal
                data = A * np.sin(omega * x_axis + phi)"""

                if param['wavelength_plot'] == 'f':    # make a fft
                    x_axis = np.fft.fftfreq(len(data))
                    data = np.abs(np.fft.fft(data))         # px -> 1/px (dT -> f = 1/dT)
                    good_values = list(range(1,int(len(data)/2)))      # 0 contains the sum of the data, and n/2+1 the negative frequencies
                    data = data * 2 / data.shape[0]
                    data, x_axis = data[good_values], 1/x_axis[good_values]
                    xlabel_text = 'period [px]'
                elif param['wavelength_plot'].find('l') != -1:       # http://joseph-long.com/writing/recovering-signals-from-unevenly-sampled-data/
                    nout = 4*len(data) # number of frequency-space points at which to calculate the signal strength (output)
                    # the posible Periods in px scale: 2 to 10% length of the data; in wave-scale: diff between 2 closest to 10% diff between furthest away:
                    period_range = np.array([ 2*np.nanmin(np.abs(x_axis[1:] - x_axis[:-1])), 0.5*(np.nanmax(x_axis)-np.nanmin(x_axis)) ])
                    freq_range = 1.0 / period_range
                    """model = LombScargleFast().fit(x_axis, data-np.nanmedian(data))
                    df = (max(freq_range) - min(freq_range)) / nout
                    freqs = min(freq_range) + df * np.arange(nout)
                    periods = 1/freqs
                    power = model.score_frequency_grid(min(freq_range), df, nout)
                    #periods, power = model.periodogram_auto(nyquist_factor=0.5)
                    data = power"""
                    freqs = np.linspace(min(freq_range), max(freq_range), nout)
                    periods = 1.0 / freqs
                    angular_freqs = 2 * np.pi * freqs
                    #print len(x_axis), len(data), angular_freqs, data.dtype, x_axis.dtype
                    pgram = scipy.signal.lombscargle(x_axis, data, angular_freqs)
                    normalized_pgram = np.sqrt(4 * (pgram / data.shape[0]))
                    data = normalized_pgram
                    x_axis = periods
                    if param['wavelength_plot'].find('w') != -1:
                        xlabel_text = 'period [Angstrom]'
                    else:
                        xlabel_text = 'period [px]'
                if len(x_axis) > 0:
                    if param['draw_legend']:
                        frame.plot(x_axis, data, label=data_label)
                    else:
                        frame.plot(x_axis, data)
                    #print sum(data[70:3000]-min(data[70:3000])), sum(data[1550:3000]-min(data[1550:3000]))
                    plot_ranges[0] = min(plot_ranges[0], min(x_axis))
                    plot_ranges[1] = max(plot_ranges[1], max(x_axis))
                    plot_ranges[2] = min(plot_ranges[2], min(data))
                    plot_ranges[3] = max(plot_ranges[3], max(data))
        
        if x_range != (0.0, 1.0) and y_range != (0.0, 1.0) and old_xlabel_text==xlabel_text and param['reset_plot']==False:
            plt.axis([x_range[0], x_range[1], y_range[0], y_range[1]])
        else:
            dx = (plot_ranges[1] - plot_ranges[0])*0.01
            dy = (plot_ranges[3] - plot_ranges[2])*0.01
            plt.axis([plot_ranges[0]-dx,plot_ranges[1]+dx, plot_ranges[2]-dy, plot_ranges[3]+dy])
            param['reset_plot']=True
        
        frame.set_xlabel(xlabel_text, fontsize=14)
        frame.set_ylabel('flux [ADU]', fontsize=14)
        if param['draw_legend']:
            frame.legend(loc='upper left', bbox_to_anchor=(adjust[4], adjust[5]))
            
    # get kwargs
    pkwargs = dict()
    pkwargs['frame'] = frame
    pkwargs['im'] = im
    
    # define valid_function
    # input is one variable (the string input)
    # return is either:
    #   True and values
    # or
    #   False and error message
    def vfunc(xs):
        try:
            xwhole = xs.replace(',', ' ')
            xwhole = xwhole.replace(': ', ':')
            new_xs = xwhole.split()
            xs = []
            for entry in new_xs:
                if entry.find(':') > 0:
                    entry = entry.split(':')
                    for j in range(int(entry[0]),int(entry[1])+1):
                        xs.append(j)
                else:
                    xs.append(int(entry))
            return True, xs
        except:
            return False, ('Error, input must consist of integers \n '
                           'separated by commas or white spaces')
    def vfunc_bool(xs):
        if xs.lower() in ['true', 'yes']:
            xs = True
        else:
            xs = False
        return True, xs
    # define widgets
    widgets = dict()
    widgets['wavelength_plot'] = dict(label='Wavelength, FFT,\nLomb-Scargle',
                                comment='w/f/l/wl' ,
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=str, start='', valid_function=None,
                                width=5)
    pkwargs['wavelength_plot'] = ''
    for i in range(iml-1):
        starta, startb = [], True
        text = [ 'Which data?', 'comma or colon\nseparated list' ]
        if (iml-1) - i - 1 == 2:
            text = [ 'Which file', 'comma or colon\nseparated list']
        if (iml-1) - i - 1 == 1:
            text = [ 'Which data', 'list']
            starta, startb = [1], False
        if (iml-1) - i - 1 == 0:
            text = [ 'Which aperture', 'list']
            starta, startb = [int(ims[-2]/2)], False
        if 'plot_sub_{0}'.format(i) in settings.keys():
            starta = settings['plot_sub_{0}'.format(i)]
        if 'exclude_{0}'.format(i) in settings.keys():
            startb = settings['exclude_{0}'.format(i)]
        startat = str(starta).replace('[','').replace(']','')               # make the list into a text
        widgets['plot_sub_{0}'.format(i)] = dict(label=text[0], comment=None, #text[1],
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=str, start=startat, valid_function=vfunc,
                                width=10)
        widgets['exclude_{0}'.format(i)] = dict(label='Exclude',
                                kind='CheckBox', start=startb)
        pkwargs['plot_sub_{0}'.format(i)] = starta
        pkwargs['exclude_{0}'.format(i)] = startb
    widgets['reset_plot'] = dict(label='reset plot',
                                kind='CheckBox', start=False)
    pkwargs['reset_plot'] = False
    widgets['draw_legend'] = dict(label='draw legend',
                                kind='CheckBox', start=True)
    pkwargs['draw_legend'] = True
    widgets['draw_filenames'] = dict(label='draw filenames',
                                kind='CheckBox', start=True)
    pkwargs['draw_filenames'] = True
    #widgets['draw_size'] = dict(label='12.8 x 9.6 inch',
    #                            kind='CheckBox', start=False)
    #pkwargs['draw_size'] = False
    
    widgets['accept'] = dict(label='Close', kind='ExitButton', position=Tk.BOTTOM)
    widgets['update'] = dict(label='Update', kind='UpdatePlot', position=Tk.BOTTOM)
    
    # run initial update plot function
    plot(**pkwargs)
    
    wprops = dict(orientation='v', position=Tk.RIGHT)

    gui3 = tkc.TkCanvas(figure=fig, ax=frame, func=plot, kwargs=pkwargs,
                        title='Plot spectral orders', widgets=widgets,
                        widgetprops=wprops)
    
    gui3.master.mainloop()

# Start of code
if __name__ == "__main__":
    #read the data
    rpath = '/data/ronny/20170630_harps_test_data_red/'
    fitsfiles=[]
    #fitsfiles.append(['MasterBias_HARPS.fits',0,'i'])		#one CCD
    #fitsfiles.append(['MasterBias_HARPS.fits',1,'i'])		#other CCD
    #fitsfiles.append(['FlatOb_HARPS.fits',0,'i'])		#object fiber
    #fitsfiles.append(['FlatOb_HARPS.fits',1,'i'])
    #fitsfiles.append(['FlatCo_HARPS.fits',0,'i'])		#other fiber
    #fitsfiles.append(['FlatCo_HARPS.fits',1,'i'])
    #fitsfiles.append(['Flat_HARPS.fits',0,'i'])			#both fibers
    #fitsfiles.append(['Flat_HARPS.fits',1,'i'])
    #fitsfiles.append(['BACR_FLAT_HARPS.fits',0,'i'])			#???
    #fitsfiles.append(['BACB_FLAT_HARPS.fits',0,'i'])
    #fitsfiles.append(['P_ob_B_HARPS.fits',0,'i'])			#traces?
    #fitsfiles.append(['P_ob_R_HARPS.fits',0,'i'])
    #fitsfiles.append(['B_flat_ob_HARPS.fits',1,'s'])			#object fiber, B chip, flux
    #fitsfiles.append(['B_flat_ob_HARPS.fits',2,'s'])			#object fiber, B chip, error?
    #fitsfiles.append(['R_flat_ob_HARPS.fits',1,'s'])			#
    #fitsfiles.append(['R_flat_ob_HARPS.fits',2,'s'])			#
    #fitsfiles.append(['BACR_FLAT_CO_HARPS.fits',0,'i'])			#
    #fitsfiles.append(['BACB_FLAT_CO_HARPS.fits',0,'i'])			#
    #fitsfiles.append(['P_co_B_HARPS.fits',0,'i'])			#traces?
    #fitsfiles.append(['P_co_R_HARPS.fits',0,'i'])			#
    #fitsfiles.append(['B_flat_co_HARPS.fits',1,'s'])			#
    #fitsfiles.append(['B_flat_co_HARPS.fits',2,'s'])			#
    #fitsfiles.append(['R_flat_co_HARPS.fits',1,'s'])			#
    #fitsfiles.append(['R_flat_co_HARPS.fits',2,'s'])			#
    #fitsfiles.append(['BACR_HARPS.2014-04-01T19:30:52.924.fits',0,'i'])			#
    #fitsfiles.append(['BACB_HARPS.2014-04-01T19:30:52.924.fits',0,'i'])			#arc lines and something
    #fitsfiles.append(['HARPS.2014-04-01T19:30:52.924.spec.ob.B.fits.S',1,'s'])			#emission
    #fitsfiles.append(['HARPS.2014-04-01T19:30:52.924.spec.ob.B.fits.S',2,'s'])
    #fitsfiles.append(['HARPS.2014-04-01T19:30:52.924.spec.ob.R.fits.S',1,'s'])			#
    #fitsfiles.append(['HARPS.2014-04-01T19:30:52.924.spec.ob.R.fits.S',2,'s'])
    #fitsfiles.append(['HARPS.2014-04-01T19:30:52.924.spec.co.B.fits.S',1,'s'])			#emission
    #fitsfiles.append(['HARPS.2014-04-01T19:30:52.924.spec.co.B.fits.S',2,'s'])
    #fitsfiles.append(['HARPS.2014-04-01T19:30:52.924.spec.co.R.fits.S',1,'s'])			#
    #fitsfiles.append(['HARPS.2014-04-01T19:30:52.924.spec.co.R.fits.S',2,'s'])
    rpath = '/data/ronny/20170630_harps_test_data/'
    #fitsfiles.append(['HARPS.2014-04-01T19:30:52.924.fits',0,'i'])			#B, both fibers: arcs
    #rpath = '/data/ronny/Reduced/20170706_test/extracted/'
    #fitsfiles.append(['master_flat_Order_13.fits',0,'s'])			#
    #rpath = 'extracted/'
    #for i in range(1,37):
    #    fitsfiles.append(['SunArc-{0}_10s.npy'.format('%4.4i'%i),0,'sg'])        #20170912
    #fitsfiles.append(['Arc_long.npy',0,'sg'])        #20170913, 20170912
    #fitsfiles.append(['Arc_short.npy',0,'sg'])        #20170913
    #fitsfiles.append(['FlatArc-0001_10s.npy',0,'sg'])        #20170913
    #fitsfiles.append(['FlatArc-0001_4s.npy',0,'sg'])        #20170922
    #for i in range(1,10):
    #     fitsfiles.append(['Fiberbend-{0}_4s.npy'.format('%4.4i'%i),0,'sg'])        #20170922
    if len(sys.argv) > 1:
        filenames = []
        for filename in sys.argv[1:]:
            if filename.find('plot_') != -1:
                filenames.append(filename)
    else:
        filenames = ['plot_files.lst']
    for filename in filenames:
        rpath = ''
        files = read_text_file(filename, no_empty_lines=True, warn_missing_file=True)
        for line in files:
            if line[0] != '#':
                fitsfiles.append([line.replace(' ',''),0,'sg'])
    
    data_sg = []
    titel_sg = ''
    for fitsfile in fitsfiles:
        if fitsfile[0].find('.npy')>0:
            im=np.load(rpath+fitsfile[0])
        else:
            im = np.array(fits.getdata(rpath+fitsfile[0],0))
        ims = im.shape
        #print ims
        if fitsfile[2] == 'i':              # For images
            if len(ims) == 3:               # select the right image, if more than one is in the data
                im = im[:,:,fitsfile[1]]
            plot_image(im, ['savepaths'], 1, True, [0.05,0.95,0.95,0.05], str(fitsfile))
        elif fitsfile[2] == 's':
            if len(ims) == 3:
                imx = []
                imy = []
                labels = []
                for i in range(ims[0]):
                    imx.append(im[i,0,:])
                    imy.append(im[i,fitsfile[1],:])
                    labels.append('Order ' + str(i+1))
                #imx = np.transpose(np.array(imx))
                #imy = np.transpose(np.array(imy))
            elif len(ims) == 1:
                ims = Table.read(rpath+fitsfile[0])
                imx = np.array(ims['pixel number'])
                imy = np.array(ims['counts'])
                imx = imx[~np.isnan(imy)]
                imy = imy[~np.isnan(imy)]
                print(imx.shape,imy.shape,imx, imy)
                labels = ''
            else:
                imx = im[0,:]
                imy = im[fitsfile[1],:]
                labels = ''
            plot_spectra(imx,imy,labels, ['savepaths'], True, [0.06,0.92,0.95,0.06, 1.0,1.01], str(fitsfile))
        elif fitsfile[2] == 'sg':
            ims = np.insert(ims, 0, 1)          # Reshape to make the insert easier
            im.shape = ims
            if len(data_sg) == 0:
                data_sg = im
            else:
                datas = data_sg.shape
                while len(ims) < len(datas):        # Different number of data-subsets
                    im = im[..., np.newaxis].reshape(np.insert(ims, 2, 1))      # add a new axis, first at the end, and then move to second place by reshaping
                    ims = im.shape
                while len(datas) < len(ims):        # Different number of data-subsets
                    data_sg = data_sg[..., np.newaxis].reshape(np.insert(datas, 2, 1))      # add a new axis, first at the end, and then move to second place by reshaping
                    print('Added a new dimension at the third place')
                    datas = data_sg.shape
                #print 2,datas, ims
                for i in range(1, len(ims)):                # Diffent sizes of the array
                    if ims[i] < datas[i]:
                        ims = list(ims)
                        ims[i] = datas[i] - ims[i]
                        temp = np.empty(ims)
                        temp.fill(np.nan)
                        im = np.append(im, temp, axis=i)
                    elif datas[i] < ims[i]:
                        datas = list(datas)
                        datas[i] = ims[i] - datas[i]
                        temp = np.empty(datas)
                        temp.fill(np.nan)
                        data_sg = np.append(data_sg, temp, axis=i)
                    ims = im.shape
                    datas = data_sg.shape
                
                data_sg = np.append(data_sg, im, axis=0)
                #print data_sg.shape
            titel_sg += ',\n{0}'.format(fitsfile[0])
        #print data_sg.shape, data_sg.dtype, data_sg.itemsize, data_sg.nbytes, sys.getsizeof(data_sg), data_sg.nbytes   # size about the size of the fits file
    if data_sg != []:
        max_orders = 0
        for im in data_sg:
            max_orders = max(max_orders, im.shape[1])
        for i in range(len(data_sg)):
            nr_orders = data_sg[i].shape[1]
            if nr_orders < max_orders:
                data_sg[i] = np.insert(data_sg[i], list(range(nr_orders,max_orders)), np.zeros(data_sg[i].shape[2]), axis=1)
        print('Data read ({0} MB), plotting data now'.format(round(data_sg.nbytes/1024./1024,1)) )
        plot_spectra_UI(np.array(data_sg), title=titel_sg[2:])
            
#no of orders, 3 different values, px


















