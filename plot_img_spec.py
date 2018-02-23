#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.cm as cmaps
import matplotlib.colors as colors
import os
import sys
from astropy.table import Table
from procedures import *
import tkcanvas as tkc
if sys.version_info[0] < 3:
    import Tkinter as Tk
    from collections import OrderedDict as dict
else:
    import tkinter as Tk
import copy
import pickle

def plot_image(image, spaths, pctile=0, show=False, adjust=[0.05,0.95,0.95,0.05], title='', return_frame=False, frame=None, autotranspose=True, colorbar=True, axis_name=['x [px]','y [px]','flux [ADU]']):
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
        fig.set_size_inches(32, 32)
    if len(ims) == 2:
        if ims[0] > ims[1] and show == True and autotranspose == True:		#transpose the image
            image = np.transpose(image)
            title = title + ' - transposed'
    plt.subplots_adjust(left=adjust[0], right=adjust[1], top=adjust[2], bottom=adjust[3])
    cframe = frame.imshow(image, cmap='gray',vmin=np.percentile(im_good,pctile), vmax=np.percentile(im_good,100-pctile))
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
        for spath in spaths:
            plt.savefig(spath, bbox_inches='tight')
        plt.close()

def plot_spectra(spectra_x, spectra_y, labels, spaths, show=False, adjust=[0.05,0.95,0.95,0.05, 1.0,1.01], title='', return_frame=False, frame=None):
    """
    Plots the Spectra to files $spaths$
    :param spectra: ???
    :param spath: path to save the plot to
    :param show: if True then plot is shown instead of saved
    :param adjust: to avoid big white areas in the figure: left, right, top, bottom, legend left, legend top
    :param title: Titel for the plot
    :return:
    """
    if frame is None:
        fig, frame = plt.subplots(1, 1)
        fig.set_size_inches(16.2, 10)
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
        for spath in spaths:
            plt.savefig(spath, bbox_inches='tight')
        plt.close()

def plot_points(data_x, data_y, labels, spaths, show=False, adjust=[0.05,0.95,0.95,0.05, 1.0,1.01], title='', return_frame=False, frame=None, x_title='x', y_title='y', linestyle="", marker="o"):
    """
    Plots the Spectra to files $spaths$
    :param spectra: ???
    :param spath: path to save the plot to
    :param show: if True then plot is shown instead of saved
    :param adjust: to avoid big white areas in the figure: left, right, top, bottom, legend left, legend top
    :param title: Titel for the plot
    :return:
    """
    if frame is None:
        fig, frame = plt.subplots(1, 1)
        fig.set_size_inches(16.2, 10)
    plt.subplots_adjust(left=adjust[0], right=adjust[1], top=adjust[2], bottom=adjust[3])
    minx, miny, maxx, maxy = 1e10, 1e10, -1e10, -1e10
    for i in range(len(data_x)):
        if data_x[i] <> []:
            minx = min(minx, np.min(data_x[i]) )
            miny = min(miny, np.min(data_y[i]) )
            maxx = max(maxx, np.max(data_x[i]) )
            maxy = max(maxy, np.max(data_y[i]) )
    #miny,maxy = -0.3,0.3
    #minx,maxx = 0,4250
    dx = max(1,maxx - minx)*0.01
    dy = max(1,maxy - miny)*0.01
    plt.axis([minx-dx,maxx+dx, miny-dy,maxy+dy])
    label = ''
    for i in range(len(data_x)):
        if len(labels) > 0:
            label = labels[i]
        if len(data_x[i]) - len(data_y[i]) == 1:        # Histogram
            width = (( 1 - i/len(data_x) ) * 0.7 + 0.1) * (data_x[i][1] - data_x[i][0])
            center = (data_x[i][:-1] + data_x[i][1:]) / 2
            plt.bar(center, data_y[i], align='center', width=width, label=label)
        elif data_x[i] <> []:
            frame.plot(data_x[i], data_y[i], label=label, linestyle=linestyle, marker=marker)
            
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
        for spath in spaths:
            plt.savefig(spath, bbox_inches='tight')
        plt.close()

def save_obj(obj, name ):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, 0)

def load_obj(name ):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

def plot_spectra_UI(im, title='', adjust=[0.07,0.93,0.94,0.06, 1.0,1.01]):
    param=dict()
    
    ims = im.shape
    iml = len(ims)
    
    settings = dict()
    if os.path.isfile('plot_settings.pkl') == True:
        settings = load_obj('plot_settings')
    
    fig, frame = plt.subplots(1, 1)
    plt.subplots_adjust(left=adjust[0], right=adjust[1], top=adjust[2], bottom=adjust[3])
    
    def plot(frame, im, settings):
        x_range, y_range, old_xlabel_text = copy.copy(frame.get_xlim()), copy.copy(frame.get_ylim()), copy.copy(frame.get_xlabel())
        frame.clear()
        for i in range(iml-1):
            try:
                param['plot_sub_{0}'.format(i)] = gui3.data['plot_sub_{0}'.format(i)]
                param['exclude_{0}'.format(i)] = gui3.data['exclude_{0}'.format(i)]
            except:
                param['plot_sub_{0}'.format(i)] = []
                param['exclude_{0}'.format(i)] = True
                if (iml-1) - i - 1 == 1:
                    param['plot_sub_{0}'.format(i)] = [1]
                    param['exclude_{0}'.format(i)] = False
                if 'plot_sub_{0}'.format(i) in settings:
                    param['plot_sub_{0}'.format(i)] = settings['plot_sub_{0}'.format(i)]
                if 'exclude_{0}'.format(i) in settings:
                    param['exclude_{0}'.format(i)] = settings['exclude_{0}'.format(i)]
        try:
            wavelength_plot = gui3.data['wavelength_plot']
        except:
            wavelength_plot =''
        try:
            reset_plot = gui3.data['reset_plot']
            gui3.ws[-3].delete(0, 200)              #clear the text field
        except:
            reset_plot = ''
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
        if 'name_data_{0}'.format(iml-1-1) in param.keys():
            for data_entry in tqdm(param['name_data_{0}'.format(iml-1-1)]):
                data_label = data_entry.replace('data_','')
                data = param[data_entry]
                if wavelength_plot == 'w' or wavelength_plot == 'wl' or wavelength_plot == 'lw':
                    data_entry_w = data_entry.rsplit('_',2)
                    data_entry_w = data_entry_w[0]+'_0_'+data_entry_w[2]        # Select the wavelength in the data file
                    if data_entry_w not in param.keys():
                        text = data_entry_w.split('_')
                        param[data_entry_w] = im[int(text[1]), int(text[2]), int(text[3]) ]     # This won't work, if it's not 3 indexes
                    x_axis = param[data_entry_w]
                else:
                    x_axis = np.arange(len(data)).astype(float)
                x_axis = x_axis[~np.isnan(data)]
                data = data[~np.isnan(data)]
                """ Create a test signal
                A = 27. # amplitude
                period = 31 # in x_axis uits
                omega = 2. * np.pi / period # radians per second
                phi = 0.5 * np.pi # radians
                # Create the full signal
                data = A * np.sin(omega * x_axis + phi)"""

                if wavelength_plot == 'f':    # make a fft
                    x_axis = np.fft.fftfreq(len(data))
                    data = np.abs(np.fft.fft(data))         # px -> 1/px (dT -> f = 1/dT)
                    good_values = range(1,len(data)/2)      # 0 contains the sum of the data, and n/2+1 the negative frequencies
                    data = data * 2 / data.shape[0]
                    data, x_axis = data[good_values], 1/x_axis[good_values]
                elif wavelength_plot == 'l' or wavelength_plot == 'wl' or wavelength_plot == 'lw':       # http://joseph-long.com/writing/recovering-signals-from-unevenly-sampled-data/
                    nout = 5000 # number of frequency-space points at which to calculate the signal strength (output)
                    # the posible Periods in px scale: 2 to 10% length of the data; in wave-scale: diff between 2 closest to 10% diff between furthest away:
                    period_range = np.array([ 2*np.nanmin(np.abs(x_axis[1:] - x_axis[:-1])), 0.5*(np.nanmax(x_axis)-np.nanmin(x_axis)) ])
                    #periods = np.linspace(min(period_range), max(period_range), nout)
                    #freqs = 1.0 / periods
                    freqs = np.linspace(min(1/period_range), max(1/period_range), nout)
                    periods = 1.0 / freqs
                    angular_freqs = 2 * np.pi * freqs
                    pgram = scipy.signal.lombscargle(x_axis, data, angular_freqs)
                    normalized_pgram = np.sqrt(4 * (pgram / data.shape[0]))
                    x_axis = periods
                    data = normalized_pgram

                frame.plot(x_axis, data, label=data_label)
                #print sum(data[70:3000]-min(data[70:3000])), sum(data[1550:3000]-min(data[1550:3000]))
                plot_ranges[0] = min(plot_ranges[0], min(x_axis))
                plot_ranges[1] = max(plot_ranges[1], max(x_axis))
                plot_ranges[2] = min(plot_ranges[2], min(data))
                plot_ranges[3] = max(plot_ranges[3], max(data))
        if wavelength_plot == 'w':
            xlabel_text = 'wavelength [Angstrom]'
        elif wavelength_plot == 'f':
            xlabel_text = 'period [px]'
        elif wavelength_plot == 'l':
            xlabel_text = 'period [px]'
        elif wavelength_plot == 'wl' or wavelength_plot == 'lw':
            xlabel_text = 'period [Angstrom]'
        else:
            xlabel_text = 'x [px]'
        if x_range <> (0.0, 1.0) and y_range <> (0.0, 1.0) and old_xlabel_text==xlabel_text and reset_plot=='':
            plt.axis([x_range[0], x_range[1], y_range[0], y_range[1]])
        else:
            dx = (plot_ranges[1] - plot_ranges[0])*0.01
            dy = (plot_ranges[3] - plot_ranges[2])*0.01
            plt.axis([plot_ranges[0]-dx,plot_ranges[1]+dx, plot_ranges[2]-dy, plot_ranges[3]+dy])
        frame.set_xlabel(xlabel_text, fontsize=14)
        frame.set_ylabel('flux [ADU]', fontsize=14)
        frame.set_title(title, fontsize=16)
        frame.legend(loc='upper left', bbox_to_anchor=(adjust[4], adjust[5]))
    # get kwargs
    pkwargs = dict(frame=frame, im = im, settings=settings)
    # run initial update plot function
    plot(**pkwargs)
    
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
        if xs.lower() in ['true']:
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
    for i in range(iml-1):
        starta, startb = '', 'True'
        text = [ 'Which data?', 'comma or colon\nseparated list' ]
        if (iml-1) - i - 1 == 2:
            text = [ 'Which file', 'comma or colon\nseparated list']
        if (iml-1) - i - 1 == 1:
            text = [ 'Which data', 'list']
            starta, startb = '1', ''
        if (iml-1) - i - 1 == 0:
            text = [ 'Which order', 'list']
        
        if 'plot_sub_{0}'.format(i) in settings.keys():
            starta = (('{0}'.format(settings['plot_sub_{0}'.format(i)])).replace('[','')).replace(']','')
        if 'exclude_{0}'.format(i) in settings.keys():
            startb = '{0}'.format(settings['exclude_{0}'.format(i)])
        widgets['plot_sub_{0}'.format(i)] = dict(label=text[0], comment=None, #text[1],
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=str, start=starta, valid_function=vfunc,
                                width=10)
        widgets['exclude_{0}'.format(i)] = dict(label='Exclude',
                                comment=None, #'True for exclude' ,
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=str, start=startb, valid_function=vfunc_bool,
                                width=5)
    widgets['reset_plot'] = dict(label='reset plot',
                                comment=None, #'anything to reset' ,
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=str, start='', valid_function=None,
                                width=5)
    widgets['accept'] = dict(label='Close', kind='ExitButton', position=Tk.BOTTOM)
    widgets['update'] = dict(label='Update', kind='UpdatePlot', position=Tk.BOTTOM)

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
    if True:
        rpath = ''
        files = read_text_file('plot_files.lst')
        for line in files:
            if line[0] <> '#' and line.find('Stray') == -1:
                fitsfiles.append([line,0,'sg'])
    
    data_sg = []
    titel_sg = ''
    for fitsfile in fitsfiles:
        if fitsfile[0].find('.npy')>0:
            im=np.load(rpath+fitsfile[0])
        else:
            im = np.array(fits.getdata(rpath+fitsfile[0],0))
        ims = im.shape
        if fitsfile[2] == 'i':
            if len(ims) == 3:
                im = im[:,:,fitsfile[1]]
            #ims2 = im.shape
            #if len(ims) == 2:
            #    if ims2[0] > ims2[1]:		#transpose the image
            #        im = np.transpose(im)
        if fitsfile[2] == 'i':
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
                print imx.shape,imy.shape,imx, imy
                labels = ''
            else:
                imx = im[0,:]
                imy = im[fitsfile[1],:]
                labels = ''
            plot_spectra(imx,imy,labels, ['savepaths'], True, [0.06,0.92,0.95,0.06, 1.0,1.01], str(fitsfile))
        elif fitsfile[2] == 'sg':
            ims = np.insert(ims, 0, 1)
            im.shape = ims
            if len(data_sg) == 0:
                data_sg = im
            else:
                data_sg = np.append(data_sg, im, axis=0)
            titel_sg += ',\n{0}'.format(fitsfile[0])
        #print data_sg.shape, data_sg.dtype, data_sg.itemsize, data_sg.nbytes, sys.getsizeof(data_sg), data_sg.nbytes   # size about the size of the fits file
    if data_sg <> []:
        max_orders = 0
        for im in data_sg:
            max_orders = max(max_orders, im.shape[1])
        for i in range(len(data_sg)):
            nr_orders = data_sg[i].shape[1]
            if nr_orders < max_orders:
                data_sg[i] = np.insert(data_sg[i], range(nr_orders,max_orders),np.zeros(data_sg[i].shape[2]), axis=1)
        print 'Data read ({0} MB), plotting data now'.format(round(data_sg.itemsize/1024./1024))
        plot_spectra_UI(np.array(data_sg), title=titel_sg[2:])
            
#no of orders, 3 different values, px


















