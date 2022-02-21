import os
import sys
import copy
import datetime
import operator
import warnings
import numpy as np
import matplotlib.pyplot as plt
import tkcanvas as tkc
if sys.version_info[0] < 3:     # Python 2
    import Tkinter as Tk
    from collections import OrderedDict as dict
else:                           # Python 3
    import tkinter as Tk
import plot_img_spec
from proc_general import *


def file_list_UI(file_list, CONFIGFILE):
    from procedures import logger
    # Extract a common path
    maxlen_files = 0
    for entry in file_list:
        maxlen_files = max(maxlen_files, len(entry[0].replace('#','').replace(' ','')) )
    common_path = ''
    if file_list[0][0][2:].find(os.sep) != -1:         # it contains subfolders
        fnames = []
        for entry in file_list:
            fnames.append( entry[0].replace('#','').replace(' ','') )
        common_path = fnames[0].rsplit(os.sep,1)[0]    # just remove the filename
        for dummy in range(len(fnames[0].split(os.sep))-1):
            found_common_path = True
            for entry in fnames:
                if entry.find(common_path) == -1:   # This path is not the same
                    found_common_path = False
                    break
            if found_common_path:
                common_path += os.sep
                break
            else:
                common_path = common_path.rsplit(os.sep,1)     # remove the next subfolder
                if len(common_path) == 2:                   # Another subfolder in the path
                    common_path = common_path[0]
                else:                                       # no more subfolders
                    common_path = ''
                    break
    # Split the common_path to make the GUI shorter (in x), if necessary
    maxlen_files -= len(common_path)
    if len(common_path) == 0:
        common_path_text = '{0}{1}(no common path)'.format(' '*int(maxlen_files*1.5), os.linesep)
    else:
        common_path_text = common_path
    if len(common_path) > maxlen_files:
        common_path_temp = common_path[:-1].split(os.sep)          # without the last '/'
        common_path_text = [common_path_temp[0]+os.sep]            # start with the first entry
        ii = 0
        for entry in common_path_temp[1:]:
            entry += os.sep
            if len(common_path_text[ii]+entry) > maxlen_files and len(common_path_text[ii]) > maxlen_files/2.:     # If too long to add
                common_path_text.append(entry)
                ii += 1
            else:
                common_path_text[ii] += entry
        common_path_text = os.linesep.join(common_path_text)
    
    # define widgets
    pkwargs = dict()
    widgets = dict()
    widgets['comment'] = dict(label='Comm-{0}ented{0}out '.format(os.linesep), kind='Label', row=0, column=0, columnspan=1, orientation=Tk.W)
    widgets['mid_exp'] = dict(label='   Observation Time    {0}(mid exposure){0}UTC'.format(os.linesep), kind='Label', row=0, column=1, orientation=Tk.W)
    widgets['exp']     = dict(label='Expo-{0}sure {0}[s]  '.format(os.linesep), kind='Label', row=0, column=2, orientation=Tk.W)
    widgets['name']    = dict(label='Path and folder{1}{0}'.format(common_path_text, os.linesep), kind='Label', row=0, column=3, orientation=Tk.W)
    #widgets['fib1']    = dict(label='Science\nfiber', kind='Label', row=0, column=5)
    #widgets['fib2']    = dict(label='Calibration\nfiber', kind='Label', row=0, column=6)
    widgets['b']       = dict(label='Bias', kind='Label', row=0, column=6)
    widgets['d']       = dict(label='Dark', kind='Label', row=0, column=7)
    widgets['a']       = dict(label='Real{0}Flat'.format(os.linesep), kind='Label', row=0, column=8)
    widgets['t1']      = dict(label='Sci.{0}tra-{0}ce'.format(os.linesep), kind='Label', row=0, column=9)
    widgets['t2']      = dict(label='Cal.{0}tra-{0}ce'.format(os.linesep), kind='Label', row=0, column=10)
    widgets['z']       = dict(label='Blaze', kind='Label', row=0, column=11)
    widgets['w1']      = dict(label='Wave{0}Sci.'.format(os.linesep), kind='Label', row=0, column=12)
    #widgets['w1l']     = dict(label='Wave\nSci.\nlong', kind='Label', row=0, column=12)
    #widgets['w1s']     = dict(label='Wave\nSci.\nshort', kind='Label', row=0, column=13)
    widgets['w2']     = dict(label='Wave{0}Cal.'.format(os.linesep), kind='Label', row=0, column=14)
    #widgets['w2l']     = dict(label='Wave\nCal.\nlong', kind='Label', row=0, column=14)
    #widgets['w2s']     = dict(label='Wave\nCal.\nshort', kind='Label', row=0, column=15)
    widgets['ws']      = dict(label='Wave{0}shft{0}Sci'.format(os.linesep), kind='Label', row=0, column=16)
    widgets['wc']      = dict(label='Wave{0}shft{0}Cal'.format(os.linesep), kind='Label', row=0, column=17)
    widgets['e']       = dict(label='Ex-{0}tract'.format(os.linesep), kind='Label', row=0, column=18)
    widgets['obname']  = dict(label='Objectname{0}(comma{0}separated)'.format(os.linesep), kind='Label', row=0, column=19)
    widgets['extra']   = dict(label='Further usage{0}of files{0}(comma sep.)'.format(os.linesep), kind='Label', row=0, column=20)
    for ii, entry in enumerate(file_list):
        pkwargs['comment_{0}'.format(ii)] = ( 0 <= entry[0].find('#') < 20 )            # Commented out
        widgets['comment_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['comment_{0}'.format(ii)], row=ii+1, column=0)
        expstr = '%1.2f'%entry[3]       # '%9.2f'%entry[3] has not long enough spaces
        if len(expstr) > 7:
            expstr = '  '+expstr[:-3]        # without fractions at the end
        else:
            #expstr = ' '*{4:7, 5:5, 6:4, 7:2, 8:1, 9:1, 10:1}[len(expstr)] + expstr            # This has the 100s not perfectly alligned
            expstr = ' '*{4:6, 5:4, 6:2, 7:1, 8:1, 9:1, 10:1}[len(expstr)] + expstr             # This has the 100s and 100s not perfectly alligned
        text = '{0}{1}   {2}'.format( datetime.datetime.utcfromtimestamp(entry[4]).strftime('%Y-%m-%d %H:%M:%S'), 
                                     expstr, entry[0].replace('#','').replace(' ','').replace(common_path,'') )
        widgets['name_{0}'.format(ii)] = dict(label=text, kind='Label', row=ii+1, column=1, columnspan=3, orientation=Tk.W)
        #fname = entry[0].replace('#','').replace(' ','').replace(common_path,'')
        #widgets['mid_exp_{0}'.format(ii)] = dict(label=datetime.datetime.utcfromtimestamp(entry[4]).strftime('%Y-%m-%dT%H:%M:%S'), kind='Label', row=ii+1, column=5)
        #widgets['exp_{0}'.format(ii)]  = dict(label=entry[3], kind='Label', row=ii+1, column=4)
        #widgets['name_{0}'.format(ii)] = dict(label=fname, kind='Label', row=ii+1, column=1, columnspan=1, orientation=Tk.W)
        #pkwargs['fib1_{0}'.format(ii)] = entry[1]     # Allows modification of the fiber content
        #widgets['fib1_{0}'.format(ii)] = dict(kind='TextEntry', minval=None, maxval=None, fmt=str, start=pkwargs['fib1_{0}'.format(ii)], width=8, row=ii+1, column=5)     # Allows modification of the fiber content
        #pkwargs['fib2_{0}'.format(ii)] = entry[2]     # Allows modification of the fiber content
        #widgets['fib2_{0}'.format(ii)] = dict(kind='TextEntry', minval=None, maxval=None, fmt=str, start=pkwargs['fib2_{0}'.format(ii)], width=8, row=ii+1, column=6)     # Allows modification of the fiber content
        flags = entry[6].replace(' ','').split(',')
        pkwargs['b_{0}'.format(ii)] = ( 'b' in flags )
        widgets['b_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['b_{0}'.format(ii)], row=ii+1, column=6)
        pkwargs['d_{0}'.format(ii)] = ( 'd' in flags )
        widgets['d_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['d_{0}'.format(ii)], row=ii+1, column=7)
        pkwargs['a_{0}'.format(ii)] = ( 'a' in flags )
        widgets['a_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['a_{0}'.format(ii)], row=ii+1, column=8)
        pkwargs['t1_{0}'.format(ii)] = ( 't1' in flags )
        widgets['t1_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['t1_{0}'.format(ii)], row=ii+1, column=9)
        pkwargs['t2_{0}'.format(ii)] = ( 't2' in flags )
        widgets['t2_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['t2_{0}'.format(ii)], row=ii+1, column=10)
        pkwargs['z_{0}'.format(ii)] = ( 'z' in flags )
        widgets['z_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['z_{0}'.format(ii)], row=ii+1, column=11)
        pkwargs['w1_{0}'.format(ii)] = ( ('w1' in flags) | ('w1l' in flags) | ('w1s' in flags) )
        widgets['w1_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['w1_{0}'.format(ii)], row=ii+1, column=12)
        #pkwargs['w1l_{0}'.format(ii)] = ( 'w1l' in flags )
        #widgets['w1l_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['w1l_{0}'.format(ii)], row=ii+1, column=12)
        #pkwargs['w1s_{0}'.format(ii)] = ( 'w1s' in flags )
        #widgets['w1s_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['w1s_{0}'.format(ii)], row=ii+1, column=13)
        pkwargs['w2_{0}'.format(ii)] = ( ('w2' in flags) | ('w2l' in flags) | ('w2s' in flags) )
        widgets['w2_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['w2_{0}'.format(ii)], row=ii+1, column=14)
        #pkwargs['w2l_{0}'.format(ii)] = ( 'w2l' in flags )
        #widgets['w2l_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['w2l_{0}'.format(ii)], row=ii+1, column=14)
        #pkwargs['w2s_{0}'.format(ii)] = ( 'w2s' in flags )
        #widgets['w2s_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['w2s_{0}'.format(ii)], row=ii+1, column=15)
        pkwargs['ws_{0}'.format(ii)] = ( 'ws' in flags )
        widgets['ws_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['ws_{0}'.format(ii)], row=ii+1, column=16)
        pkwargs['wc_{0}'.format(ii)] = ( 'wc' in flags )
        widgets['wc_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['wc_{0}'.format(ii)], row=ii+1, column=17)
        pkwargs['e_{0}'.format(ii)] = ( 'e' in flags )
        widgets['e_{0}'.format(ii)] = dict(label=None,  kind='CheckBox', start=pkwargs['e_{0}'.format(ii)], row=ii+1, column=18)
        pkwargs['obname_{0}'.format(ii)] = entry[5]
        widgets['obname_{0}'.format(ii)] = dict(kind='TextEntry', minval=None, maxval=None, fmt=str, start=pkwargs['obname_{0}'.format(ii)], width=15, row=ii+1, column=19)
        extra = ''
        for flag in flags:          # Add extra flags, if necessary
            if flag not in ['b', 'd', 'a', 't1', 't2', 'z', 'w2l', 'w2s', 'w2', 'w1l', 'w1s', 'w1', 'ws', 'wc', 'e']:
                extra += ','+flag
        if len(extra) > 0:
            extra = extra[1:]
        pkwargs['extra_{0}'.format(ii)] = extra
        widgets['extra_{0}'.format(ii)] = dict(kind='TextEntry', minval=None, maxval=None, fmt=str, start=pkwargs['extra_{0}'.format(ii)], width=15, row=ii+1, column=20)
    
    explain = ('Explanation of the columns:{1}'+\
               '- Tick first column to not use some files at all{1}'+\
               '- Mark the files to be used for calibration:{1}'+\
               '-- Bias: These file are combined into a master bias{1}'+\
               '-- Dark: exposure time will be automatically taken{1}   into account{1}'+\
               '-- Real Flat: Evenly exposed detector to calibrate{1}   pixel-to-pixel sensitivity variation{1}'+\
               '-- Science trace: To trace the science orders{1}'+\
               '-- Calibration trace: To trace the calibration orders{1}'+\
               '-- Blaze: To derive the blaze function{1}'+\
               '-- Wavelength solition for science fiber {1}   (long and short expsoure time){1}'+\
               '-- Wavelength solution for calibration fiber (*){1}'+\
               '-- Wavelength offset Science fiber (**) to correct for{1}   wavelength drit{1}'+\
               '-- Wavelength offset between the Science fiber and{1}   Calibration fiber (*){1}'+\
               '-- Object name: List of Names for Object{1}'+\
               '-- Extract: Extract these files on an individual basis{1}'+\
               '-- Further settings (manual): e.g. to combine files{1}   before or after extraction{1}'+\
               '(*) not for single fiber spectrographs{1}'+\
               '(**) important for unstabilised (single) fiber{1}     spectrographs{1}{1}'+\
               'The automatic assignment is based on the parameters{1} raw_data_* in {0} (and in procedure{1} add_new_rawfiles_file_list). ').format(CONFIGFILE, os.linesep)
              #'- Type of Science and Calibration\n  fibers are derived from header or\n  filename and can be changed here\n  (optional)\n'+\     # Allows modification of the fiber content
    for ii, commentii in enumerate(explain.split(os.linesep)):
        if len(commentii) > 0:
            widgets['explain_{0}'.format(ii)] = dict(label=commentii, kind='Label', row=ii, column=21, rowspan=1, orientation=Tk.W )#, wraplength=100 )      
    widgets['accept'] = dict(label='Accept', kind='ExitButton', row=ii+1, column=21, rowspan=2)
                              
    wprops = dict(fullscreen=False )
    #wprops['width_data'] = 800   # not neccssary, as uses the automatic width
    
    if len(file_list) > 100:
        logger('Info: A GUI with {0} elements will be created, on some machines that can take up to a few minutes.{1}'.format(len(widgets), os.linesep)+\
               '\tIf you want to use an editor instead of the GUI, please kill the process (kill {0}){1}'.format(os.getpid(), os.linesep)+\
               '\tand run the script with parameter "nogui", e.g.{1}\t\tpython {0} nogui'.format(sys.argv[0], os.linesep) )

    gui3 = tkc.TkCanvasGrid(title='HiFLEx: Asigning Files to extraction steps (What file contains what data)', 
                            kwargs=pkwargs,widgets=widgets, widgetprops=wprops )
    gui3.master.mainloop()
    
    # Get the information from the GUI
    file_list_new, file_list_full = [], []
    for ii in range(len(file_list)):
        text = ( {True:'b',False:''}[gui3.data['b_{0}'.format(ii)]] +','+
                 {True:'d',False:''}[gui3.data['d_{0}'.format(ii)]] +','+
                 {True:'a',False:''}[gui3.data['a_{0}'.format(ii)]] +','+
                 {True:'t1',False:''}[gui3.data['t1_{0}'.format(ii)]] +','+
                 {True:'t2',False:''}[gui3.data['t2_{0}'.format(ii)]] +','+
                 {True:'z',False:''}[gui3.data['z_{0}'.format(ii)]] +','+
                 {True:'w2',False:''}[gui3.data['w2_{0}'.format(ii)]] +','+
                 {True:'w1',False:''}[gui3.data['w1_{0}'.format(ii)]] +','+
                 {True:'ws',False:''}[gui3.data['ws_{0}'.format(ii)]] +','+
                 {True:'wc',False:''}[gui3.data['wc_{0}'.format(ii)]] +','+
                 {True:'e',False:''}[gui3.data['e_{0}'.format(ii)]] +','+
                 gui3.data['extra_{0}'.format(ii)]
               ).replace(',,',',').replace(',,',',').replace(',,',',').replace(',,',',').replace(',,',',')
        if text[0] == ',':
            text = text[1:]
        if len(text) > 0:
            if text[-1] == ',':
                text = text[:-1]
        #file_list_full.append([ {True:'#',False:''}[gui3.data['comment_{0}'.format(ii)]] + file_list[ii][0].replace('#','').replace(' ',''),
        #                       gui3.data['fib1_{0}'.format(ii)] , gui3.data['fib2_{0}'.format(ii)], file_list[ii][3], file_list[ii][4], 
        #                       text ])     # Allows modification of the fiber content
        file_list_full.append([ {True:'#',False:''}[gui3.data['comment_{0}'.format(ii)]] + file_list[ii][0].replace('#','').replace(' ',''),
                               file_list[ii][1] , file_list[ii][2], file_list[ii][3], file_list[ii][4], gui3.data['obname_{0}'.format(ii)], text ])
        if not gui3.data['comment_{0}'.format(ii)]:     # Ignore the commented out lines
            file_list_new.append(file_list_full[-1])
            
    return file_list_new, file_list_full

def calibration_parameters_coordinates_UI(params, conf_data, object_information_full, object_information_head, CONFIGFILE):
    max_length = 1024
    def add_widgets(widgets, pkwargs, ii, ftype, doneftype):
        if '{0}_rawfiles'.format(ftype) not in conf_data.keys() or ftype in doneftype:
            return widgets, pkwargs, ii, doneftype
            
        widgets['type_{0}'.format(ii)] = dict(label=ftype, kind='Label', row=ii, column=0, columnspan=1, orientation=Tk.W)
        if 'master_{0}_filename'.format(ftype) in conf_data.keys() and (ftype.find('extract') == -1 or ftype.find('extract_combine') != -1) and ftype.find('wavelenoffset') == -1:
            elname = 'master_{0}_filename'.format(ftype)
            pkwargs[elname] = conf_data[elname]
            widgets[elname] = dict(kind='TextEntry', fmt=str, start=pkwargs[elname], width=30, row=ii, column=1, columnspan=3, orientation=Tk.W)
        elname = '{0}_calibs_create'.format(ftype)
        pkwargs[elname] = conf_data[elname]
        widgets[elname] = dict(kind='TextEntry', fmt=str, start=pkwargs[elname], width=50, row=ii, column=4, columnspan=5, orientation=Tk.W)
        label = conf_data['{0}_rawfiles'.format(ftype)]
        if len(label) > max_length:
            label = label[:max_length-3] + '...'        # Otherwise: X Error of failed request:  BadAlloc (insufficient resources for operation)
        widgets['files_{0}'.format(ii)]  = dict(label=label, kind='Label', row=ii, column=9, columnspan=1000, orientation=Tk.W)
        doneftype.append(ftype)
        ii += 1
        return widgets, pkwargs, ii, doneftype
    
    def vfunc_float(xs):
        try:
            value = float(xs)
            return True, value
        except:
            return False, ('Error, input must be float')
    
    # define widgets
    pkwargs = dict()
    widgets = dict()
    
    # Add widgets for fits_conf.txt
    widgets['type']        = dict(label='Type', kind='Label', row=0, column=0, columnspan=1, orientation=Tk.W)
    widgets['name_master'] = dict(label='Name of Master file', kind='Label', row=0, column=1, columnspan=3, orientation=Tk.W)
    widgets['calibs']      = dict(label='Calibrations to be applied', kind='Label', row=0, column=4, columnspan=5, orientation=Tk.W)
    widgets['files_used']  = dict(label='Files included (see last GUI or {0} to change the assigned files)'.format(params['raw_data_file_list']), kind='Label', row=0, column=9, columnspan=1000, orientation=Tk.W)
    darks = []
    for entry in sorted(params['exp_darks']):
        darks.append('dark{0}'.format(entry))
    doneftype = []
    ii = 1
    for ftype in ['bias'] + darks + ['rflat', 'trace1', 'trace2', 'blazecor', 'cal2_l', 'cal2_s', 'cal1_l', 'cal1_s', 'waveoffsetsci', 'waveoffsetcal']:
            widgets, pkwargs, ii, doneftype = add_widgets(widgets, pkwargs, ii, ftype, doneftype)
    for entry in sorted(conf_data.keys()):
        if entry.find('_rawfiles') >= 0:
            ftype = entry.replace('_rawfiles','')
            widgets, pkwargs, ii, doneftype = add_widgets(widgets, pkwargs, ii, ftype, doneftype)
    
    widgets['accept1'] = dict(label='Accept', kind='ExitButton', row=ii+1, column=6, rowspan=2, columnspan=2)
    ii += 3
    
    # Add widgets for object_list.txt
    widgets['name'] = dict(label='Object{0}name'.format(os.linesep), kind='Label', row=ii, column=0, rowspan=2, columnspan=1, orientation=Tk.E)
    widgets['header_info1'] = dict(label='Header information', kind='Label', row=ii, column=1, columnspan=2, orientation=Tk.E)
    widgets['header_ra'] = dict(label='RA', kind='Label', row=ii+1, column=1)
    widgets['header_dec'] = dict(label='DEC', kind='Label', row=ii+1, column=2)
    widgets['header_info2'] = dict(label='Header information', kind='Label', row=ii, column=3, columnspan=3)
    widgets['header_pmra'] = dict(label='PMRA', kind='Label', row=ii+1, column=3)
    widgets['header_pmdec'] = dict(label='PMRA', kind='Label', row=ii+1, column=4)
    widgets['header_epoch'] = dict(label='Epoch', kind='Label', row=ii+1, column=5)
    widgets['use_header'] = dict(label='Use{0}head'.format(os.linesep), kind='Label', row=ii, column=6, rowspan=2, columnspan=1)
    widgets['ra'] = dict(label='RA{0}<:>,< >,float'.format(os.linesep), kind='Label', row=ii, column=7, rowspan=2, columnspan=1)
    widgets['dec'] = dict(label='DEC{0}<:>,< >,float'.format(os.linesep), kind='Label', row=ii, column=8, rowspan=2, columnspan=1)
    widgets['pmra'] = dict(label='PMRA{0}[mas/yr]'.format(os.linesep), kind='Label', row=ii, column=9, rowspan=2, columnspan=1)
    widgets['pmdec'] = dict(label='PMDEC{0}[mas/yr]'.format(os.linesep), kind='Label', row=ii, column=10, rowspan=2, columnspan=1)
    widgets['epoch'] = dict(label='Epoch{0}(number)'.format(os.linesep), kind='Label', row=ii, column=11, rowspan=2, columnspan=1)
    widgets['simbad_name'] = dict(label='Simbad Name{0}(check that correct)'.format(os.linesep), kind='Label', row=ii, column=12, rowspan=2, columnspan=1)
    widgets['mask'] = dict(label='Mask (optional){0}G2,K5,M2'.format(os.linesep), kind='Label', row=ii, column=13, rowspan=2, columnspan=1)
    widgets['rot'] = dict(label='rotation (optional){0}[km/s]'.format(os.linesep), kind='Label', row=ii, column=14, rowspan=2, columnspan=1)
    ii += 2
    jj = 0                              # if no object to extract, jj will not exist, but is needed later
    for jj, entry in enumerate(object_information_full):
        widgets['name_{0}'.format(jj)]  = dict(label=entry[0], kind='Label', row=ii+jj, column=0, columnspan=1, orientation=Tk.E)
        state = None
        if entry[0].lower().find('sun') == 0 or entry[0].lower().find('moon') == 0 or entry[0].lower().find('jupiter') == 0:     # Check with get_object_site_from_header()
            for i in range(1,7):
                entry[i] = ''
            state = Tk.DISABLED
        pkwargs['ra_{0}'.format(jj)]    = entry[2]
        widgets['ra_{0}'.format(jj)]    = dict(kind='TextEntry', fmt=str, start=pkwargs['ra_{0}'.format(jj)], width=13, row=ii+jj, column=7, columnspan=1, state=state)
        pkwargs['dec_{0}'.format(jj)]   = entry[3]
        widgets['dec_{0}'.format(jj)]   = dict(kind='TextEntry', fmt=str, start=pkwargs['dec_{0}'.format(jj)], width=13, row=ii+jj, column=8, columnspan=1, state=state)
        pkwargs['pmra_{0}'.format(jj)]  = entry[4]
        widgets['pmra_{0}'.format(jj)]  = dict(kind='TextEntry', fmt=str, valid_function=vfunc_float, start=pkwargs['pmra_{0}'.format(jj)], width=7, row=ii+jj, column=9, columnspan=1, state=state)
        pkwargs['pmdec_{0}'.format(jj)] = entry[5]
        widgets['pmdec_{0}'.format(jj)] = dict(kind='TextEntry', fmt=str, valid_function=vfunc_float, start=pkwargs['pmdec_{0}'.format(jj)], width=7, row=ii+jj, column=10, columnspan=1, state=state)
        pkwargs['epoch_{0}'.format(jj)] = entry[9]
        widgets['epoch_{0}'.format(jj)] = dict(kind='TextEntry', fmt=str, valid_function=vfunc_float, start=pkwargs['epoch_{0}'.format(jj)], width=6, row=ii+jj, column=11, columnspan=1, state=state)
        # simbad name at end, next to button (button to search again is missing, therefore state=Tk.DISABLED
        pkwargs['simbad_name_{0}'.format(jj)] = entry[1]
        widgets['simbad_name_{0}'.format(jj)] = dict(kind='TextEntry', fmt=str, start=pkwargs['simbad_name_{0}'.format(jj)], width=13, row=ii+jj, column=12, columnspan=1, orientation=Tk.W, state=Tk.DISABLED)

        pkwargs['mask_{0}'.format(jj)]  = entry[7]
        widgets['mask_{0}'.format(jj)]  = dict(kind='TextEntry', fmt=str, start=pkwargs['mask_{0}'.format(jj)], width=5, row=ii+jj, column=13, columnspan=1)
        pkwargs['rot_{0}'.format(jj)]   = entry[8]
        widgets['rot_{0}'.format(jj)]   = dict(kind='TextEntry', fmt=str, valid_function=vfunc_float, start=pkwargs['rot_{0}'.format(jj)], width=5, row=ii+jj, column=14, columnspan=1)
        
        if object_information_head[jj][0] == '':
            pkwargs['use_header_{0}'.format(jj)]   = False
            widgets['use_header_{0}'.format(jj)]   = dict(label=None,  kind='CheckBox', start=pkwargs['use_header_{0}'.format(jj)], row=ii+jj, column=6, state=Tk.DISABLED)
            continue
        entry = object_information_head[jj]
        widgets['header_ra_{0}'.format(jj)]    = dict(label=entry[0], kind='Label', row=ii+jj, column=1)
        widgets['header_dec_{0}'.format(jj)]   = dict(label=entry[1], kind='Label', row=ii+jj, column=2)
        widgets['header_pmra_{0}'.format(jj)]  = dict(label=entry[2], kind='Label', row=ii+jj, column=3)
        widgets['header_pmdec_{0}'.format(jj)] = dict(label=entry[3], kind='Label', row=ii+jj, column=4)
        widgets['header_epoch_{0}'.format(jj)] = dict(label=entry[4], kind='Label', row=ii+jj, column=5)
        pkwargs['use_header_{0}'.format(jj)]   = True
        widgets['use_header_{0}'.format(jj)]   = dict(label=None,  kind='CheckBox', start=pkwargs['use_header_{0}'.format(jj)], row=ii+jj, column=6)
    
    ii += jj+1
    widgets['accept2'] = dict(label='Accept', kind='ExitButton', row=ii, column=6, rowspan=2, columnspan=2)
    ii += 2
    
    explain = 'Explanation of upper half:'+os.linesep+\
              '  This assigns the calibration that will be applied to files of the different types before creating the master file or before extracting the spectra. The following comma separated options are possible:'+os.linesep+\
              '     subframe, badpx_mask, bias, dark, flat, background, normalise, combine_mean, combine_sum'+os.linesep+\
              '     Please check the manual for more information on these options.'+os.linesep+\
              '  The assigned calibration steps are read from {0} (*_calibs_create_g){1}'.format(CONFIGFILE, os.linesep)+\
              '  The information from this part of the GUI will be stored in {0}.{1}{1}'.format(params['configfile_fitsfiles'], os.linesep)+\
              '  Explanation of the lower half:'+os.linesep+\
              '  For each object to be extracted the coordinates are derived/displayed here. If the header information is available then the information is shown and the user can decide if this information should be used.'+os.linesep+\
              '  The editable coordinates are taken from the file {0}, for which is also checked in the result and raw data path. If the object does not exist in the file, then Simbad is searched using the Object name.{1}'.format(params['object_file'], os.linesep)+\
              '  The object name for the results from Simbad is shown and should be the same as the object name derived from header/filename'+os.linesep+\
              '  The user can modify this information, which should be correct to perform correct barycentric correction. The information is then stored in {0}, overwriting the previous information.{1}'.format(params['object_file'], os.linesep)+\
              '     The RA and DEC can be given in hour or degree format. The optional parameters "mask" and "rotation speed" are used if RV analysis with the CERES pipeline is performed.'
    for jj, commentjj in enumerate(explain.split(os.linesep)):
        if len(commentjj) > 0:
            widgets['explain_{0}'.format(jj)] = dict(label=commentjj, kind='Label', row=ii+jj, column=0, columnspan=20, orientation=Tk.W )#, wraplength=100 )      
    
    wprops = dict(fullscreen=False )
    #wprops['width_GUI'] = 800   # not neccssary, as uses the automatic width
    
    gui3 = tkc.TkCanvasGrid(title='HiFLEx: Asigning Calibration and Coordinates', 
                            kwargs=pkwargs,widgets=widgets, widgetprops=wprops )
    gui3.master.mainloop()
    
    # Get the information from the GUI
    for entry in conf_data.keys():      # For the conf data
        if entry in gui3.data.keys():
            conf_data[entry] = gui3.data[entry]
            if entry.find('_calibs_create') != -1:
                if len(conf_data[entry]) > 0:
                    if conf_data[entry].find('[') + conf_data[entry].find(']') == -2:
                        conf_data[entry] = '[{0}]'.format(conf_data[entry])                     # add brackets
    
    object_information = []
    for ii in range(len(object_information_full)):
        enabled = {True:0,False:1}[gui3.data['use_header_{0}'.format(ii)]]
        if len(object_information_full[ii][1]) >= 3:            # The coordinates were read from a file
            try:
                object_information_full[ii][6] = float(object_information_full[ii][6])      # see if we can make it to a number
            except:
                True
            if type(object_information_full[ii][6]).__name__ != 'str':          # it's a number
                if abs(object_information_full[ii][6]) < 0.9:
                    enabled = False                                 # Object was disabled before
        if ( len(gui3.data['mask_{0}'.format(ii)]) < 1 and type(gui3.data['rot_{0}'.format(ii)]).__name__ == 'str' ) and ( len(gui3.data['ra_{0}'.format(ii)]) < 3 or len(gui3.data['ra_{0}'.format(ii)]) < 3 ):
            continue                                            # Both coordinates are not complete and Mask or rotation is not there -> no useful information
        if (len(gui3.data['ra_{0}'.format(ii)]) >= 3 and len(gui3.data['dec_{0}'.format(ii)]) >= 3):     # Fill pm* and epoch if coordinates are full
            for entry in ['pmra', 'pmdec']:
                if type(gui3.data['{1}_{0}'.format(ii, entry)]).__name__ == 'str':
                    gui3.data['{1}_{0}'.format(ii, entry)] = 0
            if type(gui3.data['epoch_{0}'.format(ii)]).__name__ == 'str':
                gui3.data['epoch_{0}'.format(ii)] = 2000
        object_information.append([ object_information_full[ii][0], gui3.data['ra_{0}'.format(ii)], gui3.data['dec_{0}'.format(ii)], 
                                    gui3.data['pmra_{0}'.format(ii)], gui3.data['pmdec_{0}'.format(ii)], gui3.data['epoch_{0}'.format(ii)],
                                    enabled, gui3.data['mask_{0}'.format(ii)], gui3.data['rot_{0}'.format(ii)] ])
    
    return conf_data, object_information

def adjust_binning_UI(im1, binxy, searchlimit_brightness=None, min_separation=None, userinput=True):
    """
    Adjusts the binning
    """
    if not userinput:
        return binxy, searchlimit_brightness, min_separation
    from procedures import bin_im
    
    #sim_sflat, dummy, dummy = bin_im(im_sflat, params['bin_search_apertures'] )
    # set up plot
    fig, frame = plt.subplots(1, 1)
    fig.set_size_inches(10, 7.5, forward=True)

    # get kwargs
    binx, biny = binxy
    pkwargs = dict(frame=frame, im1=im1, binx=binx, biny=biny, searchlimit_brightness=searchlimit_brightness, min_separation=min_separation,
                   im_scale_min=round_sig(np.percentile(im1,10),2), im_scale_max=round_sig(np.percentile(im1,90),2) )

    # define update plot function
    def plot(frame, im1, binx, biny, searchlimit_brightness, min_separation, im_scale_min, im_scale_max):
        frame.clear()
        title = ('Adjusting the binning')
        im_bin, dummy, dummy = bin_im(im1, [binx,biny] )
        frame = plot_img_spec.plot_image(im_bin, 'dummy_filename', pctile=[im_scale_min,im_scale_max], show=False, adjust=[0.07,0.95,0.95,0.07], title=title, return_frame=True, frame=frame, autotranspose=False, colorbar=False, axis_name=['Cross-dispersion axis [px]','Dispersion axis [px]','flux [ADU]'])
    
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
    def vfunc_float(xs):
        try:
            value = float(xs)
            return True, value
        except:
            return False, ('Error, input must be float')
    
    # define widgets
    widgets = dict()
    widgets['im_scale_min'] = dict(label='Scale the image (black)', 
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=str, start=pkwargs['im_scale_min'], valid_function=vfunc_float,
                                width=10)
    widgets['im_scale_max'] = dict(label='Scale the image (white)',
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=str, start=pkwargs['im_scale_max'], valid_function=vfunc_float,
                                width=10)
    widgets['binx'] = dict(label='Binning in{0}Dispersion axis'.format(os.linesep),
                           kind='TextEntry', minval=None, maxval=None,
                           fmt=str, start=pkwargs['binx'], valid_function=vfunc_int,
                           width=10)
    widgets['biny'] = dict(label='Binning in{0}Cross-dispersion axis'.format(os.linesep),
                           kind='TextEntry', minval=None, maxval=None,
                           fmt=str, start=pkwargs['biny'], valid_function=vfunc_int,
                           width=10)
    if searchlimit_brightness is not None:
        widgets['searchlimit_brightness'] = dict(label='Brightest pixel', comment='of the faintest order{0}in binned image to be traced.'.format(os.linesep), 
                           kind='TextEntry', minval=None, maxval=None,
                           fmt=str, start=pkwargs['searchlimit_brightness'], valid_function=vfunc_float,
                           width=10)
    if min_separation is not None:
        widgets['min_separation'] = dict(label='Minumum separation{0}of orders'.format(os.linesep), comment='in unbinned image.{0}Divide by Binning in{0}Cross-dispersion to{0}compare with image.'.format(os.linesep),
                           kind='TextEntry', minval=None, maxval=None,
                           fmt=str, start=pkwargs['min_separation'], valid_function=vfunc_float,
                           width=10)
    widgets['accept'] = dict(label='Accept Values',
                             kind='ExitButton',
                             position=Tk.BOTTOM)
    widgets['update'] = dict(label='Update', kind='UpdatePlot',
                             position=Tk.BOTTOM)
                             
    wprops = dict(orientation='v', position=Tk.RIGHT)

    gui3 = tkc.TkCanvas(figure=fig, ax=frame, func=plot, kwargs=pkwargs,
                        title='HiFlEx: Adjusting the binning', widgets=widgets,
                        widgetprops=wprops)
    gui3.master.mainloop()
    
    #binxy = pkwargs['binxy']       # This doesn't work, pkwargs are only updated within the plot function
    binxy = [ gui3.data['binx'], gui3.data['biny'] ]
    if 'searchlimit_brightness' in gui3.data.keys():
        searchlimit_brightness = gui3.data['searchlimit_brightness']
    if 'min_separation' in gui3.data.keys():
        min_separation = gui3.data['min_separation']
    
    plt.close()
    return binxy, searchlimit_brightness, min_separation

def remove_adjust_orders_UI(im1, pfits, xlows, xhighs, widths=[], shift=0, userinput=True, do_rm=False, do_adj=False, do_shft=False, do_add=False):
    """
    Removes orders from the pfits array in a GUI, allows to adjust the width of th extracted area
    :param do_adj: Adjust the width of the traces
    return fmask: array of bool with the same length as original orders, True for the orders to keep
    return pfits: same format as pfits, with adjusted parameters for the polynomial for the left and right border
    return widths: same format as widths, with adjusted left and right stop of the trace
    """
    from procedures import remove_orders, plot_traces_over_image, scale_image_plot, sort_traces_along_detector, update_tr_poly_width_multiplicate
    if not userinput or (not do_rm and not do_adj and not do_shft and not do_add):
        return remove_orders(pfits, []), pfits, widths, xlows, xhighs         # No orders to remove, pfits stays the same, widths stays the same
    
    # convert to numpy arrays
    pfits = np.array(pfits)
    xlows, xhighs = np.array(xlows), np.array(xhighs)
    im_log = scale_image_plot(im1,'log10')
    # Create a list with all data points to refit the orders
    data_x = np.linspace(0, im1.shape[0], int(0.01*im1.shape[0]), dtype=int) # about every 100 pixel in dispersion direction
    data_order_position = []
    for order in range(pfits.shape[0]):
        inorder = ( data_x > xlows[order] ) & (data_x < xhighs[order]-1)
        xarr = np.hstack(( xlows[order], data_x[inorder], xhighs[order]-1 ))           # add the boundaries as points
        if len(pfits.shape) == 3:
            yarr = np.polyval(pfits[order,0,1:], xarr-pfits[order,0,0])    
            yarrl = np.polyval(pfits[order,1,1:], xarr-pfits[order,1,0])
            yarrr = np.polyval(pfits[order,2,1:], xarr-pfits[order,2,0])
        else:
            yarr = np.polyval(pfits[order,1:], xarr-pfits[order,0]) 
            yarrl = yarr
            yarrr = yarr
        data_order_position.append( np.vstack((xarr, yarr, yarrl, yarrr)).T )             # contains the x and y values in one row

    # set up plot
    fig, frame = plt.subplots(1, 1)
    gui3 = tkc.TkCanvas(figure=None, ax=None, func=None, title='HiFLEx: Test', kwargs=dict(), widgets=dict(), widgetprops=dict(orientation='v', position=Tk.RIGHT) )      # Only to get the display size
    dpi=100
    fig.set_dpi(dpi)
    fig.set_size_inches( (int(gui3.screen_w_h[0]) - 400)/dpi, (int(gui3.screen_w_h[1])-90)/dpi, forward=True)     # Make the plot as big as possible
    tkc.TkCanvas.end(gui3)
    
    # im_log = scale_image_plot(im1,'log10')   # Don't do it here
    # get kwargs; add more to pkwargs, as otherwise what is changed in for example settings_modify_end() is lost
    pkwargs = dict(frame=frame, update_frame=True, im1=im1, pfits=pfits, xlows=xlows, xhighs=xhighs, widths=widths,
                   w_mult=1.0, rm_orders=[], shift=shift, im_scale_min=round_sig(np.percentile(im_log,10),2),
                   im_scale_max=round_sig(np.percentile(im_log,90),2), order=np.nan, addpointsorder_lcr=1, 
                   data_order_position=data_order_position, addpoints=False, removepoints=False)

    # define update plot function
    def plot(frame, update_frame, im1, pfits, xlows, xhighs, widths, w_mult, rm_orders, shift, im_scale_min, im_scale_max, order, addpointsorder_lcr, data_order_position, **kwargs):
        #print('addpointsorder_lcr',addpointsorder_lcr)
        x_range, y_range = copy.copy(frame.get_xlim()), copy.copy(frame.get_ylim())
        #if update_frame:        # This will speed up plotting, but copying the old frame doesn't work, hence it always adds plots to the frame, and frame.clear() before loading creates an clear frame, more in https://stackoverflow.com/questions/57351212/matplotlib-how-to-copy-a-contour-plot-to-another-figure
        if True:
            frame.clear()
            title = ''
            if do_adj:
                title += 'Defining the width of the traces. '
            if do_shft:
                title += 'Finding the shift of the traces. '
            if do_rm:
                title += 'Removing bad orders. '
            if do_add:
                title += 'Adding/modifying orders. '
            if do_rm or do_add:
                title += '{1}(Largest order number = {0}){1}'.format(pfits.shape[0]-1, os.linesep)
            mask = remove_orders(pfits, rm_orders)
            pfits_shift = copy.deepcopy(pfits)
            if len(pfits_shift.shape) == 3:
                pfits_shift[:,:,-1] += shift        # shift all traces
            else:
                pfits_shift[:,-1] += shift          # shift all traces
            frame = plot_traces_over_image(scale_image_plot(im1,'log10'), 'dummy_filename', pfits_shift, xlows, xhighs, widths, w_mult=w_mult, mask=mask, frame=frame, return_frame=True, imscale=[im_scale_min,im_scale_max])
            frame.set_title(title[:-1])
        #    pkwargs['backup_frame'] = copy.deepcopy(frame)# This will speed up plotting, but copying the old frame doesn't work, hence it always adds plots to the frame, and frame.clear() before loading creates an clear frame
        #else:# This will speed up plotting, but copying the old frame doesn't work, hence it always adds plots to the frame, and frame.clear() before loading creates an clear frame
        #    frame = copy.deepcopy(pkwargs['backup_frame'])# This will speed up plotting, but copying the old frame doesn't work, hence it always adds plots to the frame, and frame.clear() before loading creates an clear frame

        if do_add and not np.isnan(pkwargs['order']):
            yarr = data_order_position[order][:,addpointsorder_lcr ] + shift
            xarr = data_order_position[order][:,0]
            frame.plot(yarr, xarr, color='c', linewidth=1, linestyle='dotted', marker="x", markersize=5)
        if x_range != (0.0, 1.0) and y_range != (0.0, 1.0):
            frame.axis([x_range[0], x_range[1], y_range[0], y_range[1]])
    
    def settings_modify_order():
        pkwargs['backup_data_order_position_order'] = copy.copy(pkwargs['data_order_position'][pkwargs['order']])
        gui3.ws['addpointsorder'].config(state=Tk.DISABLED)
        if 'addpointsorder_lcr' in gui3.ws.keys():
            gui3.ws['addpointsorder_lcr'].config(state=Tk.DISABLED)
        gui3.ws['addpoints'].config(state=Tk.DISABLED)
        gui3.ws['removepoints'].config(state=Tk.DISABLED)
        gui3.ws['update'].config(state=Tk.DISABLED)
        gui3.ws['stoppoints'].config(state="normal")
        gui3.ws['cancelpoints'].config(state="normal")
        tkc.TkCanvas.update(gui3)
        pkwargs['update_frame'] = False     # When clicking only the points are redrawn
    
    def settings_modify_end():
        pkwargs['order'] = np.nan
        gui3.ws['addpointsorder'].config(state="normal")
        if 'addpointsorder_lcr' in gui3.ws.keys():
            gui3.ws['addpointsorder_lcr'].config(state="normal")
        gui3.ws['addpoints'].config(state="normal")
        gui3.ws['removepoints'].config(state="normal")
        gui3.ws['update'].config(state="normal")
        gui3.ws['stoppoints'].config(state=Tk.DISABLED)
        gui3.ws['cancelpoints'].config(state=Tk.DISABLED)
        pkwargs['addpoints'] = False
        pkwargs['removepoints'] = False
        pkwargs['update_frame'] = True
        tkc.TkCanvas.update(gui3)
    
    def add_an_order():
        pkwargs['xlows'] = np.hstack(( pkwargs['xlows'], im1.shape[0] ))
        pkwargs['xhighs'] = np.hstack(( pkwargs['xhighs'], 0 ))
        pfits = pkwargs['pfits']
        if len(pfits.shape) == 3:
            pfnew = np.zeros((1,pfits.shape[1],pfits.shape[2]))
            for ii in range(pfits.shape[1]):
                pfnew[0,ii,0] = np.median(pfits[:,ii,0])        # cen_px
                pfnew[0,ii,-1] = im1.shape[1]/2.                # put it in the middle
        else:
            pfnew = np.zeros((1, pfits.shape[1]))
            pfnew[0,0] = np.median(pfits[:,0])        # cen_px
            pfnew[0,-1] = im1.shape[1]/2.             # put it in the middle
        pkwargs['pfits'] = np.concatenate((pfits,pfnew), axis=0)
        pkwargs['data_order_position'].append( np.zeros((0,4)) )
        widths = pkwargs['widths']
        if len(widths) > 0:
            widthnew = np.median(widths,axis=0)
            #print(widths, widthnew)
            if type(widths).__name__ == 'ndarray':
                pkwargs['widths'] = np.concatenate( ( widths, widthnew.reshape((1,widths.shape[1])) ), axis=0)
            elif type(widths).__name__ == 'list':
                 pkwargs['widths'].append(list(widthnew))
        gui3.prompt_info('Added order with index {0}. Now you need to add points to it.'.format(pkwargs['pfits'].shape[0]-1), width=300)
        #print(pfits.shape, pfnew.shape, pkwargs['pfits'].shape, len(pkwargs['widths']), pkwargs['widths']  )
    
    def addpoints_action():
        tkc.TkCanvas.update(gui3, updateplot=False)
        pkwargs['order'] = gui3.data['addpointsorder']
        if type(pkwargs['order']).__name__ not in ['int']:
            return
        pkwargs['addpoints'] = True
        settings_modify_order()
    
    def removepoints_action():
        tkc.TkCanvas.update(gui3, updateplot=False)
        pkwargs['order'] = gui3.data['addpointsorder']
        if type(pkwargs['order']).__name__ not in ['int']:
            return
        pkwargs['removepoints'] = True
        settings_modify_order()
    
    def cancelpoints_action():
        # Restore the original points
        pkwargs['data_order_position'][pkwargs['order']] = pkwargs['backup_data_order_position_order']
        settings_modify_end()
    
    def stoppoints_action():
        order = pkwargs['order']
        data = pkwargs['data_order_position'][order]
        addpointsorder_lcr = pkwargs['addpointsorder_lcr']
        nonnan = ~np.isnan(data[:,addpointsorder_lcr])
        if len(pkwargs['pfits'].shape) == 3:
            poly_orders = pkwargs['pfits'].shape[2]-2
            cen_px = pkwargs['pfits'][order,addpointsorder_lcr-1,0]
        else:       # addpointsorder_lcr will be always 1
            poly_orders = pkwargs['pfits'].shape[1]-2
            cen_px = pkwargs['pfits'][order,0]
        if np.sum(nonnan) < 2*poly_orders:
            gui3.prompt_info('Warn: Needs at least {} points to fit the order. Please add more points or cancel.'.format(2*poly_orders), width=300)
            return
        pkwargs['xlows'][order] = int(np.min(data[:,0]))
        pkwargs['xhighs'][order] = int(np.ceil(np.max(data[:,0])+1))
        # Fit the center
        pfs = np.polyfit(data[nonnan,0] - cen_px, data[nonnan,addpointsorder_lcr], poly_orders)
        if len(pkwargs['pfits'].shape) == 3:
            pkwargs['pfits'][order,addpointsorder_lcr-1,1:] = pfs
        else:       # addpointsorder_lcr will be always 1
            pkwargs['pfits'][order,1:] = pfs
        settings_modify_end()
    
    def onclick(event):
        #print('clicked',event.xdata, event.ydata, event.inaxes, im1.shape)
        if event.xdata is None or event.ydata is None:
            return
        if not pkwargs['addpoints'] and not pkwargs['removepoints']:            # Not set
            return
        yy, xx = event.xdata-shift, event.ydata         # event.ydata: dispersion axis, xx in dispersion axis
        uncertainty = 15
        data = pkwargs['data_order_position'][pkwargs['order']]     # data[:,0]: dispersion axis
        addpointsorder_lcr = pkwargs['addpointsorder_lcr']
        if pkwargs['addpoints']:
            if xx < 0 or xx > im1.shape[0]-1 or yy < 0 or yy > im1.shape[1]-1:    #im1.shape[0]: disersion axis
                return                  # click was outside of the image
            datanew = [xx, np.nan, np.nan, np.nan]
            datanew[addpointsorder_lcr] = yy
            data = np.vstack((data, datanew))
            datasortindex = np.argsort(data[:,0])
            pkwargs['data_order_position'][pkwargs['order']] = data[datasortindex,:]
        elif pkwargs['removepoints']:
            yarr = copy.copy(data[:,addpointsorder_lcr])
            xarr = copy.copy(data[:,0])
            badvalues = np.isnan(yarr)
            xarr[badvalues] = -1E6
            yarr[badvalues] = -1E6
            diff = (xarr-xx)**2 + (yarr-yy)**2
            posmin = np.argmin(diff)
            #print('smallest distance', np.sqrt(diff[posmin]), data, posmin)
            if np.sqrt(diff[posmin]) > uncertainty:       # Too many pixel away
                return
            data[posmin,addpointsorder_lcr] = np.nan
            pkwargs['data_order_position'][pkwargs['order']] = data
        
        tkc.TkCanvas.update(gui3, updateplot=True)
        # !!!!!!!!!!!! do it for pressing lrc to work on boundaries of the orders? 
    
    # define valid_function
    # input is one variable (the string input)
    # return is either:
    #   True and values
    # or
    #   False and error message
    def vfunc(xs):
        try:
            xwhole = xs.replace(',', ' ')
            new_xs = xwhole.split()
            xs = []
            for nxs in new_xs:
                xs.append(int(nxs))
            return True, xs
        except:
            return False, ('Error, input must consist of integers'+os.linesep+\
                           'separated by commas or white spaces')
    def vfunc_int(xs):
        try:
            value = int(xs)
            return True, value
        except:
            return False, ('Error, input must be integer'+os.linesep)
    def vfunc_float(xs):
        try:
            value = float(xs)
            return True, value
        except:
            return False, ('Error, input must be float')
    def vfunc_lcr(xs):
        if xs == 'c':
            return True, 1
        elif xs == 'l':
            return True, 2
        elif xs == 'r':
            return True, 3
        else:
            return False, ('  Error, input must be either l, c, or r  ')
    
    # run initial update plot function
    plot(**pkwargs)
    
    # define widgets
    widgets = dict()
    widgets['im_scale_min'] = dict(label='Scale the image (black)', comment='The image will be shown in log10 scale',
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=str, start=pkwargs['im_scale_min'], valid_function=vfunc_float,
                                width=10)
    widgets['im_scale_max'] = dict(label='Scale the image (white)',
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=str, start=pkwargs['im_scale_max'], valid_function=vfunc_float,
                                width=10)
    if do_rm:
        widgets['rm_orders'] = dict(label='Select orders to remove',
                                comment='Enter all order numbers to remove'+os.linesep+\
                                        'separated by a whitespace or comma'+os.linesep+\
                                        'to undo just delete the entered number',
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=str, start=" ", valid_function=vfunc,
                                width=40)
    if do_adj:
        widgets['w_mult'] = dict(label='Multiplier for the'+os.linesep+'width of the traces',
                            comment='If the results are not as wished,'+os.linesep+\
                                    'a modification of the parameter "width_percentile"'+os.linesep+\
                                    'might help. To do this'+os.linesep+\
                                    'Cancel the script with CTRL+C in the terminal'+os.linesep+\
                                    'and then restart',
                            kind='TextEntry', minval=None, maxval=None,
                            fmt=float, start=pkwargs['w_mult'], valid_function=vfunc_float,
                            width=10)
    if do_shft:
        widgets['shift'] = dict(label='Shift the traces by:',
                            #comment='',
                            kind='TextEntry', minval=None, maxval=None,
                            fmt=float, start=0.0, valid_function=vfunc_float,
                            width=10)
    if do_add:
        widgets['addpointsorder'] = dict(label='Add/remove points from order:',
                            comment='It might be necessary to zoom before'+os.linesep+\
                                    'adding/removing points will work',
                            kind='TextEntry', minval=None, maxval=None,
                            fmt=int, start='', valid_function=vfunc_int,
                            width=10)
        if len(pkwargs['pfits'].shape) == 3:
            widgets['addpointsorder_lcr'] = dict(label='Which part of the order?',
                            comment='c: center (brightest part),'+os.linesep+\
                                    'l or r: left or right boundary of the order',
                            kind='TextEntry', minval=None, maxval=None,
                            fmt=str, start='c', valid_function=vfunc_lcr,
                            width=10)
        #widgets['addpoints'] = dict(label='Add or remove points from'+os.linesep+\
        #                               'the order given below',
        #                        kind='CheckBox', start=False)
    #widgets['spacer1'] = dict(kind='Spacer', position=Tk.BOTTOM)
    widgets['accept'] = dict(label='Accept', kind='ExitButton',
                             position=Tk.BOTTOM)
    widgets['update'] = dict(label='Update', kind='UpdatePlot',
                             position=Tk.BOTTOM)
    if do_add:
        widgets['addorder'] = dict(label='Add a new Order', kind='CommandButton', command=add_an_order,
                             position=Tk.BOTTOM, width=20)
        widgets['stoppoints'] = dict(label='Use added/removed points', kind='CommandButton', command=stoppoints_action,
                             position=Tk.BOTTOM, state=Tk.DISABLED, width=30)
        widgets['cancelpoints'] = dict(label='Cancel adding/removing points', kind='CommandButton', command=cancelpoints_action,
                             position=Tk.BOTTOM, state=Tk.DISABLED, width=30)
        widgets['removepoints'] = dict(label='Remove points from an Order', kind='CommandButton', command=removepoints_action,
                             position=Tk.BOTTOM, width=30)
        widgets['addpoints'] = dict(label='Add points to an Order', kind='CommandButton', command=addpoints_action,
                             position=Tk.BOTTOM, width=30)
    wprops = dict(orientation='v', position=Tk.RIGHT)
    
    gui3 = tkc.TkCanvas(figure=fig, ax=frame, func=plot, kwargs=pkwargs,
                        title='HiFlEx: Locating orders', widgets=widgets,
                        widgetprops=wprops)
    cid = fig.canvas.callbacks.connect('button_press_event', onclick)       # It works with GUI
    
    gui3.master.mainloop()

    pfits = pkwargs['pfits']
    xlows = pkwargs['xlows']
    xhighs = pkwargs['xhighs']
    widths = pkwargs['widths']
    w_mult = pkwargs['w_mult']
    shift = pkwargs['shift']
    rm_orders = pkwargs['rm_orders']
    
    # Remove bad orders
    good_data = (xhighs - xlows >= 5)      # enough data added to the order
    if len(pfits.shape) == 3:
        pfits = pfits[good_data,:,:]
    else:
        pfits = pfits[good_data,:]
    xlows = xlows[good_data]
    xhighs = xhighs[good_data]
    if len(widths) > 0:
        widths = widths[good_data,:]
    
    # Mask removed orders   
    fmask = remove_orders(pfits, rm_orders)
    
    # Rescale the widths
    if 'w_mult' != 1.0 and len(pfits.shape) == 3 and len(widths) > 0:
        pfits, widths = update_tr_poly_width_multiplicate(pfits, widths, [w_mult, w_mult], xlows, xhighs)
        
    # Shift orders
    if shift != 0.0:
        if len(pfits.shape) == 3:
            pfits[:,:,-1] += shift        # shift all traces
        else:
            pfits[:,-1] += shift          # shift all traces
    
    # Sort the orders
    if len(pfits.shape) == 3:
        sort_index = sort_traces_along_detector(pfits[:,0,:], xlows, xhighs)
        pfits = pfits[sort_index,:,:]
    else:
        sort_index = sort_traces_along_detector(pfits, xlows, xhighs)
        pfits = pfits[sort_index,:]
    xlows = xlows[sort_index]
    xhighs = xhighs[sort_index]
    if len(widths) > 0:
        widths = widths[sort_index,:]
    fmask = fmask[sort_index]
    
    plt.close()
    return fmask, pfits, widths, xlows, xhighs


def create_new_wavelength_UI( params, cal_l_spec, cal_s_spec, arc_lines_px, reference_lines_dict, adjust=[0.12,0.88,0.85,0.12, 1.0,1.01] ):   
    """
    :param arc_lines_px: numpy arry of floats: order, pixel, width, height of the line
    """
    from procedures import logger, read_text_file, convert_readfile, fit_basic_wavelength_solution, adjust_data_log
    reference_catalog, reference_names = reference_lines_dict['reference_catalog'][0], reference_lines_dict['reference_names'][0]
    px_to_wave_txt = read_text_file(params['px_to_wavelength_file'], no_empty_lines=True, warn_missing_file=False)              # list of strings
    if len(px_to_wave_txt) == 0:
        px_to_wave_txt = read_text_file(params['logging_path']+'tmp_'+params['px_to_wavelength_file'], no_empty_lines=True, warn_missing_file=False)              # list of strings
    px_to_wave = np.array( convert_readfile(px_to_wave_txt, [int, int, float, float], delimiter='\t', replaces=[['\n',''],[os.linesep,'']], shorten_input=True, replacewithnan=True ))     # order, real order, px, wave
    if px_to_wave.shape[0] > 0:
        px_to_wave = px_to_wave[~np.isnan(px_to_wave[:,0]),:]   # clear the values that don't have order number

    nr_entries = len(px_to_wave)
    if nr_entries > 0:
        order_offset = np.nan
        if np.prod( np.isnan(px_to_wave[:,1]) ) == 0:
            order_offset = int(np.nanmedian(px_to_wave[:,1] - px_to_wave[:,0]))
        px_to_wave = np.hstack(( px_to_wave, np.zeros((nr_entries,2))*np.nan, np.expand_dims(np.arange(nr_entries), axis=1), np.zeros((nr_entries,1))*np.nan  ))     # order, real order, px, wave, width, height of line, index, nan
        order = -1E6
        for entry in arc_lines_px:
            if entry[0] != order:       # So this step is only done when neccessary 
                order = entry[0]
                px_to_wave_sub = px_to_wave[ px_to_wave[:,0] == order, : ]
            diff = np.abs(px_to_wave_sub[:,2] - entry[1])
            posi = np.where(diff < 2)[0]
            if len(posi) >= 1:
                index = int(px_to_wave_sub[posi[0], 6])      # first entry as highest signal, not min(diff)
                px_to_wave[index,4] = entry[2]
                px_to_wave[index,5] = entry[3]
            else:
                px_to_wave = np.vstack(( px_to_wave, [order, order+order_offset, entry[1], np.nan, entry[2], entry[3], len(px_to_wave), np.nan ] ))
    else:
        order_offset = np.nan
        tmp = arc_lines_px[:,0]*np.nan
        px_to_wave = np.vstack(( arc_lines_px[:,0], tmp, arc_lines_px[:,1], tmp, arc_lines_px[:,2], arc_lines_px[:,3], tmp, tmp )).T
    
    px_to_wave[np.isnan(px_to_wave[:,4]),4] = -1
    px_to_wave[np.isnan(px_to_wave[:,5]),5] = -1
    px_to_wave[:,6] = np.nan
    order = int(len(cal_l_spec)/2)
    
    fig, ax = plt.subplots(3, 1, gridspec_kw={'height_ratios': [2, 5, 2]})
    gui3 = tkc.TkCanvasGrid(figure=None, ax=None, func=None, title='HiFLEx: Test', kwargs=dict(), widgets=dict(), widgetprops=dict() )      # Only to get the display size
    dpi=100
    fig.set_dpi(dpi)
    fig.set_size_inches( (int(gui3.screen_w_h[0]) - 400)/dpi, (int(gui3.screen_w_h[1])-90)/dpi, forward=True)     # Make the plot as big as possible
    tkc.TkCanvasGrid.end(gui3)
    plt.subplots_adjust(left=adjust[0], right=adjust[1], top=adjust[2], bottom=adjust[3])
    pkwargs = dict()
    
    def mark_line_in_GUI(px):
        gui3.repopulate()
        px_to_wave = pkwargs['px_to_wave']
        px_to_wave_sub = px_to_wave[ px_to_wave[:,0] == pkwargs['order'], : ]
        match = np.where( np.abs(px_to_wave_sub[:,2] - px) < 20)[0]
        if len(match) == 0:
            unmark_line_in_GUI()
        else:
            index = match[0]
            for widget in gui3.ws.keys():
                if widget.find('wave_') > -1:
                    try:
                        number = int(widget.replace('wave_',''))
                    except:
                        continue
                    if number != index:
                        gui3.ws[widget].config(state=Tk.DISABLED)
        
    def unmark_line_in_GUI():
        gui3.repopulate()
    
    def onclick(event):
        #print('clicked',event.xdata, event.ydata, event.inaxes)
        if event.xdata is None or event.ydata is None:
            unmark_line_in_GUI()
        else:
            px = event.xdata
            mark_line_in_GUI(px)

    def plot(ax, cal_l_spec, cal_s_spec, px_to_wave, order, **pkwarks):
        for ii in range(3):
            ax[ii].clear()
            ordii = order+ii-1
            if ordii < 0 or ordii >= len(cal_l_spec):         # outside of useful range
                continue
            y1 = adjust_data_log( cal_l_spec[ordii,:] )
            y2 = adjust_data_log( cal_s_spec[ordii,:] )
            x = list(range(len(y1)))
            title   = None
            labels  = []
            x_title = ''
            y_title = 'log ( Flux [ADU] )'
            if ii == 0:
                title = 'Plot of the spectrum for previous ({0}), actual ({1}), and next order ({2})'.format(order-1, order, order+1)
                labels = ['long'+os.linesep+'exp', 'short'+os.linesep+'exp']
            if ii == 2:
                x_title = 'Dispersion direction [px]'
            ax[ii] = plot_img_spec.plot_points([x,x], [y1,y2], labels, [], show=False, adjust=[0.05,0.99,0.97,0.05, 0.92,1.2], title='', return_frame=True, frame=ax[ii], x_title=x_title, y_title=y_title, linestyle="-", marker="")
            
            pctl = 2./len(y1)*100                   # Exclude the two highest and two lowest points
            xmin, xmax = np.nanmin(x), np.nanmax(x)
            ymin1, ymax1 = np.nanpercentile(y1, pctl), np.nanpercentile(y1, 100-pctl)
            ymin2, ymax2 = np.nanpercentile(y2, pctl), np.nanpercentile(y2, 100-pctl)
            ymin, ymax = min(ymin1, ymin2), max(ymax1, ymax2)
            y_range = ymax - ymin
            x_range = xmax - xmin
            if np.isnan(y_range) or np.isnan(x_range):
                continue
            ax[ii].set_ylim([ymin-y_range*0.01, ymax+y_range*0.01])
            ax[ii].set_xlim([xmin-x_range*0.02, xmax+x_range*0.02])
            title = {0:'Previous', 1:'', 2:'Next'}[ii]
            ax[ii].set_title('{0} Order: {1}'.format(title, ordii), fontsize=11)
            
            # Plot the wavelengths from the reference list
            px_to_wave_sub = px_to_wave[ px_to_wave[:,0] == ordii, : ]
            if ii == 1:
                goodentries = 0
            else:
                goodentries = 1E6           # to only plot the identified lines
            for jj, entry in enumerate(px_to_wave_sub):
                if entry[3] is np.nan or str(entry[3]) == 'nan':
                    if goodentries > pkwargs['max_lines_list']:
                        continue
                    text = '{0}'.format(round(entry[2],1))
                    color = 'g'
                    goodentries += 1
                else:
                    if ii == 1:
                        text = '{0} - {1}'.format(round(entry[2],1), round(entry[3],4))     # px and wave
                    else:
                        text = '{0}'.format(round(entry[3],4))                              # wave
                    color = 'r'
                ax[ii].plot( [entry[2],entry[2]], [ymin,ymin+y_range*0.01], color=color )
                ax[ii].text( entry[2], ymin+y_range*0.015, text, fontsize=8, 
                                    horizontalalignment='center', verticalalignment='bottom', rotation=90, color=color, zorder=5 )

            x_f = np.arange(xmin-x_range*0.02, xmax+x_range*0.02, 0.01)
            w_fa, w_fl = [], []
            if 'poly_lin_{0}'.format(ordii) in pkwargs.keys() and ii == 1:
                    w_fl = np.polyval(pkwargs['poly_lin_{0}'.format(ordii)], x_f)           # Linear fit
            if 'wavelength_solution' in pkwargs.keys():       # when wavelength solution is not available
                    wavelength_solution = pkwargs['wavelength_solution']
                    w_fa = np.polyval(wavelength_solution[ordii,2:], x_f-wavelength_solution[ordii,1])  # Fit all wavelengths
            if len(reference_catalog) == 0:
                    continue
            pos_index = 0
            ytop = ymax+y_range*0.01
            if np.sum(~np.isnan(px_to_wave_sub[:,3])) == 0 and px_to_wave_sub.shape[0] > 0:
                px_to_wave_sub[0,3] = -10000                    # Dummy wavelength
            for kk, w_f in enumerate([w_fl, w_fa]):
                if len(w_f) < 10:
                    continue
                inorder = (reference_catalog[:,0] >= np.min(w_f)) & (reference_catalog[:,0] <= np.max(w_f))
                reference_catalog_sub = np.array(sorted(reference_catalog[inorder,:], key=operator.itemgetter(1), reverse=True ))
                if len(reference_catalog_sub) == 0:
                    continue
                y_scale = y_range / (np.max(reference_catalog_sub[:,1])+0.0)
                pos_index += 1
                num_notident = 1
                for color in {0:['b', 'r'], 1:['g', 'r']}[kk]:        # plot the green/blue lines before the red ones
                    for jj,refline in enumerate(reference_catalog_sub):
                        if px_to_wave_sub[:,3].shape[0] > 0:    temp = px_to_wave_sub[:,3]
                        else:                                   temp = [-1000]
                        if np.nanmin(np.abs(refline[0] - temp)) < 1E-2:
                            if color != 'r':
                                continue        # don't plot a matching line when it's not red lines to be plotted
                        else:
                            if color == 'r':
                                continue        # don't plot the non matching lines when red lines to be plotted
                            if num_notident > pkwargs['max_reflines']:
                                break           # stop plot non-matching lines
                            num_notident += 1   
                        index = np.argmin(np.abs(refline[0] - w_f))
                        x_pos = x_f[index]
                        #if np.isnan(spec1[order, x_position]):
                        #    continue
                        #y_position = np.nanmax(spec1[order, max(0,x_position-2):min(len(spec1[order,:]),x_position+2)])
                        #??? = max(y_pos, y_position+0.23*y_range)
                        text = '{0} - {1}'.format( round(refline[0],2+2*kk), reference_names[int(refline[2])] )
                        if kk == 0:                                         # Text below marker
                            y1 = [ymax, ymax+y_range*0.005+y_scale*0.03*refline[1]]
                            y2 = ymax
                            align = 'top'
                            ytop = max(y1[1], ytop)
                        else:                                               # Text above marker
                            y_pos = ymax - y_range*0.47                        # where is the second line
                            y1 = [y_pos-y_scale*0.03*refline[1]-y_range*0.005, y_pos]
                            y2 = y_pos+y_range*0.01
                            align = 'bottom'
                        ax[ii].plot( [x_pos,x_pos], y1, color=color )
                        ax[ii].text( x_pos, y2, text, fontsize=8, horizontalalignment='center', verticalalignment=align, rotation=90, color=color, zorder=5 )
            ax[ii].set_ylim([ymin-y_range*0.01, ytop])
    
    def calculate_linear_solution(pkwargs, order=None):
        px_to_wave = pkwargs['px_to_wave']
        if order is None:
            orders = np.unique(px_to_wave[:,0]).astype(int)
        else:
            orders = [order]
        for order in orders:
            inorder = px_to_wave[:,0] == order
            px_to_wave_sub = px_to_wave[ inorder, : ]
            good_values = ~np.isnan(px_to_wave_sub[:,3])
            if np.sum(good_values) <= 1:
                px_to_wave[inorder,7] = np.nan
                if 'poly_lin_{0}'.format(order) in pkwargs.keys():
                    del pkwargs['poly_lin_{0}'.format(order)]
                continue
            poly = np.polyfit( px_to_wave_sub[good_values,2], px_to_wave_sub[good_values,3], 1)     # 1 for linear fit
            px_to_wave[inorder,7] = np.polyval(poly, px_to_wave[inorder,2])
            pkwargs['poly_lin_{0}'.format(order)] = poly
        pkwargs['px_to_wave'] = px_to_wave
        return pkwargs
    
    def calculate_wavesolution_calc(px_to_wave, cal_l_spec):
        px_to_wave_sub = px_to_wave[ ~np.isnan(px_to_wave[:,3]), : ]

        wavelength_solution, wavelength_solution_arclines, wavesol_result_txt = fit_basic_wavelength_solution(params, px_to_wave_sub, cal_l_spec, 'GUI')
        # can be [], [], in case there are not enough lines
        
        for order in range(len(wavelength_solution)):
            values_order = px_to_wave[:,0] == order
            px_to_wave_sub = px_to_wave[ values_order, : ]
            px_to_wave[ values_order, 6 ] = np.polyval(wavelength_solution[order,2:], px_to_wave[ values_order, 2 ]-wavelength_solution[order,1])
        return px_to_wave, wavelength_solution, wavelength_solution_arclines, wavesol_result_txt
    
    def calculate_wavesolution():
        update_order(updateplot=False)                      # To make sure the data is read again
        px_to_wave = pkwargs['px_to_wave']
        if np.sum(np.isnan(px_to_wave[:,1])) > 0:
            gui3.prompt('The order offset is not yet defined.{0}{0}Please insert the correct number.'.format(os.linesep))
            return
        pkwargs['wavelength_solution_info'] = True
        px_to_wave, wavelength_solution, wavelength_solution_arclines, wavesol_result_txt = calculate_wavesolution_calc(px_to_wave, cal_l_spec)
        if len(wavelength_solution) == 0:
            return
        pkwargs['px_to_wave']                   = px_to_wave
        pkwargs['wavelength_solution']          = wavelength_solution
        pkwargs['wavelength_solution_arclines'] = wavelength_solution_arclines
        # Replot
        """gui3.funkwargs['notclose'] = True
        tkc.TkCanvasGrid.end(gui3)          # will update also the plot"""
        if pkwargs['wavelength_solution_info']:
            gui3.prompt_info(wavesol_result_txt, width=600)
            pkwargs['wavelength_solution_info'] = False
        update_order()                      # To update the wavelengths
        
    def remove_widgets():
        variables = [gui3.validation, gui3.fmts, gui3.entries, gui3.funcs, gui3.onclickvalue, gui3.data]
        for entry in pkwargs['widgets_change'].keys():
            if entry in pkwargs.keys():
                del pkwargs[entry]
            gui3.ws[entry].destroy()
            for variable in variables:
                if entry in variable.keys():
                    del variable[entry]
    
    def add_widgets(pkwargs):
        order = pkwargs['order']
        px_to_wave = pkwargs['px_to_wave']
        px_to_wave_sub = px_to_wave[ px_to_wave[:,0] == order, : ]
        if px_to_wave_sub.shape[0] > pkwargs['max_lines_list']:        # Limit to 80 lines
            px_to_wave_sub = px_to_wave_sub[:pkwargs['max_lines_list'],:]
        pkwargs['oldmax_lines_list'] = copy.copy(pkwargs['max_lines_list'])
        wid_sub = dict()
        offset = 5
        ii = 0              # if the order doesn't contain any entries
        for ii, entry in enumerate(px_to_wave_sub):
            pkwargs['delete_{0}'.format(ii)]   = False
            wid_sub['delete_{0}'.format(ii)]   = dict(label=None,  kind='CheckBox', start=pkwargs['delete_{0}'.format(ii)], row=ii+offset, column=0)
            wid_sub['px_pos_{0}'.format(ii)] = dict(label='%6.1f'%entry[2], kind='Label', row=ii+offset, column=1, rowspan=1, columnspan=1, orientation=Tk.W)
            wave = {True:'', False:str(round(entry[3],6))}[entry[3] is np.nan or str(entry[3]) == 'nan']
            pkwargs['wave_{0}'.format(ii)] = wave
            wid_sub['wave_{0}'.format(ii)] = dict(kind='TextEntry', fmt=str, start=pkwargs['wave_{0}'.format(ii)], width=10, row=ii+offset, column=2, columnspan=2)
            text = '{0}  {1}'.format('%4.1f'%entry[4], '%5.1i'%entry[5])
            wid_sub['lineprop_{0}'.format(ii)] = dict(label=text, kind='Label', row=ii+offset, column=4, rowspan=1, columnspan=2, orientation=Tk.W)
            text = ''
            if np.isnan(entry[6]) == False and entry[6] > 1000 and entry[6] < 1E6:
                # Add button to copy the value
                text += '%8.3f'%entry[6]
            text += ' - '
            if np.isnan(entry[7]) == False and entry[7] > 1000 and entry[7] < 1E6:
                text += '%8.3f'%entry[7]
            wid_sub['wave_sol_{0}'.format(ii)] = dict(label=text, kind='Label', row=ii+offset, column=6, rowspan=1, columnspan=2, orientation=Tk.W)
        for jj in range(3):
            pkwargs['px_pos_extra_{0}'.format(jj)] = ''
            wid_sub['px_pos_extra_{0}'.format(jj)] = dict(kind='TextEntry', fmt=str, start=pkwargs['px_pos_extra_{0}'.format(jj)], width=6, row=ii+1+jj+offset, column=1, columnspan=1)
            pkwargs['wave_extra_{0}'.format(jj)] = ''
            wid_sub['wave_extra_{0}'.format(jj)] = dict(kind='TextEntry', fmt=str, start=pkwargs['wave_extra_{0}'.format(jj)], width=10, row=ii+1+jj+offset, column=2, columnspan=2)
        
        pkwargs['widgets_change'] = wid_sub
        return pkwargs, wid_sub
    
    def read_data(gui3):
        oldorder = gui3.funkwargs['oldorder']
        px_to_wave = gui3.funkwargs['px_to_wave']
        params['tmp_polynom_order_traces'] = [gui3.funkwargs['disporder']]
        params['tmp_polynom_order_intertraces'] = [gui3.funkwargs['crossdisporder']]
        keep = np.ones(( len(px_to_wave) )).astype(bool)
        # Read the results
        indexes = np.where(px_to_wave[:,0] == oldorder)[0]
        if indexes.shape[0] > pkwargs['oldmax_lines_list']:        # Limit to 80 lines
            indexes = indexes[:pkwargs['oldmax_lines_list']]
        for ii,index in enumerate(indexes):
            #if 'wave_{0}'.format(ii) not in gui3.data.keys():       # can be empty when adding the extra stuff
            #    continue
            if gui3.data['delete_{0}'.format(ii)]:                  # delete
                keep[index] = False
            data = gui3.data['wave_{0}'.format(ii)]
            good_data = False
            if len(data) > 0:
                try:
                    data = float(data)
                    good_data = True
                except:
                    print('Warn: Can not convert entry into a number: {0}'.format(data))
            if good_data:
                px_to_wave[index,3] = data
            else:
                px_to_wave[index,3] = np.nan
        px_to_wave = px_to_wave[keep,:]
        for jj in range(3):
            good_data = False
            px = gui3.data['px_pos_extra_{0}'.format(jj)]
            wave = gui3.data['wave_extra_{0}'.format(jj)]
            if len(px) > 0 and len(wave) > 0:
                try:
                    px = float(px)
                    good_data = True
                except:
                    print('Warn: Can not convert entry into a number: {0}'.format(data))
                try:
                    wave = float(wave)
                except:
                    good_data = False
                    print('Warn: Can not convert entry into a number: {0}'.format(data))
            if good_data:
                px_to_wave = np.vstack(( px_to_wave, [oldorder, np.nan, px, wave, -1, -1, np.nan, np.nan ] ))
        if type(gui3.funkwargs['order_offset']).__name__ == 'int':
            px_to_wave[:,1] = px_to_wave[:,0] + gui3.funkwargs['order_offset']
        save_px_to_wave(px_to_wave, params['logging_path']+'tmp_'+params['px_to_wavelength_file'])
        gui3.funkwargs['px_to_wave'] = px_to_wave
        gui3.funkwargs = calculate_linear_solution(gui3.funkwargs, order=oldorder)
        
        #print 'new', gui3.data['order'], 'old', gui3.funkwargs['oldorder'], gui3.funkwargs['px_to_wave'][:,3]
    
    def update_order(updateplot=False):
        tkc.TkCanvasGrid.update(gui3, updateplot=updateplot)               # Read the values: Update order, get the new wavelengths
        read_data(gui3)              # read the new wavelengths, necessary here, as switched off at closing the gui
        remove_widgets()
        gui3.funkwargs['oldorder'] = copy.copy(gui3.data['order'])
        gui3.funkwargs, wid_sub = add_widgets(gui3.funkwargs)
        gui3.widgets = wid_sub
        gui3.add_widgets(parent=gui3.scrollable_frame)
        tkc.TkCanvasGrid.update(gui3)               # Update to get the plot with the new fitted wavelengths
    
    def save_px_to_wave(px_to_wave, fname):
        with open(fname,'w') as file:
            for entry in px_to_wave:
                if np.isnan(entry[0]):
                    continue
                entry = list(entry)
                for i in [1,3]:
                    if np.isnan(entry[i]) or str(entry[i]) == 'nan':
                        entry[i] = ''
                    elif i == 1:
                        entry[i] = int(entry[i])
                file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}'.format( int(entry[0]), entry[1], round(entry[2],2), entry[3], round(entry[4],2), round(entry[5],2), os.linesep ))
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
            return False, ('Error, input must be integer'+os.linesep)
    def vfunc_float(xs):
        try:
            value = float(xs)
            return True, value
        except:
            return False, ('Error, input must be float'+os.linesep)
    
    # get kwargs
    pkwargs = dict(ax=ax, cal_l_spec=cal_l_spec, cal_s_spec=cal_s_spec, px_to_wave=px_to_wave, order=order, order_offset=order_offset, max_reflines=30, max_lines_list=80)
    pkwargs['wavelength_solution_info'] = False
    if np.sum(np.isnan(px_to_wave[:,1])) == 0:
        px_to_wave, wavelength_solution, wavelength_solution_arclines, wavesol_result_txt = calculate_wavesolution_calc(px_to_wave, cal_l_spec)
        if len(wavelength_solution) > 0:
            pkwargs['wavelength_solution'] = wavelength_solution
            pkwargs['px_to_wave']          = px_to_wave
        pkwargs = calculate_linear_solution(pkwargs)
    # run initial update plot function
    plot(**pkwargs)
    # define widgets
    widgets = dict()
    widgets['ordertext']  = dict(label='Order:', kind='Label', row=0, column=0, rowspan=1, columnspan=4, orientation=Tk.W)
    pkwargs['order'] = order
    pkwargs['oldorder'] = order
    widgets['order'] = dict(kind='TextEntry', fmt=int, start=pkwargs['order'], valid_function=vfunc_int, minval=-1E-9, maxval=len(cal_l_spec)-1+1E-9, width=5, row=0, column=4, columnspan=1)
    widgets['update'] = dict(label='Update', kind='CommandButton', command=update_order, row=0, column=5, columnspan=2, width=9)
    widgets['accept'] = dict(label='Accept', kind='ExitButton', row=0, column=7, width=4)     # Replace by new fuction, as issue will occur if use changed the order but didn't updata
    widgets['order_offsettext']  = dict(label='Order offset between arbitary{0}numbering (starting at 0){0}and real physical order:'.format(os.linesep), kind='Label', row=1, column=0, rowspan=1, columnspan=4, orientation=Tk.W)
    pkwargs['order_offset'] = {True:pkwargs['order_offset'], False:''}[pkwargs['order_offset'] is not np.nan]
    widgets['order_offset'] = dict(kind='TextEntry', fmt=int, start=pkwargs['order_offset'], valid_function=vfunc_int, minval=-1000, maxval=1000, width=5, row=1, column=4, columnspan=1)
    widgets['calculate_wavesolution'] = dict(label='Calculate{0}Wavelength{0}solution'.format(os.linesep), kind='CommandButton', command=calculate_wavesolution, row=1, column=5, columnspan=2, width=9)
    widgets['disptxt']   = dict(label='Polynom-order for{0}dispersion direction:'.format(os.linesep), kind='Label', row=2, column=0, rowspan=1, columnspan=3, orientation=Tk.W)
    pkwargs['disporder'] = max(params.get('tmp_polynom_order_traces' ,params['polynom_order_traces']))
    widgets['disporder'] = dict(kind='TextEntry', fmt=int, start=pkwargs['disporder'], valid_function=vfunc_int, minval=1-1E-9, maxval=100, width=4, row=2, column=3, columnspan=1)
    widgets['crossdisptxt']   = dict(label='Polynom-order for cross-{0}dispersion direction:'.format(os.linesep), kind='Label', row=2, column=4, rowspan=1, columnspan=3, orientation=Tk.W)
    pkwargs['crossdisporder'] = max(params.get('tmp_polynom_order_intertraces' ,params['polynom_order_intertraces']))
    widgets['crossdisporder'] = dict(kind='TextEntry', fmt=int, start=pkwargs['crossdisporder'], valid_function=vfunc_int, minval=1-1E-9, maxval=100, width=4, row=2, column=7, columnspan=1)
    
    widgets['txtmax_reflines'] = dict(label='Max reflines in plot:', kind='Label', row=3, column=0, rowspan=1, columnspan=3, orientation=Tk.W)
    widgets['max_reflines'] = dict(kind='TextEntry', fmt=int, start=pkwargs['max_reflines'], valid_function=vfunc_int, minval=0, maxval=10000, width=4, row=3, column=3, columnspan=1)
    widgets['txtmax_lines_list'] = dict(label='Max identified lines in{0}list below/plot'.format(os.linesep), kind='Label', row=3, column=4, rowspan=1, columnspan=3, orientation=Tk.W)
    widgets['max_lines_list'] = dict(kind='TextEntry', fmt=int, start=pkwargs['max_lines_list'], valid_function=vfunc_int, minval=0, maxval=10000, width=4, row=3, column=7, columnspan=1)
    
    widgets['txtdelete']   = dict(label='Del-{0}ete'.format(os.linesep), kind='Label', row=4, column=0, rowspan=1, columnspan=1, orientation=Tk.W)
    widgets['txtpx_pos']   = dict(label='Pixel', kind='Label', row=4, column=1, rowspan=1, columnspan=1, orientation=Tk.W)
    widgets['txtwave']     = dict(label='Wavelength', kind='Label', row=4, column=2, rowspan=1, columnspan=2, orientation=Tk.W)
    widgets['txtlineprop'] = dict(label='Line width{0}+ height'.format(os.linesep), kind='Label', row=4, column=4, rowspan=1, columnspan=2, orientation=Tk.W)   # Can't be float, as otherwise deleting an entry wouldn't work
    widgets['txtwave_sol'] = dict(label='Wavelength-solution{0}global - linear'.format(os.linesep), kind='Label', row=4, column=6, rowspan=1, columnspan=2, orientation=Tk.W)   # Can't be float, as otherwise deleting an entry wouldn't work
    pkwargs, wid_sub = add_widgets(pkwargs)
    widgets.update(wid_sub)
    
    wprops = dict(fullscreen=False )
    gui3 = tkc.TkCanvasGrid(figure=fig, ax=ax, func=plot, title='HiFLEx: Create a new wavelength solution', 
                            kwargs=pkwargs, widgets=widgets, widgetprops=wprops )
    cid = fig.canvas.callbacks.connect('button_press_event', onclick)       # It works with GUI
    fig.set_size_inches( (int(gui3.screen_w_h[0]) - gui3.width_GUI-50)/dpi, (int(gui3.screen_w_h[1])-90)/dpi, forward=True)     # Make the plot as big as possible
        
    gui3.master.mainloop()
        
    read_data(gui3)             # read the new wavelengths only when closing the window
    plt.close()
    
    px_to_wave = gui3.funkwargs['px_to_wave']
    if os.path.isfile(params['px_to_wavelength_file']):
        os.system('mv {0} old_{0}'.format(params['px_to_wavelength_file']))
    save_px_to_wave(px_to_wave, params['px_to_wavelength_file'])
    del params['tmp_polynom_order_traces']
    del params['tmp_polynom_order_intertraces']
    
    if np.sum(np.isnan(px_to_wave[:,1])) == 0:
        px_to_wave, wavelength_solution, wavelength_solution_arclines, wavesol_result_txt = calculate_wavesolution_calc(px_to_wave, cal_l_spec)
    else:
        logger('Error: Not able to find a wavelength solution with the information provided by the user. Please rerun the program.', params=params)
    
    return dict(wavesol=wavelength_solution, reflines=wavelength_solution_arclines)


def correlate_px_wave_result_UI(params, cal_l_spec, cal_s_spec, wavelength_solution, wavelength_solution_arclines, reference_catalog, reference_names, plot_log=False, adjust=[0.12,0.88,0.85,0.12, 1.0,1.01]):
    
    from procedures import logger, adjust_data_log, create_wavelengths_from_solution, plot_overlapping_orders
    
    plot_pos = [0.01, 0.035, 0.04]      # begin line, end line, beginn of text; all relative to the the flux at the position and proportional to the range of the flux

    y_title = 'Flux [ADU]'
    
    if plot_log:
        y_title = 'log10 ( {0} )'.format(y_title)
        cal_l_spec = adjust_data_log(cal_l_spec)
        cal_s_spec = adjust_data_log(cal_s_spec)
    reference_catalog = np.array(sorted(reference_catalog, key=operator.itemgetter(1), reverse=True ))          # Sort by intensity in reverse order
    xarr = np.arange(cal_l_spec.shape[1])
    
    order = int(cal_l_spec.shape[0]/2)
    
    fig, ax = plt.subplots(3, 1, gridspec_kw={'height_ratios': [2, 5, 2]})
    gui3 = tkc.TkCanvasGrid(figure=None, ax=None, func=None, title='HiFLEx: Test', kwargs=dict(), widgets=dict(), widgetprops=dict() )      # Only to get the display size
    dpi=100
    fig.set_dpi(dpi)
    fig.set_size_inches( (int(gui3.screen_w_h[0]) - 400)/dpi, (int(gui3.screen_w_h[1])-90)/dpi, forward=True)     # Make the plot as big as possible
    tkc.TkCanvasGrid.end(gui3)
    plt.subplots_adjust(left=adjust[0], right=adjust[1], top=adjust[2], bottom=adjust[3])
    pkwargs = dict()
    
    def plot(ax, cal_l_spec, cal_s_spec, order, max_reflines, max_lines_list):
        for ii in range(3):
            ax[ii].clear()
            ordii = order+ii-1
            if ordii < 0 or ordii >= cal_l_spec.shape[0]:         # outside of useful range
                continue
            title   = None
            labels  = []
            x_title = ''
            if ii == 0:
                title = 'Plot of the spectrum for previous ({0}), actual ({1}), and next order ({2})'.format(order-1, order, order+1)
            if ii == 1:
                labels = ['long exp', 'short exp']
            if ii == 2:
                x_title = 'Wavelength [Angstroms]'
                    
            wavelengths, dummy = create_wavelengths_from_solution(params, wavelength_solution, cal_l_spec)
            x_range = wavelengths[ ordii, ~np.isnan(cal_l_spec[ordii,:]) ]
            if x_range.shape[0] < 10:                                       # the reference trace cod fall completely out of the CCD
                continue
            good_values = np.where((reference_catalog[:,0] >= np.min(x_range)) & (reference_catalog[:,0] <= np.max(x_range)) )[0]     # lines in the right wavelength range
            reference_lines = reference_catalog[good_values,:]                                                                  # possible reference lines in the order
            
            num_ident = len( np.where(np.array(wavelength_solution_arclines[ordii]) > 0)[0] )                                            # np.where only works on arrays
            #if num_ident == 0:                                          # Add a dummy value in case no lines are available for this order, as otherwise it will crash ; !! Not necessary after making everything into an numpy array?
            #    wavelength_solution_arclines[order].append(-1E5)
            num_notident = len(reference_lines)
            
            # Create a the plotting data
            x_data, y_data, label, color = plot_overlapping_orders(ordii, wavelengths, cal_l_spec, cal_s_spec, labels)        # create the overlapping pixel
            ax[ii] = plot_img_spec.plot_points(x_data, y_data, label, [], show=False, adjust=[0.05,0.95,0.97,0.05, 0.92,1.2], title=title, 
                                      return_frame=True, frame=ax[ii], x_title=x_title, y_title=y_title, linestyle="-", marker="",
                                      color=color)
            #ymin, ymax = axes.get_ylim()
            #y_range = ymax - ymin
            #xmin, xmax = axes.get_xlim()        # wavelengths contains also values for nans
            # Redifine the axis using only the current order
            pctl = 5./len(cal_l_spec[ordii,:])*100                   # Exclude the five highest and five lowest points
            xmin, xmax = np.nanmin(wavelengths[ordii,:]), np.nanmax(wavelengths[ordii,:])
            ymin1, ymax1 = np.nanpercentile(cal_l_spec[ordii,:], pctl), np.nanpercentile(cal_l_spec[ordii,:], 100-pctl)
            ymin2, ymax2 = np.nanpercentile(cal_s_spec[ordii,:], pctl), np.nanpercentile(cal_s_spec[ordii,:], 100-pctl)
            ymin, ymax = min(ymin1, ymin2), max(ymax1, ymax2)
            ymin, ymax = min(ymin1, ymin2), max(ymax1, ymax2)
            y_range = ymax - ymin
            x_range = xmax - xmin
            
            if len(reference_lines) > 0:
                y_scale = (plot_pos[1] - plot_pos[0]) / (max(reference_lines[:,1])+0.0)
                num_notident = 1
                for color in ['g', 'r']:        # plot the green lines before the red ones
                    for refline in reference_lines:
                        if np.min(np.abs(refline[0] - wavelength_solution_arclines[ordii])) < 1E-2:
                            if color != 'r':
                                continue        # don't plot a matching line when it's not red lines to be plotted
                        else:
                            if color == 'r':
                                continue        # don't plot the non matching lines when red lines to be plotted
                            if num_notident > max_reflines:
                                break           # stop plot non-matching lines
                            num_notident += 1   
                        x_position = np.argmin(np.abs(refline[0] - wavelengths[ordii,:]))
                        if np.isnan(cal_l_spec[ordii, x_position]):
                            continue
                        y_position = np.nanmax(cal_l_spec[ordii, max(0,x_position-2):min(len(cal_l_spec[ordii,:]),x_position+2)])
                        ymax = max(ymax, y_position+0.23*y_range)
                        ax[ii].plot( [refline[0],refline[0]], [y_position+plot_pos[0]*y_range,y_position+(plot_pos[0]+y_scale*refline[1])*y_range], color=color )
                        ax[ii].text( refline[0], y_position+plot_pos[2]*y_range, '{0} - {1}'.format(round(refline[0],4),reference_names[int(refline[2])]), 
                                    fontsize=8, horizontalalignment='center', verticalalignment='bottom', rotation=90, color=color, zorder=5 )
                    
            ax[ii].set_ylim(ymin,ymax)
            ax[ii].set_xlim(xmin-0.01*x_range,xmax+0.01*x_range)
            
            """wave_to_px = [ [], [], [] ]
            px_to_wave = [ [], [], [] ]
            secax_x = [ [], [], [] ]
            wave_to_px[ii] = lambda wave: np.polyval(poly_inv, wave)
            px_to_wave[ii] = lambda px: np.polyval(wavelength_solution[ordii, 2:], px - wavelength_solution[ordii,1])"""
            
            if ii == 1:
                # Add top x-axis
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', np.RankWarning)
                    poly_inv = np.polyfit(wavelengths[ordii,:], xarr, wavelength_solution.shape[1]-2-1)        # fit x(y), might cause rank warning, 
                def wave_to_px(wave):
                    """Convert wavelength to pixel"""
                    px = np.polyval(poly_inv, wave)
                    return px
                def px_to_wave(px):
                    """Conver pixel to wavelength"""
                    if ordii < cal_l_spec.shape[0]:
                        wave = np.polyval(wavelength_solution[ordii, 2:], px - wavelength_solution[ordii,1])
                        return wave
                    else:
                        return px
                try:
                    secax_x = ax[ii].secondary_xaxis('top', functions=(wave_to_px, px_to_wave))
                    #secax_x.set_xlabel('x [px]')
                except:
                    logger('Info: you are using an old matplotlib version ({0}). Consider updating it.'.format(matplotlib.__version__))

               
            """    
            y1 = adjust_data_log( cal_l_spec[ordii,:] )
            y2 = adjust_data_log( cal_s_spec[ordii,:] )
            x = list(range(len(y1)))
            title   = None
            labels  = []
            x_title = ''
            y_title = 'log ( Flux [ADU] )'
            if ii == 0:
                title = 'Plot of the spectrum for previous ({0}), actual ({1}), and next order ({2})'.format(order-1, order, order+1)
                labels = ['long'+os.linesep+'exp', 'short'+os.linesep+'exp']
            if ii == 2:
                x_title = 'Dispersion direction [px]'
            ax[ii] = plot_img_spec.plot_points([x,x], [y1,y2], labels, [], show=False, adjust=[0.05,0.99,0.97,0.05, 0.92,1.2], title='', return_frame=True, frame=ax[ii], x_title=x_title, y_title=y_title, linestyle="-", marker="")
            
            pctl = 2./len(y1)*100                   # Exclude the two highest and two lowest points
            xmin, xmax = np.nanmin(x), np.nanmax(x)
            ymin1, ymax1 = np.nanpercentile(y1, pctl), np.nanpercentile(y1, 100-pctl)
            ymin2, ymax2 = np.nanpercentile(y2, pctl), np.nanpercentile(y2, 100-pctl)
            ymin, ymax = min(ymin1, ymin2), max(ymax1, ymax2)
            y_range = ymax - ymin
            x_range = xmax - xmin
            if np.isnan(y_range) or np.isnan(x_range):
                continue
            ax[ii].set_ylim([ymin-y_range*0.01, ymax+y_range*0.01])
            ax[ii].set_xlim([xmin-x_range*0.02, xmax+x_range*0.02])
            title = {0:'Previous', 1:'', 2:'Next'}[ii]
            ax[ii].set_title('{0} Order: {1}'.format(title, ordii), fontsize=11)
            
            # Plot the wavelengths from the reference list
            px_to_wave_sub = px_to_wave[ px_to_wave[:,0] == ordii, : ]
            if ii == 1:
                goodentries = 0
            else:
                goodentries = 1E6           # to only plot the identified lines
            for jj, entry in enumerate(px_to_wave_sub):
                if entry[3] is np.nan or str(entry[3]) == 'nan':
                    if goodentries > pkwargs['max_lines_list']:
                        continue
                    text = '{0}'.format(round(entry[2],1))
                    color = 'g'
                    goodentries += 1
                else:
                    if ii == 1:
                        text = '{0} - {1}'.format(round(entry[2],1), round(entry[3],4))     # px and wave
                    else:
                        text = '{0}'.format(round(entry[3],4))                              # wave
                    color = 'r'
                ax[ii].plot( [entry[2],entry[2]], [ymin,ymin+y_range*0.01], color=color )
                ax[ii].text( entry[2], ymin+y_range*0.015, text, fontsize=8, 
                                    horizontalalignment='center', verticalalignment='bottom', rotation=90, color=color, zorder=5 )

            x_f = np.arange(xmin-x_range*0.02, xmax+x_range*0.02, 0.01)
            w_fa, w_fl = [], []
            if 'poly_lin_{0}'.format(ordii) in pkwargs.keys() and ii == 1:
                    w_fl = np.polyval(pkwargs['poly_lin_{0}'.format(ordii)], x_f)           # Linear fit
            if 'wavelength_solution' in pkwargs.keys():       # when wavelength solution is not available
                    wavelength_solution = pkwargs['wavelength_solution']
                    w_fa = np.polyval(wavelength_solution[ordii,2:], x_f-wavelength_solution[ordii,1])  # Fit all wavelengths
            if len(reference_catalog) == 0:
                    continue
            pos_index = 0
            ytop = ymax+y_range*0.01
            if np.sum(~np.isnan(px_to_wave_sub[:,3])) == 0 and px_to_wave_sub.shape[0] > 0:
                px_to_wave_sub[0,3] = -10000                    # Dummy wavelength
            for kk, w_f in enumerate([w_fl, w_fa]):
                if len(w_f) < 10:
                    continue
                inorder = (reference_catalog[:,0] >= np.min(w_f)) & (reference_catalog[:,0] <= np.max(w_f))
                reference_catalog_sub = np.array(sorted(reference_catalog[inorder,:], key=operator.itemgetter(1), reverse=True ))
                if len(reference_catalog_sub) == 0:
                    continue
                y_scale = y_range / (np.max(reference_catalog_sub[:,1])+0.0)
                pos_index += 1
                num_notident = 1
                for color in {0:['b', 'r'], 1:['g', 'r']}[kk]:        # plot the green/blue lines before the red ones
                    for jj,refline in enumerate(reference_catalog_sub):
                        if np.nanmin(np.abs(refline[0] - px_to_wave_sub[:,3])) < 1E-2:
                            if color != 'r':
                                continue        # don't plot a matching line when it's not red lines to be plotted
                        else:
                            if color == 'r':
                                continue        # don't plot the non matching lines when red lines to be plotted
                            if num_notident > pkwargs['max_reflines']:
                                break           # stop plot non-matching lines
                            num_notident += 1   
                        index = np.argmin(np.abs(refline[0] - w_f))
                        x_pos = x_f[index]
                        #if np.isnan(spec1[order, x_position]):
                        #    continue
                        #y_position = np.nanmax(spec1[order, max(0,x_position-2):min(len(spec1[order,:]),x_position+2)])
                        #??? = max(y_pos, y_position+0.23*y_range)
                        text = '{0} - {1}'.format( round(refline[0],2+2*kk), reference_names[int(refline[2])] )
                        if kk == 0:                                         # Text below marker
                            y1 = [ymax, ymax+y_range*0.005+y_scale*0.03*refline[1]]
                            y2 = ymax
                            align = 'top'
                            ytop = max(y1[1], ytop)
                        else:                                               # Text above marker
                            y_pos = ymax - y_range*0.47                        # where is the second line
                            y1 = [y_pos-y_scale*0.03*refline[1]-y_range*0.005, y_pos]
                            y2 = y_pos+y_range*0.01
                            align = 'bottom'
                        ax[ii].plot( [x_pos,x_pos], y1, color=color )
                        ax[ii].text( x_pos, y2, text, fontsize=8, horizontalalignment='center', verticalalignment=align, rotation=90, color=color, zorder=5 )
            ax[ii].set_ylim([ymin-y_range*0.01, ytop])
            """
            
    
    
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
            return False, ('Error, input must be integer'+os.linesep)
    def vfunc_float(xs):
        try:
            value = float(xs)
            return True, value
        except:
            return False, ('Error, input must be float'+os.linesep)
    
    # get kwargs
    pkwargs = dict(ax=ax, cal_l_spec=cal_l_spec, cal_s_spec=cal_s_spec, order=order, max_reflines=50, max_lines_list=100)
    
    # run initial update plot function
    plot(**pkwargs)
    # define widgets
    widgets = dict()
    widgets['colorexplain'] = dict(label='Red:{0}identified{0}catalogue lines{0}{0}Green:{0}unidentified lines'.format(os.linesep),
                                kind='Label')
    widgets['order'] = dict(label='Order',
                                comment='which order{0}in the center'.format(os.linesep),
                                kind='TextEntry', minval=-1E-9, maxval=cal_l_spec.shape[0]-1+1E-9,
                                fmt=int, start=pkwargs['order'], valid_function=vfunc_int,
                                width=6)
    widgets['max_reflines'] = dict(label='Maximum number of{0}reference lines'.format(os.linesep),
                                comment='(unidentified)' ,
                                kind='TextEntry', minval=0, maxval=1E8,
                                fmt=int, start=pkwargs['max_reflines'], valid_function=vfunc_int,
                                width=10)
    #widgets['max_lines_list'] = dict(label='Max number of{0}identified lines'.format(os.linesep),
    #                            kind='TextEntry', minval=0, maxval=1E8,
    #                            fmt=int, start=pkwargs['max_lines_list'], valid_function=vfunc_int,
    #                            width=10)
    widgets['accept'] = dict(label='Close', kind='ExitButton', position=Tk.BOTTOM)
    widgets['update'] = dict(label='Update Plot', kind='UpdatePlot', position=Tk.BOTTOM)

    wprops = dict(orientation='v', position=Tk.RIGHT)

    gui3 = tkc.TkCanvas(figure=fig, ax=ax, func=plot, kwargs=pkwargs,
                        title='HiFlEx: Check the identified and correlated emission lines', widgets=widgets,
                        widgetprops=wprops)
    
    fig.set_size_inches( (int(gui3.screen_w_h[0]) - gui3.width_GUI-50)/dpi, (int(gui3.screen_w_h[1])-90)/dpi, forward=True)     # Make the plot as big as possible
        
    gui3.master.mainloop()

    plt.close()
    
    


def correlate_px_wave_result_UI_old(im, arc_lines_wavelength, reference_catalog, arc_lines_px, reference_names, wavelength_solution, wavelength_solution_arclines, adjust=[0.12,0.88,0.85,0.12, 1.0,1.01]): # replace by create_new_wavelength_UI
    ims = im.shape
            
    fig, frame = plt.subplots(1, 1)
    plt.subplots_adjust(left=adjust[0], right=adjust[1], top=adjust[2], bottom=adjust[3])
    
    order = 10
    high_limit = 98
    arc_stretch = 0.5
    
    def plot(frame, im, order, high_limit, arc_stretch):
        x_range, y_range = copy.copy(frame.get_xlim()), copy.copy(frame.get_ylim())
        frame.clear()
        """try:
            order = gui3.data['order']
        except:
            order = order
        try:
            high_limit = gui3.data['high_limit']
        except:
            high_limit = high_limit
        try:
            arc_stretch = gui3.data['arc_stretch']
        except:
            arc_stretch = arc_stretch"""
        
        #print 'order, arc_setting', order, arc_setting
        title = ('Order {0}: identified emission lines: blue [px],{1}'+\
                 'catalog lines: red (0.1px precission, name and wavelength at the botom),{1}'+\
                 'corellated lines: green (px, wavelength at top)').format(order, os.linesep)
        
        # Plot the extracted arc spectrum
        xarr = np.arange(len(im[order,:]))
        yarr = im[order,:]
        notnans = ~np.isnan(yarr)
        yarr = yarr[notnans]
        xarr = xarr[notnans]
        yarr_max = max(0,np.percentile(yarr,high_limit))
        yarr[yarr > yarr_max] = yarr_max
        frame = plot_img_spec.plot_points([xarr], [yarr], [], [], show=False, adjust=adjust, title=title, 
                                      return_frame=True, frame=frame, x_title='x [px]', y_title='flux [ADU]', linestyle="-", marker="")
        
        # Plot the correlated lines in the extracted spectrum
        arc_line_order = arc_lines_px[(arc_lines_px[:,0] == order),:]
        xarr = np.vstack(( arc_line_order[:, 1], arc_line_order[:, 1])).T
        yplot = np.repeat( [[np.min(yarr), np.max(yarr)]], len(arc_line_order), axis=0)
        frame = plot_img_spec.plot_points(xarr, yplot, [], [], show=False, adjust=adjust, title=title, 
                                      return_frame=True, frame=frame, x_title='x [px]', y_title='flux [ADU]', linestyle="-", marker="", color='b')

        # Plot the catalog lines using the fit
        xarr_catalog = np.arange(-100,ims[1]+100, 0.1)
        xarr_wave = np.polyval(wavelength_solution[order,2:], xarr_catalog-wavelength_solution[order,1])
        inorder = ( reference_catalog[:,0] >= np.min(xarr_wave) ) & ( reference_catalog[:,0] <= np.max(xarr_wave) )
        reference_catalog_sub = reference_catalog[inorder, :]
        reference_names_sub = np.array(reference_names)[inorder]
        arc_stretch *= (np.max(yarr)-np.min(yarr)) / np.max(reference_catalog_sub[:,1])      # Rescale the arc lines so that arc_stretch=1 fills it just once
        for line_index in range(reference_catalog_sub.shape[0]):
            diff = np.abs(reference_catalog_sub[line_index,0] - xarr_wave)
            xarc = xarr_catalog[np.argmin(diff)]
            frame.plot( [xarc,xarc],[np.min(yarr),np.min(yarr)+reference_catalog_sub[line_index,1]*arc_stretch], color='r' )
            frame.text(xarc, np.min(yarr)+reference_catalog_sub[line_index,1]*arc_stretch, '{0} {1}'.format(reference_names_sub[line_index],reference_catalog_sub[line_index,0]),
                            horizontalalignment='center', verticalalignment='bottom', rotation=90, color='k', zorder=5)
           
        # Plot the identified catalog lines
        arc_line_order = arc_lines_wavelength[(arc_lines_wavelength[:,0] == order),:]
        xarr = np.vstack(( arc_line_order[:, 1], arc_line_order[:, 1])).T
        yplot = np.repeat( [[np.min(yarr), np.max(yarr)]], len(arc_line_order), axis=0)
        frame = plot_img_spec.plot_points(xarr, yplot, [], [], show=False, adjust=adjust, title=title, 
                                      return_frame=True, frame=frame, x_title='x [px]', y_title='flux [ADU]', linestyle="-", marker="", color='g')
        for arc_line in arc_line_order:
            frame.text(arc_line[1], np.max(yarr), r'{0} $\pm$ {1}'.format(arc_line[2], round(arc_line[3],3) ), horizontalalignment='center', verticalalignment='top', rotation=90, color='k', zorder=5)

    # get kwargs
    pkwargs = dict(frame=frame, im=im, order=order, high_limit=high_limit, arc_stretch=arc_stretch)
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
    def vfunc_float(xs):
        try:
            value = float(xs)
            return True, value
        except:
            return False, ('Error, input must be float')
            
    # define widgets
    widgets = dict()
    widgets['order'] = dict(label='Order',
                                comment='which order to show' ,
                                kind='TextEntry', minval=0, maxval=ims[0]-1,
                                fmt=int, start=order, valid_function=vfunc_int,
                                width=10)
    widgets['high_limit'] = dict(label='High limit [%]',
                                comment='Reject highest pixels' ,
                                kind='TextEntry', minval=1, maxval=105,
                                fmt=float, start=high_limit, valid_function=vfunc_float,
                                width=10)
    widgets['arc_stretch'] = dict(label='Scale of{0}catalogue lines'.format(os.linesep),
                                comment='float value' ,
                                kind='TextEntry', minval=None, maxval=None,
                                fmt=float, start=arc_stretch, valid_function=vfunc_float,
                                width=10)
    widgets['accept'] = dict(label='Close', kind='ExitButton', position=Tk.BOTTOM)
    widgets['update'] = dict(label='Update', kind='UpdatePlot', position=Tk.BOTTOM)

    wprops = dict(orientation='v', position=Tk.RIGHT)

    gui3 = tkc.TkCanvas(figure=fig, ax=frame, func=plot, kwargs=pkwargs,
                        title='HiFlEx: Check the identified and correlated emission lines', widgets=widgets,
                        widgetprops=wprops)
    gui3.master.mainloop()
    plt.close()





