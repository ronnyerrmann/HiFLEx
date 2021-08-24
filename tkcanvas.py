#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 06/07/17 at 4:56 PM

Last Modified on March 31 2020

@author: Neil Cook
@modified extensively: Ronny Errmann

Version 0.0.11
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
backend = 'TkAgg'           # The NavigationToolbar doesn't contain the settings for x,y zoom or line styles
# backend = 'Qt4Agg'        # This requires more work
if backend == 'TkAgg':
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    try:        # after version 2.2
        from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk as NavigationToolbar2TkAgg
    except:     # before version 2.2
        from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
    from matplotlib.backend_bases import cursors
    import matplotlib.backends.backend_tkagg as tkagg
else:
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvasTkAgg
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar2TkAgg
    from matplotlib.backend_bases import cursors
    import matplotlib.backends.backend_qt4agg as tkagg

import time

if sys.version_info[0] < 3:
    import Tkinter as Tk
    from collections import OrderedDict as dict
    import ttk
    import Queue as queue
else:
    import tkinter as Tk
    import tkinter.ttk as ttk
    import queue

# =============================================================================
# Define variables
# =============================================================================

def get_curr_screen_geometry():
    """
    Workaround to get the size of the current screen in a multi-screen setup.

    Returns:
        geometry (str): The standard Tk geometry string.
            [width]x[height]+[left]+[top]
    """
    root = Tk.Tk()
    root.update_idletasks()
    root.attributes('-fullscreen', True)
    root.state('iconic')
    geometry = root.winfo_geometry()
    root.destroy()
    return geometry

class TkCanvas:
    def __init__(self, figure=None, ax=None, **kwargs):         # "=None" added by Ronny; to use it just as text field
        try:
            screeen_geometry = get_curr_screen_geometry()
        except:
            print('Warn: Could not get the Screen Geometry: This terminal needs to be able to acces a graphical display.')
            screeen_geometry = '1024x768+1+1'
        self.master = Tk.Tk()
        self.screen_w_h = (str(screeen_geometry).replace('-','+').split('+')[0]).split('x')       # list of width and height

        # Adopted from https://blog.tecladocode.com/tkinter-scrollable-frames/
        container = ttk.Frame(self.master)
        canvas = Tk.Canvas( container)
        scrollbarx = ttk.Scrollbar(container, orient="horizontal", command=canvas.xview)
        scrollbary = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind("<Configure>",lambda e: canvas.configure( scrollregion=canvas.bbox("all") ) )
        # canvas.bbox("all") gives a 4-value tuple describing the position of two corners of a rectangle which is the scroll region
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(xscrollcommand=scrollbarx.set, yscrollcommand=scrollbary.set)
        
        container.pack(side=Tk.RIGHT)
        canvas.pack(side="left", fill="both", expand=True)
        scrollbarx.pack(side="bottom", fill="x")
        scrollbary.pack(side="right", fill="y")
        
        self.scrollable_frame = scrollable_frame
        
        # extract keyword arguments
        title = kwargs.get('title', 'Matplotlib window')
        self.master.title(title)

        # get figure
        self.figure = figure
        self.ax = ax
        self.func = kwargs.get('func', None)
        self.funkwargs = kwargs.get('kwargs', dict())
        self.widgets = kwargs.get('widgets', dict())
        self.wprops = kwargs.get('widgetprops',
                                 dict(orientation='h', position='bottom'))
        self.bardir = self.wprops['orientation']

        # default validation command args
        self.vcmdargs = ['%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W']
        # validation variables
        self.current_problem = ""

        # label variables
        self.statusmessage = Tk.StringVar(value="Ready")

        # empty data
        self.data = dict()
        self.entries = dict()
        self.funcs = dict()
        self.validation = dict()
        self.fmts = dict()
        self.onclickvalue = dict()

        # populate matplotlib drawing area and buttons
        self.populate()
        
        #self.width_GUI  = max(300, min( self.wprops.get('width_GUI' , 1E9), int(self.screen_w_h[0])-50 ) )
        #self.height_GUI = max(200, min( self.wprops.get('height_GUI', 1E9), int(self.screen_w_h[1])-90 ) )
        self.width_GUI  = max(300, min( self.wprops.get('width_GUI' , 1E9), int(self.screen_w_h[0])-50, scrollable_frame.winfo_width() ) )
        self.height_GUI = max(200, min( self.wprops.get('height_GUI', 1E9), int(self.screen_w_h[1])-90, scrollable_frame.winfo_height() ) )
        canvas.config(width=self.width_GUI, height=self.height_GUI)
        canvas.update()     # adds the winfo information


    def populate(self):

        self.button_bar(self.master, orientation=self.bardir)
        if self.figure is not None:                  # Added by Ronny
            self.matplotlib_window()
            self.add_progressbar(self.master, orientation='h')
        self.add_widgets(parent=self.buttonbar)

    def depopulate(self):
        for widget in self.ws.keys():
            self.ws[widget].config(state='disable')

    def repopulate(self):
        for widget in self.ws.keys():
            self.ws[widget].config(state="normal")

    def matplotlib_window(self):
        self.mwindow = FigureCanvasTkAgg(self.figure, master=self.master)
        self.mwindow.get_tk_widget().config(bg='white')
        # add matplotlib toolbar
        toolbar = NavigationToolbar2TkAgg(self.mwindow, self.master)
        toolbar.update()
        self.mwindow.get_tk_widget().pack(side=Tk.TOP, fill='both',
                                          expand='YES')
        self.mwindow.draw()         # 20/12 replaced show with draw
        #self.mwindow.show()        20/12 replaced show with draw
        try:
            # fix native matplotlib cursor bug  # fails in version 2.2
            tkagg.cursord[cursors.POINTER] = ""
        except:
            "do nothing"
        #extra room at top
        #plt.subplots_adjust(top=0.75)      # commented out by ronny
        self.mwindow.draw()

    def button_bar(self, master, orientation='h'):
        self.buttonbar = Tk.Frame(master)
        if 'h' in orientation.lower():
            self.buttonbar.pack(side=Tk.BOTTOM, fill='y')
        else:
            self.buttonbar.pack(side=Tk.RIGHT, fill='x')

    def add_widgets(self, parent, orientation='h'):
        widgets = self.widgets
        self.ws = dict()
        if 'h' in orientation.lower():
            dpos = Tk.TOP
        else:
            dpos = Tk.LEFT

        for widget in widgets:
            kind = widgets[widget].get('kind', 'TextEntry')
            position = widgets[widget].get('position', dpos)
            name = widget
            label = widgets[widget].get('label', 'Enter {0}'.format(name))
            start = widgets[widget].get('start', 0)
            fmt = widgets[widget].get('fmt', int)
            minval = widgets[widget].get('minval', None)
            maxval = widgets[widget].get('maxval', None)
            comment = widgets[widget].get('comment', None)
            result = widgets[widget].get('result', None)
            onclick = widgets[widget].get('onclick', None)
            vfunc = widgets[widget].get('valid_function', None)
            width = widgets[widget].get('width', 10)
            orientation = widgets[widget].get('orientation', None)
            wraplength = widgets[widget].get('wraplength', None)
            state = widgets[widget].get('state', None)
            command = widgets[widget].get('command', None)
            
            kwargs = dict(width=width, minval=minval, maxval=maxval, valid_function=vfunc, comment=comment, state=state)
            
            if kind == 'Label':                     #added by ronny
                w = Tk.Label(parent, text=label, state=state)    #added by ronny
                w.pack(side=position, anchor=Tk.NW)                    #added by ronny
            elif kind == 'Spacer':                  #added by ronny
                tframe = Tk.Frame(parent)
                tframe.pack(side=Tk.LEFT, anchor=Tk.NW, padx=100)
                parent = tframe
                w = None
            elif kind == 'TextEntry':
                if name not in self.data:
                    self.data[name] = start
                w = self.add_TextEntry(parent, name=name, desc=label, fmt=fmt,
                                       position=position, minval=minval,
                                       maxval=maxval, comment=comment,
                                       valid_function=vfunc, width=width)
            elif kind == 'CheckBox':                    #added by ronny
                if name not in self.data:                #added by ronny
                    self.data[name] = start                #added by ronny
                w = self.add_CheckBox(parent, name=name, desc=label,
                                      position=position, comment=comment)                #added by ronny
            elif kind == 'ExitButton':
                if result is not None:
                    self.data[name] = result
                    self.onclickvalue[name] = onclick
                if onclick is None:
                    w = self.add_button(parent, label, self.end, side=position)
                else:
                    w = self.add_button(parent, label,
                                        command=lambda name=name:
                                                self.valueend(name),
                                        side=position)
            elif kind == 'UpdatePlot':
                w = self.add_button(parent, label, self.update,
                                    side=position)
                if result is not None:
                    self.data[name] = result
                    self.fmts[name] = type(result)
            elif kind == 'CommandButton':
                w = self.add_button(parent, label, command, **kwargs)
                if result is not None:
                    self.data[name] = result
                    self.fmts[name] = type(result)
            else:
                w = None
            # append to storage
            if w is not None:
                #self.ws.append(w)
                self.ws[name] = w

    #---------------------------------------------------------------------------
    #define sub widgets
    #---------------------------------------------------------------------------
    def add_TextEntry(self, master, name, desc, fmt, position, **kwargs):

        background = kwargs.get('background', 'white')
        foreground = kwargs.get('foreground', 'black')
        width = kwargs.get('width', 10)
        padx = kwargs.get('padx', 5)
        pady = kwargs.get('pady', 2)
        minval = kwargs.get('minval', None)
        maxval = kwargs.get('maxval', None)
        comment = kwargs.get('comment', None)
        vfunc = kwargs.get('valid_function', None)
        state = kwargs.get('state', None)

        # first make a frame to contain all elements
        tframe = Tk.Frame(master)
        tframe.pack(side=position, pady=10)

        # set initial value
        labelv = Tk.StringVar(value=desc)
        label = Tk.Label(tframe, textvariable=labelv)
        label.pack(side=Tk.TOP)
        label.config(font=("Times", 12))
        if comment is not None:
            commentv = Tk.StringVar(value=comment)
            comment = Tk.Label(tframe, textvariable=commentv)
            comment.pack(side=Tk.TOP)
            comment.config(font=("Times", 10))

        if fmt == int:
            entryv = Tk.StringVar(value=self.data[name])
        elif fmt == float:
            entryv = Tk.StringVar(value=self.data[name])
        elif fmt == bool:
            entryv = Tk.BooleanVar(value=self.data[name])
        else:
            entryv = Tk.StringVar(value=self.data[name])

        vcmd = [master.register(self.validate)]
        vcmd += self.vcmdargs + [minval, maxval, fmt, name, self.data[name]]
        entry = Tk.Entry(tframe, textvariable=entryv, bg=background,
                         fg=foreground, width=width,
                         validate="all", validatecommand=tuple(vcmd), state=state)
        entry.pack(side=Tk.BOTTOM, padx=padx, pady=pady)
        entry.focus_set()
        self.validation[name] = True
        self.add_entry(name, entry, vfunc)
        self.fmts[name] = fmt
        # self.allbind(entry, self.update)
        return entry

    def add_CheckBox(self, master, name, desc, position, **kwargs):     # added by ronny
        width = kwargs.get('width', 10)
        side = kwargs.get('side', None)
        padx = kwargs.get('padx', 5)
        pady = kwargs.get('pady', 2)
        state = kwargs.get('state', None)
        self.entries[name] = Tk.IntVar()                                # added by ronny
        self.funcs[name] = lambda x: x                                  # added by ronny
        checkbox = Tk.Checkbutton(master, text=desc, variable=self.entries[name], onvalue=True, offvalue=False, state=state)     # added by ronny, variable has to be set before, different to TextEntry
        checkbox.pack(in_=master, side=side, padx=padx, pady=pady)
        self.entries[name].set(self.data[name])                         # added by ronny, set the startvalue
        self.validation[name] = True
        self.fmts[name] = bool
        return checkbox     # added by ronny
        
    def add_button(self, master, text, command, **kwargs):
        state = kwargs.get('state', None)
        background = kwargs.get('background', 'white')
        foreground = kwargs.get('foreground', 'black')
        width = kwargs.get('width', 10)
        side = kwargs.get('side', Tk.BOTTOM)
        padx = kwargs.get('padx', 5)
        pady = kwargs.get('pady', 2)
        
        button = Tk.Button(master=self.master, text=text,
                           command=command, bg=background, fg=foreground,
                           width=width, state=state)
        button.pack(in_=master, side=side, padx=padx, pady=pady)
        return button
    
    def prompt(self, message):
        top = Tk.Toplevel(master=self.master)
        # top.grab_set()
        self.depopulate()
        self.center_window_on_root(top, width=300, height=100)
        top.title('Error')
        self.current_problem = message
        labelv = Tk.StringVar(value=message)
        label = Tk.Label(top, textvariable=labelv)
        label.pack(side=Tk.TOP, padx=30, pady=10)
        button = Tk.Button(master=top, text='Close',
                           command=lambda top=top: self.end_prompt(top))
        button.pack(in_=top, side=Tk.BOTTOM)

    def add_progressbar(self, master, orientation='h'):

        self.progresbar = Tk.Frame(master)
        if 'h' in orientation.lower():
            self.progresbar.pack(side=Tk.BOTTOM, fill='x')
        else:
            self.progresbar.pack(side=Tk.RIGHT, fill='y')

        # set initial value

        lavel1v = Tk.StringVar(value="Status: ")
        label1 = Tk.Label(self.progresbar, textvariable=lavel1v)
        label1.pack(side=Tk.LEFT)
        label1.config(font=("Times", 12, "bold"))
        label1.pack()


        label2 = Tk.Label(self.progresbar, textvariable=self.statusmessage)
        label2.config(font=("Times", 12))
        label2.pack(side=Tk.LEFT)

    # -------------------------------------------------------------------------
    # Events, Updates and Validation
    # -------------------------------------------------------------------------
    def end(self, event=None):          # Run wenn clicking the ExitButton
        """
        Event for clicking the finish button - closes the graph

        :return:
        """
        self.update()   # Update everything (first gui3.data, afterwards kwargs (self.update() calls self.update_plot()) before exiting - Added by Ronny
        
        self.master.quit()
        time.sleep(1)
        self.master.destroy()

    def valueend(self, name):
        """
        Event for clicking the finish button - closes the graph

        :return:
        """
        if name in self.onclickvalue:
            if self.onclickvalue is not None:
                self.data[name] = self.onclickvalue[name]
        self.update(updateplot=False)   # Update everything (first gui3.data, afterwards kwargs (self.update() calls self.update_plot()) before exiting - Added by Ronny
        
        self.master.quit()
        self.master.destroy()

    def allbind(self, entry, update):
        entry.bind('<FocusOut>', update)
        entry.bind('<Leave>', update)
        entry.bind('<Return>', update)

    def add_entry(self, name, entry, vfunc=None):
        self.entries[name] = entry
        if vfunc is not None:
            def proxy(x):
                check, result = vfunc(x)
                if check:
                    return vfunc(x)[1]
                else:
                    self.prompt(result)
                    return None
            self.funcs[name] = lambda x: proxy(x)
        else:
            self.funcs[name] = lambda x: x

    def update(self, event=None, updateplot=True):
        for name in self.entries:
            value = self.entries[name].get()
            F = self.funcs[name]
            if name in self.fmts:
                fmt = self.fmts[name]
                if not self.validation[name]:
                    continue
                if self.validate_fmt(value, fmt):
                    new = F(fmt(value))             # Added by Ronny
                    if new is not None:             # Added by Ronny
                        self.data[name] = new       # Modified by Ronny
        if self.figure is not None and updateplot:  # Added by Ronny
            self.update_plot()

    def update_plot(self):
        # update kwargs with values in self.data
        for name in self.data:
            value = self.data[name]
            if name in self.fmts:
                fmt = self.fmts[name]
                if not self.validation[name]:
                    continue
                if self.data[name] is None:
                    continue
                if self.validate_fmt(value, fmt) and name in self.funkwargs:
                    self.funkwargs[name] = self.data[name]
        # apply function
        if self.func is not None:
            self.update_func()
        if type(self.ax).__name__ in ['numpy.ndarray', 'ndarray']:
            for ii in range(self.ax.shape[0]):
                self.ax[ii].figure.canvas.draw()
        else:
            self.ax.figure.canvas.draw()

    def update_func(self):
        # set working message
        self.master.config(cursor="watch")
        self.statusmessage.set("Loading...")
        self.master.update_idletasks()
        # run function
        self.func(**self.funkwargs)
        # set working message to blank
        self.master.config(cursor="")
        self.statusmessage.set("Ready")
        self.master.update_idletasks()

    def validate_fmt(self, value, fmt):
        try:
            _ = fmt(value)
            return True
        except ValueError:
            return False
        except TypeError:
            return False

    def validate(self, *args):
        # %d = Type of action (1=insert, 0=delete, -1 for others)
        # %i = index of char string to be inserted/deleted, or -1
        # %P = value of the entry if the edit is allowed
        # %s = value of entry prior to editing
        # %S = the text string being inserted or deleted, if any
        # %v = the type of validation that is currently set
        # %V = the type of validation that triggered the callback
        #      (key, focusin, focusout, forced)
        # %W = the tk name of the widget
        # extract params
        d, i, P, s, S, v, V, W, minval, maxval, fmt, name, oldval = args

        # print('Check {0} is a {1}'.format(P, fmt))
        # argnames = ['d', 'i', 'P', 's', 'S', 'v', 'V', 'W', 'minval', 'maxval',
        #             'fmt', 'name', 'oldval']
        # for a, arg in enumerate(args):
        #     print(argnames[a], arg)
        # print('\n\n\n')


        # convert fmt back to a type operator
        if 'float' in fmt:
            fmt, strfmt = float, 'float'
        elif 'int' in fmt:
            fmt, strfmt = int, 'int'
        elif 'bool' in fmt:
            fmt, strfmt = bool, 'bool'
        else:
            fmt, strfmt = str, 'string'
        # convert minval and maxval to floats or None
        minval = None if minval == 'None' else float(minval)
        maxval = None if maxval == 'None' else float(maxval)

        # if v == 'key' or V == 'key':
        #      return True
        if P == '-' or P == '':
            return True

        # Check if P is an int
        if not self.validate_fmt(P, fmt):
            message = '{0} is not a valid {1}'.format(name, strfmt)
            message += '\n Current value = "{0}"'.format(P)
            self.prompt(message)
            self.validation[name] = False
            return False

        P = fmt(P)
        # Check that P is in bounds (if minval and maxval are not None)
        if minval is not None:
            if P <= minval:
                message = '{0} must be greater than {1}'.format(name, minval)
                message += '\n Current value = "{0}"'.format(P)
                self.prompt(message)
                self.validation[name] = False
                return False

        if maxval is not None:
            if P > maxval:
                message = '{0} must be less than {1}'.format(name, maxval)
                message += '\n Current value = "{0}"'.format(P)
                self.prompt(message)
                self.validation[name] = False
                return False

        # else we have validated P
        self.validation[name] = True
        return True

    def end_prompt(self, master):
        master.destroy()
        self.repopulate()
        # master.grab_release()
        
    def prompt_info(self, message, **kwargs):
        width = kwargs.get('width', 400)
        height = kwargs.get('height', 200)
        top = Tk.Toplevel(master=self.master)
        # top.grab_set()
        # Split the message into several lines to comply with width
        charlen = 8.0       # 6.1 Tested 20201030; 8.0 Tested 20210312
        message = message.split('\n')
        for ii,entry in enumerate(message):
            ltxt = len(entry)
            index = 0
            start = int(width/charlen)
            while ltxt*charlen > width:
                pos = entry[start:].find(' ')
                if pos > -1:
                    entry = entry[:start+pos]+'\n'+entry[start+pos+1:]
                    ltxt = len(entry) - (start+pos)
                    start = start+pos + int(width/charlen)
                else:
                    ltxt = 0
            message[ii] = entry
        message = '\n'.join(message)
        height = max( height, (len(message.split('\n'))+5)*12 )     # Tested 20201030

        self.center_window_on_root(top, width=width, height=height)
        top.title('Info')
        labelv = Tk.StringVar(value=message)
        label = Tk.Label(top, textvariable=labelv)
        button = Tk.Button(master=top, text='Close',
                           command=lambda top=top: top.destroy())
        label.pack()
        button.pack()

    # -------------------------------------------------------------------------
    # Other
    # -------------------------------------------------------------------------
    def center_on_screen(self, win, width=None, height=None):
        if width is None:
            width=self.master.winfo_width()
        if height is None:
            height=self.master.winfo_height()
        # get screen size
        x0 = 0
        y0 = 0
        w0 = self.master.winfo_screenwidth()
        h0 = self.master.winfo_screenheight()
        # work out center of self.master
        cx0 = x0 + w0//2
        cy0 = y0 + h0//2
        # work out position of win
        x1 = cx0 - width//2
        y1 = cy0 - height//2
        win.geometry('%dx%d+%d+%d' % (width, height, x1, y1))

    def center_window_on_root(self, win, width=300, height=200):
        x0 = self.master.winfo_rootx()
        y0 = self.master.winfo_rooty()
        w0 = self.master.winfo_width()
        h0 = self.master.winfo_height()
        # work out center of self.master
        cx0 = x0 + w0//2
        cy0 = y0 + h0//2
        # work out position of win
        x1 = cx0 - width//2
        y1 = cy0 - height//2
        win.geometry('%dx%d+%d+%d' % (width, height, x1, y1))
        
# =============================================================================
# Grid class
# ============================================================================= 
class TkCanvasGrid:                                             # Completely new by Ronny
    def __init__(self, figure=None, ax=None, **kwargs):         # "=None" added by Ronny; to use it just as text field
        try:
            screeen_geometry = get_curr_screen_geometry()
        except:
            print('Warn: Could not get the Screen Geometry: This terminal needs to be able to acces a graphical display.')
            screeen_geometry = '1024x768+1+1'
        self.master = Tk.Tk()
        self.screen_w_h = (str(screeen_geometry).replace('-','+').split('+')[0]).split('x')       # list of width and height
        
        # Adopted from https://blog.tecladocode.com/tkinter-scrollable-frames/
        container = ttk.Frame(self.master)
        canvas = Tk.Canvas( container)
        scrollbarx = ttk.Scrollbar(container, orient="horizontal", command=canvas.xview)
        scrollbary = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind("<Configure>",lambda e: canvas.configure( scrollregion=canvas.bbox("all") ) )
        # canvas.bbox("all") gives a 4-value tuple describing the position of two corners of a rectangle which is the scroll region
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(xscrollcommand=scrollbarx.set, yscrollcommand=scrollbary.set)
        
        container.pack(side=Tk.RIGHT)
        canvas.pack(side="left", fill="both", expand=True)
        scrollbarx.pack(side="bottom", fill="x")
        scrollbary.pack(side="right", fill="y")
        
        self.scrollable_frame = scrollable_frame
        
        # extract keyword arguments
        title = kwargs.get('title', 'Matplotlib window')
        self.master.title(title)

        # get figure
        self.figure = figure
        self.ax = ax
        self.func = kwargs.get('func', None)
        self.funkwargs = kwargs.get('kwargs', dict())
        self.widgets = kwargs.get('widgets', dict())
        self.wprops = kwargs.get('widgetprops', dict())
        self.bardir = self.wprops.get('orientation', 'h')

        # default validation command args
        self.vcmdargs = ['%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W']
        # validation variables
        self.current_problem = ""

        # label variables
        self.statusmessage = Tk.StringVar(value="Ready")

        # empty data
        self.data = dict()
        self.entries = dict()
        self.funcs = dict()
        self.validation = dict()
        self.fmts = dict()
        self.onclickvalue = dict()

        # populate matplotlib drawing area and buttons
        self.canvas = canvas    # To allow regular updating and the use of tqdm
        self.populate()
        #print(time.time(),'before canvas.update()')
        canvas.update()     # adds the winfo information
        #print(time.time(),'after canvas.update()')
        self.width_GUI  = max(300, min( self.wprops.get('width_GUI' , 1E9), int(self.screen_w_h[0])-50, scrollable_frame.winfo_width() ) )
        self.height_GUI = max(200, min( self.wprops.get('height_GUI', 1E9), int(self.screen_w_h[1])-90, scrollable_frame.winfo_height() ) )
        canvas.config(width=self.width_GUI, height=self.height_GUI)
        canvas.update()     # adds the winfo information

        #print(self.master.winfo_height(), self.master.winfo_width(), self.master.winfo_screenwidth(), self.master.winfo_screenheight(), self.master.winfo_geometry(), self.master.geometry())   # all of these only show a size of 1x1
        if self.wprops.get('fullscreen', False):
            self.master.geometry(screeen_geometry)


    def populate(self):
        if self.figure is not None:                  # Added by Ronny
            self.matplotlib_window()
            self.add_progressbar(self.master, orientation='h')
        self.add_widgets(parent=self.scrollable_frame)

    def depopulate(self):
        for widget in self.ws.keys():
            self.ws[widget].config(state='disable')

    def repopulate(self):
        for widget in self.ws.keys():
            self.ws[widget].config(state="normal")

    def matplotlib_window(self):
        self.mwindow = FigureCanvasTkAgg(self.figure, master=self.master)
        self.mwindow.get_tk_widget().config(bg='white')
        # add matplotlib toolbar
        toolbar = NavigationToolbar2TkAgg(self.mwindow, self.master)
        toolbar.update()
        self.mwindow.get_tk_widget().pack(side=Tk.TOP, fill='both',
                                          expand='YES')
        #self.mwindow.get_tk_widget().grid()        # This is a pack subframe
        try:
            # fix native matplotlib cursor bug  # fails in version 2.2
            tkagg.cursord[cursors.POINTER] = ""
        except:
            "do nothing"
        self.mwindow.draw()

    def add_widgets(self, parent):
        widgets = self.widgets
        self.ws = dict()
        
        window_size = [0, 0]
        
        # Used tqdm, if installed and if the number of data is worth it to use
        def no_tqdm(input, desc=''):
            """
            In order to switch between tqdm on and off, this is needed
            """
            return input
        try:
            from tqdm import tqdm
        except:
            tqdm = no_tqdm
        if len(widgets) > 1500:
            plot_tqdm = tqdm
        else:
            plot_tqdm = no_tqdm
            
        for ii,widget in enumerate(plot_tqdm(widgets, desc='The statusbar might stop updating for 60s')):
            #print(time.time(),widget)
            kind = widgets[widget].get('kind', 'TextEntry')
            name = widget
            label = widgets[widget].get('label', 'Enter {0}'.format(name))
            start = widgets[widget].get('start', 0)
            fmt = widgets[widget].get('fmt', int)
            minval = widgets[widget].get('minval', None)
            maxval = widgets[widget].get('maxval', None)
            comment = widgets[widget].get('comment', None)
            result = widgets[widget].get('result', None)
            onclick = widgets[widget].get('onclick', None)
            vfunc = widgets[widget].get('valid_function', None)
            width = widgets[widget].get('width', 10)
            row = widgets[widget].get('row', None)
            column = widgets[widget].get('column', None)
            rowspan = widgets[widget].get('rowspan', 1)
            columnspan = widgets[widget].get('columnspan', 1)
            orientation = widgets[widget].get('orientation', None)
            wraplength = widgets[widget].get('wraplength', None)
            state = widgets[widget].get('state', None)
            command = widgets[widget].get('command', None)
            
            kwargs = dict(width=width, minval=minval, maxval=maxval, valid_function=vfunc, comment=comment, state=state)

            if kind == 'Label':
                w = Tk.Label(parent, text=label, wraplength=wraplength)
            elif kind == 'TextEntry':
                if name not in self.data:
                    self.data[name] = start
                w = self.add_TextEntry(parent, name=name, fmt=fmt, command=command, **kwargs)
            elif kind == 'CheckBox':
                if name not in self.data:
                    self.data[name] = start
                w = self.add_CheckBox(parent, name=name, desc=label, **kwargs)
            elif kind == 'ExitButton':
                if result is not None:
                    self.data[name] = result
                    self.onclickvalue[name] = onclick
                if onclick is None:
                    w = self.add_button(parent, label, self.end, **kwargs)
                else:
                    w = self.add_button(parent, label,
                                        command=lambda name=name:
                                                self.valueend(name), **kwargs)
            elif kind == 'UpdatePlot':
                w = self.add_button(parent, label, self.update, **kwargs)
                if result is not None:
                    self.data[name] = result
                    self.fmts[name] = type(result)
            elif kind == 'CommandButton':
                w = self.add_button(parent, label, command, **kwargs)
                if result is not None:
                    self.data[name] = result
                    self.fmts[name] = type(result)
            else:
                w = None
            # append to storage
            if w is not None:
                #if row is not None and column is not None:
                w.grid(row=row, column=column, rowspan=rowspan, columnspan=columnspan, sticky=orientation)
                #print('row,column',row,column)
                self.ws[name] = w
                if len(widgets) > 1500:
                    if np.remainder(ii, 128) == 0:
                        self.canvas.update()

    #---------------------------------------------------------------------------
    #define sub widgets
    #---------------------------------------------------------------------------
    def add_TextEntry(self, master, name, fmt, **kwargs):

        background = kwargs.get('background', 'white')
        foreground = kwargs.get('foreground', 'black')
        width = kwargs.get('width', 10)
        minval = kwargs.get('minval', None)
        maxval = kwargs.get('maxval', None)
        vfunc = kwargs.get('valid_function', None)
        state = kwargs.get('state', None)
        command = kwargs.get('command', None)

        if fmt == int:
            entryv = Tk.StringVar(value=self.data[name])
        elif fmt == float:
            entryv = Tk.StringVar(value=self.data[name])
        elif fmt == bool:
            entryv = Tk.BooleanVar(value=self.data[name])
        else:
            entryv = Tk.StringVar(value=self.data[name])

        vcmd = [master.register(self.validate)]
        vcmd += self.vcmdargs + [minval, maxval, fmt, name, self.data[name]]
        entry = Tk.Entry(master, textvariable=entryv, bg=background,
                         fg=foreground, width=width,
                         validate="all", validatecommand=tuple(vcmd), state=state)
        entry.focus_set()
        self.validation[name] = True
        self.add_entry(name, entry, vfunc)
        self.fmts[name] = fmt
        # self.allbind(entry, self.update)
        return entry

    def add_CheckBox(self, master, name, desc, **kwargs):     # added by ronny
        state = kwargs.get('state', None)
        self.entries[name] = Tk.IntVar()                                # added by ronny
        self.funcs[name] = lambda x: x                                  # added by ronny
        checkbox = Tk.Checkbutton(master, text=desc, variable=self.entries[name], onvalue=True, offvalue=False, state=state)     # added by ronny, variable has to be set before, different to TextEntry
        self.entries[name].set(self.data[name])                         # added by ronny, set the startvalue
        self.validation[name] = True
        self.fmts[name] = bool
        return checkbox     # added by ronny
        
    def add_button(self, master, text, command, **kwargs):

        background = kwargs.get('background', 'white')
        foreground = kwargs.get('foreground', 'black')
        width = kwargs.get('width', 10)
        
        button = Tk.Button(master=master, text=text,
                           command=command, bg=background, fg=foreground,
                           width=width)
        return button
    
    def prompt(self, message):
        top = Tk.Toplevel(master=self.master)
        # top.grab_set()
        self.depopulate()
        self.center_window_on_root(top, width=300, height=100)
        top.title('Error')
        self.current_problem = message
        labelv = Tk.StringVar(value=message)
        label = Tk.Label(top, textvariable=labelv)
        button = Tk.Button(master=top, text='Close',
                           command=lambda top=top: self.end_prompt(top))
        label.pack()
        button.pack()
        
    def prompt_info(self, message, **kwargs):
        width = kwargs.get('width', 400)
        height = kwargs.get('height', 200)
        top = Tk.Toplevel(master=self.master)
        # top.grab_set()
        # Split the message into several lines to comply with width
        charlen = 8.0       # 6.1 Tested 20201030; 8.0 Tested 20210312
        message = message.split('\n')
        for ii,entry in enumerate(message):
            ltxt = len(entry)
            index = 0
            start = int(width/charlen)
            while ltxt*charlen > width:
                pos = entry[start:].find(' ')
                if pos > -1:
                    entry = entry[:start+pos]+'\n'+entry[start+pos+1:]
                    ltxt = len(entry) - (start+pos)
                    start = start+pos + int(width/charlen)
                else:
                    ltxt = 0
            message[ii] = entry
        message = '\n'.join(message)
        height = max( height, (len(message.split('\n'))+5)*12 )     # Tested 20201030

        self.center_window_on_root(top, width=width, height=height)
        top.title('Info')
        labelv = Tk.StringVar(value=message)
        label = Tk.Label(top, textvariable=labelv)
        button = Tk.Button(master=top, text='Close',
                           command=lambda top=top: top.destroy())
        label.pack()
        button.pack()

    def add_progressbar(self, master, orientation='h'):

        self.progresbar = Tk.Frame(master)
        if 'h' in orientation.lower():
            self.progresbar.pack(side=Tk.BOTTOM, fill='x')      # This is a pack subframe
        else:
            self.progresbar.pack(side=Tk.RIGHT, fill='y')
       
        # set initial value
        lavel1v = Tk.StringVar(value="Status: ")
        label1 = Tk.Label(self.progresbar, textvariable=lavel1v)
        label1.grid()
        label1.config(font=("Times", 12, "bold"))
        label1.grid()


        label2 = Tk.Label(self.progresbar, textvariable=self.statusmessage)
        label2.config(font=("Times", 12))
        label2.grid()

    # -------------------------------------------------------------------------
    # Events, Updates and Validation
    # -------------------------------------------------------------------------
    def end(self, event=None):          # Run wenn clicking the ExitButton
        """
        Event for clicking the finish button - closes the graph

        :return:
        """
        self.update()   # Update everything (first gui3.data, afterwards kwargs (self.update() calls self.update_plot()) before exiting - Added by Ronny
        
        self.master.quit()
        time.sleep(1)
        self.master.destroy()

    def valueend(self, name):
        """
        Event for clicking the finish button - closes the graph

        :return:
        """
        if name in self.onclickvalue:
            if self.onclickvalue is not None:
                self.data[name] = self.onclickvalue[name]
        self.update(updateplot=False)   # Update everything (first gui3.data, afterwards kwargs (self.update() calls self.update_plot()) before exiting - Added by Ronny
        
        self.master.quit()
        time.sleep(1)
        self.master.destroy()

    def allbind(self, entry, update):
        entry.bind('<FocusOut>', update)
        entry.bind('<Leave>', update)
        entry.bind('<Return>', update)

    def add_entry(self, name, entry, vfunc=None):
        self.entries[name] = entry
        if vfunc is not None:
            def proxy(x):
                check, result = vfunc(x)
                if check:
                    return vfunc(x)[1]
                else:
                    self.prompt(result)
                    return None
            self.funcs[name] = lambda x: proxy(x)
        else:
            self.funcs[name] = lambda x: x

    def update(self, event=None, updateplot=True):
        # update kwargs with values in self.data
        for name in self.entries:
            value = self.entries[name].get()
            F = self.funcs[name]
            if name in self.fmts:
                fmt = self.fmts[name]
                if not self.validation[name]:
                    continue
                if self.validate_fmt(value, fmt):
                    new = F(fmt(value))             # Added by Ronny
                    if new is not None:             # Added by Ronny
                        self.data[name] = new       # Modified by Ronny
                        self.funkwargs[name] = new  # Moved this from update_plot to here
            #print('update', name, value, self.data[name])
            #if name in self.fmts: print('fmts', self.fmts[name], self.validation[name], self.validate_fmt(value, fmt))
            #if name in self.fmts and name in self.funkwargs: print('funkawgs', self.funkwargs[name])
            
        if self.figure is not None and updateplot:  # Added by Ronny
            self.update_plot()

    def update_plot(self):
        """# update kwargs with values in self.data
        for name in self.data:
            value = self.data[name]
            if name in self.fmts:
                fmt = self.fmts[name]
                if not self.validation[name]:
                    continue
                if self.data[name] is None:
                    continue
                if self.validate_fmt(value, fmt) and name in self.funkwargs:
                    self.funkwargs[name] = self.data[name]
            print('update_plot', name, value, self.data[name])
            if name in self.fmts: print('fmts', self.fmts[name], self.validation[name], self.validate_fmt(value, fmt))
            if name in self.fmts and name in self.funkwargs: print('funkawgs', self.funkwargs[name])"""
        # apply function
        if self.func is not None:
            self.update_func()
        self.figure.canvas.draw()               # The GUI might hang on this step if another GUI is open/something is messed up? However, only have the widgets have been changed

    def update_func(self):
        # set working message
        self.master.config(cursor="watch")
        self.statusmessage.set("Loading...")
        #self.master.update_idletasks()          # The GUI might hang on this step if another GUI is open?
        # run function
        self.func(**self.funkwargs)
        # set working message to blank
        self.master.config(cursor="")
        self.statusmessage.set("Ready")
        #self.master.update_idletasks()

    def validate_fmt(self, value, fmt):
        try:
            _ = fmt(value)
            return True
        except ValueError:
            return False
        except TypeError:
            return False

    def validate(self, *args):
        # %d = Type of action (1=insert, 0=delete, -1 for others)
        # %i = index of char string to be inserted/deleted, or -1
        # %P = value of the entry if the edit is allowed
        # %s = value of entry prior to editing
        # %S = the text string being inserted or deleted, if any
        # %v = the type of validation that is currently set
        # %V = the type of validation that triggered the callback
        #      (key, focusin, focusout, forced)
        # %W = the tk name of the widget
        # extract params
        d, i, P, s, S, v, V, W, minval, maxval, fmt, name, oldval = args

        # print('Check {0} is a {1}'.format(P, fmt))
        # argnames = ['d', 'i', 'P', 's', 'S', 'v', 'V', 'W', 'minval', 'maxval',
        #             'fmt', 'name', 'oldval']
        # for a, arg in enumerate(args):
        #     print(argnames[a], arg)
        # print('\n\n\n')


        # convert fmt back to a type operator
        if 'float' in fmt:
            fmt, strfmt = float, 'float'
        elif 'int' in fmt:
            fmt, strfmt = int, 'int'
        elif 'bool' in fmt:
            fmt, strfmt = bool, 'bool'
        else:
            fmt, strfmt = str, 'string'
        # convert minval and maxval to floats or None
        minval = None if minval == 'None' else float(minval)
        maxval = None if maxval == 'None' else float(maxval)

        # if v == 'key' or V == 'key':
        #      return True
        if P == '-' or P == '':
            return True

        # Check if P is an int
        if not self.validate_fmt(P, fmt):
            message = '{0} is not a valid {1}'.format(name, strfmt)
            message += '\n Current value = "{0}"'.format(P)
            self.prompt(message)
            self.validation[name] = False
            return False

        P = fmt(P)
        # Check that P is in bounds (if minval and maxval are not None)
        if minval is not None:
            if P <= minval:
                message = '{0} must be greater than {1}'.format(name, minval)
                message += '\n Current value = "{0}"'.format(P)
                self.prompt(message)
                self.validation[name] = False
                return False

        if maxval is not None:
            if P > maxval:
                message = '{0} must be less than {1}'.format(name, maxval)
                message += '\n Current value = "{0}"'.format(P)
                self.prompt(message)
                self.validation[name] = False
                return False

        # else we have validated P
        self.validation[name] = True
        return True

    def end_prompt(self, master):
        master.destroy()
        self.repopulate()
        # master.grab_release()

    # -------------------------------------------------------------------------
    # Other
    # -------------------------------------------------------------------------
    def center_on_screen(self, win, width=None, height=None):
        if width is None:
            width=self.master.winfo_width()
        if height is None:
            height=self.master.winfo_height()
        # get screen size
        x0 = 0
        y0 = 0
        w0 = self.master.winfo_screenwidth()
        h0 = self.master.winfo_screenheight()
        # work out center of self.master
        cx0 = x0 + w0//2
        cy0 = y0 + h0//2
        # work out position of win
        x1 = cx0 - width//2
        y1 = cy0 - height//2
        win.geometry('%dx%d+%d+%d' % (width, height, x1, y1))

    def center_window_on_root(self, win, width=300, height=200):
        x0 = self.master.winfo_rootx()
        y0 = self.master.winfo_rooty()
        w0 = self.master.winfo_width()
        h0 = self.master.winfo_height()
        # work out center of self.master
        cx0 = x0 + w0//2
        cy0 = y0 + h0//2
        # work out position of win
        x1 = cx0 - width//2
        y1 = cy0 - height//2
        win.geometry('%dx%d+%d+%d' % (width, height, x1, y1))


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":

    fig, frame = plt.subplots(nrows=1, ncols=1)
    xlow, xhigh = 0.0, 3.0
    testronny, testronny2 = False, True

    def plot(ax, xlow, xhigh, testronny, testronny2):
        ax.clear()
        x = np.arange(xlow, xhigh, 0.01)
        y = np.sin(2 * np.pi * x)
        ax.plot(x, y)
        ax.set(xlabel='x', ylabel='y')
        print(testronny,testronny2)

    pkwargs = dict(ax=frame, xlow=xlow, xhigh=xhigh, testronny=testronny, testronny2=testronny2)

    plot(**pkwargs)


    widgets = dict()
    widgets['xlow'] = dict(label='Enter x low', comment='-20 < x low < 20',
                           kind='TextEntry',
                           minval=-20, maxval=+20, fmt=float, start=xlow)
    widgets['testronny'] = dict(label='Test for Checkbox',
                           kind='CheckBox', start=testronny)
    widgets['xhigh'] = dict(label='Enter x high', comment='-20 < x low < 20',
                            kind='TextEntry',
                            minval=-20, maxval=+20, fmt=float, start=xhigh)
    widgets['testronny2'] = dict(label='Test for Checkbox2',
                           kind='CheckBox', start=testronny2)
    widgets['close'] = dict(label='Next', kind='ExitButton',
                            position=Tk.BOTTOM, result='Close')
    widgets['update'] = dict(label='Update', kind='UpdatePlot',
                             position=Tk.BOTTOM)
    

    wprops = dict(orientation='v')

    gui = TkCanvas(figure=fig, ax=frame, func=plot, kwargs=pkwargs,
                   title='sin test run', widgets=widgets, widgetprops=wprops)
    # root.eval('tk::PlaceWindow {0} center'.format(root.winfo_pathname(root.winfo_id())))
    gui.master.mainloop()

# =============================================================================
# End of code
# =============================================================================
