# To use this file, best copy it to your result folder and run it there
import sys
sys.path.append('/home/ronny/Scripts/hiflex/')
sys.path.append('/hiflex/')
from procedures import *
import csv

params = dict()

savefile = 'measurement_table.pdf'
datafile = 'measurement_table.csv'      # assumes that the objectnames are somewhat in column one or two
sensorsfile = 'TSP01_20210624-1822.csv'      # optional
objectfile = 'objects_reduced.lst'
x_range = None
size_inches = [16.2, 10]
objects = ['*']                     # Create a plot with everything
objects += read_text_file(objectfile, no_empty_lines=False)
if len(objects) == 1:
    # Code defined what objects to use
    objects = ['*', 'Sun', ['alpOri','alphori','alfori','alphaori'], 'AlphaAri', 'alpcmi', ['gamVirA','gamVir2'], ['gamVirB','gamVir1'], 'TauBoo', 'Tung' ]

select = []     
if True:
    selectd = []    # 'data'/'text', index x, index y, error index, label (empty to use header), color, linestyle, marker, markersize, [x_title, y_title]
    #selectd.append([ 'y_range', -0.7, 0.7])
    selectd.append([ 'data', 5, 8, None, '', 'g', '', 'o', 2, '', 'Offset [px]' ])
    selectd.append([ 'data', 5, 10, 11, '', 'b', '', 'o', 2, '', 'Offset [px]' ])
    selectd.append([ 'data', 5, 15, None, '', 'k', '', 'o', 2, '', 'Offset [px]' ])
    select.append(selectd)
    selectd = []
    selectd.append([ 'data', 5, 9, None, '', 'g', '', 'o', 2, '', 'Width [px]' ])
    selectd.append([ 'data', 5, 13, None, '', 'b', '', 'o', 2, '', 'Width [px]' ])
    #selectd.append([ 'y_range', 0.5, 1.5])
    select.append(selectd)
    selectd = []
    selectd.append([ 'data', 5, 12, None, '', 'k', '', 'o', 2, '', 'Number' ])
    select.append(selectd)
    selectd = []
    selectd.append([ 'sensordata', 0, 1, None, 'Behind camera (water)', 'tab:blue', '', 'o', 2, '', '' ])
    #selectd.append([ 'sensordata', 0, 2, None, 'Room - front right', 'tab:orange', '', 'o', 2, '', '' ])
    selectd.append([ 'sensordata', 0, 3, None, 'Breadboard - next to grating holder', 'tab:green', '', 'o', 2, '', '' ])
    #selectd.append([ 'sensordata', 0, 5, None, 'Spectrograph - low between fiber and AO', 'tab:red', '', 'o', 2, '', '' ])
    selectd.append([ 'sensordata', 0, 6, None, 'Spectrograph - above fiber', 'tab:purple', '', 'o', 2, '', '' ])
    #selectd.append([ 'sensordata', 0, 7, None, 'Spectrograph - front of grating', 'tab:brown', '', 'o', 2, '', '' ])
    #selectd.append([ 'sensordata', 0, 9, None, 'On Breadboard - behind camera', 'tab:pink', '', 'o', 2, '', '' ])
    selectd.append([ 'sensordata', 0, 10, None, 'Room/Water cooling', 'tab:gray', '', 'o', 2, '', '' ])
    #selectd.append([ 'sensordata', 0, 11, None, 'Breadboard - next to grating holder', 'tab:olive', '', 'o', 2, 'JD [d]', 'Temperature' ])
    selectd.append([ 'text', 2459342.011, 23, 'Helium off', 'k', 90, 'left', 'bottom' ])
    selectd.append([ 'grid' ])
    #selectd.append([ 'y_range', 20, 25])
    #selectd.append([ 'text', 294, 40, 'Helium 1 l/min', 'k', 90, 'left', 'bottom' ])
    select.append(selectd)
    selectd = []
    selectd.append([ 'data', 5, 6, None, '', 'k', '', 'o', 2, '', 'Flux (ADU)' ])
    select.append(selectd)
    selectd = []
    #selectd.append([ 'y_range', -20, 70])
    #selectd.append([ 'data', 17, 18, 20, '', 'grey', '', 'o', 2, '', 'RV [km/s]' ])
    selectd.append([ 'data', 17, 19, 20, '', 'k', '', 'o', 2, '', 'RV [km/s]' ])
    select.append(selectd)
    selectd = []
    selectd.append([ 'data', 17, 24, 25, '', 'r', '', 'o', 2, '', 'RV [m/s]' ])
    selectd.append([ 'data', 17, 26, 27, '', 'c', '', 'o', 2, 'JD/BJD [d]', 'RV [m/s]' ])
    #selectd.append([ 'text', 294, 40, 'Helium 1 l/min', 'k', 90, 'left', 'bottom' ])
    #selectd.append([ 'y_range', -800, 400])
    select.append(selectd)
    #x_range = [-10, 470]
    #size_inches = [10, 6]

if datafile.find(os.sep) == -1:
    datafile = os.getcwd() + os.sep + datafile      # To select settings depending on datafile

if os.path.isfile(sensorsfile):
    with open(sensorsfile, 'r') as f:
        header = list(csv.reader(f, delimiter=','))
        data = np.array(header[1:])
        times1 = [ get_julian_datetime( datetime.datetime.strptime(dat, '%Y%m%d-%H%M%S') ) for dat in data[:,0]]
        times2 = [ get_julian_datetime( datetime.datetime.strptime(dat, '%Y%m%d-%H%M%S') ) for dat in data[:,-1]]
        sensordata = np.vstack((times1, data[:,1:-1].astype(float).T, times2)).T
else:
    sensordata = np.zeros((10,100))*np.nan

data = read_text_file(datafile, no_empty_lines=False)
if len(data) == 0:
    print('Error: no datapoints found')
    exit()
len_entries = len( data[0].split('\t') )
data = convert_readfile(data, [str]*len_entries, delimiter='\t', replaces=['\n', os.linesep])
data = np.array(data).T
header = data[:,0:3]
data = data[:,3:]
data[data==''] = 'nan'
converted = dict()

with PdfPages(savefile) as pdf:
    for star in objects:        # New page for each opject
        x_range_data = [1E10, -1E10]
        nr_g = len(select)      # How many graphs
        if type(star).__name__ in ['str']:
            star = [star]
        starindata = ( data[1,:] == star[0] )
        for entry in copy.copy(star):
            if entry == '*':
                starindata[:] = True
                break
            starindata_e1 = [ (ii.lower().find(entry.lower()) > -1) for ii in data[0,:] ]
            starindata_e2 = [ (ii.lower().find(entry.lower()) > -1) for ii in data[1,:] ]
            starindata = (starindata + starindata_e1 + starindata_e2) > 0
        if not np.any(starindata):      # Only plot objects, that are in the data
            continue
        fig, frame = plt.subplots(nr_g,1, sharex=True)
        fig.set_size_inches(size_inches[0], size_inches[1])
        plt.subplots_adjust(left=0.05, right=0.85, top=0.95, bottom=0.14)
        
        x_titles, y_titles = [], []
        for ii_g in range(nr_g):
            selectd = select[ii_g]
            nr_d = len(selectd) # How many data sets
            if nr_g == 1:
                frame_ii_g = frame
            else:
                frame_ii_g = frame[ii_g]
            if ii_g == 0:
                frame_ii_g.set_title('{0}'.format(star))
            x_title, y_title, y_range = '', '', None
            for ii_d in range(nr_d):
                if selectd[ii_d][0] in ['data','sensordata']:
                    [dset, indexx, indexy, indexe, label, color, linestyle, marker, markersize ] = selectd[ii_d][0:9]
                    if dset == 'sensordata':
                        x_range_tmp = copy.copy(frame_ii_g.get_xlim())
                    if indexy < 0:
                        indexy = abs(indexy)
                        sign = -1
                    else:
                        sign = 1
                    for index in [indexx, indexy, indexe]:
                        if index is None:   continue
                        if '{0}_{1}'.format(index, dset) not in converted.keys():
                            if dset == 'data':
                                converted['{0}_{1}'.format(index, dset)] = data[index,:].astype(float)
                            elif dset == 'sensordata':
                                converted['{0}_{1}'.format(index, dset)] = sensordata[:,index]
                    #x_title = header[1,indexx]+' '+header[2,indexx]
                    #if x_title not in x_titles:
                    #    x_titles.append(x_title)
                    if label == '' and dset == 'data':
                        label = header[indexy, 1]
                        if len(label) > 25:
                            for ii in list(range(25, 18, -1)) + list(range(25,min(len(label),35))):
                                if label[ii] == ' ':
                                    label = label[:ii] + os.linesep + label[ii+1:]
                                    break
                    if dset == 'data':
                        good_data = ~np.isnan(converted['{0}_{1}'.format(indexx, dset)]) & ~np.isnan(converted['{0}_{1}'.format(indexy, dset)]) & starindata
                    else:
                        good_data = np.ones_like(converted['{0}_{1}'.format(indexx, dset)], dtype=bool)
                    if np.sum(good_data) > 0:       # Only plot if data is available
                        xx = converted['{0}_{1}'.format(indexx, dset)][good_data]
                        yy = sign*converted['{0}_{1}'.format(indexy, dset)][good_data]
                        if y_range is not None:
                            yy[ yy > y_range[1] ] = y_range[1]
                            yy[ yy < y_range[0] ] = y_range[0]
                        if dset == 'data':
                            x_range_data = [ min(np.min(xx), x_range_data[0]), max(np.max(xx), x_range_data[1]) ]
                        if indexe is None:          # No error bars
                            frame_ii_g.plot(xx, yy, label=label, color=color, linestyle=linestyle, marker=marker, markersize=markersize)
                        else:                       # error bars
                            frame_ii_g.errorbar(xx, yy, yerr=converted['{0}_{1}'.format(indexe, dset)][good_data], label=label, color=color, linestyle=linestyle, marker=marker, markersize=markersize)
                    if len(selectd[ii_d]) == 11:
                        if len(selectd[ii_d][9]) > 0: x_title = selectd[ii_d][9]
                        if len(selectd[ii_d][10]) > 0: y_title = selectd[ii_d][10]
                    if dset == 'sensordata':
                        frame_ii_g.set_xlim(x_range_tmp[0],x_range_tmp[1])
                elif selectd[ii_d][0] == 'text':
                    [posx, posy, text, color, rotation, hpos, vpos ] = selectd[ii_d][1:8]
                    frame_ii_g.text( posx, posy, text, horizontalalignment=hpos, verticalalignment=vpos, rotation=rotation, color=color, zorder=5 )
                elif selectd[ii_d][0] == 'y_range':
                    y_range = selectd[ii_d][1:3]
                elif selectd[ii_d][0] == 'grid':
                    frame_ii_g.grid(which='both', axis='both')
            frame_ii_g.legend(loc='upper left', bbox_to_anchor=(    1.0,1.01), fontsize=11, handletextpad=0.2)
            frame_ii_g.set_xlabel(x_title, fontsize=11)
            frame_ii_g.set_ylabel(y_title, fontsize=11)
            frame_ii_g.yaxis.set_ticks_position('both')
            if y_range is not None:
                frame_ii_g.set_ylim(y_range[0],y_range[1])
        if x_range is not None:
            frame_ii_g.set_xlim(x_range[0],x_range[1])
        elif x_range_data[0] < 1E9 and x_range_data[1] > -1E9:
            b = 0.01 * (x_range_data[1] - x_range_data[0])
            frame_ii_g.set_xlim(x_range_data[0]-b,x_range_data[1]+b)
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()






