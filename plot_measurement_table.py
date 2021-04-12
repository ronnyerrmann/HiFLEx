# To use this file, best copy it to your result folder and run it there
import sys
sys.path.append('/home/ronny/Scripts/hiflex/')
from procedures import *

params = dict()

savefile = 'measurement_table.pdf'
datafile = 'measurement_table.csv'      # assumes that the objectnames are somewhat in column one or two
objectfile = 'objects_reduced.lst'
x_range = None
size_inches = [16.2, 10]
objects = ['*']                     # Create a plot with everything
objects += read_text_file(objectfile, no_empty_lines=False)
if len(objects) == 1:
    # Code defined what objects to use
    objects = ['*', 'Sun', ['alpOri','alphori','alfori','alphaori'], 'AlphaAri', 'alpcmi', ['gamVirA','gamVir2'], ['gamVirB','gamVir1'], 'TauBoo', 'Tung' ]

if datafile.find(os.sep) == -1:
    datafile = os.getcwd() + os.sep + datafile      # To select settings depending on datafile

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

data = read_text_file(datafile, no_empty_lines=False)
len_entries = len( data[0].split('\t') )
data = convert_readfile(data, [str]*len_entries, delimiter='\t', replaces=['\n'])
data = np.array(data).T
header = data[:,0:3]
data = data[:,3:]
data[data==''] = 'nan'
converted = dict()

with PdfPages(savefile) as pdf:
    for star in objects:        # New page for each opject
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
                if selectd[ii_d][0] == 'data':
                    [indexx, indexy, indexe, label, color, linestyle, marker, markersize ] = selectd[ii_d][1:9]
                    if indexy < 0:
                        indexy = abs(indexy)
                        sign = -1
                    else:
                        sign = 1
                    for index in [indexx, indexy, indexe]:
                        if index is None:   continue
                        if '{0}'.format(index) not in converted:
                            converted['{0}'.format(index)] = data[index,:].astype(float)
                    #x_title = header[1,indexx]+' '+header[2,indexx]
                    #if x_title not in x_titles:
                    #    x_titles.append(x_title)
                    if label == '':
                        label = header[indexy, 1]
                        if len(label) > 25:
                            for ii in list(range(25, 18, -1)) + list(range(25,min(len(label),35))):
                                if label[ii] == ' ':
                                    label = label[:ii] + '\n' + label[ii+1:]
                                    break
                    good_data = ~np.isnan(converted['{0}'.format(indexx)]) & ~np.isnan(converted['{0}'.format(indexy)]) & starindata
                    if np.sum(good_data) > 0:       # Only plot if data is available
                        xx = converted['{0}'.format(indexx)][good_data]
                        yy = sign*converted['{0}'.format(indexy)][good_data]
                        if y_range is not None:
                            yy[ yy > y_range[1] ] = y_range[1]
                            yy[ yy < y_range[0] ] = y_range[0]
                        if indexe is None:          # No error bars
                            frame_ii_g.plot(xx, yy, label=label, color=color, linestyle=linestyle, marker=marker, markersize=markersize)
                        else:                       # error bars
                            frame_ii_g.errorbar(xx, yy, yerr=converted['{0}'.format(indexe)][good_data], label=label, color=color, linestyle=linestyle, marker=marker, markersize=markersize)
                    if len(selectd[ii_d]) == 11:
                        [x_title, y_title] = selectd[ii_d][9:11]
                elif selectd[ii_d][0] == 'text':
                    [posx, posy, text, color, rotation, hpos, vpos ] = selectd[ii_d][1:8]
                    frame_ii_g.text( posx, posy, text, horizontalalignment=hpos, verticalalignment=vpos, rotation=rotation, color=color, zorder=5 )
                elif selectd[ii_d][0] == 'y_range':
                    y_range = selectd[ii_d][1:3]
            frame_ii_g.legend(loc='upper left', bbox_to_anchor=(    1.0,1.01), fontsize=11, handletextpad=0.2)
            frame_ii_g.set_xlabel(x_title, fontsize=11)
            frame_ii_g.set_ylabel(y_title, fontsize=11)
            if y_range is not None:
                frame_ii_g.set_ylim(y_range[0],y_range[1])
        if x_range is not None:
            frame_ii_g.set_xlim(x_range[0],x_range[1])
        
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()






