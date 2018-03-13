import numpy as np
import os
from procedures import *

# =============================================================================
# Define variables
# =============================================================================
params = dict()     # default param dictionary
params['path_extraction'] = 'converted/'

if __name__ == "__main__":
    
    # Prepare folder
    entry = 'path_extraction'
    if not os.path.exists(params[entry]):
        try:
            os.makedirs(params[entry])
        except:
            logger('Warn: Cant create directory {0}'.format(params[entry]))
    # Read the file which contains all the file to be converted
    fitsfiles = []
    if not os.path.isfile('convert_files.lst'):
        logger('Error: convert_files.lst is missing')
    files = read_text_file('convert_files.lst')
    for line in files:
        if line[0] <> '#':
            fname = line.rsplit('/',1)
            fitsfiles.append([line, fname[1]])
    # Perform the convertion
    for fname in fitsfiles:
        if not os.path.isfile(fname[0]):
            logger('Error: {0} is missing'.format(fname[0]))
        hdul = fits.open(fname[0])     # Read the file
        im_head = hdul[0].header
        spectra = np.array(hdul[0].data)
        hdul.close()
        spectras = spectra.shape
        # Get the arc spectrum
        aspecname = fname[0].replace('e2ds_A.fits','e2ds_B.fits')
        if os.path.isfile(aspecname):
            aspectra = np.array(fits.getdata(aspecname))
            aspectras = aspectra.shape
            if aspectras[0] == spectras[0] - 1:                # Order 45 (starting at 0) is missing in the data
                aspectra = np.insert(aspectra, 45, 0, axis=0)   # Fills line 45 with 0
            elif aspectras <> spectras:
                logger('Error: The size of the science spectrum ({0}: {2}) and the calibration spectrum ({1}: {3}) do not match.'.format(fname[0], fname[0].replace('e2ds_A.fits','e2ds_B.fits'), spectras, aspectras))
        else:
            aspectra = np.zeros(spectras)
        # Create the Wavelength Array
        deg = im_head['HIERARCH ESO DRS CAL TH DEG LL']             # Degrees of Freedom
        x0 = np.arange(spectras[1])                                  # x is probably 0 based
        #x1 = np.arange(spectras[1]) + 1                                 # x is probably 0 based
        wavelengths = np.zeros(spectras)                                    # Array for the wavelengths
        for order in range(spectras[0]):
            wave0 = np.zeros(spectras[1])
            #wave1 = np.zeros(spectras[1])
            for i in range(deg+1):
                A = im_head['HIERARCH ESO DRS CAL TH COEFF LL{0}'.format( i + order * ( deg + 1 ) )]
                wave0 +=  A * x0**i
                #wave1 +=  A * x1**i
            #print min(wave0), min(wave1), max(wave0), max(wave1), max(np.abs(wave0[1:]-wave1[:-1]))     # No difference between wave0 and wave1
            wavelengths[order,:] = wave0
        # Flip up/down to have the reddest orders first:
        spectra = np.flipud(spectra)
        aspectra = np.flipud(aspectra)
        wavelengths = np.flipud(wavelengths)
        # Create some more data to be compatible with EXOhSPEC
        cspectra, sn_cont = normalise_continuum(spectra, wavelengths)
        save_multispec([wavelengths, spectra, spectra*0, spectra*0, spectra*0, cspectra, sn_cont, spectra*0, aspectra], params['path_extraction']+fname[1], im_head, bitpix=im_head['BITPIX'])
        add_text_to_file(params['path_extraction']+fname[1], 'plot_files.lst')
        logger('Info: File created: {0}'.format(params['path_extraction']+fname[1]))
        
