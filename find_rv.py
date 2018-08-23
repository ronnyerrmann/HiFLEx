import numpy as np
import os
from procedures import *
import sys
base = '/home/ronny/software/ceres-master/'
sys.path.append(base+"utils/Correlation")
#sys.path.append(base+"utils/OptExtract")
sys.path.append(base+"utils/GLOBALutils")
sys.path.append(base+"utils/OptExtract")
sys.path.append(base+"utils/CCF")
#sys.path.append(base+'utils/SSEphem/')
baryc_dir= base+'utils/SSEphem/'
sys.path.append(baryc_dir)
ephemeris='DEc403'

import GLOBALutils
import correlation
import statsmodels.api as sm
lowess = sm.nonparametric.lowess
import pickle
import ephem
from math import radians as rad
import jplephem

# =============================================================================
# Define variables
# =============================================================================
params = dict()     # default param dictionary
# location of config file
#CONFIGFILE = 'conf.txt'

# rm analysis_files.lst ; ls extracted/SunArc*2s.fits > analysis_files.lst ; python ~/Scripts/exohspec/find_rv.py | tee RV_logfile ; grep 'RV = ' RV_logfile | cut -d' ' -f 3,5

def mjd_fromheader(head):       # from CERES, slightly modified
    """
    return modified Julian date from header
    """
    secinday = 24*3600.0

    datetu   = head['DATE-OBS'][:10] 
    ut       = head['DATE-OBS'][11:]

    mjd0,mjd,i = GLOBALutils.iau_cal2jd(int(datetu[0:4]),int(datetu[5:7]),int(datetu[8:10]))

    ut        = (float(ut[:2])*3600. + float(ut[3:5])*60. + float(ut[6:]))
    mjd_start = mjd + ut / secinday
    
    if 'HIERARCH ESO INS DET1 TMMEAN' in head:
        fraction = head['HIERARCH ESO INS DET1 TMMEAN']
    elif 'ESO INS DET1 TMMEAN' in head:
        fraction = head['ESO INS DET1 TMMEAN']
    else:
        fraction = 0.5
    texp     = head['EXPTIME'] #sec

    mjd = mjd_start + (fraction * texp) / secinday

    return mjd, mjd0

def logger(message, show=True, printarrayformat=[], printarray=[], logfile='logfile'):
    """
    Saves the status information to a logfile
    :param message: Text to log
    :param show: if False than it will be logged only to the logfile but not printed on the screen
    :param printarrayformat: Format of the columns in printarray, e.g. '%3.1f'
    :param printarray: Array with values to log
    :param logfile: filename in which the information will be written
    """
    if show:
        print message
    file = open(logfile, 'a')
    file.write(time.strftime("%Y%m%d%H%M%S - ", time.localtime()) + message + '\n')
    if printarrayformat <> [] and printarray <> []:
        for line in printarray:
            text = '\t'
            for i,printformat in enumerate(printarrayformat):
                #print printformat,line[i]
                text += '\t'+printformat%line[i]
            file.write(text+'\n')
            if show:
                print text
    file.close()
    if message.find('Error') == 0:
        print '\t-> exiting'
        exit(1)
  
def get_spec(spectra, waves, cen_wave, range_px):
    """
    :param spectra: 1d or 2d array of floats with the spectral data
    :param waves: same format as spectra, corresponding wavelengths
    :param cen_wave: wave_length to get the spectra around
    :param range_px: how many pixel around that wavelength
    :return spec, waves: returns the spectra and the wavelengths of the area specified
    """
    if len(spectra.shape) < 2:                      # make into a 2d spectrum, if necessary
        spectra = spectra.reshape(1, spectra.shape[0])
        waves = waves.reshape(1, waves.shape[0])
    wavediff = np.nanmax( np.abs( waves[:,1:] - waves[:,:-1] ))     # maximum wavelengths difference between consecutive pixel
    diff = waves - cen_wave
    pos = np.where( np.abs(diff) <= wavediff )      # gives a 2d array for the positions (order, px)
    diffp = np.abs(pos[1] - spectra.shape[1]/2)     # what is the most central px
    for dummy in np.unique(pos[0]):                 # unique ignores if wavelength is not covered in the spectrum
        pos2 = np.argmin(diffp)                     # what is the most central pixel
        order, px = pos[0][pos2], pos[1][pos2]
        spec = spectra[order, max(0,px-range_px) : min(spectra.shape[1], px+range_px+1) ]
        if len( spec[~np.isnan(spec)] ) == 0:           # no spectrum for this data
            pos[1][pos[0] == order] = 1E9               # -> check the others orders
            diffp = np.abs(pos[1] - spectra.shape[1]/2)
        else:
            #print spec, waves[order, max(0,px-range_px) : min(spectra.shape[1], px+range_px+1) ]
            return spec, waves[order, max(0,px-range_px) : min(spectra.shape[1], px+range_px+1) ]

    return [], []                                      # everything went wrong

def rv_analysis(params, spec, im_head, fitsfile):
    specs = spec.shape
    spec = np.vstack(( spec, np.zeros([11-specs[0], specs[1], specs[2]]) ))
    spec[7:11,:,:] = np.nan
    reffile = 'reffile.txt'
    obname = fitsfile.split('/')
    obname = obname[-1].split('-')
    obname = obname[0]
    lbary_ltopo = 1
    npools = 1
    refvel = 0
    know_moon = False
    here_moon = False
    models_path = params['rv_models_ceres'] # "/home/ronny/software/ceres-master/data/COELHO_MODELS/R_40000b/"
    #dirout = './'          # replaced by params['path_rv']
    fsim = fitsfile
    avoid_plot = False
    RESI = 120000.
    force_stellar_pars = False
    
    
    #ron1,gain1 = h[1].header['HIERARCH ESO DET OUT1 RON'],h[1].header['HIERARCH ESO DET OUT1 GAIN']
    #ron2,gain2 = h[2].header['HIERARCH ESO DET OUT1 RON'],h[2].header['HIERARCH ESO DET OUT1 GAIN']
    #halfcounts = h[0].header['HIERARCH ESO INS DET1 TMMEAN']

    # Find lambda_bary/lambda_topo using baryc
    site_keys       = ['TELESCOP']
    altitude_keys   = ['HIERARCH ESO TEL GEOELEV', 'ESO TEL GEOELEV']
    latitude_keys   = ['HIERARCH ESO TEL GEOLAT', 'ESO TEL GEOLAT']
    longitude_keys  = ['HIERARCH ESO TEL GEOLON', 'ESO TEL GEOLON']
    ra_keys         = ['RA']
    dec_keys        = ['DEC']
    epoch_keys      = ['HIERARCH ESO TEL TARG EQUINOX', 'ESO TEL TARG EQUINOX']
    params['site']        = 'UH'
    params['altitude']    = 30
    params['latitude']    = 51.7534
    params['longitude']   = -0.2401
    params['epoch']       = 2000.0
    settings = []
    settings.append( [0, site_keys, 'site'] )
    settings.append( [0, altitude_keys, 'altitude'] )
    settings.append( [0, latitude_keys, 'latitude'] )
    settings.append( [0, longitude_keys, 'longitude'] )
    settings.append( [0, ra_keys, 'ra'] )
    settings.append( [0, dec_keys, 'dec'] )
    settings.append( [0, epoch_keys, 'epoch'] )
    for [i, keys, parentr] in settings:
        for entry in keys:
            if entry in im_head.keys():
                params[parentr] = im_head[entry]    # Get the information from the header
                if entry == 'ra':                   # Assume that dec is coming from the same source
                    source_radec = 'The object coordinates are derived from the image header'
        if parentr == 'longitude':                  # Enough information to calculate the ephemerides of sun and moon
            gobs = ephem.Observer()  
            gobs.name = params['site']  
            gobs.lat  = rad(params['latitude'])     # lat/long in decimal degrees  
            gobs.long = rad(params['longitude'])
            DDATE     = im_head['DATE-OBS'][:10]
            HHOUR     = im_head['DATE-OBS'][11:]
            gobs.date = str(DDATE[:4]) + '-' +  str(DDATE[5:7]) + '-' + str(DDATE[8:]) + ' ' +  HHOUR[:2] + ':' + HHOUR[3:5] +':' + str(float(HHOUR[6:]) + 0.5*im_head[params['raw_data_exptim_keyword']] )
            mephem    = ephem.Moon()
            mephem.compute(gobs)
            sephem    = ephem.Sun()
            sephem.compute(gobs)
            if obname.find('Sun') == 0:
                params['ra']          = sephem.ra
                params['dec']         = sephem.dec
                source_radec = 'The object coordinates are derived from the calculated solar ephermeris'
            elif obname.find('Moon') == 0:
                params['ra']          = mephem.ra
                params['dec']         = mephem.dec
                source_radec = 'The object coordinates are derived from the calculated lunar ephermeris'
            else:
                params['ra']          = mephem.ra       # To fill in the parameter
                params['dec']         = mephem.dec      # To fill in the parameter
                source_radec = 'The object coordinates were made up'
    
    mjd,mjd0 = mjd_fromheader(im_head)
    ra2,dec2 = GLOBALutils.getcoords(obname,mjd,filen=reffile)
    if ra2 !=0 and dec2 != 0:
        params['ra']  = ra2
        params['dec'] = dec2
        source_radec = 'The object coordinates are derived from the reference file {0}'.format(reffile)
    
    logger(('Info: Using the following data: altitude = {0}, latitude = {1}, longitude = {2}, ra = {3}, dec = {4}, epoch = {5}. '+\
           '{6}').format(params['altitude'], params['latitude'], params['longitude'], 
                                                                params['ra'], params['dec'], params['epoch'], source_radec ))
    ra, dec = params['ra'], params['dec']

    iers          = GLOBALutils.JPLiers( baryc_dir, mjd-999.0, mjd+999.0 )
    obsradius, R0 = GLOBALutils.JPLR0( params['latitude'], params['altitude'])
    obpos         = GLOBALutils.obspos( params['longitude'], obsradius, R0 )
    jplephem.set_ephemeris_dir( baryc_dir , ephemeris )
    jplephem.set_observer_coordinates( obpos[0], obpos[1], obpos[2] )

    res         = jplephem.doppler_fraction(ra/15.0, dec, int(mjd), mjd%1, 1, 0.0)
    lbary_ltopo = 1.0 + res['frac'][0]
    bcvel_baryc = ( lbary_ltopo - 1.0 ) * 2.99792458E5

    #logger("Info: Barycentric velocity: {0}, MJD: {1}".format( bcvel_baryc, mjd))
    """ Data description of the file
    0: wavelength for each order and pixel'
    1: extracted spectrum'
    2: measure of error (photon noise, read noise)'
    3: flat corrected spectrum'
    4: error of the flat corrected spectrum (residuals to a {0} order polynomial)'.format(measure_noise_orders)
    5: continuum normalised spectrum'
    6: error in continuum (fit to residuals of {0} order polynomial)'.format(measure_noise_orders)
    7: Mask with good areas of the spectrum: 0.1=saturated_px, 0.2=badpx'
    8: spectrum of the emission line lamp'
    """
    for order in range(spec.shape[1]):
        L  = np.where( (spec[1,order,:] != 0) & (~np.isnan(spec[1,order,:])) )              # good values
        #ratio              = np.polyval(ccoefs[order],spec[0,order,:][L])*Rnorms[order]
        ratio = spec[1,order,:][L]/  spec[5,order,:][L]                                     # ratio between extracted spectrum and continuum normalised spectrum -> blaze function, cancels absorption lines
        spec[7,order,:][L] = ratio
        #spec[8,order,:][L] = spec[6,order,:][L]                                             # error continuum (first guess), but not good. #sn_order=8
        spec[8,order,:][L] = spec[2,order,:][L]                                             # error of the extracted data, depending on what is used, the RV changes by few 100 km/s -> > several \AA
        #spec[8,order,:][L] = ratio * R_flat_ob_n[order,1,:][L] / np.sqrt( ratio * R_flat_ob_n[order,1,:][L] / gain2 + (ron2/gain2)**2 )        # something like S/N -> used as this by XCor
        #spl           = scipy.interpolate.splrep(np.arange(WavSol.shape[0]), WavSol,k=3)
        #dlambda_dx    = scipy.interpolate.splev(np.arange(WavSol.shape[0]), spl, der=1)
        #NN            = np.average(dlambda_dx)
        #dlambda_dx    /= NN
        LL = np.where(spec[5,order,:] > 1 + 10. / scipy.signal.medfilt(spec[8,order,:],21))[0]          # remove emission lines and cosmics
        spec[5,order,LL] = 1.
        spec[9,order,:][L] = spec[5,order,:][L]# * (dlambda_dx[L] ** 1)         # used for the analysis in XCor (spec_order=9, iv_order=10)
        spec[10,order,:][L] = spec[2,order,:][L]# / (dlambda_dx[L] ** 2)        # used for the analysis in XCor (spec_order=9, iv_order=10)
    #plot_img_spec.plot_spectra_UI(np.array([spec]))
    T_eff, logg, Z, vsini, vel0 = 5777, 4.4374, 0.0134, 2, 0
    if True:  
        if False:       # is not run
            pars_file = params['path_rv'] + fsim.split('/')[-1][:-4]+'_stellar_pars.txt'

            if os.access(pars_file,os.F_OK) == False or force_stellar_pars:
                print "\t\t\tEstimating atmospheric parameters:"
                Rx = np.around(1./np.sqrt(1./40000.**2 - 1./RESI**2))
                spec2 = spec.copy()
                for i in range(spec.shape[1]):
                    IJ = np.where(spec[5,i]!=0.)[0]
                    spec2[5,i,IJ] = GLOBALutils.convolve(spec[0,i,IJ],spec[5,i,IJ],Rx)
                T_eff, logg, Z, vsini, vel0, ccf = correlation.CCF(spec2,model_path=models_path,npools=npools, base='/home/ronny/software/ceres-master/utils/Correlation/')     # Fails, because our spectrum doesn't cover the hard coded wavelength
                line = "%6d %4.1f %4.1f %8.1f %8.1f\n" % (T_eff,logg, Z, vsini, vel0)
                f = open(pars_file,'w')
                f.write(line)
                f.close()
               
            else:
                print "\t\t\tAtmospheric parameters loaded from file:"
                T_eff, logg, Z, vsini, vel0 = np.loadtxt(pars_file,unpack=True)
        
        print "\t\tRadial Velocity analysis:"
        # assign mask
        #   obname is the name of the object
        #   reffile is a reference file: reffile = dirin+'reffile.txt' -> this file doesn't exist
        sp_type, mask = GLOBALutils.get_mask_reffile(obname,reffile=reffile,base='/home/ronny/software/ceres-master/data/xc_masks/')
        print "\t\t\tWill use",sp_type,"mask for CCF."

        # Read in mask
        ml, mh, weight = np.loadtxt(mask,unpack=True)
        ml_v = GLOBALutils.ToVacuum( ml )
        mh_v = GLOBALutils.ToVacuum( mh )
        av_m = 0.5*( ml_v + mh_v )
        mask_hw_kms = (GLOBALutils.Constants.c/1e3) * 0.5*(mh_v - ml_v) / av_m

        disp = GLOBALutils.get_disp(obname, reffile=reffile)
        if disp == 0:
            known_sigma = False
            if vsini != -999 and vsini != 0.:
                disp = vsini
            else:
                disp = 3.
        else:
            known_sigma = True

        mask_hw_wide = av_m * disp / (GLOBALutils.Constants.c/1.0e3)
        ml_v = av_m - mask_hw_wide
        mh_v = av_m + mask_hw_wide 

        print '\t\t\tComputing the CCF...'
        cond = True

        if sp_type == 'M5':
            moon_sig = 4.5
        elif sp_type == 'K5':
            moon_sig = 4.2
        else:
            moon_sig = 4.0

    
        while (cond):
            # first rough correlation to find the minimum
            #   spec: spectrum in the form [data type, orders, pixel]
            #   ml_v, mh_v: masks 
            #   lbary_ltopo = 1.0 + res['frac'][0]
            vels, xc_full, sn, nlines_ccf, W_ccf = \
                    GLOBALutils.XCor(spec, ml_v, mh_v, weight,\
                    0, lbary_ltopo, vel_width=300, vel_step=3,\
                    spec_order=9, iv_order=10, sn_order=8, max_vel_rough=300)

            xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=3.0, Simple=True, W=W_ccf)

            # Normalize the continuum of the CCF robustly with lowess     
            yy = scipy.signal.medfilt(xc_av,11)
            pred = lowess(yy, vels,frac=0.4,it=10,return_sorted=False)
            tck1 = scipy.interpolate.splrep(vels,pred,k=1)
            xc_av_orig = xc_av.copy()
            xc_av /= pred
            vel0_xc = vels[ np.argmin( xc_av ) ] 
                
            rvels, rxc_av, rpred, rxc_av_orig, rvel0_xc = \
                    vels.copy(), xc_av.copy(), pred.copy(),\
                    xc_av_orig.copy(), vel0_xc

            xc_av_rough = xc_av
            vels_rough  = vels
                
            vel_width = np.maximum( 20.0, 6*disp )                      # Adjusted in order to avoid crashes because of unphysical disp
            #print vel_width, disp, vsini       # problem with vel_width, due to disp, due to p1gau below
            vels, xc_full, sn, nlines_ccf, W_ccf =\
                    GLOBALutils.XCor(spec, ml_v, mh_v, weight,\
                    vel0_xc, lbary_ltopo, vel_width=vel_width,\
                    vel_step=0.1, spec_order=9, iv_order=10, sn_order=8,max_vel_rough=300)

            xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=3.0, Simple=True, W=W_ccf)
            pred = scipy.interpolate.splev(vels,tck1)
            xc_av /= pred

            p1,XCmodel,p1gau,XCmodelgau,Ls2 = \
                    GLOBALutils.XC_Final_Fit( vels, xc_av, sigma_res = 4,\
                     horder=8, moonv=refvel, moons=moon_sig, moon=False)
            print 'plgau 345', p1gau

            moonmatters = False
            if (know_moon and here_moon):
                moonmatters = True
                ismoon = True
                confused = False
                p1_m,XCmodel_m,p1gau_m,XCmodelgau_m,Ls2_m = GLOBALutils.XC_Final_Fit( vels, xc_av, \
                sigma_res = 4, horder=8, moonv = refvel, moons = moon_sig, moon = True)
                moon_flag = 1
            else:
                confused = False
                ismoon = False
                p1_m,XCmodel_m,p1gau_m,XCmodelgau_m,Ls2_m = p1,XCmodel,p1gau,XCmodelgau,Ls2
                moon_flag = 0

            bspan = GLOBALutils.calc_bss(vels,xc_av)
            SP = bspan[0]
            
            if (not known_sigma):
                disp = np.floor(p1gau[2])
                if (disp < 3.0): 
                    disp = 3.0
                mask_hw_wide = av_m * disp / (GLOBALutils.Constants.c/1.0e3)
                ml_v = av_m - mask_hw_wide
                mh_v = av_m + mask_hw_wide            
                known_sigma = True
            else:
                cond = False
                
            if p1gau[2] > 1E3:
                cond = False
                
        xc_dict = {'vels':vels,'xc_av':xc_av,'XCmodelgau':XCmodelgau,'Ls2':Ls2,'refvel':refvel,\
               'rvels':rvels,'rxc_av':rxc_av,'rpred':rpred,'rxc_av_orig':rxc_av_orig,\
               'rvel0_xc':rvel0_xc,'xc_full':xc_full, 'p1':p1, 'sn':sn, 'p1gau':p1gau,\
               'p1_m':p1_m,'XCmodel_m':XCmodel_m,'p1gau_m':p1gau_m,'Ls2_m':Ls2_m,\
               'XCmodelgau_m':XCmodelgau_m}

        #moon_dict = {'moonmatters':moonmatters,'moon_state':moon_state,'moonsep':moonsep,\
        #     'lunation':lunation,'mephem':mephem,'texp':im_head['EXPTIME']}
        moon_dict = {'moonmatters':moonmatters,'moon_state':'dummy','moonsep':0,\
             'lunation':0,'mephem':mephem,'texp':0}

        pkl_xc = params['path_rv'] + fsim.split('/')[-1][:-4]+obname+'_XC_'+sp_type+'.pkl'
        pickle.dump( xc_dict, open( pkl_xc, 'w' ) )

        ccf_pdf = params['logging_path'] + fsim.split('/')[-1][:-4] + obname + '_XCs_' + sp_type + '.pdf'       # dirout + 'logging/'

        if not avoid_plot:
            GLOBALutils.plot_CCF(xc_dict,moon_dict,path=ccf_pdf)

        #SNR_5130 = np.median(spec[8,28,1900:2101] ) This wavelength is not covered in exohspec
        #SNR_5130 = np.nanmedian(spec[8,0,1900:2101] ) # Set to order 20 artificially
        SNR_5130 = np.nan
        for cen_wave in [5130,6200,4000,7300,3000,8400,9500]:       # try different wavelengths in case one is not covered
            subspec, wave = get_spec(spec[8,:,:], spec[0,:,:], cen_wave, 100)
            if len(subspec) > 0:
                SNR_5130 = np.nanmedian(subspec)       # 5130 shows significant absorption lines
                break
            logger('Warn: The spectra around the wavelength of {0} Angstrom is not covered'.format(cen_wave))
        
        B,A = -0.00257864,0.07765779            # from Monte Carlo Simulation, different for each instrument
        B,A = 0.010, 5.0
        #print SNR_5130
        RVerr  =  B + ( 1.6 + 0.2 * p1gau[2] ) * A / np.round(SNR_5130)
        depth_fact = 1. + p1gau[0]/(p1gau[2]*np.sqrt(2*np.pi))
        if depth_fact < 0.6:
            depth_fact = 0.6
        depth_fact = (1 - 0.6) / (1 - depth_fact)
        RVerr *= depth_fact
        #print RVerr, depth_fact, p1gau, SNR_5130
        if RVerr < 0.002:
            RVerr = .002

        B,A = -0.00348879, 0.10220848
        BSerr = B + ( 1.6 + 0.2 * p1gau[2] ) * A / np.round(SNR_5130)
        if BSerr<0.002:
            BSerr = .002

        RV     = np.around(p1gau_m[1],4)  
        BS     = np.around(SP,4)   
        RVerr2 = np.around(RVerr,4)
        BSerr  = np.around(BSerr,4)
        bcvel_baryc = np.around(bcvel_baryc,4)

        return RV+bcvel_baryc, RVerr2, BS, BSerr, bcvel_baryc
        
        
# Start of code
if __name__ == "__main__":
    #if not os.path.exists('proc'):
    #    os.makedirs('proc')
    params['raw_data_exptim_keyword'] = 'EXPTIME'
    fitsfiles = []
    files = read_text_file('analysis_files.lst')
    for line in files:
        if line[0] <> '#':
            spec = np.array(fits.getdata(line))
            #spec = np.delete(spec, 0, axis=1)
            #spec = np.delete(spec, -1, axis=1)
            #print('!!! first and last order are deleted')
            #spec[0] = spec[0,:,:]+3
            im_head = fits.getheader(line)
            RV, RVerr2, BS, BSerr, bcvel_baryc = rv_analysis(params, spec, im_head, line)
            print '\t\t\tRV = '+str(RV)+' +- '+str(RVerr2)+' ,\tBS = '+str(BS)+' +- '+str(BSerr)+' ,\t'+str(bcvel_baryc)+' ,RV contains barycentric velocity correction'
        

