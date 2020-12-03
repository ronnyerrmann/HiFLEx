from procedures import *

# =============================================================================
# Define variables
# =============================================================================
#global calimages
params = dict()     # default param dictionary
#calimages = dict()  # dictionary for all calibration images
# location of config file
CONFIGFILE = 'conf.txt'

# Start of code
# deal with arguments from a text file
params = textfileargs(params, CONFIGFILE)
global calimages    # dictionary for all calibration images

for calib in [ ['_cal','cal2','cal'], ['_sci','cal1', 'sci'] ]:
    if 'master_wavelensolution'+calib[0]+'_filename' not in params.keys():
        params['master_wavelensolution'+calib[0]+'_filename'] = params['master_wavelensolution_filename'].replace('.fit',calib[0]+'.fit')
    if os.path.isfile(params['master_wavelensolution'+calib[0]+'_filename']):
        wave_sol_dict = read_wavelength_solution_from_fits(params['master_wavelensolution'+calib[0]+'_filename'])
        calimages['wave_sol_dict_'+calib[2]] = copy.deepcopy( wave_sol_dict )
if 'wave_sol_dict_sci' not in calimages.keys():                             # Use the calibration fiber solution for the science solution
        calimages['wave_sol_dict_sci'] = calimages['wave_sol_dict_cal']

rv_results_to_hiflex(params)

header_keywords = []
if False:       # Change to True to use below list
        header_keywords.append(['HIERARCH HiFLEx OBJNAME',      'Object name',                  ''  ])
        header_keywords.append(['HIERARCH HiFLEx EXPOSURE',     'Exposure time',                '[s]'])
        header_keywords.append(['HIERARCH HiFLEx DATE-OBS',     'UTC, Begin of expose',         ''  ])
        header_keywords.append(['HIERARCH HiFLEx DATE-MID',     'Middle of exposure',           ''  ])
        header_keywords.append(['HIERARCH HiFLEx JD',           'JD at middle of exposure',     '[d]'])
        header_keywords.append(['HIERARCH HiFLEx fsum_all',     'Total extracted flux',         '[ADU]'])
        header_keywords.append(['HIERARCH HiFLEx BCKNOISE',     'Background noise',             '[ADU]'])
        header_keywords.append(['HIERARCH HiFLEx CD_SHIFT',     'Offset in Cross-Dispersion',   '[px]'])
        header_keywords.append(['HIERARCH HiFLEx CD_S_WDTH',    'Width of the offset in CD',    '[px]'])
        header_keywords.append(['HIERARCH HiFLEx CD_SHIFT_GAUSS',     'Offset in Cross-Dispersion (Gauss)',   'px'  ])
        header_keywords.append(['HIERARCH HiFLEx CD_SHIFT_POLY',      'Offset in Cross-Dispersion (Poly)',    'px'  ])
        header_keywords.append(['HIERARCH HiFLEx CD_SHIFT_MAXFLUX',   'Offset in Cross-Dispersion (Maxflux)', 'px'  ])
        header_keywords.append(['HIERARCH HiFLEx D_SHIFT',      'Offset in Dispersion',         '[px]'])
        header_keywords.append(['HIERARCH HiFLEx D_SHIFT_ERR',  'Uncertainty of the offset in D',               '[px]'])
        header_keywords.append(['HIERARCH HiFLEx D_SHIFT_NUMBER_LINES', 'Number of lines used to calculate the shift', ''])
        header_keywords.append(['HIERARCH HiFLEx D_SHIFT_GAUSS',      'Offset in Dispersion (Gauss)',         'px'  ])
        header_keywords.append(['HIERARCH HiFLEx D_WIDTH',      'Gaussian width of the calibration lines',      '[px]'])
        header_keywords.append(['HIERARCH HiFLEx D_SHIFT_KMS',  'Offset in Dispersion',                         '[km/s]'])
        header_keywords.append(['HIERARCH HiFLEx DT_SHIFT',     'Offset between wavelength solutions',          '[px]'])
        header_keywords.append(['HIERARCH HiFLEx DT_SHIFT_KMS', 'Offset between wavelength solutions',          '[km/s]'])
        header_keywords.append(['HIERARCH HiFLEx BJDTDB',       'Barycentric correct JD (incl. leap seconds)',  '[d]'])
        header_keywords.append(['HIERARCH HiFLEx CERES RV',     'CERES RV (not corrected for BVC)',             '[km/s]'])
        header_keywords.append(['HIERARCH HiFLEx CERES RV_BARY','CERES RV (corrected for BVC)', '[km/s]'])
        header_keywords.append(['HIERARCH HiFLEx CERES RV_ERR', 'CERES RV error',               '[km/s]'])
        header_keywords.append(['HIERARCH HiFLEx BCV',          'Barycentric Velocity',         '[km/s]'])
        header_keywords.append(['HIERARCH HiFLEx CERES BS',     'CERES Bisector',               ''  ])
        header_keywords.append(['HIERARCH HiFLEx CERES BS_ERR', 'CERES Bisector error',         ''  ])
        
header_results_to_texfile(params, header_keywords)


