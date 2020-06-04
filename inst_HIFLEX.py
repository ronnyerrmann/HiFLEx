from read_spec import *

# Instrument parameters
name = inst = 'HIFLEX'
pat = '*.fits'              # Pattern of the file name
# obsname = (28.75728, -17.88508, 2382) from wiki
obsname = "hiflex" # for barycorrpy, not import if the BJD is in the header

hiflex_file = 'conf_hiflex_serval.txt'

params = dict()
if os.path.isfile(hiflex_file):
    try:
        keys, values = np.genfromtxt(hiflex_file, dtype=str, comments='#', delimiter='=', filling_values='', autostrip=True, unpack=True)
        params = dict(zip(keys, values))
    except ValueError as error:
        print('Problems when reading {0}'.format(hiflex_file))
else:
    print('Warn: No {0} found.'.format(hiflex_file))
for key in params.keys():
    params[key] = int(params[key])
iomax = params.get('orders', 1000)           # 1000: stupidly big number for orders. Real control is done with -oset
hdatf = params.get('data_dataset', 5)        # 1: higher scatter of resulting RVs; 5: no information about flux (continuum corrected)
hdate = params.get('error_dataset', 6)
hdatw = params.get('wave_dataset', 9)        # 9: Wavelength, drift correcte; 0: + barycentric correction
hdats = params.get('mask_dataset', 7)
print('Info: Using dataset {0} for the spectral data, {1} for the uncertainty, {2} for the wavelength, and {3} for the mastk'.format(hdatf, hdate, hdatw, hdats))

maskfile = 'telluric_mask_atlas.dat'

# Instrument read functions
def scan(self, s, pfits=True):
   """
   Returns
   -------
   namedtuple('spectrum', 'w f berv bjd blaze drift timeid sn55 ')
           w    - wavelength
           f    - flux
           berv - Barycentric Earth Radial Velocity
           bjd  - Barycentric Julian Day
           blaze - Blaze filename
           drift - Used RV Drift
           sn55  - S_N order center55
   Example
   -------

   """
   HIERARCH = 'HIERARCH '
   HIERHIFLEX = HIERARCH + 'HiFLEx '
   hdulist = self.hdulist = pyfits.open(s) # slow 30 ms
   if 1:
      self.header = self.hdr = hdr = hdulist[0].header
      iomax = hdr['NAXIS2']
      self.instname = 'HIFLEX'  # hdr['INSTRUME']
      self.drsberv = hdr.get(HIERHIFLEX + 'BCV', np.nan)
      self.drsbjd = hdr.get(HIERHIFLEX + 'BJDTDB', np.nan)
      self.dateobs = hdr[HIERHIFLEX + 'DATE-OBS']
      self.mjd = hdr.get(HIERHIFLEX + 'MJD')
      # for HiFLEx spectra the drift is already included in the wavelength solution
      #self.drift = hdr.get(HIERHIFLEX+'D_SHIFT_KMS', np.nan)    
      #self.e_drift = hdr.get(HIERHIFLEX+'CARACAL DRIFT FP E_RV', np.nan)                         # Adapt HiFLEx so this value could be given
      self.sn55 = hdr.get(HIERHIFLEX + 'SN_order%2.2i'%(int(iomax/2)), 50)                        # Modfied, flexible, maybe fix it to one order in HiFLEx
      self.fileid = self.timeid = hdr.get(HIERHIFLEX + 'DATE-OBS', 0)
      self.calmode = "%s,%s,%s" % (hdr.get('SCI-OBJ', ''), hdr.get('CAL-OBJ', ''), hdr.get('SKY-OBJ', ''))
   #calmodedict = {'objcal':'OBJ,CAL','objsky':'OBJ,SKY'}
   #if calmode in calmodedict: calmode = calmodedict[calmode]

      self.ccf.rvc = hdr.get(HIERHIFLEX+'RV_BARY___', np.nan)
      self.ccf.err_rvc = hdr.get(HIERHIFLEX+'RV_ERR____', np.nan)

      self.ra = hdr[HIERHIFLEX + 'RA']
      self.de = hdr[HIERHIFLEX + 'DEC']
      self.airmass = hdr.get(HIERHIFLEX + 'AIRMASS', np.nan)
      self.exptime = hdr[HIERHIFLEX+'EXPOSURE']
      self.tmmean = hdr.get(HIERHIFLEX+'EXP_FRAC', 0.5)
      if 'OBJECT' not in hdr.keys():
         hdr['OBJECT'] = hdr.get(HIERHIFLEX+'OBJNAME', 'Dummy')

def data(self, orders=None, pfits=True):
   if 1:  # read order data
      if hasattr(self, 'hdu'):   # read directly
         print 1111
         data = self.hdu.getdata()
      else:
         if not hasattr(self, 'hdulist'):
            scan(self, self.filename)
         data = self.hdulist[0].section[:]
      data = data.astype(dtype=np.float64)
      f = data[hdatf,orders,:]          # higher scatter of resulting RVs
      #f = data[5,orders,:]          # no information about flux (continuum corrected)
      e = data[hdate,orders,:]
      #w = data[0,orders,:]          # Wavelength, drift corrected, barycentric correction
      w = data[hdatw,orders,:]          # Wavelength, drift corrected
      s = data[hdats,orders,:]          # Mask

      bpmap = np.isnan(f).astype(int)   # flag 1 for nan
      e[np.isnan(e)] = 0
      e[e<0] = np.median(e[e>0],axis=None)

      with np.errstate(invalid='ignore'):
         bpmap[f <= 0.001] |= flag.neg      # flag 2 for zero and negative flux
         bpmap[f < -3*e] |= flag.neg      # flag 2 for zero and negative flux
         bpmap[s == 0.2] |= flag.neg      # bad-pixel
         bpmap[s == 0.1] |= flag.sat    # saturated
         bpmap[e==0] |= flag.nan

      w = airtovac(w)
      
      return w, f, e, bpmap
    



 
