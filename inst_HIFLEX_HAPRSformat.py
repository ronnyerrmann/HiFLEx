from read_spec import *
from read_spec import Inst

# Instrument parameters
name = inst = 'HIFLEX'
pat = '*.fits'              # Pattern of the file name
# obsname = (28.75728, -17.88508, 2382) from wiki
obsname = "testobsname" # for barycorrpy, not import if the BJD is in the header

iomax = 1000          # stupidly big number for orders. Real control is done with -oset

maskfile = 'telluric_mask_atlas.dat'

def scan(self, s, pfits=True, verb=False):
   """
   SYNTAX: read_harps(filename)
   OUTPUT: namedtuple('spectrum', 'w f berv bjd blaze drift timeid sn55 ')
           w    - wavelength
           f    - flux
           berv - Barycentric Earth Radial Velocity
           bjd  - Barycentric Julian Day
           blaze - Blaze filename
           drift - Used RV Drift
           sn55  - S_N order center55

   """
   drs = self.drs
   if 1:
      HIERARCH = 'HIERARCH '
      HIERINST = HIERARCH + 'ESO '
      k_tmmean = HIERINST + 'INS DET1 TMMEAN'
      # In old HARPN the keyword is different and the value absolute
      #k_tmmean = {'HARPS': HIERINST + 'INS DET1 TMMEAN', 'HARPN': HIERINST + 'EXP1 TMMEAN'}[inst]
      if 1:
         self.HIERDRS = HIERDRS = HIERINST + 'DRS '
         iomax = hdr['NAXIS2']
         k_sn55 = HIERDRS + 'SPE EXT SN{0}'.format(int(iomax/2))            # Was 30 and changes the results
         k_berv = HIERDRS + 'BERV'
         k_bjd = HIERDRS + 'BJD'
      
      if pfits is True:  # header with pyfits
         self.hdulist = hdulist = pyfits.open(s) # slow 30 ms
         hdr = self.hdulist[0].header # pyfits.getheader(s)
      elif pfits==2:     # a faster version
         args = ('INSTRUME', 'OBJECT', 'MJD-OBS', 'DATE-OBS', 'OBS TARG NAME', 'EXPTIME',
                 'MJD-OBS', 'FILENAME', 'RA', 'DEC', k_tmmean, HIERINST+'DPR TYPE',
                 HIERINST+'DPR TECH', HIERINST+'INS MODE', HIERINST+'OBS TARG NAME')
         args += (k_bjd, k_berv, k_sn55)
         # args += (HIERINST+'OBS PI-COI NAME', HIERINST+'OBS PROG ID')
         if 1:
            args += (HIERDRS+'BLAZE FILE', HIERDRS+'DRIFT RV USED',
                     HIERDRS+'CAL TH DEG LL', HIERDRS+'CAL LOC NBO',
                     HIERDRS+'CAL TH COEFF LL')

         hdr = imhead(s, *args)
         self.hdu = getext(s)
      
      #self.drs = 'DRS CAL LOC NBO' in "".join(hdr.keys())  # check DRS or FOX
      self.instname = 'HIFLEX' #hdr['INSTRUME'] #if self.drs else 'HARPS'
      self.HIERARCH = HIERARCH
      
      self.airmass = hdr.get('AIRMASS', np.nan)
      self.exptime = hdr['EXPTIME']
      self.mjd = hdr.get('MJD-OBS')
      self.dateobs = self.fileid = self.timeid = fileid = hdr['DATE-OBS']
      self.ra = hdr['RA']
      self.de = hdr['DEC']
      self.utc = datetime.datetime.strptime(self.dateobs, '%Y-%m-%dT%H:%M:%S.%f')

      self.obs.lon = -70.7345               # Maybe adapt HiFLEx to get from header, shouldn't be used
      self.obs.lat = -29.2584               # Maybe adapt HiFLEx to get from header, shouldn't be used

      self.tmmean = hdr.get(k_tmmean, 0.5)

      self.drsbjd = hdr.get(k_bjd)
      if self.drsbjd is None:
         self.drsbjd = hdr.get('MJD-OBS')
      #if self.drsbjd: # comment out because the sa cannot be calculated with str
         #self.drsbjd = repr(self.drsbjd)
      self.drsberv = hdr.get(k_berv, np.nan)
      self.sn55 = hdr.get(k_sn55, np.nan)
      self.blaze = hdr.get(HIERDRS+'BLAZE FILE', 0)
      self.drift = hdr.get(HIERDRS+'DRIFT RV USED', np.nan)
      if abs(self.drift) > 1000:
         # sometimes there are crazy drift values ~2147491.59911, e.g. 2011-06-15T08:11:13.465
         self.drift = np.nan

      
      #calmode = hdr.get('IMAGETYP',0).split(",")[:2]
      calmode = hdr.get(HIERINST+'DPR TYPE','NOTFOUND').split(',')[:2]
      self.calmode = ','.join(calmode)
      calmodedict = {'STAR,WAVE': 'OBJ,CAL', 'STAR,DARK': 'OBJ,SKY'}
      if self.calmode in calmodedict: self.calmode = calmodedict[self.calmode]

      if hdr[HIERINST+'DPR TECH'] == 'ECHELLE,ABSORPTION-CELL':
         self.flag |= sflag.iod
      if hdr[HIERINST+'INS MODE'] == 'EGGS':
         self.flag |= sflag.eggs

      hdr['OBJECT'] = hdr.get('OBJECT', 'FOX')
      self.header = self.hdr = hdr # self.header will be set to None


def data(self, orders, pfits=True):
   hdr = self.hdr
   drs = self.drs
   if 1:  # read order data
      if hasattr(self, 'hdu'):   # read directly
         f = self.hdu.getdata(o=orders)
      else:
         if not hasattr(self, 'hdulist'):
            scan(self, self.filename)
         f = self.hdulist[0].section[orders]
         
      bpmap = np.isnan(f).astype(int)   # flag 1 for nan
      if drs:
         # print " applying wavelength solution ", file
         # omax = self.hdu['SPEC'].NAXIS1
         omax = hdr[self.HIERDRS+'CAL LOC NBO'] # 72 for A and 71 for B
         d = hdr[self.HIERDRS+'CAL TH DEG LL']
         xmax = hdr['NAXIS1']
         x = np.empty((d+1, xmax), 'int64')
         x[0].fill(1)                               # x[0,*] = x^0 = 1,1,1,1,1,...
         x[1] = np.arange(xmax)                     #        = x^1 = 0,1,2,3,4,...
         for i in range(1,d): x[i+1] = x[i] * x[1]  #        = x^i
         if not hasattr(self, 'A'):
            #A = np.array([hdr[self.HIERDRS+'CAL TH COEFF LL'+str(i)] for i in range(omax*(d+1))],dtype='float64').reshape(omax,d+1) #slow 30 ms
            self.A = np.reshape([hdr[self.HIERDRS+'CAL TH COEFF LL'+str(i)] for i in range(omax*(d+1))], (omax,d+1)) #slow 30 ms
         w = np.dot(self.A[orders], x)  # wavelength lambda
         e = np.sqrt(np.where(bpmap, 0., 5**2 * 6 + np.abs(f, dtype=float)))


      """with np.errstate(invalid='ignore'):
         bpmap[f < -3*e] |= flag.neg      # flag 2 for zero and negative flux
         bpmap[f > 300000] |= flag.sat    # estimate for saturation level:
                                       # HARPS.2004-10-03T01:30:44.506.fits:
                                       # last order: e2ds_B: 346930 (x=2158) raw: 62263 (y=1939)"""

      w = airtovac(w)

      return w, f, e, bpmap


