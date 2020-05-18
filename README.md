# HiFLEx
A highly flexible package to reduce echelle data (taken with a bifurcated or single fiber)

The Package is described in https://ui.adsabs.harvard.edu/abs/2020PASP..132f4504E/abstract [[pdf](https://ui.adsabs.harvard.edu/link_gateway/2020PASP..132f4504E/PUB_PDF)].

For information on how to use the software please check the Manual: [[HiFLEx_UserManual.pdf](https://github.com/ronnyerrmann/HiFLEx/blob/master/HiFLEx_UserManual.pdf)].

If you publish data using the barycentric correction, please cite (https://github.com/shbhuk/barycorrpy#citation).

If you publish radial velocities using [[TERRA](https://adsabs.harvard.edu/abs/2012ApJS..200...15A)], [[SERVAL](http://adsabs.harvard.edu/abs/2017A&A...609A..12Z)], or [[CERES](https://adsabs.harvard.edu/abs/2017PASP..129c4002B)], please cite the packages.


## Recent changes
* Bugfixing

## Install instruction

Requirements:
- python 2.7 + numpy, scipy, pyfits, astropy
- (only tested under linux)

Create a new [[Anaconda](https://www.anaconda.com/distribution/#linux)] environment:
```conda create --name hiflex python=2.7 numpy scipy matplotlib astropy pycurl ephem rpy2 tqdm psutil statsmodels 
conda activate hiflex
pip install gatspy barycorrpy==0.2.2.1 PyAstronomy multiprocessing
```

### TERRA (optional RV analysis)

### SERVAL (optional RV analysis)

### CERES (optional RV analysis)

## First steps
