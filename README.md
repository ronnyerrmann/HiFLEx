# HiFLEx
A highly flexible package to reduce echelle data (taken with a bifurcated or single fiber)

The Package is described in https://ui.adsabs.harvard.edu/abs/2020PASP..132f4504E/abstract [pdf](https://ui.adsabs.harvard.edu/link_gateway/2020PASP..132f4504E/PUB_PDF).

For information on how to use the software please check the Manual: [HiFLEx_UserManual.pdf](https://github.com/ronnyerrmann/HiFLEx/blob/master/HiFLEx_UserManual.pdf).

If you publish data using the barycentric correction, please cite (https://github.com/shbhuk/barycorrpy#citation).

If you publish radial velocities using [TERRA](https://adsabs.harvard.edu/abs/2012ApJS..200...15A), [SERVAL](http://adsabs.harvard.edu/abs/2017A&A...609A..12Z), or [CERES](https://adsabs.harvard.edu/abs/2017PASP..129c4002B), please cite the packages.


## Recent changes
* Bugfixing

## Install instruction

Requirements:
- python 2.7 + numpy, scipy, pyfits, astropy
- (only tested under linux)

Create a new [Anaconda](https://www.anaconda.com/distribution/#linux) environment:
```
conda create --name hiflex python=2.7 numpy scipy matplotlib astropy pycurl ephem rpy2 tqdm psutil statsmodels 
conda activate hiflex
pip install gatspy barycorrpy==0.2.2.1 PyAstronomy multiprocessing
```

Download and extract the [latest relase](https://github.com/ronnyerrmann/HiFLEx/releases) or clone the repository
git clone https://github.com/ronnyerrmann/hiflex.git

## First steps
Create a new folder and copy the `conf.txt` file from your HiFLEx installation path into this folder. Edit the following entries:


## Optional radial velocity packages
### TERRA (optional RV analysis)
Download the [TERRA achive](https://drive.google.com/file/d/1xK-lYghFwpwtdXG9b4IbryYRd102q7So/view) and extract. Set the variable **terra_jar_file** in `conf.txt` to the full path of the *PRV.jar* file.
To check that all dependencies are installed on the system one can run.
```
java -jar <full/path/to/terra>/terra/PRV.jar
```
This should produce the entry \textit{*** TERRA v1.8 ***} before failing. If not, please check that java is installed
```
java --version
```

### SERVAL (optional RV analysis)
To install Serval please follow the instructions as given on (https://github.com/mzechmeister/serval). It is worth to check that the software is working by running the package on the provided test data. Set the variable **path_serval** in `conf.txt` to the full path of the *mzechmeister* folder, e.g. the path given in *$SERVALHOME*.

Serval can be controlled with an instrument file. This has to be provided in the serval **src** folder. Below is a way to create a soft link it when following the installing instructions.
```
ln -s <full/path/to/hiflex>/inst_HIFLEX*.py $SERVALHOME/serval/src/
```
Alternatively, the following line can be used:
```
ln -s <full/path/to/hiflex>/inst_HIFLEX*.py <path/to/mzechmeister>/serval/src/
```

### CERES (optional RV analysis)
To install Serval please follow the instructions as given on (https://github.com/rabrahm/ceres). Please note that the compilers in the CERES install process ignore any anconda settings and go straight to the system packages (e.g. for SSEPhem), hence you might want to install it when anaconda is not loaded.
Set the variable **path_ceres** in `conf.txt` to the full path of the *ceres* folder.






