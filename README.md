# HiFLEx
A highly flexible package to reduce echelle data (taken with a bifurcated or single fiber)

The Package is described in [Errmann et al.](https://ui.adsabs.harvard.edu/abs/2020PASP..132f4504E) (https://ui.adsabs.harvard.edu/abs/2020PASP..132f4504E).

For information on how to use the software please check the Manual: [HiFLEx_UserManual.pdf](https://github.com/ronnyerrmann/HiFLEx/blob/master/HiFLEx_UserManual.pdf).

If you publish data using the barycentric correction, please cite (https://github.com/shbhuk/barycorrpy#citation). If you publish radial velocities using [TERRA](https://ui.adsabs.harvard.edu/abs/2012ApJS..200...15A), [SERVAL](http://ui.adsabs.harvard.edu/abs/2017A&A...609A..12Z), or [CERES](https://ui.adsabs.harvard.edu/abs/2017PASP..129c4002B), please cite the packages (paper behind the links).


## Recent changes
* Python 3 ready (Pipeline will also keep working under Python 2)
* Blaze correction with a fitted function (polynomial of user defined order)
* Extra logging information for wavelength solution and drift of the wavelength solution
* Improved measurement of the offset between emission line spectra and the wavelength solution
* Bugfixing

## Install instruction

Requirements:
- python 3.8 (and 2.7) + numpy, scipy, pyfits, astropy
- (only tested under linux)

Create two new [Anaconda](https://www.anaconda.com/distribution/#linux) environments (python 3 required for cosmic ray removal, python 2.7 required for SERVAL and CERES):
```
conda create --name hiflex_p2 python=2.7 numpy scipy matplotlib astropy pycurl ephem rpy2 tqdm psutil statsmodels 
conda activate hiflex_p2
pip install gatspy barycorrpy==0.2.2.1 PyAstronomy multiprocessing
conda deactivate

conda create --name hiflex python=3 numpy scipy matplotlib astropy pycurl ephem rpy2 tqdm psutil statsmodels
conda activate hiflex
pip install gatspy barycorrpy deepCR PyAstronomy
```

Download and extract the [latest relase](https://github.com/ronnyerrmann/HiFLEx/releases) or clone the repository
git clone https://github.com/ronnyerrmann/hiflex.git

## First steps
Create a new folder and copy the `conf.txt` file from your HiFLEx installation path into this folder. Edit the following entries:
- Change parameter `raw_data_paths` to point to the folder or folders with your raw or pre-reduced fits-files. The data can be stored in sub-folders.
- Set the parameter `badpx_mask_filename` if a bad-pixel-mask exists. Otherwise leave empty or 'NA'.
- Adjust `rotate_frame` and `flip_frame` so that the orders are up to down (blue wavelength on top) and red orders are on the left. The pipeline will first rotate and then flip, the flip is swapping left and right.
- Adjust parameter `subframe` to a subframe, if only part of the CCD should be used. Otherwise leave empty or set to the full detector size.
- Set parameter `original_master_traces_filename = ` (empty), or to a non-existing file.
- To increase signal to noise and to speed up the search for the traces of the orders, binning as given in the parameters `bin_search_apertures` (rough estimate) and `bin_adjust_apertures` (fine tuning) can be adjusted. Please note that after binning, the orders should still be well distinguished from each other.
- The side of the calibration traces compared to the science traces is given in parameter `arcshift_side` (left or right). Set it to center for a single fiber spectrograph.
- Change the path to the catalogue of the reference lines (`reference_catalog`), if necessary. The file must contain one entry per line, each entry consists of tab-separated wavelength, line strength, and element. Line strength can be empty.
- Set parameter `original_master_wavelensolution_filename = master_wavelength_manual.fits`. To find a wavelength solution a file which provides the wavelength for some (extracted) pixel and wavelengths can be given optionally (called *pixel_to_wavelength.txt*).
- The number of degrees of freedom for the wavelength solution (2-dimensional polynomial fit) can be adjusted with `polynom_order_traces` (polynomial orders along the dispersion direction) and `polynom_order_intertraces` (polynomial orders between the traces).
- Adjust the parameters starting with `raw_data_*`. The example configuration file shows the settings for data taken with MaximDL and for data from the HARPS spectrograph.
- Adjust the parameters `*_calibs_create_g` to define the reduction steps which should be applied. Some of the `*_calibs_create_g` parameters might not be relevant, e.g. if no *dark* or *rflat* (real flatfield) corrections should be applied.
- The parameters `site`, `altitude`, `latitude` (negative for west), and `longitude` need to be set, if the information is not stored in the fits-header
- For optional RV analysis:  Set path `terra_jar_file` to the TERRA *PRV.jar* file. Set paths `path_serval` and `path_ceres` to the home folders of the SERVAL and CERES pipelines.
- Run the scrip `file_assignment.py`:
  * Define what files should be used for what calibration.
  * Define what files to extract (and in which RVs will be measured).
  * Please note that you can use already reduced images. In this case the reduction steps as defined in `conf.txt` should be empty for the file type (e.g. dark or real flats) to avoid a second application of the correction.
- Run the scrip `hiflex.py`
- Afterwards: check the output in the `logfile`, the images in the logging path, or in results in the extracted files.

More information can be found in the [manual](https://github.com/ronnyerrmann/HiFLEx/blob/master/HiFLEx_UserManual.pdf).

## Optional radial velocity packages
### TERRA (optional RV analysis)
Download the [TERRA achive](https://drive.google.com/file/d/1xK-lYghFwpwtdXG9b4IbryYRd102q7So/view) and extract. Set the variable **terra_jar_file** in `conf.txt` to the full path of the *PRV.jar* file.
To check that all dependencies are installed on the system one can run.
```
java -jar <full/path/to/terra>/terra/PRV.jar
```
This should produce the entry _*** TERRA v1.8 ***_ before failing. If not, please check that java is installed
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






