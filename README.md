# likelihood_ratio
A python code for likelihood ratio matching between radio and multi-wavelength surveys, based on the methods described in, e.g. McAlpine et al. (2012).

**Software requirements**

The software in this repository also makes use of the Starlink Tables Infrastructure Library Tool Set (STILTS), which is a command line scriptable version of topcat.  Please see http://www.star.bris.ac.uk/~mbt/stilts/sun256/sun256.html for installation instructions. 

**python libraries**

This software uses the following python libraries:

argparse, numpy, pandas, os, astropy, scipy

**Preparing to run**

1. You will need the following, in binary fits table foramt:  

- Multi-wavelength catalogue with magnitudes for the bands you wish to match  
- Radio catalogue  
- Mask image for the multi-wavelength data

2. Set up the configuration file. Please look at the example configuration file and update the parameters as appropriate. This is meant for you to be able to use the software flexibly for catalogues where the column names might be slightly different (e.g., RA vs. ALPHA_J2000, etc.). Please surround strings with quotation marks. The following parameters are set in the configuration file:

**outdir**          *the out directory where data will be written*  
**bands**             *a comma-separated list of bands for which you want to run the likelihood ratio (case sensitive)*  
**id_col**            *the column name for the ID column*  
**ra_col**            *the column name for right ascension*  
**dec_col**           *the column name for declination*  
**halo_col**          *the column name for halo masks (for stars -- assumed 1 = flag, 0 = no flag)*  
**sg_col**             *the column name for star / galaxy separation (assumed 1 = galaxy, 0 = star)*  
**mag_col**           *an example column name for the magnitude column, with 'X' instead of the band -- e.g.: mag_X*  
**mag_err_col**       *same as mag_col but for the errors on the magnitude column -- e.g.: mag_err_X*  
**flux_col**          *Total flux column name from radio catalogue*  
**flux_err_col**      *Column name for errors on total flux*  
**beam_size**         *the beam size in arcsec*  
**cal_errors**        *calibration errors (e.g., 0.1 is 10 percent)*


**Running the software**

run  
> likelihood_ratio_matching.py -h  

to see the optional and required inputs.
