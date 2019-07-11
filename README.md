# R_likelihood_ratio
An R code for likelihood ratio matching between radio and multi-wavelength surveys, based on the methods described in, e.g. McAlpine et al. (2012).

**Software requirements**

R is a free, widely distributed software that is usually available on UNIX systems.  Please see https://www.r-project.org/ for more details. 

The software in this repository also makes use of the Starlink Tables Infrastructure Library Tool Set (STILTS), which is a command line scriptable version of topcat.  Please see http://www.star.bris.ac.uk/~mbt/stilts/sun256/sun256.html for installation instructions. 

**R packages**

Before running, several packages need to be installed on R.  This is very easy to do -- simply start R anywhere on your system and run the *install.packages()* command. Here is an example:

> $R  
> \> install.packages('argparser')


... and then follow the onscreen instructions. 

Please do this for the following packages:

*argparser*

**Preparing to run**

1. Mask the radio data based on the multi-wavelength data. I have not included this step since there are many ways to do this (e.g., Aladin's MOC, using a fits mask, ds9 regions, etc.). Your radio catalogue should not include any sources which are outside of or straddle the boundary for the multi-wavelength data, *including halo regions*. 

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


**Running the software**

1. Set up the configuration file with the correct parameters.
2. Run *make_master_cat.r*  to generate a slimmer version of the catalogue with added SNR information that will be used.
3. 
