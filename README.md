# R_likelihood_ratio
An R code for likelihood ratio matching between radio and multi-wavelength surveys.

**Software requirements**

R is a free, widely distributed software that is usually available on UNIX systems.  Please see https://www.r-project.org/ for more details. 

The software in this repository also makes use of the Starlink Tables Infrastructure Library Tool Set (STILTS), which is a command line scriptable version of topcat.  Please see http://www.star.bris.ac.uk/~mbt/stilts/sun256/sun256.html for installation instructions. 

**R packages**

Before running, several packages need to be installed on R.  This is very easy to do -- simply start R anywhere on your system and run the *install.packages()* command. Here is an example:

$R

>install.packages('argparser')


... and then follow the onscreen instructions. 

Please do this for the following packages:

*argparser*


**Running the software**

1. Set up the configuration file with the correct parameters.
2. Run *make_master_cat.r* to generate a slimmer version of the catalogue with added SNR information that will be used.
3. ....
