#!/usr/bin/python
from astropy.io import fits
import numpy as np
import pandas
import argparse
import matplotlib
from matplotlib import pyplot as plt
from astropy.wcs import WCS
import scipy.stats as st
from scipy.optimize import curve_fit
import os

## define the model
def quad_model( x, a, b, c ):
    return a * x + b * x * x + c

def outlier_cut( x, a, b ):
    return a * x - b

def selection_function( x, x1, x2, y1, y2, a, b, c ):
    y = quad_model( x, a, b, c )
    y[np.where(x < x1)[0]] = y1
    y[np.where(x > x2)[0]] = y2
    return y

def fit_stellar_locus( xvals, yvals, iso_slope=0.6, iso_intersect=0.475, x1=0.15, x2=1.5, xmin=-0.5, xmax=2, ymin=-1.5, ymax=2, offset=0.14 ):

    ## limit the outliers in x and y
    x_idx = np.where( (xvals >= x1) & (xvals <= x2) )[0]
    y_idx = np.where( (yvals >= ymin) & (yvals <= ymax) )[0]    
    good_idx = np.intersect1d( x_idx, y_idx )

    ## apply selection to values to fit
    fit_xvals = xvals[good_idx]
    fit_yvals = yvals[good_idx]

    ## also limit the galaxy contribution
    x_div = np.arange( np.floor(xmin), np.ceil(xmax), 0.001 )
    y_div = iso_slope * x_div - iso_intersect

    ## find where the points are below this
    proj_y = outlier_cut( fit_xvals, iso_slope, iso_intersect )
    good_idx = np.where( fit_yvals <= proj_y )[0] 

    ## apply the selection
    fit_xvals = fit_xvals[good_idx]
    fit_yvals = fit_yvals[good_idx]

    result, rescov = curve_fit( quad_model, fit_xvals, fit_yvals )

    y_pred = quad_model( x_div, result[0], result[1], result[2] )
   
    ## final model
    level_idx = np.where( x_div < x1 )[0]
    y_pred[level_idx] = y_pred[np.max(level_idx)]
    start_idx = np.where( y_pred == np.max( y_pred ) )[0]
    y_pred[np.arange(start_idx,len(y_pred))] = np.max(y_pred)
    xmax_val = x_div[start_idx][0]

    ## division
    y_cut = y_pred + offset

    ## make a plot
    os.system( 'rm Stellar_locus_fit.png' )
    plt.plot( xvals, yvals, '.', markersize=1 )
    proj_y = outlier_cut( x_div, iso_slope, iso_intersect )
    plt.plot( x_div, proj_y, label='Initial cut' )
    plt.plot( x_div, y_pred, label='Fit to stellar locus' )
    plt.plot( x_div, y_cut, '-', label='Star/Galaxy separation' )
    plt.xlim( xmin,xmax )
    plt.ylim( ymin,ymax )
    plt.legend()
    plt.savefig('Stellar_locus_fit.png')
    plt.close(fig='all')
    print('Please check Stellar_locus_fit.png and run again if fit is not suitable.')
    cut_result = result
    cut_result[2] = cut_result[2] + offset

    print( 'Fit results to select galaxies: ' )
    print( ' y > {:f} for x < {:f}'.format( y_pred[np.max(level_idx)], x1 ) )
    print( ' y > {:f} * x + {:f} * x^2 + {:f} for {:f} < x < {:f}'.format( cut_result[0], cut_result[1], cut_result[2], x1, xmax_val ) )
    print( ' y > {:f} for x > {:f}'.format( np.max(y_cut), xmax_val ) )


    result = np.array([x1, xmax_val, y_pred[np.max(level_idx)], np.max(y_cut), cut_result[0], cut_result[1], cut_result[2]])

    return( result  )


def main( fits_cat, mask_name ):

    print( 'Reading in configuration file.' )
    config_params = pandas.read_table( config_file, delim_whitespace=True ).replace("'","",regex=True)
    outdir = config_params['value'][np.where( config_params['parameter'] == 'outdir' )[0][0]]
    bands = config_params['value'][np.where( config_params['parameter'] == 'bands' )[0][0]]
    my_bands = bands.split(',')
    id_col = config_params['value'][np.where( config_params['parameter'] == 'id_col' )[0][0]]
    ra_col = config_params['value'][np.where( config_params['parameter'] == 'ra_col' )[0][0]]
    dec_col = config_params['value'][np.where( config_params['parameter'] == 'dec_col' )[0][0]]
    mask_col = config_params['value'][np.where( config_params['parameter'] == 'mask_col' )[0][0]]
    sg_col = config_params['value'][np.where( config_params['parameter'] == 'sg_col' )[0][0]]
    mag_col = config_params['value'][np.where( config_params['parameter'] == 'mag_col' )[0][0]]
    mag_err_col = config_params['value'][np.where( config_params['parameter'] == 'mag_err_col' )[0][0]]
    flux_col = config_params['value'][np.where( config_params['parameter'] == 'flux_col' )[0][0]]
    flux_err_col = config_params['value'][np.where( config_params['parameter'] == 'flux_err_col' )[0][0]]
    beam_size = np.float(config_params['value'][np.where( config_params['parameter'] == 'beam_size' )[0][0]])
    cal_errors = np.float(config_params['value'][np.where( config_params['parameter'] == 'cal_errors' )[0][0]])

    hdul = fits.open( fits_cat )
    mytab = hdul[1].data

    tab_cols = mytab.columns
    colnames =  tab_cols.names

    ## open the mask    
    mask_fits = fits.open( mask_name )
    mask_wcs = WCS( mask_name )
    mask_wcs.all_world2pix(  )

    ## get RA, DEC pairs for the catalogue
    RA = mytab[ra_col]
    DEC = mytab[dec_col]
    pixx, pixy = mask_wcs.wcs_world2pix( RA, DEC, 1 )
    mask_vals = mask_fits[0].data
    mask_col_vals = [ mask_vals[np.int(y),np.int(x)] for x,y in zip(pixx, pixy) ]

    ## add this to the table
    tab_cols = tab_cols + fits.Column( name='Mask', format='I', array=mask_col_vals )


    ## get a list of bands
    mag_pref = mag_col.split('_')[0]
    mag_err_pref = mag_err_col.split('_')[0] + '_' + mag_err_col.split('_')[1]
    mybands = [ a for a in colnames if mag_pref in a and mag_err_pref not in a ]

    ## make SNR cuts for each
    for band in mybands:
        err_col = band.replace(mag_pref,mag_err_pref)
        snr_col = 1.08574 / mytab[err_col] ## for 5-sigma cut
        new_col = fits.Column( name = band.replace(mag_pref,'SNR'), format='D', array=snr_col )
        tab_cols = tab_cols + new_col

    ## add the various colour-colour columns we need    
    g_col = 'mag_CFHT-g'
    i_col = 'mag_CFHT-iy'
    r_col = 'mag_CFHT-r'
    j_col = 'mag_J'
    ks_col = 'mag_Ks'
    g_minus_i = mytab[g_col] - mytab[i_col]
    g_minus_r = mytab[g_col] - mytab[r_col]
    J_minus_Ks = mytab[j_col] - mytab[ks_col]

    ## add them to the table
    tab_cols = tab_cols + fits.Column( name='g_minus_i', format='D', array=g_minus_i )
    tab_cols = tab_cols + fits.Column( name='g_minus_r', format='D', array=g_minus_r )
    tab_cols = tab_cols + fits.Column( name='J_minus_Ks', format='D', array=J_minus_Ks )

    hdu = fits.BinTableHDU.from_columns(tab_cols)
    ## write the table
    hdu.writeto('tmp.fits')

    new_hdul = fits.open('tmp.fits')
    
    new_tab = new_hdul[1].data
    new_tab_cols = new_tab.columns

    ## apply mask
    unmasked_idx = np.where( new_tab[mask_col] == 1 )[0]
    unmasked_tab = new_tab[unmasked_idx]

    ## make 5-sigma selection
    g_idx = np.where( unmasked_tab[g_col.replace(mag_pref,'SNR')] >= 5. )[0]
    r_idx = np.where( unmasked_tab[r_col.replace(mag_pref,'SNR')] >= 5. )[0]
    J_idx = np.where( unmasked_tab[j_col.replace(mag_pref,'SNR')] >= 5. )[0]
    Ks_idx = np.where( unmasked_tab[ks_col.replace(mag_pref,'SNR')] >= 5. )[0]
    ## combine them
    idx = np.intersect1d( g_idx, r_idx )
    idx = np.intersect1d( idx, J_idx )
    idx = np.intersect1d( idx, Ks_idx )

    final_tab = unmasked_tab[idx] 

    ## make a plot
    xmin, xmax = -0.5, 2.
    ymin, ymax = -1.5, 2.

    xvals = final_tab['g_minus_r']
    yvals = final_tab['J_minus_Ks']

    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([xvals,yvals])
    kernel = st.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)

    mylevels = np.linspace(np.min(f),np.max(f),20)
    mylevels = mylevels[np.array([0,1,2,3,4,5,6,7,11,15,19])]

    fig = plt.figure()
    ax = fig.gca()
    ax.scatter(xvals, yvals, s=0.5, marker='.')
    ## overlay the kernel
    kde2d = ax.contour( xx, yy, f, mylevels )
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    fig.savefig('sg-sep.png')


    ## fit the stellar locus
    ## PLEASE NOTE: YOU MAY HAVE TO ADJUST THE VALUES GIVEN TO THE INITIAL OUTLIER CUTS AND FITS
    result = fit_stellar_locus( final_tab['g_minus_r'], yvals = final_tab['J_minus_Ks'], iso_slope=0.59, iso_intersect=0.55, x1=0.175, x2=1.5, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, offset=0.14 )

    ## make a star/galaxy separation column
    ## 1 = galaxy, 0 = star
    predict_y = selection_function( new_tab['g_minus_r'], result[0], result[1], result[2], result[3], result[4], result[5], result[6] )
    galaxy_idx = np.where( new_tab['J_minus_Ks'] >= predict_y )
    sg_sep = np.zeros(len(new_tab))
    sg_sep[galaxy_idx] = 1

    new_tab_cols = new_tab_cols + fits.Column( name='sg_sep', format='I', array=sg_sep )

    final_hdu = fits.BinTableHDU.from_columns(new_tab_cols)
    ## write the table
    final_hdu.writeto( fits_cat.replace('.fits','_SG_sep.fits') )
    os.system('rm tmp.fits')

        




