import argparse
import glob
from astropy.table import Table, Column, vstack
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

def check_duplicates( mytab ):

    ## check for duplicates
    unique_vals, unique_idx, dup_count = np.unique( mytab['radio_ID'], return_index=True, return_counts=True )
    if len( unique_vals ) == len( mytab ):
        print( 'Each match is unique.' )
        return mytab
    else:
        print( 'Some radio sources have more than one valid match, selecting based on highest Rel value.' )
        unique_tab = mytab[unique_idx]
        dup_idx = np.where( dup_count > 1 )[0]
        ## get the rel col
        rel_col = [ xx for xx in tmp_vals.colnames if 'Rel' in xx ]
        rel_col = rel_col[0]
        for xx in dup_idx:
            dup_val = unique_vals[xx]
            tmp_vals = mytab[np.where( mytab['radio_ID'] == dup_val)[0]]
            best_val = tmp_vals[np.where(tmp_vals[rel_col] == np.max(tmp_vals[rel_col]))[0]]
            unique_tab[np.where(unique_tab['radio_ID'] == dup_val)[0]] = best_val
        return unique_tab


def main( ordered_bands='Ks,H,J', LR_threshold = 0.8 ):

    my_bands = ordered_bands.split(',')
    myfiles =glob.glob( my_bands[0]+'*LR_matches.dat' )
    ## read in the primary band
    if len( myfiles ) > 0:
        mytab = Table.read( glob.glob( my_bands[0]+'*LR_matches.dat' )[0], format='ascii' )
    else:
        print( 'Files not found, please check you are running in the right directory.' )
        return

    ## filter by LR_threshold
    rel_col = my_bands[0] + '_Rel'
    LR_idx = np.where( mytab[rel_col] > LR_threshold )[0]
    print( 'There are %s matches in the primary band, %s.'%(str(len(LR_idx)), my_bands[0] ) )

    good_matches = mytab[LR_idx]
    ## for more than one match, keep one with highest rel
    good_matches = check_duplicates( good_matches )

    ## clean up the table
    for tmp_col in good_matches.colnames:
        new_col = tmp_col.replace( my_bands[0]+'_', '' )
        good_matches[tmp_col].name = new_col
        
    ## add a column for the match band
    good_matches.add_column( Column( name='Band', data=np.repeat(my_bands[0], len(good_matches) ) ) )

    ## for the rest of the bands, read in the data and look for unique IDs
    for my_band in my_bands[1:len(my_bands)]:

        print( 'Checking band %s ...'%my_band )
        tmp_tab = Table.read( glob.glob( my_band+'*LR_matches.dat' )[0], format='ascii' )
        rel_col = my_band + '_Rel'
        LR_idx = np.where( tmp_tab[rel_col] > LR_threshold )[0] 
        tmp_matches = tmp_tab[LR_idx]
        ## for more than one match, keep one with highest rel
        tmp_matches = check_duplicates( tmp_matches )
        ## find the radio sources which do not already have a match
        rad_idx = [ xx for xx, id in enumerate(tmp_matches['radio_ID']) if id not in good_matches['radio_ID'] ]
        print( 'There are %s unique matches found in band %s.'%(str(len(rad_idx)), my_band ) )
        new_tab = tmp_matches[rad_idx]

        ## clean up the table
        for tmp_col in new_tab.colnames:
            new_col = tmp_col.replace( my_band+'_', '' )
            new_tab[tmp_col].name = new_col

        new_tab.add_column( Column( name = 'Band', data=np.repeat( my_band, len(new_tab) ) ) )
        
        good_matches = vstack( [good_matches, new_tab] )

    ## save the table
    print( 'Writing final matches to Final_matches.dat and plotting' )
    good_matches.write( 'Final_matches.dat', format='ascii' )

    ## make a plot
    ## load a cmap
    cmap = matplotlib.cm.get_cmap('magma')
    cvals = np.linspace( 0, 1, num=len(my_bands)+2 )
    fig, axs = plt.subplots( 2, sharex=True, gridspec_kw={'hspace':0} )
    for xx in np.arange(len(my_bands)):
        band_idx = np.where( good_matches['Band'] == my_bands[xx] )[0]
        axs[0].scatter( good_matches['separation'][band_idx], good_matches['Rel'][band_idx], marker='.', facecolor=cmap(cvals[xx+1]), edgecolor='0.3', label=my_bands[xx] )
        axs[1].scatter( good_matches['separation'][band_idx], good_matches['LR'][band_idx], marker='.', facecolor=cmap(cvals[xx+1]), edgecolor='0.3', label=my_bands[xx] )

    axs[0].set( ylabel='Reliability' )
    axs[1].set( xlabel='separation [arcsec]', ylabel='LR value' )
    axs[0].legend()
    axs[1].legend()
    axs[1].set( ylim=(0,3.*np.median(good_matches['LR']) ) )
    fig.savefig( 'Final_matches.png' )
    fig.clear()

    print( 'Done.' )

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--ordered_bands', type=str, help='string of comma separated bands in order of priority' )
    parser.add_argument('--LR_threshold', type=float, help='Value of LR reliability threshold [default 0.8]', default=0.8 )
    args = parser.parse_args()

    main( ordered_bands = args.ordered_bands, LR_threshold = args.LR_threshold )
