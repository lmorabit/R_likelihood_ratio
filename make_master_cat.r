library('argparser')


main <- function( master_cat, outdir='/mnt/vardy/leah/likelihood_ratio/', my_bands=c('J','H','Ks'), band_suffix='', id_col='UID', ra_col='RAcen', dec_col='Deccen', halo_col='Mask', sgcol='SGsep', mag_col_to_use='mag', mag_err_col_to_use='mag_err' ){


	## basic columns to keep
	basic_cols <- c( id_col, ra_col, dec_col, halo_col, sgcol )

	## magnitudes -- may need to paste these in a different order
	keep_bands <- paste( mag_col_to_use, '_', my_bands, band_suffix, sep='' )
	keep_band_errs <- paste( mag_err_col_to_use, '_', my_bands, band_suffix, sep='' )

	## prepare stilts commands
	## calculate SNR
	snr_cmd <- c()
	snr_cols <- c()
	for ( xx in 1:length(keep_bands) ){
		
		colname <- paste( keep_bands[xx], '_SNR', sep='' )
		snr_cols <- c( snr_cols, colname )
		snr_cmd <- c( snr_cmd, paste( " cmd=\'addcol ", colname, " 1.08574/", keep_band_errs[xx], "\'",  sep='' ) )
	}

	## update the keep columns
	keep_cols <- c( basic_cols, keep_bands, keep_band_errs, snr_cols )

	outtmp <- strsplit( master_cat, '/' )[[1]]
	outcat <- paste( outdir, outtmp[length(outtmp)], sep='' )
	outcat <- gsub('.fits', '_master.fits', outcat )
	
	## run stilts
	ss <- paste( "stilts tpipe in=", master_cat, " out=", outcat, " omode=out ", paste( snr_cmd, collapse='' ), ' cmd=\'keepcols "', paste( keep_cols, collapse=' ' ), '"', "\' ", sep='' )
	system( ss )

}

p <- arg_parser("prepare catalogue for likelihood ratio fitting")
p <- add_argument( p, "master_cat", help="base catalogue to start with" )
p <- add_argument( p, "--outdir", help="out directory to store the processed catalogue", default="/mnt/vardy/leah/likelihood_ratio/", type="string" )
p <- add_argument( p, "--mybands", help="bands to use for LR matching", default=c('J','H','Ks'))
p <- add_argument( p, "--band_suffix", help="suffix for bands -- implies aperture", default="", type="string" )
p <- add_argument( p, "--idcol", help="ID column name", default="UID", type="string" )
p <- add_argument( p, "--racol", help="RA column name", default="RAcen", type="string" )
p <- add_argument( p, "--deccol", help="Dec column name", default="Deccen", type="string" )
p <- add_argument( p, "--halo", help="Halo Mask Column", default="Mask", type="string" )
p <- add_argument( p, "--sgsep", help="star galaxy separation Column", default="SGsep", type="string" )
p <- add_argument( p, "--magcol", help="mag column name", default="mag", type="string" )
p <- add_argument( p, "--magerrcol", help="mag error column name", default="mag_err", type="string" )		  

argv <- parse_args( p )

main( argv$master_cat, outdir=argv$outdir, my_bands=argv$mybands,  band_suffix=argv$band_suffix, id_col=argv$idcol, ra_col=argv$racol, dec_col=argv$deccol, halo_col=argv$halo, sgcol=argv$sgsep, mag_col_to_use=argv$magcol, mag_err_col_to_use=argv$magerrcol )

