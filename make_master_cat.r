library('argparser')


main <- function( master_cat, config_file='lr_config.txt', stilts_cmd='stilts' ){

	## read the configuration file
	cat( 'Reading configuration file ...' )
	A <- read.table( config_file, stringsAsFactors=FALSE, header=FALSE )
	outdir <- A$V2[which(A$V1 == 'outdir')]
	bands <- A$V2[which(A$V1 == 'bands')]
	my_bands <- strsplit(bands, ',')[[1]]
	id_col <- A$V2[which(A$V1 == 'id_col')]
	ra_col <- A$V2[which(A$V1 == 'ra_col')]
	dec_col <- A$V2[which(A$V1 == 'dec_col')]
	halo_col <- A$V2[which(A$V1 == 'halo_col')]
	sg_col <- A$V2[which(A$V1 == 'sg_col')]
	mag_col <- A$V2[which(A$V1 == 'mag_col')]
	mag_err_col <- A$V2[which(A$V1 == 'mag_err_col')]
	cat( ' configuration complete.\n' )


	## basic columns to keep
	basic_cols <- c( id_col, ra_col, dec_col, halo_col, sg_col )

	## magnitude columns
	keep_bands <- c()
	keep_bands_err <- c()
	for ( my_band in my_bands ){
		keep_bands <- c( keep_bands, gsub( 'X', my_band, mag_col ) )
		keep_bands_err <- c( keep_bands_err, gsub( 'X', my_band, mag_err_col ) )
	}

	## prepare stilts commands
	cat( 'Preparing STILTS command ... \n' )
	## calculate SNR
	snr_cmd <- c()
	snr_cols <- c()
	for ( xx in 1:length(keep_bands) ){
		
		colname <- paste( keep_bands[xx], '_SNR', sep='' )
		snr_cols <- c( snr_cols, colname )
		snr_cmd <- c( snr_cmd, paste( " cmd=\'addcol ", colname, " 1.08574/", keep_bands_err[xx], "\'",  sep='' ) )
	}

	## update the keep columns
	keep_cols <- c( basic_cols, keep_bands, keep_bands_err, snr_cols )

	outtmp <- strsplit( master_cat, '/' )[[1]]
	outcat <- paste( outdir, outtmp[length(outtmp)], sep='/' )
	outcat <- gsub('.fits', '_master.fits', outcat )
	
	## run stilts
	ss <- paste( stilts_cmd, " tpipe in=", master_cat, " out=", outcat, " omode=out ", paste( snr_cmd, collapse='' ), ' cmd=\'keepcols "', paste( keep_cols, collapse=' ' ), '"', "\' ", sep='' )
	cat( 'Running the following command:\n' )
	cat( ss )
	cat( '\n' )
	system( ss )
	cat( 'Master cat finished: ', outcat, '\n' )

}

p <- arg_parser("prepare catalogue for likelihood ratio fitting")
p <- add_argument( p, "master_cat", help="base catalogue to start with" )
p <- add_argument( p, "--config_file", help="configuration file (see example)" )
p <- add_argument( p, "--stilts_cmd", help="executable for STILTS if not stilts", default="stilts" )

argv <- parse_args( p )

main( argv$master_cat, config_file=argv$config_file, stilts_cmd=argv$stilts_cmd )

