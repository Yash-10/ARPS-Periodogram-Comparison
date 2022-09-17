# tcffit.R
# Part of KARPS pipeline, flow.control bit 2=P
# Construct Transit Comb Filter periodogram of ARMA-type
# residuals for Kepler lightcurves

# Slightly adapted by Eric Feigelson for KARPS pipeline from
# karps_analysis.R by Gabriel Caceres
# Oct 14 2015


# The following may be needed if the pipeline modules are run out of order

#if(exists('Nkid') == FALSE | exists('results_file_name')  == FALSE) {
	good_list_name <- paste0(dir.path, pipeline.name, '_good.lst')
	kid <- as.character(read.table(good_list_name)$V1)
	Nkid <- length(kid)
	leading_zeros <- NULL
	for (i in 1:Nkid) {
	    leading_zeros <- paste(rep(0,(9 - nchar(kid[i]))), sep="", collapse="")
	    kid[i] <- paste0(leading_zeros, kid[i])
	}
	results_file_name <- paste0(dir.path, 'results/', kid, '_results.Rdat')
	if(verbose > 0) cat('Recovered Nkid = ', Nkid, ' KID names from file ', 
		good_list_name, '\n')
#} 

# Preliminaries:
#	- Skip this pipeline module when 'results' list has the 'tcffit' element
#	- Get ~500,000 test periods for periodogram
#	- Get R function 'tcf' that contains .Fortran call to 'main_tcf' for TCF calculation
#	NOTE: At present, the R function 'tcf' is available only within the astro.psu.edu domain

skip_tcffit_ind <- vector(length=Nkid)
for(i in 1:Nkid) {
    load(file=results_file_name[i]) # contents in 'results'
    if(length(results$tcffit) != 0) skip_tcffit_ind[i] <- TRUE
}
test_periods <- readRDS("/bulk/binah2/edf/KARPS/KARPS_pipeline/src/periods_to_test.rds")
source("/bulk/daat2/gac158/kepler_AR/analysis/TCF3.0/intf_libtcf.R") ## TCF 3.0

# Loop over each KID star: load ARIMA & ARFIMA residuals and call tcf function

for (i in 1:Nkid) {
	if(skip_tcffit_ind[i] == TRUE) next
        load(file=results_file_name[i]) # contents in 'results'

	for (k in 1:2) {            # loop over ARIMA and ARFIMA type of TCF periodogram
	    	if(k==1) {
			resids <- results$arfit$arima_fit$residuals
	    		tcf_arima <- tcf(resids, p.try = test_periods, print.output <- FALSE)
		}
		if(k==2) {
			resids <- results$arfit$arfima_fit$residuals
			lc_nas <- results$arfit$arfima_fit$x == 0.0
			resids[lc_nas] <- NA   ### REPLACE WITH NA FOR ARFIMA 
	    		tcf_arfima <- tcf(resids, p.try = test_periods, print.output <- FALSE)
		}
	}
	results[['tcf']] <- list(test_periods = test_periods, tcf_arima = tcf_arima, tcf_arfima = tcf_arfima)
	save(results, file=results_file_name[i])
	if(verbose>0) cat('TCF calculation for KID ', kid[i], ' number ', i, ' is completed ', date(),  '\n') 
}
