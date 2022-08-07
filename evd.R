#################################################
############## IMPORTANT REFERENCES #############
# http://quantdevel.com/BootstrappingTimeSeriesData/BootstrappingTimeSeriesData.pdf (About bootstrapping in time-series data).
# http://www.ccpo.odu.edu/~klinck/Reprints/PDF/omeyHUB2009.pdf (suggested by Suveges, 2014).

########### Resources for extreme value statistics ##########
# (1) http://personal.cityu.edu.hk/xizhou/first-draft-report.pdf
# (2) Playlist on Extreme Value Statistics: https://youtube.com/playlist?list=PLh35GyCXlQaTJtTq4OQGzMblwEcVIWW9n

#################################################

library('extRemes')
library('boot')
# library('GoFKernel')
library('cobs')
library('EnvStats')
source('BLS/bls.R')
source('TCF3.0/intf_libtcf.R')
source('standardize_periodogram.R')
source('test_periodograms.R')

statFunc <- function(a) { return (a) }

evd <- function(
    period,
    depth,
    duration,
    noiseType=1,  # Noise model present in y. Either 1 (white gaussian noise) or 2 (autoregressive noise). Resampling technique is dependent on this, see http://quantdevel.com/BootstrappingTimeSeriesData/BootstrappingTimeSeriesData.pdf
    # Note: noiseType is not used for adding noise to series, but instead used for deciding the way of resampling.
    algo="BLS",  # TODO: Add support for TCF
    ntransits=10
) {

    R <- 1000  # No. of bootstrap resamples of the original time series.
    K <- 100    # No. of distinct frequencies in a frequency bin.
    L <- 50    # No. of distinct frequency bins.

    # Generate light curve using the parameters.
    yt <- getLightCurve(period, depth, duration, noiseType=noiseType, ntransits=ntransits)
    y <- unlist(yt[1])
    t <- unlist(yt[2])

    # (1) Bootstrap the time series.
    if (noiseType == 1) {
        bootTS <- replicate(R, replicate(length(y), sample(y, 1, replace=TRUE)))
        bootTS <- aperm(bootTS)
        # bootTS will be of shape (R, length(y)).
    }
    else {  # TODO: (Q) Am I resampling with replacement in this case?
        bootTS <- tsboot(ts(y), statistic=statFunc, R=R, sim="fixed", l=length(y), n.sim=length(y))  # Moving-block bootstrap.
    }

    ### Create a frequency grid.
    perMin <- t[3] - t[1]
    perMax <- t[length(t)] - t[1]
    freqMin <- 1 / perMax
    freqMax <- 1 / perMin
    nfreq <- length(t) * 10
    freqStep <- (freqMax - freqMin) / nfreq
    freqGrid <- seq(freqMin, by=freqStep, length.out=nfreq)  # Goes from ~0.001 to 0.4999 (NOTE: Since delta_t = 1, fmax must be <= Nyquist frequency = 1/(2*delta_t) = 0.5 -- from Suveges, 2014).

    stopifnot(exprs={
        any(freqGrid) > 0.5  # No frequency must be greater than the Nyquist frequency.
    })

    # Divide the frequency into L bins, each with K datapoints.
    ## From https://stackoverflow.com/questions/57889573/how-to-randomly-divide-interval-into-non-overlapping-spaced-bins-of-equal-lengt
    intervalLength <- length(freqGrid)
    nBins <- L
    binWidth <- K
    binMinDistance <- 1
    spaceToDistribute <- intervalLength - (nBins * binWidth + (nBins - 1) * binMinDistance)
    distances <- diff(floor(c(0, sort(runif(nBins))) * spaceToDistribute))
    startOfBin <- cumsum(distances) + (0:(nBins-1)) * 101
    KLinds<- data.frame(bin = 1:nBins, startOfBin = startOfBin, endOfBin = startOfBin + binWidth - 1)

    print(KLinds)

    # TODO: Currently, both BLS and TCF use uniform frequency sampling. Ofir, 2014 suggested to use the optimal frequency sampling, so try to use it instead?

    # (2) Max of each partial periodogram
    # TODO: Need to use standardization/normalization somewhere! See astropy _statistics module under LombScargle to know at which step to normalize -- should we normalize these bootstrap periodograms or only the final full periodogram.
    maxOverAll_R_samples <- c()
    for (j in 1:R) {
        partialPeriodograms <- c()
        for (ll in 1:L) {
            freqs <- freqGrid[KLinds[ll,2]:KLinds[ll,3]]
            # TODO: Decide whether to use standardize or normal periodogram only for the partial ones?
            if (algo == "BLS") {
                partialPeriodogram <- bls(bootTS[j,], t, per.min=min(1/freqs), per.max=max(1/freqs), nper=K, bls.plot = FALSE)$spec
                # partialPeriodogram <- unlist(standardPeriodogram(bootTS[j,], t, perMin=min(1/freqs), perMax=max(1/freqs), nper=K, plot = FALSE, noiseType=noiseType)[1])
            }
            else {
                # For TCF, select K frequencies from freqs.
                freqsTCF <- seq(min(freqs), max(freqs), length.out=K)
                # print(freqsTCF)
                partialPeriodogram <- tcf(bootTS[j,], p.try = 1 / freqsTCF, print.output = FALSE)$outpow
                # partialPeriodogram <- unlist(standardPeriodogram(bootTS[j,], t, perMin=min(1/freqs), perMax=max(1/freqs), nper=K, plot = FALSE, noiseType=noiseType, algo="TCF")[1])
            }

            partialPeriodograms <- append(partialPeriodograms, partialPeriodogram)
        }
        maxOverAll_R_samples <- append(maxOverAll_R_samples, max(unlist(partialPeriodograms)))
    }
    print("Done calculating maxima...")
    print(maxOverAll_R_samples)
    # plot(maxOverAll_R_samples, type='l')

    # Decluster the peaks: https://search.r-project.org/CRAN/refmans/extRemes/html/decluster.html
    # Some intution on how to choose the threshold: https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.996.914&rep=rep1&type=pdf (search for threshold - Ctrl+F - in the paper)
    threshold <- quantile(maxOverAll_R_samples, 0.75)  # TODO: How to choose best threshold?
    maxOverAll_R_samples <- decluster(maxOverAll_R_samples, threshold = threshold)

    # (3) GEV modelling of partial maxima
    fitEVD <- fevd(maxOverAll_R_samples, type='GEV')
    # See https://www.dataanalysisclassroom.com/lesson60/ for discussion on the fevd function.
    print(summary(fitEVD))
    distill(fitEVD)

    # Diagnostic plots.
    # TODO: Why ci fails sometimes?
    try(plot(fitEVD))
    # plot(fitEVD, "trace")
    # return.level(fitEVD)
    # return.level(fitEVD, do.ci = TRUE)
    # ci(fitEVD, return.period = c(2, 20, 100))
    # See some description on how ci's are calculated: https://reliability.readthedocs.io/en/latest/How%20are%20the%20confidence%20intervals%20calculated.html

    # (4) Extrapolation to full periodogram
    print("Extrapolating to full periodogram...")
    print("Calculating FAP...")
    ## Get the parameters
    location <- findpars(fitEVD)$location[1]  # For some reason, the parameter values repeat 10 times, and all are same. So extract the first.
    scale <- findpars(fitEVD)$scale[1]
    shape <- findpars(fitEVD)$shape[1]

    ## Important note: It would be better to find an automatic way to judge whether we want to select a GEV model or not, instead of manually looking at the diagnostic plots. This is because we want to apply this method on several periodograms.

    # Compute full periodogram (note: standardized periodogram is used).
    # TODO: Since standardized periodogram's scale has changed (due to scatter-removal), it lies at the end of gev cdf, thus always giving fap=0.000 -- fix this: either remove the scatter or do some hackery to prevent this from happening.
    if (algo == "BLS") {
        op <- getStandardPeriodogram(period, depth, duration, noiseType=noiseType, algo=algo, ntransits=ntransits)
        output <- unlist(op[1])
        periodsTested <- unlist(op[2])
        periodEstimate <- periodsTested[which.max(output)]

        fullPeriodogramReturnLevel <- pgevd(max(output), location=location, scale=scale, shape=shape)
        print(sprintf("FAP (standardized periodogram): %f", nfreq * (1 - fullPeriodogramReturnLevel) / (K * L)))  # This formula is from Suveges, 2014.

        output <- bls(y, t, bls.plot = FALSE)$spec

        # Verify that the period corresponding to the largest peak in standardized periodogram is the same as in original periodogram.
        stopifnot(exprs={
            all.equal(periodEstimate, periodsTested[which.max(output)], tolerance = sqrt(.Machine$double.eps))
        })

        fullPeriodogramReturnLevel <- pgevd(max(output), location=location, scale=scale, shape=shape)
        print(sprintf("FAP (original periodogram): %f", nfreq * (1 - fullPeriodogramReturnLevel) / (K * L)))
    }
    else {
        op <- getStandardPeriodogram(period, depth, duration, noiseType=noiseType, algo=algo, ntransits=ntransits)
        output <- unlist(op[1])
        periodsTested <- unlist(op[2])
        periodEstimate <- periodsTested[which.max(output)]

        fullPeriodogramReturnLevel <- pgevd(max(output), location=location, scale=scale, shape=shape)
        print(sprintf("FAP (standardized periodogram): %f", nfreq * (1 - fullPeriodogramReturnLevel) / (K * L)))  # This formula is from Suveges, 2014.

        perMin = t[3] - t[1]
        perMax = t[length(t)] - t[1]
        freqMax = 1 / perMin
        freqMin = 1 / perMax
        nfreq = length(y) * 10
        freqStep = (freqMax - freqMin) / nfreq
        f = seq(freqMin, by=freqStep, length.out=nfreq)
        periodsToTry = 1 / f
        output <- tcf(y, p.try = periodsToTry, print.output = FALSE)$outpow

        stopifnot(exprs={
            all.equal(periodEstimate, periodsTested[which.max(output)], tolerance = sqrt(.Machine$double.eps))
        })

        fullPeriodogramReturnLevel <- pgevd(max(output), location=location, scale=scale, shape=shape)
        print(sprintf("FAP (original periodogram): %f", nfreq * (1 - fullPeriodogramReturnLevel) / (K * L)))
    }

    ###### Interpreting what FAP is good (from Baluev: https://academic.oup.com/mnras/article/385/3/1279/1010111):
    # > Given some small critical value FAP* (usually between 10−3 and 0.1), we can claim that the candidate signal is statistically
    # significant (if FAP < FAP*) or is not (if FAP > FAP*)
}

validate1_evd <- function(
    y,
    t,
    bootTS,
    R
) {
    for (j in 1:R) {
        for (i in 1:length(y)) {
            myVec <- c(bootTS$t[j,])
            stopifnot(exprs = {
                y[i] %in% myVec  # Obviously, values in the bootstrap sample must be there in the original time series since we are sampling from it.
            })
            any(duplicated(myVec))  # Fine if observations in the bootstrap resamples series duplicates.
        }
    }
}


################# Questions not yet understood by me ##################
# (1) What is "high quantiles of a distribution"? See online where mainly talk about heavy-tailed distributions..