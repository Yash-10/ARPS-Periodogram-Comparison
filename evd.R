#################################################
############## IMPORTANT REFERENCES #############
# http://quantdevel.com/BootstrappingTimeSeriesData/BootstrappingTimeSeriesData.pdf (About bootstrapping in time-series data).
# http://www.ccpo.odu.edu/~klinck/Reprints/PDF/omeyHUB2009.pdf (suggested by Suveges, 2014).

library('extRemes')
library('boot')
source('BLS/bls.R')

statFunc <- function(a) { return (a) }

evd <- function(
    y,
    t
) {

    R <- 10  # No. of bootstrap resamples of the original time series.
    K <- 100    # No. of distinct frequencies in a frequency bin.
    L <- 80    # No. of distinct frequency bins.

    # (1) Bootstrap the time series.
    # TODO: (Q) Am I resampling with replacement??
    bootTS <- tsboot(y, statistic=statFunc, R=R, sim="fixed", l=length(y), n.sim=length(y))

    # plot(bootTS$t[1,], col='red', type='l')
    # lines(bootTS$t[2,], col='black', type='l')

    ### Create a frequency grid.
    perMin <- t[3] - t[1]
    perMax <- t[length(t)] - t[1]
    freqMin <- 1 / perMax
    freqMax <- 1 / perMin
    nfreq <- length(t) * 10
    freqStep <- (freqMax - freqMin) / nfreq
    freqGrid <- seq(freqMin, by=freqStep, length.out=nfreq)  # Goes from ~0.001 to 0.4999 (NOTE: Since delta_t = 1, fmax must be <= Nyquist frequency = 1/(2*delta_t) = 0.5 -- from Suveges, 2014).

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

    # (2) Max of each partial periodogram
    maxOverAll_R_samples <- c()
    for (j in 1:R) {
        maxOfPartialPeriodograms <- c()
        for (ll in 1:L) {
            freqs <- freqGrid[KLinds[ll,2]:KLinds[ll,3]]
            partialPeriodogram <- bls(bootTS$t[j,], t, per.min=max(1/freqs), per.max=min(1/freqs), bls.plot = FALSE)
            maxOfPartialPeriodograms <- append(maxOfPartialPeriodograms, max(partialPeriodogram$spec))
        }
        maxOverAll_R_samples <- append(maxOverAll_R_samples, max(maxOfPartialPeriodograms))
    }
    print("Done calculating maxima...")
    print(maxOverAll_R_samples)
    # plot(maxOverAll_R_samples, type='l')

    # Decluster the peaks: https://search.r-project.org/CRAN/refmans/extRemes/html/decluster.html
    # TODO

    # (3) GEV modelling of partial maxima
    fitEVD <- fevd(maxOverAll_R_samples)
    distill(fitEVD)

    # Diagnostic plots.
    plot(fitEVD)
    # plot(fitEVD, "trace")
    # return.level(fitEVD)
    # return.level(fitEVD, do.ci = TRUE)
    # ci(fitEVD, return.period = c(2, 20, 100))
}
