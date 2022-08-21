#################################################
############## IMPORTANT REFERENCES #############
# http://quantdevel.com/BootstrappingTimeSeriesData/BootstrappingTimeSeriesData.pdf (About bootstrapping in time-series data).
# http://www.ccpo.odu.edu/~klinck/Reprints/PDF/omeyHUB2009.pdf (suggested by Suveges, 2014).
# See about aliasing at the end of this page, for example: https://docs.gammapy.org/0.8/time/period.html and this also: https://hea-www.harvard.edu/~swolk/thesis/period/node5.html
# See discussion on period/frequency spacing considerations for BLS: https://johnh2o2.github.io/cuvarbase/bls.html#period-spacing-considerations
# Mathematical description of the Anderson-Darling test: https://bookdown.org/egarpor/NP-UC3M/nptests-dist.html

########### Resources for extreme value statistics ##########
# (1) http://personal.cityu.edu.hk/xizhou/first-draft-report.pdf
# (2) Playlist on Extreme Value Statistics: https://youtube.com/playlist?list=PLh35GyCXlQaTJtTq4OQGzMblwEcVIWW9n
# https://www.lmd.ens.fr/E2C2/class/naveauRomaniaE2C207.pdf

#############################################################
# Good set of papers: https://arxiv.org/pdf/1712.00734.pdf

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
library('goftest')  # install.packages("goftest")

statFunc <- function(a) { return (a) }

evd <- function(
    period,
    depth,
    duration,
    noiseType=1,  # Noise model present in y. Either 1 (white gaussian noise) or 2 (autoregressive noise). Resampling technique is dependent on this, see http://quantdevel.com/BootstrappingTimeSeriesData/BootstrappingTimeSeriesData.pdf
    # Note: noiseType is not used for adding noise to series, but instead used for deciding the way of resampling.
    algo="BLS",
    ntransits=10,
    plot = TRUE,
    ofac=2,  # ofac is also called as "samples per peak" sometimes.
    useOptimalFreqSampling = FALSE,  # If want to use the optimal frequency sampling from Ofir, 2014: delta_freq = q / (s * os), where s is whole time series duration, os is oversampling factor and q is the duty cycle (time in single transit / total time series duration).
    alpha=0.05  # Significance level for hypothesis testing on the GEV fit on periodogram maxima. TODO: How to choose a significance level beforehand - any heuristics to follow?
) {

    R <- 1000  # No. of bootstrap resamples of the original time series.
    K <- 50   # No. of distinct frequencies in a frequency bin.  # TODO: Note that in Suveges, 2014, K = 16 is used and K is called as the oversampling factor - here we are not doing that, i.e. K is not the oversampling factor, ofac.
    L <- 100    # No. of distinct frequency bins.

    # Generate light curve using the parameters.
    yt <- getLightCurve(period, depth, duration, noiseType=noiseType, ntransits=ntransits)
    y <- unlist(yt[1])
    t <- unlist(yt[2])

    # (1) Bootstrap the time series.
    # The reason why we first bootstrap the time series and then take block maxima rather than simply bootstrapping block maxima of original series is mentioned in first paragraph in https://personal.eur.nl/zhou/Research/WP/bootstrap_revision.pdf
    # Non-parametric bootstrap with replacement of blocks.
    # if (noiseType == 1) {
    #     bootTS <- replicate(R, replicate(length(y), sample(y, 1, replace=TRUE)))
    #     bootTS <- aperm(bootTS)  # This just permutes the dimension of bootTS - rows become columns and columns become rows - just done for easier indexing further in the code.
    #     # At this point, bootTS will be of shape (R, length(y)).
    # }
    # Note that bootstrapping, by definition, is resampling "with replacement": https://en.wikipedia.org/wiki/Bootstrapping_(statistics)
    # We use block resampling irrespective of the noise (i.e. block resampling even if noise is uncorrelated in white Gaussian noise), because the underlying time-series is in the form of repeated box-like shapes, and we would like to preserve that in order to look more like the original time-series.
    # Note: A general option for time-series with independent data is to use: random, uniform bootstrapping, but that can distort the repeating box-like shapes in the time series.
    bootTS <- tsboot(ts(y), statistic=statFunc, R=R, sim="fixed", l=period*25, n.sim=length(y))$t  # block resampling with fixed block lengths
    # Here, the block length is chosen to be slightly larger than the period (so that each block atleast contains a period -- a heuristic).

    # This answer: https://stats.stackexchange.com/a/317724 seems to say that block resampling resamples blocks with replacement.

    ### Create a frequency grid.
    ################## Ofir, 2014 - optimal frequency sampling - notes #################
    # (1) "It is now easy to see that the frequency resolution ∆f is no longer constant - it depends on f itself due to the physics of the problem."
    # (2) Section 3.2 also says that by using a very fine frequency grid (suitable for long-period signals), we (a) increase computation time a lot, and (b) it will be too sensitive to noise and less to actual real signals.

    # In this code, there is also an option to use Ofir, 2014's suggestion to use the optimal frequency sampling rather than the default uniform frequency sampling.
    # Note that the fact that we uniformly sample in "frequency" rather than "period" is itself a good choice: see last para in sec 7.1 in https://iopscience.iop.org/article/10.3847/1538-4365/aab766/pdf
    # Note that while using min frequency as zero is often not a problem (does not add suprious peaks - as described in 7.1 in https://iopscience.iop.org/article/10.3847/1538-4365/aab766/pdf), here we start with min_freq = 1 / (duration of time series).
    # One motivation for oversampling (from https://iopscience.iop.org/article/10.3847/1538-4365/aab766/pdf): "...it is important to choose grid spacings smaller than the expected widths of the periodogram peaks...To ensure that our grid sufficiently samples each peak, it is prudent to oversample by some factor—say, n0 samples per peak--and use a grid of size 1 / (n0 * T)"
    # The above paper also says that n0 = 5 to 10 is common.
    perMin <- t[3] - t[1]
    perMax <- t[length(t)] - t[1]
    freqMin <- 1 / perMax
    freqMax <- 1 / perMin
    # nfreq <- length(t) * 10

    if (useOptimalFreqSampling) {
        if (algo == "BLS") {
            q = duration  # single transit duration / light curve duration.
        }
        else if (algo == "TCF") {  # TODO: This actually yields a very bad GEV fit - so something is wrong in calculating the duty cycle for TCF.
            # Duty cycle for TCF taken from Caceres, 2019 methodology paper: https://iopscience.iop.org/article/10.3847/1538-3881/ab26b8
            q = 1 / (period * 24)
        }
        s = length(t)
        freqStep = q / (s * ofac)
    }
    else {
        # Note: When we oversample, we are essentially imposing no constraints on the frequencies to be tested - that is helpful in general: see https://arxiv.org/pdf/1712.00734.pdf
        # Note that too much oversampling can lead to artifacts. These artifacts can be wrongly interpreted as a true periodic component in the periodogram.
        freqStep <- (freqMax - freqMin) / (nfreq * ofac)  # Oversampled by a factor, `ofac`.
    }

    freqGrid <- seq(from = freqMin, to = freqMax, by=freqStep)  # Goes from ~0.001 to 0.5 (NOTE: Since delta_t = 1, fmax must be <= Nyquist frequency = 1/(2*delta_t) = 0.5 -- from Suveges, 2014).
    sprintf("No. of frequencies in grid: %f", length(freqGrid))

    stopifnot(exprs={
        all(freqGrid <= 0.5)  # No frequency must be greater than the Nyquist frequency.
        length(freqGrid) >= K * L  # If this is not true, the FAP values could be greater than one, which is not realistic (this is my interpretation).
    })

    # Divide the frequency into L bins, each with K datapoints.
    ## From https://stackoverflow.com/questions/57889573/how-to-randomly-divide-interval-into-non-overlapping-spaced-bins-of-equal-lengt
    intervalLength <- length(freqGrid)
    nBins <- L
    binWidth <- K
    binMinDistance <- 1
    spaceToDistribute <- intervalLength - (nBins * binWidth + (nBins - 1) * binMinDistance)
    distances <- diff(floor(c(0, sort(runif(nBins))) * spaceToDistribute))
    startOfBin <- cumsum(distances) + (0:(nBins-1)) * 10
    KLinds <- data.frame(bin = 1:nBins, startOfBin = startOfBin, endOfBin = startOfBin + binWidth - 1)

    stopifnot(exprs={  # Check if the no. of frequencies in a bin is in fact equal to the desired number.
        length(freqGrid[KLinds[1, 2]:KLinds[1, 3]]) == binWidth
    })

    # (2) Max of each partial periodogram
    # TODO: Need to use standardization/normalization somewhere! See astropy _statistics module under LombScargle to know at which step to normalize -- should we normalize these bootstrap periodograms or only the final full periodogram.
    maxima_R <- c()
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
                partialPeriodogram <- tcf(diff(bootTS[j,]), p.try = 1 / freqsTCF, print.output = FALSE)$outpow
                # partialPeriodogram <- unlist(standardPeriodogram(bootTS[j,], t, perMin=min(1/freqs), perMax=max(1/freqs), nper=K, plot = FALSE, noiseType=noiseType, algo="TCF")[1])
            }

            # Note: If we use oversampling, then while it increases the flexibility to choose frequencies in the frequency grid, it also has important issues as noted in https://academic.oup.com/mnras/article/388/4/1693/981666:
            # (1) "if we oversample the periodogram, the powers at the sampled frequencies are no longer independent..."
            # To solve the above problem, we decluster the partial periodograms. Even without oversampling, the peaks tend to be clustered and we need to decluster the peaks.
            # TODO: See performance with and without declustering.
            partialPeriodogram <- decluster(partialPeriodogram, threshold = quantile(partialPeriodogram, probs=c(0.75)))

            partialPeriodograms <- append(partialPeriodograms, partialPeriodogram)
        }
        maxima_R <- append(maxima_R, max(unlist(partialPeriodograms)))
    }
    print("Done calculating maxima...")
    print(maxima_R)
    # plot(maxima_R, type='l')

    # Decluster the peaks: https://search.r-project.org/CRAN/refmans/extRemes/html/decluster.html
    # Some intution on how to choose the threshold: https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.996.914&rep=rep1&type=pdf (search for threshold - Ctrl+F - in the paper)
    # See section 5.3.2 in Coles, 2001 to see why declustering is needed: Extremes tend to cluster themselves and tend to occur in groups. Note that log-likelihood can be decomposed into a product of individual marginal distribution functions only under iid. So declustering "tries" to make them independent to prevent the violation of the iid assumption while fitting the GEV model below.
    # In short, declustering (approximately) solves the dependence issue of extremes.
    # TODO: Ensure where exactly to perform declustering -- over the periodogram itself or only on the extracted maxima? I think most likely the former.
    # TODO: We might not need declustering in all cases -- we can calculate the extremel index and do declustering only if index < 1...
    # Due to our way of extracting maxima of periodograms (i.e. not from whole periodogram but only from partial periodogram), maybe we do not even need declustering.
    # threshold <- quantile(maxima_R, 0.75)  # TODO: How to choose best threshold?
    # maxima_R <- decluster(maxima_R, threshold = threshold)

    # (3) GEV modelling of partial periodograms' maxima
    fitEVD <- fevd(maxima_R, type='GEV')
    # See https://www.dataanalysisclassroom.com/lesson60/ for discussion on the fevd function.
    print(summary(fitEVD))
    distill(fitEVD)

    ## Get the fitted GEV parameters
    location <- findpars(fitEVD)$location[1]  # For some reason, the parameter values repeat 10 times, and all are same. So extract the first.
    scale <- findpars(fitEVD)$scale[1]
    shape <- findpars(fitEVD)$shape[1]

    # Diagnostic goodness-of-fit tests (we use the Anderson-Darling (AD) test: https://search.r-project.org/CRAN/refmans/DescTools/html/AndersonDarlingTest.html)
    # A simple reason why we use the Anderson–Darling (AD) test rather than Komogorov-Smirnov (KS) is that AD is able to detect better the situations in which F0 and F differ on the tails (that is, for extreme data), where H0: F = F0 and H1: F \neq F0.
    result <- ad.test(maxima_R, null = "pgevd", location=location, scale=scale, shape=shape, nullname = "pgevd", estimated = FALSE)  # estimated = TRUE would have been fine as well since the gevd parameters (location, scale, shape) are estimated using the data itself - those three parameters are not data-agnostic. But here we use estimated = FALSE because using TRUE uses a different variant of AD test using the Braun's method which we do not want.
    print(result)
    print(sprintf("p-value for Anderson-Darling goodness-of-fit test of the periodogram maxima: %f", result$p.value))

    # Check if AD fit is good enough. If not, return a dummy fap value.
    # This check serves as a way to "automatically" find if the GEV fit is good and if it can be extrapolated to the full periodogram.
    # Suveges, 2014 suggests looking at the diagnostic plots before extrapolating to full periodogram, but that is cumbersome for large-scale simulations. Hence, this is a simple way to overcome manual fit quality inspection.
    # TODO: Make this work.
    # if (result$p.value < alpha) {  # Reject null hypothesis
    #     dummyFap <- -999
    #     return (dummyFap)
    # }

    # Diagnostic plots.
    if (plot) {
        # TODO: Why ci fails sometimes?
        try(plot(fitEVD))
        # plot(fitEVD, "trace")
        # return.level(fitEVD)
        # return.level(fitEVD, do.ci = TRUE)
        # ci(fitEVD, return.period = c(2, 20, 100))
        # See some description on how ci's are calculated: https://reliability.readthedocs.io/en/latest/How%20are%20the%20confidence%20intervals%20calculated.html
    }

    # # (4) Extrapolation to full periodogram
    # print("Extrapolating to full periodogram...")
    # print("Calculating FAP...")

    # ## Important note: It would be better to find an automatic way to judge whether we want to select a GEV model or not, instead of manually looking at the diagnostic plots. This is because we want to apply this method on several periodograms.

    # # Compute full periodogram (note: standardized periodogram is used).
    # # TODO: Since standardized periodogram's scale has changed (due to scatter-removal), it lies at the end of gev cdf, thus always giving fap=0.000 -- fix this: either remove the scatter or do some hackery to prevent this from happening.
    # if (algo == "BLS") {
    #     ## On standardized periodogram
    #     # op <- getStandardPeriodogram(period, depth, duration, noiseType=noiseType, algo=algo, ntransits=ntransits)
    #     # output <- unlist(op[1])
    #     # periodsTested <- unlist(op[2])
    #     # periodEstimate <- periodsTested[which.max(output)]

    #     # fullPeriodogramReturnLevel <- pgevd(max(output), location=location, scale=scale, shape=shape)
    #     # print(sprintf("FAP (standardized periodogram): %f", nfreq * (1 - fullPeriodogramReturnLevel) / (K * L)))  # This formula is from Suveges, 2014.

    #     ## On original periodogram
    #     output <- bls(y, t, bls.plot = FALSE)$spec
    #     print("max of output")
    #     print(max(output))

    #     print("Calculating return level corresponding to FAP = 0.01")
    #     ppgevd <- function(x) pgevd(x, location=location, scale=scale, shape=shape)
    #     ppgevdInv <- inverse(ppgevd)
    #     returnLevel <- ppgevdInv(1 - ((0.01 * K * L) / length(freqGrid)))
    #     print(returnLevel)

    #     # Verify that the period corresponding to the largest peak in standardized periodogram is the same as in original periodogram.
    #     # stopifnot(exprs={
    #     #     all.equal(periodEstimate, periodsTested[which.max(output)], tolerance = sqrt(.Machine$double.eps))
    #     # })
    #     ##### TODO: FAP calculation here is most likely not correct since here it is assumed that return level is same as full periodogram maxima, but it does not seem true.
    #     quantInv <- function(distr, value) ecdf(distr)(value)
    #     maxOutputIsWhichQuantile <- quantInv(maxima_R, max(output))

    #     fullPeriodogramValue <- pgevd(maxOutputIsWhichQuantile, location=location, scale=scale, shape=shape)
    #     fap <- nfreq * (1 - fullPeriodogramValue) / (K * L)
    #     print(sprintf("FAP (original periodogram): %f", fap))
    # }
    # else {
    #     # op <- getStandardPeriodogram(period, depth, duration, noiseType=noiseType, algo=algo, ntransits=ntransits)
    #     # output <- unlist(op[1])
    #     # periodsTested <- unlist(op[2])
    #     # periodEstimate <- periodsTested[which.max(output)]

    #     # fullPeriodogramValue <- pgevd(max(output), location=location, scale=scale, shape=shape)
    #     # print(sprintf("FAP (standardized periodogram): %f", nfreq * (1 - fullPeriodogramValue) / (K * L)))  # This formula is from Suveges, 2014.

    #     perMin = t[3] - t[1]
    #     perMax = t[length(t)] - t[1]
    #     freqMax = 1 / perMin
    #     freqMin = 1 / perMax
    #     nfreq = length(y) * 10
    #     freqStep = (freqMax - freqMin) / nfreq
    #     f = seq(freqMin, by=freqStep, length.out=nfreq)
    #     periodsToTry = 1 / f
    #     output <- tcf(diff(y), p.try = periodsToTry, print.output = FALSE)$outpow
    #     print("max of output")
    #     print(max(output))

    #     print("Calculating return level corresponding to FAP = 0.01")
    #     ppgevd <- function(x) pgevd(x, location=location, scale=scale, shape=shape)
    #     ppgevdInv <- inverse(ppgevd)
    #     print(1 - ((0.01 * K * L) / length(freqGrid)))
    #     returnLevel <- ppgevdInv(1 - ((0.01 * K * L) / length(freqGrid)))
    #     print(returnLevel)

    #     # stopifnot(exprs={
    #     #     all.equal(periodEstimate, periodsTested[which.max(output)], tolerance = sqrt(.Machine$double.eps))
    #     # })

    #     fullPeriodogramValue <- pgevd(max(output), location=location, scale=scale, shape=shape)
    #     fap <- nfreq * (1 - fullPeriodogramValue) / (K * L)
    #     print(sprintf("FAP (original periodogram): %f", fap))
    # }

    # return (fap);

    ###### Interpreting what FAP is good (from Baluev: https://academic.oup.com/mnras/article/385/3/1279/1010111):
    # (1) > Given some small critical value FAP* (usually between 10−3 and 0.1), we can claim that the candidate signal is statistically
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

smallestPlanetDetectableTest <- function(
    period,  # in days
    depths,  # in %
    duration,  # in hours
    algo,  # either BLS or TCF
    noiseType=1  # 1 for Gaussian and 2 for autoregressive noise
) {
    faps <- c()
    for (depth in depths) {
        fap <- evd(period, depth, duration, algo=algo, plot=FALSE)
        sprintf("depth (ppm): %f, fap: %f", depth*1e4, fap)
        faps <- append(faps, fap)
    }

    plot(depths*1e4, faps, log='y', xlab='Depth (ppm)', ylab='FAP', type='o')
    if (noiseType == 1) {
        abline(h=0.003, col='black', lty=2)  # FAP=0.003 corresponds to 3-sigma criterion for Gaussian -- commonly used in astronomy.
    }
    else {
        abline(h=0.002, col='black', lty=2)  # TODO: Decide what threshold FAP to use for the autoregressive case.
    }
}

################# Questions not yet understood by me ##################
# (1) What is "high quantiles of a distribution"? See online where mainly talk about heavy-tailed distributions..