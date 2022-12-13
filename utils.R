## Utility functions.
library(forecast)
library(reticulate)

source_python("python_utils.py")

getResidForTCF <- function(
    y  # Time series (must not be differenced because it is done internally).
) {
    max.p = 5
    max.q = 5
    max.d = 0
    ARIMA.fit = auto.arima(diff(y), stepwise=FALSE, approximation=FALSE, seasonal=FALSE, max.p=max.p, max.q=max.q, max.d=max.d) #leave d as 0. 
    # Simple statistics of ARIMA residuals 
    ARIMA.resid = residuals(ARIMA.fit)
    return (ARIMA.resid)
}

plot_folded_lc <- function(bestper, harmonic=1, phase_offset=0, info=T, depth=NA, p_se =F)
{  
    TCF.fold.Phase = (lc[,1] - min(lc[,1], na.rm=TRUE)) %% bestper / bestper
    Range_lc = 1.3*diff(range(lc[,2], na.rm=TRUE)) + min(lc[,2], na.rm=TRUE)
    Text0.pos = 1.35*diff(range(lc[,2], na.rm=TRUE)) + min(lc[,2], na.rm=TRUE)
    plot(TCF.fold.Phase, lc[,2], xlab = "Phase", ylab ="Normalized Flux", cex.axis=1.0, cex.lab=1.2,
        cex = 0.5, pch = 20, col="#00000080", ylim=c(min(lc[,2], na.rm=TRUE), Range_lc), xaxs='i')
    axis(1,at=seq(from=0, to=1, by=.1), labels = F, tcl=0.25)
    axis(1,at=seq(from=0, to=1, by=.2), labels = F, tcl=0.5)
    text(0.5, Text0.pos, "Phase-folded light curve", pos=1, cex=1.2)
}

getFreqGridToTest <- function(
    t,  # Observation epochs.
    period, duration,
    res=2,  # This the resolution in the time series that controls the candence. res=2 means cadence of 30 min and res=1 means 1hr, for example.
    ofac=1,  # Oversampling factor for frequency selection.
    useOptimalFreqSampling=FALSE,
    algo="BLS"
) {
    ### Create a frequency grid.
    ## Ofir, 2014 - optimal frequency sampling - notes ##
    # (1) "It is now easy to see that the frequency resolution ∆f is no longer constant - it depends on f itself due to the physics of the problem."
    # (2) Section 3.2 also says that by using a very fine frequency grid (suitable for long-period signals), we (a) increase computation time a lot, and (b) it will be too sensitive to noise and less to actual real signals.

    # In this code, there is also an option to use Ofir, 2014's suggestion to use the optimal frequency sampling rather than the default uniform frequency sampling.
    # Note that the fact that we uniformly sample in "frequency" rather than "period" is itself a good choice: see last para in sec 7.1 in https://iopscience.iop.org/article/10.3847/1538-4365/aab766/pdf
    # Note that while using min frequency as zero is often not a problem (does not add suprious peaks - as described in 7.1 in https://iopscience.iop.org/article/10.3847/1538-4365/aab766/pdf), here we start with min_freq = 1 / (duration of time series).
    # One motivation for oversampling (from https://iopscience.iop.org/article/10.3847/1538-4365/aab766/pdf): "...it is important to choose grid spacings smaller than the expected widths of the periodogram peaks...To ensure that our grid sufficiently samples each peak, it is prudent to oversample by some factor—say, n0 samples per peak--and use a grid of size 1 / (n0 * T)"
    # The above paper also says that n0 = 5 to 10 is common.

    # Observation: TCF peak changes as res changes. In fact, the estimated period scales by `res`. This is expected as per our conversation. So we need to scale the estimated period with res.
    # I think the reason for that could be that bls takes `t` (apart from `y`, the observation values), so it knows the time-spacing internally whereas TCF does not take t.
    perMin <- t[2*res+1] - t[1]  # `res` is needed because t[3] - t[1] = 2 time units for res=1 but 1 for res=2. So t[3] - t[1] would become smaller and smaller for larger res, which should not happen since it will give errors at some point.
    perMax <- t[length(t)] - t[1]
    freqMin <- 1 / perMax
    freqMax <- 1 / perMin
    nfreq <- length(t) * 10  # This particular value is taken from BLS - see bls.R. Here, it is used for both BLS and TCF.

    if (useOptimalFreqSampling) {
        if (algo == "BLS") {
            q = duration / (period * 24)  # single transit duration / period, both in same units.
        }
        else if (algo == "TCF") {  # TODO: Confirm duty cycle scaling by res.
            # Duty cycle for TCF taken from Caceres, 2019 methodology paper: https://iopscience.iop.org/article/10.3847/1538-3881/ab26b8
            # Again, scaling by `res` is performed because TCF calculates periodogram in terms of cadence rather than absolute times.
            q = res * 1 / (period * 24)
        }
        s = length(t)
        freqStep = q / (s * ofac)
    }
    else {
        # Note: When we oversample, we are essentially trying to impose no constraints on the frequencies to be tested - that is helpful in general: see https://arxiv.org/pdf/1712.00734.pdf
        # Note that too much oversampling can lead to artifacts. These artifacts can be wrongly interpreted as a true periodic component in the periodogram.
        freqStep <- (freqMax - freqMin) / (nfreq * ofac)  # Oversampled by a factor, `ofac`.
    }

    freqGrid <- seq(from = freqMin, to = freqMax, by = freqStep)  # Goes from ~0.001 to max freq set by the time spacing (NOTE: fmax must be <= Nyquist frequency = 1/(2*delta_t) -- from Suveges, 2014), where delta_t here is res.
    freqGrid <- freqGrid[-length(freqGrid)]  # Remove the last frequency. This is done to prevent periodogram returning nan power at frequency = Nyquist frequency.
    print(sprintf("No. of frequencies in grid: %f", length(freqGrid)))

    return (freqGrid);
}

calculateSNR <- function(  # TODO: For making this more efficient, compute trend fit only for region around the periodogram peak - will save some time.
    periods,  # The periods at which the periodogram is computed.
    periodogramPower,  # This power must not be detrended since it is done internally.
    lambdaTrend=1
) {
    cobsTrend <- cobs(log10(periods), periodogramPower, ic='BIC', tau=0.5, lambda=lambdaTrend)
    detrended <- cobsTrend$resid
    # return (list(periods, detrended, cobsTrend$fitted, periodogramPower))
    lowerInd <- which.max(detrended) - 100
    upperInd <- which.max(detrended) + 100
    if (which.max(detrended) - 100 < 1) {
        lowerInd <- 1
    }
    if (which.max(detrended) + 100 > length(detrended)) {
        upperInd <- length(detrended)
    }
    consider <- detrended[lowerInd:upperInd]
    snr <- max(consider) / IQR(consider)
    return (snr)
}

# divide <- function(vec, n, min_spacing = 1) {  # Below approach to ensure the selected L frequencies are spaced by atleast the oversampling factor is taken from https://stackoverflow.com/a/66036847:
#     unname(tapply(vec, ceiling(seq_along(vec) / min_spacing), sample, size = 1))
#     # head(u[seq(1, length(u), by = 2)], n)
# }

# divide <- function(vec, n, min_spacing = 2) {
#   idx <- seq_along(vec)
#   repeat {
#     k <- sort(sample(idx,n))
#     if (all(diff(k)>=min_spacing)) break
#   }
#   vec[k]
# }

freqdivideFreqGrid <- function(freqGrid, L, K, seedValue=1) {
    set.seed(seedValue)  # Set seed for reproducibility.

    if ((K %% 2) == 0) {
        safeDist <- 1 + K/2  # 1 is added just to be more safe at the edges of the frequency grid. This is just a hackery.
    }
    else {
        safeDist <- 1 + (K-1)/2  # 1 is added just to be more safe at the edges of the frequency grid. This is just a hackery.
    }
    endIndex <- length(freqGrid) - safeDist
    freqConsider <- freqGrid[safeDist:endIndex]
    KLfreqs <- rand_parts(freqConsider, n=L, l=K)
    return (unlist(KLfreqs))
}


findAGoodOfac <- function() {  # This function is just a quick, rough method, not intended to use for sophisticated analyses. 
    hwhms <- c()
    things<-list(c(3, 0.01, 1/36), c(3, 0.005, 1/36), c(5, 0.01, 1/60), c(5, 0.005, 1/60), c(7, 0.01, 1/84), c(7, 0.005, 1/84), c(11, 0.005, 1/132))
    for (thing in things) {
        yt <- getLightCurve(thing[1], thing[2], thing[3], noiseType=1, ntransits=10, gaussStd=gaussStd, ar=ar, ma=ma, order=order, res=res, checkConditions=TRUE)
        y <- unlist(yt[1])
        t <- unlist(yt[2])
        freqGrid <- getFreqGridToTest(t, res=2, ofac=1, useOptimalFreqSampling=FALSE, algo="TCF")
        output <- bls(y, t, bls.plot = FALSE, per.min=min(1/freqGrid), per.max=max(1/freqGrid), nper=length(freqGrid))

        # fstep <- (max(freqGrid) - min(freqGrid)) / length(freqGrid)
        # freqs <- seq(from = min(freqGrid), by = fstep, length.out = length(freqGrid))
        # periodsToTry <- 1 / freqs
        # residTCF <- getResidForTCF(y)
        # output <- tcf(residTCF, p.try = periodsToTry * res, print.output = TRUE)

        out<-output$spec
        ptest<-output$periodsTested
        i<-which.max(out)
        l<-i-100
        u<-i+100
        plot(1/ptest[l:u], out[l:u])
        hwhm_ <- fwhm(1/ptest[l:u], out[l:u])/2
        hwhms <- append(hwhms, hwhm_)
        print("hwhm")
        print(hwhm_)
    }
    print(hwhms)
    print("mean")
    print(mean(hwhms))
}

### Store ###
    # This is the earlier code (that did not work well), to divide the frequency grid into K*L frequencies.
    # # Divide the frequency into L bins, each with K datapoints.
    # ## From https://stackoverflow.com/questions/57889573/how-to-randomly-divide-interval-into-non-overlapping-spaced-bins-of-equal-lengt
    # intervalLength <- length(freqGrid)
    # nBins <- L
    # binWidth <- K
    # binMinDistance <- 1
    # spaceToDistribute <- intervalLength - (nBins * binWidth + (nBins - 1) * binMinDistance)
    # distances <- diff(floor(c(0, sort(runif(nBins))) * spaceToDistribute))
    # startOfBin <- cumsum(distances) + (0:(nBins-1)) * 101
    # KLinds <- data.frame(bin = 1:nBins, startOfBin = startOfBin, endOfBin = startOfBin + binWidth - 1)

    # stopifnot(exprs={  # Check if the no. of frequencies in a bin is in fact equal to the desired number.
    #     length(freqGrid[KLinds[1, 2]:KLinds[1, 3]]) == binWidth
    # })

    # KLfreqs <- c()
    # for (i in 1:nrow(KLinds)) {
    #     Kfreqs <- freqGrid[KLinds[i, 2]:KLinds[i, 3]]
    #     KLfreqs <- append(KLfreqs, Kfreqs)
    # }
    # return (KLfreqs);

    # Store (2):
    # freqdivideFreqGrid <- function(freqGrid, L, K) {
    # # set.seed(1)  # Set seed for reproducibility.

    # if ((K %% 2) == 0) {
    #     safeDist <- 1 + K/2  # 1 is added just to be more safe at the edges of the frequency grid. This is just a hackery.
    # }
    # else {
    #     safeDist <- 1 + (K-1)/2  # 1 is added just to be more safe at the edges of the frequency grid. This is just a hackery.
    # }
    # endIndex <- length(freqGrid) - safeDist
    # freqConsider <- freqGrid[safeDist:endIndex]
    # # LcentralFreqs <- divide(freqConsider, n=L, min_spacing=K+1)
    # LcentralFreqs <- sample(freqConsider, L, replace=FALSE, prob=rep(1/length(freqConsider), length(freqConsider)))  # replace=FALSE to prevent sampling the same frequency again. According to Suveges, each of the L central freqeuencies is selected with equal probability, so we pass an equal probability vector.
    # KLfreqs <- c()
    # for (i in 1:length(LcentralFreqs)) {
    #     index <- match(LcentralFreqs[i], freqGrid)
    #     if ((K %% 2) == 0) {
    #         k_ <- as.integer(K/2)
    #         lowerIndx <- index-k_
    #         upperIndx <- index+(K-k_-1)
    #         KLfreqs <- append(KLfreqs, freqGrid[lowerIndx:upperIndx])
    #     }
    #     else {
    #         kminusonehalf <- as.integer((K-1) / 2)
    #         lowerIndx <- index-kminusonehalf
    #         upperIndx <- index+kminusonehalf
    #         KLfreqs <- append(KLfreqs, freqGrid[lowerIndx:upperIndx])
    #     }
    # }
    # return (KLfreqs);
# }
