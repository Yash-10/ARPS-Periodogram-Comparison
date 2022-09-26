## Utility functions.
library(forecast)

getResidForTCF <- function(
    y  # Time series (must not be differenced because it is done internally).
) {
    max.p = 5
    max.q = 5
    max.d = 0
    ARIMA.fit = auto.arima(diff(y), stepwise=FALSE, approximation=FALSE, seasonal=FALSE, max.p=max.p, max.q=max.q, max.d=max.d) #leave d as 0. 
    # Simple statistics of ARIMA residuals 
    ARIMA.resid = residuals(ARIMA.fit)
    return (ARIMA.resid);
}

getFreqGridToTest <- function(
    t,  # Observation epochs.
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
        # Note: When we oversample, we are essentially trying to impose no constraints on the frequencies to be tested - that is helpful in general: see https://arxiv.org/pdf/1712.00734.pdf
        # Note that too much oversampling can lead to artifacts. These artifacts can be wrongly interpreted as a true periodic component in the periodogram.
        freqStep <- (freqMax - freqMin) / (nfreq * ofac)  # Oversampled by a factor, `ofac`.
    }

    freqGrid <- seq(from = freqMin, to = freqMax, by = freqStep)  # Goes from ~0.001 to max freq set by the time spacing (NOTE: fmax must be <= Nyquist frequency = 1/(2*delta_t) -- from Suveges, 2014), where delta_t here is res.
    freqGrid <- freqGrid[-length(freqGrid)]  # Remove the last frequency. This is done to prevent periodogram returning nan power at frequency = Nyquist frequency.
    print(sprintf("No. of frequencies in grid: %f", length(freqGrid)))

    return (freqGrid);
}

# TODO: The below freq grid division does not ensure non-overlap, which is needed as per suveges....

freqdivideFreqGrid <- function(freqGrid, L, K) {
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

    # TODO: Need to verify this works as expected for large oversampling factors. I think that the hackery below might still yield errors for larger ofac values than 2.
    if ((K %% 2) == 0) {
        safeDist <- 1 + K/2  # 1 is added just to be more safe at the edges of the frequency grid. This is just a hackery.
    }
    else {
        safeDist <- 1 + (K-1)/2  # 1 is added just to be more safe at the edges of the frequency grid. This is just a hackery.
    }
    endIndex <- length(freqGrid) - safeDist
    freqConsider <- freqGrid[safeDist:endIndex]
    LcentralFreqs <- sample(freqConsider, L, replace=FALSE, prob=rep(1/length(freqConsider), length(freqConsider)))  # replace=FALSE to prevent sampling the same frequency again. According to Suveges, each of the L central freqeuencies is selected with equal probability, so we pass an equal probability vector.
    KLfreqs <- c()
    for (i in 1:length(LcentralFreqs)) {
        index <- match(LcentralFreqs[i], freqGrid)
        if ((K %% 2) == 0) {
            k_ <- as.integer(K/2)
            lowerIndx <- index-k_
            upperIndx <- index+(K-k_-1)
            KLfreqs <- append(KLfreqs, freqGrid[lowerIndx:upperIndx])
        }
        else {
            kminusonehalf <- as.integer((K-1) / 2)
            lowerIndx <- index-kminusonehalf
            upperIndx <- index+kminusonehalf
            KLfreqs <- append(KLfreqs, freqGrid[lowerIndx:upperIndx])
        }
    }
    return (KLfreqs);
}
