## Utility functions.
library(forecast)
library(reticulate)
library('microbenchmark')
library('kernlab')
library(moments)

source_python("python_utils.py")


###### Note: Multiply the depth value by 100 to get the depth in % ######

getLightCurve <- function(
    period,  # What period (in days) do you want to have in your light curve, will be a single value. eg: 1/3/5/7/9.
    depth,  # What depth (in % of the star's presumed constant level which is 1) do you want to have in your light curve, will be a single value. eg: 0.01/0.05/0.1/0.15/0.2.
    duration,  # What transit duration (as a fraction of the period) do you want to have in your light curve. eg: 1/24.
    res=2,  # Resolution to create the light curve. If res = 1, then the light curve cannot handle things that happen on scales less than an hour. Hence, a 1hr transit duration and 0.75hr transit duration would be treated same. If you want to avoid this, set a higher resolution using this parameter.
    # Note: The definition of transit duration used in the code is how many points there are at the in-transit level whereas in astronomy it is, how many points are there before you reach the constant value again taking into account points in going from 1 --> inTransitValue.
    noiseType=0,  # 1 --> Gaussian noise, 2 --> Autoregressive noise. If autoregressive noise, (1, 0, 1) model is used. To change it, need to change the source code.
    ntransits=10,  # No. of transits in the whole time series. Note: It must be >=3, otherwise BLS/TCF matching filter periodograms might not work.
    ### VIMP note: While giving inputs, never give a period like 1 since duration would be a fraction of 1 which is a float number and since the code rounds the result, results might not be correct.
    ### To prevent such issues, if your `duration` is, say, 1/24 times the period (eg: 1 hr duration for a 1 day period), then pass period = 24 and duration = 1/24 instead of period = 1 and duration = 1/24.
    gaussStd=1e-4,
    ar=0.2,
    ma=0.2,
    order=c(1, 0, 1),  # ARMA, no differencing.
    checkConditions=TRUE,  # Whether to perform small unit tests within the code to ensure the output is as expected. Recommended: Set to TRUE for almost all cases, but added this option so that `uniroot` does not throw error when finding the limiting depth - see below functions.
    seedValue=1
) {
    # *** IMPORTANT LIMITATION OF THIS FUNCTION ***
    # -> Both the duration (in hours) and the period (in days) must be either integer or half-integer because the code requires two times the duration and period to be an integer.
    # *********************************************
    if (checkConditions) {
        stopifnot(exprs = {
            period > 0
            depth >= 0
            depth <= 100
            duration >= 0
            ntransits >= 0
        })
    }

    # Create a simulated planet transit time series based on period, depth, and transit duration.
    inTransitValue = 1 - (depth / 100) * 1
    inTransitTime = duration  # inTransitTime is the actual absolute in-transit time (in hours).
    constTime = (period * 24 - 2 / res - inTransitTime)

    if (checkConditions) {
        stopifnot(exprs = {
            inTransitValue < 1
            inTransitValue >= 0
            inTransitTime > 0  # NOTE: 0 transit duration is not allowed by the code as of now.
            constTime > 0
        })
    }

    # The deemed constant value will be 1 and the transits would be scaled with respect to 1. For eg: 1 --> 0.998 --> 1 --> 0.998 ...
    # Note one caveat here is that if constTime * res is not an integer, this light curve simulation might fail.
    y <- rep(1.0, each = as.integer(constTime*res))  # Start with some constant level.

    for (n in 1:ntransits) {
        y <- append(y, seq(from = 1, to = inTransitValue, length.out = res))  # Decrement from 1 to depth.
        endt <- as.integer(inTransitTime*res) - 1
        if (endt != 0) {  # Special case when 0 < duration < 1 (assuming res=2).
            for (j in 1:endt) {
                y <- append(y, inTransitValue)
            }
        }
        y <- append(y, seq(from = inTransitValue, to = 1, length.out = res))  # Increment from depth to 1.
        constt <- as.integer(constTime*res) + 1
        for (j in 1:constt) {
            y <- append(y, 1)
        }
    }
    y <- append(y, 1)

    # Generate the time epochs.
    tIncrement <- 1 / res
    t <- seq(from = 0, by = tIncrement, length.out = length(y))

    if (noiseType == 1) {
        set.seed(seedValue)
        noise <- rnorm(length(y), mean = 0, sd = gaussStd)
        y <- y + noise  # 0.01% Gaussian noise.
        noiseStd <- sd(noise)
        noiseIQR <- IQR(noise)
    }
    else if (noiseType == 2) {
        set.seed(seedValue)
        # Note that autoregresive noise has been scaled by `gaussStd` so that the range of values in both Gaussian and autoregressive case look similar.
        # We can as well remove that scaling, but using it allows us to compare Gaussian and autoregressive cases much easier since then the only difference is that Gaussian is uncorrelated and autoregressive is correlated noise. And we don't need to worry about autoregressive and Gaussian noises being on different scales.
        autoRegNoise <- arima.sim(list(order = c(3,0,3), ar = c(0.2,0.3,0.2), ma=c(0.2,0.2,0.3)), n = length(y)) * gaussStd  # It has only AR and MA components, no differencing component. So it is ARMA and not ARIMA. Note: Keep ar and ma < 0.5.
        # Note: Use the above example (immediate above) or modifications to it to use more stronger autoregressive noises. 
        # autoRegNoise <- arima.sim(list(order=c(1, 0, 1), ar=0.2, ma=0.2), n = length(y)) * gaussStd
        y <- y + autoRegNoise
        noiseStd <- sd(autoRegNoise)
        noiseIQR <- IQR(autoRegNoise)
    }
    else if (noiseType == 0) {
        noiseStd <- 0
        noiseIQR <- 0
    }

    if (checkConditions) {
        stopifnot(exprs={
            length(y) == length(t)
        })
    }

    print(sprintf('Length of time series = %d', length(y)))

    # Print some things
    print("==== Simulated light curve parameters ====")
    print(noquote(paste("Period (hours) = ", sprintf("%.3f", period * 24))))
    print(noquote(paste("Depth = ", sprintf("%.6f",  depth))))
    print(noquote(paste("Transit duration (hours) = ", sprintf("%.3f", inTransitTime))))
    print("==========================================")

    return (list(y, t, noiseStd, noiseIQR))
}

getStandardPeriodogram <- function(
    period,  # What period (in days) do you want to have in your light curve, will be a single value. eg: 1/3/5/7/9.
    depth,  # What depth (in % of the star's presumed constant level which is 1) do you want to have in your light curve, will be a vector. eg: 0.01/0.05/0.1/0.15/0.2.
    duration,  # What transit duration (in hours) do you want to have in your light curve. eg: 1/24.
    # Note: The definition of transit duration used in the code is how many points there are at the in-transit level whereas in astronomy it is, how many points are there before you reach the constant value again taking into account points in going from 1 --> inTransitValue.
    noiseType=0,  # 1 --> Gaussian noise, 2 --> Autoregressive noise. If autoregressive noise, (1, 0, 1) model is used. To change it, need to change the source code.
    ntransits=10,  # No. of transits in the whole time series. Note: It must be >=3, otherwise BLS/TCF matching filter periodograms might not work.
    ### VIMP note: While giving inputs, never give a period like 1 since duration would be a fraction of 1 which is a float number and since the code rounds the result, results might not be correct.
    ### To prevent such issues, if your `duration` is, say, 1/24 times the period (eg: 1 hr duration for a 1 day period), then pass period = 24 and duration = 1/24 instead of period = 1 and duration = 1/24.
    algo="BLS"  # or "TCF"
){

    yt <- getLightCurve(period, depth, duration, noiseType=noiseType, ntransits=ntransits)

    output <- standardPeriodogram(unlist(yt[1]), unlist(yt[2]), period, depth, duration, algo=algo, noiseType = 0, plot=FALSE)  # 0 noise type we add custom noise in this function, so we should not add it again.
    return (output);
}

getResidForTCF <- function(
    y  # Time series (must not be differenced because it is done internally).
) {
    max.p = 5
    max.q = 5
    max.d = 0
    ARIMA.fit = auto.arima(diff(y), stepwise=FALSE, approximation=FALSE, seasonal=FALSE, max.p=max.p, max.q=max.q, max.d=max.d) #leave d as 0. 
    print(ARIMA.fit)
    # Simple statistics of ARIMA residuals
    # plot(diff(y), col='black', type='l')
    # lines(fitted(ARIMA.fit), col='red')
    # return (1)
    ARIMA.resid = residuals(ARIMA.fit)
    return (ARIMA.resid)
}

getARMAresid <- function(
    y
) {
    max.p = 5
    max.q = 5
    max.d = 0
    ARMA.fit = auto.arima(y, stepwise=FALSE, approximation=FALSE, seasonal=FALSE, max.p=max.p, max.q=max.q, max.d=max.d, d=0) #leave d as 0. 
    ARMA.resid = residuals(ARMA.fit)
    # plot(y, col='black', type='l')
    # lines(fitted(ARMA.fit), col='red')
    # plot(ARMA.resid)
    return (ARMA.resid)
}

getGPRResid <- function(
    t, y
) {
    gp <- gausspr(t, y)
    predicted_y <- predict(gp, t)
    y <- y - predicted_y
    return (y)
}

plot_folded_lc <- function(lc, bestper)
{  
    Phase = (lc[,1] - min(lc[,1], na.rm=TRUE)) %% bestper / bestper
    Range_lc = 1.3*diff(range(lc[,2], na.rm=TRUE)) + min(lc[,2], na.rm=TRUE)
    Text0.pos = 1.35*diff(range(lc[,2], na.rm=TRUE)) + min(lc[,2], na.rm=TRUE)
    plot(
        Phase, lc[,2], xlab = "Phase", ylab ="Normalized Flux", cex.axis=1.0, cex.lab=1.2, cex = 0.5,
        pch = 20, col="#00000080", ylim=c(min(lc[,2], na.rm=TRUE), Range_lc), xaxs='i'
    )
    axis(1,at=seq(from=0, to=1, by=.1), labels = F, tcl=0.25)
    axis(1,at=seq(from=0, to=1, by=.2), labels = F, tcl=0.5)
    text(0.5, Text0.pos, "Phase-folded light curve", pos=1, cex=1.2)
}

getFreqGridToTest <- function(
    t,  # Observation epochs.
    period=NULL, duration=NULL,
    res=2,  # This the resolution in the time series that controls the candence. res=2 means cadence of 30 min and res=1 means 1hr, for example.
    ofac=1,  # Oversampling factor for frequency selection.
    useOptimalFreqSampling=FALSE,
    algo="BLS",
    lctype="sim"
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

    if (useOptimalFreqSampling & (lctype == "sim")) {  # TODO: This is done because as of now I don't know how to use optimal freq samplin when the true period/frequency is NOT known.
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
    print(sprintf("No. of frequencies in grid [algo = %s]: %f", algo, length(freqGrid)))

    return (freqGrid);
}

calculateSNR <- function(  # TODO: For making this more efficient, compute trend fit only for region around the periodogram peak - will save some time.
    periods,  # The periods at which the periodogram is computed.
    periodogramPower,  # This power must not be detrended since it is done internally.
    lambdaTrend=1,
    oneSideWindowLength=1500  # NOTE: This is one side, the actual window length will be 2*oneSideWindowLength.
) {
    cobsTrend <- cobs(log10(periods), periodogramPower, ic='BIC', tau=0.5, lambda=lambdaTrend)
    detrended <- cobsTrend$resid

    # Uncomment the below five lines to use standardization while SNR calculation.
    # scatterWindowLength <- 100
    # Scatter <- computeScatter(detrended, windowLength=scatterWindowLength)
    # lambdaScatter <- 1
    # cobsScatter <- cobs(periods, Scatter, ic='BIC', tau=0.5, lambda=lambdaScatter)
    # detrended <- detrended / cobsScatter$fitted

    # return (list(periods, detrended, cobsTrend$fitted, periodogramPower))
    lowerInd <- which.max(detrended) - oneSideWindowLength
    upperInd <- which.max(detrended) + oneSideWindowLength
    if (which.max(detrended) - oneSideWindowLength < 1) {
        lowerInd <- 1
    }
    if (which.max(detrended) + oneSideWindowLength > length(detrended)) {
        upperInd <- length(detrended)
    }
    consider <- detrended[lowerInd:upperInd]
    snr <- max(consider) / mad(consider)
    return (snr)
}

timeAnalysis <- function(
    period, depth=0.01, duration=2, noiseType=1, ntransits=10,
    gaussStd=1e-4, ar=0.2, ma=0.2, res=2, order=c(1, 0, 1), algo="BLS",
    ofac=2, useOptimalFreqSampling=TRUE, times=1, lctype="sim"
) {
    yt <- getLightCurve(period, depth, duration, noiseType=noiseType, ntransits=ntransits, gaussStd=gaussStd, ar=ar, ma=ma, order=order, res=res, checkConditions=TRUE, seedValue=42)
    y <- unlist(yt[1])
    t <- unlist(yt[2])

    freqGrid <- getFreqGridToTest(t, period, duration, res=res, ofac=ofac, useOptimalFreqSampling=useOptimalFreqSampling, algo=algo, lctype=lctype)

    stopifnot(exprs={
        all(freqGrid <= res / 2)  # No frequency must be greater than the Nyquist frequency.
    })

    print(sprintf("Max frequency: %f, Min frequency: %f", max(freqGrid), min(freqGrid)))

    if (algo == "BLS") {
        print("BLS periodogram time benchmark...")
        if (isTRUE(noiseType == 2) | applyGPRforBLS) {
            y <- getGPRResid(t, y)  # Run Gaussian Processes Regression on light curve if autoregressive noise is present.
        }
        if (applyARMAforBLS) {
            y <- getARMAresid(y)
        }
        microbenchmark(
            bls(y, t, bls.plot = FALSE, per.min=min(1/freqGrid), per.max=max(1/freqGrid), nper=length(freqGrid)),
            times=times
        )
    }
    else if (algo == "TCF") {
        fstep <- (max(freqGrid) - min(freqGrid)) / length(freqGrid)
        freqs <- seq(from = min(freqGrid), by = fstep, length.out = length(freqGrid))
        periodsToTry <- 1 / freqs
        residTCF <- getResidForTCF(y)
        print("TCF periodogram time benchmark...")
        microbenchmark(
            tcf(residTCF, p.try = periodsToTry * res, print.output = FALSE),
            times=times
        )
    }
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