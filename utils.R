# ***************************************************************************
# Author: Yash Gondhalekar  Last updated: March, 2023

# Description: Contains utility functions used in other R scripts.

#              Notes:
#              Limitation of the function getLightCurve:
#              Both the duration (in hours) and the period (in days) must be
#              either an integer or a half-integer because the code requires
#              two times the duration and period to be an integer.
#              IMPORTANT: The arguments ar, ma, and order are obsolete. The
#              getLightCurve function, to which these were intended to passed
#              have these values harcoded inside that function. So passing
#              different values will not have any difference. Future versions
#              would try solving for this caveat.

# ***************************************************************************


library(forecast)
library(reticulate)
library('microbenchmark')
library('kernlab')
library(moments)

source_python("python_utils.py")

# This function simulates a light curve given some parameters of the light curve.
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
    # Ensure the user specifices legitimate values.
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

    # Add noise to the simulated transits.
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

    # Print information about the light curve created.
    print(sprintf('Length of time series = %d', length(y)))

    # Print some things
    print("==== Simulated light curve parameters ====")
    print(noquote(paste("Period (hours) = ", sprintf("%.3f", period * 24))))
    print(noquote(paste("Depth = ", sprintf("%.6f",  depth))))
    print(noquote(paste("Transit duration (hours) = ", sprintf("%.3f", inTransitTime))))
    print("==========================================")

    return (list(y, t, noiseStd, noiseIQR))
}

# Function to get ARIMA residuals ready for use with TCF.
getResidForTCF <- function(
    y  # Time series (must not be differenced because it is done internally).
) {
    max.p = 5
    max.q = 5
    max.d = 0
    ARIMA.fit = auto.arima(diff(y), stepwise=FALSE, approximation=FALSE, seasonal=FALSE, max.p=max.p, max.q=max.q, max.d=max.d) # leave d as 0. 
    print(ARIMA.fit)
    # return (ARIMA.fit)
    ARIMA.resid = residuals(ARIMA.fit)
    return (ARIMA.resid)
}

# Function to get ARMA residuals.
getARMAresid <- function(
    y
) {
    max.p = 5
    max.q = 5
    max.d = 0
    ARMA.fit = auto.arima(y, stepwise=FALSE, approximation=FALSE, seasonal=FALSE, max.p=max.p, max.q=max.q, max.d=max.d, d=0) #leave d as 0. 
    ARMA.resid = residuals(ARMA.fit)
    return (ARMA.resid)
}

# Function to get Gaussian Processes Regression residuals.
getGPRResid <- function(
    t, y
) {
    gp <- gausspr(t, y)
    predicted_y <- predict(gp, t)
    y <- y - predicted_y
    return (y)
}

# THIS IS A TRAIL FUNCTION, NOT USED ANYWHERE : It can be used to undifference a light curve.
# It is not guaranteed to work as expected.
undifferenceATimeSeries <- function(
    differenced,  # The differenced time series.
    orig_first_observation  # The first observation value in the original time series.
    # This original time series is the one that we want to predict. But this function assumes the user
    # themselves differenced the original light curve, in which case the user would have access to the first value.
) {
    undifferenced <- c()
    undifferenced <- c(undifferenced, orig_first_observation)
    undiff_value <- undifferenced[1]
    for (i in 1:length(differenced)) {
        undiff_value <- undiff_value + differenced[i]
        undifferenced <- c(undifferenced, undiff_value)
    }
    return (undifferenced)
}

# Function for folding a light curve at a given period.
plot_folded_lc <- function(lc, bestper)  # bestper is the period at which you want to fold the light curve.
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

# Function to create the frequency grid used for computing the BLS and TCF periodograms.
getFreqGridToTest <- function(
    t,  # Observation epochs (or the time values).
    period=NULL, duration=NULL,
    res=2,  # This the resolution in the time series that controls the candence. res=2 means cadence of 30 min and res=1 means 1hr, for example.
    ofac=1,  # Oversampling factor for frequency selection.
    useOptimalFreqSampling=FALSE,
    algo="BLS",  # Either BLS or TCF.
    lctype="sim"  # Either real or sim.
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

# Function to calculate the SNR of the periodogram peak. This automatically selects the highest peak and calculates SNR
# in a region around the peak. That region is defined by the argument oneSideWindowLength, which is set to 1500 by default.
# This means a window of 2*1500 = 3000 periods will be used as the region around the peak.
calculateSNR <- function(  # TODO: For making this more efficient, compute trend fit only for region around the periodogram peak - will save some time.
    periods,  # The periods at which the periodogram is computed.
    periodogramPower,  # This power must not be detrended since it is done internally.
    lambdaTrend=1,  # Parameter passed to cobs.
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

# Using the microbenchmark CRAN package, this allows timing the periodogram run - BLS or TCF.
timeAnalysis <- function(
    period, depth=0.01, duration=2, noiseType=1, ntransits=10,
    gaussStd=1e-4, ar=0.2, ma=0.2, res=2, order=c(1, 0, 1), algo="BLS",
    ofac=2, useOptimalFreqSampling=TRUE, times=1, lctype="sim", applyGPRforBLS=FALSE, applyARMAforBLS=FALSE
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

# Function to divide the frequency grid (freqGrid) into K*L frequencies used in the extreme value application.
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
