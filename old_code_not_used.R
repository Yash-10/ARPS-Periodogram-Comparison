###########################################################################
# Aim: Collection of functions written sometime ago, but not used anywhere.
#    - Kept only for demonstration.
###########################################################################

# Below cv function from https://rpubs.com/mengxu/loess_cv
library(bootstrap)
loess_wrapper_extrapolate <- function (x, y, span.vals = seq(0.5, 1, by = 0.05), folds = 5){
    # Do model selection using mean absolute error, which is more robust than squared error.
    mean.abs.error <- numeric(length(span.vals))

    # Quantify error for each span, using CV
    loess.model <- function(x, y, span){
    loess(y ~ x, span = span, control=loess.control(surface="direct"))
    }

    loess.predict <- function(fit, newdata) {
    predict(fit, newdata = newdata)
    }

    span.index <- 0
    for (each.span in span.vals) {
    span.index <- span.index + 1
    y.hat.cv <- crossval(x, y, theta.fit = loess.model, theta.predict = loess.predict, span = each.span, ngroup = folds)$cv.fit
    non.empty.indices <- !is.na(y.hat.cv)
    mean.abs.error[span.index] <- mean(abs(y[non.empty.indices] - y.hat.cv[non.empty.indices]))
    }

    # find the span which minimizes error
    best.span <- span.vals[which.min(mean.abs.error)]

    # fit and return the best model
    best.model <- loess(y ~ x, span = best.span, control=loess.control(surface="direct"))
    return(best.model)
}

standardizeAPeriodogram <- function(
    output,
    periodsToTry=NULL,  # This argument is only needed when algo="TCF" and not needed for algo="BLS".
    algo="BLS",
    mode='detrend'  # Other option is 'detrend' in which case only detrending is performed, no normalization using scatter is performed.
) {
    lambdaTrend <- 1
    lambdaScatter <- 1

    # (1) Remove trend.
    if (algo == "BLS") {
        lambdaTrend <- 1
        cobsxy50 <- cobs(output$periodsTested, output$spec, ic='BIC', tau=0.5, lambda=lambdaTrend, constraint="increase")  # If tau = 0.5 and lambda = 0 => Median regression fit.
    }
    else if (algo == "TCF") {
        lambdaTrend <- 1
        cobsxy50 <- cobs(periodsToTry, output$outpow, ic='BIC', tau=0.5, lambda=lambdaTrend, constraint="increase")
    }

    periodogramTrendRemoved <- cobsxy50$resid

    if (mode == 'detrend_normalize') {
        # (2) Remove local scatter in periodogram.
        scatterWindowLength <- 100
        Scatter <- computeScatter(periodogramTrendRemoved, windowLength=scatterWindowLength)
        if (algo == "BLS") {
            lambdaScatter <- 1
            cobsScatter <- cobs(output$periodsTested, Scatter, ic='BIC', tau=0.5, lambda=lambdaScatter)
        }
        else if (algo == "TCF") {
            lambdaScatter <- 1
            cobsScatter <- cobs(periodsToTry, Scatter, ic='BIC', tau=0.5, lambda=lambdaScatter)
        }

        normalizedPeriodogram <- periodogramTrendRemoved / cobsScatter$fitted
        return (normalizedPeriodogram);
    }
    else {  # mode == 'detrend'
        if (min(periodogramTrendRemoved) < 0) {
            return (periodogramTrendRemoved + abs(min(periodogramTrendRemoved)) + .Machine$double.eps)
        }
        else {
            return (periodogramTrendRemoved);
        }
    }
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


divide <- function(vec, n, min_spacing = 1) {  # Below approach to ensure the selected L frequencies are spaced by atleast the oversampling factor is taken from https://stackoverflow.com/a/66036847:
    unname(tapply(vec, ceiling(seq_along(vec) / min_spacing), sample, size = 1))
    # head(u[seq(1, length(u), by = 2)], n)
}

divide <- function(vec, n, min_spacing = 2) {
  idx <- seq_along(vec)
  repeat {
    k <- sort(sample(idx,n))
    if (all(diff(k)>=min_spacing)) break
  }
  vec[k]
}

### Store ###
# This is the earlier code (that did not work well), to divide the frequency grid into K*L frequencies.
# Divide the frequency into L bins, each with K datapoints.
## From https://stackoverflow.com/questions/57889573/how-to-randomly-divide-interval-into-non-overlapping-spaced-bins-of-equal-lengt
intervalLength <- length(freqGrid)
nBins <- L
binWidth <- K
binMinDistance <- 1
spaceToDistribute <- intervalLength - (nBins * binWidth + (nBins - 1) * binMinDistance)
distances <- diff(floor(c(0, sort(runif(nBins))) * spaceToDistribute))
startOfBin <- cumsum(distances) + (0:(nBins-1)) * 101
KLinds <- data.frame(bin = 1:nBins, startOfBin = startOfBin, endOfBin = startOfBin + binWidth - 1)

stopifnot(exprs={  # Check if the no. of frequencies in a bin is in fact equal to the desired number.
    length(freqGrid[KLinds[1, 2]:KLinds[1, 3]]) == binWidth
})

KLfreqs <- c()
for (i in 1:nrow(KLinds)) {
    Kfreqs <- freqGrid[KLinds[i, 2]:KLinds[i, 3]]
    KLfreqs <- append(KLfreqs, Kfreqs)
}
return (KLfreqs);

# Store (2):
freqdivideFreqGrid <- function(freqGrid, L, K) {
    # set.seed(1)  # Set seed for reproducibility.

    if ((K %% 2) == 0) {
        safeDist <- 1 + K/2  # 1 is added just to be more safe at the edges of the frequency grid. This is just a hackery.
    }
    else {
        safeDist <- 1 + (K-1)/2  # 1 is added just to be more safe at the edges of the frequency grid. This is just a hackery.
    }
    endIndex <- length(freqGrid) - safeDist
    freqConsider <- freqGrid[safeDist:endIndex]
    # LcentralFreqs <- divide(freqConsider, n=L, min_spacing=K+1)
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


validate1_evd <- function(  # Checks whether the values in the bootstrapped resample are actually from the original time series, which is a must.
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


findbestLandR <- function(  # Finds the optimal L and R values via grid search. It uses the AIC for finding the best {L, R} pair.
    Ls,
    Rs,
    period,
    depth,
    duration,
    ...
) {
    # *** CAUTION: Do not use this code with large Ls and Rs lengths. It is only meant to compare a few L and R pairs and not for large scale tuning ***

    stopifnot(exprs={
        length(Ls) == length(Rs)
    })
    minAIC <- Inf
    bestLR <- NULL
    for (i in 1:length(Ls)) {
        result <- evd(period, depth, duration, Ls[i], Rs[i], ...)
        aic <- result[2]
        if (aic < minAIC) {
            minAIC = aic
            bestLR <- c(Ls[i], Rs[i])
        }
    }
    return (bestLR)
}