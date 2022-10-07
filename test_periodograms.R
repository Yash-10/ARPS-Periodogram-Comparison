library(moments)

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
    checkConditions=TRUE  # Whether to perform small unit tests within the code to ensure the output is as expected. Recommended: Set to TRUE for almost all cases, but added this option so that `uniroot` does not throw error when finding the limiting depth - see below functions.
) {
    if (checkConditions) {
        stopifnot(exprs = {
            period > 0
            depth >= 0
            depth <= 100
            duration > 0
            duration <= 1
            ntransits >= 0
        })   
    }

    # Create a simulated planet transit time series based on period, depth, and transit duration.
    inTransitValue = 1 - (depth / 100) * 1
    # TODO: do not do round() below since then 1/48 and 1/36 for a 3-day period is the same, which MUST NOT be. UPDATE: STILL does not work.
    inTransitTime = duration * period * 24  # inTransitTime is the actual absolute in-transit time (in hours).
    constTime = (period * 24 - 2 - inTransitTime)

    if (checkConditions) {
        stopifnot(exprs = {
            inTransitValue < 1
            inTransitValue >= 0
            inTransitTime > 0
            constTime > 0
        })
    }

    # First generate the time epochs.
    tIncrement <- 1 / res
    tEnd <- constTime + ntransits * (2 + inTransitTime + constTime)
    t <- seq(from = 0, by = tIncrement, to = tEnd)

    # The deemed constant value will be 1 and the transits would be scaled with respect to 1. For eg: 1 --> 0.998 --> 1 --> 0.998 ...
    y <- rep(1.0, each = as.integer(constTime*res)) # Start with some constant level.

    for (n in 1:ntransits) {
        y <- append(y, seq(from = 1, to = inTransitValue, length.out = res))  # Decrement from 1 to depth.
        endt <- as.integer(inTransitTime*res) - 1
        for (j in 1:endt) {
            y <- append(y, inTransitValue)
        }
        y <- append(y, seq(from = inTransitValue, to = 1, length.out = res))  # Increment from depth to 1.
        constt <- as.integer(constTime*res) + 1
        for (j in 1:constt) {
            y <- append(y, 1)
        }
    }
    y <- append(y, 1)

    if (noiseType == 1) {
        set.seed(1)
        noise <- rnorm(length(y), mean = 0, sd = gaussStd)
        y <- y + noise  # 0.01% Gaussian noise.
        noiseStd <- sd(noise)
        noiseIQR <- IQR(noise)
    }
    else if (noiseType == 2) {
        set.seed(1)
        # Note that autoregresive noise has been scaled by `gaussStd` so that the range of values in both Gaussian and autoregressive case look similar.
        # We can as well remove that scaling, but using it allows us to compare Gaussian and autoregressive cases much easier since then the only difference is that Gaussian is uncorrelated and autoregressive is correlated noise. And we don't need to worry about autoregressive and Gaussian noises being on different scales.
        # autoRegNoise <- arima.sim(list(order = c(3,0,3), ar = c(0.2,0.3,0.2), ma=c(0.2,0.2,0.3)), n = length(y)) * gaussStd  # It has only AR and MA components, no differencing component. So it is ARMA and not ARIMA. Note: Keep ar and ma < 0.5.
        # Note: Use the above example (immediate above) or modifications to it to use more stronger autoregressive noises. 
        autoRegNoise <- arima.sim(list(order=c(1, 0, 1), ar=0.2, ma=0.2), n = length(y)) * gaussStd
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

    # Print some things
    print("==== Simulated light curve parameters ====")
    print(noquote(paste("Period (hours) = ", sprintf("%.3f", period * 24))))
    print(noquote(paste("Depth = ", sprintf("%.6f",  depth))))
    print(noquote(paste("Transit duration (hours) = ", sprintf("%.3f", inTransitTime))))

    return (list(y, t, noiseStd, noiseIQR))

    # # Add noise to time series, if any.
    # if (noiseType == 1) {
    #     set.seed(1)
    #     y <- y + 0.3 * rnorm(length(y))  # TODO: Remove scaling and incorporate it in stddev?
    # }
    # else if (noiseType == 2) {
    #     set.seed(1)
    #     autoRegNoise <- arima.sim(model = list(order=c(1, 0, 1), ar=0.2, ma=0.2), n = length(y))  # It has only AR and MA components, no differencing component. So it is ARMA and not ARIMA. Note: Keep ar and ma < 0.5.
    #     y <- y + autoRegNoise
    # }
}

getStandardPeriodogram <- function(
    period,  # What period (in days) do you want to have in your light curve, will be a single value. eg: 1/3/5/7/9.
    depth,  # What depth (in % of the star's presumed constant level which is 1) do you want to have in your light curve, will be a vector. eg: 0.01/0.05/0.1/0.15/0.2.
    duration,  # What transit duration (as a fraction of the period) do you want to have in your light curve. eg: 1/24.
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
