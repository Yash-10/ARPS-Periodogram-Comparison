# ***************************************************************************
# Author: Yash Gondhalekar  Last updated: March, 2023

# Description: This script contains code for standardizing a periodogram.
#              The function standardPeriodogram can be used for getting the
#              standardized periodogram and (optionally) plotting the results.
#              In this code, tcf periods are scaled by `res` before passing
#              into tcf() since tcf works with cadences rather than absolute
#              time values, unlike BLS. The period tcf prints on terminal
#              would differ from the period at max(power) shown in the plot.
#              Such a scaling is also done in eva_periodogram.R since TCF
#              calculates periodogram in units of cadences rather than
#              absolute time values.
#              IMPORTANT: The arguments ar, ma, and order are obsolete. The
#              getLightCurve function, to which these were intended to passed
#              have these values harcoded inside that function. So passing
#              different values will not have any difference. Future versions
#              would try solving for this caveat.

# Example:
# An example run is: standardPeriodogram(1, 0.01, 1) - which simulates a
# period = 1 day, depth = 0.01%, and duration = 1 hr planet. By default BLS is
# run but can be changed using the `algo` argument.

# ***************************************************************************

library(moments)
library('cobs')
source('TCF3.0/intf_libtcf.R')
source('BLS/bls.R')
source('utils.R')
source('eva_periodogram.R')

computeScatter <- function(
    cobsTrendResid,     # The data for which the scatter needs to be computed. This is generally the detrended periodogram.
    windowLength=1000,  # Window length on one side of the focal point. The actual window length used is 2*windowLength.
    algo="BLS"          # Periodogram algorithm to use. BLS or TCF.
){
    scatterVals = c()
    for (i in 1:length(cobsTrendResid)) {
        if (i <= windowLength) {
            u = 2 * windowLength
            l = i + 1
            scatterVal = IQR(c(cobsTrendResid[1:i], cobsTrendResid[l:u]))
            stopifnot(exprs={
                length(c(cobsTrendResid[1:i], cobsTrendResid[l:u])) == 2 * windowLength
            })
            # print(length(c(cobsTrendResid[1:i], cobsTrendResid[l:u])))
            scatterVals <- append(scatterVals, scatterVal)
        }
        else if (i + windowLength >= length(cobsTrendResid)) {
            ll = i - 2 * windowLength - i + length(cobsTrendResid)
            ul = i - 1
            ls = i
            us = length(cobsTrendResid)-1
            # stopifnot(exprs={
            #     length(c(cobsTrendResid[ll:ul], cobsTrendResid[ls:us])) == 2 * windowLength
            # })
            # print(length(c(cobsTrendResid[ll:ul], cobsTrendResid[ls:us])))
            scatterVal = IQR(c(cobsTrendResid[ll:ul], cobsTrendResid[ls:us]))
            scatterVals <- append(scatterVals, scatterVal)
        }
        else {
            l_ = i-windowLength
            u_ = i+windowLength-1
            stopifnot(exprs={
                length(cobsTrendResid[l_:u_]) == 2 * windowLength
            })
            scatterVal <- IQR(cobsTrendResid[l_:u_])
            scatterVals <- append(scatterVals, scatterVal)
        }
    }
    stopifnot(exprs={
        length(scatterVals) == length(cobsTrendResid)
    })
    return (scatterVals)
}

standardPeriodogram <- function(
    period=NULL,  # in days.
    depth=NULL,  # in %
    duration=NULL,  # in hours.
    y=NULL, t=NULL,
    noiseType = 0,  # 0 for no noise, 1 for white Gaussian noise, and 2 for autoregressive noise (if 2, the ARMA model is internally fixed; also no differencing is used, so it assumes the time-series is stationary or already differenced.)
    algo = "BLS",  # or "TCF"
    windowLength=100,  # Window length for SNR estimation.
    plot = TRUE,  # Whether to plot result.
    ntransits=10,  # no. of transits
    ofac=2,  # oversampling factor to compute the periodograms.
    res=2,  # Light curve resolution, see getLightCurve().
    showFAP = FALSE,  # Whether to show the calculated false alarm probability in the plot. If TRUE, it will take much more time since internally the evd() function is run.
    # Below four arguments set the noise parameters. See getLightCurve in utils.R
    gaussStd=1e-4,
    ar=0.2,
    ma=0.2,
    order=c(1, 0, 1),
    L=500, R=500,  # Parameters used for extreme value application.
    useOptimalFreqSampling=FALSE,  # Whether to use Ofir's optimal frequency sampling.
    seedValue=1,  # Seed value to use. Can be used for reproducibility.
    lctype="sim",  # Either sim or real.
    applyARMAforBLS=FALSE  # Whether to apply ARMA before BLS.
){
    # Perform some checks.
    if (lctype == "sim" && (is.null(period) | is.null(depth) | is.null(duration))) {
        stop("type is set to `sim`, but at least one of {period, depth, or duration} is not specified!")
    }
    if (lctype == "real" && (is.null(y) | is.null(t))) {
        stop("type is set to `real`, but at least one of {y, t} is not specified!")
    }
    if (lctype == "real") {
        period <- depth <- duration <- noiseType <- ntransits <- ar <- ma <- order <- gaussStd <- NULL
        significanceMode <- 'max'  # Since for real light curves, passing `expected_peak` is not possible.
        res <- 2
    }

    if (lctype == "sim") {
        # Generate light curve using the parameters.
        yt <- getLightCurve(period, depth, duration, noiseType=noiseType, ntransits=ntransits, res=res, gaussStd=gaussStd, ar=ar, ma=ma, order=order, seedValue=seedValue)
        y <- unlist(yt[1])
        t <- unlist(yt[2])
        noiseStd <- unlist(yt[3])
        noiseIQR <- unlist(yt[4])
    }

    if (lctype == "sim") {
        # Special case (TCF fails if absolutely no noise -- so add a very small amount of noise just to prevent any errors).
        if (noiseType == 0 && algo == "TCF") {
            y <- y + 10^-10 * rnorm(length(y))
        }
    }

    # Create frequency grid.
    freqGrid <- getFreqGridToTest(t, period, duration, res=res, ofac=ofac, useOptimalFreqSampling=useOptimalFreqSampling, algo=algo, lctype=lctype)

    print(min(freqGrid))
    print(max(freqGrid))

    stopifnot(exprs={
        all(freqGrid <= res / 2)
    })

    if (algo == "BLS") {
        if (isTRUE(noiseType) == 2) {
            # Run Gaussian Process Regression and get the fit residuals.
            y <- getGPRResid(t, y)
        }
        if (applyARMAforBLS) {
            # Run ARMA fit and get the residuals.
            y <- getARMAresid(y)
        }
        # Run BLS.
        output <- bls(y, t, bls.plot = FALSE, per.min=min(1/freqGrid), per.max=max(1/freqGrid), nper=length(freqGrid))
    }
    else if (algo == "TCF") {
        # Setup for setting the frequencies (or periods) to use to compute the TCF periodogram.
        # This setup ensures that the frequencies used for BLS and TCF are the same. This has been tested.
        fstep <- (max(freqGrid) - min(freqGrid)) / length(freqGrid)
        freqs <- seq(from = min(freqGrid), by = fstep, length.out = length(freqGrid))
        periodsToTry <- 1 / freqs
        # Get ARIMA residual.
        residTCF <- getResidForTCF(y)
        # Run TCF on ARIMA residual.
        output <- tcf(residTCF, p.try = periodsToTry*res, print.output = TRUE)
    }

    # (1) Remove trend in periodogram
    # We use constraint='increase'. But that condition may be removed.
    if (algo == "BLS") {
        lambdaTrend <- 1
        cobsxy50 <- cobs(output$periodsTested, output$spec, ic='BIC', tau=0.5, lambda=lambdaTrend, constraint="increase")  # If tau = 0.5 and lambda = 0 => Median regression fit.
        cobsxy501 <- cobs(output$periodsTested, output$spec, ic='BIC', tau=0.9, lambda=lambdaTrend)
        cobsxy502 <- cobs(output$periodsTested, output$spec, ic='BIC', tau=0.99, lambda=lambdaTrend)
    }
    else if (algo == "TCF") {
        lambdaTrend <- 1
        cobsxy50 <- cobs(periodsToTry, output$outpow, ic='BIC', tau=0.5, lambda=lambdaTrend, constraint="increase")
        cobsxy501 <- cobs(periodsToTry, output$outpow, ic='BIC', tau=0.9, lambda=lambdaTrend)
        cobsxy502 <- cobs(periodsToTry, output$outpow, ic='BIC', tau=0.99, lambda=lambdaTrend)
    }

    periodogramTrendRemoved <- cobsxy50$resid

    # (2) Remove local scatter in periodogram
    scatterWindowLength <- 100
    Scatter <- computeScatter(periodogramTrendRemoved, windowLength=scatterWindowLength)
    # print("Scatter")
    if (algo == "BLS") {
        lambdaScatter <- 1
        cobsScatter <- cobs(output$periodsTested, Scatter, ic='BIC', tau=0.5, lambda=lambdaScatter)
        # cobss50 <- cobs(output$periodsTested, periodogramTrendRemoved, ic='BIC', tau=0.5, lambda=lambdaTrend, constraint="increase")  # If tau = 0.5 and lambda = 0 => Median regression fit.
        # cobss501 <- cobs(output$periodsTested, periodogramTrendRemoved, ic='BIC', tau=0.9, lambda=lambdaTrend, constraint="increase")
        # cobss502 <- cobs(output$periodsTested, periodogramTrendRemoved, ic='BIC', tau=0.99, lambda=lambdaTrend)
    }
    else if (algo == "TCF") {
        lambdaScatter <- 1
        cobsScatter <- cobs(periodsToTry, Scatter, ic='BIC', tau=0.5, lambda=lambdaScatter)
    }

    # Get the normalized periodogram.
    normalizedPeriodogram <- periodogramTrendRemoved / cobsScatter$fitted

    if (algo == "BLS") {
        returnVals <- list(normalizedPeriodogram, output$periodsTested)
    }
    else {
        returnVals <- list(normalizedPeriodogram, periodsToTry)
    }

    if (showFAP) {
        # Call extreme value analysis code.
        result <- evd(period, depth, duration, noiseType=noiseType, algo=algo, ofac=ofac, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, ar=ar, ma=ma, order=order, seedValue=seedValue)
        print(sprintf("FAP = %.10f", result[1]))
        fap <- result[1]
    }

    # Plotting things.
    if (plot) {
        dev.new(width=20, height=10)

        par("mar" = c(5, 6, 4, 2))

        cexVal = 1.7
        mat1 <- matrix(c(
            1, 3, 4,
            2, 3, 4,
            5, 6, 7,
            5, 6, 7), nrow = 4, ncol = 3, byrow = TRUE
        )
        layout(mat = mat1,
            heights = c(1),    # Heights of the two rows
            widths = c(1)
        )     # Widths of the two columns

        if (algo == "BLS") {
            pergram <- output$spec
        }
        else if (algo == "TCF") {
            pergram <- output$outpow
        }

        plot(t, y, type='l', main=sprintf("Period: %.1f days, depth: %.3f (pct), duration: %.1f (hrs)", period, depth, duration), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xlab='time (hrs)', ylab='normalized flux')
        acfEstimate <- acf(y, plot = FALSE, na.action = na.pass)
        lJStats <- Box.test(y, lag = 1, type = "Ljung")  # We want to see autocorrelation with each lag, hence pass lag = 1.
        plot(acfEstimate, main=sprintf("P(Ljung-Box) = %s, lag-1 acf = %s", lJStats[3], acfEstimate$acf[[2]]), cex=2)

        plot(cobsxy50$x, pergram, type = 'l', main=sprintf("Original %s periodogram", algo), log='x', xlab='Period (hrs) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
        lines(cobsxy50$x, cobsxy50$fitted, type = 'l', col='red', lwd=3.0)
        # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
        # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
        rug(cobsxy50$knots)
        legend(x = "topleft", lty = 1, text.font = 6, 
            col= c("red"), text.col = "black", 
            legend=c("trend fit")
        )

        plot(cobsxy50$x, normalizedPeriodogram, type = 'l', main=sprintf('Standardized %s periodogram', algo), log='x', xlab='Period (hrs) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)

        # plot.new()  # Just show an empty plot.
        lc <- data.frame(t=t, y=y)
        if (lctype == "sim") {  # TODO: In future, we can still allow plotting the folded lc but with the predicted period.
            plot_folded_lc(lc, bestper=period*24*res)
        }
        else {
            if (algo == "BLS") {
                plot_folded_lc(lc, bestper=output$per*24*res)
            }
            else if (algo == "TCF") {
                plot_folded_lc(lc, bestper=output$inper[which.max(output$outpow)])
            }
        }

        # Plot histogram of periodogram. Shows log-frequency on y-axis in histogram for better visualization.
        hist.data = hist(pergram, breaks=50, plot = FALSE)
        # Compute skewness and kurtosis of the original and standardized histograms.
        ### Refer https://brownmath.com/stat/shape.htm for more information ###
        ### Note: R does NOT compute the "excess kurtosis".
        # The kurtosis is calculated as follows:
        # ```
        # n <- length(x)
        # n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
        # ``` Taken from https://stackoverflow.com/a/21484052
        SkewnessBefore <- skewness(pergram)
        KurtosisBefore <- kurtosis(pergram)

        if (showFAP) {
            plot(hist.data$count, type='h', log='y', main=sprintf('Original %s periodogram histogram, FAP: %s', algo, formatC(fap, format = "e", digits = 5)), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
        }
        else {
            plot(hist.data$count, type='h', log='y', main=sprintf('Original %s periodogram histogram', algo, SkewnessBefore, KurtosisBefore), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
        }
        axis(1, at=1:length(hist.data$mids), labels=hist.data$mids)

        SkewnessAfter <- skewness(normalizedPeriodogram)
        KurtosisAfter <- kurtosis(normalizedPeriodogram)

        hist.data = hist(normalizedPeriodogram, breaks=50, plot = FALSE)
        plot(hist.data$count, type='h', log='y', main=sprintf("Standardized %s periodogram histogram", algo), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
        axis(1, at=1:length(hist.data$mids), labels=hist.data$mids)
    }
}
