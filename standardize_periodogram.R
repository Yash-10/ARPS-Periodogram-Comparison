##########################################################################################################
## To run, do source(bls.R) and source(norm_bls_periodogram.R) if using BLS, and similarly if using TCF ##
##########################################################################################################

######################################### NOTE #################################################
# 1. In this code, tcf periods are scaled by `res` before passing into tcf() since tcf works with cadences rather than absolute time values, unlike BLS.
#    - The period tcf prints on terminal would differ from the period at max(power) shown in the plot.
#    - Such scaling is also done in eva_periodogram.R since TCF calculates periodogram in units of cadences rather than absolute time values.
################################################################################################

# # Below cv function from https://rpubs.com/mengxu/loess_cv
# library(bootstrap)
# loess_wrapper_extrapolate <- function (x, y, span.vals = seq(0.5, 1, by = 0.05), folds = 5){
#   # Do model selection using mean absolute error, which is more robust than squared error.
#   mean.abs.error <- numeric(length(span.vals))
  
#   # Quantify error for each span, using CV
#   loess.model <- function(x, y, span){
#     loess(y ~ x, span = span, control=loess.control(surface="direct"))
#   }
  
#   loess.predict <- function(fit, newdata) {
#     predict(fit, newdata = newdata)
#   }

#   span.index <- 0
#   for (each.span in span.vals) {
#     span.index <- span.index + 1
#     y.hat.cv <- crossval(x, y, theta.fit = loess.model, theta.predict = loess.predict, span = each.span, ngroup = folds)$cv.fit
#     non.empty.indices <- !is.na(y.hat.cv)
#     mean.abs.error[span.index] <- mean(abs(y[non.empty.indices] - y.hat.cv[non.empty.indices]))
#   }
  
#   # find the span which minimizes error
#   best.span <- span.vals[which.min(mean.abs.error)]
  
#   # fit and return the best model
#   best.model <- loess(y ~ x, span = best.span, control=loess.control(surface="direct"))
#   return(best.model)
# }

library(moments)
library('cobs')
source('TCF3.0/intf_libtcf.R')
source('BLS/bls.R')
source('utils.R')
source('eva_periodogram.R')

computeScatter <- function(
    cobsTrendResid,
    windowLength=1000,  # window length on one side of the focal point.
    algo="BLS"
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
            # print(length(cobsTrendResid[l_:u_]))
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
    windowLength=100,
    plot = TRUE,
    ntransits=10,
    ofac=2,
    res=2,  # Light curve resolution, see getLightCurve().
    showFAP = FALSE,  # Whether to show the calculated false alarm probability in the plot. If TRUE, it will take much more time since internally the evd() function is run.
    gaussStd=1e-4,
    ar=0.2,
    ma=0.2,
    order=c(1, 0, 1),
    L=500, R=500,
    useOptimalFreqSampling=FALSE,
    seedValue=1,
    lctype="sim"
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
        # TODO: Once confirmed with Prof. about the cadence of the real light curves, overwrite res value here.
    }

    if (lctype == "sim") {
        # Generate light curve using the parameters.
        yt <- getLightCurve(period, depth, duration, noiseType=noiseType, ntransits=ntransits, res=res, gaussStd=gaussStd, ar=ar, ma=ma, order=order)
        y <- unlist(yt[1])
        t <- unlist(yt[2])
    }
    noiseStd <- unlist(yt[3])
    noiseIQR <- unlist(yt[4])

    if (lctype == "sim") {
        # Special case (TCF fails if absolutely no noise -- so add a very small amount of noise just to prevent any errors).
        if (noiseType == 0 && algo == "TCF") {
            y <- y + 10^-10 * rnorm(length(y))
        }
    }

    # Create frequency grid.
    freqGrid <- getFreqGridToTest(t, period, duration, res=res, ofac=ofac, useOptimalFreqSampling=useOptimalFreqSampling, algo=algo, lctype=lctype)

    stopifnot(exprs={
        all(freqGrid <= res / 2)
    })

    if (algo == "BLS") {
        if (isTRUE(noiseType) == 2) {
            y <- getGPRResid(t, y)
        }
        output <- bls(y, t, bls.plot = FALSE, per.min=min(1/freqGrid), per.max=max(1/freqGrid), nper=length(freqGrid))
    }
    else if (algo == "TCF") {
        fstep <- (max(freqGrid) - min(freqGrid)) / length(freqGrid)
        freqs <- seq(from = min(freqGrid), by = fstep, length.out = length(freqGrid))
        periodsToTry <- 1 / freqs
        residTCF <- getResidForTCF(y)
        output <- tcf(residTCF, p.try = periodsToTry*res, print.output = TRUE)
        print(str(output))
        # output$inper = output$inper / 2
    }

    # (1) Remove trend in periodogram
    # TODO: Is constraint='increase' really needed??
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
        acfEstimate <- acf(y, plot = FALSE)
        lJStats <- Box.test(y, lag = 1, type = "Ljung")  # We want to see autocorrelation with each lag, hence pass lag = 1.
        plot(acfEstimate, main=sprintf("P(Ljung-Box) = %s, lag-1 acf = %s", lJStats[3], acfEstimate$acf[[2]]), cex=2)

        plot(cobsxy50$x, pergram, type = 'l', main=sprintf("Original %s periodogram", algo), log='x', xlab='Period (hrs) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
        lines(cobsxy50$x, cobsxy50$fitted, type = 'l', col='red')
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