##########################################################################################################
## To run, do source(bls.R) and source(norm_bls_periodogram.R) if using BLS, and similarly if using TCF ##
##########################################################################################################

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

computeScatter <- function(
    cobsTrendResid,
    windowLength=1000,  # window length on one side of the focal point.
    algo="BLS"
){
    scatterVals = c()
    for (i in 1:length(cobsTrendResid)) {
        if (i <= windowLength) {
            u = 2 * windowLength + 1
            l = i + 1
            scatterVal = IQR(c(cobsTrendResid[1:i], cobsTrendResid[l:u]))
            # print(length(c(cobsTrendResid[1:i], cobsTrendResid[l:u])))
            scatterVals <- append(scatterVals, scatterVal)
        }
        else if (i + windowLength >= length(cobsTrendResid)) {
            ll = i - 2 * windowLength - i + length(cobsTrendResid)
            ul = i - 1
            ls = i
            us = length(cobsTrendResid)
            # print(length(c(cobsTrendResid[ll:ul], cobsTrendResid[ls:us])))
            scatterVal = IQR(c(cobsTrendResid[ll:ul], cobsTrendResid[ls:us]))
            scatterVals <- append(scatterVals, scatterVal)
        }
        else {
            l_ = i-windowLength
            u_ = i+windowLength
            # print(length(cobsTrendResid[l_:u_]))
            scatterVal <- IQR(cobsTrendResid[l_:u_])
            scatterVals <- append(scatterVals, scatterVal)
        }
    }
    return (scatterVals)
}

standardPeriodogram <- function(
    y,  # time series values
    t,  # time
    ### NOTE: period, depth, and duration are only taken as inputs for plotting -- not used anywhere else. Remove them if you do not want, but then change the "main" argument to "plot" below.
    period,
    depth,
    duration,
    noiseType = 0,  # 0 for no noise, 1 for white Gaussian noise, and 2 for autoregressive noise (if 2, the ARMA model is internally fixed; also no differencing is used, so it assumes the time-series is stationary or already differenced.)
    algo = "BLS",  # or "TCF"
    windowLength=100,
    plot = TRUE,
    perMin=t[3]-t[1],
    perMax=t[length(t)]-t[1],
    nper=length(t)*10
){
    if (noiseType == 1) {
        set.seed(1)
        y <- y + 0.3 * rnorm(length(y))  # scale by 0.3 to keep noise levels low.
    }
    else if (noiseType == 2) {
        set.seed(1)
        autoRegNoise <- arima.sim(model = list(order=c(1, 0, 1), ar=0.2, ma=0.2), n = length(y))  # It has only AR and MA components, no differencing component. So it is ARMA and not ARIMA. Note: Keep ar and ma < 0.5.
        y <- y + autoRegNoise
    }

    # Run periodogram algorithm.

    # By default, BLS tests these periods:
    # per.min = t[3]-t[1], ## min period to test
    # per.max = t[length(t)]-t[1]
    # f = seq(freq.min,by=freq.step,length.out=nfreq)
    # per = 1/f
    if (algo == "BLS") {
        output <- bls(y, t, bls.plot = FALSE, per.min=perMin, per.max=perMax, nper=nper)
    }
    else {
        perMin = t[3] - t[1]
        perMax = t[length(t)] - t[1]
        freqMax = 1 / perMin
        freqMin = 1 / perMax
        nfreq = length(y) * 10
        freqStep = (freqMax - freqMin) / nfreq
        f = seq(freqMin, by=freqStep, length.out=nfreq)
        periodsToTry = 1 / f
        output <- tcf(y, p.try = periodsToTry)
    }

    # (1) Remove trend in periodogram
    if (algo == "BLS") {
        lambdaTrend <- 1
        cobsxy50 <- cobs(output$periodsTested, output$spec, ic='BIC', tau=0.5, lambda=lambdaTrend, constraint="increase")  # If tau = 0.5 and lambda = 0 => Median regression fit.
        cobsxy501 <- cobs(output$periodsTested, output$spec, ic='BIC', tau=0.9, lambda=lambdaTrend, constraint="increase")
        cobsxy502 <- cobs(output$periodsTested, output$spec, ic='BIC', tau=0.99, lambda=lambdaTrend)
    }
    else if (algo == "TCF") {
        lambdaTrend <- 1
        cobsxy50 <- cobs(periodsToTry, output$outpow, ic='BIC', tau=0.5, lambda=lambdaTrend, constraint="increase")
        cobsxy501 <- cobs(periodsToTry, output$outpow, ic='BIC', tau=0.9, lambda=lambdaTrend, constraint="increase")
        cobsxy502 <- cobs(periodsToTry, output$outpow, ic='BIC', tau=0.99, lambda=lambdaTrend)
    }

    periodogramTrendRemoved <- cobsxy50$resid

    # (2) Remove local scatter in periodogram
    scatterWindowLength <- 3000
    Scatter <- computeScatter(periodogramTrendRemoved, windowLength=scatterWindowLength)
    # print("Scatter")
    if (algo == "BLS") {
        lambdaScatter <- 1
        cobsScatter <- cobs(output$periodsTested, Scatter, ic='BIC', tau=0.5, lambda=lambdaScatter)
    }
    else if (algo == "TCF") {
        lambdaScatter <- 1
        cobsScatter <- cobs(periodsToTry, Scatter, ic='BIC', tau=0.5, lambda=lambdaScatter)
    }

    # plot(output$periodsTested[1:200], Scatter[1:200], type='l')
    # lines(output$periodsTested[1:200], cobsScatter$fitted[1:200], col='red', type='l')
    # return (1);

    normalizedPeriodogram <- periodogramTrendRemoved / cobsScatter$fitted

    if (algo == "BLS") {
        returnVals <- list(normalizedPeriodogram, output$periodsTested)
    }
    else {
        returnVals <- list(normalizedPeriodogram, periodsToTry)
    }

    if (plot) {
        if (algo == "BLS") {
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

            plot(t, y, type='l', main=sprintf("period: %.3f days, depth: %.3f (pct), duration: %.3f (hrs)", period, depth, period * 24 * duration), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xlab='time (hrs)')
            acfEstimate <- acf(y, plot = FALSE)
            lJStats <- Box.test(y, lag = 1, type = "Ljung")  # We want to see autocorrelation with each lag, hence pass lag = 1.
            plot(acfEstimate, main=sprintf("P(Ljung-Box) = %s", lJStats[3]), cex=2)

            plot(cobsxy50$x, output$spec, type = 'l', main=sprintf("Original BLS periodogram and cobs fit, lambdaTrend=%d", lambdaTrend), log='x', xlab='log Period (hrs)', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
            lines(cobsxy50$x, cobsxy50$fitted, type = 'l', col='red')
            lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
            lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
            rug(cobsxy50$knots)
            legend(x = "topleft", lty = 1, text.font = 6, 
                col= c("red", "cyan", "magenta"), text.col = "black", 
                legend=c("trend fit", "tau=0.9", "tau=0.99")
            )

            plot(cobsxy50$x, periodogramTrendRemoved, type = 'l', main=sprintf("BLS periodogram (detrended), lamScatter=%d, windowLength=%d", lambdaScatter, scatterWindowLength), log='x', xlab='log Period', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
            lines(cobsxy50$x, cobsScatter$fitted, type = 'l', col='blue')
            # lines(cobsxy50$x, Scatter, type = 'l', col='red')
            legend(x="topleft", lty=1, text.font = 6,
                col="blue", text.col="black",
                legend="local scatter fit"
            )

            plot(cobsxy50$x, normalizedPeriodogram, type = 'l', main='Standardized BLS periodogram (detrended and local scatter removed)', log='x', xlab='log Period', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)

            # Plot histogram of periodogram. Shows log-frequency on y-axis in histogram for better visualization.
            hist.data = hist(output$spec, breaks=50, plot = FALSE)

            # Compute skewness and kurtosis of the original and standardized histograms.
            ### Refer https://brownmath.com/stat/shape.htm for more information ###
            ### Note: R does NOT compute the "excess kurtosis".
            # The kurtosis is calculated as follows:
            # ```
            # n <- length(x)
            # n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
            # ``` Taken from https://stackoverflow.com/a/21484052
            SkewnessBefore <- skewness(output$spec)
            KurtosisBefore <- kurtosis(output$spec)

            plot(hist.data$count, type='h', log='y', main=sprintf('Original BLS periodogram histogram:\nskewness: %.3f, kurtosis: %.3f', SkewnessBefore, KurtosisBefore), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
            axis(1, at=1:length(hist.data$mids), labels=hist.data$mids)
            
            SkewnessAfter <- skewness(normalizedPeriodogram)
            KurtosisAfter <- kurtosis(normalizedPeriodogram)

            hist.data = hist(normalizedPeriodogram, breaks=50, plot = FALSE)
            plot(hist.data$count, type='h', log='y', main=sprintf("Standardized BLS periodogram histogram:\nskewness: %.3f, kurtosis: %.3f", SkewnessAfter, KurtosisAfter), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
            axis(1, at=1:length(hist.data$mids), labels=hist.data$mids)
            # hist.data = hist(output$spec, breaks=50, plot=FALSE)
            # plot(hist.data$count, type='h', log='y', main='Original BLS periodogram Histogram')
            # hist.data = hist(normalizedBLS, breaks=50, plot=FALSE)
            # plot(hist.data$count, type='h', log='y', main='Standardized BLS periodogram Histogram')
        }
        else {
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

            plot(t, y, type='l', main=sprintf("period: %.3f days, depth: %.3f (pct), duration: %.3f (hrs)", period, depth, period * 24 * duration), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xlab='time (hrs)')
            acfEstimate <- acf(y, plot = FALSE)
            lJStats <- Box.test(y, lag = 1, type = "Ljung")  # We want to see autocorrelation with each lag, hence pass lag = 1.
            plot(acfEstimate, main=sprintf("P(Ljung-Box) = %s", lJStats[3]), cex=2)

            # Note that when I write log period, I mean period on log-scale rather than log of the period.
            plot(cobsxy50$x, output$outpow, type = 'l', main=sprintf("Original TCF periodogram and cobs fit, lambdaTrend=%d", lambdaTrend), log='x', xlab='log Period (hrs)', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
            lines(cobsxy50$x, cobsxy50$fitted, type = 'l', col='red')
            lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
            lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
            rug(cobsxy50$knots)
            legend(x = "topleft", lty = 1, text.font = 6, 
                col= c("red", "cyan", "magenta"), text.col = "black", 
                legend=c("trend fit", "tau=0.9", "tau=0.99")
            )

            plot(cobsxy50$x, periodogramTrendRemoved, type = 'l', main=sprintf("TCF periodogram (detrended), lamScatter=%d, windowLength=%d", lambdaScatter, scatterWindowLength), log='x', xlab='log Period', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
            lines(cobsxy50$x, cobsScatter$fitted, type = 'l', col='blue')
            legend(x="topleft", lty=1, text.font = 6,
                col="blue", text.col="black",
                legend="local scatter fit"
            )

            plot(cobsxy50$x, normalizedPeriodogram, type = 'l', main='Standardized TCF periodogram (detrended and local scatter removed)', log='x', xlab='log Period', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)

            # Plot histogram of periodogram. Shows log-frequency on y-axis in histogram for better visualization.
            hist.data = hist(output$spec, breaks=50, plot = FALSE)

            # Compute skewness and kurtosis of the original and standardized histograms.
            ### Refer https://brownmath.com/stat/shape.htm for more information ###
            sampleSkewnessBefore <- skewness(output$spec)  # Note: I assume, like Excel, R also computes the sample skewness rather than population skewness.
            sampleKurtosisBefore <- kurtosis(output$spec)

            plot(hist.data$count, type='h', log='y', main=sprintf('Original TCF periodogram histogram:\nskewness: %.3f, kurtosis: %.3f', sampleSkewnessBefore, sampleKurtosisBefore), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
            axis(1, at=1:length(hist.data$mids), labels=hist.data$mids)
            
            SkewnessAfter <- skewness(normalizedPeriodogram)
            KurtosisAfter <- kurtosis(normalizedPeriodogram)

            hist.data = hist(normalizedPeriodogram, breaks=50, plot = FALSE)
            plot(hist.data$count, type='h', log='y', main=sprintf("Standardized TCF periodogram histogram:\nskewness: %.3f, kurtosis: %.3f", SkewnessAfter, KurtosisAfter), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
            axis(1, at=1:length(hist.data$mids), labels=hist.data$mids)
            
            # Plot histogram of periodogram. Shows log-frequency on y-axis in histogram for better visualization.
            # hist.data = hist(output$outpow, breaks=50, plot=FALSE)
            # plot(hist.data$count, type='h', log='y', main='Original TCF periodogram Histogram')
            # hist.data = hist(normalizedTCF, breaks=50, plot=FALSE)
            # plot(hist.data$count, type='h', log='y', main='Standardized TCF periodogram Histogram')
        }
    }
    return (returnVals)
}
