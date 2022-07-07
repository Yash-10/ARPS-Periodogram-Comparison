#################################################################
## To run, do source(bls.R) and source(norm_bls_periodogram.R) if using BLS, and similarly if using TCF ##
#################################################################

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

computeScatter <- function(
    cobsTrendResid,
    windowLength=1000,  # window length on one side of the focal point.
    algo="BLS"
){
    scatterVals = c()
    if (algo == "BLS") {
        for (i in 1:length(cobsTrendResid)) {
            if (i <= windowLength) {
                u = 2 * windowLength + 1
                l = i + 1
                scatterVal = IQR(c(cobsTrendResid[1:i], cobsTrendResid[l:u]))
                scatterVals <- append(scatterVals, scatterVal)
                next
            }
            else if (i + windowLength >= length(cobsTrendResid)) {
                ll = 2 * windowLength - length(cobsTrendResid) + i - 1
                ls = i
                us = length(cobsTrendResid)
                scatterVal = IQR(c(cobsTrendResid[ll:i-1], cobsTrendResid[ls:us]))
                scatterVals <- append(scatterVals, scatterVal)
                next
            }

            scatterVal <- IQR(cobsTrendResid[i-windowLength:i+windowLength])
            scatterVals <- append(scatterVals, scatterVal)
        }
    }
    else if (algo == "TCF") {
        for (i in 1:length(cobsTrendResid)) {
            if (i <= windowLength) {
                scatterVal = IQR(cobsTrendResid[1:windowLength])
                scatterVals <- append(scatterVals, scatterVal)
                next
            }
            else if (i + windowLength >= length(cobsTrendResid)) {
                scatterVal = IQR(cobsTrendResid[length(cobsTrendResid)-windowLength:length(cobsTrendResid)])
                scatterVals <- append(scatterVals, scatterVal)
                next
            }

            scatterVal <- IQR(cobsTrendResid[i-windowLength:i+windowLength])
            scatterVals <- append(scatterVals, scatterVal)
        }
    }
    return (scatterVals)
}

standardPeriodogram <- function(
    y,  # time series values
    t,  # time
    noiseType = 0,  # 0 for no noise, 1 for white Gaussian noise, and 2 for autoregressive noise (if 2, the ARMA model is internally fixed; also no differencing is used, so it assumes the time-series is stationary or already differenced.)
    algo = "BLS",  # or "TCF"
    windowLength=100,
    plot = TRUE
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
        output <- bls(y, t, bls.plot = FALSE)
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
        cobsxy50 <- cobs(output$periodsTested, output$spec, ic='BIC', tau=0.5)  # Median regression fit.
    }
    else {
        cobsxy50 <- cobs(periodsToTry, output$outpow, ic='BIC', tau=0.5)
    }

    if (algo == "BLS") {
        periodogramTrendRemoved <- cobsxy50$resid
    }
    else {
        periodogramTrendRemoved <- cobsxy50$resid
    }

    # (2) Remove local scatter in periodogram
    Scatter <- computeScatter(cobsxy50$resid, windowLength=100)

    normalizedPeriodogram <- periodogramTrendRemoved / Scatter

    if (algo == "BLS") {
        returnVals <- c(normalizedPeriodogram, output$periodsTested)
    }
    else {
        returnVals <- c(normalizedPeriodogram, periodsToTry)
    }

    if (plot) {
        if (algo == "BLS") {
            cexVal = 1.7
            par(mfrow=c(2,3))
            plot(t, y, type='l', main='Original time series', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)

            plot(cobsxy50$x, output$spec, type = 'l', main='Original BLS periodogram and cobs fit', log='x', xlab='log Period', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
            lines(cobsxy50$x, cobsxy50$fitted, type = 'l', col='red')
            print(length(cobsxy50$x))
            print(length(Scatter))
            print(length(periodogramTrendRemoved))
            lines(cobsxy50$x, Scatter, type = 'l', col='blue')
            legend(x = "topright", lty = c(1, 2), text.font = 6, 
                col= c("red", "blue"), text.col = "black", 
                legend=c("trend fit", "local scatter fit")
            )

            plot(cobsxy50$x, periodogramTrendRemoved, type = 'l', main='BLS periodogram (detrended)', log='x', xlab='log Period', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)

            plot(cobsxy50$x, normalizedPeriodogram, type = 'l', main='Standardized BLS periodogram (detrended and local scatter removed)', log='x', xlab='log Period', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)

            # Plot histogram of periodogram. Shows log-frequency on y-axis in histogram for better visualization.
            hist.data = hist(output$spec, breaks=50, plot = FALSE)
            plot(hist.data$count, type='h', log='y', main='Original BLS periodogram Histogram', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
            hist.data = hist(normalizedPeriodogram, breaks=50, plot = FALSE)
            plot(hist.data$count, type='h', log='y', main="Standardized BLS periodogram Histogram", cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)

            # hist.data = hist(output$spec, breaks=50, plot=FALSE)
            # plot(hist.data$count, type='h', log='y', main='Original BLS periodogram Histogram')
            # hist.data = hist(normalizedBLS, breaks=50, plot=FALSE)
            # plot(hist.data$count, type='h', log='y', main='Standardized BLS periodogram Histogram')
        }
        else {
            cexVal = 1.7  # TODO: Make it work
            par(mfrow=c(2,3))
            plot(t, y, type='l', main='Original time series', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)

            plot(cobsxy50$x, output$outpow, type = 'l', main='Original TCF periodogram and cobs fit', log='x', xlab='log Period', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
            lines(cobsxy50$x, cobsxy50$fitted, type = 'l', col='red')
            lines(cobsxy50$x, Scatter, type = 'l', col='blue')
            legend(x = "topright", lty = c(1, 1), text.font = 6,
                col= c("red", "blue"), text.col = "black",
                legend=c("trend fit", "local scatter fit"),
            )

            plot(cobsxy50$x, periodogramTrendRemoved, type = 'l', main='TCF periodogram (detrended)', log='x', xlab='log Period', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)

            plot(cobsxy50$x, normalizedPeriodogram, type = 'l', main='Standardized TCF periodogram (detrended and local scatter removed)', log='x', xlab='log Period', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)

            # Plot histogram of periodogram. Shows log-frequency on y-axis in histogram for better visualization.
            hist.data = hist(output$outpow, breaks=50, plot = FALSE)
            plot(hist.data$count, type='h', log='y', main='Original TCF periodogram Histogram', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
            hist.data = hist(normalizedPeriodogram, breaks=50, plot = FALSE)
            plot(hist.data$count, type='h', log='y', main="Standardized TCF periodogram Histogram", cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)

            # hist.data = hist(output$outpow, breaks=50, plot=FALSE)
            # plot(hist.data$count, type='h', log='y', main='Original TCF periodogram Histogram')
            # hist.data = hist(normalizedTCF, breaks=50, plot=FALSE)
            # plot(hist.data$count, type='h', log='y', main='Standardized TCF periodogram Histogram')
        }
    }
    return (returnVals)
}
