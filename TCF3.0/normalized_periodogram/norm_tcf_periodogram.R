#####################################################################
## To run, do source(intf_libtcf.R) and source(norm_tcf_periodogram.R) ##
#####################################################################

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

norm_tcf <- function(
    y,  # time series values
    t,  # time
    noiseType = 0,  # 0 for no noise, 1 for white Gaussian noise, and 2 for autoregressive noise (if 2, the ARMA model is internally fixed; also no differencing is used, so it assumes the time-series is stationary or already differenced.)
    tcfLoessSpan = 0.9,  # degree of smoothing: lower values mean lesser smoothing.
    scatterSpan = 0.5,  # degree of smoothing: lower values mean lesser smoothing.
    useCrossValLoess = FALSE,
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

    # By default, BLS tests these periods, so we currently use the same test periods in TCF:
    # per.min = t[3]-t[1], ## min period to test
    # per.max = t[length(t)]-t[1]
    # f = seq(freq.min,by=freq.step,length.out=nfreq)
    # per = 1/f
    perMin = t[3] - t[1]
    perMax = t[length(t)] - t[1]
    freqMax = 1 / perMin
    freqMin = 1 / perMax
    nfreq = length(y) * 10
    freqStep = (freqMax - freqMin) / nfreq
    f = seq(freqMin, by=freqStep, length.out=nfreq)
    periodsToTry = 1 / f

    ################################################################################
    ## Note: TCF assumes light curve is already differenced using, say, ARIMA/ARFIMA
    ## So, `y` must not have any non-stationarity like trends/volatility.
    ################################################################################
    tcf_output <- tcf(y, p.try = periodsToTry)

    # (1) Remove trend in periodogram
    if (useCrossValLoess) {
        tcfLoess <- loess_wrapper_extrapolate(periodsToTry, tcf_output$outpow)
    }
    else {
        tcfLoess <- loess(tcf_output$outpow ~ periodsToTry, span = tcfLoessSpan)
    }
    tcfTrendRemoved <- tcf_output$outpow - tcfLoess$fitted

    # (2) Remove local scatter in periodogram
    tcfScatter <- abs(diff(tcf_output$outpow))
    # After taking diff, we have one less value, i.e. if time series was of length 100, now diff will be 99 in length. So add a value to match length.
    tcfScatter <- append(tcf_output$outpow[1], tcfScatter)

    if (useCrossValLoess) {
        scatterLoess <- loess_wrapper_extrapolate(tcfLoess$x, tcfScatter)
    }
    else {
        scatterLoess <- loess(tcfScatter ~ tcfLoess$x, span = scatterSpan)
    }

    normalizedTCF <- tcfTrendRemoved / scatterLoess$fitted
    if (plot) {
        cexVal = 1.7  # TODO: Make it work
        par(mfrow=c(2,3))
        plot(t, y, type='l', main='Original time series', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)

        plot(tcfLoess$x, tcf_output$outpow, type = 'l', main='Original TCF periodogram and loess fit', log='x', xlab='log Period', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
        lines(tcfLoess$x, tcfLoess$fitted, type = 'l', col='red')
        lines(tcfLoess$x, scatterLoess$fitted, type = 'l', col='blue')
        legend(x = "topright", lty = c(1, 2), text.font = 6,
            col= c("red", "blue"), text.col = "black",
            legend=c("trend fit", "local scatter fit")
        )

        plot(tcfLoess$x, tcfTrendRemoved, type = 'l', main='TCF periodogram (detrended)', log='x', xlab='log Period', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)

        plot(tcfLoess$x, normalizedTCF, type = 'l', main='Standardized TCF periodogram (detrended and local scatter removed)', log='x', xlab='log Period', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)

        # Plot histogram of periodogram. Shows log-frequency on y-axis in histogram for better visualization.
        hist.data = hist(tcf_output$outpow, breaks=50, plot = FALSE)
        plot(hist.data$count, type='h', log='y', main='Original TCF periodogram Histogram', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
        hist.data = hist(normalizedTCF, breaks=50, plot = FALSE)
        plot(hist.data$count, type='h', log='y', main="Standardized TCF periodogram Histogram", cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)

        # hist.data = hist(tcf_output$outpow, breaks=50, plot=FALSE)
        # plot(hist.data$count, type='h', log='y', main='Original TCF periodogram Histogram')
        # hist.data = hist(normalizedTCF, breaks=50, plot=FALSE)
        # plot(hist.data$count, type='h', log='y', main='Standardized TCF periodogram Histogram')
    }
    return (normalizedTCF)
}
