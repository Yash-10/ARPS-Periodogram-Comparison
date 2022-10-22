blsAndTCF <- function(
    period, depth, duration, noiseType=1, ntransits=10, res=2, ofac=2,
    gaussStd=1e-4, ar=0.2, ma=0.2, order=c(1, 0, 1), showFAP=FALSE, useOptimalFreqSampling=TRUE
) {
    # Generate light curve using the parameters.
    yt <- getLightCurve(period, depth, duration, noiseType=noiseType, ntransits=ntransits, res=res, gaussStd=gaussStd, ar=ar, ma=ma, order=order)
    y <- unlist(yt[1])
    t <- unlist(yt[2])
    noiseStd <- unlist(yt[3])
    noiseIQR <- unlist(yt[4])

    # Special case (TCF fails if absolutely no noise -- so add a very small amount of noise just to prevent any errors).
    if (noiseType == 0) {
        y <- y + 10^-10 * rnorm(length(y))
    }

    # Create frequency grid.
    bfreqGrid <- getFreqGridToTest(t, period, duration, res=res, ofac=ofac, useOptimalFreqSampling=useOptimalFreqSampling, algo="BLS")
    tfreqGrid <- getFreqGridToTest(t, period, duration, res=res, ofac=ofac, useOptimalFreqSampling=useOptimalFreqSampling, algo="TCF")

    boutput <- bls(y, t, bls.plot = FALSE, per.min=min(1/bfreqGrid), per.max=max(1/bfreqGrid), nper=length(bfreqGrid))
    tfstep <- (max(tfreqGrid) - min(tfreqGrid)) / length(tfreqGrid)
    tfreqs <- seq(from = min(tfreqGrid), by = tfstep, length.out = length(tfreqGrid))
    tperiodsToTry <- 1 / tfreqs
    tresidTCF <- getResidForTCF(y)
    toutput <- tcf(residTCF, p.try = tperiodsToTry*res, print.output = TRUE)
    # output$inper = output$inper / 2

    # (1) Remove trend in periodogram
    # TODO: Is constraint='increase' really needed??
    blambdaTrend <- 1
    bcobsxy50 <- cobs(boutput$periodsTested, boutput$spec, ic='BIC', tau=0.5, lambda=blambdaTrend, constraint="increase")  # If tau = 0.5 and lambda = 0 => Median regression fit.
    bcobsxy501 <- cobs(boutput$periodsTested, boutput$spec, ic='BIC', tau=0.9, lambda=blambdaTrend)
    bcobsxy502 <- cobs(boutput$periodsTested, boutput$spec, ic='BIC', tau=0.99, lambda=blambdaTrend)
    tlambdaTrend <- 1
    tcobsxy50 <- cobs(tperiodsToTry, toutput$outpow, ic='BIC', tau=0.5, lambda=tlambdaTrend, constraint="increase")
    tcobsxy501 <- cobs(tperiodsToTry, toutput$outpow, ic='BIC', tau=0.9, lambda=tlambdaTrend)
    tcobsxy502 <- cobs(tperiodsToTry, toutput$outpow, ic='BIC', tau=0.99, lambda=tlambdaTrend)

    bperiodogramTrendRemoved <- bcobsxy50$resid
    tperiodogramTrendRemoved <- tcobsxy50$resid

    # (2) Remove local scatter in periodogram
    scatterWindowLength <- 100
    bScatter <- computeScatter(bperiodogramTrendRemoved, windowLength=scatterWindowLength)
    tScatter <- computeScatter(tperiodogramTrendRemoved, windowLength=scatterWindowLength)

    # print("Scatter")
    blambdaScatter <- 1
    bcobsScatter <- cobs(boutput$periodsTested, bScatter, ic='BIC', tau=0.5, lambda=blambdaScatter)
    # cobss50 <- cobs(output$periodsTested, periodogramTrendRemoved, ic='BIC', tau=0.5, lambda=lambdaTrend, constraint="increase")  # If tau = 0.5 and lambda = 0 => Median regression fit.
    # cobss501 <- cobs(output$periodsTested, periodogramTrendRemoved, ic='BIC', tau=0.9, lambda=lambdaTrend, constraint="increase")
    # cobss502 <- cobs(output$periodsTested, periodogramTrendRemoved, ic='BIC', tau=0.99, lambda=lambdaTrend)
    tlambdaScatter <- 1
    tcobsScatter <- cobs(tperiodsToTry, tScatter, ic='BIC', tau=0.5, lambda=tlambdaScatter)

    bnormalizedPeriodogram <- bperiodogramTrendRemoved / bcobsScatter$fitted
    tnormalizedPeriodogram <- tperiodogramTrendRemoved / tcobsScatter$fitted

    if (showFAP) {
        # Call extreme value analysis code.
        result <- evd(period, depth, duration, noiseType=noiseType, algo=algo, ofac=ofac, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, ar=ar, ma=ma, order=order)
        print(sprintf("FAP = %.10f", result[1]))
        fap <- result[1]
    }

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

    bpergram <- boutput$spec
    tpergram <- toutput$outpow

    plot(t, y, type='l', main=sprintf("Period: %.1f days, depth: %.3f (pct), duration: %.1f (hrs)", period, depth, duration), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xlab='time (hrs)', ylab='normalized flux')
    acfEstimate <- acf(y, plot = FALSE)
    lJStats <- Box.test(y, lag = 1, type = "Ljung")  # We want to see autocorrelation with each lag, hence pass lag = 1.
    plot(acfEstimate, main=sprintf("P(Ljung-Box) = %.3f, lag-1 acf = %.3f", lJStats[3], acfEstimate$acf[[2]]), cex=2)

    plot(bcobsxy50$x, bpergram, type = 'l', main="Original BLS periodogram", log='x', xlab='Period (hrs) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
    lines(bcobsxy50$x, bcobsxy50$fitted, type = 'l', col='red')
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug(bcobsxy50$knots)
    legend(x = "topleft", lty = 1, text.font = 6, 
        col= c("red"), text.col = "black", 
        legend=c("trend fit")
    )
    plot(tcobsxy50$x, tpergram, type = 'l', main="Original TCF periodogram", log='x', xlab='Period (hrs) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
    lines(tcobsxy50$x, tcobsxy50$fitted, type = 'l', col='red')
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug(tcobsxy50$knots)
    legend(x = "topleft", lty = 1, text.font = 6, 
        col= c("red"), text.col = "black", 
        legend=c("trend fit")
    )

    plot.new()  # Just show an empty plot.

    # Plot histogram of periodogram. Shows log-frequency on y-axis in histogram for better visualization.
    bhist.data = hist(bpergram, breaks=50, plot = FALSE)
    # Compute skewness and kurtosis of the original and standardized histograms.
    ### Refer https://brownmath.com/stat/shape.htm for more information ###
    ### Note: R does NOT compute the "excess kurtosis".
    # The kurtosis is calculated as follows:
    # ```
    # n <- length(x)
    # n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
    # ``` Taken from https://stackoverflow.com/a/21484052
    bSkewnessBefore <- skewness(bpergram)
    bKurtosisBefore <- kurtosis(bpergram)

    thist.data = hist(tpergram, breaks=50, plot = FALSE)
    # Compute skewness and kurtosis of the original and standardized histograms.
    ### Refer https://brownmath.com/stat/shape.htm for more information ###
    ### Note: R does NOT compute the "excess kurtosis".
    # The kurtosis is calculated as follows:
    # ```
    # n <- length(x)
    # n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
    # ``` Taken from https://stackoverflow.com/a/21484052
    tSkewnessBefore <- skewness(tpergram)
    tKurtosisBefore <- kurtosis(tpergram)

    if (showFAP) {
        plot(bhist.data$count, type='h', log='y', main=sprintf('Original BLS periodogram histogram, FAP: %s', formatC(fap, format = "e", digits = 5)), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
        axis(1, at=1:length(bhist.data$mids), labels=bhist.data$mids)
        plot(thist.data$count, type='h', log='y', main=sprintf('Original TCF periodogram histogram, FAP: %s', formatC(fap, format = "e", digits = 5)), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
        axis(1, at=1:length(thist.data$mids), labels=thist.data$mids)
    }
    else {
        plot(bhist.data$count, type='h', log='y', main=sprintf('Original BLS periodogram histogram'), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
        axis(1, at=1:length(bhist.data$mids), labels=bhist.data$mids)
        plot(thist.data$count, type='h', log='y', main=sprintf('Original TCF periodogram histogram'), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
        axis(1, at=1:length(thist.data$mids), labels=thist.data$mids)
    }
}