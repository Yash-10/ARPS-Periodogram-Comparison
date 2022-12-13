blsAndTCF <- function(
    depth, noiseType=1, ntransits=10, res=2, ofac=2,
    gaussStd=1e-4, ar=0.2, ma=0.2, order=c(1, 0, 1), useOptimalFreqSampling=TRUE
) {

    period <- 1
    duration <- 2

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
    toutput <- tcf(tresidTCF, p.try = tperiodsToTry*res, print.output = TRUE)
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

    # Call extreme value analysis code.
    resultBLS <- evd(period, depth, duration, noiseType=noiseType, algo='BLS', ofac=ofac, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, ar=ar, ma=ma, order=order, FAPSNR_mode=0)
    resultTCF <- evd(period, depth, duration, noiseType=noiseType, algo='TCF', ofac=ofac, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, ar=ar, ma=ma, order=order, FAPSNR_mode=0)
    fapBLS <- resultBLS[1]
    fapTCF <- resultTCF[1]

    par("mar" = c(5, 6, 4, 2), cex=15)

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
    n <- length(acfEstimate$acf)
    plot(
        acfEstimate$acf[2:n], main=sprintf("P(Ljung-Box) = %.3f, lag-1 acf = %.3f", lJStats[3], acfEstimate$acf[[2]]), cex=2, type="h", 
        xlab="Lag",     
        ylab="ACF", 
        ylim=c(-0.2,0.2), # this sets the y scale to -0.2 to 0.2
        las=1,
        xaxt="n"
    )
    abline(h=0)
    # Add labels to the x-axis
    x <- c(1:n)
    y <- c(1:n)
    axis(1, at=x, labels=y)

    plot(bcobsxy50$x, bpergram, type = 'l', main="BLS periodogram", log='x', xlab='Period (hrs) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
    lines(bcobsxy50$x, bcobsxy50$fitted, type = 'l', col='red')
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug(bcobsxy50$knots)
    legend("topleft", lty = 1, 
        col= c("red"), text.col = "black", 
        legend=c("trend fit"), bty="n"
    )
    text(100, 2.65e-5, sprintf("SNR = %.2f, FAP = %.2e", calculateSNR(boutput$periodsTested, bpergram), fapBLS), cex=1.2)
    plot(tcobsxy50$x, tpergram, type = 'l', main="TCF periodogram", log='x', xlab='Period (hrs) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
    lines(tcobsxy50$x, tcobsxy50$fitted, type = 'l', col='red')
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug(tcobsxy50$knots)
    legend("topleft", lty = 1, 
        col= c("red"), text.col = "black", 
        legend=c("trend fit"), bty="n"
    )
    text(100, 62, sprintf("SNR = %.2f, FAP = %.2e", calculateSNR(tperiodsToTry * res, tpergram), fapTCF), cex=1.2)

    print(calculateSNR(tperiodsToTry * res, tpergram))
    print(calculateSNR(boutput$periodsTested, bpergram))
}