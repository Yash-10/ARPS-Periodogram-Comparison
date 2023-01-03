source('test_periodograms.R')
source('standardize_periodogram.R')

blsAndTCF <- function(
    depth, ntransits=10, res=2, ofac=2,
    gaussStd=1e-4, ar=0.2, ma=0.2, order=c(1, 0, 1), useOptimalFreqSampling=TRUE
) {

    period <- 3
    duration <- 2
    L <- 300
    R <- 300

    # Generate light curve using the parameters.
    gyt <- getLightCurve(period, depth, duration, noiseType=1, ntransits=ntransits, res=res, gaussStd=gaussStd, ar=ar, ma=ma, order=order)
    gy <- unlist(gyt[1])
    gt <- unlist(gyt[2])
    noiseStd <- unlist(gyt[3])
    noiseIQR <- unlist(gyt[4])

    # # Special case (TCF fails if absolutely no noise -- so add a very small amount of noise just to prevent any errors).
    # if (noiseType == 0) {
    #     gy <- gy + 10^-10 * rnorm(length(gy))
    # }

    # Create frequency grid.
    gbfreqGrid <- getFreqGridToTest(gt, period, duration, res=res, ofac=ofac, useOptimalFreqSampling=useOptimalFreqSampling, algo="BLS")
    gtfreqGrid <- getFreqGridToTest(gt, period, duration, res=res, ofac=ofac, useOptimalFreqSampling=useOptimalFreqSampling, algo="TCF")

    gboutput <- bls(gy, gt, bls.plot = FALSE, per.min=min(1/gbfreqGrid), per.max=max(1/gbfreqGrid), nper=length(gbfreqGrid))
    gtfstep <- (max(gtfreqGrid) - min(gtfreqGrid)) / length(gtfreqGrid)
    gtfreqs <- seq(from = min(gtfreqGrid), by = gtfstep, length.out = length(gtfreqGrid))
    gtperiodsToTry <- 1 / gtfreqs
    gtresidTCF <- getResidForTCF(gy)
    gtoutput <- tcf(gtresidTCF, p.try = gtperiodsToTry*res, print.output = TRUE)
    # output$inper = output$inper / 2

    # (1) Remove trend in periodogram
    # TODO: Is constraint='increase' really needed??
    gblambdaTrend <- 1
    gbcobsxy50 <- cobs(gboutput$periodsTested, gboutput$spec, ic='BIC', tau=0.5, lambda=gblambdaTrend, constraint="increase")  # If tau = 0.5 and lambda = 0 => Median regression fit.
    gbcobsxy501 <- cobs(gboutput$periodsTested, gboutput$spec, ic='BIC', tau=0.9, lambda=gblambdaTrend)
    gbcobsxy502 <- cobs(gboutput$periodsTested, gboutput$spec, ic='BIC', tau=0.99, lambda=gblambdaTrend)
    gtlambdaTrend <- 1
    gtcobsxy50 <- cobs(gtperiodsToTry, gtoutput$outpow, ic='BIC', tau=0.5, lambda=gtlambdaTrend, constraint="increase")
    gtcobsxy501 <- cobs(gtperiodsToTry, gtoutput$outpow, ic='BIC', tau=0.9, lambda=gtlambdaTrend)
    gtcobsxy502 <- cobs(gtperiodsToTry, gtoutput$outpow, ic='BIC', tau=0.99, lambda=gtlambdaTrend)

    gbperiodogramTrendRemoved <- gbcobsxy50$resid
    gtperiodogramTrendRemoved <- gtcobsxy50$resid

    # (2) Remove local scatter in periodogram
    gscatterWindowLength <- 100
    gbScatter <- computeScatter(gbperiodogramTrendRemoved, windowLength=gscatterWindowLength)
    gtScatter <- computeScatter(gtperiodogramTrendRemoved, windowLength=gscatterWindowLength)

    # print("Scatter")
    gblambdaScatter <- 1
    gbcobsScatter <- cobs(gboutput$periodsTested, gbScatter, ic='BIC', tau=0.5, lambda=gblambdaScatter)
    # cobss50 <- cobs(output$periodsTested, periodogramTrendRemoved, ic='BIC', tau=0.5, lambda=lambdaTrend, constraint="increase")  # If tau = 0.5 and lambda = 0 => Median regression fit.
    # cobss501 <- cobs(output$periodsTested, periodogramTrendRemoved, ic='BIC', tau=0.9, lambda=lambdaTrend, constraint="increase")
    # cobss502 <- cobs(output$periodsTested, periodogramTrendRemoved, ic='BIC', tau=0.99, lambda=lambdaTrend)
    gtlambdaScatter <- 1
    gtcobsScatter <- cobs(gtperiodsToTry, gtScatter, ic='BIC', tau=0.5, lambda=gtlambdaScatter)

    gbnormalizedPeriodogram <- gbperiodogramTrendRemoved / gbcobsScatter$fitted
    gtnormalizedPeriodogram <- gtperiodogramTrendRemoved / gtcobsScatter$fitted

    ###### For AR ######
    # Generate light curve using the parameters.
    ayt <- getLightCurve(period, depth, duration, noiseType=2, ntransits=ntransits, res=res, gaussStd=gaussStd, ar=ar, ma=ma, order=order)
    ay <- unlist(ayt[1])
    at <- unlist(ayt[2])
    noiseStd <- unlist(ayt[3])
    noiseIQR <- unlist(ayt[4])

    # # Special case (TCF fails if absolutely no noise -- so add a very small amount of noise just to prevent any errors).
    # if (noiseType == 0) {
    #     y <- y + 10^-10 * rnorm(length(y))
    # }

    # Create frequency grid.
    abfreqGrid <- getFreqGridToTest(at, period, duration, res=res, ofac=ofac, useOptimalFreqSampling=useOptimalFreqSampling, algo="BLS")
    atfreqGrid <- getFreqGridToTest(at, period, duration, res=res, ofac=ofac, useOptimalFreqSampling=useOptimalFreqSampling, algo="TCF")

    aboutput <- bls(ay, at, bls.plot = FALSE, per.min=min(1/abfreqGrid), per.max=max(1/abfreqGrid), nper=length(abfreqGrid))
    atfstep <- (max(atfreqGrid) - min(atfreqGrid)) / length(atfreqGrid)
    atfreqs <- seq(from = min(atfreqGrid), by = atfstep, length.out = length(atfreqGrid))
    atperiodsToTry <- 1 / atfreqs
    atresidTCF <- getResidForTCF(ay)
    atoutput <- tcf(atresidTCF, p.try = atperiodsToTry*res, print.output = TRUE)
    # output$inper = output$inper / 2

    # (1) Remove trend in periodogram
    # TODO: Is constraint='increase' really needed??
    ablambdaTrend <- 1
    abcobsxy50 <- cobs(aboutput$periodsTested, aboutput$spec, ic='BIC', tau=0.5, lambda=ablambdaTrend, constraint="increase")  # If tau = 0.5 and lambda = 0 => Median regression fit.
    abcobsxy501 <- cobs(aboutput$periodsTested, aboutput$spec, ic='BIC', tau=0.9, lambda=ablambdaTrend)
    abcobsxy502 <- cobs(aboutput$periodsTested, aboutput$spec, ic='BIC', tau=0.99, lambda=ablambdaTrend)
    atlambdaTrend <- 1
    atcobsxy50 <- cobs(atperiodsToTry, atoutput$outpow, ic='BIC', tau=0.5, lambda=atlambdaTrend, constraint="increase")
    atcobsxy501 <- cobs(atperiodsToTry, atoutput$outpow, ic='BIC', tau=0.9, lambda=atlambdaTrend)
    atcobsxy502 <- cobs(atperiodsToTry, atoutput$outpow, ic='BIC', tau=0.99, lambda=atlambdaTrend)

    abperiodogramTrendRemoved <- abcobsxy50$resid
    atperiodogramTrendRemoved <- atcobsxy50$resid

    # (2) Remove local scatter in periodogram
    ascatterWindowLength <- 100
    abScatter <- computeScatter(abperiodogramTrendRemoved, windowLength=ascatterWindowLength)
    atScatter <- computeScatter(atperiodogramTrendRemoved, windowLength=ascatterWindowLength)

    # print("Scatter")
    ablambdaScatter <- 1
    abcobsScatter <- cobs(aboutput$periodsTested, abScatter, ic='BIC', tau=0.5, lambda=ablambdaScatter)
    # cobss50 <- cobs(output$periodsTested, periodogramTrendRemoved, ic='BIC', tau=0.5, lambda=lambdaTrend, constraint="increase")  # If tau = 0.5 and lambda = 0 => Median regression fit.
    # cobss501 <- cobs(output$periodsTested, periodogramTrendRemoved, ic='BIC', tau=0.9, lambda=lambdaTrend, constraint="increase")
    # cobss502 <- cobs(output$periodsTested, periodogramTrendRemoved, ic='BIC', tau=0.99, lambda=lambdaTrend)
    atlambdaScatter <- 1
    atcobsScatter <- cobs(atperiodsToTry, atScatter, ic='BIC', tau=0.5, lambda=atlambdaScatter)

    abnormalizedPeriodogram <- abperiodogramTrendRemoved / abcobsScatter$fitted
    atnormalizedPeriodogram <- atperiodogramTrendRemoved / atcobsScatter$fitted

    ########################################################################

    # Call extreme value analysis code.
    gresultBLS <- evd(period, depth, duration, noiseType=1, algo='BLS', ofac=ofac, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, ar=ar, ma=ma, order=order, FAPSNR_mode=0)
    gresultTCF <- evd(period, depth, duration, noiseType=1, algo='TCF', ofac=ofac, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, ar=ar, ma=ma, order=order, FAPSNR_mode=0)
    aresultBLS <- evd(period, depth, duration, noiseType=2, algo='BLS', ofac=ofac, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, ar=ar, ma=ma, order=order, FAPSNR_mode=0)
    aresultTCF <- evd(period, depth, duration, noiseType=2, algo='TCF', ofac=ofac, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, ar=ar, ma=ma, order=order, FAPSNR_mode=0)
    gfapBLS <- gresultBLS[1]
    gfapTCF <- gresultTCF[1]
    afapBLS <- aresultBLS[1]
    afapTCF <- aresultTCF[1]

    png("combined.png", width=600, height=200, units='mm', res=300)

    par("mar" = c(5, 6, 3, 2))

    cexVal = 2.5
    mat1 <- matrix(c(
        1, 3, 4,
        2, 3, 4,
        5, 7, 8,
        6, 7, 8), nrow = 4, ncol = 3, byrow = TRUE
    )
    layout(mat = mat1,
        heights = c(1),    # Heights of the two rows
        widths = c(1)
    )     # Widths of the two columns

    gbpergram <- gboutput$spec
    gtpergram <- gtoutput$outpow

    plot(gt, gy, type='l', main=sprintf("Period: %.1f days, depth: %.3f (pct), duration: %.1f hrs", period, depth, duration), cex.main=cexVal, cex.lab=2.1, cex.axis=2.0, xlab='time (hrs)', ylab='normalized flux', font.main = 1)
    text(400, 1.0003, "Gaussian", cex=2.0, col='#3182bd')
    # mtext("Gaussian", side=2)
    gacfEstimate <- acf(gy, plot = FALSE)
    glJStats <- Box.test(gy, lag = 1, type = "Ljung")  # We want to see autocorrelation with each lag, hence pass lag = 1.
    gn <- length(gacfEstimate$acf)
    plot(
        gacfEstimate$acf[2:gn], main=sprintf("P(Ljung-Box) = %.3f, lag-1 acf = %.3f", glJStats[3], gacfEstimate$acf[[2]]), cex=1.5, type="h", 
        xlab="Lag",     
        ylab="ACF", 
        ylim=c(-0.2,0.2), # this sets the y scale to -0.2 to 0.2
        las=1,
        xaxt="n",
        cex.main=cexVal, cex.lab=2.1, cex.axis=1.5,
        font.main = 1
    )
    # text(gn/2, 0.15, "Gaussian")
    abline(h=0)
    # Add labels to the x-axis
    x <- c(1:gn)
    y <- c(1:gn)
    axis(1, at=x, labels=y)

    plot(gbcobsxy50$x, gbpergram, type = 'l', main="BLS periodogram", log='x', xlab='Period (hrs) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=2.1, cex.axis=2.0)
    lines(gbcobsxy50$x, gbcobsxy50$fitted, type = 'l', col='red')
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug(gbcobsxy50$knots)
    legend("topleft", lty = 1, 
        col= c("red"), text.col = "black", 
        legend=c("trend fit"), bty="n", cex=2.0, pt.cex = 1
    )
    text(11.5, 5e-5, sprintf("SNR = %.2f, FAP = %.2e", calculateSNR(gboutput$periodsTested, gbpergram), gfapBLS), cex=2.0)
    plot(gtcobsxy50$x, gtpergram, type = 'l', main="TCF periodogram", log='x', xlab='Period (hrs) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=2.1, cex.axis=2.0)
    lines(gtcobsxy50$x, gtcobsxy50$fitted, type = 'l', col='red')
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug(gtcobsxy50$knots)
    legend("topleft", lty = 1, 
        col= c("red"), text.col = "black", 
        legend=c("trend fit"), bty="n", cex=2.0, pt.cex = 1
    )
    text(11.5, 425, sprintf("SNR = %.2f, FAP = %.2e", calculateSNR(gtperiodsToTry * res, gtpergram), gfapTCF), cex=2.0)

    print(calculateSNR(gtperiodsToTry * res, gtpergram))
    print(calculateSNR(gboutput$periodsTested, gbpergram))

    ######################################################
    abpergram <- aboutput$spec
    atpergram <- atoutput$outpow

    plot(at, ay, type='l', main=sprintf("Period: %.1f days, depth: %.3f (pct), duration: %.1f hrs", period, depth, duration), cex.main=cexVal, cex.lab=2.1, cex.axis=2.0, xlab='time (hrs)', ylab='normalized flux', font.main = 1)
    text(400, 1.0005, "Autoregressive", cex=2.0, col='#3182bd')
    aacfEstimate <- acf(ay, plot = FALSE)
    alJStats <- Box.test(ay, lag = 1, type = "Ljung")  # We want to see autocorrelation with each lag, hence pass lag = 1.
    an <- length(aacfEstimate$acf)
    plot(
        aacfEstimate$acf[2:an], main=sprintf("P(Ljung-Box) = %.3f, lag-1 acf = %.3f", alJStats[3], aacfEstimate$acf[[2]]), cex=1.5, type="h", 
        xlab="Lag",     
        ylab="ACF", 
        ylim=c(-0.2,0.2), # this sets the y scale to -0.5 to 0.5
        las=1,
        xaxt="n",
        cex.main=cexVal, cex.lab=2.0, cex.axis=1.2,
        font.main = 1
    )
    abline(h=0)
    # Add labels to the x-axis
    x <- c(1:an)
    y <- c(1:an)
    axis(1, at=x, labels=y)

    plot(abcobsxy50$x, abpergram, type = 'l', main="BLS periodogram", log='x', xlab='Period (hrs) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=2.1, cex.axis=2.0)
    lines(abcobsxy50$x, abcobsxy50$fitted, type = 'l', col='red')
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug(abcobsxy50$knots)
    legend("topleft", lty = 1, 
        col= c("red"), text.col = "black", 
        legend=c("trend fit"), bty="n", cex=2.0, pt.cex = 1
    )
    text(10.5, 4.4e-5, sprintf("SNR = %.2f, FAP = %.2e", calculateSNR(aboutput$periodsTested, abpergram), afapBLS), cex=2.0)
    plot(atcobsxy50$x, atpergram, type = 'l', main="TCF periodogram", log='x', xlab='Period (hrs) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=2.1, cex.axis=2.0)
    lines(atcobsxy50$x, atcobsxy50$fitted, type = 'l', col='red')
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug(atcobsxy50$knots)
    legend("topleft", lty = 1, 
        col= c("red"), text.col = "black", 
        legend=c("trend fit"), bty="n", cex=2.0, pt.cex = 1
    )
    text(11, 137, sprintf("SNR = %.2f, FAP = %.2e", calculateSNR(atperiodsToTry * res, atpergram), afapTCF), cex=2.0)

    dev.off()

    print(calculateSNR(atperiodsToTry * res, atpergram))
    print(calculateSNR(aboutput$periodsTested, abpergram))
}