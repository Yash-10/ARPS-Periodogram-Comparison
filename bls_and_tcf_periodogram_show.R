#################################################################################################
# bls_and_tcf_periodogram_show.R : This plots the BLS and TCF periodograms in a single plot for
# comparison.
#
# Notes:
# 1. **Important**
#    - The arguments ar, ma, and order are obsolete. The getLightCurve function, to which these
#      were intended to passed have these values harcoded inside that function. So passing
#      different values will not have any difference.
#################################################################################################

source('BLS/bls.R')
source('TCF3.0/intf_libtcf.R')
source('utils.R')

blsAndTCF <- function(
    period=NULL, depth=NULL, duration=NULL, y=NULL, t=NULL, noiseType=1, ntransits=10, res=2, ofac=2,
    gaussStd=1e-4, ar=0.2, ma=0.2, order=c(1, 0, 1), useOptimalFreqSampling=TRUE, lctype="sim",
    applyGPRforBLS=FALSE, seedValue=42, L=300, R=300
) {
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
        yt <- getLightCurve(period, depth, duration, noiseType=noiseType, ntransits=ntransits, res=res, gaussStd=gaussStd, ar=ar, ma=ma, order=order)
        y <- unlist(yt[1])
        t <- unlist(yt[2])
        noiseStd <- unlist(yt[3])
        noiseIQR <- unlist(yt[4])
    }

    if (lctype == "sim") {
        # Special case (TCF fails if absolutely no noise -- so add a very small amount of noise just to prevent any errors).
        if (noiseType == 0) {
            y <- y + 10^-10 * rnorm(length(y))
        }
    }

    # The BLS code does not handle non-numeric/missing values (NA, NaN, Inf).
    # So we need to get rid of all such observations and the corresponding times.
    # If there are non-numeric values, the y and t vectors will be different for BLS and TCF.
    # This also means the frequency grid (that is based on the observation times) will also be different.
    na_check <- any(is.na(y))
    y_BLS <- y[!is.na(y)]
    t_BLS <- t[!is.na(y)]

    if (isTRUE(noiseType == 2) | applyGPRforBLS) {
        y_BLS <- getGPRResid(t_BLS, y_BLS)  # Run Gaussian Processes Regression on light curve if autoregressive noise is present. Only in the case of BLS.
    }

    # Create frequency grid.
    bfreqGrid <- getFreqGridToTest(
        t_BLS, period, duration, res=res, ofac=ofac,
        useOptimalFreqSampling=useOptimalFreqSampling, algo="BLS", lctype=lctype
    )
    tfreqGrid <- getFreqGridToTest(
        t, period, duration, res=res, ofac=ofac,
        useOptimalFreqSampling=useOptimalFreqSampling, algo="TCF", lctype=lctype
    )

    # if (!na_check) {
    #     stopifnot(exprs={
    #         identical(bfreqGrid, tfreqGrid)
    #     })
    # }
    
    boutput <- bls(y_BLS, t_BLS, bls.plot = FALSE, per.min=min(1/bfreqGrid), per.max=max(1/bfreqGrid), nper=length(bfreqGrid))
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
    if (lctype == 'sim') {
        resultBLS <- evd(period, depth, duration, noiseType=noiseType, algo='BLS', ofac=ofac, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, ar=ar, ma=ma, order=order, FAPSNR_mode=0, seedValue=seedValue, lctype=lctype)
        resultTCF <- evd(period, depth, duration, noiseType=noiseType, algo='TCF', ofac=ofac, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, ar=ar, ma=ma, order=order, FAPSNR_mode=0, seedValue=seedValue, lctype=lctype)
        fapBLS <- resultBLS[1]
        fapTCF <- resultTCF[1]
    }
    else if (lctype == 'real') {
        resultBLS <- evd(y=y_BLS, t=t_BLS, noiseType=noiseType, algo='BLS', ofac=ofac, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, ar=ar, ma=ma, order=order, FAPSNR_mode=0, seedValue=seedValue, lctype=lctype)
        resultTCF <- evd(y=y, t=t, noiseType=noiseType, algo='TCF', ofac=ofac, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, ar=ar, ma=ma, order=order, FAPSNR_mode=0, seedValue=seedValue, lctype=lctype)
        fapBLS <- resultBLS[1]
        fapTCF <- resultTCF[1]
    }

    png(filename="real_4.png", width = 500, height = 250, units='mm', res = 300)

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

    plot(t, y, type='l', main= if (lctype == 'sim') sprintf("Period: %.1f days, depth: %.3f (pct), duration: %.1f (hrs)", period, depth, duration) else 'TESS 4', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xlab='time (hrs)', ylab='normalized flux')
    acfEstimate <- acf(y, plot = FALSE, na.action = na.pass)
    lJStats <- Box.test(y, lag = 1, type = "Ljung")  # We want to see autocorrelation with each lag, hence pass lag = 1.
    n <- length(acfEstimate$acf)
    # plot(
    #     acfEstimate$acf[2:n], main=sprintf("P(Ljung-Box) = %.3f, lag-1 acf = %.3f", lJStats[3], acfEstimate$acf[[2]]), cex=2, type="h", 
    #     xlab="Lag",     
    #     ylab="ACF", 
    #     ylim=c(-0.2,0.2), # this sets the y scale to -0.2 to 0.2
    #     las=1,
    #     xaxt="n"
    # )
    # abline(h=0)
    # # Add labels to the x-axis
    # x <- c(1:n)
    # yy <- c(1:n)
    # axis(1, at=x, labels=yy)

    plot(acfEstimate, main=sprintf("P(Ljung-Box) = %.3f, lag-1 acf = %.3f", lJStats[3], acfEstimate$acf[[2]]), cex=2)

    plot(bcobsxy50$x, bpergram, type = 'l', main="BLS periodogram", log='x', xlab='Period (hrs) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
    lines(bcobsxy50$x, bcobsxy50$fitted, type = 'l', col='red', lwd=3.0)
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug(bcobsxy50$knots)
    # legend("topleft", lty = 1, 
    #     col= c("red"), text.col = "black", 
    #     legend=c("trend fit"), bty="n"
    # )
    # text(100, 100, "Phase-folded light curve", cex=1.2)

    text(7, 1.85e-4, sprintf("SNR = %.2f, FAP = %.2e", calculateSNR(boutput$periodsTested, bpergram), fapBLS), cex=1.5)

    plot(tcobsxy50$x, tpergram, type = 'l', main="TCF periodogram", log='x', xlab='Period (hrs) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
    lines(tcobsxy50$x, tcobsxy50$fitted, type = 'l', col='red', lwd=3.0)
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug(tcobsxy50$knots)
    # legend("topleft", lty = 1, 
    #     col= c("red"), text.col = "black", 
    #     legend=c("trend fit"), bty="n"
    # )

    text(7, 125, sprintf("SNR = %.2f, FAP = %.2e", calculateSNR(tperiodsToTry * res, tpergram), fapTCF), cex=1.5)

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

    plot(bhist.data$count, type='h', log='y', main=sprintf('BLS periodogram histogram'), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
    axis(1, at=1:length(bhist.data$mids), labels=bhist.data$mids)
    plot(thist.data$count, type='h', log='y', main=sprintf('TCF periodogram histogram'), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
    axis(1, at=1:length(thist.data$mids), labels=thist.data$mids)
    tKurtosisBefore <- kurtosis(tpergram)

    print(calculateSNR(tperiodsToTry * res, tpergram))
    print(calculateSNR(boutput$periodsTested, bpergram))

    dev.off()
}

# 7, 5.2e-4; 7, 142
# 12, 4.8e-5; 7, 370
# 300, 1.6e-4; 300, 283
# 