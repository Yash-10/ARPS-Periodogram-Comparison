# ***************************************************************************
# Author: Yash Gondhalekar  Last updated: March, 2023

# Description: Given a transit depth, run the extreme value code and related
#              procedures and then plots the results.

#       Notes:
#          1. This script can only be used if simulation is the intention.
#          2. This script was used for making Fig. 1 and 2 in our paper.
# ***************************************************************************

library(forecast)
library('cobs')
library('boot')
library('extRemes')
library('kernlab')
library('goftest')  # install.packages("goftest")

source('BLS/bls.R')
source('TCF3.0/intf_libtcf.R')
source('utils.R')


# Use this function for simulated data.
# It requires only the depth. Period and transit duration are currently hardcoded inside the function.
# This can be helpful if one wants to look at the periodograms by only varying the depth, keeping everything else the same.
# In such a case, this sunction can be run multiple times with different depths.
# Much of the code inside this and the blsAndTCFDepthChangeReal functions use from the code inside `evd` from `eva_periodogram.R` script.
# A future TODO is to factor out the common code into a separate function to prevent duplicacy.
blsAndTCFDepthChange <- function(
    depth, noiseType=1, ntransits=10, res=2, ofac=2,
    gaussStd=1e-4, ar=0.2, ma=0.2, order=c(1, 0, 1), useOptimalFreqSampling=TRUE,
    L=300, R=300, K=2, seedValue=465, applyGPRforBLS=FALSE
) {

    set.seed(seedValue)

    period <- 2
    duration <- 2

    # Generate light curve using the parameters.
    yt <- getLightCurve(period, depth, duration, noiseType=noiseType, ntransits=ntransits, res=res, gaussStd=gaussStd, ar=ar, ma=ma, order=order, seedValue=seedValue)
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

    boutput <- bls(if (noiseType == 2 | applyGPRforBLS) getGPRResid(t, y) else y, t, bls.plot = FALSE, per.min=min(1/bfreqGrid), per.max=max(1/bfreqGrid), nper=length(bfreqGrid))
    bperResults <- c(boutput$per, boutput$depth, boutput$dur)
    tfstep <- (max(tfreqGrid) - min(tfreqGrid)) / length(tfreqGrid)
    tfreqs <- seq(from = min(tfreqGrid), by = tfstep, length.out = length(tfreqGrid))
    tperiodsToTry <- 1 / tfreqs
    tresidTCF <- getResidForTCF(y)
    toutput <- tcf(tresidTCF, p.try = tperiodsToTry * res, print.output = TRUE)
    tpowmax.loc = which.max(toutput$outpow)
    tperResults <- c(toutput$inper[tpowmax.loc]/res, toutput$outdepth[tpowmax.loc], toutput$outdur[tpowmax.loc]/res)
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
    resultBLS <- evd(period, depth, duration, noiseType=noiseType, algo='BLS', ofac=ofac, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, ar=ar, ma=ma, order=order, FAPSNR_mode=0, seedValue=seedValue)
    resultTCF <- evd(period, depth, duration, noiseType=noiseType, algo='TCF', ofac=ofac, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, ar=ar, ma=ma, order=order, FAPSNR_mode=0, seedValue=seedValue)
    fapBLS <- resultBLS[1]
    fapTCF <- resultTCF[1]
    perResultsBLS <- resultBLS[2:4]
    perResultsTCF <- resultTCF[2:4]

    # par("mar" = c(5, 6, 4, 2), cex=15)
    png(filename="depth_change_ar_0.026.png", width = 420, height = 210, units='mm', res = 300)
    par(mar=c(5,6,4,2), cex=15)

    cexVal <- 2.0
    mat1 <- matrix(c(
        1, 3, 4,
        2, 3, 4,
        5, 6, 7,
        5, 6, 7
        ), nrow = 4, ncol = 3, byrow = TRUE
    )
    layout(matrix(c(1,1,1,1,2,2, 3,3,3,4,4,4, 5,5,5,6,6,6), nrow=3, ncol=6, byrow=TRUE))
    # layout(mat = mat1,
    #     heights = c(1, 1, 1, 1),    # Heights of the two rows
    #     widths = c(1.5, 2, 2)
    # )     # Widths of the two columns

    bpergram <- boutput$spec
    tpergram <- toutput$outpow

    plot(t/24, y, type='l', main=sprintf("Period: %.1f days, Depth: %.1f ppm, Duration: %.1f hrs | Noise: Autoregressive", period, depth*1e4, duration), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xlab='Time (days)', ylab='Flux')
    # acfEstimate <- acf(y, plot = FALSE)
    # lJStats <- Box.test(y, lag = 1, type = "Ljung")  # We want to see autocorrelation with each lag, hence pass lag = 1.
    # n <- length(acfEstimate$acf)
    # plot(
    #     acfEstimate$acf[2:n], main=sprintf("P(Ljung-Box) = %.2f, ACF(1) = %.2f", lJStats[3], acfEstimate$acf[[2]]), cex=2, type="h", 
    #     xlab="Lag",     
    #     ylab="ACF", 
    #     ylim=c(-0.2,0.2), # this sets the y scale to -0.2 to 0.2
    #     las=1,
    #     xaxt="n"
    # )
    # abline(h=0)
    # # Add labels to the x-axis
    # x <- c(1:n)
    # y <- c(1:n)
    # axis(1, at=x, labels=y)

    acfEstimate <- acf(y, plot = FALSE, na.action = na.pass)
    lJStats <- Box.test(y, lag = 1, type = "Ljung")  # We want to see autocorrelation with each lag, hence pass lag = 1.
    plot(acfEstimate, main="", cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xlim=c(1, 20), ylim=c(-0.2, +0.5))
    text(10, 0.36, sprintf("P(Ljung-Box) = %.2f, ACF(1) = %.2f\n", lJStats[3], acfEstimate$acf[[2]]), cex=1.9)
    # plot(acfEstimate, main=sprintf("P(Ljung-Box) = %.2f, ACF(1) = %.2f", lJStats[3], acfEstimate$acf[[2]]), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, cex=cexVal, xlim=c(1, 20), ylim=(-0.2, +0.5))

    plot(bcobsxy50$x/24, bpergram, type = 'l', main="BLS periodogram", log='x', xlab='Period (days) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
    lines(bcobsxy50$x/24, bcobsxy50$fitted, type = 'l', col='red', lwd=3.0)
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    print(bcobsxy50$knots)
    rug(bcobsxy50$knots/24)
    # legend("topleft", lty = 1, 
    #     col= c("red"), text.col = "black", 
    #     legend=c("trend fit"), bty="n", cex=1.5, pt.cex = 1
    # )
    text(0.08, 5.4e-5, paste0(sprintf("SNR = %.1f, FAP = %.1e", calculateSNR(boutput$periodsTested, bpergram), fapBLS)), cex=1.9, adj=0)
    text(3.5, 1e-5, paste0(sprintf("Period = %.5f days\nDepth = %.1f", perResultsBLS[1]/24, perResultsBLS[2]*1e6), " ppm"), cex=1.9, adj=0)
    plot(tcobsxy50$x/24, tpergram, type = 'l', main="TCF periodogram", log='x', xlab='Period (days) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
    lines(tcobsxy50$x/24, tcobsxy50$fitted, type = 'l', col='red', lwd=3.0)
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug(tcobsxy50$knots/24)
    # legend("topleft", lty = 1, 
    #     col= c("red"), text.col = "black", 
    #     legend=c("trend fit"), bty="n", cex=1.5, pt.cex = 1
    # )
    text(0.08, 105, paste0(sprintf("SNR = %.1f, FAP = %.1e\nPeriod = %.5f days, Depth = %.1f", calculateSNR(tperiodsToTry * res, tpergram), fapTCF, perResultsTCF[1]/24, perResultsTCF[2]*1e6), " ppm"), cex=1.9, adj=0)

    print(calculateSNR(tperiodsToTry * res, tpergram))
    print(calculateSNR(boutput$periodsTested, bpergram))

    # Plot histogram of detrended periodogram. Shows log-frequency on y-axis in histogram for better visualization.
    bhist.data = hist(bperiodogramTrendRemoved, breaks=50, plot = FALSE)
    # Compute skewness and kurtosis of the original and standardized histograms.
    ### Refer https://brownmath.com/stat/shape.htm for more information ###
    ### Note: R does NOT compute the "excess kurtosis".
    # The kurtosis is calculated as follows:
    # ```
    # n <- length(x)
    # n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
    # ``` Taken from https://stackoverflow.com/a/21484052
    bSkewnessBefore <- skewness(bperiodogramTrendRemoved)
    bKurtosisBefore <- kurtosis(bperiodogramTrendRemoved)

    thist.data = hist(tperiodogramTrendRemoved, breaks=50, plot = FALSE)
    # Compute skewness and kurtosis of the original and standardized histograms.
    ### Refer https://brownmath.com/stat/shape.htm for more information ###
    ### Note: R does NOT compute the "excess kurtosis".
    # The kurtosis is calculated as follows:
    # ```
    # n <- length(x)
    # n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
    # ``` Taken from https://stackoverflow.com/a/21484052
    tSkewnessBefore <- skewness(tperiodogramTrendRemoved)

    plot(bhist.data$count, type='h', log='y', main=sprintf('BLS periodogram (detrended) histogram'), cex=cexVal, cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
    axis(1, at=1:length(bhist.data$mids), labels=sprintf(bhist.data$mids, fmt="%.1e"), cex.axis=cexVal)
    plot(thist.data$count, type='h', log='y', main=sprintf('TCF periodogram (detrended) histogram'), cex=cexVal, cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
    axis(1, at=1:length(thist.data$mids), labels=sprintf(thist.data$mids, fmt="%.1e"), cex.axis=cexVal)
    tKurtosisBefore <- kurtosis(tperiodogramTrendRemoved)

    # dev.print(png, 'depth_change_gaussian_0.008.png')
    dev.off()
}


# MDD for Gaussian case (465 seed value, period=2days, dur=2hrs):
# 0.006836735% for BLS and 0.009591837 for TCF
# MDD for Autoregressive case (544 seed value, period=2days, dur=2hrs):
# 0.02336735 for BLS and 0.02612245 for TCF

# Text coordinates used:
# depth_change_gaussian_0.006.png (BLS, TCF): (22, 1.7155e-5), (22, 45)  # BUT USE adj=1 in text().
# depth_change_gaussian_0.009.png: (0.08, 2.1e-5), (0.08, 60)
# depth_change_gaussian_0.02.png: (0.08, 3.8e-5), (0.08, 153)

# depth_change_ar_0.023.png: (BLS, TCF): (0.08, 4.85e-5 + 3.5, 1e-5), (0.08, 86)
# depth_change_ar_0.026.png: (BLS, TCF): (6.2/24, 5.3e-5), (6.2/24, 115)
# depth_change_ar_0.04.png: (BLS, TCF): (0.08, 7.5e-5), (0.08, 215)


# Rough:
# 1. Output with Autoregressive noise case: [Used seed 544]

#"""
#[1] "Starting with period (FAPSNR_mode = 0) = 2.000000 days"
#[1] "Seed 797.000000: Limiting depths for BLS and TCF:"
#[1] 0.03438776 0.01785714
#[1] "Seed 91.000000: Limiting depths for BLS and TCF:"
#[1] 0.02428571 0.02244898
#[1] "Seed 164.000000: Limiting depths for BLS and TCF:"
#[1] 0.03071429 0.03071429
#[1] "Seed 544.000000: Limiting depths for BLS and TCF:"
#[1] 0.02336735 0.02612245
#[1] "Seed 479.000000: Limiting depths for BLS and TCF:"
#[1] 0.02061224 0.02612245
#[1] "Seed 51.000000: Limiting depths for BLS and TCF:"
#[1] 0.02153061 0.02887755
#[1] "Seed 119.000000: Limiting depths for BLS and TCF:"
#[1] 0.02244898 0.02153061
#[1] "Seed 465.000000: Limiting depths for BLS and TCF:"
#[1] 0.01602041 0.04540816
#[1] "Seed 229.000000: Limiting depths for BLS and TCF:"
#[1] 0.01969388 0.02612245
#[1] "Seed 995.000000: Limiting depths for BLS and TCF:"
#[1] 0.02428571 0.01693878
#"""

# Average: 0.023734694, 0.026214286

# 2. Output with Gaussian noise: [Used seed 465]

#[1] "Starting with period (FAPSNR_mode = 0) = 2.000000 days"
#[1] "Seed 797.000000: Limiting depths for BLS and TCF:"
#[1] 0.005918367 0.007755102
#[1] "Seed 91.000000: Limiting depths for BLS and TCF:"
#[1] 0.005918367 0.007755102
#[1] "Seed 164.000000: Limiting depths for BLS and TCF:"
#[1] 0.005918367 0.009591837
#[1] "Seed 544.000000: Limiting depths for BLS and TCF:"
#[1] 0.007755102 0.009591837
#[1] "Seed 479.000000: Limiting depths for BLS and TCF:"
#[1] 0.007755102 0.014183673
#[1] "Seed 51.000000: Limiting depths for BLS and TCF:"
#[1] 0.01051020 0.01785714
#[1] "Seed 119.000000: Limiting depths for BLS and TCF:"
#[1] 0.009591837 0.013265306
#[1] "Seed 465.000000: Limiting depths for BLS and TCF:"
#[1] 0.006836735 0.009591837
#[1] "Seed 229.000000: Limiting depths for BLS and TCF:"
#[1] 0.007755102 0.013265306
#[1] "Seed 995.000000: Limiting depths for BLS and TCF:"
#[1] 0.005000000 0.007755102

# Average: 0.00729591836734694, 0.0110612244897959


# Use this function on real light curves.
# It requires the argument `table`, that contains fluxes and times as columns.
# See real_light_curve_application.R for an example to setup the table.
blsAndTCFDepthChangeReal <- function(
    table, noiseType=1, ntransits=10, res=2, ofac=2,
    gaussStd=1e-4, ar=0.2, ma=0.2, order=c(1, 0, 1), useOptimalFreqSampling=TRUE,
    L=300, R=300, K=2, seedValue=465, applyGPRforBLS=TRUE
) {
    period <- depth <- duration <- noiseType <- ntransits <- ar <- ma <- order <- gaussStd <- NULL
    significanceMode <- 'max'  # Since for real light curves, passing `expected_peak` is not possible.
    res <- 2

    rtBLS <- c()
    rtTCF <- c()

    table_BLS <- table[!is.na(table$Flux),]

    resultBLS <- evd(y=table_BLS$Flux, t=table_BLS$times, algo="BLS", FAPSNR_mode=0, lctype="real", applyGPRforBLS=TRUE, noiseType=noiseType, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, ar=ar, ma=ma, order=order, seedValue=seedValue)
    fapBLS <- resultBLS[1]
    perResultsBLS <- resultBLS[2:4]
    resultBLS <- evd(y=table_BLS$Flux, t=table_BLS$times, algo="BLS", FAPSNR_mode=1, lctype="real", applyGPRforBLS=TRUE, noiseType=noiseType, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, ar=ar, ma=ma, order=order, seedValue=seedValue)
    snrBLS <- 1 / resultBLS[1]
    rtBLS <- c(rtBLS, c(fapBLS, snrBLS, perResultsBLS))

    resultTCF <- evd(y=table$Flux, t=table$times, algo="TCF", FAPSNR_mode=0, lctype="real", noiseType=noiseType, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, ar=ar, ma=ma, order=order, seedValue=seedValue)
    fapTCF <- resultTCF[1]
    perResultsTCF <- resultTCF[2:4]
    resultTCF <- evd(y=table$Flux, t=table$times, algo="TCF", FAPSNR_mode=1, lctype="real", noiseType=noiseType, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, ar=ar, ma=ma, order=order, seedValue=seedValue)
    snrTCF <- 1 / resultTCF[1]
    rtTCF <- c(rtTCF, c(fapTCF, snrTCF, perResultsTCF))

    y_BLS <- table_BLS$Flux
    t_BLS <- table_BLS$times
    y <- table$Flux
    t <- table$times

    # Create a frequency grid.
    # Note: Optimal frequency sampling is NOT used for real light curves since the period is unknown.
    bfreqGrid <- getFreqGridToTest(t_BLS, period, duration, res=res, ofac=ofac, algo=algo, lctype="real")
    if (any(is.na(bfreqGrid))) {
        stop("Atleast one frequency in the frequency grid is NaN!")
    }

    stopifnot(exprs={
        all(bfreqGrid <= res / 2)  # No frequency must be greater than the Nyquist frequency.
        length(bfreqGrid) >= K * L  # K*L is ideally going to be less than N, otherwise the bootstrap has no benefit in terms of compuation time.
        length(bfreqGrid) / (K * L) <= length(t_BLS) / 2  # This condition is mentioned in https://ui.adsabs.harvard.edu/abs/2012ada..confE..16S.
    })

    print(sprintf("Max frequency: %f, Min frequency: %f", max(bfreqGrid), min(bfreqGrid)))

    tfreqGrid <- getFreqGridToTest(t, period, duration, res=res, ofac=ofac, algo=algo, lctype="real")
    if (any(is.na(bfreqGrid))) {
        stop("Atleast one frequency in the frequency grid is NaN!")
    }

    stopifnot(exprs={
        all(tfreqGrid <= res / 2)  # No frequency must be greater than the Nyquist frequency.
        length(tfreqGrid) >= K * L  # K*L is ideally going to be less than N, otherwise the bootstrap has no benefit in terms of compuation time.
        length(tfreqGrid) / (K * L) <= length(t) / 2  # This condition is mentioned in https://ui.adsabs.harvard.edu/abs/2012ada..confE..16S.
    })

    print(sprintf("Max frequency: %f, Min frequency: %f", max(tfreqGrid), min(tfreqGrid)))

    # Get periodograms
    if (isTRUE(noiseType == 2) | applyGPRforBLS) {
        y <- getGPRResid(t, y)  # Run Gaussian Processes Regression on light curve if autoregressive noise is present.
    }
    boutput <- bls(y_BLS, t_BLS, bls.plot = FALSE, per.min=min(1/bfreqGrid), per.max=max(1/bfreqGrid), nper=length(bfreqGrid))

    fstep <- (max(tfreqGrid) - min(tfreqGrid)) / length(tfreqGrid)
    freqs <- seq(from = min(tfreqGrid), by = fstep, length.out = length(tfreqGrid))
    tperiodsToTry <- 1 / freqs
    # Empirical observation: In reality, applying ARMA when Gaussian noise is present will give another Gaussian, so is not that helpful.
    # However, for some reason, it was found to be important to keep ARMA irrespective of the noise to get expected results.
    # Hence we apply ARMA even if Gaussian noise is present.
    tresidTCF <- getResidForTCF(y)
    toutput <- tcf(tresidTCF, p.try = tperiodsToTry * res, print.output = TRUE)

    # Normalize, detrend, etc...
    # (1) Remove trend in periodogram
    # constraint='increase' is used, but may not be needed.
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

    blambdaScatter <- 1
    bcobsScatter <- cobs(boutput$periodsTested, bScatter, ic='BIC', tau=0.5, lambda=blambdaScatter)
    # cobss50 <- cobs(output$periodsTested, periodogramTrendRemoved, ic='BIC', tau=0.5, lambda=lambdaTrend, constraint="increase")  # If tau = 0.5 and lambda = 0 => Median regression fit.
    # cobss501 <- cobs(output$periodsTested, periodogramTrendRemoved, ic='BIC', tau=0.9, lambda=lambdaTrend, constraint="increase")
    # cobss502 <- cobs(output$periodsTested, periodogramTrendRemoved, ic='BIC', tau=0.99, lambda=lambdaTrend)
    tlambdaScatter <- 1
    tcobsScatter <- cobs(tperiodsToTry, tScatter, ic='BIC', tau=0.5, lambda=tlambdaScatter)

    bnormalizedPeriodogram <- bperiodogramTrendRemoved / bcobsScatter$fitted
    tnormalizedPeriodogram <- tperiodogramTrendRemoved / tcobsScatter$fitted

    # Plotting.
    png(filename="real_4.png", width = 420, height = 210, units='mm', res = 300)
    par(mar=c(5,6,4,2), cex=15)

    cexVal <- 2.0
    mat1 <- matrix(c(
        1, 3, 4,
        2, 3, 4,
        5, 6, 7,
        5, 6, 7
        ), nrow = 4, ncol = 3, byrow = TRUE
    )
    layout(matrix(c(1,1,1,1,2,2, 3,3,3,4,4,4, 5,5,5,6,6,6), nrow=3, ncol=6, byrow=TRUE))
    # layout(mat = mat1,
    #     heights = c(1, 1, 1, 1),    # Heights of the two rows
    #     widths = c(1.5, 2, 2)
    # )     # Widths of the two columns

    bpergram <- boutput$spec
    tpergram <- toutput$outpow

    plot(t/24, y, type='l', main="DTARPS 103 = TIC 89020549", cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xlab='Time (days)', ylab='Flux')
    # acfEstimate <- acf(y, plot = FALSE)
    # lJStats <- Box.test(y, lag = 1, type = "Ljung")  # We want to see autocorrelation with each lag, hence pass lag = 1.
    # n <- length(acfEstimate$acf)
    # plot(
    #     acfEstimate$acf[2:n], main=sprintf("P(Ljung-Box) = %.2f, ACF(1) = %.2f", lJStats[3], acfEstimate$acf[[2]]), cex=2, type="h", 
    #     xlab="Lag",     
    #     ylab="ACF", 
    #     ylim=c(-0.2,0.2), # this sets the y scale to -0.2 to 0.2
    #     las=1,
    #     xaxt="n"
    # )
    # abline(h=0)
    # # Add labels to the x-axis
    # x <- c(1:n)
    # y <- c(1:n)
    # axis(1, at=x, labels=y)

    # Note: We show the ACF of the time series as it is. In reality, we slightly modify the time series for BLS (see above).
    acfEstimate <- acf(y, plot = FALSE, na.action = na.pass)
    lJStats <- Box.test(y, lag = 1, type = "Ljung")  # We want to see autocorrelation with each lag, hence pass lag = 1.
    plot(acfEstimate, main="", cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xlim=c(1, 20), ylim=c(-0.2, +0.5))
    text(10, 0.36, sprintf("P(Ljung-Box) = %.2f, ACF(1) = %.2f\n", lJStats[3], acfEstimate$acf[[2]]), cex=1.9)
    # plot(acfEstimate, main=sprintf("P(Ljung-Box) = %.2f, ACF(1) = %.2f", lJStats[3], acfEstimate$acf[[2]]), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, cex=cexVal, ylim=(-0.2, +0.5))

    plot(bcobsxy50$x/24, bpergram, type = 'l', main="BLS periodogram", log='x', xlab='Period (days) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
    lines(bcobsxy50$x/24, bcobsxy50$fitted, type = 'l', col='red', lwd=3.0)
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug(bcobsxy50$knots/24)
    # legend("topleft", lty = 1, 
    #     col= c("red"), text.col = "black", 
    #     legend=c("trend fit"), bty="n", cex=1.5, pt.cex = 1
    # )
    text(2/24, 1.62e-4, paste0(sprintf("SNR = %.1f, FAP = %.1e\nPeriod = %.5f days, Depth = %.3f", calculateSNR(boutput$periodsTested, bpergram), fapBLS, perResultsBLS[1]/24, perResultsBLS[2]), "%"), cex=1.8, adj=0)
    plot(tcobsxy50$x/24, tpergram, type = 'l', main="TCF periodogram", log='x', xlab='Period (days) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
    lines(tcobsxy50$x/24, tcobsxy50$fitted, type = 'l', col='red', lwd=3.0)
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug(tcobsxy50$knots/24)
    # legend("topleft", lty = 1, 
    #     col= c("red"), text.col = "black", 
    #     legend=c("trend fit"), bty="n", cex=1.5, pt.cex = 1
    # )
    text(2/24, 107, paste0(sprintf("SNR = %.1f, FAP = %.1e\nPeriod = %.5f days, Depth = %.3f", calculateSNR(tperiodsToTry * res, tpergram), fapTCF, perResultsTCF[1]/24, perResultsTCF[2]), "%"), cex=1.8, adj=0)

    print(calculateSNR(tperiodsToTry * res, tpergram))
    print(calculateSNR(boutput$periodsTested, bpergram))

    # Plot histogram of detrended periodogram. Shows log-frequency on y-axis in histogram for better visualization.
    bhist.data = hist(bperiodogramTrendRemoved, breaks=50, plot = FALSE)
    # Compute skewness and kurtosis of the original and standardized histograms.
    ### Refer https://brownmath.com/stat/shape.htm for more information ###
    ### Note: R does NOT compute the "excess kurtosis".
    # The kurtosis is calculated as follows:
    # ```
    # n <- length(x)
    # n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
    # ``` Taken from https://stackoverflow.com/a/21484052
    bSkewnessBefore <- skewness(bperiodogramTrendRemoved)
    bKurtosisBefore <- kurtosis(bperiodogramTrendRemoved)

    thist.data = hist(tperiodogramTrendRemoved, breaks=50, plot = FALSE)
    # Compute skewness and kurtosis of the original and standardized histograms.
    ### Refer https://brownmath.com/stat/shape.htm for more information ###
    ### Note: R does NOT compute the "excess kurtosis".
    # The kurtosis is calculated as follows:
    # ```
    # n <- length(x)
    # n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
    # ``` Taken from https://stackoverflow.com/a/21484052
    tSkewnessBefore <- skewness(tperiodogramTrendRemoved)

    plot(bhist.data$count, type='h', log='y', main=sprintf('BLS periodogram (detrended) histogram'), cex=cexVal, cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
    axis(1, at=1:length(bhist.data$mids), labels=sprintf(bhist.data$mids, fmt="%.1e"), cex.axis=cexVal)
    plot(thist.data$count, type='h', log='y', main=sprintf('TCF periodogram (detrended) histogram'), cex=cexVal, cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
    axis(1, at=1:length(thist.data$mids), labels=sprintf(thist.data$mids, fmt="%.1e"), cex.axis=cexVal)
    tKurtosisBefore <- kurtosis(tperiodogramTrendRemoved)

    dev.off()
}

# Coordinates
# real_1.png: 2/24, 4.5e-4; 2/24, 123
# real_2.png: 4/24, 4.15e-5; 2/24, 123
# real_3.png: 35/24, 1.4e-4; 35/24, 240
# real_4.png: 2/24, 1.62e-4; 2/24, 107
