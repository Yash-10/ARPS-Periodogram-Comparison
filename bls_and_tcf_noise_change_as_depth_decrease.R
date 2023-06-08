#################################################################################
# Aim: Given a transit depth, run the extreme value code and related procedures
# and then plot the results.

# Notes:
# 1. This script can only be used if simulation is the intention.
# 2. This script was used for making some of the plots in the paper.

#################################################################################

library(forecast)
library('cobs')
library('boot')
library('extRemes')
library('kernlab')
library('goftest')  # install.packages("goftest")

source('BLS/bls.R')
source('TCF3.0/intf_libtcf.R')
source('utils.R')
source('eva_periodogram.R')

helperGetFAP <- function(
    depth, noiseType=1, ntransits=10, res=2, ofac=2,
    gaussStd=1e-4, useOptimalFreqSampling=TRUE,
    L=300, R=300, K=2, seedValue=465, applyGPRforBLS=FALSE
) {
    set.seed(seedValue)

    period <- 2
    duration <- 2

    # Generate light curve using the parameters.
    yt <- getLightCurve(period, depth, duration, noiseType=noiseType, ntransits=ntransits, res=res, gaussStd=gaussStd, seedValue=seedValue)
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

    stopifnot(exprs={
        identical(tfreqGrid, bfreqGrid)
    })

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

    # Standardize periodograms
    nboutput <- standardizeAPeriodogram(
        boutput,
        periodsToTry=NULL,  # This argument is only needed when algo="TCF" and not needed for algo="BLS".
        algo="BLS",
        mode='detrend_normalize',  # Other option is 'detrend' in which case only detrending is performed, no normalization using scatter is performed.
        scatterWindowLength=length(boutput$periodsTested)/10
    )
    bcobsxy50 <- nboutput[[2]]
    ntoutput <- standardizeAPeriodogram(
        toutput,
        periodsToTry=tperiodsToTry,  # This argument is only needed when algo="TCF" and not needed for algo="BLS".
        algo="TCF",
        mode='detrend_normalize',  # Other option is 'detrend' in which case only detrending is performed, no normalization using scatter is performed.
        scatterWindowLength=length(tperiodsToTry)/10
    )
    tcobsxy50 <- ntoutput[[2]]

    # Call extreme value analysis code.
    resultBLS <- evd(period, depth, duration, noiseType=noiseType, algo='BLS', ofac=ofac, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, FAPSNR_mode=0, seedValue=seedValue, useStandardization=TRUE, mode='detrend_normalize')
    resultTCF <- evd(period, depth, duration, noiseType=noiseType, algo='TCF', ofac=ofac, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, FAPSNR_mode=0, seedValue=seedValue, useStandardization=TRUE, mode='detrend_normalize')
    fapBLS <- resultBLS[1]
    fapTCF <- resultTCF[1]
    perResultsBLS <- resultBLS[2:4]
    perResultsTCF <- resultTCF[2:4]

    return (list(perResultsBLS, perResultsTCF, fapBLS, fapTCF))
}

helperGetFAPMultipleSeeds <- function(
    ...
) {
    set.seed(42)
    seeds <- sample(1:1e5, 50)

    faps_BLS <- c()
    faps_TCF <- c()

    for (seed in seeds) {
        x <- helperGetFAP(..., seedValue=seed)
        fapBLS <- x[[3]]
        fapTCF <- x[[4]]
        faps_BLS <- c(faps_BLS, fapBLS)
        faps_TCF <- c(faps_TCF, fapTCF)
    }
    return (c(median(faps_BLS), median(faps_TCF)))
}

blsAndTCFDepthChange <- function(
    depth, noiseType=1, ntransits=10, res=2, ofac=2,
    gaussStd=1e-4, useOptimalFreqSampling=TRUE,
    L=300, R=300, K=2, seedValue=465, applyGPRforBLS=FALSE,
    fapBLS=NaN, fapTCF=NaN
) {
    set.seed(seedValue)

    period <- 2
    duration <- 2

    # Generate light curve using the parameters.
    yt <- getLightCurve(period, depth, duration, noiseType=noiseType, ntransits=ntransits, res=res, gaussStd=gaussStd, seedValue=seedValue)
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

    stopifnot(exprs={
        identical(tfreqGrid, bfreqGrid)
    })

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

    # Standardize periodograms
    nboutput <- standardizeAPeriodogram(
        boutput,
        periodsToTry=NULL,  # This argument is only needed when algo="TCF" and not needed for algo="BLS".
        algo="BLS",
        mode='detrend_normalize',  # Other option is 'detrend' in which case only detrending is performed, no normalization using scatter is performed.
        scatterWindowLength=length(boutput$periodsTested)/10
    )
    bcobsxy50 <- nboutput[[2]]
    ntoutput <- standardizeAPeriodogram(
        toutput,
        periodsToTry=tperiodsToTry,  # This argument is only needed when algo="TCF" and not needed for algo="BLS".
        algo="TCF",
        mode='detrend_normalize',  # Other option is 'detrend' in which case only detrending is performed, no normalization using scatter is performed.
        scatterWindowLength=length(tperiodsToTry)/10
    )
    tcobsxy50 <- ntoutput[[2]]

    # Call extreme value analysis code.
    resultBLS <- evd(period, depth, duration, noiseType=noiseType, algo='BLS', ofac=ofac, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, FAPSNR_mode=0, seedValue=seedValue, useStandardization=TRUE, mode='detrend_normalize')
    resultTCF <- evd(period, depth, duration, noiseType=noiseType, algo='TCF', ofac=ofac, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, FAPSNR_mode=0, seedValue=seedValue, useStandardization=TRUE, mode='detrend_normalize')
    fapBLS <- resultBLS[1]
    fapTCF <- resultTCF[1]
    perResultsBLS <- resultBLS[2:4]
    perResultsTCF <- resultTCF[2:4]

    png(filename="depth_change_gaussian_0.02.png", width = 430, height = 320, units='mm', res = 300)
    par(mar=c(5,6,4,2), cex=15)

    cexVal <- 2.0
    layout(matrix(c(1,1,1,1,2,2, 3,3,3,4,4,4, 5,5,5,6,6,6, 7,7,7,8,8,8, 9,9,9,10,10,10), nrow=5, ncol=6, byrow=TRUE))
    # layout(mat = mat1,
    #     heights = c(1, 1, 1, 1),    # Heights of the two rows
    #     widths = c(1.5, 2, 2)
    # )     # Widths of the two columns

    bpergram <- boutput$spec
    tpergram <- toutput$outpow

    plot(t/24, y, type='l', main=sprintf("Period: %.1f days, Depth: %.1f ppm, Duration: %.1f hrs | Noise: Gaussian", period, depth*1e4, duration), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xlab='Time (days)', ylab='Flux')

    acfEstimate <- acf(y, plot = FALSE, na.action = na.pass)
    lJStats <- Box.test(y, lag = 1, type = "Ljung")  # We want to see autocorrelation with each lag, hence pass lag = 1.
    plot(acfEstimate, main="", cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xlim=c(1, 20), ylim=c(-0.2, +0.5))
    text(10, 0.36, sprintf("P(Ljung-Box) = %.2f, ACF(1) = %.2f\n", lJStats[3], acfEstimate$acf[[2]]), cex=1.9)
    # plot(acfEstimate, main=sprintf("P(Ljung-Box) = %.2f, ACF(1) = %.2f", lJStats[3], acfEstimate$acf[[2]]), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, cex=cexVal, xlim=c(1, 20), ylim=(-0.2, +0.5))

    # ROW 2 #######################################################################################################
    plot((10**bcobsxy50$x)/24, bpergram, type = 'l', main="BLS periodogram", log='x', xlab='Period (days) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
    lines((10**bcobsxy50$x)/24, bcobsxy50$fitted, type = 'l', col='red', lwd=3.0)
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug((10**bcobsxy50$knots)/24)
    # legend("topleft", lty = 1, 
    #     col= c("red"), text.col = "black", 
    #     legend=c("trend fit"), bty="n", cex=1.5, pt.cex = 1
    # )
    text(0.08, 3.8e-5, paste0(sprintf("Period = %.5f days, Depth = %.1f ppm\nSNR = %.1f", perResultsBLS[1]/24, perResultsBLS[2]*1e6, calculateSNR(boutput$periodsTested, bpergram))), cex=1.9, adj=0)
    plot(10**(tcobsxy50$x)/24, tpergram, type = 'l', main="TCF periodogram", log='x', xlab='Period (days) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
    lines((10**tcobsxy50$x)/24, tcobsxy50$fitted, type = 'l', col='red', lwd=3.0)
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug((10**tcobsxy50$knots)/24)
    # legend("topleft", lty = 1, 
    #     col= c("red"), text.col = "black", 
    #     legend=c("trend fit"), bty="n", cex=1.5, pt.cex = 1
    # )
    text(0.08, 153, paste0(sprintf("Period = %.5f days, Depth = %.1f ppm\nSNR = %.1f", perResultsTCF[1]/24, perResultsTCF[2]*1e6, calculateSNR(tperiodsToTry * res, tpergram))), cex=1.9, adj=0)
    ###############################################################################################################

    # print(calculateSNR(tperiodsToTry * res, tpergram))
    # print(calculateSNR(boutput$periodsTested, bpergram))

    # Row 3 ##############################################################################################################
    # Plot histogram of original periodograms. Shows log-frequency on y-axis in histogram for better visualization.
    bhist.data = hist(boutput$spec, breaks=50, plot = FALSE)
    # bSkewnessBefore <- skewness(boutput$spec)
    # bKurtosisBefore <- kurtosis(boutput$spec)

    thist.data = hist(toutput$outpow, breaks=50, plot = FALSE)
    # Compute skewness and kurtosis of the original and standardized histograms.
    ### Refer https://brownmath.com/stat/shape.htm for more information ###
    ### Note: R does NOT compute the "excess kurtosis".
    # The kurtosis is calculated as follows:
    # ```
    # n <- length(x)
    # n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
    # ``` Taken from https://stackoverflow.com/a/21484052
    # tSkewnessBefore <- skewness(toutput$outpow)

    plot(bhist.data$count, type='h', log='y', main=sprintf('BLS periodogram histogram'), cex=cexVal, cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
    axis(1, at=1:length(bhist.data$mids), labels=sprintf(bhist.data$mids, fmt="%.1e"), cex.axis=cexVal)
    plot(thist.data$count, type='h', log='y', main=sprintf('TCF periodogram histogram'), cex=cexVal, cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
    axis(1, at=1:length(thist.data$mids), labels=sprintf(thist.data$mids, fmt="%.1e"), cex.axis=cexVal)
    # tKurtosisBefore <- kurtosis(toutput$outpow)
    ###############################################################################################################

    # Row 4 ##############################################################################################################
    plot((10**bcobsxy50$x)/24, nboutput[[1]], type = 'l', main="Standardized BLS periodogram", log='x', xlab='Period (days) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug((10**bcobsxy50$knots)/24)
    # legend("topleft", lty = 1, 
    #     col= c("red"), text.col = "black", 
    #     legend=c("trend fit"), bty="n", cex=1.5, pt.cex = 1
    # )
    text(0.08, 15, paste0(sprintf("FAP = %s", formatC(fapBLS, format = "e", digits = 0))), cex=1.9, adj=0)
    plot(10**(tcobsxy50$x)/24, ntoutput[[1]], type = 'l', main="Standardized TCF periodogram", log='x', xlab='Period (days) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug((10**tcobsxy50$knots)/24)
    # legend("topleft", lty = 1, 
    #     col= c("red"), text.col = "black", 
    #     legend=c("trend fit"), bty="n", cex=1.5, pt.cex = 1
    # )
    text(0.08, 30, paste0(sprintf("FAP = %s", formatC(fapTCF, format = "e", digits = 0))), cex=1.9, adj=0)
    ###############################################################################################################

    # Plot histogram of standardized periodograms. Shows log-frequency on y-axis in histogram for better visualization.
    bhist.data = hist(nboutput[[1]], breaks=50, plot = FALSE)
    # bSkewnessBefore <- skewness(boutput$spec)
    # bKurtosisBefore <- kurtosis(boutput$spec)

    thist.data = hist(ntoutput[[1]], breaks=50, plot = FALSE)
    # Compute skewness and kurtosis of the original and standardized histograms.
    ### Refer https://brownmath.com/stat/shape.htm for more information ###
    ### Note: R does NOT compute the "excess kurtosis".
    # The kurtosis is calculated as follows:
    # ```
    # n <- length(x)
    # n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
    # ``` Taken from https://stackoverflow.com/a/21484052
    # tSkewnessBefore <- skewness(toutput$outpow)

    plot(bhist.data$count, type='h', log='y', main=sprintf('Standardized BLS periodogram histogram'), cex=cexVal, cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
    axis(1, at=1:length(bhist.data$mids), labels=sprintf(bhist.data$mids, fmt="%.1e"), cex.axis=cexVal)
    plot(thist.data$count, type='h', log='y', main=sprintf('Standardized TCF periodogram histogram'), cex=cexVal, cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
    axis(1, at=1:length(thist.data$mids), labels=sprintf(thist.data$mids, fmt="%.1e"), cex.axis=cexVal)
    # tKurtosisBefore <- kurtosis(toutput$outpow)

    # dev.print(png, 'depth_change_gaussian_0.008.png')
    dev.off()
}


###### 1. MDD for Gaussian case (465 seed value, period=2days, dur=2hrs):
# 0.006836735% for BLS and 0.009591837 for TCF
# MDD for Autoregressive case (544 seed value, period=2days, dur=2hrs):
# 0.02336735 for BLS and 0.02612245 for TCF

###### Text coordinates used:
# depth_change_gaussian_0.006.png (BLS, TCF): [adj=1]
    # For period, etc: (22, 1.67e-5), (22, 43)
    # For FAP: (22, 6), (22, 9.5)
# depth_change_gaussian_0.02.png: [adj=0]
    # For period, etc: (0.08, 3.8e-5), (0.08, 153)
    # For FAP: (0.08, 15), (0.08, 30)

# depth_change_ar_0.023.png: (BLS, TCF): (0.08, 4.6e-5), (0.08, 86); For FAP: 
# depth_change_ar_0.04.png: (BLS, TCF):
    # For period, etc: (0.08, 7.2e-5), (0.08, 210) -- both adj=0
    # For FAP: (22, 9), (22, 39)

# Data
# 0.02%         FAP_BLS = , FAP_TCF = 
# 0.006836735%: FAP_BLS = , FAP_TCF = 



########### Showing results on real light curves ###########
blsAndTCFDepthChangeReal <- function(
    table, noiseType=1, ntransits=10, res=2, ofac=2, useOptimalFreqSampling=TRUE,
    L=300, R=300, applyGPRforBLS=TRUE
) {
    period <- depth <- duration <- noiseType <- ntransits <- ar <- ma <- order <- gaussStd <- NULL
    significanceMode <- 'max'  # Since for real light curves, passing `expected_peak` is not possible.
    res <- 2
    K <- ofac

    rtBLS <- c()
    rtTCF <- c()

    table_BLS <- table[!is.na(table$Flux),]

    resultBLS <- evd(y=table_BLS$Flux, t=table_BLS$times, algo="BLS", FAPSNR_mode=0, lctype="real", applyGPRforBLS=TRUE, noiseType=noiseType, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, seedValue=seedValue, useStandardization=TRUE, mode='detrend_normalize')
    fapBLS <- resultBLS[1]
    perResultsBLS <- resultBLS[2:4]
    resultBLS <- evd(y=table_BLS$Flux, t=table_BLS$times, algo="BLS", FAPSNR_mode=1, lctype="real", applyGPRforBLS=TRUE, noiseType=noiseType, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, seedValue=seedValue, useStandardization=TRUE, mode='detrend_normalize')
    snrBLS <- 1 / resultBLS[1]
    rtBLS <- c(rtBLS, c(fapBLS, snrBLS, perResultsBLS))

    resultTCF <- evd(y=table$Flux, t=table$times, algo="TCF", FAPSNR_mode=0, lctype="real", noiseType=noiseType, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, seedValue=seedValue, useStandardization=TRUE, mode='detrend_normalize')
    fapTCF <- resultTCF[1]
    perResultsTCF <- resultTCF[2:4]
    resultTCF <- evd(y=table$Flux, t=table$times, algo="TCF", FAPSNR_mode=1, lctype="real", noiseType=noiseType, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, seedValue=seedValue, useStandardization=TRUE, mode='detrend_normalize')
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

    # Note that if NaN's are present in the light curve, then the BLS and TCF frequencies will not be the same.
    # Since all the four real TESS curves have missing values, BLS and TCF frequencies will be different.
    # Hence, it doesn't make sense to to check the below condition.
    # stopifnot(exprs={
    #     identical(tfreqGrid, bfreqGrid)
    # })

    # Normalized periodograms
    # Standardize periodograms
    nboutput <- standardizeAPeriodogram(
        boutput,
        periodsToTry=NULL,  # This argument is only needed when algo="TCF" and not needed for algo="BLS".
        algo="BLS",
        mode='detrend_normalize',  # Other option is 'detrend' in which case only detrending is performed, no normalization using scatter is performed.
        scatterWindowLength=length(boutput$periodsTested)/10
    )
    bcobsxy50 <- nboutput[[2]]
    ntoutput <- standardizeAPeriodogram(
        toutput,
        periodsToTry=tperiodsToTry,  # This argument is only needed when algo="TCF" and not needed for algo="BLS".
        algo="TCF",
        mode='detrend_normalize',  # Other option is 'detrend' in which case only detrending is performed, no normalization using scatter is performed.
        scatterWindowLength=length(tperiodsToTry)/10
    )
    tcobsxy50 <- ntoutput[[2]]

    # Plotting.
    png(filename="real_4.png", width = 430, height = 320, units='mm', res = 300)
    par(mar=c(5,6,4,2), cex=15)

    cexVal <- 2.0
    layout(matrix(c(1,1,1,1,2,2, 3,3,3,4,4,4, 5,5,5,6,6,6, 7,7,7,8,8,8, 9,9,9,10,10,10), nrow=5, ncol=6, byrow=TRUE))

    # layout(mat = mat1,
    #     heights = c(1, 1, 1, 1),    # Heights of the two rows
    #     widths = c(1.5, 2, 2)
    # )     # Widths of the two columns

    bpergram <- boutput$spec
    tpergram <- toutput$outpow

    plot(t/24, y, type='l', main=sprintf("DTARPS 103 = TIC 89020549"), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xlab='Time (days)', ylab='Flux')

    acfEstimate <- acf(y, plot = FALSE, na.action = na.pass)
    lJStats <- Box.test(y, lag = 1, type = "Ljung")  # We want to see autocorrelation with each lag, hence pass lag = 1.
    plot(acfEstimate, main="", cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xlim=c(1, 20), ylim=c(-0.2, +0.5))
    text(10, 0.36, sprintf("P(Ljung-Box) = %.2f, ACF(1) = %.2f\n", lJStats[3], acfEstimate$acf[[2]]), cex=1.9)
    # plot(acfEstimate, main=sprintf("P(Ljung-Box) = %.2f, ACF(1) = %.2f", lJStats[3], acfEstimate$acf[[2]]), cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, cex=cexVal, xlim=c(1, 20), ylim=(-0.2, +0.5))

    # ROW 2 #######################################################################################################
    plot((10**bcobsxy50$x)/24, bpergram, type = 'l', main="BLS periodogram", log='x', xlab='Period (days) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
    lines((10**bcobsxy50$x)/24, bcobsxy50$fitted, type = 'l', col='red', lwd=3.0)
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug((10**bcobsxy50$knots)/24)
    # legend("topleft", lty = 1, 
    #     col= c("red"), text.col = "black", 
    #     legend=c("trend fit"), bty="n", cex=1.5, pt.cex = 1
    # )
    text(2/24, 1.62e-4, paste0(sprintf("Period = %.5f days, Depth = %.1f\nSNR = %.1f", perResultsBLS[1]/24, perResultsBLS[2]*1e6, calculateSNR(boutput$periodsTested, bpergram))), cex=1.9, adj=0)
    plot(10**(tcobsxy50$x)/24, tpergram, type = 'l', main="TCF periodogram", log='x', xlab='Period (days) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
    lines((10**tcobsxy50$x)/24, tcobsxy50$fitted, type = 'l', col='red', lwd=3.0)
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug((10**tcobsxy50$knots)/24)
    # legend("topleft", lty = 1, 
    #     col= c("red"), text.col = "black", 
    #     legend=c("trend fit"), bty="n", cex=1.5, pt.cex = 1
    # )
    text(2/24, 107, paste0(sprintf("Period = %.5f days, Depth = %.1f\nSNR = %.1f", perResultsTCF[1]/24, perResultsTCF[2]*1e6, calculateSNR(tperiodsToTry * res, tpergram))), cex=1.9, adj=0)
    ###############################################################################################################

    # print(calculateSNR(tperiodsToTry * res, tpergram))
    # print(calculateSNR(boutput$periodsTested, bpergram))

    # Row 3 ##############################################################################################################
    # Plot histogram of original periodograms. Shows log-frequency on y-axis in histogram for better visualization.
    bhist.data = hist(boutput$spec, breaks=50, plot = FALSE)
    # bSkewnessBefore <- skewness(boutput$spec)
    # bKurtosisBefore <- kurtosis(boutput$spec)

    thist.data = hist(toutput$outpow, breaks=50, plot = FALSE)
    # Compute skewness and kurtosis of the original and standardized histograms.
    ### Refer https://brownmath.com/stat/shape.htm for more information ###
    ### Note: R does NOT compute the "excess kurtosis".
    # The kurtosis is calculated as follows:
    # ```
    # n <- length(x)
    # n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
    # ``` Taken from https://stackoverflow.com/a/21484052
    # tSkewnessBefore <- skewness(toutput$outpow)

    plot(bhist.data$count, type='h', log='y', main=sprintf('BLS periodogram histogram'), cex=cexVal, cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
    axis(1, at=1:length(bhist.data$mids), labels=sprintf(bhist.data$mids, fmt="%.1e"), cex.axis=cexVal)
    plot(thist.data$count, type='h', log='y', main=sprintf('TCF periodogram histogram'), cex=cexVal, cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
    axis(1, at=1:length(thist.data$mids), labels=sprintf(thist.data$mids, fmt="%.1e"), cex.axis=cexVal)
    # tKurtosisBefore <- kurtosis(toutput$outpow)
    ###############################################################################################################

    # Row 4 ##############################################################################################################
    plot((10**bcobsxy50$x)/24, nboutput[[1]], type = 'l', main="Standardized BLS periodogram", log='x', xlab='Period (days) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug((10**bcobsxy50$knots)/24)
    # legend("topleft", lty = 1, 
    #     col= c("red"), text.col = "black", 
    #     legend=c("trend fit"), bty="n", cex=1.5, pt.cex = 1
    # )
    text(22, 16, paste0(sprintf("FAP = %s", formatC(fapBLS, format = "e", digits = 0))), cex=1.9, adj=1)
    plot(10**(tcobsxy50$x)/24, ntoutput[[1]], type = 'l', main="Standardized TCF periodogram", log='x', xlab='Period (days) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug((10**tcobsxy50$knots)/24)
    # legend("topleft", lty = 1, 
    #     col= c("red"), text.col = "black", 
    #     legend=c("trend fit"), bty="n", cex=1.5, pt.cex = 1
    # )
    text(22, 19.5, paste0(sprintf("FAP = %s", formatC(fapTCF, format = "e", digits = 0))), cex=1.9, adj=1)
    ###############################################################################################################

    # Plot histogram of standardized periodograms. Shows log-frequency on y-axis in histogram for better visualization.
    bhist.data = hist(nboutput[[1]], breaks=50, plot = FALSE)
    # bSkewnessBefore <- skewness(boutput$spec)
    # bKurtosisBefore <- kurtosis(boutput$spec)

    thist.data = hist(ntoutput[[1]], breaks=50, plot = FALSE)
    # Compute skewness and kurtosis of the original and standardized histograms.
    ### Refer https://brownmath.com/stat/shape.htm for more information ###
    ### Note: R does NOT compute the "excess kurtosis".
    # The kurtosis is calculated as follows:
    # ```
    # n <- length(x)
    # n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
    # ``` Taken from https://stackoverflow.com/a/21484052
    # tSkewnessBefore <- skewness(toutput$outpow)

    plot(bhist.data$count, type='h', log='y', main=sprintf('Standardized BLS periodogram histogram'), cex=cexVal, cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
    axis(1, at=1:length(bhist.data$mids), labels=sprintf(bhist.data$mids, fmt="%.1e"), cex.axis=cexVal)
    plot(thist.data$count, type='h', log='y', main=sprintf('Standardized TCF periodogram histogram'), cex=cexVal, cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", lwd=10, lend=2, col='grey61', xlab='Power', ylab='Count')
    axis(1, at=1:length(thist.data$mids), labels=sprintf(thist.data$mids, fmt="%.1e"), cex.axis=cexVal)
    # tKurtosisBefore <- kurtosis(toutput$outpow)

    # dev.print(png, 'depth_change_gaussian_0.008.png')
    dev.off()
}

# Coordinates
# real_1.png: 2/24, 4.5e-4; 2/24, 123 ;; 22, 7; 22, 15.5
# real_2.png: 4/24, 4.15e-5; 2/24, 123 ;; 22, 14.5; 22, 29
# real_3.png: 35/24, 1.4e-4; 35/24, 240 ;; 
# real_4.png: 2/24, 1.62e-4; 2/24, 107 ;; 22, 16; 22, 19.5


# showFitOverlayed <- function(
#     table
# ) {
#     # TODO: For BLS, I think we first remove rows with Na flux values and only then run GPR, so try doing that here.
#     gp <- gausspr(table$times, table$Flux)
#     predicted_y <- predict(gp, table$times)

#     max.p = 5
#     max.q = 5
#     max.d = 0
#     ARIMA.fit = auto.arima(diff(table$Flux), stepwise=FALSE, approximation=FALSE, seasonal=FALSE, max.p=max.p, max.q=max.q, max.d=max.d, d=0)
#     # print(ARIMA.fit)

#     png(filename="gpr_fit.png", width = 420, height = 150, units='mm', res = 300)
#     par(mar=c(5,6,4,2), cex=15)

#     cexVal <- 1.7
#     # mat1 <- matrix(c(
#     #     1, 1, 1, 1, 2, 2,
#     #     3, 3, 3, 3, 4, 4
#     #     ), nrow = 2, ncol = 6, byrow = TRUE
#     # )
#     # layout(mat1)
#     layout(matrix(c(
#         1, 1, 1, 1, 2, 2,
#         3, 3, 3, 3, 4, 4
#     ), ncol = 6, nrow=2, byrow=TRUE))

#     ################# CODE FOR GPR ################
#     par(mar=c(0,6,4,2))
#     plot(table$times/24, table$Flux, col='black', type='l', main="", cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, ylab="flux", xlab="", xaxt="n")
#     lines(table$times/24, predicted_y, col='red', lwd=1.0)
#     text(24, 1.0012, "Gaussian Processes Regression fit", cex=cexVal, adj=1)

#     acfEstimate <- acf(table$Flux, plot = FALSE, na.action = na.pass)
#     lJStats <- Box.test(y, lag = 1, type = "Ljung")  # We want to see autocorrelation with each lag, hence pass lag = 1.
#     plot(acfEstimate, main="", cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", yaxt="n", xlim=c(1, 20), ylim=c(-0.2, +0.5))
#     text(15, 0.8, sprintf("P(Ljung-Box) = %.2f, ACF(1) = %.2f\n", lJStats[3], acfEstimate$acf[[2]]), cex=1.5)

#     axis(2, at=c(-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis=cexVal)

#     par(mar=c(5,6,0,2))
#     plot(table$times/24, table$Flux-predicted_y, col='black', type='l', main="", cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, ylab="Residual flux", xlab="Time (days)", yaxt="n")
#     text(24, 0.00135, "Gaussian Processes Regression residuals", cex=cexVal, adj=1)

#     axis(2, at=c(-0.001, 0.000, 0.001), cex.axis=cexVal)

#     acfEstimate <- acf(table$Flux - predicted_y, plot = FALSE, na.action = na.pass)
#     lJStats <- Box.test(y, lag = 1, type = "Ljung")  # We want to see autocorrelation with each lag, hence pass lag = 1.
#     plot(acfEstimate, main="", cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, yaxt="n", xlim=c(1, 20), ylim=c(-0.2, +0.5))
#     text(15, 0.8, sprintf("P(Ljung-Box) = %.2f, ACF(1) = %.2f\n", lJStats[3], acfEstimate$acf[[2]]), cex=1.5)
#     axis(2, at=c(-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis=cexVal)

#     ###############################################


#     ################ CODE FOR ARIMA ###############
#     # par(mar=c(0,6,4,2))
#     # plot(head(table$times/24, -1), diff(table$Flux), col='black', type='l', main="", cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, ylab="Differenced flux", xlab="", xaxt="n")
#     # lines(head(table$times/24, -1), fitted(ARIMA.fit), col='red', lwd=1.0)
#     # text(17, 0.0015, paste0("ARIMA(", ARIMA.fit$arma[[1]], ",1,", ARIMA.fit$arma[[2]], ") fit"), cex=cexVal, adj=1)

#     # acfEstimate <- acf(diff(table$Flux), plot = FALSE, na.action = na.pass)
#     # lJStats <- Box.test(y, lag = 1, type = "Ljung")  # We want to see autocorrelation with each lag, hence pass lag = 1.
#     # plot(acfEstimate, main="", cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, xaxt="n", yaxt="n", xlim=c(1, 20), ylim=c(-0.2, +0.5))
#     # text(15, 0.8, sprintf("P(Ljung-Box) = %.2f, ACF(1) = %.2f\n", lJStats[3], acfEstimate$acf[[2]]), cex=1.5)

#     # axis(2, at=c(-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis=cexVal)

#     # par(mar=c(5,6,0,2))
#     # plot(head(table$times/24, -1), diff(table$Flux)-fitted(ARIMA.fit), col='black', type='l', main="", cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, ylab="Residual flux", xlab="Time (days)", yaxt="n")
#     # text(17, 0.00125, "ARIMA residuals", cex=cexVal, adj=1)

#     # axis(2, at=c(-0.001, 0.000, 0.001), cex.axis=cexVal)

#     # acfEstimate <- acf(residuals(ARIMA.fit), plot = FALSE, na.action = na.pass)
#     # lJStats <- Box.test(y, lag = 1, type = "Ljung")  # We want to see autocorrelation with each lag, hence pass lag = 1.
#     # plot(acfEstimate, main="", cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal, yaxt="n", xlim=c(1, 20), ylim=c(-0.2, +0.5))
#     # text(15, 0.8, sprintf("P(Ljung-Box) = %.2f, ACF(1) = %.2f\n", lJStats[3], acfEstimate$acf[[2]]), cex=1.5)
#     # axis(2, at=c(-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis=cexVal)
#     ###############################################

#     dev.off()
# }

centerPieceFigure <- function(
    period=2, depth=0.0265, noiseType=2, duration=2, ntransits=10, res=2, ofac=2, useOptimalFreqSampling=FALSE, lctype="sim",
    applyGPRforBLS=TRUE, gaussStd=1e-4, seedValue=465, L=300, R=300
) {
    # Generate light curve using the parameters.
    yt <- getLightCurve(period, depth, duration, noiseType=noiseType, ntransits=ntransits, res=res, gaussStd=gaussStd, seedValue=seedValue)
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
    scatterWindowLength <- length(tfreqGrid) / 10
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
    resultBLS <- evd(period, depth, duration, noiseType=noiseType, algo='BLS', ofac=ofac, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, FAPSNR_mode=0, seedValue=seedValue)
    resultTCF <- evd(period, depth, duration, noiseType=noiseType, algo='TCF', ofac=ofac, L=L, R=R, res=res, ntransits=ntransits, gaussStd=gaussStd, FAPSNR_mode=0, seedValue=seedValue)
    fapBLS <- resultBLS[1]
    fapTCF <- resultTCF[1]

    bpergram <- boutput$spec
    tpergram <- toutput$outpow

    gp <- gausspr(t, y)
    predicted_y <- predict(gp, t)

    max.p = 5
    max.q = 5
    max.d = 0
    ARIMA.fit = auto.arima(diff(y), stepwise=FALSE, approximation=FALSE, seasonal=FALSE, max.p=max.p, max.q=max.q, max.d=max.d, d=0) #leave d as 0. 

    cexVal <-1.5
    layout(matrix(c(
        1, 1, 1, 2, 2, 2,
        3, 3, 3, 4, 4, 4
    ), ncol = 6, nrow=2, byrow=TRUE))

    plot(t, y, type='l')
    lines(t, predicted_y, col='red')

    plot(diff(y), col='black', type='l')
    lines(fitted(ARIMA.fit), col='red')

    plot(bcobsxy50$x/24, bpergram, type = 'l', main="BLS periodogram", log='x', xlab='Period (days) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
    lines(bcobsxy50$x/24, bcobsxy50$fitted, type = 'l', col='red', lwd=3.0)
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug(bcobsxy50$knots/24)
    # legend("topleft", lty = 1, 
    #     col= c("red"), text.col = "black", 
    #     legend=c("trend fit"), bty="n", cex=1.5, pt.cex = 1
    # )
    text(2/24, 1.62e-4, paste0(sprintf("SNR = %.1f, FAP = %.1e", calculateSNR(boutput$periodsTested, bpergram), fapBLS), "%"), cex=1.8, adj=0)
    plot(tcobsxy50$x/24, tpergram, type = 'l', main="TCF periodogram", log='x', xlab='Period (days) [log scale]', ylab='Power', cex.main=cexVal, cex.lab=cexVal, cex.axis=cexVal)
    lines(tcobsxy50$x/24, tcobsxy50$fitted, type = 'l', col='red', lwd=3.0)
    # lines(cobsxy501$x, cobsxy501$fitted, type = 'l', col='cyan')
    # lines(cobsxy502$x, cobsxy502$fitted, type = 'l', col='magenta')
    rug(tcobsxy50$knots/24)
    # legend("topleft", lty = 1, 
    #     col= c("red"), text.col = "black", 
    #     legend=c("trend fit"), bty="n", cex=1.5, pt.cex = 1
    # )
    text(2/24, 107, paste0(sprintf("SNR = %.1f, FAP = %.1e", calculateSNR(tperiodsToTry * res, tpergram), fapTCF), "%"), cex=1.8, adj=0)

    blsFalsePeakAInd <- which(rev(bpergram) > 3e-5)[7]
    blsFalsePeakBInd <- head(which(bpergram > 4e-5), n=1)

    x <- bcobsxy50$x/24
    print(sprintf("BLS false peak A at period = %f days with power = %f", rev(x)[blsFalsePeakAInd], rev(bpergram)[blsFalsePeakAInd]))
    print(sprintf("BLS false peak B at period = %f days with power = %f", x[blsFalsePeakBInd], bpergram[blsFalsePeakBInd]))
}