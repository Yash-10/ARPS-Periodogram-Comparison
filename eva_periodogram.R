# ***************************************************************************
# Author: Yash Gondhalekar  Last updated: March, 2023

# Description: This script contains the functions used for the extreme value
#              application. The `evd` function is the driver code.

#        Note:
#           1. Here, TCF periods are scaled by `res` since TCF calculates the
#              periodogram in units of cadence rather than absolute time values.
#           2. So this means the periods passed to TCF and BLS, in terms of
#              values, could be different, but they are using the same periods.
#           3. All of the period values printed on terminal are in hours,
#              unless explicitly mentioned anywhere.

# ***************************************************************************


library('extRemes')
library('boot')
library('cobs')
source('BLS/bls.R')
source('TCF3.0/intf_libtcf.R')
source('utils.R')
library('goftest')  # install.packages("goftest")
library('gbutils')  # https://search.r-project.org/CRAN/refmans/gbutils/html/cdf2quantile.html

# Resample a given vector.
boot_stat <- function(original_vector, resample_vector) {
    original_vector[resample_vector]
}

# Calculates the return level corresponding to a given FAP.
calculateReturnLevel <- function(
    fap,   # Requested FAP, e.g. 0.01.
    # Below three lines are the parameters of the fitted GEV model.
    location,
    scale,
    shape,
    K, L,  # These are parameters used for bootstrapping time series.
    n  # Length of the full frequency grid.
) {
    returnLevel <- qevd(
        1 - ((fap * K * L) / n),
        loc=location, scale=scale, shape=shape, type="GEV"
    )
    return (returnLevel);
}

# Calculates the FAP for the observed periodogram maxima.
calculateFAP <- function(
    location,
    scale,
    shape,
    K, L,
    n,  # Length of the full frequency grid.
    periodogramMaxima
) {
    # See equation 5 in https://academic.oup.com/mnras/article/450/2/2052/983840
    calculatedFAP <- 1 - pevd(periodogramMaxima, loc=location, scale=scale, shape=shape, type="GEV")
    # To manually calculate the FAP, use: 1 - exp(-(1 + shape * ((periodogramMaxima - location) / scale)) ^ (-1 / shape))
    return (calculatedFAP)
}

# The main extreme value function.
evd <- function(
    period=NULL,  # in days.
    depth=NULL,  # in %
    duration=NULL,  # in hours.
    y=NULL, t=NULL,  # If type == "real" (see below), then y (observations) and t (time) are needed.
    L=300,  # No. of distinct frequency bins.
    R=300,  # No. of bootstrap resamples of the original time series.
    noiseType=1,  # Noise model present in y. Either 1 (white gaussian noise) or 2 (autoregressive noise).
    # Note: noiseType is passed as argument to the `getLightCurve` function.
    useStandardization=FALSE,  # If true, uses standardized periodograms for GEV fitting and extrapolation.
    algo="BLS",  # Either BLS or TCF.
    ntransits=10,  # No. of transits to simulate.
    plot=TRUE,  # Whether to plot.
    ofac=2,  # ofac is also called as "samples per peak" sometimes.
    # It is recommended to use `useOptimalFreqSampling` for transit periodograms. For periodograms of sine-like signals, 1/s is the frequency resolution, so it is not necessary to use optimal frequency sampling and instead better use uniform frequency sampling.
    useOptimalFreqSampling=TRUE,  # If want to use the optimal frequency sampling from Ofir, 2014: delta_freq = q / (s * os), where s is whole time series duration, os is oversampling factor and q is the duty cycle (time in single transit / total time series duration).
    alpha=0.01,  # Significance level for hypothesis testing on the GEV fit on periodogram maxima. TODO: How to choose a significance level beforehand - any heuristics to follow?
    # Parameters for controlling noise. NOTE: `gaussStd` is considered for noiseType=2 as well for appropriate scaling of the autoregressive noise.
    # HENCE IT IS IMPORTANT TO KEEP THE SAME gaussStd value WHEN COMPARING BETWEEN AUTOREGRESSIVE AND GAUSSIAN NOISE CASES.
    gaussStd=1e-4,  # 0.01% Gaussian noise
    # Imp: ar, ma, and order arguments are now obsolete. They are not used anywhere.
    # ar=0.2,
    # ma=0.2,
    # order=c(1, 0, 1),
    res=2,  # Resolution for creating the time series. Refer getLightCurve from utils.R
    mode='detrend',  # Standardization mode: either detrend_normalize or detrend, see the function `standardizeAPeriodogram`. Only used if useStandardization=TRUE.
    checkConditions=TRUE,  # Mostly passed to light curve generation code. Also used in ad.test p-value check in this function.
    significanceMode='max',  # This tells whose significance (in terms of FAP) should be reported. 'max' means report significance of max(output), where output is the original periodogram (note output can be detrended/standardized based on `useStandardization` and `mode`).
    # (cont...) Other option is 'expected_peak' which tells to calculate significance of not the maximum power but the power corresponding to the expected period. 'expected_peak' option can only be used in simulations.
    seedValue=1,  # Seed to sue for reproducibility.
    FAPSNR_mode=0,  # 0 means only FAP, 1 means only SNR, and 2 means a linear combination of FAP and SNR.
    lctype="sim",  # Light curve type. Allowed values: sim or real. This parameter controls whether the light curve needs to be simulated or is a real light curve. In the former case, period, depth, and duration is needed at the least. In the latter case, y and t are needed as input.
    applyGPRforBLS=FALSE,  # This controls whether Gaussian Process Regression needs to be run before BLS.
    applyARMAforBLS=FALSE  # This controls whethter an ARMA model must be fit to the original light curve before applying BLS.
    # It is not recommended to use `applyARMAforBLS` since the ARMA model was shown to significantly model the transits as well.
    # This issue does not occur when ARMA is used on the differenced light curve before TCF since the no. of points of the transit in the differenced ligt curve is very small.
) {
    # TODO: Add comment if any assumption about Nan is assumed by this code.
    # Perform some checks.
    if (lctype == "sim" && (is.null(period) | is.null(depth) | is.null(duration))) {
        stop("type is set to `sim`, but at least one of {period, depth, or duration} is not specified!")
    }
    if (lctype == "real" && (is.null(y) | is.null(t))) {
        stop("type is set to `real`, but at least one of {y, t} is not specified!")
    }
    if (lctype == "real") {
        period <- depth <- duration <- noiseType <- ntransits <- ar <- ma <- order <- gaussStd <- NULL
        significanceMode <- 'max'  # Since for real light curves, passing `expected_peak` is not possible since the period is not known.
        res <- 2
    }

    # L, R, noiseType, ntransits must be integers.
    stopifnot(exprs={
        L %% 1 == 0
        R %% 1 == 0
        noiseType %% 1 == 0
        ntransits %% 1 == 0
    })

    # Oversampling factor must be an integer greater than or equal to zero.
    stopifnot(exprs={
        ofac %% 1 == 0  # Checks for integer.
        ofac >= 1
    })

    K <- ofac  # No. of distinct frequencies in a frequency bin.  # Note that in Suveges, 2014, K = 16 is used and K is called as the oversampling factor. So we also do that.
    # In short, L allows capturing long-range dependence while K prevents spectral leakage -- from Suveges.
    # FAP and SNR estimates should remain largely unaffected by L and R. That's what was shown in our paper.

    if (lctype == "sim") {
        # Generate light curve using the parameters.
        yt <- getLightCurve(period, depth, duration, noiseType=noiseType, ntransits=ntransits, gaussStd=gaussStd, res=res, checkConditions=checkConditions, seedValue=seedValue)
        y <- unlist(yt[1])
        t <- unlist(yt[2])
    }

    # Error handling.
    if (lctype == 'sim' & (any(is.na(y)) | any(is.na(t)))) {
        stop("Atleast one value in the observations or the time epochs is NaN! while lctype is sim!")
    }

    # Get the ACF estimate.
    acfEstimate <- acf(y, plot = FALSE, na.action = na.pass)
    print(sprintf("ACF at lag-1: %s", acfEstimate$acf[[2]]))

    # Special case (TCF fails if absolutely no noise -- so add a very small amount of noise just to prevent any errors).
    if (lctype == "sim") {
        if (noiseType == 0 && algo == "TCF") {
            y <- y + 10^-10 * rnorm(length(y))
        }
    }

    # (1) Bootstrap the time series.
    # Note that bootstrapping, by definition, is resampling "with replacement": https://en.wikipedia.org/wiki/Bootstrapping_(statistics)
    # Also note that Suveges says that the marginal distribution of each of the bootstrapped resample time series must approximately be the same as the original time series.

    # bootTS <- replicate(R, sample(y, length(y), replace=TRUE))
    # bootTS <- aperm(bootTS)  # This just permutes the dimension of bootTS - rows become columns and columns become rows - just done for easier indexing further in the code.
    set.seed(seedValue)
    bootTS <- boot(y, statistic=boot_stat, R=R)$t

    # NOTE: The below commented out lines of code can be used to test if the bootstrap resampling suggests the data is white noise or not.
    # However, it was observed in around 1 in 500 cases, the box test suggested the bootstrapped time series is not white noise. But we ignore that since 1/500 is a low probability.
    # for (j in 1:R) {
    #     lJStats <- Box.test(bootTS[j,], lag = 1, type = "Ljung")
    #     print(lJStats[3])
    #     stopifnot(exprs={
    #         lJStats[3] > 0.01
    #     })
    # }

    # Ensure dimensions of the bootstrapped set is as expected.
    # There must be R bootstrapped light curves, each with length = length(y).
    stopifnot(exprs={
        dim(bootTS) == c(R, length(y))
    })

    # Create a frequency grid.
    freqGrid <- getFreqGridToTest(t, period, duration, res=res, ofac=ofac, useOptimalFreqSampling=useOptimalFreqSampling, algo=algo, lctype=lctype)
    if (any(is.na(freqGrid))) {
        stop("Atleast one frequency in the frequency grid is NaN!")
    }

    # Some error checking.
    stopifnot(exprs={
        all(freqGrid <= res / 2)  # No frequency must be greater than the Nyquist frequency.
        length(freqGrid) >= K * L  # K*L is ideally going to be less than N, otherwise the bootstrap has no benefit in terms of compuation time.
        length(freqGrid) / (K * L) <= length(t) / 2  # This condition is mentioned in https://ui.adsabs.harvard.edu/abs/2012ada..confE..16S.
    })

    print(sprintf("Max frequency: %f, Min frequency: %f", max(freqGrid), min(freqGrid)))

    # Compute full periodogram (to be used afterwards when using FAP (mode=0 or 2), and terminate after this only if using mode=1).
    if (algo == "BLS") {
        if (isTRUE(noiseType == 2) | applyGPRforBLS) {
            y <- getGPRResid(t, y)  # Run Gaussian Processes Regression on light curve if autoregressive noise is present.
        }
        if (applyARMAforBLS) {
            y <- getARMAresid(y)  # Run ARMA before BLS is applied. This was shown to deteriorate performance and is NOT recommended.
        }
        output <- bls(y, t, bls.plot = FALSE, per.min=min(1/freqGrid), per.max=max(1/freqGrid), nper=length(freqGrid))
        ptested <- output$periodsTested

        # Calculate SNR of the periodogram peak.
        snr <- calculateSNR(ptested, output$spec, lambdaTrend=1, oneSideWindowLength=1500)

        scatterWindowLength <- length(ptested) / 10
        print(sprintf('Scatter window length for standardization used: %d', as.integer(scatterWindowLength)))
        perResults <- c(output$per, output$depth, output$dur)
        if (useStandardization) {
            output <- standardizeAPeriodogram(output, periodsToTry=NULL, algo="BLS", mode=mode, scatterWindowLength=scatterWindowLength)[[1]]
        }
        else {
            output <- output$spec
        }
        periodAtMaxOutput <- ptested[which.max(output)]
    }
    else if (algo == "TCF") {
        fstep <- (max(freqGrid) - min(freqGrid)) / length(freqGrid)
        freqs <- seq(from = min(freqGrid), by = fstep, length.out = length(freqGrid))
        periodsToTry <- 1 / freqs
        # Empirical observation: In reality, applying ARMA when Gaussian noise is present will give another Gaussian, so is not that helpful.
        # However, for some reason, it was found to be important to keep ARMA irrespective of the noise to get expected results.
        # Hence we apply ARMA even if Gaussian noise is present.
        tresidTCF <- getResidForTCF(y)
        output <- tcf(tresidTCF, p.try = periodsToTry * res, print.output = TRUE)

        snr <- calculateSNR(periodsToTry * res, output$outpow, lambdaTrend=1, oneSideWindowLength=1500)

        powmax.loc = which.max(output$outpow)
        perResults <- c(output$inper[powmax.loc]/res, output$outdepth[powmax.loc], output$outdur[powmax.loc]/res)
        if (useStandardization) {
            scatterWindowLength <- length(periodsToTry) / 10
            print(sprintf('Scatter window length for standardization used: %d', as.integer(scatterWindowLength)))
            output <- standardizeAPeriodogram(output, periodsToTry=periodsToTry, algo="TCF", mode=mode, scatterWindowLength=scatterWindowLength)[[1]]
        }
        else {
            output <- output$outpow
        }
        periodAtMaxOutput <- periodsToTry[which.max(output)]
    }
    if (snr < 0) {  # Ideally this will not occur because periodogram peaks will never be negative (unless some strong detrending has been applied) and IQR, by definition, cannot be negative.
        print('Negative SNR, returning NA score.')
        score <- NA
    }
    print(sprintf("Signal-to-noise ratio of periodogram peak = %f", snr))

    # (2) Max of each partial periodogram
    # Note that from the Suveges, 2014 paper, the reason for doing block maxima is: "The principal goal is to decrease the computational load due to a bootstrap. At the same time, the reduced frequency set should reflect the fundamental characteristics of a full periodogram: ..."
    maxima_R <- c()
    for (j in 1:R) {
        KLfreqs <- freqdivideFreqGrid(freqGrid, L, K, seedValue=seedValue)

        if (any(is.na(KLfreqs))) {
            stop("Atleast one frequency in the subset grid is NaN!")
        }

        stopifnot(exprs={
            length(KLfreqs) == K * L  # Obviously, length of the K*L frequencies must be K*L.
            length(bootTS[j,]) == length(y)  # Bootstrapped time series must be of same length as original time series.
            all(unique(KLfreqs) == KLfreqs)  # Note that the K*L frequencies must be non-overlapping. So all must be unique.
        })

        if (algo == "BLS") {
            out <- bls(bootTS[j,], t, per.min=min(1/KLfreqs), per.max=max(1/KLfreqs), nper=K*L, bls.plot = FALSE, print.output = FALSE)
            if (useStandardization) {
                ptestedPartial <- out$periodsTested
                # Compute scatter in windows only if the length of the periodogram on bootstrapped light curve is suficiently high.
                # By sufficiently high, we use 2000 as the limit.
                # The reason for doing this is because we have found the scatter estimate to be unreliable/noisy if the window size is too small.
                # Since these partial periodograms are much smaller than the original periodogram, we check this condition so that windows are used only when really needed.
                # *Important* : The length of partial periodograms = K * L, so it is entirely fixed by the parameters K and L. So unless those are changed, the partial periodograms would have the same length.
                if (length(ptestedPartial) > 2000) {
                    scatterWindowLength <- length(ptestedPartial) / 10
                }
                else {
                    scatterWindowLength <- length(ptestedPartial)
                }
                partialPeriodogram <- standardizeAPeriodogram(out, periodsToTry=NULL, algo="BLS", mode=mode, scatterWindowLength=scatterWindowLength)[[1]]  # For BLS, the periods tested is gotten from the R object `out`, hence we do not need to pass periodsToTy.
                # print(sprintf('Scatter window length for standardization used: %d', scatterWindowLength))
            }
            else {
                partialPeriodogram <- out$spec
            }
        }
        else if (algo == "TCF") {
            # Note: We do not need auto.arima here since the bootstrapped time series corresponds to white noise, and so ARIMA is of no use here.
            # Note: The reason for defining freqsPartial like this is to ensure that TCF and BLS use same set of frequencies for calculating periodogram.
            freqStepPartial <- (max(KLfreqs) - min(KLfreqs)) / (K * L)
            # freqsPartial are the same frequencies as used in BLS (verified).
            freqsPartial <- seq(from = min(KLfreqs), by = freqStepPartial, length.out = K*L)
            pToTry <- 1 / freqsPartial
            # If the original light curve contains NA values, then TCF on bootstrapped light curve can result in
            # TCF powers to blow up to very large numbers. Hence, we remove NA values only in this bootstrapped case.
            # NOTE: It is important not to remove NA in this way for the original TCF periodogram already computed above, but should be done for the bootstrapped ligth curve's periodogram for TCF..
            out <- tcf(diff(bootTS[j,][!is.na(bootTS[j,])]), p.try = pToTry * res, print.output = FALSE)  # Multiplying by res because TCF works according to cadence rather than actual time values, unlike BLS. Doing this ensures, TCF still peaks at 72 hr for different res values, for example.
            if (useStandardization) {
                ptestedPartial <- pToTry
                if (length(ptestedPartial) > 2000) {
                    scatterWindowLength <- length(ptestedPartial) / 10
                }
                else {
                    scatterWindowLength <- length(ptestedPartial)
                }
                partialPeriodogram <- standardizeAPeriodogram(out, periodsToTry = pToTry, algo="TCF", mode=mode, scatterWindowLength=scatterWindowLength)[[1]]
            }
            else {
                partialPeriodogram <- out$outpow
            }
        }

        # *** Some notes on declustering below -- DECLUSTERING IS NOT USED IN THIS CODE ***
        # Note: If we use oversampling, then while it increases the flexibility to choose frequencies in the frequency grid, it also has important issues as noted in https://academic.oup.com/mnras/article/388/4/1693/981666:
        # (1) "if we oversample the periodogram, the powers at the sampled frequencies are no longer independent..."
        # To solve the above problem, we decluster the partial periodograms. Even without oversampling, the peaks tend to be clustered and we need to decluster the peaks. One can see performance with and without declustering.
        # Decluster the peaks: https://search.r-project.org/CRAN/refmans/extRemes/html/decluster.html
        # Some intution on how to choose the threshold: https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.996.914&rep=rep1&type=pdf (search for threshold - Ctrl+F - in the paper)
        # See section 5.3.2 in Coles, 2001 to see why declustering is needed: Extremes tend to cluster themselves and tend to occur in groups. Note that log-likelihood can be decomposed into a product of individual marginal distribution functions only under iid. So declustering "tries" to make them independent to prevent the violation of the iid assumption while fitting the GEV model below.
        # In short, declustering (approximately) solves the dependence issue of extremes.
        # We might not need declustering in all cases -- we can calculate the extremel index and do declustering only if index < 1...
        # **NOTE**: Due to our way of extracting maxima of periodograms (i.e. not from whole periodogram but only from partial periodogram), maybe we do not even need declustering because the peaks, while fitting the GEV, do not have any physical sense, so declustering is not required.
        # How to choose the best threshold for declustering is an important aspect in this.
        # partialPeriodogram <- decluster(partialPeriodogram, threshold = quantile(partialPeriodogram, probs=c(0.75)))
        # ****************************************
        maxima_R <- append(maxima_R, max(partialPeriodogram))
    }

    print(sprintf("Maxima of the R maxima: %f", max(maxima_R)))
    print(sprintf("Maximum of the original periodogram: %f", max(output)))

    if (FAPSNR_mode == 1) {
        score <- 1 / snr
        return (c(score, perResults))
    }

    print("Done calculating maxima...")

    # (3) GEV modelling of partial periodograms' maxima
    fitEVD <- fevd(maxima_R, type='GEV')
    # See https://www.dataanalysisclassroom.com/lesson60/ for discussion on the fevd function.
    print(summary(fitEVD))
    distill(fitEVD)

    ## Get the fitted GEV parameters
    location <- findpars(fitEVD)$location[1]  # In extRemes, the parameter values repeat R times (for stationary models), and all are same. So extract the first.
    scale <- findpars(fitEVD)$scale[1]
    shape <- findpars(fitEVD)$shape[1]

    print(sprintf("location: %f, scale: %f, shape: %f", location, scale, shape))

    ## Important note: It would be better to find an automatic way to judge whether we want to select a GEV model or not, instead of manually looking at the diagnostic plots. This is because we want to apply this method on several periodograms. Hence we perform the A-D test.
    # Diagnostic goodness-of-fit tests (we use the Anderson-Darling (AD) test: https://search.r-project.org/CRAN/refmans/DescTools/html/AndersonDarlingTest.html)
    # A simple reason why we use the Anderson–Darling (AD) test rather than Komogorov-Smirnov (KS) is that AD is able to detect better the situations in which F0 and F differ on the tails (that is, for extreme data), where H0: F = F0 and H1: F \neq F0.
    # A NOTE: Take a look at https://academic.oup.com/mnras/article/449/1/1098/1314897: where it is written "[It should be admitted that a shortcut was taken in these tests: estimated parameters were treated as known...". It means that paper also used something similar to estimated = FALSE. See that paragraph in the paper to know why this shortcut is justified.
    result <- ad.test(maxima_R, null = "pevd", loc=location, scale=scale, shape=shape, nullname = "pevd", estimated = FALSE)  # estimated = TRUE would have been fine as well since the gevd parameters (location, scale, shape) are estimated using the data itself - those three parameters are not data-agnostic. But here we use estimated = FALSE because using TRUE uses a different variant of AD test using the Braun's method which we do not want.
    print(result)
    print(sprintf("p-value for Anderson-Darling goodness-of-fit test of the periodogram maxima: %f", result$p.value))
    # Check if AD fit is good enough. If not, returns NULL.
    # This check serves as a way to "automatically" find if the GEV fit is good and if it can be extrapolated to the full periodogram.
    # Suveges, 2014 suggests looking at the diagnostic plots before extrapolating to full periodogram, but that is cumbersome for large-scale simulations. Hence, this is a simple way to overcome manual fit quality inspection.
    if (checkConditions) {
        if (result$p.value < alpha) {  # Reject null hypothesis: the maxima sample is in the favor of alternate hypothesis (that the sample comes from a different distribution than GEV).
            warning("Anderson-Darling test failed while fitting GEV to the sample periodogram maxima. This means the GEV fit was sub-optimal and the FAP estimate may not be reliable.")
        }
    }

    # Diagnostic plots.
    if (plot) {
        # TODO: Why ci fails sometimes? See some discussion here: https://www.facebook.com/groups/254966138199366/posts/1167467316949239/
        try(plot(fitEVD))
        # try(plot(fitEVD, "trace"))
        # return.level(fitEVD)
        # return.level(fitEVD, do.ci = TRUE)
        # ci(fitEVD, return.period = c(2, 20, 100))
        # See some description on how ci's are calculated: https://reliability.readthedocs.io/en/latest/How%20are%20the%20confidence%20intervals%20calculated.html
    }

    # (4) Extrapolation to full periodogram
    # Note: The full periodogram was already computed towards the start of this function.
    print("Extrapolating to full periodogram...")

    print("Calculating return level...")
    returnLevel <- calculateReturnLevel(0.01, location, scale, shape, K, L, length(freqGrid))
    print(sprintf("Return level corresponding to FAP = %f: %f", 0.01, returnLevel))

    # For interpretation, we would like to get FAP given a return level rather than giving return level from a given FAP.
    if (significanceMode == 'max') {
        toCheck <- max(output)    
    }
    else if ((significanceMode == 'expected_peak') & (lctype == "sim")) {
        if (algo == "BLS") {
            x <- which.min(abs(ptested - period * 24))
            # The below is done in reality, there would be differences between estimated period and actual period, hence we need to select max from an interval around the expected period.
            indsLow <- x - 100
            indsUp <- x + 100
            toCheck <- max(output[indsLow:indsUp])
        }
        else if (algo == "TCF") {
            x <- which.min(abs(periodsToTry - period * 24))
            indsLow <- x - 100
            indsUp <- x + 100
            toCheck <- max(output[indsLow:indsUp])
        }
    }

    if (FAPSNR_mode == 0 || FAPSNR_mode == 2) {
        print("Calculating FAP...")
        fap <- calculateFAP(location, scale, shape, K, L, length(freqGrid), toCheck)
        print(sprintf("FAP = %.10f", fap))
    }

    if (FAPSNR_mode == 2) {
        score <- 0.75 * fap + 0.25 * (1 / snr)   # Note if 0 <= SNR < 1, then 1/snr will be > 1 and hence the score will be heavily penalized.
    }
    else if (FAPSNR_mode == 0) {
        score <- fap
    }
    print(sprintf("Overall score for this periodogram peak = %f", score))

    return (c(score, perResults))

    ###### Interpreting what FAP is good (from Baluev, 2008 paper: https://academic.oup.com/mnras/article/385/3/1279/1010111):
    # (1) > Given some small critical value FAP* (usually between 10−3 and 0.1), we can claim that the candidate signal is statistically
    # significant (if FAP < FAP*) or is not (if FAP > FAP*)
}

smallestPlanetDetectableTest <- function(  # This function returns the smallest planet detectable (in terms of transit depth) using the FAP criterion.
    period,  # a single period, in days.
    depths,  # vector of depths for which FAP needs to be calculated, each in %. An example: c(0.1, 0.08, 0.06, 0.04, 0.02, 0.015, 0.012, 0.01, 0.005)
    duration,  # a single duration, in hours.
    ...  # Arguments passed to evd() internally.
) {
    faps <- c()
    for (depth in depths) {
        result <- evd(period, depth, duration, ...)
        fap <- result[1]
        print(sprintf("depth (ppm): %f, fap: %f", depth*1e4, fap))
        faps <- append(faps, fap)
    }

    png(filename=sprintf("%sdays_%shours.png", period, duration))
    plot(depths*1e4, faps, xlab='Depth (ppm)', ylab='FAP', type='o', ylim=c(1e-7, 0.02), log='y')  # Upper limit is set to 0.02 which is slightly larger than 0.01, the threshold FAP.
    axis(1, at=1:length(depths), labels=depths*1e4)
    abline(h=0.01, col='black', lty=2)  # Here 1% FAP is used. Another choice is to use FAP=0.003, which corresponds to 3-sigma criterion for Gaussian -- commonly used in astronomy.
    dev.off()
}

# This function finds the root of the equation: FAP(depth, **params) - 0.01 = 0, i.e., given the period and duration of a planet,
# it finds the depth corresponding to the case FAP = 0.01 called the limiting_depth. So any transit with depth < limiting_depth
# is statistically insignificant using the FAP = 0.01 criterion.
depthEquation <- function(depth, period, duration, ...) {
    result <- try(evd(period, depth, duration, ...))
    # The below hack is only performed to prevent errors in long-running.
    result <- if (inherits(result, "try-error")) 1.0 else result 
    return (result[1] - 0.01);
}

# This function is a high-level wrapper for `findLimitingDepth` that prints the limiting depth.
# Root solving is done using the Newton-Raphson iteration method via the `uniroot` function in R.
findLimitingDepth <- function(period, duration, ...) {
    print('Finding limiting depth corresponding to FAP = 0.01, the fixed threshold FAP...')
    de <- function(depth) { return (depthEquation(depth, period=period, duration=duration, ...)); }
    # TODO: In future, we can set the upper limit of interval depth (currently, 0.3) intelligently based on the IQR of noise, for example.
    # Note that extendInt = "yes" is passed so that limiting depths for already significant planets (that extend beyond the passed interval) can also be availed using this function. So extendInt is added only for allowing applications to various types of planets.
    # Below try() hack from https://stackoverflow.com/questions/70200174/how-can-i-avoid-uniroot-error-that-stops-the-loop
    findLDepth <- try(uniroot(de, interval=c(0., 1.), tol=1e-4, extendInt="yes"))
    return (if (inherits(findLDepth, "try-error")) NA else findLDepth$root)
}

# This function is only for a quick verification test. One would not expect to get the exact depth where the planet starts to become insignificant.
periodDurationDepthTest <- function(
    algo="BLS",
    depths=c(0.1, 0.08, 0.06, 0.04, 0.02, 0.015, 0.012, 0.01, 0.005),  # in %
    ofac=1
) {
    periodDurations <- list()
    periodDurations[[1]] <- c(2, 2)  # 2 days, 2 hrs
    periodDurations[[2]] <- c(3, 2)  # 3 days, 2 hrs
    periodDurations[[3]] <- c(4, 2)  # 4 days, 2 hrs
    periodDurations[[4]] <- c(5, 2)  # 5 days, 2 hrs

    for (periodDuration in periodDurations) {
        period <- periodDuration[1]
        duration <- periodDuration[2]
        smallestPlanetDetectableTest(period, depths, duration, algo=algo, ofac=ofac)
    }
}


# ********** Some important resources **********
# http://www.ccpo.odu.edu/~klinck/Reprints/PDF/omeyHUB2009.pdf (suggested by Suveges, 2014).
# See discussion on period/frequency spacing considerations for BLS: https://johnh2o2.github.io/cuvarbase/bls.html#period-spacing-considerations
# Mathematical description of the Anderson-Darling test: https://bookdown.org/egarpor/NP-UC3M/nptests-dist.html

# Resources for extreme value statistics:
    # (1) http://personal.cityu.edu.hk/xizhou/first-draft-report.pdf
    # (2) Playlist on Extreme Value Statistics: https://youtube.com/playlist?list=PLh35GyCXlQaTJtTq4OQGzMblwEcVIWW9n
    # (3) https://www.lmd.ens.fr/E2C2/class/naveauRomaniaE2C207.pdf

# Any other papers:
# Good set of papers: https://arxiv.org/pdf/1712.00734.pdf

