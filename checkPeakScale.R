checkPeakScale <- function(y, t, algo="BLS") {
        # set.seed(1)
        # y <- y + 0.3 * rnorm(length(y))  # scale by 0.3 to keep noise levels low.

        set.seed(1)
        autoRegNoise <- arima.sim(model = list(order=c(1, 0, 1), ar=0.2, ma=0.2), n = length(y))  # It has only AR and MA components, no differencing component. So it is ARMA and not ARIMA. Note: Keep ar and ma < 0.5.
        y <- y + autoRegNoise

	if (algo == "BLS") {
                output <- bls(y, t, bls.plot = FALSE)
                print(max(output$spec))
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
                print(max(output$outpow))
	}
}
