# install.packages("extRemes")
# install.packages("cobs")
# install.packages("moments")
# install.packages("goftest")
# install.packages("gbutils")
# install.packages("forecast")
# install.packages("reticulate")

source("eva_periodogram.R")

getScore <- function(depth, ...) {
  result <- try(
      suppressWarnings(evd(depth=depth, ...))
  )
  result <- if (inherits(result, "try-error")) c(NA, NA) else result 
  score <- result[1]
}

##################### LIST OF PARAMETERS #####################
L = 300
R = 300
period = 2
depth = 0.02336735
duration = 2
ntransits = 10
useStandardization = TRUE
mode = 'detrend_normalize'
noiseType = 1
useOptimalFreqSampling = TRUE
ofac = 2
res = 2
FAPSNR_mode = 0
significanceMode = 'max'
##############################################################

# start_time <- Sys.time()

# Generate random seeds.
thousandSeeds <- sample(1:1e5, 10)

# Create vector to store BLS and TCF FAPs.
BLS_FAPs <- c()
TCF_FAPs <- c()

# Run FAP procedure for each seed [BLS and TCF].
for (seedValue in thousandSeeds) {
    # BLS
    msgs <- capture.output(
        scoreBLS <- getScore(depth, L=L, R=R, period=period, duration=duration, noiseType=noiseType,
            algo="BLS", ntransits=ntransits, significanceMode=significanceMode,
            useStandardization=useStandardization, ofac=ofac, res=res, useOptimalFreqSampling=TRUE,
            checkConditions=TRUE, plot=FALSE, seedValue=seedValue, FAPSNR_mode=FAPSNR_mode, mode=mode,
        )
    )

    # TCF
    msgs <- capture.output(
        scoreTCF <- getScore(depth, L=L, R=R, period=period, duration=duration, noiseType=noiseType,
            algo="TCF", ntransits=ntransits, significanceMode=significanceMode,
            useStandardization=useStandardization, ofac=ofac, res=res, useOptimalFreqSampling=TRUE,
            checkConditions=TRUE, plot=FALSE, seedValue=seedValue, FAPSNR_mode=FAPSNR_mode, mode=mode,
        )
    )

    BLS_FAPs <- c(BLS_FAPs, scoreBLS)
    TCF_FAPs <- c(TCF_FAPs, scoreTCF)
}

# Save the FAPs.
# One can read it after saving using `readRDS(...)`.
saveRDS(BLS_FAPs, file="BLS_FAPs.rds")
saveRDS(TCF_FAPs, file="TCF_FAPs.rds")

# end_time <- Sys.time()

# print(sprintf("Time required: %.1fs", end_time - start_time))
