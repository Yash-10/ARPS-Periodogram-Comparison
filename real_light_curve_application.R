# ***************************************************************************
# Author: Yash Gondhalekar  Last updated: March, 2023

# Description: This script reads all the four DTARPS TESS light curves showed
#              in our paper, runs the extreme value code on each, and then
#              plots the results. This code also provides an example of how
#              to use our code on real light curves.

# ***************************************************************************

source('eva_periodogram.R')
source('bls_and_tcf_periodogram_show.R')

seedValue <- 465

# Read the TESS light curves stored in .txt files.
tables <- list()
filenames <- c()
list_of_files <- list.files("./", pattern=glob2rx("DTARPS*_lc.txt"))
for (file in list_of_files) {
    table <- read.table(file, header=TRUE)
    # NOTE: Our extreme value function `evd` expects time to be in hours. Hence, we modify from JD values to hours.
    table$times = (table$Time - table$Time[1]) * 24
    tables[[file]] <- table
    filenames <- c(filenames, file)
}

# Create an empty dataframe to store results.
results = data.frame(matrix(nrow=0, ncol=10))

# Assign column names
colnames(results) = c('fap_BLS', 'snr_BLS', 'period_BLS (hrs)', 'depth_BLS', 'duration_BLS (hrs)', 'fap_TCF', 'snr_TCF', 'period_TCF (hrs)', 'depth_TCF', 'duration_TCF (hrs)')
# Run extreme value code.
counter <- 1
for (i in 1:length(tables)) {
    table <- tables[[filenames[i]]]
    rt <- c()
    for (algo in c("BLS", "TCF")) {
        if (algo == "BLS") {
            # BLS does not handle non-numeric values in input (NA/NaN/Inf). So we need to manually remove all rows
            # where the observation flux had such non-numeric values.
            # NOTE: It is just that the light curves have NA, so the below condition assumes NA is present.
            # If instead of NA, one has NaN/Inf, it can be changed according to need.
            table_BLS <- table[!is.na(table$Flux),]
        }

        if (algo == "BLS") {
            # Run EVD on the light curve,
            result <- evd(y=table_BLS$Flux, t=table_BLS$times, algo=algo, FAPSNR_mode=0, lctype="real", applyGPRforBLS=TRUE, seedValue=seedValue)
            fap <- result[1]
            perResults <- result[3:5]
            result <- evd(y=table_BLS$Flux, t=table_BLS$times, algo=algo, FAPSNR_mode=1, lctype="real", applyGPRforBLS=TRUE, seedValue=seedValue)
            snr <- 1 / result[1]
            rt <- c(rt, c(fap, snr, perResults))
        }
        else {
            result <- evd(y=table$Flux, t=table$times, algo=algo, FAPSNR_mode=0, lctype="real", seedValue=seedValue)
            fap <- result[1]
            perResults <- result[3:5]
            result <- evd(y=table$Flux, t=table$times, algo=algo, FAPSNR_mode=1, lctype="real", seedValue=seedValue)
            snr <- 1 / result[1]
            rt <- c(rt, c(fap, snr, perResults))
        }
    }
    results[counter,] <- rt
    counter <- counter + 1
}

rownames(results) <- filenames

print("==================================")
print(results)
print("==================================")

write.csv(x=results, file="real_lc_bls_and_tcf_compare_gpr.csv")
