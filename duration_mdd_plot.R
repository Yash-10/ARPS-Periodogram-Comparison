getAveragedMDD <- function(
	vector
) {
	return (
		c(
			mean(unlist(lapply(vector, `[[`, 1))), mean(unlist(lapply(vector, `[[`, 2)))
		)
	)
}

gfperiod05 <- readRDS('GAUSSIAN_ONLY_FAP_MDD_duration/result_dur1_gaussian.rds')
gfperiod1 <- readRDS('GAUSSIAN_ONLY_FAP_MDD_duration/result_dur2_gaussian.rds')
gfperiod4 <- readRDS('GAUSSIAN_ONLY_FAP_MDD_duration/result_dur3_gaussian.rds')
gfperiod7 <- readRDS('GAUSSIAN_ONLY_FAP_MDD_duration/result7_dur4_gaussian.rds')

afperiod05 <- readRDS('AR_ONLY_FAP_MDD_duration/result_dur1_ar.rds')
afperiod1 <- readRDS('AR_ONLY_FAP_MDD_duration/result_dur2_ar.rds')
afperiod4 <- readRDS('AR_ONLY_FAP_MDD_duration/result_dur3_ar.rds')
afperiod7 <- readRDS('AR_ONLY_FAP_MDD_duration/result7_dur4_ar.rds')

gsperiod05 <- readRDS('GAUSSIAN_ONLY_SNR_MDD_duration/result_dur1_gaussian_SNR.rds')
gsperiod1 <- readRDS('GAUSSIAN_ONLY_SNR_MDD_duration/result_dur2_gaussian_SNR.rds')
gsperiod4 <- readRDS('GAUSSIAN_ONLY_SNR_MDD_duration/result_dur3_gaussian_SNR.rds')
gsperiod7 <- readRDS('GAUSSIAN_ONLY_SNR_MDD_duration/result7_dur4_gaussian_SNR.rds')

asperiod05 <- readRDS('AR_ONLY_SNR_MDD_duration/result_dur1_ar_SNR.rds')
asperiod1 <- readRDS('AR_ONLY_SNR_MDD_duration/result_dur2_ar_SNR.rds')
asperiod4 <- readRDS('AR_ONLY_SNR_MDD_duration/result_dur3_ar_SNR.rds')
asperiod7 <- readRDS('AR_ONLY_SNR_MDD_duration/result7_dur4_ar_SNR.rds')

gfperiod05_MDD <- getAveragedMDD(gfperiod05)
gfperiod1_MDD <- getAveragedMDD(gfperiod1)
gfperiod4_MDD <- getAveragedMDD(gfperiod4)
gfperiod7_MDD <- getAveragedMDD(gfperiod7)

afperiod05_MDD <- getAveragedMDD(afperiod05)
afperiod1_MDD <- getAveragedMDD(afperiod1)
afperiod4_MDD <- getAveragedMDD(afperiod4)
afperiod7_MDD <- getAveragedMDD(afperiod7)

gsperiod05_MDD <- getAveragedMDD(gsperiod05)
gsperiod1_MDD <- getAveragedMDD(gsperiod1)
gsperiod4_MDD <- getAveragedMDD(gsperiod4)
gsperiod7_MDD <- getAveragedMDD(gsperiod7)

asperiod05_MDD <- getAveragedMDD(asperiod05)
asperiod1_MDD <- getAveragedMDD(asperiod1)
asperiod4_MDD <- getAveragedMDD(asperiod4)
asperiod7_MDD <- getAveragedMDD(asperiod7)

########################################################

tcfColor <- '#FC9272' #'#225ea8' # '#fd8d3c' # '#fdae6b'
blsColor <-  '#1C9099' #'#d95f0e'
titleColor <- 'black'
ablineColor <- 'black'

# Gaussian
glimitingDepthBLSArr_FAPSNR_mode_0 <- c(gfperiod05_MDD[1], gfperiod1_MDD[1], gfperiod4_MDD[1], gfperiod7_MDD[1])
glimitingDepthTCFArr_FAPSNR_mode_0 <- c(gfperiod05_MDD[2], gfperiod1_MDD[2], gfperiod4_MDD[2], gfperiod7_MDD[2])
glimitingDepthBLSArr_FAPSNR_mode_1 <- c(gsperiod05_MDD[1], gsperiod1_MDD[1], gsperiod4_MDD[1], gsperiod7_MDD[1])
glimitingDepthTCFArr_FAPSNR_mode_1 <- c(gsperiod05_MDD[2], gsperiod1_MDD[2], gsperiod4_MDD[2], gsperiod7_MDD[2])

# AR
alimitingDepthBLSArr_FAPSNR_mode_0 <- c(afperiod05_MDD[1], afperiod1_MDD[1], afperiod4_MDD[1], afperiod7_MDD[1])
alimitingDepthTCFArr_FAPSNR_mode_0 <- c(afperiod05_MDD[2], afperiod1_MDD[2], afperiod4_MDD[2], afperiod7_MDD[2])
alimitingDepthBLSArr_FAPSNR_mode_1 <- c(asperiod05_MDD[1], asperiod1_MDD[1], asperiod4_MDD[1], asperiod7_MDD[1])
alimitingDepthTCFArr_FAPSNR_mode_1 <- c(asperiod05_MDD[2], asperiod1_MDD[2], asperiod4_MDD[2], asperiod7_MDD[2])

cex <- 1.3 
png(filename="duration_BLS_TCF.png", width = 250, height = 140, units='mm', res = 300) 
par("mar" = c(5, 5, 2, 2)) 
mat1 <- matrix(c(
	1, 3,
	1, 3,
	2, 4,
	2, 4), nrow = 4, ncol = 2, byrow = TRUE
)
layout(mat = mat1,
	heights = c(1),    # Heights of the two rows
	widths = c(1)
	)     # Widths of the two columns
 
plot(glimitingDepthBLSArr_FAPSNR_mode_0, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Gaussian - FAP', xlab='Duration (hrs)', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0,, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(glimitingDepthTCFArr_FAPSNR_mode_0, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex)
legend(
x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:4, labels=c(1, 2, 3, 4), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab='white', col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
abline(h=0.005, lty='dashed', col=ablineColor) 
text(1.1, 0.0025, "0.005%", col=ablineColor, cex=cex) 
plot(glimitingDepthBLSArr_FAPSNR_mode_1, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Gaussian - SNR', xlab='Duration (hrs)', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(glimitingDepthTCFArr_FAPSNR_mode_1, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
legend(
x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:4, labels=c(1, 2, 3, 4), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
abline(h=0.005, lty='dashed', col=ablineColor)

#text(1.1, 0.007, "0.005")
 
plot(alimitingDepthBLSArr_FAPSNR_mode_0, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Autoregressive - FAP', xlab='Duration (hrs)', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(alimitingDepthTCFArr_FAPSNR_mode_0, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
legend(
x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:4, labels=c(1, 2, 3, 4), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
#abline(h=0, lty='dashed') 
abline(h=0.005, lty='dashed', col=ablineColor)
 
#text(1.1, 0.007, "0.005")
 
plot(alimitingDepthBLSArr_FAPSNR_mode_1, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Autoregressive - SNR', xlab='Duration (hrs)', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(alimitingDepthTCFArr_FAPSNR_mode_1, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0)
legend(
	x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:4, labels=c(1, 2, 3, 4), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
#abline(h=0, lty='dashed')
 
abline(h=0.005, lty='dashed', col=ablineColor)
 
#text(1.1, 0.007, "0.005")
dev.off()