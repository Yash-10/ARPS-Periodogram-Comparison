getAveragedMDD <- function(
	vector
) {
	return (
		c(
			mean(unlist(lapply(vector, `[[`, 1))), mean(unlist(lapply(vector, `[[`, 2)))
		)
	)
}

gfntransits2 <- readRDS('GAUSSIAN_ONLY_FAP_MDD_ntransits/result2_gaussian.rds')
gfntransits3 <- readRDS('GAUSSIAN_ONLY_FAP_MDD_ntransits/result3_gaussian.rds')
gfntransits5 <- readRDS('GAUSSIAN_ONLY_FAP_MDD_ntransits/result5_gaussian.rds')
gfntransits10 <- readRDS('GAUSSIAN_ONLY_FAP_MDD_ntransits/result10_gaussian.rds')
gfntransits20 <- readRDS('GAUSSIAN_ONLY_FAP_MDD_ntransits/result20_gaussian.rds')
gfntransits40 <- readRDS('GAUSSIAN_ONLY_FAP_MDD_ntransits/result40_gaussian.rds')
gfntransits60 <- readRDS('GAUSSIAN_ONLY_FAP_MDD_ntransits/result60_gaussian.rds')
gfntransits80 <- readRDS('GAUSSIAN_ONLY_FAP_MDD_ntransits/result80_gaussian.rds')
gfntransits100 <- readRDS('GAUSSIAN_ONLY_FAP_MDD_ntransits/result100_gaussian.rds')

afntransits2 <- readRDS('AR_ONLY_FAP_MDD_ntransits/result2_ar.rds')
afntransits3 <- readRDS('AR_ONLY_FAP_MDD_ntransits/result3_ar.rds')
afntransits5 <- readRDS('AR_ONLY_FAP_MDD_ntransits/result5_ar.rds')
afntransits10 <- readRDS('AR_ONLY_FAP_MDD_ntransits/result10_ar.rds')
afntransits20 <- readRDS('AR_ONLY_FAP_MDD_ntransits/result20_ar.rds')
afntransits40 <- readRDS('AR_ONLY_FAP_MDD_ntransits/result40_ar.rds')
afntransits60 <- readRDS('AR_ONLY_FAP_MDD_ntransits/result60_ar.rds')
afntransits80 <- readRDS('AR_ONLY_FAP_MDD_ntransits/result80_ar.rds')
afntransits100 <- readRDS('AR_ONLY_FAP_MDD_ntransits/result100_ar.rds')

gsntransits2 <- readRDS('GAUSSIAN_ONLY_SNR_MDD_ntransits/result2_gaussian_SNR.rds')
gsntransits3 <- readRDS('GAUSSIAN_ONLY_SNR_MDD_ntransits/result3_gaussian_SNR.rds')
gsntransits5 <- readRDS('GAUSSIAN_ONLY_SNR_MDD_ntransits/result5_gaussian_SNR.rds')
gsntransits10 <- readRDS('GAUSSIAN_ONLY_SNR_MDD_ntransits/result10_gaussian_SNR.rds')
gsntransits20 <- readRDS('GAUSSIAN_ONLY_SNR_MDD_ntransits/result20_gaussian_SNR.rds')
gsntransits40 <- readRDS('GAUSSIAN_ONLY_SNR_MDD_ntransits/result40_gaussian_SNR.rds')
gsntransits60 <- readRDS('GAUSSIAN_ONLY_SNR_MDD_ntransits/result60_gaussian_SNR.rds')
gsntransits80 <- readRDS('GAUSSIAN_ONLY_SNR_MDD_ntransits/result80_gaussian_SNR.rds')
gsntransits100 <- readRDS('GAUSSIAN_ONLY_SNR_MDD_ntransits/result100_gaussian_SNR.rds')

asntransits2 <- readRDS('AR_ONLY_SNR_MDD_ntransits/result2_ar_SNR.rds')
asntransits3 <- readRDS('AR_ONLY_SNR_MDD_ntransits/result3_ar_SNR.rds')
asntransits5 <- readRDS('AR_ONLY_SNR_MDD_ntransits/result5_ar_SNR.rds')
asntransits10 <- readRDS('AR_ONLY_SNR_MDD_ntransits/result10_ar_SNR.rds')
asntransits20 <- readRDS('AR_ONLY_SNR_MDD_ntransits/result20_ar_SNR.rds')
asntransits40 <- readRDS('AR_ONLY_SNR_MDD_ntransits/result40_ar_SNR.rds')
asntransits60 <- readRDS('AR_ONLY_SNR_MDD_ntransits/result60_ar_SNR.rds')
asntransits80 <- readRDS('AR_ONLY_SNR_MDD_ntransits/result80_ar_SNR.rds')
asntransits100 <- readRDS('AR_ONLY_SNR_MDD_ntransits/result100_ar_SNR.rds')

gfntransits2_MDD <- getAveragedMDD(gfntransits2)
gfntransits3_MDD <- getAveragedMDD(gfntransits3)
gfntransits5_MDD <- getAveragedMDD(gfntransits5)
gfntransits10_MDD <- getAveragedMDD(gfntransits10)
gfntransits20_MDD <- getAveragedMDD(gfntransits20)
gfntransits40_MDD <- getAveragedMDD(gfntransits40)
gfntransits60_MDD <- getAveragedMDD(gfntransits60)
gfntransits80_MDD <- getAveragedMDD(gfntransits80)
gfntransits100_MDD <- getAveragedMDD(gfntransits100)

afntransits2_MDD <- getAveragedMDD(afntransits2)
afntransits3_MDD <- getAveragedMDD(afntransits3)
afntransits5_MDD <- getAveragedMDD(afntransits5)
afntransits10_MDD <- getAveragedMDD(afntransits10)
afntransits20_MDD <- getAveragedMDD(afntransits20)
afntransits40_MDD <- getAveragedMDD(afntransits40)
afntransits60_MDD <- getAveragedMDD(afntransits60)
afntransits80_MDD <- getAveragedMDD(afntransits80)
afntransits100_MDD <- getAveragedMDD(afntransits100)

gsntransits2_MDD <- getAveragedMDD(gsntransits2)
gsntransits3_MDD <- getAveragedMDD(gsntransits3)
gsntransits5_MDD <- getAveragedMDD(gsntransits5)
gsntransits10_MDD <- getAveragedMDD(gsntransits10)
gsntransits20_MDD <- getAveragedMDD(gsntransits20)
gsntransits40_MDD <- getAveragedMDD(gsntransits40)
gsntransits60_MDD <- getAveragedMDD(gsntransits60)
gsntransits80_MDD <- getAveragedMDD(gsntransits80)
gsntransits100_MDD <- getAveragedMDD(gsntransits100)

asntransits2_MDD <- getAveragedMDD(asntransits2)
asntransits3_MDD <- getAveragedMDD(asntransits3)
asntransits5_MDD <- getAveragedMDD(asntransits5)
asntransits10_MDD <- getAveragedMDD(asntransits10)
asntransits20_MDD <- getAveragedMDD(asntransits20)
asntransits40_MDD <- getAveragedMDD(asntransits40)
asntransits60_MDD <- getAveragedMDD(asntransits60)
asntransits80_MDD <- getAveragedMDD(asntransits80)
asntransits100_MDD <- getAveragedMDD(asntransits100)

########################################################

tcfColor <- '#FC9272' #'#225ea8' # '#fd8d3c' # '#fdae6b'
blsColor <-  '#1C9099' #'#d95f0e'
titleColor <- 'black'
ablineColor <- 'black'

# Gaussian
glimitingDepthBLSArr_FAPSNR_mode_0 <- c(gfntransits2_MDD[1], gfntransits3_MDD[1], gfntransits5_MDD[1], gfntransits10_MDD[1], gfntransits20_MDD[1], gfntransits40_MDD[1], gfntransits60_MDD[1], gfntransits80_MDD[1], gfntransits100_MDD[1])
glimitingDepthTCFArr_FAPSNR_mode_0 <- c(gfntransits2_MDD[2], gfntransits3_MDD[2], gfntransits5_MDD[2], gfntransits10_MDD[2], gfntransits20_MDD[2], gfntransits40_MDD[2], gfntransits60_MDD[2], gfntransits80_MDD[2], gfntransits100_MDD[2])
glimitingDepthBLSArr_FAPSNR_mode_1 <- c(gsntransits2_MDD[1], gsntransits3_MDD[1], gsntransits5_MDD[1], gsntransits10_MDD[1], gsntransits20_MDD[1], gsntransits40_MDD[1], gsntransits60_MDD[1], gsntransits80_MDD[1], gsntransits100_MDD[1])
glimitingDepthTCFArr_FAPSNR_mode_1 <- c(gsntransits2_MDD[2], gsntransits3_MDD[2], gsntransits5_MDD[2], gsntransits10_MDD[2], gsntransits20_MDD[2], gsntransits40_MDD[2], gsntransits60_MDD[2], gsntransits80_MDD[2], gsntransits100_MDD[2])

# AR
alimitingDepthBLSArr_FAPSNR_mode_0 <- c(afntransits2_MDD[1], afntransits3_MDD[1], afntransits5_MDD[1], afntransits10_MDD[1], afntransits20_MDD[1], afntransits40_MDD[1], afntransits60_MDD[1], afntransits80_MDD[1], afntransits100_MDD[1])
alimitingDepthTCFArr_FAPSNR_mode_0 <- c(afntransits2_MDD[2], afntransits3_MDD[2], afntransits5_MDD[2], afntransits10_MDD[2], afntransits20_MDD[2], afntransits40_MDD[2], afntransits60_MDD[2], afntransits80_MDD[2], afntransits100_MDD[2])
alimitingDepthBLSArr_FAPSNR_mode_1 <- c(asntransits2_MDD[1], asntransits3_MDD[1], asntransits5_MDD[1], asntransits10_MDD[1], asntransits20_MDD[1], asntransits40_MDD[1], asntransits60_MDD[1], asntransits80_MDD[1], asntransits100_MDD[1])
alimitingDepthTCFArr_FAPSNR_mode_1 <- c(asntransits2_MDD[2], asntransits3_MDD[2], asntransits5_MDD[2], asntransits10_MDD[2], asntransits20_MDD[2], asntransits40_MDD[2], asntransits60_MDD[2], asntransits80_MDD[2], asntransits100_MDD[2])

cex <- 1.3 
png(filename="ntransits_BLS_TCF.png", width = 250, height = 140, units='mm', res = 300) 
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
 
plot(glimitingDepthBLSArr_FAPSNR_mode_0, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Gaussian - FAP', xlab='No. of transits', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0,, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(glimitingDepthTCFArr_FAPSNR_mode_0, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex)
legend(
x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:9, labels=c(2, 3, 5, 10, 20, 40, 60, 80, 100), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab='white', col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
abline(h=0.005, lty='dashed', col=ablineColor) 
text(1.4, 0.008, "0.005%", col=ablineColor, cex=cex) 
plot(glimitingDepthBLSArr_FAPSNR_mode_1, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Gaussian - SNR', xlab='No. of transits', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(glimitingDepthTCFArr_FAPSNR_mode_1, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
legend(
x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:9, labels=c(2, 3, 5, 10, 20, 40, 60, 80, 100), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
abline(h=0.005, lty='dashed', col=ablineColor)

#text(1.1, 0.007, "0.005")
 
plot(alimitingDepthBLSArr_FAPSNR_mode_0, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Autoregressive - FAP', xlab='No. of transits', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(alimitingDepthTCFArr_FAPSNR_mode_0, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
legend(
x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:9, labels=c(2, 3, 5, 10, 20, 40, 60, 80, 100), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
#abline(h=0, lty='dashed') 
abline(h=0.005, lty='dashed', col=ablineColor)
 
#text(1.1, 0.007, "0.005")
 
plot(alimitingDepthBLSArr_FAPSNR_mode_1, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Autoregressive - SNR', xlab='No. of transits', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(alimitingDepthTCFArr_FAPSNR_mode_1, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0)
legend(
	x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:9, labels=c(2, 3, 5, 10, 20, 40, 60, 80, 100), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
#abline(h=0, lty='dashed')
 
abline(h=0.005, lty='dashed', col=ablineColor)
 
#text(1.1, 0.007, "0.005")
dev.off()