# ***************************************************************************
# Author: Yash Gondhalekar  Last updated: March, 2023

# Description: This script contains code for creating figures 3 and 4 from
#			   our paper. Values used for plotting are hardcoded inside the
# 			   script since they were separately calculated. See the Google
#			   Colab example tutorial for details on how to get these values
#			   in the first place.

# ***************************************************************************

# NTRANSITS
tcfColor <- '#FC9272' #'#225ea8' # '#fd8d3c' # '#fdae6b'
blsColor <-  '#1C9099' #'#d95f0e'
titleColor <- 'black'
ablineColor <- 'black'

# Gaussian
glimitingDepthBLSArr_FAPSNR_mode_0 <- c(0.024010204, 0.043295918, 0.024010204, 0.013448980, 0.010785714, 0.007479592, 0.005642857, 0.005091837, 0.005000000, 0.005000000, 0.005000000)
glimitingDepthTCFArr_FAPSNR_mode_0 <- c(0.050000000, 0.050000000, 0.050000000, 0.043846939, 0.024469388, 0.013540816, 0.008581633, 0.006010204, 0.005091837, 0.005091837, 0.005000000)
glimitingDepthBLSArr_FAPSNR_mode_1 <- c(0.025663265, 0.048438776, 0.025663265, 0.015377551, 0.011244898, 0.007204082, 0.005459184, 0.005091837, 0.005000000, 0.005000000, 0.005000000)
glimitingDepthTCFArr_FAPSNR_mode_1 <- c(0.006469388, 0.013448980, 0.006469388, 0.006010204, 0.005000000, 0.005000000, 0.005000000, 0.005000000, 0.005000000, 0.005000000, 0.005000000)

# AR
alimitingDepthBLSArr_FAPSNR_mode_0 <- c(0.05000000, 0.05000000, 0.05000000, 0.04210204, 0.02924490, 0.01767347, 0.01721429, 0.01436735, 0.01280612, 0.01152041, 0.01078571)
alimitingDepthTCFArr_FAPSNR_mode_0 <- c(0.050000000, 0.050000000, 0.050000000, 0.050000000, 0.044030612, 0.028234694, 0.017489796, 0.011704082, 0.009683673, 0.008030612, 0.007479592)
alimitingDepthBLSArr_FAPSNR_mode_1 <- c(0.04687755, 0.05000000, 0.04687755, 0.03135714, 0.02685714, 0.01978571, 0.01501020, 0.01409184, 0.01170408, 0.01014286, 0.00922449)
alimitingDepthTCFArr_FAPSNR_mode_1 <- c(0.009132653, 0.005551020, 0.009132653, 0.005826531, 0.005000000, 0.005000000, 0.005000000, 0.005000000, 0.005000000, 0.005, 0.005)

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
axis(1, at=1:11, labels=c(0, 1, 2,3, 5, 10, 20, 40, 60, 80, 100), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
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
axis(1, at=1:11, labels=c(0, 1, 2,3, 5, 10, 20, 40, 60, 80, 100), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
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
axis(1, at=1:11, labels=c(0, 1, 2,3, 5, 10, 20, 40, 60, 80, 100), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
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
axis(1, at=1:11, labels=c(0, 1, 2,3, 5, 10, 20, 40, 60, 80, 100), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
#abline(h=0, lty='dashed')
 
abline(h=0.005, lty='dashed', col=ablineColor)
 
#text(1.1, 0.007, "0.005")
dev.off()

# PERIOD
glimitingDepthBLSArr_FAPSNR_mode_0 <- c(0.009132653, 0.007479592, 0.007571429, 0.009132653)
glimitingDepthTCFArr_FAPSNR_mode_0 <- c(0.01960204, 0.01354082, 0.01096939, 0.01280612)
glimitingDepthBLSArr_FAPSNR_mode_1 <- c(0.009224490, 0.007204082, 0.006102041, 0.005551020)
glimitingDepthTCFArr_FAPSNR_mode_1 <- c(0.005091837, 0.005000000, 0.005000000, 0.005000000)

alimitingDepthBLSArr_FAPSNR_mode_0 <- c(0.01400000, 0.01767347, 0.02483673, 0.02997959)
alimitingDepthTCFArr_FAPSNR_mode_0 <- c(0.02896939, 0.02823469, 0.02198980, 0.02079592)
alimitingDepthBLSArr_FAPSNR_mode_1 <- c(0.01014286, 0.01978571, 0.03264286, 0.0382449)
alimitingDepthTCFArr_FAPSNR_mode_1 <- c(0.005734694, 0.005000000, 0.005000000, 0.005)

cex <- 1.3
png(filename="period_BLS_TCF.png", width = 250, height = 140, units='mm', res = 300)
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

plot(glimitingDepthBLSArr_FAPSNR_mode_0, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Gaussian - FAP', xlab='Period (days)', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0,, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(glimitingDepthTCFArr_FAPSNR_mode_0, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex)
legend(
x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:4, labels=c(0.5, 1, 4, 7), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab='white', col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
abline(h=0.005, lty='dashed', col=ablineColor)
text(1.1, 0.0025, "0.005%", col=ablineColor, cex=cex)
plot(glimitingDepthBLSArr_FAPSNR_mode_1, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Gaussian - SNR', xlab='Period (days)', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(glimitingDepthTCFArr_FAPSNR_mode_1, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
legend(
x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:4, labels=c(0.5, 1, 4, 7), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
abline(h=0.005, lty='dashed', col=ablineColor)

#text(1.1, 0.007, "0.005")

plot(alimitingDepthBLSArr_FAPSNR_mode_0, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Autoregressive - FAP', xlab='Period (days)', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(alimitingDepthTCFArr_FAPSNR_mode_0, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
legend(
x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:4, labels=c(0.5, 1, 4, 7), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
#abline(h=0, lty='dashed')
abline(h=0.005, lty='dashed', col=ablineColor)

#text(1.1, 0.007, "0.005")

plot(alimitingDepthBLSArr_FAPSNR_mode_1, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Autoregressive - SNR', xlab='Period (days)', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(alimitingDepthTCFArr_FAPSNR_mode_1, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0)
legend(
	x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:4, labels=c(0.5, 1, 4, 7), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
#abline(h=0, lty='dashed')

abline(h=0.005, lty='dashed', col=ablineColor)

#text(1.1, 0.007, "0.005")
dev.off()


# DURATION
glimitingDepthBLSArr_FAPSNR_mode_0 <- c(0.010142857, 0.007479592, 0.006469388, 0.005826531)
glimitingDepthTCFArr_FAPSNR_mode_0 <- c(0.01932653, 0.01354082, 0.01142857, 0.01188776)
glimitingDepthBLSArr_FAPSNR_mode_1 <- c(0.010418367, 0.007204082, 0.006653061, 0.006010204)
glimitingDepthTCFArr_FAPSNR_mode_1 <- c(0.005, 0.005, 0.005, 0.005)

alimitingDepthBLSArr_FAPSNR_mode_0 <- c(0.02254082, 0.01767347, 0.01804082, 0.01638776)
alimitingDepthTCFArr_FAPSNR_mode_0 <- c(0.02942857, 0.02823469, 0.02713265, 0.02979592)
alimitingDepthBLSArr_FAPSNR_mode_1 <- c(0.02520408, 0.01978571, 0.01243878, 0.01115306)
alimitingDepthTCFArr_FAPSNR_mode_1 <- c(0.005, 0.005, 0.005, 0.005)

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
axis(1, at=1:4, labels=c(1,2,3,4), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
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
axis(1, at=1:4, labels=c(1,2,3,4), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
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
axis(1, at=1:4, labels=c(1,2,3,4), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
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
axis(1, at=1:4, labels=c(1,2,3,4), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
#abline(h=0, lty='dashed')

abline(h=0.005, lty='dashed', col=ablineColor)

#text(1.1, 0.007, "0.005")
dev.off()
