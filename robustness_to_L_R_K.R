# ***************************************************************************
# Author: Yash Gondhalekar  Last updated: March, 2023

# Description: This script was used for making figures in Appendix B in our
#			   paper, which is about the robustness tests.

# ***************************************************************************

# L
glimitingDepthBLSArr_FAPSNR_mode_0 <- c(0.005918367, 0.005000000, 0.006836735, 0.005918367, 0.005918367)
glimitingDepthTCFArr_FAPSNR_mode_0 <- c(0.01142857, 0.01142857, 0.01234694, 0.01234694, 0.01142857)
glimitingDepthBLSArr_FAPSNR_mode_1 <- c(0.006836735, 0.006836735, 0.006836735, 0.006836735, 0.006836735)
glimitingDepthTCFArr_FAPSNR_mode_1 <- c(0.005, 0.005, 0.005, 0.005, 0.005)

alimitingDepthBLSArr_FAPSNR_mode_0 <- c(0.01969388, 0.01969388, 0.01969388, 0.01969388, 0.01969388)
alimitingDepthTCFArr_FAPSNR_mode_0 <- c(0.03255102, 0.02979592, 0.02979592, 0.02887755, 0.03255102)
alimitingDepthBLSArr_FAPSNR_mode_1 <- c(0.02244898, 0.02244898, 0.02244898, 0.02244898, 0.02244898)
alimitingDepthTCFArr_FAPSNR_mode_1 <- c(0.005, 0.005, 0.005, 0.005, 0.005)

cex <- 1.3 
png(filename="L_robustness.png", width = 250, height = 140, units='mm', res = 300) 
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
 
plot(glimitingDepthBLSArr_FAPSNR_mode_0, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Gaussian - Only FAP', xlab='L', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0,, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(glimitingDepthTCFArr_FAPSNR_mode_0, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex)
legend(
x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:5, labels=c(100,200,300,400,500), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab='white', col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
abline(h=0.005, lty='dashed', col=ablineColor) 
text(1.15, 0.0021, "0.005%", col=ablineColor, cex=cex) 
plot(glimitingDepthBLSArr_FAPSNR_mode_1, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Gaussian - Only SNR', xlab='L', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(glimitingDepthTCFArr_FAPSNR_mode_1, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
legend(
x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:5, labels=c(100,200,300,400,500), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
abline(h=0.005, lty='dashed', col=ablineColor)

#text(1.1, 0.007, "0.005")
 
plot(alimitingDepthBLSArr_FAPSNR_mode_0, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Autoregressive - Only FAP', xlab='L', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(alimitingDepthTCFArr_FAPSNR_mode_0, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
legend(
x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:5, labels=c(100,200,300,400,500), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
#abline(h=0, lty='dashed') 
abline(h=0.005, lty='dashed', col=ablineColor)
 
#text(1.1, 0.007, "0.005")
 
plot(alimitingDepthBLSArr_FAPSNR_mode_1, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Autoregressive - Only SNR', xlab='L', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(alimitingDepthTCFArr_FAPSNR_mode_1, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0)
legend(
	x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:5, labels=c(100,200,300,400,500), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
#abline(h=0, lty='dashed')
 
abline(h=0.005, lty='dashed', col=ablineColor)
 
#text(1.1, 0.007, "0.005")
dev.off()


# R
glimitingDepthBLSArr_FAPSNR_mode_0 <- c(0.005918367, 0.005000000, 0.005000000, 0.005000000, 0.005000000)
glimitingDepthTCFArr_FAPSNR_mode_0 <- c(0.01142857, 0.01142857, 0.01142857, 0.01142857, 0.01142857)
glimitingDepthBLSArr_FAPSNR_mode_1 <- c(0.006836735, 0.006836735, 0.006836735, 0.006836735, 0.006836735)
glimitingDepthTCFArr_FAPSNR_mode_1 <- c(0.005, 0.005, 0.005, 0.005, 0.005)

alimitingDepthBLSArr_FAPSNR_mode_0 <- c(0.01969388, 0.01969388, 0.02061224, 0.01969388, 0.01969388)
alimitingDepthTCFArr_FAPSNR_mode_0 <- c(0.03255102, 0.02979592, 0.03255102, 0.03255102, 0.02979592)
alimitingDepthBLSArr_FAPSNR_mode_1 <- c(0.02244898, 0.02244898, 0.02244898, 0.02244898, 0.02244898)
alimitingDepthTCFArr_FAPSNR_mode_1 <- c(0.005, 0.005, 0.005, 0.005, 0.005)

cex <- 1.3
png(filename="R_robustness.png", width = 250, height = 140, units='mm', res = 300)
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

plot(glimitingDepthBLSArr_FAPSNR_mode_0, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Gaussian - Only FAP', xlab='R', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0,, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(glimitingDepthTCFArr_FAPSNR_mode_0, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex)
legend(
x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:5, labels=c(100,200,300,400,500), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab='white', col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
abline(h=0.005, lty='dashed', col=ablineColor)
text(1.15, 0.0021, "0.005%", col=ablineColor, cex=cex)
plot(glimitingDepthBLSArr_FAPSNR_mode_1, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Gaussian - Only SNR', xlab='R', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(glimitingDepthTCFArr_FAPSNR_mode_1, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
legend(
x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:5, labels=c(100,200,300,400,500), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
abline(h=0.005, lty='dashed', col=ablineColor)

#text(1.1, 0.007, "0.005")

plot(alimitingDepthBLSArr_FAPSNR_mode_0, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Autoregressive - Only FAP', xlab='R', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(alimitingDepthTCFArr_FAPSNR_mode_0, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
legend(
x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:5, labels=c(100,200,300,400,500), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
#abline(h=0, lty='dashed')
abline(h=0.005, lty='dashed', col=ablineColor)

#text(1.1, 0.007, "0.005")

plot(alimitingDepthBLSArr_FAPSNR_mode_1, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Autoregressive - Only SNR', xlab='R', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(alimitingDepthTCFArr_FAPSNR_mode_1, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0)
legend(
	x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:5, labels=c(100,200,300,400,500), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
#abline(h=0, lty='dashed')

abline(h=0.005, lty='dashed', col=ablineColor)

#text(1.1, 0.007, "0.005")
dev.off()

# K
glimitingDepthBLSArr_FAPSNR_mode_0 <- c(0.005000000, 0.005918367, 0.005000000, 0.005000000, 0.005000000)
glimitingDepthTCFArr_FAPSNR_mode_0 <- c(0.01418367, 0.01142857, 0.01142857, 0.01142857, 0.01142857)
glimitingDepthBLSArr_FAPSNR_mode_1 <- c(0.006836735, 0.006836735, 0.006836735, 0.006836735, 0.006836735)
glimitingDepthTCFArr_FAPSNR_mode_1 <- c(0.005, 0.005, 0.005, 0.005, 0.005)

alimitingDepthBLSArr_FAPSNR_mode_0 <- c(0.01969388, 0.01969388, 0.01969388, 0.02061224, 0.01969388)
alimitingDepthTCFArr_FAPSNR_mode_0 <- c(0.03255102, 0.03255102, 0.03255102, 0.03071429, 0.03071429)
alimitingDepthBLSArr_FAPSNR_mode_1 <- c(0.00500000, 0.02244898, 0.02336735, 0.02520408, 0.02612245)
alimitingDepthTCFArr_FAPSNR_mode_1 <- c(0.005, 0.005, 0.005, 0.005, 0.005)

cex <- 1.3
png(filename="K_robustness.png", width = 250, height = 140, units='mm', res = 300)
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

plot(glimitingDepthBLSArr_FAPSNR_mode_0, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Gaussian - Only FAP', xlab='K', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0,, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(glimitingDepthTCFArr_FAPSNR_mode_0, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex)
legend(
x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:5, labels=c(1,2,3,4,5), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab='white', col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
abline(h=0.005, lty='dashed', col=ablineColor)
text(1.15, 0.0021, "0.005%", col=ablineColor, cex=cex)
plot(glimitingDepthBLSArr_FAPSNR_mode_1, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Gaussian - Only SNR', xlab='K', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(glimitingDepthTCFArr_FAPSNR_mode_1, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
legend(
x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:5, labels=c(1,2,3,4,5), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
abline(h=0.005, lty='dashed', col=ablineColor)

#text(1.1, 0.007, "0.005")

plot(alimitingDepthBLSArr_FAPSNR_mode_0, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Autoregressive - Only FAP', xlab='K', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(alimitingDepthTCFArr_FAPSNR_mode_0, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
legend(
x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:5, labels=c(1,2,3,4,5), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
#abline(h=0, lty='dashed')
abline(h=0.005, lty='dashed', col=ablineColor)

#text(1.1, 0.007, "0.005")

plot(alimitingDepthBLSArr_FAPSNR_mode_1, type='o', xaxt='n', col=blsColor, ylab='Minimum detectable depth (%)', ylim=c(0.0, 0.055), main='Autoregressive - Only SNR', xlab='K', col.main=titleColor, col.lab=titleColor, col.axis=titleColor, bty='l', lwd=2.0, cex.axis=cex, cex.lab=cex, cex=cex, cex.main=cex)
lines(alimitingDepthTCFArr_FAPSNR_mode_1, type='o', col=tcfColor, col.lab=titleColor, col.axis=titleColor, col.axis=titleColor, lwd=2.0)
legend(
	x = "topright", lty = 1, text.font = 6,
col= c(blsColor, tcfColor),
legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0, cex=cex
)
axis(1, at=1:5, labels=c(1,2,3,4,5), col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
axis(2, col.lab=titleColor, col.axis=titleColor, col.ticks=titleColor, cex.axis=cex)
#abline(h=0, lty='dashed')

abline(h=0.005, lty='dashed', col=ablineColor)

#text(1.1, 0.007, "0.005")
dev.off()
