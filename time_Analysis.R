# in seconds

tcfColor <- '#FC9272' #'#225ea8' # '#fd8d3c' # '#fdae6b'
blsColor <-  '#1C9099' #'#d95f0e'
titleColor <- 'black'

timeSeriesLengths <- c(279, 543, 1599, 2655, 5295, 10575)
numPeriodsTested <- c(1649.000000, 6467.000000, 57419.000000, 159059.000000, 634919.000000, 2537039.000000)

bls_0.5 <- 17.2232e-3
bls_1 <- 133.8394e-3
bls_3 <- 4.464666
bls_5 <- 25.36005
bls_10 <- 287.7993
bls_20 <- 3693.079
bls_40 <- NA
bls_80 <- NA

tcf_0.5 <- 12.13852e-3
tcf_1 <- 80.76386e-3
tcf_3 <- 1.85525
tcf_5 <- 8.227688
tcf_10 <- 63.28432
tcf_20 <- 495.2101
tcf_40 <- 3938.99
tcf_80 <- 31777.66

blsTimes <- c(bls_0.5, bls_1, bls_3, bls_5, bls_10, bls_20)
tcfTimes <- c(tcf_0.5, tcf_1, tcf_3, tcf_5, tcf_10, tcf_20)

options(scipen=999)
png(filename="time.png", width = 125, height = 125, units='mm', res = 300)

par("mar" = c(5, 5, 2, 2), bg = 'white')

plot(timeSeriesLengths, blsTimes, log='xy', bty='l', type='o', col=blsColor, ylab='Execution time (sec)', xlab='Number of points', lwd=2.0, main='')  # main='Timing Comparison'
axis(1, at=1:6, labels=c(200, 500, 1000, 5000, 10000, 160000))
lines(timeSeriesLengths, tcfTimes, type='o', col=tcfColor, lwd=2.0)

abline(0, 1, lty=1)
abline(a=10, b=0, lty=1)


times <- data.frame(timeSeriesLengths=timeSeriesLengths, tcfTimes=tcfTimes)
model = lm(tcfTimes~timeSeriesLengths, data=times)

p<-data.frame(timeSeriesLengths=timeSeriesLengths)
p$out<-predict(model, p)

#lines(x=p$timeSeriesLengths, y=p$out, col='red')

legend(
	x = "topleft", lty = 1, text.font = 6,
	col= c(blsColor, tcfColor),
	legend=c("BLS", "TCF"), text.col=titleColor, bty = "n", lwd=2.0
)
dev.off()
