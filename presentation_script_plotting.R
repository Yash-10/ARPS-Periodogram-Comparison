png("foo.png", width=15, height=10)

par("mar" = c(5, 6, 4, 2))
cexVal = 1.7
mat1 <- matrix(c(
    1, 2, 3,
    4, 5, 6,
    7, 8, 9), nrow = 3, ncol = 3, byrow = TRUE
)
layout(mat = mat1,
    heights = c(1),    # Heights of the two rows
    widths = c(1)
)
depths <- c(
  seq(from = 0.005, to = 0.03, length.out=30),
  seq(from = 0.031, to = 0.1, length.out=10)
)
plot(depths, bls1, type='o', main='Period = 1 day', ylab='FAP, duration = 2 hr', xlab='depths (%)', ylim=c(1e-10, 0.03), col='#1C9099', log='x', cex.main=2, cex.axis=1.5, cex.lab=2)
lines(depths, tcf1, type='o', col='#FC9272', cex.main=2, cex.axis=1.5, cex.lab=2)
abline(h=0.01, col='black', lty=2)
legend("topright", lty = 1, text.font = 6, col = c('#1C9099', '#FC9272'), text.col='black', legend=c('BLS', 'TCF'), cex=1.5)
plot(depths, bls7, type='o', main='Period = 7 days', ylab='FAP, duration = 2 hr', xlab='depths (%)', ylim=c(1e-10, 0.03), col='#1C9099', log='x', cex.main=2, cex.axis=1.5, cex.lab=2)
lines(depths, tcf7, type='o', col='#FC9272', cex.main=2, cex.axis=1.5, cex.lab=2)
abline(h=0.01, col='black', lty=2)
plot(depths, bls13, type='o', main='Period = 13 days', ylab='FAP, duration = 2 hr', xlab='depths (%)', ylim=c(1e-10, 0.03), col='#1C9099', log='x', cex.main=2, cex.axis=1.5, cex.lab=2)
lines(depths, tcf13, type='o', col='#FC9272', cex.main=2, cex.axis=1.5, cex.lab=2)
abline(h=0.01, col='black', lty=2)

plot(depths, bls13, type='o', main='Period = 1 days', ylab='FAP, duration = 3 hr', xlab='depths (%)', ylim=c(1e-10, 0.03), col='#1C9099', log='x', cex.main=2, cex.axis=1.5, cex.lab=2)
lines(depths, tcf13, type='o', col='#FC9272', cex.main=2, cex.axis=1.5, cex.lab=2)
abline(h=0.01, col='black', lty=2)
legend("topright", lty = 1, text.font = 6, col = c('#1C9099', '#FC9272'), text.col='black', legend=c('BLS', 'TCF'), cex=1.5)
plot(depths, bls73, type='o', main='Period = 7 days', ylab='FAP, duration = 3 hr', xlab='depths (%)', ylim=c(1e-10, 0.03), col='#1C9099', log='x', cex.main=2, cex.axis=1.5, cex.lab=2)
lines(depths, tcf73, type='o', col='#FC9272', cex.main=2, cex.axis=1.5, cex.lab=2)
abline(h=0.01, col='black', lty=2)
plot(depths, bls133, type='o', main='Period = 13 days', ylab='FAP, duration = 3 hr', xlab='depths (%)', ylim=c(1e-10, 0.03), col='#1C9099', log='x', cex.main=2, cex.axis=1.5, cex.lab=2)
lines(depths, tcf133, type='o', col='#FC9272', cex.main=2, cex.axis=1.5, cex.lab=2)
abline(h=0.01, col='black', lty=2)


dev.off()
