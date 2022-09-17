
library(forecast)

tot_points <- 2000
### gaussian noise
t.noise <- rnorm(tot_points,0,1)
### simulated AR process
t.ar <- arima.sim(n=tot_points,model=list(ar=c(.2,.3,.1),ma=c(.8,.15,.4)))
#t.ar <- arima.sim(n=tot_points,model=list(ar=c(.2)))
### fake background
t <- t.noise #+ t.ar
par(mfrow = c(3,1))
plot(t,type="l")
title(paste("Total Simulated Background --","SD =",round(sd(t),2)," MAD =",round(mad(t),2)))
plot(t.ar,type="l")
title(paste("Autoregressive Progress --","SD =",round(sd(t.ar),2)," MAD =",round(mad(t.ar),2)))
plot(t.noise,type="l")
title(paste("Gaussian Noise --","SD =",round(sd(t.noise),2)," MAD =",round(mad(t.noise),2)))

#t <- t+10
### fake planet transit
p <- t
P <- 311  ## period
A <- 1 ## amplitude
l <- 19   ## duration
p.signal.start <- seq(P-sample(1:(P-1),1),length(t),by=P)#(P+l))
p.signal <- NULL
for (i in p.signal.start){
  p.signal <- c(p.signal,i:min((i+l-1),length(t)))
}
for (i in p.signal.start){
  p.signal <- c(p.signal,i:(i+l-1))
}
p[p.signal] <- p[p.signal]-A
par(mfrow = c(3,1))
plot(t,type="l")#,ylim=y_range)
title("Simulated Background")
plot(p-t,type="l")#,ylim=y_range)
title("Simulated Planet Transit")
plot(p,type="l")#,ylim=y_range)
title("Simulated Signal")
t.d <- c(0,diff(t))
p.d <- c(0,diff(p))

p.fit <- auto.arima(p,d=1,max.p = 3,max.q = 0)
p.ar <- fitted.values(p.fit)
par(mfrow = c(3,1))
plot(t,type="l")#,ylim=y_range)
plot(p,type="l")#,ylim=y_range)
plot(p.ar,type="l")#,ylim=y_range)
###############################
###############################

p.bls = bls(p,1:length(p))

p.ar.bls = bls(p.ar,1:length(p.ar))


dyn.load("eebls.so")
#.Fortran("bar", n=as.integer(5), x=as.double(rnorm(5)))

##eebls(n,t,x,u,v,nf,fmin,df,nb,qmi,qma,
##      p,bper,bpow,depth,qtran,in1,in2)
bls2 <- .Fortran("eebls",
                n = as.integer(length(p)),
                t = as.numeric(1:length(p)),
                x = as.numeric(p),
                u = as.numeric(1:length(p)),
                v = as.numeric(1:length(p)),
                nf = as.integer(1000),
                fmin = as.numeric(0.001),
                df = as.numeric(1e-5),
                nb = as.integer(200),
                qmi = as.numeric(0.01),
                qma = as.numeric(0.1),
                ##
                p = as.numeric(1:1000),
                bper = as.numeric(1),
                bpow = as.numeric(1),
                depth = as.numeric(1),
                qtran = as.numeric(1),
                in1 = as.integer(1),
                in2 = as.integer(1)
                )

f = seq(0.001,by=1e-5,length.out=1000)
per = 1/f
plot(per,bls2$p,type="l")
abline(v=P,col="red")
abline(v=P*2,col="blue")
abline(v=P/2,col="blue")
abline(v=P*3,col="blue")
abline(v=P*2/3,col="blue")
