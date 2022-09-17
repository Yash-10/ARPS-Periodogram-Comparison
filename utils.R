## Utility functions.
library(forecast)

getResidForTCF <- function(
    y  # Time series (must not be differenced because it is done internally).
) {
    max.p = 5
    max.q = 5
    max.d = 0
    ARIMA.fit = auto.arima(diff(y), stepwise=FALSE, approximation=FALSE, seasonal=FALSE, max.p=max.p, max.q=max.q, max.d=max.d) #leave d as 0. 
    # Simple statistics of ARIMA residuals 
    ARIMA.resid = residuals(ARIMA.fit)
    return (ARIMA.resid);
}
