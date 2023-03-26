# ARPS - Optimizing future transiting exoplanet surveys

This repository contains R code for comparing the [BLS](https://www.aanda.org/articles/aa/abs/2002/31/aa2422/aa2422.html) and [TCF](https://iopscience.iop.org/article/10.3847/1538-3881/ab26b8) transit periodogram algorithms. It can be extended for comparing any set of periodograms. It uses extreme value theory and/or signal-to-noise ratio (SNR) for comparison.

## Installation

1. Clone the repository
2. Edit the shared library paths inside `BLS/bls.R` and `TCF3.0/intf_libtcf.R` (locate the line `dyn.load(...)` at the top of these files) to match the path to the `.so` files.

## Usage

The function `evd` inside `eva_periodogram.R` contains the main functionality for calculating the False Alarm Probability and/or the SNR of periodogram peaks. It can handle both simulated and real observational data. No constraints on the observations, such as evenly spaced, no gaps, etc. are imposed.

Here is a basic example to use this function:

```R
source("eva_periodogram.R")  # To source the R script that contains the evd function.
evd(
    2, 0.01, 2, algo="BLS", noiseType=1, ntransits=10,
    ofac=2, L=300, R=300, FAPSNR_mode=0, lctype="sim"
)
```

- This simulates a planet with period = 2 days (first argument), a depth of 0.01% (second argument), and transit duration = 2 hours (third argument).
- `algo` can be either `BLS` or `TCF`.
- `noiseType=1` simulates Gaussian noise with a fixed mean and standard deviation
    - Currently, the user cannot tune it via this function interface, but can be done by editing the source code inside `utils.R` - we recommend the user instead manipulate the `depth` argument since for simulations, only the relative difference between the noise standard deviation and depth matters, and not the actual values.
- `ntransits` controls how many transits must be contained inside the ligt curve.
- `ofac` is the oversampling factor for computing the periodogram. A value around 2-5 generally suffices.
- `L` and `R` are parameters used for the extreme value calculation.
- `FAPSNR_mode` controls what metric must be used to get the periodogram score. 0 means only FAP of the peak is computed. 1 means only the SNR of the peak is computed. 2 means both are computed, in which both FAP and SNR are weighted by fixed factors.
- `lctype` can take values "sim" and "real". The former is to be used for simulations (in which case the period, depth, and duration are needed as input). The latter is to be used for calculations on real ligt curves (in which case, the observations (i.e. fluxes) and time epochs need to be passed. See below).

For passing custom flux values and time epochs (for eg: in the case of real observational data), this can be done by:

```R
evd(y=y, t=t, ...)
```
where `y` and `t` denote the fluxes and time epochs, respectively. Observational fluxes are generally associated with errors, but this cannot be used in the code as of now.

## Example application
The `evd` function can be run on a set of transit depths for different periodogram algorithms independently, to yield plots like the below:

![ntransits_comparison_BLS_and_TCF](images/ntransits_BLS_TCF.png)

## Code motivation
The extreme value part of the code is a replication of the approach described in [Süveges (2014)](https://academic.oup.com/mnras/article/440/3/2099/1077179). This code is not an official implementation of that paper.

## Data availability
Most of our analysis used simulations. The applicability of our method was described on four TESS light curves. These are present as `DTARPS*.txt` files in this repository.

## Bugs or issues
If you find something not working as expected or anything weird, we would like to know and improve it! Please feel free to open an issue in the issue tracker or [send an email](yashgondhalekar567@gmail.com).

## License
[MIT](https://github.com/Yash-10/arps/blob/main/LICENSE)
