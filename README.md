# Periodogram Comparison for Optimizing Small Transiting Planet Detection

[![DOI](https://zenodo.org/badge/509952217.svg)](https://zenodo.org/badge/latestdoi/509952217) [![Notebook example](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1vTNkruatH5Uk4mGoT1WT1K5EU8ej4bca?usp=sharing)


An approach for optimizing the detection of small transiting planets in future transiting exoplanet surveys ([arXiv:2308.04282](https://arxiv.org/abs/2308.04282))

---

**Motivation for this work**\
A major goal of space-based exoplanet transit survey missions (Kepler, K2, TESS, Plato, Roman) is to detect small Earth-sized planets that produce very faint (~0.01%) periodic dips in stellar brightness. This requires an extremely effective reduction of aperiodic stellar and instrumental variations and a very sensitive periodogram applied to the detrended data. Our recent paper (Gondhalekar, Feigelson, Montalto, Caceres & Saha, 2023, submitted to the AAS Journals) provides a detailed analysis of two statistical approaches to this problem. We find that ARIMA detrending and the Transit Comb Filter periodogram is substantially more effective for small planet detection than commonly used detrenders (e.g., splines, Gaussian Processes regression) followed by Box Least Squares periodogram. This AutoRegressive Planet Search (ARPS) approach, developed by Caceres et al. 2019a, has been used to find small transiting planets from the 4-year Kepler dataset (Caceres et al. 2019b) and the Year 1 TESS Full Field Image dataset (Melton et al. 2023a,b,c). 

**Methodology and coding**\
This GitHub repository contains R code for comparing the BLS and TCF transit periodogram algorithms. It can be extended for comparing any set of periodograms. It uses extreme value theory and/or signal-to-noise ratio (SNR) for comparison. The R language (with computationally intensive periodogram computations in Fortran) is preferred over Python because of its comprehensive packages for time series analysis and other statistical methodologies.

Several detrending approaches have been used in the past for detrending stellar light curves prior to the search for transiting exoplanets: Gaussian Processes Regression or local regression methods such as Splines and LOWESS. These approaches do an excellent job of removing (generally systematic but can also be stellar) long-term trends from the light curve. The class of parametric AutoRegressive Moving Average (ARMA) models is effective in removing (generally stellar) short-term trends. One well-known issue is that detrending methods often alter the planet's depth since, in its naive usage, the approach has no information about whether and where the transit is present (transit masking, for example, can be used to resolve this). Hence, ARMA used before BLS was not helpful since ARMA advertently modeled the planetary signal.

Our study observed that the ARMA model is beneficial if fitted on a differenced light curve. Once the light curve is differenced, BLS can no longer be used since the box-shaped transit is changed to a double spike. Hence, one can use the TCF periodogram after differencing + ARMA modeling. The resulting pipeline has proven more beneficial than Gaussian Processes Regression/Splines + BLS for detecting small planets. The TCF periodogram also possesses better properties (noise characteristics, peak width) than BLS.

Please see our [paper](https://arxiv.org/abs/2308.04282) for more details.

## Repository description

- The `BLS` and `TCF3.0` folders contain the underlying implementation of these algorithms.
	- The BLS Fortran code is taken from the [original BLS website](https://konkoly.hu/staff/kovacs/bls_code.html) (the modification to handle edge effects, `eebls.f` is used). The folder `BLS` also contains `bls.R`, which is the R interface for the Fortran code. We recommend users use this script for running BLS.
	- The `TCF3.0` folder contains the TCF implementation. The lower-level implementation is again in Fortran (Fortran was used due to its swift computation). The `tcf.f95` code is the cleaned version after consulting with the original code author, Gabriel A. Caceres. We recommend using this version of TCF for better compatibility with the rest of the code in this repository.
	- Note that both these implementations are not newly introduced in this repository.
- The `images` folder contains `.png` figures used in our paper.
- The `results` folder contains the data used for plotting Figures 5 and 6 in our paper.
- `DTARPS*.txt` contain four real TESS light curves used in the study (see our paper for details).
- The rest of the repository contains various scripts (mainly in R, but one in Python). The most important script is `eva_periodogram.R` since it contains the code for running the significance analysis procedure. Note that it depends on the following scripts: `BLS/bls.R`, `TCF3.0/intf_libtcf.R`, and `utils.R`. `utils.R` contains various utility functions. The code contains comments, so it might be easy to understand each function's usage.
- The notebook `extra_figure_why_tcf_better.ipynb` contains Python code for generating Figure 11 from our paper.
- The script `bls_and_tcf_noise_change_as_depth_decrease.R` was used to generate Figures 1-4 and 7-10 in the paper.
- `python_utils.py` contains a standalone Python utility function (it was created since I was not able to implement it in R)
- `real_light_curve_application.R` is kept in this repository to showcase how this code can be used on real light curves. For a complete example, see the function `blsAndTCFDepthChangeReal` in the script `bls_and_tcf_noise_change_as_depth_decrease.R` here.
- `time_Analysis.R` contains the script used for plotting Figure 13 in the paper. The data in the script was obtained by running the (long) computations on the [Kaggle](https://www.kaggle.com/) kernel.
- Any other remaining scripts not described here are perhaps not necessary but have been kept for storage and more context.

---

## Installation

1. Clone the repository.
2. Create shared object files by compiling the fortran source code:


This can be done by using:

```bash
# For BLS
R CMD SHLIB eebls.f -o a.out  # will create a.out
# For TCF
R CMD SHLIB main_tcf.f95 median.f90 rand_tools.f95 tcf.f95 -o a.out  # will create a.out
```

Alternatively, one can achieve the same using a set of terminal commands ([`gfortran`](https://fortran-lang.org/en/learn/os_setup/install_gfortran/) must be installed for this):
```bash
# For BLS
cd BLS
gfortran -c eebls.f  # will create eebls.o
gfortran -shared eebls.o  # will create a.out

# For TCF
cd ../TCF3.0
gfortran -c median.f90
gfortran -c rand_tools.f95
gfortran -c tcf.f95
gfortran -c main_tcf.f95

# Combine all object files into a single shared object file
gfortran -shared median.o rand_tools.o tcf.o main_tcf.o  # will create a.out
```

3. Edit the shared library paths inside `BLS/bls.R` and `TCF3.0/intf_libtcf.R` (locate the line `dyn.load(...)` at the top of these files) to match the path to the BLS's `a.out` and TCF's `a.out` files.

**Note**:

1. While BLS's and TCF's (backend) implementation is written in Fortran, their usage is facilitated in R (frontend). That is why the above procedure is used.
2. The installation instructions here have been tested only on the Ubuntu OS.


## Usage

### Simulations

The function `evd` inside `eva_periodogram.R` contains the main functionality for calculating the False Alarm Probability and/or the SNR of periodogram peaks. It can handle both simulated and real observational data. No constraints on the observations, such as evenly spaced, no gaps, etc., are imposed.

Here is a basic example of using this function:

```R
source("eva_periodogram.R")  # To source the R script that contains the evd function.
result <- evd(
    2, 0.01, 2, algo="BLS", noiseType=1, ntransits=10,
    ofac=2, L=300, R=300, FAPSNR_mode=0, lctype="sim"
)
score <- result[1]  # score is the FAP if FAPSNR_mode == 0, SNR if FAPSNR_mode == 1, and a weighted sum of FAP and SNR if FAPSNR_mode == 2.
```

- This simulates a planet with a period = 2 days (first argument), a depth of 0.01% (second argument), and transit duration = 2 hours (third argument).
- `algo` can be either `BLS` or `TCF`.
- `noiseType=1` simulates Gaussian noise with a fixed mean and standard deviation.
    - Currently, the user cannot tune the Gaussian noise parameters via this function interface, but can be done by editing the source code inside `utils.R` - we recommend the user instead manipulate the `depth` argument since only the relative difference between the noise standard deviation and depth matters, and not the actual values.
- `ntransits` controls how many transits must be contained inside the light curve.
- `ofac` is the oversampling factor for computing the periodogram. A value of around 2-5 generally suffices.
- `L` and `R` are parameters used for the extreme value calculation.
- `FAPSNR_mode` controls what metric must be used to get the periodogram score. 0 means only the FAP of the peak is computed. 1 means only the SNR of the peak is computed. 2 means both are computed, in which both FAP and SNR are weighted by fixed factors.
- `lctype` can take values "sim" and "real". The former is to be used for simulations (in which case the period, depth, and duration are needed as input). The latter is to be used for calculations on real light curves (in which case, the observations (i.e., fluxes) and time epochs need to be passed. See below).

### Custom datasets

For passing custom flux values and time epochs (e.g., in the case of real observational data), this can be done by:

```R
result <- evd(y=y, t=t, ...)
```
where `y` and `t` denote the fluxes and time epochs, respectively. `...` denote any other arguments to `evd` as stated above. Observational fluxes are generally associated with errors, but this cannot be used in the code as of now.

**Important notes**:

1. The BLS code present in this repository (which is extracted directly from the [original source code](https://ui.adsabs.harvard.edu/abs/2016ascl.soft07008K/abstract)), fails when non-finite (e.g., NaN/Inf) values are present in the input. Hence, one needs to manually remove the flux values and the corresponding observation epoch that are non-finite before running the extreme value code above. If you think the extreme value code should handle non-finite value check internally so that one need not do it beforehand, please open an issue, and we can can discuss it further.
2. Noting the above point, [`real_light_curve_application.R`](https://github.com/Yash-10/arps/blob/main/real_light_curve_application.R) contains a code example of how to use the extreme value code for BLS without any errors. See [line 41 in `real_light_curve_application.R`](https://github.com/Yash-10/arps/blob/main/real_light_curve_application.R#L41).

## Example application
The `evd` function can be run on a set of transit depths for different periodogram algorithms independently to yield the minimum detectable depth. These can be used to make plots like the below:

![ntransits_comparison_BLS_and_TCF](images/ntransits_BLS_TCF.png)

For more details, please see our paper.

## Code motivation
The extreme value part of the code is a replication of the approach described in [Süveges (2014)](https://academic.oup.com/mnras/article/440/3/2099/1077179). This code is not an official implementation of that paper.

## Data availability
Most of our analysis used simulations. The applicability of our method was described on four TESS light curves. These are present as `DTARPS*.txt` files in this repository.

## Bugs or issues
If you find something not working as expected or anything weird, we would like to know and improve it! Please feel free to open an issue in the issue tracker or [send an email](yashgondhalekar567@gmail.com).

## Current status (as of March 2023)
The code was tested using some internal quick tests, which were designed while writing the code. These tests are not in this repository. It does not have a dedicated `tests/` folder where ideally all the test scripts can be written. This is planned for the future. The code, as of now, is not optimized for speed. However, it is expected that it runs without any troubles in most use cases. Future versions of the code might focus on these issues.

A Python version of the comparison code is aimed to be released sometime in the future for wider accessibility. If you are interested in this, please let us know your thoughts at https://github.com/Yash-10/ARPS-Periodogram-Comparison/issues/2.

## License
[MIT](https://github.com/Yash-10/arps/blob/main/LICENSE)
