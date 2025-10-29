# Counterfactual Forecasting for Panel Data (FOCUS)

This Github respository maintains the codes and experiments presented in our paper [**"Counterfactual Forecasting for Panel Data"**](https://github.com/navonildeb/navonildeb.github.io/blob/master/files/focus_preprint.pdf)(under review at AISTATS 2026).
 
---

## Overview

**FOCUS** is a forecasting framework for counterfactual prediction in panel data with temporally correlated latent factors.
This repository provides implementations of FOCUS alongside benchmark methods (**mSSA** and **SyNBEATS**), simulation scripts, and real-data experiments using the HeartSteps study.
Comprehensive experimental details and figure reproduction instructions are provided in the reproducibility vignette (`vignette.Rmd`), from which the accompanying `Vignette.pdf` is generated.  
To regenerate all results and figures from scratch, set the option  
`knitr::opts_chunk$set(echo = TRUE, cache = FALSE)`  
in line 11 of `vignette.Rmd` before knitting.

---

## Setup

### Python dependencies

Some benchmark algorithms require Python. Example installations:
```bash
pip install numpy scipy matplotlib
```
Benchmark repositories:
- [mSSA](https://github.com/AbdullahO/mSSA)
- [SyNBEATS](https://github.com/Crabtain959/SyNBEATS)

### R dependencies
Install necessary R packages and source helper functions:
```r
install.packages("<package_name>")
source("library_causalTS.R")
```

---
## Running Experiments

### Simulation studies
```r
source("gen_data_DGP0.R")     # data generation
source("forecast_DGP0.R")     # forecasting (FOCUS and mSSA)
```

### SyNBEATS results
Run `SyNBEATS_forecast.ipynb` and process outputs with:
```r
source("syn_error_preprocess.R")
```

### HeartSteps study
```r
source("HeartSteps_Experiments/suggestions_clean.R")
source("HeartSteps_Experiments/HeartSteps_forecast.R")
source("HeartSteps_Experiments/HeartSteps_forecast_plot.R")
```

---

## Citation
If you use this code, please cite:
> **Deb, Dwivedi, and Basu (2025).** *Counterfactual Forecasting for Panel Data.* (Under review at AISTATS 2026)

---

## Authors
**Navonil Deb**, Raaz Dwivedi, and Sumanta Basu





