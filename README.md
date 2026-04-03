# har-xgb-block-conformal
Hybrid HAR + XGBoost + block-conformal pipeline for daily variance-proxy forecasting
## Overview

This repository provides an R implementation of a conservative hybrid forecasting framework for daily variance-proxy prediction. The pipeline combines:

- a heterogeneous autoregressive (HAR) benchmark,
- an XGBoost nonlinear forecasting model,
- a convex stacking layer that combines HAR and XGBoost forecasts,
- standard split-conformal prediction bands,
- block-conformal prediction bands based on block maxima.

The empirical target in the main specification is a **daily variance proxy** constructed from squared daily returns. The framework is designed to be transparent, reproducible, and suitable for both point forecasting and uncertainty quantification under volatility clustering.

## Asset universe

The main empirical analysis uses the following six assets:

- XOM
- EFA
- XLB
- FDX
- BTC-USD
- BRK-B

## Sample period

- Start date: 2015-01-01
- End date: 2024-12-31

## Main outputs

### Empirical tables
- **T0** Descriptive statistics
- **T1** RMSE comparison
- **T2** Coverage comparison
- **T3** Band-width comparison
- **T4** Tail-region performance
- **T5** XGBoost tuning summary
- **T6** Block-length sensitivity

### Empirical figures
- **F1** Test-period variance proxy and forecasts
- **F2** Standard vs block-conformal bands
- **F4** Calibration nonconformity score histograms

### Simulation outputs
- **S1** Simulation coverage table
- **S2** Simulation width table
- **F5** Simulation trade-off figure

### Optional appendix robustness
- **T8** Alternative variance-proxy robustness

## Repository contents

Recommended repository structure:

```text
.
├── final_har_xgb_block_conformal.R
├── README.md
├── .gitignore
└── outml_final/
    ├── plots/
    ├── tables/
    └── rds/
