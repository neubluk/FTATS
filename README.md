# FTATS

Forecasting Temporally Aggregated Time Series

[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

This R package features the essential functions from [Neubauer and Filzmoser (2024)](arxiv). Exported functions are

- `fit_agg_arima`: Fit fixed or automatically selected ARIMA models to a temporal hierarchy.
- `reconcile_forecasts`: Reconcile base forecasts with reconciliation method of choice.
- `reconcile_forecasts.cv`: Time-Series Cross-Validation function reconciliation methods with hyperparameters.
- `forecast_and_reconcile`: Combines fitting ARIMA models and reconciling resulting base forecasts.

## Citation

## Installation

``` s
# install.packages("remotes")
remotes::install_github("neubluk/FTATS")
```

## Data Examples

Data examples are available in `demo/data_example.R`, `demo/data_example2.R`, `demo/data_examples_all.R`, as well as simulation studies in `demo/simulation.R`

## Data

Four datasets are available.

- Daily Food Demand Data from Schrankel GmbH (```food_demand_daily```)
- Australian Electricity Generation (```nem_generation_by_source```), taken from *Panagiotelis, A., Gamakumara, P., Athanasopoulos, G., and Hyndman, R. J. (2023). Probabilistic forecast reconciliation: Properties, evaluation and score optimisation. European Journal of Operational Research, 306(2):693–706.*
- Prison Population (```prison_population```), taken from *Hyndman, R. J., & Athanasopoulos, G. (2018). Forecasting: principles and practice. OTexts.*
- Nightly Visitors in Australia (```VNdata```), taken from *Wickramasuriya, S. L., Athanasopoulos, G., and Hyndman, R. J. (2019). Optimal forecast reconciliation for hierarchical and grouped time series through trace minimization. Journal of the American Statistical Association, 114(526):804–819.* and *Girolimetto, D., Athanasopoulos, G., Di Fonzo, T., and Hyndman, R. J. (2023). Cross-temporal probabilistic forecast reconciliation: Methodological and practical issues. International Journal of Forecasting.*

## License

This package is free and open source software, licensed under GPL-3.
