# FTATS

Forecasting Temporally Aggregated Time Series

[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

This R package features the essential functions from [Neubauer and Filzmoser (2024a)](http://arxiv.org/abs/2407.02367). Exported functions are

- `fit_agg_arima`: Fit fixed or automatically selected ARIMA models to a temporal hierarchy.
- `reconcile_forecasts`: Reconcile base forecasts with reconciliation method of choice.
- `reconcile_forecasts.cv`: Time-Series Cross-Validation function reconciliation methods with hyperparameters.
- `forecast_and_reconcile`: Combines fitting ARIMA models and reconciling resulting base forecasts.

An extension in v1.1.0 features the functionality from [Neubauer and Filzmoser (2024b)](http://arxiv.org/). Exported functions are
- `reconcile_forecasts_partly': Reconcile partly observed base forecasts in order to perform forecast updating

## Citation
``` s
@Manual{,
    title = {Forecasting Temporally Aggregated Time Series},
    author = {Lukas Neubauer},
    year = {2024},
    note = {R package version 1.0.0 to 'Rediscovering Bottom-Up: Effective Forecasting in Temporal Hierarchies'},
    url = {https://github.com/neubluk/FTATS},
    doi = {https://doi.org/10.48550/arXiv.2407.02367},
  }
```

## Installation

``` s
# install.packages("remotes")
remotes::install_github("neubluk/FTATS")
```

## Data Examples

Data examples are available in `demo/data_example.R`, `demo/data_example2.R`, `demo/data_examples_all.R`, as well as simulation studies in `demo/simulation.R`

Similar examples for the updating extension are available in `demo/updating`.

## Data

Five datasets are available.

- Daily Food Demand Data from Schrankel GmbH (```food_demand_daily```)
- Australian Electricity Generation (```nem_generation_by_source```), taken from *Panagiotelis, A., Gamakumara, P., Athanasopoulos, G., and Hyndman, R. J. (2023). Probabilistic forecast reconciliation: Properties, evaluation and score optimisation. European Journal of Operational Research, 306(2):693–706.*
- Prison Population (```prison_population```), taken from *Hyndman, R. J., & Athanasopoulos, G. (2018). Forecasting: principles and practice. OTexts.*
- Nightly Visitors in Australia (```VNdata```), taken from *Wickramasuriya, S. L., Athanasopoulos, G., and Hyndman, R. J. (2019). Optimal forecast reconciliation for hierarchical and grouped time series through trace minimization. Journal of the American Statistical Association, 114(526):804–819.* and *Girolimetto, D., Athanasopoulos, G., Di Fonzo, T., and Hyndman, R. J. (2023). Cross-temporal probabilistic forecast reconciliation: Methodological and practical issues. International Journal of Forecasting.*
- Solar Power Generation (```spower```), taken from *Panamtash, H. and Zhou, Q. (2018). Coherent probabilistic solar power forecasting. In 2018 IEEE International Conference on Probabilistic Methods Applied to Power Systems (PMAPS), pages 1–6. IEEE.*

## License

This package is free and open source software, licensed under GPL-3.
