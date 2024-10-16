#' Australian Electricity Generation
#'
#' Daily time series of electricity generation in Australia per source.
#'
#' @format
#' A tibble of size 8,418 × 3 containing
#' \itemize{
#' \item date
#' \item Source
#' \item Generation
#' }
#' @references Panagiotelis, A., Gamakumara, P., Athanasopoulos, G., and Hyndman, R. J. (2023). Probabilistic forecast reconciliation: Properties, evaluation and score optimisation. European Journal of Operational Research, 306(2):693–706.
#' @source <https://github.com/PuwasalaG/Probabilistic-Forecast-Reconciliation>, <http://opennem.org.au>
"nem_generation_by_source"

#' Australian Prison Population
#'
#' Quarterly time series of prison population in Australia per various categories.
#'
#' @format
#' A tibble of size 3,072 x 6 containing
#' \itemize{
#' \item Date
#' \item State
#' \item Gender
#' \item Legal
#' \item Indigenous
#' \item Count
#' }
#' @references Hyndman, R. J., & Athanasopoulos, G. (2018). Forecasting: principles and practice. OTexts.
#' @source <https://otexts.com/fpp3/extrafiles/prison_population.csv>
"prison_population"

#' Australian Nightly Visitors
#'
#' Monthly time series of nightly visitors in various regions in Australia.
#'
#' @format
#' A time series of size 228 x 525 containing with each column representing a region in Australia.
#'
#' @references Wickramasuriya, S. L., Athanasopoulos, G., and Hyndman, R. J. (2019). Optimal forecast reconciliation for hierarchical and grouped time series through trace minimization. Journal of the American Statistical Association, 114(526):804–819.
#' @references  Girolimetto, D., Athanasopoulos, G., Di Fonzo, T., and Hyndman, R. J. (2023). Cross-temporal probabilistic forecast reconciliation: Methodological and practical issues. International Journal of Forecasting.
#' @source <https://robjhyndman.com/publications/mint/>, <https://github.com/daniGiro/ctprob>
"VNdata"

#' Daily Food Demand Dataset
#'
#' Daily time series of smart fridges, containing delivered and sold items.
#'
#' @format
#' A tibble of size 40,464 x 5 containing
#' \itemize{
#' \item fridge_id
#' \item date
#' \item day
#' \item delivered
#' \item sold
#' }
"food_demand_daily"

#' Solar Power Generation
#'
#' Simulated photovoltaic generation data per state and capacity in MW, on a 5-min basis for the year 2006.
#'
#' @format
#' A data.frame of size 8,198,502 x 4 containing
#' \itemize{
#' \item State
#' \item CapacityMW
#' \item LocalTime
#' \item Power
#' }
#' @source <https://www.nrel.gov/grid/solar-power-data.html>
"spower"
