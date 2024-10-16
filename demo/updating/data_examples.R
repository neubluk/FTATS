
library(FTATS)
library(tidyverse)
source("demo/updating/misc.R")

# energy generation ----

data("nem_generation_by_source")

res_energy <- nem_generation_by_source %>%
  mutate(week_date = as.Date(cut(date, "weeks"))) %>%
  group_by(Source, week_date) %>%
  mutate(n = n()) %>%
  filter(week_date >= min(week_date[n == 7]) &
           week_date <= max(week_date[n == 7])) %>%
  group_by(Source) %>%
  summarise(y = list({
    y <- ts(Generation[1:336], frequency = 28)

    y_agg <- thief::tsaggregates(y, aggregatelist = c(1, 7, 28))

    y_df <- data.frame(index = seq_along(y),
                       x = as.numeric(y))

    y_agg_df <- rev(sapply(c(1, 2, 3), function(i) {
      tsa <- y_agg[[i]]
      data.frame(index = rev(rev(y_df$index)[seq(1,
                                                 by = 28 / frequency(tsa),
                                                 length.out = length(tsa))]),
                 x = tsa)
    }, simplify = FALSE))
    y_agg_df
  })) %>%
  expand_grid(s = 1:27) %>%
  group_by(Source, s) %>%
  rowwise() %>%
  mutate(res = list({
    t0 <-
      max(head(y[[1]]$index, n = -1)) # model train index
    y_train <- lapply(y, function(d)
      d[d$index <= t0,])
    y_test <- lapply(y, function(d)
      d[d$index > t0,])

    model <- list(fit = fit_agg_arima)

    fits <- fit_agg_arima(x ~ 1, y_train, c(1, 4, 28))

    parres <- tryCatch({
      recon_res <-
        reconcile_forecasts_partly(
          fits,
          cur_group()$s,
          var_type = c("var", "bu"),
          cov_type = c("full", "bu", "cov_shrink","scaling"),
          ytest = y_test,
          fix = "order",
          err_measure = function(y, yhat=0) sqrt(mean( (y-yhat)^2, na.rm=TRUE)),
          lambda=seq(0,1,by=0.1),
          verbose = TRUE
        )
      recon_res %>%
        # select(-c("hat","test_hat", "tilde", "test_tilde")) %>%
        unnest(errors)
    },
    error = function(e) {
      print(e)
      return(NULL)
    })

    parres
  }))

p <- plot_data_result(res_energy, c(1,4,28), steps=seq(3,27, by=4), est_mapping = c("bu_bu"="Bottom-Up","var_cov_shrink"="Cov. Shrinkage"))
# Coal is the bad outlier
p + xlab("Additional Data in Days")

p <- plot_sample_ts(res_energy %>% filter(Source == "non-Renewable"),
                    k = c(1,4,28),
                    idx_shown = 4,
                    steps=c(7,14,21),
                    xvar="index",
                    est_mapping = c("bu_bu"="Bottom-Up","var_cov_shrink"="Cov. Shrinkage"))

p +
  labs(y="MW", x="", color="Type") +
  theme(legend.position = "bottom", axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_continuous(name = "",
                     labels=unique(nem_generation_by_source$date)[tail(res_energy$y[[1]][[1]]$index,n=4)],
                     breaks=tail(res_energy$y[[1]][[1]]$index,n=4),
                     limits=range(tail(res_energy$y[[1]][[1]]$index,n=4)))
# food demand ----

data("food_demand_daily")
source("demo/updating/fit_agg_arima2.R")
res_food <- food_demand_daily %>%
  filter(!(day %in% c("saturday","sunday"))) %>%
  mutate(week_date = as.Date(cut(date,"weeks"))) %>%
  group_by(fridge_id, week_date) %>%
  mutate(n=n()) %>%
  filter(week_date >= min(week_date[n==5]) & week_date <= max(week_date[n==5]) & week_date >= as.Date("2022-01-01")) %>%
  group_by(fridge_id) %>%
  filter(length(sold)>25)%>%
  summarise(y = list({
    y <- ts(sold, frequency = 5)

    y_agg <- thief::tsaggregates(y)

    y_df <- data.frame(index = seq_along(y),
                       x = as.numeric(y))

    y_agg_df <- rev(sapply(c(1, 2), function(i) {
      tsa <- y_agg[[i]]
      data.frame(index = rev(rev(y_df$index)[seq(1,
                                                 by = 5 / frequency(tsa),
                                                 length.out = length(tsa))]),
                 x = tsa)
    }, simplify = FALSE))
    y_agg_df
  })) %>%
  expand_grid(s = 1:4) %>%
  group_by(fridge_id, s) %>%
  rowwise() %>%
  mutate(res = list({
    t0 <-
      max(head(y[[1]]$index, n = -1)) # model train index
    y_train <- lapply(y, function(d)
      d[d$index <= t0,])
    y_test <- lapply(y, function(d)
      d[d$index > t0,])

    model <- list(
      fit = function(formula,
                     x,
                     k,
                     orders = rep(NA, length(x)),
                     seasonal = rep(NA, length(x)),
                     m = rep(NA, length(x)),
                     mean = TRUE,
                     fixed = NULL,
                     ...) {
        fit_agg_arima2(
          formula,
          x,
          k,
          orders,
          seasonal,
          m,
          mean,
          fixed,
          trans = c(list(list(
            x = function(y)
              log(y + 1)
          )), list(list(
            x = function(y)
              log(y + 1)
          ))),
          trans_back = c(list(list(
            x = function(y)
              exp(y) - 1
          )), list(list(
            x = function(y)
              exp(y) - 1
          ))),
          ...
        )
      }
      # params = list(order= c(list(c(1,1,2)),
      #                        list(c(1,1,1))))
    )

    fits <- model$fit(x ~ 1, y_train, c(1, 5))
    parres <- tryCatch({
      recon_res <-
        reconcile_forecasts_partly(
          fits,
          cur_group()$s,
          var_type = c("var", "bu"),
          cov_type = c("full", "bu", "scaling"),
          ytest = y_test,
          fix = "order",
          err_measure = function(y, yhat=0) sqrt(mean( (y-yhat)^2, na.rm=TRUE)),
          sntz = TRUE,
          x = y_train,
          verbose=TRUE
        )
      recon_res %>%
        # select(-c("hat","test_hat", "tilde", "test_tilde")) %>%
        unnest(errors)
    },
    error = function(e) {
      print(e)
      return(NULL)
    })

    parres
  }))

plot_data_result(res_food, c(1,5))
plot_sample_ts(res_food, c(1,5), 7)
