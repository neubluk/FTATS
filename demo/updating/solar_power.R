library(FTATS)
library(tidyverse)
source("demo/updating/misc.R")
source("demo/updating/fit_agg_arima2.R")

data("spower")

spower_agg <- spower %>%
  group_by(State,CapacityMW) %>%
  filter(month(LocalTime) <= 1, DataType == "Actual") %>%
  drop_na() %>%
  arrange(LocalTime, .by_group = TRUE) %>%
  summarise(n(),
            y = list({
    y <- ts(Power, frequency = 288)

    # print(head(LocalTime))
    # print(tail(LocalTime))

    y_agg <- thief::tsaggregates(y, m = 288, aggregatelist = c(1,12,12*24))

    # print(lapply(y_agg,tail))

    y_df <- data.frame(index = seq_along(y),
                       x = as.numeric(y))

    y_agg_df <- rev(sapply(c(1, 2, 3), function(i) {
      tsa <- y_agg[[i]]
      data.frame(index = rev(rev(y_df$index)[seq(1,
                                                 by = 288 / frequency(tsa),
                                                 length.out = length(tsa))]),
                 # time = rev(rev(LocalTime)[seq(1,
                 #                                by = 288 / frequency(tsa),
                 #                                length.out = length(tsa))]),
                 x = as.numeric(tsa))
    }, simplify = FALSE))

    class(y_agg_df) <- "agg_ts"
    attr(y_agg_df, "k") <- c(1,24,12*24)

    y_agg_df
  }))

# one test rep

y <- spower_agg$y[[1]]
y[[1]]$time <- y[[1]]$time - hours(23) - minutes(55)
y[[2]]$time <- y[[2]]$time - minutes(55)
y2 <- lapply(y, function(yi){
  yi[,c("index","x")]
})
attributes(y2) <- attributes(y)
plot(y2)

y3 <- lapply(y2, function(yi){
  yi$x <- log(yi$x+1)
  return(yi)
})
attributes(y3) <- attributes(y)
plot(y3)


plot(y[[1]]$x, type="l", xlim = c(0,20))
lines(fitted(forecast::auto.arima(y[[1]]$x)), col="blue")
#
plot(y[[2]]$x, type="l", xlim = c(0,1000))
lines(fitted(forecast::auto.arima(y[[2]]$x)), col="blue")
lines(fitted(forecast::auto.arima(y[[2]]$x, xreg=forecast::fourier(ts(y[[2]]$x,frequency = 12*24), K=5))), col="red")

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
      )),
      list(list(
        x = function(y)
          log(y + 1)
      ))),
      trans_back = c(list(list(
        x = function(y)
          exp(y) - 1
      )), list(list(
        x = function(y)
          exp(y) - 1
      )),
      list(list(
        x = function(y)
          exp(y) - 1
      ))),
      ...
    )
  }
    # params = list(order= c(list(c(1,1,2)),
    #                        list(c(1,1,1))))
)

fs <- c(0, 24, 12*24)
k <- c(1,24,288)
y_reg <- lapply(seq_along(y), function(i){
  yi <- y[[i]]
  yi <- as_tibble(yi)
  if (fs[i] > 0) {
    # print(tail(forecast::fourier(ts(yi$x,frequency = fs[i]), K=5), n=5))
    # print(tail(rbind(forecast::fourier(ts(head(yi$x,n=-1),frequency = fs[i]), K=5),
    #             forecast::fourier(ts(head(yi$x,n=-1),frequency = fs[i]), h=1, K=5)), n=5))
    yi$reg <- rbind(forecast::fourier(ts(head(yi$x,n=-k[i]),frequency = fs[i]), K=5),
                    forecast::fourier(ts(head(yi$x,n=-k[i]),frequency = fs[i]), h=k[i], K=5))
  }
  yi
})

fits <- model$fit(
  x ~ reg,
  x = lapply(seq_along(y_reg), function(i) {
    head(y_reg[[i]], n = -c(1,24,288)[i])
  }),
  k = k
)

# fits
# update_and_predict(fits, tail(FTATS:::hierarchy_list_to_matrix(y_reg, c(1,4)),n=2),
#                    newxreg = lapply(seq_along(y_reg), function(i) {
#                      yi <- y_reg[[i]]
#                      tail(yi[,"reg"], n = c(1, 4)[i])
#
                   # }))

# fv <- fitted(fits)
# attributes(fv) <- attributes(y2)
#
# y_resi <- lapply(seq_along(y2), function(i){
#   y2[[i]]$x <- y2[[i]]$x-fv[[i]]$x
#   return(y2[[i]])
# })
# attributes(y_resi) <- attributes(y2)
# plot(y_resi)
# lapply(y_resi, function(ri) sqrt(mean(ri$x^2, na.rm=TRUE)))
#
# forecast_and_reconcile(x~reg, x=y_reg, k=c(1,12), model=model,
#                        var_type = c("bu"), cov_type=c("bu"), h_test = 1,
#                        err_measure = function(y, yhat = 0) sqrt(mean((y - yhat)^2, na.rm = TRUE)))

xregs <- lapply(seq_along(y_reg), function(i) {
  yi <- y_reg[[i]]
  if (is.null(yi$reg)) return(NULL)
  tail(yi[,"reg"], n = k[i])
})
t0 <-
  max(head(y_reg[[1]]$index, n = -1))
y_test <- lapply(y_reg, function(d)
  d[d$index > t0,])

reconcile_forecasts_partly(
  fits,
  12,
  var_type = c("var"),
  cov_type = c("spectral"),
  ytest = y_test,
  fix = "order",
  err_measure = function(y, yhat=0) sqrt(mean( (y-yhat)^2, na.rm=TRUE)),
  xreg = xregs,
  sntz = TRUE,
  lambda=seq(0,1,by=0.1),
  neig=28
)

library(multidplyr)

cl <- new_cluster(4)
cluster_library(cl, c("tidyverse","FTATS"))
cluster_send(cl, {
  source("reconcile_partly.R")
  source("fit_agg_arima2.R")
  fs <- c(0, 24, 12*24)
  k <- c(1,24,288)
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
        )),
        list(list(
          x = function(y)
            log(y + 1)
        ))),
        trans_back = c(list(list(
          x = function(y)
            exp(y) - 1
        )), list(list(
          x = function(y)
            exp(y) - 1
        )),
        list(list(
          x = function(y)
            exp(y) - 1
        ))),
        ...
      )
    }
    # params = list(order= c(list(c(1,1,2)),
    #                        list(c(1,1,1))))
  )
})

res_solar <- spower_agg %>%
  # slice_head(n=1) %>%
  # filter(CapacityMW=="38MW") %>%
  expand_grid(s = seq(12,180,by=12)) %>%
  group_by(CapacityMW, s) %>%
  rowwise() %>%
  partition(cl) %>%
  mutate(res = list({

    y <- lapply(seq_along(y), function(i){
      yi <- y[[i]]
      yi <- as_tibble(yi)
      if (fs[i] > 0) {
        # print(tail(forecast::fourier(ts(yi$x,frequency = fs[i]), K=5), n=5))
        # print(tail(rbind(forecast::fourier(ts(head(yi$x,n=-1),frequency = fs[i]), K=5),
        #             forecast::fourier(ts(head(yi$x,n=-1),frequency = fs[i]), h=1, K=5)), n=5))
        yi$reg <- rbind(forecast::fourier(ts(head(yi$x,n=-k[i]),frequency = fs[i]), K=3),
                        forecast::fourier(ts(head(yi$x,n=-k[i]),frequency = fs[i]), h=k[i], K=3))
      }
      yi
    })

    t0 <-
      max(head(y[[1]]$index, n = -1)) # model train index
    y_train <- lapply(y, function(d)
      d[d$index <= t0,])
    y_test <- lapply(y, function(d)
      d[d$index > t0,])

    fits <- model$fit(x ~ reg, y_train, k)

    print(fits)

    parres <- tryCatch({
      recon_res <-
        reconcile_forecasts_partly(
          fits,
          cur_group()$s,
          var_type = c("var", "bu"),
          cov_type = c("cov_shrink", "bu"),
          ytest = y_test,
          fix = "order",
          err_measure = function(y, yhat=0) sqrt(mean( (y-yhat)^2, na.rm=TRUE)),
          sntz = TRUE,
          lambda=seq(0,1,by=.1),
          x = y_train
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
  })) %>%
  collect()

p <- plot_data_result(res_solar, c(1,24,288),
                      est_mapping = c("bu_bu"="Bottom-Up", "var_cov_shrink"="Cov. Shrinkage"),
                      plot_lims = c(-1,1))

p +
  ylab("rRMSE") +
  scale_x_discrete(name="Additional Data in Hours", labels=0:15, breaks = seq(0,180,by=12)) +
  theme(legend.position = "bottom", axis.text.x=element_text(angle=45, hjust=1))

p <- plot_sample_ts(res_solar %>% filter(s %in% c(96, 144, 180)), c(1,24,288), 4)
p +
  labs(y="MW",color="") +
  scale_x_continuous(name = "", labels=as.Date(unique(spower$LocalTime)[tail(res_solar$y[[1]][[1]]$index,n=4)]), breaks=tail(res_solar$y[[1]][[1]]$index,n=4), limits=range(tail(res_solar$y[[1]][[1]]$index,n=4)))+
  theme(axis.text.x = element_text(angle=45,hjust=1), legend.position = "bottom")

