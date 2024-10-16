library(tidyverse)

ncluster <- 4

var_type <- c("var","bu")
cov_type <- c("bu", "full", "cov_shrink")

grid <- expand_grid(
  Nrep = 1:50,
  n = c(20, 50, 100),
  k = list(c(1, 4), c(1, 12), c(1, 4, 12)),
  p = 0:2,
  q = 0:2,
  auto = c(FALSE, TRUE)
) %>%
  rowwise() %>%
  mutate(
    y = list(
      FTATS:::sim_agg_arma(
        n,
        k,
        phi = FTATS:::generate_ar_params(p),
        theta = FTATS:::generate_ar_params(q)
      )
    ),
    s = list(c(1:(k[length(k)]-1))),
    full_cov = c("adv")
  ) %>%
  unnest(s)

cl <- parallel::makeCluster(ncluster, "FORK", outfile="./log.txt")
parallel::clusterExport(cl, c("var_type","cov_type","grid"))
parallel::clusterEvalQ(cl,{
  library(FTATS)
  library(tidyverse)
  source("reconcile_partly.R")

  model <- list(fit = FTATS::fit_agg_arima,
                params = list()
  )


})

# parallel::clusterExport(cl, c("var_type","cov_type","grid"))
sim <-
  parallel::parApply(cl, grid, 1, function(r) {
  # apply(grid, 1, function(r){
    p <- as.numeric(r["p"])
    q <- as.numeric(r["q"])
    k <- r["k"][[1]]
    s <- as.numeric(r["s"])
    y <- r["y"][[1]]
    auto <- as.logical(r["auto"])
    # y <- FTATS:::sim_agg_arma(50, k, phi=0.8)

    t0 <-
      max(head(y[[1]]$index,n=-1)) # model train index
    y_train <- lapply(y, function(d)
      d[d$index <= t0, ,drop=FALSE])
    y_test <- lapply(y, function(d)
      d[d$index > t0, ,drop=FALSE])

    # y_test[[2]] <- y_test[[2]][1:s, , drop=FALSE]

    if (!auto) {
      orders <- list(c(p,0,q))
      p_old <- p
      q_old <- q
      k2 <- k[-1]/k[-length(k)]
      for (ki in rev(k2)){

        agg_order <-
          FTATS:::get_agg_order(
            p = p_old,
            k = ki,
            d = 0,
            q = q_old,
            P = 0,
            D = 0,
            Q = 0,
            m = 1
          )
        p_agg <- agg_order[1]
        q_agg <- agg_order[3]

        orders <- c(list(c(p_agg,0,q_agg)),
                    orders)

        p_old <- p_agg
        q_old <- q_agg
      }

      model$params$orders <- orders
    }


    # res <- forecast_and_reconcile(x~1, data=y, k=k, model,
    #                        var_type = c("var","bu"),
    #                        cov_type = c("full","spectral","ols","bu"))

    fits <- model$fit(x~1, y_train, k)

    if (r["full_cov"] == "simple") {
      k2 <- k[-1]/k[-length(k)]
      steps <- c(floor(s/k2)[-1], s)
      ytmp <- hierarchy_matrix_reduce_obs(hierarchy_list_to_matrix(y,k), k, steps)
      full_cov <- cov(hierarchy_list_to_matrix(fitted(fits), k), use="na.or.complete")[colnames(ytmp), colnames(ytmp)]
    } else {
      full_cov <- NULL
    }


    parres <- tryCatch({
      recon_res <- reconcile_forecasts_partly(
        fits,
        s,
        var_type = var_type,
        cov_type = cov_type,
        ytest = y_test,
        full_cov = full_cov,
        err_measure = function(y, yhat)
          sqrt(mean((y - yhat) ^ 2)),
        lambda = seq(0,1,by=0.1),
        verbose = FALSE
      )
      recon_res %>%
        #select(-c("hat","test_hat", "tilde", "test_tilde")) %>%
        unnest(errors)
    },
    error = function(e) {
      print(e)
      return(NULL)
    })

    parres
  })

parallel::stopCluster(cl)

tmp <- tibble(grid) %>%
  mutate(sim = sim) %>%
  unnest(sim) %>%
  mutate(type = paste(var_type,cov_type,sep="_"))

saveRDS(tmp, file="./sim_res.rds")
