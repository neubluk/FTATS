library(parallel)
library("tidyverse")

# grid for AR(1) simulations
grid <- expand_grid(Nrep = 50,
                    phi = seq(-0.9,0.9,by=.1),
                    sigma.sq = c(1),
                    h = 1:2,
                    n = c(20,50,100),
                    k = list(c(1,4),c(1,12),c(1,4,12), c(1,5))) %>%
  rowwise() %>%
  mutate(yid = list(1:Nrep)) %>%
  unnest(c(yid)) %>%
  expand_grid(test_recon = "non_cum",
              auto = c(TRUE, FALSE),
              full_cov = c(FALSE))

# grid for random ARMA simulations

# grid <- expand_grid(Nrep = 1:100,
#                p = 1:2,
#                q = 0:2,
#                sigma.sq = 1,
#                h=1:2,
#                n=c(20, 50, 100),
#                k=list(c(1,4),c(1,4,12),c(1,5))) %>%
#   rowwise() %>%
#   mutate(phi = list(generate_ar_params(p)),
#          theta = list(generate_ar_params(q)),
#          yid = list(1:20)) %>%
#   unnest(yid) %>%
#   expand_grid(test_recon = "non_cum",
#               auto = c(FALSE, TRUE),
#               full_cov = c(FALSE))

cl <- makeCluster(5)

clusterEvalQ(cl,{
  library("tidyverse")
  library("FTATS")
  model <- list(fit = fit_agg_arima,
                params = list(
                  mean=FALSE
                ))

  err_fn <- function(y, yhat){
    mean(abs(y-yhat))/mean(abs(y))
  }
})

clusterExport(cl, "grid")

result <-
  #parApply(cl, grid, 1, function(r){
  parSapply(cl = cl, 1:10, simplify = FALSE, FUN = function(i) {
    #apply(grid[1:10,], 1, function(r){
    #sapply(1:nrow(grid[1:10,]), simplify=FALSE, FUN=function(i){
    r <- grid[i,]

    k <- r[["k"]][[1]]
    h <- r[["h"]]

    p <- length(r[["phi"]][[1]])
    q <- length(r[["theta"]][[1]])

    test_recon <- r[["test_recon"]]
    full_cov <- r[["full_cov"]]
    auto <- r[["auto"]]

    # y <- r[["y"]]
    st <- Sys.time()
    y <- FTATS:::sim_agg_arma(n = as.numeric(r[["n"]]),
                      k = k,
                      phi = r[["phi"]][[1]],
                      theta = r[["theta"]][[1]],
                      sigma.sq = as.numeric(r[["sigma.sq"]]))

    # if (full_cov) cov_th <- agg_ar_cov(0.8, 5)
    # else
    cov_th <- NULL

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

    result <- tryCatch(forecast_and_reconcile(
      formula=x~1,
      x=y,
      k=k,
      model=model,
      var_type=c("var", "bu"),
      cov_type=c("ols",
                 "scaling",
                 "full",
                 "cov_shrink",
                 #"autocov_scaling",
                 #"markov_scaling",
                 #"inv_glasso",
                 "corr_shrink",
                 "spectral",
                 "bu"),
      h=h,
      model_train_prop = 0.75,
      lambda = seq(0, 1, by = 0.1),
      neig = 1:sum(k),
      rho = c(0, 0.01, 0.1, 0.25, 1, 5),
      test_recon = test_recon,
      full_cov = cov_th,
      err_measure = err_fn
    ),
    error=function(e)return(NULL))

    et <- Sys.time()

    if (!is.null(result)){
      par_res <- tibble(result$result) %>%
        unnest_wider(1) %>%
        unnest_wider(error) %>%
        unnest(c(level,contains("err"))) %>%
        select(-recon) %>%
        pivot_longer(contains("err"))
    } else {
      par_res <- NULL
    }

    cat(sprintf("rep %d of %d finished in %1.0f s. pars: %s \n",
                i,
                nrow(grid),
                as.numeric(et-st, units="secs"),
                paste(names(r),r, sep="=",collapse=", ")))

    return(par_res)
  })

stopCluster(cl)

grid[1:10, ] %>%
  mutate(res=result) %>%
  unnest(res) %>%
  filter(!full_cov) %>%
  select(Nrep, n, sigma.sq, h, k, auto, any_of(c("p","q", "phi","theta")), yid, var_type, cov_type, level, name, value) -> sim_res



