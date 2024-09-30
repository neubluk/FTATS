library(parallel)
library("tidyverse")
library(FTATS)

# grid for AR(1) simulations
grid <- expand_grid(Nrep = 100,
                    # phi = 0.8,
                    # phi_fit = seq(-0.9, 0.9, by=0.1),
                    phi = replicate(20, FTATS:::generate_ar_params(1), simplify=FALSE),
                    theta = replicate(20, FTATS:::generate_ar_params(1), simplify=FALSE),
                    p_fit = 0:2,
                    q_fit = 0:2,
                    sigma.sq = c(1),
                    h = 1,
                    n = c(100),
                    k = list(c(1,4),c(1,10))) %>%
  rowwise() %>%
  mutate(yid = list(1:Nrep)) %>%
  unnest(c(yid)) %>%
  expand_grid(test_recon = "non_cum",
              auto = c(FALSE),
              full_cov = c(FALSE))

cl <- makeCluster(100, outfile="./simlog.txt")

clusterEvalQ(cl,{
  library("tidyverse")
  library("FTATS")
  model <- list(fit = fit_agg_arima,
                params = list(
                  mean=FALSE
                ))

  err_fn <- function(y, yhat){
    sqrt(mean((y-yhat)^2))
  }
})

clusterExport(cl, "grid")

result <-
  #parApply(cl, grid, 1, function(r){
  parSapply(cl = cl, 1:nrow(grid), simplify = FALSE, FUN = function(i) {
    # apply(grid[1:10,], 1, function(r){
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

      if (!is.null(r[["p_fit"]]) & !is.null(r[["q_fit"]])){
        model$params$orders[[length(model$params$orders)]] <- c(r[["p_fit"]], 0, r[["q_fit"]])
      }

      if (!is.null(r[["phi_fit"]])){
        model$params$fixed <- lapply(orders, function(oi) rep(NA, sum(oi[-2]!=0)))
        model$params$fixed[[length(model$params$fixed)]] <- c(r[["phi_fit"]][[1]])
      }
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
                 #"corr_shrink",
                 #"spectral",
                 "bu"),
      h=h,
      model_train_prop = 0.75,
      lambda = seq(0, 1, by = 0.1),
      #neig = 1:sum(k),
      #rho = c(0, 0.01, 0.1, 0.25, 1, 5),
      test_recon = test_recon,
      full_cov = cov_th,
      err_measure = err_fn
    ),
    error=function(e){
      print(e)
      return(NULL)
    })

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

grid %>%
  mutate(res=result) %>%
  unnest(res) %>%
  filter(!full_cov) %>%
  select(Nrep, n, sigma.sq, h, k, auto, any_of(c("p","p_fit","q","q_fit", "phi","phi_fit", "theta")), yid, var_type, cov_type, level, name, value) -> sim_res

saveRDS(sim_res, file="./sim_misspec.rds")

# tmp <- sim_res %>% rowwise() %>% mutate(k=paste(k,collapse="-"),pq_fit = paste0("p=",p_fit,", q=",q_fit))
#
# plot_bounds <- tmp %>%
#   filter(k=="1-10", grepl("rel_err", name)) %>%
#   group_by(level, name) %>%
#   mutate(value = ifelse(value==0,NA,value),
#          level = factor(level, levels=c(1:(length(k)), "overall")) )%>%
#   summarise(lower_bound = max(-1.1,quantile(value,0.25,na.rm=TRUE)-1.5*(quantile(value,0.75,na.rm=TRUE)-quantile(value,0.25,na.rm=TRUE))),
#             upper_bound = quantile(value,0.75,na.rm=TRUE)+1.5*(quantile(value,0.75,na.rm=TRUE)-quantile(value,0.25,na.rm=TRUE))) %>%
#   group_by(level) %>%
#   summarise(bounds = list(c(min(lower_bound), max(upper_bound))))
#
# plot_lims <- c(NA,NA)
# plot_bounds <- lapply(plot_bounds$bounds, function(x) {
#   scale_y_continuous(limits = c(if(is.na(plot_lims[1])) x[1] else plot_lims[1],
#                                 if(is.na(plot_lims[2])) x[2] else plot_lims[2]))
# })
#
#
# tmp %>%
#   filter(k=="1-10", grepl("rel_err", name)) %>%
#   ggplot(aes(x=pq_fit, y=value, color=interaction(var_type, cov_type))) +
#   geom_boxplot() +
#   facet_grid(level~name, scales="free_y")
#   # ggh4x::facetted_pos_scales(y=plot_bounds)
sim_res %>%
  distinct(phi, theta) %>%
  unnest_wider(everything(), names_sep="_") %>%
  filter(abs(phi_1) > 0.7,
         abs(theta_1) > 0.7) %>%
  slice_sample(n=20) %>%
  rowwise() %>%
  summarise(phi = list(phi_1),
            theta = list(theta_1)) -> one_model

one_model %>%
  unnest_wider(everything(), names_sep="_")

sim_res %>%
  semi_join(one_model) -> sim_res1

sim_res1 %>%
  filter(grepl("hat_err", name), p_fit == 1, q_fit == 1) %>%
  distinct(phi, theta, yid, name, level, value) %>%
  pivot_wider() -> opt_hat_errs

sim_res_new <- sim_res1 %>%
  filter(!grepl("rel_err", name)) %>%
  pivot_wider() %>%
  left_join(opt_hat_errs, by=c("yid","phi","theta", "level")) %>%
  mutate(train_hat_rel_err = train_hat_err.x/train_hat_err.y - 1,
         train_tilde_rel_err = train_tilde_err/train_hat_err.y - 1,
         test_hat_rel_err = test_hat_err.x/test_hat_err.y - 1,
         test_tilde_rel_err = test_tilde_err/test_hat_err.y - 1) %>%
  select(yid, p_fit, q_fit, var_type, cov_type, level, contains("rel_err"))

sim_res_new %>%
  distinct(yid, p_fit, q_fit, level, train_hat_rel_err, test_hat_rel_err) %>%
  rename(train_tilde_rel_err = train_hat_rel_err,
         test_tilde_rel_err = test_hat_rel_err) %>%
  bind_rows(sim_res_new) %>%
  select(-contains("hat_rel_err")) %>%
  pivot_longer(contains("rel_err")) -> sim_res_new2

plot_bounds <- sim_res_new2 %>%
  group_by(level, name) %>%
  mutate(value = ifelse(value==0,NA,value),
         level = factor(level, levels=c(1:2, "overall")))%>%
  summarise(lower_bound = max(-1.1,quantile(value,0.25,na.rm=TRUE)-1.5*(quantile(value,0.75,na.rm=TRUE)-quantile(value,0.25,na.rm=TRUE))),
            upper_bound = quantile(value,0.75,na.rm=TRUE)+1.5*(quantile(value,0.75,na.rm=TRUE)-quantile(value,0.25,na.rm=TRUE))) %>%
  group_by(level) %>%
  summarise(bounds = list(c(min(lower_bound), max(upper_bound))))

plot_bounds <- lapply(plot_bounds$bounds, function(x) {
  scale_y_continuous(limits = c(x[1], x[2]))
})


sim_res_new2 %>%
  filter(cov_type %in% c("bu","full","spectral", NA)) %>%
  mutate(pq_fit = paste0("p=",p_fit, "\nq=",q_fit),
         name = ifelse(grepl("train",name), 1, 2)) %>%
  ggplot(aes(x=pq_fit, y=value, color=interaction(var_type, cov_type,sep="_"), alpha=!is.na(var_type))) +
  geom_boxplot() +
  geom_hline(yintercept=0) +
  facet_grid(level~name, scales="free_y", labeller = labeller(name=c("1"="Training","2"="Test"),
                                                              level=c("1"="Level 1","2"="Level 2","overall"="Overall"))) +
  ggh4x::facetted_pos_scales(y=plot_bounds) +
  labs(x="Bottom Level Model Fits", y="rMSE", color="Recon. Type") +
  scale_color_discrete(breaks = c("bu_bu","var_full","var_spectral"),
                     labels = est_mapping[c("bu_bu","var_full","var_spectral")])+
  theme(legend.position = "bottom")+
  guides(alpha="none")

ggsave("./arma11_misspec.pdf", width=7, height=4)

# all misspecified
sim_misspec_full %>%
  distinct(phi, theta) %>%
  filter(abs(phi) > 0.7,
         abs(theta) > 0.7) %>%
  slice_sample(n=20) -> one_model

sim_misspec_full %>%
  semi_join(one_model) -> sim_res1

sim_res1 %>%
  filter(grepl("hat_err", name), p_fit1 == 1, p_fit2 == 1, q_fit1 == 1, q_fit2 == 1) %>%
  distinct(phi, theta, yid, name, level, value) %>%
  pivot_wider() -> opt_hat_errs

sim_res_new <- sim_res1 %>%
  filter(!grepl("rel_err", name)) %>%
  pivot_wider() %>%
  left_join(opt_hat_errs, by=c("yid","phi","theta", "level")) %>%
  mutate(train_hat_rel_err = train_hat_err.x/train_hat_err.y - 1,
         train_tilde_rel_err = train_tilde_err/train_hat_err.y - 1,
         test_hat_rel_err = test_hat_err.x/test_hat_err.y - 1,
         test_tilde_rel_err = test_tilde_err/test_hat_err.y - 1) %>%
  select(yid, p_fit1, p_fit2, q_fit1, q_fit2, var_type, cov_type, level, contains("rel_err"))

sim_res_new %>%
  distinct(yid, p_fit1, p_fit2, q_fit1, q_fit2, level, train_hat_rel_err, test_hat_rel_err) %>%
  rename(train_tilde_rel_err = train_hat_rel_err,
         test_tilde_rel_err = test_hat_rel_err) %>%
  bind_rows(sim_res_new) %>%
  select(-contains("hat_rel_err")) %>%
  pivot_longer(contains("rel_err")) -> sim_res_new2

plot_bounds <- sim_res_new2 %>%
  filter(!grepl("train",name)) %>%
  group_by(level, name) %>%
  mutate(value = ifelse(value==0,NA,value),
         level = factor(level, levels=c(1:2, "overall")))%>%
  summarise(lower_bound = max(-1.1,quantile(value,0.25,na.rm=TRUE)-1.5*(quantile(value,0.75,na.rm=TRUE)-quantile(value,0.25,na.rm=TRUE))),
            upper_bound = quantile(value,0.75,na.rm=TRUE)+1.5*(quantile(value,0.75,na.rm=TRUE)-quantile(value,0.25,na.rm=TRUE))) %>%
  group_by(level) %>%
  summarise(bounds = list(c(min(lower_bound), max(upper_bound))))


plot_bounds <- lapply(plot_bounds$bounds, function(x) {
  scale_y_continuous(limits = c(x[1], x[2]))
})



sim_res_new2 %>%
  filter(cov_type %in% c("bu","full", NA), !grepl("train",name)) %>%
  mutate(pq_fit1 = paste0("p=",p_fit1, "\nq=",q_fit1),
         pq_fit2 = paste0("p=",p_fit2, "\nq=",q_fit2),
         name = ifelse(grepl("train",name), 1, 2)) %>%
  filter(pq_fit2 %in% c("p=0\nq=0","p=1\nq=1","p=2\nq=2")) %>%
  # group_by(pq_fit1, pq_fit2, level, var_type, cov_type, name) %>%
  # summarise(value = median(value)) %>%
  ggplot(aes(x=pq_fit1, y=value, color=interaction(var_type, cov_type,sep="_"), alpha=!is.na(var_type))) +
  geom_boxplot() +
  geom_hline(yintercept=0) +
  facet_grid(level~pq_fit2, scales="free_y", labeller = labeller(pq_fit2=function(x)sprintf("Bottom Level Model Fit\n%s",x),
                                                              level=c("1"="Level 1","2"="Level 2","overall"="Overall"))) +
  ggh4x::facetted_pos_scales(y=plot_bounds) +
  labs(x="Top Level Model Fits", y="rMSE", color="Recon. Type") +
  scale_color_discrete(breaks = c("bu_bu","var_full","var_spectral"),
                       labels = est_mapping[c("bu_bu","var_full","var_spectral")])+
  theme(legend.position = "bottom", axis.text.x=element_text(size=7))+
  guides(alpha="none")

ggsave("./arma11_misspec2.pdf", width=7, height=4)

