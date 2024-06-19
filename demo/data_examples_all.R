library("tidyverse")
library("tsutils")

est_mapping <- c("bu_bu"="Bottom-Up",
                 "var_full"="Full Cov.",
                 "var_spectral"="Spectral",
                 "var_ols"="OLS",
                 "var_scaling"="Var. Scaling")

err_fn <- function(y, yhat){
  mean( (y-yhat)^2)
}

# M3 data
res_M3_14 <- lapply(subset(Mcomp::M3, "quarterly"), function(ts){

  sq <- start(ts$x)[2]
  eq <- end(ts$x)[2]
  y <- c(ts$x, ts$xx)
  if (sq > 1){
    y <- y[-c(1:(5-sq))]
  }

  if (eq < 4 ){
    y <- y[-c((length(y)-eq+1):length(y))]
  }

  y <- ts(y, frequency = 4)
  print(length(y))
  if (length(y)/4 > 0){

    y_agg <- thief::tsaggregates(y, m=4)

    y_df <- data.frame(index=seq_along(y),
                       x=as.numeric(y))

    y_agg_df <- rev(sapply(c(1,3), function(i){
      tsa <- y_agg[[i]]
      # tsa_decomp <- stl(tsa, s.window="periodic", robust=TRUE)
      # tsa_deseason <- as.numeric(tsa_decomp$time.series[,2]+tsa_decomp$time.series[,3])
      data.frame(index=rev(rev(y_df$index)[seq(1, by=4/frequency(tsa), length.out=length(tsa))]),
                 x=tsa - mean(tsa))
    }, simplify=FALSE))

    model <- list(fit = fit_agg_arima
                  # params = list(order= c(list(c(1,1,2)),
                  #                        list(c(1,1,1))))
    )

    res <- tryCatch(forecast_and_reconcile(x~1,
                                           x=y_agg_df,
                                           k=c(1,4),
                                           model=model,
                                           var_type=c("var", "bu"),
                                           cov_type=c("full","bu", "spectral", "ols"),
                                           lambda = seq(0, 1, by = 0.1),
                                           neig = 1:5,
                                           #h_test = 10,
                                           model_train_prop = 0.8,
                                           err_measure=err_fn),
                    error=function(e){
                      print(e)
                      return(NULL)
                    })
  } else{
    res <- NULL
  }

  if (is.null(res)) return(NULL)

  return(tibble(res$result) %>%
           unnest_wider(1) %>%
           unnest(error) %>%
           select(-recon) %>%
           mutate(n = length(y)))
})

## tourism data
data("VNdata")
res_VNdata <- apply(VNdata, 2, function(ts){

  y <- ts(ts, frequency = 12)

  if (length(y)/3 > 0){

    y_agg <- thief::tsaggregates(y, m=12)

    y_df <- data.frame(index=seq_along(y),
                       x=as.numeric(y))

    y_agg_df <- rev(sapply(c(1,3), function(i){
      tsa <- y_agg[[i]]
      data.frame(index=rev(rev(y_df$index)[seq(1, by=12/frequency(tsa), length.out=length(tsa))]),
                 x=tsa - mean(tsa))
    }, simplify=FALSE))

    model <- list(fit = fit_agg_arima
                  # params = list(order= c(list(c(1,1,2)),
                  #                        list(c(1,1,1))))
    )

    res <- tryCatch(forecast_and_reconcile(x~1,
                                           x=y_agg_df,
                                           k=c(1,3),
                                           model=model,
                                           var_type=c("var", "bu"),
                                           cov_type=c("full","bu", "spectral", "ols"),
                                           lambda = seq(0, 1, by = 0.1),
                                           neig = 1:4,
                                           # h_test = 10,
                                           model_train_prop = 0.8,
                                           err_measure=err_fn),
                    error=function(e){
                      print(e)
                      return(NULL)
                    })
  } else {
    res <- NULL
  }
  if (is.null(res)) return(NULL)

  return(tibble(res$result) %>%
           unnest_wider(1) %>%
           unnest(error) %>%
           select(-recon) %>%
           mutate(n = length(y)))
})

# daily demand
data("food_demand_daily")
res_food_demand_daily <- food_demand_daily %>%
  filter(!(day %in% c("saturday","sunday"))) %>%
  mutate(week_date = as.Date(cut(date,"weeks"))) %>%
  group_by(fridge_id, week_date) %>%
  mutate(n=n()) %>%
  filter(week_date >= min(week_date[n==5]) & week_date <= max(week_date[n==5]) & week_date >= as.Date("2022-01-01")) %>%
  group_by(fridge_id) %>%
  filter(length(sold)>25) %>%
  summarise(n=length(sold),
            res = list({
              #sold <- log(sold+1)
              y <- ts(sold, frequency = 5)

              if (length(y)/5 > 0){

                y_agg <- thief::tsaggregates(y, m=5)

                y_df <- data.frame(index=seq_along(y),
                                   x=as.numeric(y))

                y_agg_df <- rev(sapply(c(1,2), function(i){
                  tsa <- y_agg[[i]]
                  data.frame(index=rev(rev(y_df$index)[seq(1, by=5/frequency(tsa), length.out=length(tsa))]),
                             x=tsa - mean(tsa))
                }, simplify=FALSE))

                model <- list(fit = fit_agg_arima
                              # params = list(order= c(list(c(1,1,2)),
                              #                        list(c(1,1,1))))
                )

                res <- tryCatch(forecast_and_reconcile(x~1,
                                                       x=y_agg_df,
                                                       k=c(1,5),
                                                       model=model,
                                                       var_type=c("var", "bu"),
                                                       cov_type=c("full","bu", "spectral", "ols"),
                                                       lambda = seq(0, 1, by = 0.1),
                                                       neig = 1:6,
                                                       #h_test = 10,
                                                       model_train_prop = 0.8,
                                                       err_measure=err_fn),
                                error=function(e){
                                  print(e)
                                  return(NULL)
                                })

              } else {
                res <- NULL
              }

              if (is.null(res)){
                pres <- NULL
              } else {
                tibble(res$result) %>%
                  unnest_wider(1) %>%
                  unnest(error) %>%
                  select(-recon)
              }
            }))

# electricity
data("nem_generation_by_source")
res_energy <- nem_generation_by_source %>%
  mutate(week_date = as.Date(cut(date,"weeks"))) %>%
  group_by(Source, week_date) %>%
  mutate(n=n()) %>%
  filter(week_date >= min(week_date[n==7]) & week_date <= max(week_date[n==7])) %>%
  group_by(Source) %>%
  summarise(n = length(Generation),
            res = list({
              y <- ts(Generation, frequency = 7)

              if (length(y)/7 > 0){

                y_agg <- thief::tsaggregates(y, m=7)

                y_df <- data.frame(index=seq_along(y),
                                   x=as.numeric(y))

                y_agg_df <- rev(sapply(c(1,2), function(i){
                  tsa <- y_agg[[i]]
                  data.frame(index=rev(rev(y_df$index)[seq(1, by=7/frequency(tsa), length.out=length(tsa))]),
                             x=tsa - mean(tsa))
                }, simplify=FALSE))

                model <- list(fit = fit_agg_arima
                              # params = list(order= c(list(c(1,1,2)),
                              #                        list(c(1,1,1))))
                )

                res <- tryCatch(forecast_and_reconcile(x~1,
                                                       x=y_agg_df,
                                                       k=c(1,7),
                                                       model=model,
                                                       var_type=c("var", "bu"),
                                                       cov_type=c("full","bu", "spectral", "ols"),
                                                       lambda = seq(0, 1, by = 0.1),
                                                       neig = 1:8,
                                                       # h_test = 1,
                                                       model_train_prop = 0.8,
                                                       err_measure=err_fn),
                                error=function(e){
                                  print(e)
                                  return(NULL)
                                })
              } else {
                res <- NULL
              }

              if (is.null(res)) sres <-  NULL
              else {
                sres <- tibble(res$result) %>%
                  unnest_wider(1) %>%
                  unnest(error) %>%
                  select(-recon)
              }
              sres
            }))

# prison population
data("prison_population")
res_prison <- prison_population %>%
  group_by(Date, State) %>%
  summarise(Count=sum(Count)) %>%
  group_by(State) %>%
  summarise(n = length(Count),
            res = list({
              y <- ts(Count, frequency = 4)

              if (length(y)/4 > 0){

                y_agg <- thief::tsaggregates(y, m=4)

                y_df <- data.frame(index=seq_along(y),
                                   x=as.numeric(y))

                y_agg_df <- rev(sapply(c(1,3), function(i){
                  tsa <- y_agg[[i]]
                  data.frame(index=rev(rev(y_df$index)[seq(1, by=4/frequency(tsa), length.out=length(tsa))]),
                             x=tsa - mean(tsa))
                }, simplify=FALSE))

                model <- list(fit = fit_agg_arima
                              # params = list(order= c(list(c(1,1,2)),
                              #                        list(c(1,1,1))))
                )

                res <- tryCatch(forecast_and_reconcile(x~1,
                                                       x=y_agg_df,
                                                       k=c(1,4),
                                                       model=model,
                                                       var_type=c("var", "bu"),
                                                       cov_type= c("full","bu", "spectral", "ols"),
                                                       lambda = seq(0, 1, by = 0.1),
                                                       neig = 1:5,
                                                       # h_test = 10,
                                                       model_train_prop = 0.8,
                                                       err_measure=err_fn),
                                error=function(e){
                                  print(e)
                                  return(NULL)
                                })
              } else {
                res <- NULL
              }

              if (is.null(res)) sres <-  NULL
              else {
                sres <- tibble(res$result) %>%
                  unnest_wider(1) %>%
                  unnest(error) %>%
                  select(-recon)
              }
              sres
            }))

# overall table
res_all <- bind_rows(
  do.call(rbind, res_M3_14) %>% mutate(dataset="M3", k=4, N=length(res_M3_14)),
  do.call(rbind,res_VNdata) %>% mutate(dataset="Tourism", k=3, N=length(res_VNdata)),
  res_food_demand_daily %>% unnest(res) %>% mutate(dataset="Food", k=5, N=nrow(res_food_demand_daily)),
  res_energy %>% unnest(res) %>% mutate(dataset="Energy", k=7, N=nrow(res_energy)),
  res_prison %>% unnest(res) %>% mutate(dataset="Prison", k=4, N=nrow(res_prison))
)

res_all_kab <- res_all %>%
  group_by(dataset, var_type, cov_type, level) %>%
  summarise(N=unique(N),
            k=unique(k),
            nt=paste(range(n/unique(k)), collapse="-"),
            nb=paste(range(n), collapse="-"),
            across(contains("rel_err"), list(mean = function(x)mean(x,trim=0.1, na.rm=TRUE),
                                             se = function(x)chemometrics::sd_trim(x[!is.na(x)],trim=0.1)/(1-0.1)/sqrt(sum(!is.na(x)))))) %>%
  filter(level=="overall") %>%
  pivot_longer(contains("rel_err"), names_to=c("name","value_type"), names_pattern="(.*_.*_.*)_(.*)") %>%
  pivot_wider(names_from=value_type) %>%
  mutate(type = est_mapping[paste(var_type, cov_type, sep="_")],
         value = paste0(format(round(mean,2),nsmall=2)," (",format(round(se,2),nsmall=2), ")")) %>%
  ungroup() %>%
  select(dataset, N, k, nt, nb, type, name, value) %>%
  arrange(dataset,type) %>%
  pivot_wider(names_from=c(name,type)) %>%
  mutate(hierarchy = paste0("$\\{",k,",1\\}$")) %>%
  select(dataset, N, nb, nt, hierarchy, contains("train"), contains("test")) %>%
  arrange(dataset)

kab_tex <- res_all_kab %>%
  select(dataset, N, nt, nb, hierarchy) %>%
  knitr::kable(format="latex",
               booktabs=TRUE,
               align="lllll",
               caption="Dataset properties. $N$ denotes the number of total time series in the dataset, and $n_\\text{bottom},n_\\text{top}$ give the range of the available lengths in the hierarchy given by $k$.",
               label = "dat_props",
               escape = FALSE,
               col.names = c("Dataset", "$N$", "$n_{\\text{bottom}}$", "$n_{\\text{top}}$", "$k\\in$")) #%>%  kableExtra::kable_styling(latex_options = "scale_down")

kab_tex <- gsub("\\raggedleft\\arraybackslash", "", kab_tex, fixed=TRUE)
kab_tex <- gsub("\\raggedright\\arraybackslash", "", kab_tex, fixed=TRUE)

print(kab_tex)

res_all_kab %>%
  select(-c(N, nt, nb, hierarchy)) %>%
  knitr::kable(format="latex",
               booktabs=TRUE,
               align="lrrrrrrrr",
               caption="$10\\%$-trimmed overall means for $5$ datasets and selected reconciliation methods. The standard errors are available in parenthesis.",
               label = "dat_res",
               escape = FALSE,
               col.names = c("Dataset",  rep(c("Bottom-Up","Full Cov.","OLS","Spectral"),2))) %>%
  #kableExtra::group_rows(index=lcount) %>%
  # kableExtra::collapse_rows(columns=1, latex_hline = "custom", custom_latex_hline = 1:2) %>%
  kableExtra::add_header_above(c(" "=1, "Training rMSE"=4, "Test rMSE"=4)) %>%
  kableExtra::kable_styling(latex_options = "scale_down") -> kab_tex

kab_tex <- gsub("\\raggedleft\\arraybackslash", "", kab_tex, fixed=TRUE)
kab_tex <- gsub("\\raggedright\\arraybackslash", "", kab_tex, fixed=TRUE)

res_all %>%
  filter(level=="overall") %>%
  select(dataset, var_type, cov_type, train_tilde_err, test_tilde_err) %>%
  pivot_longer(contains("err")) %>%
  mutate(name = factor(ifelse(name=="train_tilde_err","Train","Test"), levels=c("Train","Test"))) %>%
  pivot_wider(names_from=c(var_type,cov_type),values_from=value) %>%
  unnest() %>%
  group_by(dataset, name) %>%
  summarise(nem_test = list({
    tsutils::nemenyi( pick(everything()), plottype="none")[c("means","intervals")]
  })) %>%
  unnest_wider(nem_test) %>%
  rowwise() %>%
  mutate(type=list(est_mapping[names(means)]),
         intervals=list(as.data.frame(t(intervals)))) %>%
  unnest(c(type,means,intervals)) %>%
  group_by(dataset, name) %>%
  mutate(min_mean_right = V2[which(means==min(means))],
         best = min_mean_right > V1) -> res_all_nem

res_all_nem_plots <- res_all_nem %>%
  arrange(dataset) %>%
  group_by(name) %>%
  summarise(p=list({

    if (all(cur_data()$best)) colours <- c(2)
    else colours <- c(1,2)

    cur_data() %>%
      ggplot(aes(x=means, y=tidytext::reorder_within(type, means, list(dataset))))+
      geom_point(aes(color=best))+
      geom_errorbarh(aes(xmin=V1, xmax=V2))+
      tidytext::scale_y_reordered()+
      facet_grid(dataset~., scales="free_y",) +
      scale_color_manual(values=colours)+
      guides(color="none")+
      labs(x="Mean ranks", y=if (cur_group()$name=="Train") "Recon. Type" else "")+
      ggtitle(sprintf("Relative %s Reconciliation Errors",cur_group()$name))+
      theme(plot.title = element_text(size = 10, face = "bold"))
  }))

p <- do.call(gridExtra::grid.arrange, c(res_all_nem_plots$p, nrow=ceiling(nrow(res_all_nem_plots)/2), as.table=FALSE))

print(p)

p <- res_all %>%
  filter(level=="overall") %>% #, cov_type %in% c("full","bu")) %>%
  select(dataset, var_type, cov_type, train_rel_err, test_rel_err) %>%
  pivot_longer(contains("err")) %>%
  mutate(Type = est_mapping[paste(var_type, cov_type, sep="_")],
         name = factor(ifelse(name=="train_rel_err","Train","Test"), levels=c("Train","Test"))) %>%
  group_by(dataset, var_type, cov_type, name) %>%
  mutate(rank = rank(value)/length(value),
         intcept = max(rank[value<=0], na.rm=TRUE)) %>%
  ggplot(aes(x=rank, y=value, color=Type))+
  geom_step()+
  geom_vline(aes(xintercept=intcept, color=Type), linetype="dashed")+
  geom_hline(yintercept=0)+
  facet_grid(name~dataset, scales="free")+
  ggh4x::facetted_pos_scales(
    y = list(name == "Train" ~ scale_y_continuous(limits = c(-0.5,0.5)),
             name == "Test" ~ scale_y_continuous(limits = c(-0.5,0.5)))
  )+
  theme(legend.position = "bottom")+
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), labels = c(0,0.25,0.5,0.75,1))+
  labs(x="Percentiles", y="rMSE")

print(p)
