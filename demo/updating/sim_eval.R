tmp <- readRDS("./sim_res.rds")

path <- "sims_basic"
library(tidyverse)

est_mapping <- c("bu_bu"="Bottom-Up",
                 "var_full"="Full Cov.",
                 "var_spectral"="Spectral",
                 "var_cov_shrink"="Cov. Shrink",
                 "no_recon"="No Recon.")

# basic boxplots for each scenario, both train and test errors
plot_type <- "med_line" # or "boxplot"
tmp %>%
  group_by(across(c(Nrep, any_of(c("n", "k", "p", "q", "auto", "sigma.sq", "m")), type, level))) %>%
  distinct(train_hat_err, train_tilde_err, test_hat_err, test_tilde_err) %>%
  summarise(s=0,
            train_tilde_now_base_rel_err = train_tilde_err/train_hat_err - 1,
            test_tilde_now_base_rel_err = test_tilde_err/test_hat_err - 1) %>%
  bind_rows(tmp %>%
              group_by(across(c(Nrep, any_of(c("n", "k", "p", "q", "auto", "sigma.sq", "m")), s, type, level))) %>%
              distinct(train_hat_err, train_hat_now_err, test_hat_err, test_hat_now_err) %>%
              summarise(train_tilde_now_base_rel_err = train_hat_now_err/train_hat_err - 1,
                        test_tilde_now_base_rel_err = test_hat_now_err/test_hat_err - 1,
                        type = "no_recon")) %>%
  bind_rows(tmp %>%
              select(Nrep, any_of(c("n", "k", "p", "q", "auto", "s", "sigma.sq", "m")), type, level, train_tilde_now_base_rel_err, test_tilde_now_base_rel_err)) %>%
  rename(train=train_tilde_now_base_rel_err,
         test=test_tilde_now_base_rel_err) %>%
  # filter(s %in% seq(0,length.out=20, by=15)) %>%
  filter(n == 100, # choose here the filter for appropriate plots
         # level == 2,
          type %in% c("var_cov_shrink","no_recon"),
         p == q) %>%
  mutate(pq = paste0("p=",p,", q=",q)) %>%
  ungroup() %>%
  select(-c(p,q)) %>%
  pivot_longer(c(train,test)) %>%
  group_by( name,
           across(any_of(c("k", "p", "q", "auto", "sigma.sq", "m")))) %>%
  group_walk(~{
     k <- unlist(.y$k)
    .y$k <- paste0(rev(unlist(.y$k)), collapse=",")
    title <- paste(names(.y),
                   .y,
                   collapse=",",
                   sep="=")
    .x <- .x %>%
      mutate(level = factor(level, levels=c(1:(length(k)), "overall")),
             Type=est_mapping[type],
             # name = factor(name, levels=c("train","test"))
             )

    lbler <- labeller(name = function(x) ifelse(x=="train","Training","Test"),
                      level=function(x) ifelse(x=="overall","Overall",sprintf("Level %s",x)),
                      n=function(x)sprintf("n=%s",x),
                      pq=label_value)

    if (plot_type == "boxplot"){
          plot_bounds <- .x %>%
      group_by(level, s, pq, Type) %>% # choose appropriate grouping variables depending on the filter above
      mutate(value = ifelse(value==0,NA,value))%>%
      summarise(lower_bound = max(-1,quantile(value,0.25,na.rm=TRUE)-1.5*(quantile(value,0.75,na.rm=TRUE)-quantile(value,0.25,na.rm=TRUE))),
                upper_bound = min(1.5,quantile(value,0.75,na.rm=TRUE)+1.5*(quantile(value,0.75,na.rm=TRUE)-quantile(value,0.25,na.rm=TRUE)))) %>%
      group_by(level) %>%
      summarise(bounds = list(c(min(lower_bound, na.rm=TRUE), max(upper_bound, na.rm=TRUE))))
       print(plot_bounds %>% unnest_wider(bounds, names_sep="_"))
    plot_bounds <- lapply(plot_bounds$bounds, function(x) scale_y_continuous(limits = c(x[1],x[2])))


    p <-  .x %>%
      ggplot(aes(color=Type,y=value,x=factor(s)))+
      geom_boxplot(outlier.size = 1)+
      facet_grid(level~pq, scales="free", labeller = lbler)+ # choose appropriate faceting variables depending on the filter above
      ggh4x::facetted_pos_scales(y=plot_bounds)+
      geom_hline(yintercept=0, linetype="dashed")+
      theme(legend.position = "bottom")+
      labs(x="New data available", y="rRMSE")

      print(p)

    } else if (plot_type == "med_line") {
      p <- .x %>%
        group_by(Type, level, s, across(any_of(c("n","pq","name")))) %>%
        summarise(value = quantile(value, na.rm=TRUE, 0.5)) %>%
        ggplot(aes(color=Type,y=value,x=s,fill=Type))+
        geom_point(alpha=0.1, size=0.1)+
        # geom_smooth(method="lm", se=FALSE)+
        geom_line()+
        facet_grid(level~pq, scales="free", labeller = lbler)+ # choose appropriate faceting variables depending on the filter above
        # ggh4x::facetted_pos_scales(y=plot_bounds)+
        geom_hline(yintercept=0, linetype="dashed")+
        theme(legend.position = "bottom")+
        labs(x="New data available", y="rRMSE")+
        ylim(c(-1,1))

      print(p)
    }
})
