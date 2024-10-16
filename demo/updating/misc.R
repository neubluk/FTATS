plot_data_result <- function(res, k, steps = NULL, est_mapping=c("bu_bu"="Bottom-Up",
                                                               "var_full"="Full Cov.",
                                                               "var_spectral"="Spectral",
                                                               "var_cov_shrink"="Cov. Shrink",
                                                                "var_corr_shrink"="Corr. Shrink",
                                                                "var_scaling"="Var. Scaling",
                                                               "no_recon"="No Recon."),
                             plot_lims = c(NA,NA)){
  res <- res %>%
    unnest(res)

  if (!is.null(steps)){
    res <- res %>%
      filter(s %in% steps)
  }

  if (!is.null(est_mapping)){
    res <- res %>%
      mutate(type = est_mapping[paste(var_type,cov_type,sep="_")]) %>%
      filter(type %in% est_mapping)
  } else {
    res <- res %>%
      mutate(type = paste(var_type,cov_type,sep="_"))
  }


  res <- res %>%
    group_by(y, type, level) %>%
    distinct(train_hat_err, train_tilde_err, test_hat_err, test_tilde_err) %>%
    summarise(s=0,
              train_tilde_now_base_rel_err = train_tilde_err/train_hat_err - 1,
              test_tilde_now_base_rel_err = test_tilde_err/test_hat_err - 1) %>%
    bind_rows(res %>%
                group_by(y, s, type, level) %>%
                distinct(train_hat_err, train_hat_now_err, test_hat_err, test_hat_now_err) %>%
                summarise(train_tilde_now_base_rel_err = train_hat_now_err/train_hat_err - 1,
                          test_tilde_now_base_rel_err = test_hat_now_err/test_hat_err - 1,
                          type = "No Recon.")) %>%
    bind_rows(res)

  plot_bounds <- res %>%
    pivot_longer(c(train_tilde_now_base_rel_err, test_tilde_now_base_rel_err)) %>%
    group_by(level, s, name) %>%
    mutate(value = ifelse(value==0,NA,value),
           level = factor(level, levels=c(1:(length(k)), "overall")) )%>%
    summarise(lower_bound = max(-1.1,quantile(value,0.25,na.rm=TRUE)-1.5*(quantile(value,0.75,na.rm=TRUE)-quantile(value,0.25,na.rm=TRUE))),
              upper_bound = quantile(value,0.75,na.rm=TRUE)+1.5*(quantile(value,0.75,na.rm=TRUE)-quantile(value,0.25,na.rm=TRUE))) %>%
    group_by(level) %>%
    summarise(bounds = list(c(min(lower_bound), max(upper_bound))))

  plot_bounds <- lapply(plot_bounds$bounds, function(x) {
      scale_y_continuous(limits = c(if(is.na(plot_lims[1])) x[1] else plot_lims[1],
                                    if(is.na(plot_lims[2])) x[2] else plot_lims[2]))
  })

  lblr <- labeller(level = function(x) ifelse(x=="overall", "Overall",sprintf("Level %s",x)),
                   name = function(x) ifelse(grepl("train",x), "Training","Test"))

  print(
    res %>%
    pivot_longer(c(train_tilde_now_base_rel_err, test_tilde_now_base_rel_err)) %>%
    mutate(name = factor(name, levels=c("train_tilde_now_base_rel_err","test_tilde_now_base_rel_err"))) %>%
    ggplot(aes(x=factor(s), y=value, color=type))+
    geom_boxplot(outlier.fill = "white", outlier.shape = 21)+
    geom_hline(yintercept=0, linetype="dashed")+
    facet_grid(level~name, scales="free_y", labeller = lblr)+
    ggh4x::facetted_pos_scales(y=plot_bounds)+
    labs(x="New data available", y="rRMSE", color="Type")+
    theme(legend.position = "bottom")
  )
}

plot_sample_ts <-
  function(res,
           k,
           idx_shown,
           steps = NULL,
           xvar = "index",
           est_mapping = c(
             "bu_bu" = "Bottom-Up",
             "var_full" =
               "Full Cov.",
             "var_spectral" =
               "Spectral",
             "var_cov_shrink" =
               "Cov. Shrink",
             "var_corr_shrink" =
               "Corr. Shrink",
             "var_scaling" =
               "Var. Scaling",
             "no_recon" =
               "No Recon."
           ),
           ...) {

  res <- res %>%
    unnest(res)

  if (!is.null(steps)){
    res <- res %>%
      filter(s %in% steps)
  }

  if (!is.null(est_mapping)){
    res <- res %>%
      mutate(type = est_mapping[paste(var_type,cov_type,sep="_")]) %>%
      filter(type %in% est_mapping)
  } else {
    res <- res %>%
      mutate(type = paste(var_type,cov_type,sep="_"))
  }
  sample_ts <- res %>% ungroup() %>% distinct(y) %>% slice_sample(n=1)


  res_ts <- semi_join(res, sample_ts)
  # print(res_ts)
  pdata <-
    res %>%
    semi_join(sample_ts) %>%
    distinct(y, s, type, hat, test_hat, tilde_now, test_tilde_now) %>%
    bind_rows(res %>%
                semi_join(sample_ts) %>%
                distinct(y, type, hat, test_hat, tilde, test_tilde) %>%
                rename(tilde_now=tilde,
                       test_tilde_now=test_tilde) %>%
                mutate(s=0)) %>%
    rowwise() %>%
    mutate(y=list(FTATS:::hierarchy_list_to_matrix(y, k))) %>%
    pivot_longer(c(hat,test_hat,tilde_now,test_tilde_now,y)) %>%
    rowwise() %>%
    mutate(value={
      tryCatch(list(FTATS:::hierarchy_list_to_plot_data(FTATS:::hierarchy_matrix_to_list(value,k))),
               error=function(e)return(NA)
      )}) %>%
    unnest(value) %>%
    mutate(name = ifelse(name=="y","obs",c("tilde","base")[name %in% c("hat","test_hat") + 1]),
           name2 = ifelse(name %in% c("obs","base"), name, type)) %>%
    mutate(name2 = ifelse(name2 == "base", "Base", ifelse(name2 == "obs", "Observed", name2)))

  pdata_obs <- res %>%
    semi_join(sample_ts) %>%
    distinct(s,y,obs_dat) %>%
    unnest(obs_dat) %>%
    rowwise() %>%
    mutate(y_obs = list(
      FTATS:::hierarchy_list_to_plot_data(
        FTATS:::hierarchy_matrix_to_list({
          m <- FTATS:::hierarchy_list_to_matrix(y, k)
          m2 <- matrix(nrow=nrow(m), ncol=ncol(m),dimnames = dimnames(m))
          m2[ ,obs_dat] <- m[,obs_dat,drop=FALSE]
          m2
        },k,drop_na=FALSE)
      )),
      obs_dat = paste(obs_dat,collapse="/")) %>%
    unnest(y_obs)

  leg_colors <- 1:(length(unique(res_ts$type))+2)
  names(leg_colors) <- c("Observed","Base",unique(res_ts$type))

  lblr <- labeller(s = function(x)sprintf("%s Additional Hours",x),
                   level = function(x) sprintf("Level %s",x))

  plot_bounds <- pdata %>%
    filter(get(xvar) >= min(tail(res_ts$y[[1]][[1]]$index,n=idx_shown))) %>%
    group_by(level) %>%
    summarise(b = list(scale_y_continuous(limits=range(x))))

  print(
    ggplot(mapping=aes(x=.data[[xvar]], y=x, color=name2, linetype = name2 %in% c("Observed"))) +
    geom_line(data = pdata[pdata$level == max(pdata$level),]) +
    geom_step(data = pdata[pdata$level != max(pdata$level),], direction="hv")+
    # scale_color_manual(values=c("obs"=1,"base"=2,"Full Cov."=3,"Spectral"=4,"Var. Scaling"=5,  "Bottom-Up"=6))+
    scale_color_manual(values=leg_colors,
                       # labels = leg_colors[leg_colors != "Observed"],
                       limits = names(leg_colors)[-1])+
    geom_point(data = pdata_obs, mapping=aes(x=index, y=x), inherit.aes = FALSE, size=0.5)+
    geom_vline(xintercept=min(as.numeric(rownames(res_ts$test_hat[[1]])))-k[length(k)], linetype="dotted")+
    # geom_vline(xintercept=as.numeric(rownames(res_ts$test_hat[[1]])), linetype="dotted",alpha=0.2)+
    facet_grid(level~s, scales="free_y", labeller = lblr)+
    ggh4x::facetted_pos_scales(y=plot_bounds$b)+
    xlim(range(tail(res_ts$y[[1]][[1]]$index,n=idx_shown)))+
    guides(linetype = "none")
  )
}
