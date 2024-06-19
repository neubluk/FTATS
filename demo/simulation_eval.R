est_mapping <- c("bu_bu"="Bottom-Up",
                 "var_full"="Full Cov.",
                 "var_spectral"="Spectral",
                 "var_ols"="OLS",
                 "var_scaling"="Var. Scaling")

# ar1 ----

sim_res %>%
  filter(cov_type %in% c("full","bu")) %>%
  rowwise() %>%
  mutate(k = paste0("(",paste(k,collapse = ","),")")) %>%
  group_by(name, h, k, auto, sigma.sq) %>%
  group_walk(~{
    title <- paste(
      names(.y),
      unlist(.y),
      sep = "=",
      collapse = ", "
    )

    levels <- unique(.x$level)
    levels <- c(1:(length(levels)-1), "overall")

    dat <- .x %>%
      group_by(phi, level, n) %>%
      mutate(value_bu = rep(value[cov_type=="bu"],2),
             levels = factor(level, levels=levels)) %>%
      # group_by(phi, var_type, cov_type, level, n, sigma.sq) %>%
      # summarise(mean_value = mean(value-value_bu),
      #           med_value = median(value-value_bu),
      #           se_value = sd(value-value_bu)/sqrt(n()),
      #           mad_value = mad(value-value_bu)) %>%
      # group_by(phi, level, n, sigma.sq) %>%
      # mutate(diff_med_value = med_value - med_value[cov_type=="bu"]) %>%
      filter(cov_type != "bu") %>%
      drop_na()

    level_labeller <- setNames(c(paste("Level", levels[-length(levels)]), "Overall"), levels)

    if (nrow(dat) == 0) {
      p <- ggplot()
    } else {
      p <- dat %>%
        ggplot(aes(x=phi,
                   y=value-value_bu,
                   #color=interaction(var_type,cov_type),
                   group=interaction(phi)))+
        geom_boxplot()+
        #geom_point()+
        # geom_errorbar(aes(ymin=med_value-mad_value,
        #                   ymax=med_value+mad_value))+
        geom_hline(yintercept=0)+
        facet_grid(level~n, scales="free", labeller = labeller(n = function(x) sprintf("n=%s",x),
                                                               level = level_labeller))+
        labs(x=expression(phi), y="rMSE(full) - rMSE(bu)")+
        theme(legend.position = "bottom")#axis.text.x=element_text(angle=90))
      #ggtitle(title)

      if (max(dat$value-dat$value_bu) > 2 & min(dat$value-dat$value_bu) < -2){
        p <- p + ylim(c(-0.5,0.5))
      } else {
        if (max(dat$value-dat$value_bu) > 2) p <- p + ylim(c(NA,0.5))
        else p <- p + ylim(c(-0.5,NA))
      }
    }

    print(p)
  })


sim_res %>%
  filter(cov_type %in% c("bu","full","spectral","ols"), grepl("rel", name)) %>%
  mutate(phi = cut(phi, breaks=c(-0.9, -0.5, 0.5, 0.9), include.lowest = TRUE),
         train = c("1train","2test")[grepl("test",name)+1]) %>%
  # rowwise() %>%
  # mutate(k = paste0("(",paste(k,collapse = ","),")")) %>%
  group_by(n, sigma.sq, h, k, auto, phi, var_type, cov_type, level, train) %>%
  summarise(mean=mean(value),
            se = sd(value)/sqrt(n())) %>%
  group_by(h,k,sigma.sq,auto) %>%
  group_walk(~{
    .y2 <- .y %>%
      rowwise() %>%
      mutate(k = paste0("(",paste(k,collapse = ","),")"))

    title <- paste(
      names(.y2),
      unlist(.y2),
      sep = "=",
      collapse = ", "
    )

    x2 <- .x %>%
      mutate(level = factor(level, levels = c(1:(length(unique(.x$level))-1), "overall"),
                            labels = c(paste("Level",1:(length(unique(.x$level))-1)), "Overall")
      ),
      type = paste(var_type, cov_type, sep="_"),
      type = factor(type, levels=names(est_mapping), labels=est_mapping))

    dat <- x2 %>% mutate(
      mean = ifelse(is.na(mean), "-", format(round(mean,2),nsmall=2)),
      value = ifelse(is.na(se),
                     mean,
                     paste0(mean," (",format(round(se,2),nsmall=2), ")"))
    ) %>%
      select(level, n, type, value, phi, train) %>%
      pivot_wider(names_from=c(train, phi), values_from=value, names_sort = TRUE) %>%
      arrange(level, n, type)

    dat_min_ind <- x2 %>%
      group_by(level, n, train, phi) %>%
      reframe(type,
              is_min = mean == min(mean, na.rm=TRUE)) %>%
      pivot_wider(names_from=c(train, phi), values_from=is_min, names_sort = TRUE) %>%
      arrange(level, n, type) %>%
      mutate(level=FALSE,
             n=FALSE,
             type=FALSE,
             across(where(is.logical), replace_na, FALSE)) %>%
      as.matrix()

    lcount <- dat %>% group_by(level) %>% summarise(n=n_distinct(n, type)) %>% pull(n)
    names(lcount) <- levels(dat$level)

    # if (grepl("train",.y$name)){
    #   caption <- "Training mean rMSE per buckets of $\\phi$ for"
    # } else {
    #   caption <- "Test mean rMSE per buckets of $\\phi$ for"
    # }

    caption <- "Mean rMAE per buckets of $\\phi$ for"

    math_cap <- c("h"="h","k"="k","sigma.sq"="\\sigma^2")
    seps <- c("=", "\\in", "=")

    .y2 <- .y %>%
      select(-any_of(c("auto" ,"name"))) %>%
      rowwise() %>%
      mutate(k = paste0("\\{",paste(rev(k),collapse = ","),"\\}"))

    caption <- paste(caption,
                     "$",
                     paste0(
                       math_cap[names(.y2)],
                       seps,
                       unlist(.y2),
                       collapse = ", "
                     ),
                     "$",
                     collapse=" ")

    if (.y$auto){
      caption <- paste(caption,
                       "and auto-selected models.")
    } else {
      caption <- paste(caption,
                       "and fixed order of the used models.")
    }

    caption <- paste(caption,
                     "The standard errors are given in parenthesis.")

    dat <- as.matrix(dat)

    for (i in 1:ncol(dat_min_ind)){
      dat[, i] <- kableExtra::cell_spec(dat[,i],format="latex", bold=dat_min_ind[,i])
    }

    dat %>%
      knitr::kable(format="latex",
                   booktabs=TRUE,
                   align="lllrrrrrr",
                   caption=caption,
                   label = title,
                   escape = FALSE,
                   col.names = c("Level", "n","Recon. Type", rep(levels(.x$phi),2))) %>%
      #kableExtra::group_rows(index=lcount) %>%
      kableExtra::collapse_rows(columns=1:2, latex_hline = "custom", custom_latex_hline = 1:2) %>%
      kableExtra::add_header_above(c(" "=3, "Training rMAE"=3, "Test rMAE"=3)) %>%
      kableExtra::kable_styling(latex_options = "scale_down") -> kab_tex

    kab_tex <- gsub("\\raggedright\\arraybackslash", "", kab_tex, fixed=TRUE)

    print(kab_tex)
  })
