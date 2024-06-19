est_mapping <- c("bu_bu"="Bottom-Up",
                 "var_full"="Full Cov.",
                 "var_spectral"="Spectral",
                 "var_ols"="OLS",
                 "var_scaling"="Var. Scaling")

# Quarterly production of woollen yarn in Australia data from forecast package
wool <- forecast::woolyrnq
wool <- head(wool, n=-3) # remove incomplete last year

wool_dates <- strftime(seq(as.Date("1965-03-01"),as.Date("1993-12-01"),by="3 months"), format="%Y-%m")

wool_aggs <- thief::tsaggregates(wool)
plot(wool_aggs)

wool_df <- data.frame(index=seq_along(wool),
                      x=as.numeric(wool))

wool_agg_df <- rev(sapply(c(1,2,3), function(i){
  tsa <- wool_aggs[[i]]
  if (frequency(tsa)>100){
    tsa_decomp <- stl(tsa, s.window="periodic", robust=TRUE)
    tsa_deseason <- as.numeric(tsa_decomp$time.series[,2]+tsa_decomp$time.series[,3])
  } else {
    tsa_deseason <- tsa
  }
  data.frame(index=rev(rev(wool_df$index)[seq(1, by=4/frequency(tsa), length.out=length(tsa))]),
             x=tsa_deseason - mean(tsa_deseason))
}, simplify=FALSE))

model <- list(fit = fit_agg_arima)
model2 <- list(fit = fit_agg_arima,
               params = list(order= c(list(c(3,1,4)),
                                      list(c(3,1,3)),
                                      list(c(3,1,2))))
)

res <- forecast_and_reconcile(x~1,
                              x=wool_agg_df,
                              k=c(1,2,4),
                              model=model,
                              var_type=c("var", "bu"),
                              cov_type=c("full","bu", "spectral", "ols"),
                              lambda = seq(0, 1, by = 0.1),
                              neig = 1:7)
res2 <- forecast_and_reconcile(x~1,
                               x=wool_agg_df,
                               k=c(1,2,4),
                               model=model2,
                               var_type=c("var", "bu"),
                               cov_type=c("full","bu", "spectral", "ols"),
                               lambda = seq(0, 1, by = 0.1),
                               neig = 1:7)
summary(res)
summary(res2)
x <- .Last.value

kab_df <- x %>%
  mutate(
    level = factor(
      level,
      levels = c(1:(length(unique(x$level)) - 1), "overall"),
      labels = c(paste("Level", 1:(length(unique(x$level)) - 1)), "Overall")
    ),
    type = paste(var_type, cov_type, sep = "_"),
    type = factor(type, levels = names(est_mapping), labels = est_mapping),
    across(c(train_hat_err,test_hat_err,train_tilde_err,test_tilde_err), function(x)x / 1e4),
    across(is.numeric, round, 2)
  ) %>%
  select(level, train_hat_err, test_hat_err, type, contains("err")) %>%
  rename("Training Base MSE"=train_hat_err,
         "Training Recon. MSE"=train_tilde_err,
         "Test Base MSE"=test_hat_err,
         "Test Recon. MSE"=test_tilde_err,
         "Training rMSE"=train_rel_err,
         "Test rMSE"=test_rel_err)

kab_df %>%
  arrange(level,type) %>%
  rename("Level"=level,
         "Recon. Type"=type) %>%
  knitr::kable(format="latex",
               booktabs=TRUE,
               align="lrrlrrrrr",
               caption="Results for the Wool Production Time Series in a Quarter-Biannual-Annual Hierarchy in units of $10^4$ tonnes and fixed order models.",
               label = "data:wool2",
               escape = FALSE) %>%
  #kableExtra::group_rows(index=lcount) %>%
  kableExtra::collapse_rows(columns=1:3, latex_hline = "custom", custom_latex_hline = 1) %>%
  kableExtra::kable_styling(latex_options = "scale_down") -> kab_tex

kab_tex <- gsub("\\raggedright\\arraybackslash", "", kab_tex, fixed=TRUE)
kab_tex <- gsub("\\raggedleft\\arraybackslash", "", kab_tex, fixed=TRUE)

print(kab_tex)


p <- plot(res, facet_fmla = level~type, filter= c("var_type"=list(c("bu","var")), "cov_type"=list(c("bu","full"))))
p2 <- p +
  scale_linetype_manual(values=c(1,1,1)) +
  guides(linetype="none")+
  facet_grid(level~type,
             scales="free",
             labeller=labeller(level=function(x)c("Annual","Biannual","Quarterly")[as.numeric(x)]))+
  labs(x="Year-Month", y="Wool produced in tonnes (de-meaned)")+
  scale_x_continuous(breaks=seq_along(wool_dates)[seq(1,length(wool_dates),length.out=10)],
                     labels=wool_dates[seq(1,length(wool_dates),length.out=10)])+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))

print(p2)
