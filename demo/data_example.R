est_mapping <- c("bu_bu"="Bottom-Up",
                 "var_full"="Full Cov.",
                 "var_spectral"="Spectral",
                 "var_ols"="OLS",
                 "var_scaling"="Var. Scaling")
# dataset from Athanasopoulos et al. (2017) - temporal hierarchy og paper
library(thief)
library(tidyverse)
str(AEdemand)
head(AEdemand)
totals <- AEdemand[ ,"Total Attendances"] #
#totals <- AEdemand[ ,"Total Emergency Admissions via A&E"]
totals_date <- seq(as.Date("2010-11-08"), length.out=length(totals), by="1 week")
complete_years <- totals_date >= as.Date("2011-01-01") & totals_date < as.Date("2015-01-01")
totals <- ts(totals[complete_years],
             start = c(2011,1),
             end = c(2014,52),
             frequency=52)
totals_date <- totals_date[complete_years]

totals_df <- data.frame(index=seq_along(totals),
                        x=as.numeric(totals))

totals_agg <- tsaggregates(totals)

plot(totals_agg, series = c(1,2,3))

totals_agg_df <- rev(sapply(c(1,3), function(i){
  tsa <- totals_agg[[i]]
  tsa_decomp <- stl(tsa, s.window="periodic", robust=TRUE)
  tsa_deseason <- as.numeric(tsa_decomp$time.series[,2]+tsa_decomp$time.series[,3])
  data.frame(index=rev(rev(totals_df$index)[seq(1, by=52/frequency(tsa), length.out=length(tsa))]),
             x=tsa_deseason - mean(tsa_deseason))
}, simplify=FALSE))

class(totals_agg_df) <- "agg_ts"
attr(totals_agg_df, "k") <- c(1,4)

p <- plot(totals_agg_df)

p + geom_vline(xintercept=max(totals_agg_df[[2]]$index)*0.8, color="red")

model <- list(fit = fit_agg_arima,
              params = list(order= c(list(c(1,1,2)),
                                     list(c(1,1,1))))
)

res <- forecast_and_reconcile(x~1,
                              x=totals_agg_df,
                              k=c(1,4),
                              model=model,
                              var_type=c("var", "bu"),
                              cov_type=c("full","bu", "spectral", "ols"),
                              lambda = seq(0, 1, by = 0.1),
                              neig = 1:5)

x <- summary(res)

kab_df <- x %>%
  mutate(
    level = factor(
      level,
      levels = c(1:(length(unique(x$level)) - 1), "overall"),
      labels = c(paste("Level", 1:(length(unique(x$level)) - 1)), "Overall")
    ),
    type = paste(var_type, cov_type, sep = "_"),
    type = factor(type, levels = names(est_mapping), labels = est_mapping),
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
               caption="Results for the Total Addendances Time Series in a Month-Week Hierarchy",
               label = "data:ae:ta2",
               escape = FALSE) %>%
  #kableExtra::group_rows(index=lcount) %>%
  kableExtra::collapse_rows(columns=1:3, latex_hline = "custom", custom_latex_hline = 1) %>%
  kableExtra::kable_styling(latex_options = "scale_down") -> kab_tex

kab_tex <- gsub("\\raggedright\\arraybackslash", "", kab_tex, fixed=TRUE)
kab_tex <- gsub("\\raggedleft\\arraybackslash", "", kab_tex, fixed=TRUE)

print(kab_tex)

p <- plot(res, facet_fmla = level~type, filter= c("var_type"=list(c("bu","var")), "cov_type"=list(c("bu","full"))))
p2 <- p +
  scale_linetype_manual(values=c(1,1)) +
  guides(linetype="none")+
  facet_grid(level~type,
             scales="free",
             labeller=labeller(level=function(x)c("Monthly","Weekly")[as.numeric(x)]))+
  labs(x="Week Date", y="Total Attendances (transformed)")+
  scale_x_continuous(breaks=seq_along(totals_date)[seq(1,length(totals_date),length.out=10)],
                     labels=totals_date[seq(1,length(totals_date),length.out=10)])+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))

print(p2)
