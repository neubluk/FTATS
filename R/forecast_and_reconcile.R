#' forecast_and_reconcile
#'
#' Fit ARIMA models to each level of temporal hierarchy and reconcile those base forecasts using selected reconciliation methods
#' @details
#' in contrast to reconcile_forecasts, multiple choices for var_type, cov_type are allowed.
#'
#' @inheritParams fit_agg_arima
#' @inheritParams reconcile_forecasts
#' @param model list with fit function and possible params for the fit function
#' @param h forecast horizon based on top level of the hierarchy
#' @param model_train_prop proportion of x to be used for model fitting
#' @param recon_train_prop proportion of x to be used for only reconciliation purposes (must be 0, not used)
#' @param err_measure error measure to be used in possible cross-validation as well as for calculating train/test errors
#' @param test_recon how should reconciliation be handled on the test data? non-cumulatively, ex-ante or ex-post?
#' @param h_test fixed test set size. required if model_train_prop is NULL
#' @param ... additional parameters for reconcile_forecasts()
#' @return object of class fit_recon with
#' \itemize{
#'  \item fits output of fit_agg_arima
#'  \item y training data
#'  \item ytest test data
#'  \item yhat training base forecasts
#'  \item yhat_test test base forecasts
#'  \item k aggregation parameter
#'  \item result named list of recon object for each reconciliation type
#' }
#' @import forecast
#' @import glasso
#' @import Matrix
#' @seealso [fit_agg_arima()], [reconcile_forecasts()]
#' @export
forecast_and_reconcile <- function(formula,
                                   x,
                                   k,
                                   model,
                                   var_type=c("var","struc","series","bu"),
                                   cov_type=c("ols",
                                              "scaling",
                                              "full",
                                              "cov_shrink",
                                              "autocov_scaling",
                                              "markov_scaling",
                                              "inv_glasso",
                                              "corr_shrink",
                                              "spectral",
                                              "bu"),
                                   h = 1,
                                   model_train_prop = 0.8,
                                   recon_train_prop = 0,
                                   err_measure = function(y, yhat=0) (mean( (y-yhat)^2, na.rm = TRUE)),
                                   test_recon = c("non_cum","ex_ante","ex_post"),
                                   h_test = NULL,
                                   ...) {

  stopifnot(recon_train_prop == 0)

  if (is.null(model_train_prop) & is.null(h_test)) {
    stop()
  }
  else {
    if (!is.null(model_train_prop)) stopifnot(model_train_prop > 0 & model_train_prop < 1)
    else stopifnot(h_test > 0)
  }

  var_type <- match.arg(var_type, several.ok = TRUE)
  cov_type <- match.arg(cov_type, several.ok = TRUE)

  test_recon <- match.arg(test_recon)

  if (is.null(h_test)){
    t0 <-
      max(x[[1]][1:(model_train_prop * nrow(x[[1]])), ]$index) # model train index
  } else {
    t0 <- x[[1]][nrow(x[[1]]) - h_test, ]$index
  }
  data_train <- lapply(x, function(d)
    d[d$index <= t0, ])


  fits <-
    tryCatch(
      do.call(model$fit,
              c(list(formula=formula,
                     x = data_train,
                     k = k),
                model$params)),
      error = function(e) {
        #print(e)
        return(NULL)
      }
    )

  if (is.null(fits)) return(NULL)

  if (!is.null(recon_train_prop) & recon_train_prop > 0){
    stopifnot(recon_train_prop > 0 & recon_train_prop + model_train_prop < 1)
    t1 <-
      max(x[[1]][1:( (recon_train_prop + model_train_prop)* nrow(x[[1]])), ]$index) # recon train index
  } else {
    t1 <- t0
  }
  data_recon <-
    lapply(x, function(d)
      d[d$index > t0 & d$index <= t1, ])
  data_test <- lapply(x, function(d)
    d[d$index > t1, ])

  xreg_all <-
    lapply(Map(rbind, data_recon, data_test), function(d){
      fmla_d <- extract_formula(formula, d)
      fmla_d$regs
    })
  xreg_init <-
    sapply(seq_along(xreg_all), function(i){
      x <- head(xreg_all[[i]], n = k[i]*h)
      if (!is.null(x)) as.matrix(x) else NULL
    }, simplify = FALSE)

  hats <- predict(fits, partial = TRUE, newxreg = xreg_init, h=h)$pred
  if (nrow(data_test[[1]]) > 1){
    hats <- rbind(hats,
                  update_and_predict(fits,
                                     newdata_mat = head(hierarchy_list_to_matrix(
                                       Map(rbind, data_recon, data_test), k
                                     ),
                                     n = -h),
                                     cumulative = TRUE,
                                     h = h,
                                     newxreg = xreg_all))
  }

  ytest <- hierarchy_list_to_matrix(data_test, k)

  if (t0 == t1) {
    # do recon. in-sample
    y <- hierarchy_list_to_matrix(data_train, k)
    yhat <- hierarchy_list_to_matrix(fitted(fits, h=h), k)

    yhat[rowSums(is.na(yhat))>0, ] <- NA

    yhat_test <- hats
    if (h > 1) ytest <- tail(ytest, n=-(h-1))
  } else {
    # have specific recon. subset
    y <- hierarchy_list_to_matrix(data_recon, k)
    #yhat <- hats[1:((t1 - t0) / k[length(k)]), ]
    yhat <- hats[rownames(hats) %in% rownames(y), , drop=FALSE]
    y <- y[rownames(y) %in% rownames(yhat), , drop=FALSE]

    yhat_test <- hats[rownames(hats) %in% rownames(ytest), , drop=FALSE]
  }

  types <- expand.grid(var_type, cov_type)

  types <- types[ rowSums(types == "bu") != 1, ]

  result <- apply(types, 1, function(row) {
    names(row) <- NULL
    exec_time <- system.time(
      recon_result <-
        reconcile_forecasts(
          y,
          yhat,
          k = k,
          var_type = row[1],
          cov_type = row[2],
          yhat_test = yhat_test,
          test_recon_cum = test_recon != "non_cum",
          ytest = if (test_recon == "ex_post") ytest else NULL,
          err_measure = err_measure,
          ...
        )
    )

    errors_train <-
      sapply(c(
        list(yhat),
        list(recon_result$yhat_tilde)
      ),
      function(x) {
        yml <- hierarchy_matrix_to_list(y[rowSums(is.na(x))==0, , drop=FALSE], k)
        x <- x[rowSums(is.na(x))==0, , drop=FALSE]
        #x[, !(colnames(x) %in% colnames(ym_now))] <- (ym - fvs_m)[rowSums(is.na(ym - fvs_m))==0, , drop=FALSE][, !(colnames(x) %in% colnames(ym_now))]
        xl <- hierarchy_matrix_to_list(x, k)
        e <- sapply(seq_along(yml), function(i) err_measure(yml[[i]]$x, xl[[i]]$x))
        c(sum(e), e)
      })

    errors_test <-
      sapply(c(
        list(yhat_test),
        list(recon_result$yhat_test_tilde)
      ),
      function(x) {
        yml <- hierarchy_matrix_to_list(ytest[rowSums(is.na(x))==0, , drop=FALSE], k)
        x <- x[rowSums(is.na(x))==0, , drop=FALSE]
        #x[, !(colnames(x) %in% colnames(ym_now))] <- (ym - fvs_m)[rowSums(is.na(ym - fvs_m))==0, , drop=FALSE][, !(colnames(x) %in% colnames(ym_now))]
        xl <- hierarchy_matrix_to_list(x, k)
        e <- sapply(seq_along(yml), function(i) err_measure(yml[[i]]$x, xl[[i]]$x))
        c(sum(e), e)
      })

    errors <- cbind(errors_train, errors_test)

    colnames(errors) <- c("train_hat_err", "train_tilde_err", "test_hat_err", "test_tilde_err")

    errors <- data.frame(level=factor(c("overall",seq_along(k))),
                         errors)

    errors$train_rel_err <- errors[,3]/errors[,2] - 1
    errors$test_rel_err <- errors[,5]/errors[,4] - 1

    list(var_type = row[1],
         cov_type = row[2],
         recon = recon_result,
         error = errors,
         exec_time = exec_time["user.self"])
  })

  names(result) <- paste(types[,1], types[,2], sep="_")

  final_res <- list(fits = fits,
                    y = y,
                    ytest = ytest,
                    yhat = yhat,
                    yhat_test = yhat_test,
                    k = k,
                    result = result)

  attr(final_res, "class") <- "fit_recon"

  return(final_res)
}

#' Print summary of forecast_and_reconcile output. It prints the summary of the fits as well as a selection of errors (which is also returned invisibly)
#'
#' @param object object of class fit_recon
#' @param ... not used
#' @method summary fit_recon
#' @export
summary.fit_recon <- function(object, ...){
  print(object$fits)
  obj_errs <- tibble(object$result) %>%
          unnest_wider(1) %>%
          unnest(error) %>%
          select(-recon)
  print(obj_errs)
  invisible(obj_errs)
}

#' Plot the time series, base and reconciled forecasts
#'
#' @param obj object of class fit_recon
#' @param facet_fmla control the facetting of the plot, e.g. should the reconciled forecasts across all method be displayed in one plot?
#' @param filter list with names var_type, cov_type to filter recon. methods to plot
#' @param est_mapping named vector to give recon. methods proper names
#' @return ggplot2 object
#' @method plot fit_recon
#' @export
plot.fit_recon <- function(obj, facet_fmla = ~ var_type+cov_type, filter=NULL, est_mapping=NULL){
  y_list <- hierarchy_matrix_to_list(rbind(obj$y,
                                           obj$ytest), k = obj$k)

  yhat_list <- hierarchy_matrix_to_list(rbind(obj$yhat,
                                              obj$yhat_test), k = obj$k)

  test_idx <- max(as.numeric(rownames(obj$yhat)))

  data <- hierarchy_list_to_plot_data(yhat_list)

  data_tilde <- tibble(obj$result) %>%
    unnest_wider(1) %>%
    select(var_type,cov_type, recon) %>%
    unnest_wider(recon) %>%
    rowwise() %>%
    summarise(data = {
      ytilde_list <- hierarchy_matrix_to_list(rbind(yhat_tilde,
                                                    yhat_test_tilde), k = obj$k)
      list(cbind(var_type, cov_type, hierarchy_list_to_plot_data(ytilde_list)))
    })

  data_tilde <- do.call(bind_rows,data_tilde) %>%
    mutate(type = paste(var_type, cov_type, sep="_"))

  if (!is.null(est_mapping)){
    data_tilde$type <- est_mapping[data_tilde$type]
  }

  if (!is.null(filter)){
    data_tilde <- data_tilde %>%
      filter(var_type %in% unlist(filter["var_type"]),
             cov_type %in% unlist(filter["cov_type"]))
  }

  p <- invisible(plot(y_list))

  if (any(grepl("type",all.vars(facet_fmla)))){
    p <- p +
      geom_line(data = data[data$level == max(data$level),], aes(color="Base")) +
      geom_step(data = data[data$level != max(data$level),], aes(color="Base"), direction="hv")+
      geom_vline(xintercept = test_idx)
    p <- p +
      geom_line(data = data_tilde[data_tilde$level == max(data_tilde$level),], aes(color="Recon.")) +
      geom_step(data = data_tilde[data_tilde$level != max(data_tilde$level),], aes(color="Recon."),direction="hv") +
      facet_grid(facet_fmla, scales="free")+
      labs(color="Fc. Type")
  } else {
    data <- data %>%
      mutate(type = "base")
    p <- p +
      geom_line(data = data[data$level == max(data$level),], aes(color=type)) +
      geom_step(data = data[data$level != max(data$level),], aes(color=type),direction="hv")+
      geom_vline(xintercept = test_idx)
    p <- p +
      geom_line(data = data_tilde[data_tilde$level == max(data_tilde$level),], aes(color=type)) +
      geom_step(data = data_tilde[data_tilde$level != max(data_tilde$level),], aes(color=type) ,direction="hv") +
      facet_grid(facet_fmla, scales="free")+
      labs(color="Fc. Type")
  }

  p +
    theme(legend.position = "bottom")
}
