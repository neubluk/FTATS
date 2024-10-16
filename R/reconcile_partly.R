#' reconcile_forecasts_partly
#'
#' Reconcile base forecasts using selected reconciliation method in a partly manner
#' @details
#' This function can be used to perform temporal hierarchical forecast updating by using partial information of the hierarchy.
#'
#' cov_shrink and corr_shrink types need additional lambda parameter between 0 and 1.
#' inv_glasso requires rho>0 parameter.
#' spectral scaling requires lambda and neig parameters.
#'
#' if multiple are given, then cv is performed to find optimal one.
#'
#'
#' @param fits a fitted hierarchical model with fitted, fitted2, update, predict, update_and_predict functions and fields k (and data x)
#' @param bottom_steps number of bottom level steps, must be greater or equal to 1
#' @param var_type type of variance to be used on each aggregation level
#' @param cov_type type of covariance to be used between various aggregation levels
#' @param err_measure error measure function of form f(y,yhat)
#' @param ytest list of test data to compute test errors
#' @param full_cov if not given, it is estimated using fits
#' @param fix for model updating fix either "model" or "order"
#' @param sntz use set-negative-to-zero heuristic (T/F)
#' @param x optional list containining the data used to model the fits
#' @param xreg optional regressors, if not given, should be contained in x (either by fits or given)
#' @param verbose detailed output (T/F)
#' @param ... additional parameters for selected reconciliation methods
#' @return tibble with columns
#' \itemize{
#' \item var_type
#' \item cov_type
#' \item level
#' \item \strong{train_hat_err} base forecast training errors
#' \item \strong{train_tilde_err} reconciled training errors for no additional data
#' \item \strong{train_hat_now_err} updated base forecast training errors
#' \item \strong{train_tilde_now_err} updated reconciled training errors
#' \item \strong{test_hat_err} base forecast test errors
#' \item \strong{test_tilde_err} reconciled test errors for no additional data
#' \item \strong{test_hat_now_err} updated base forecast test errors
#' \item \strong{test_tilde_now_err} updated reconciled test errors
#' \item \strong{train_tilde_now_base_rel_err} relative errors train_tilde_now_err/train_hat_err - 1
#' \item \strong{train_tilde_base_rel_err} relative errors train_tilde_err/train_hat_err - 1
#' \item \strong{train_tilde_now_rel_err} relative errors train_tilde_now_err/train_tilde_err - 1
#' \item \strong{test_tilde_now_base_rel_err} relative errors test_tilde_now_err/test_hat_err - 1
#' \item \strong{test_tilde_rel_err} relative errors test_tilde_err/test_hat_err - 1
#' \item \strong{test_tilde_now_rel_err} relative errors test_tilde_now_err/test_tilde_err - 1
#' }
#' if verbose, then also included list-columns are: hat, test_hat, tilde, test_tilde, hat_now, test_hat_now, tilde_now, test_tilde_now, obs_dat, SG, SG_now
#' @import glasso
#' @import Matrix
#' @seealso [reconcile_forecasts()]
#' @export
reconcile_forecasts_partly <- function(fits,
                                       bottom_steps=1L,
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
                                       err_measure = function(y, yhat=0) mean( (y-yhat)^2, na.rm=TRUE),
                                       ytest = NULL,
                                       full_cov = NULL,
                                       fix = c("model","order"),
                                       sntz = FALSE,
                                       x = NULL,
                                       xreg = NULL,
                                       verbose = FALSE,
                                       ...){

  var_type <- match.arg(var_type, several.ok = TRUE)
  cov_type <- match.arg(cov_type, several.ok = TRUE)

  k <- fits$k
  k2 <- k[-1]/k[-length(k)]
  if (is.null(x)) {
    x <- fits$x
  }

  stopifnot(bottom_steps >= 1 & bottom_steps < k[length(k)])

  steps <- bottom_steps
  if (length(k2) > 1){
    for (i in length(k2):2){
      hie_step <- floor(steps[1]/k2[i])
      steps <- c(hie_step, steps)
    }
  }

  # steps <- c(0,steps)

  # fvs <- fitted(fits, simplify = FALSE)
  fvs <- fitted2(fits, rep(0, length(steps)+1), h=1)

  # now_fvs <- fvs

  xreg <- lapply(x, function(xi) xi$reg)

  if (!is.null(xreg)){
    xreg_init <-
      sapply(seq_along(xreg), function(i){
        x <- head(xreg[[i]], n = k[i])
        if (!is.null(x)) as.matrix(x) else NULL
      }, simplify = FALSE)
  } else {
    xreg_init <- NULL
  }

  test_fc <- predict(fits, h=1, partial=TRUE, newxreg = xreg_init)$pred
  if (nrow(ytest[[1]]) > 1){
    test_fc <- rbind(test_fc,
                     update_and_predict(fits,
                                        newdata_mat = head(FTATS:::hierarchy_list_to_matrix(ytest, k), n = -1),
                                        cumulative = TRUE,
                                        h = 1,
                                        fix,
                                        newxreg = xreg))
  }


  ytest_now <- c(list(NULL),
                 sapply(2:(length(steps) + 1), function(i) {
                   if (steps[i-1] == 0) return(NULL)
                   ytest[[i]][1:steps[i - 1], , drop = FALSE]

                 }, simplify = FALSE))
  # test_fc_now <- predict(update(fits, newdata=ytest_now), h=1, partial=TRUE)$pred

  # sind <- min(x[[1]]$index)
  # timestep <- diff(x[[length(x)]]$index)[1]

  # for (i in seq_along(steps)){
  #   if (steps[length(steps)-i+1] == 0) next
  #   i2 <- length(x)-i+1
  #   # now_ind <- seq(sind + timestep + steps[length(steps)-i+1], length.out=nrow(x[[1]])-1, by=k[length(k)])
  #   now_ind <- seq(k[i2] + steps[length(steps)-i+1], length.out=nrow(x[[1]])-1, by=k[i2])
  #   new_fcs <- sapply(now_ind,
  #                     function(j){
  #                       c(
  #                         x[[i2]]$x[(j-steps[length(steps)-i+1]+1):j],
  #                         predict(forecast::Arima(x[[i2]]$x[1:j],
  #                                                 model = fits$models[[i2]],
  #                                 xreg = x[[i2]]$reg[1:j, , drop=FALSE]),
  #                                 n.ahead = k[i2] - steps[length(steps) - i + 1],
  #                                 newxreg = x[[i2]]$reg[(j+1):(j+k[i2] - steps[length(steps) - i + 1]), , drop=FALSE])$pred
  #                       )
  #                     })
  #
  #   now_fvs[[i2]]$x[-(1:k[i2])] <- as.numeric(new_fcs)
  # }

  now_fvs <- fitted2(fits, c(0,steps), h=1)

  test_fc_now <- sapply(1:nrow(ytest[[1]]), function(i){
    ytest_now_tmp <- c(
      sapply(1:(length(steps) + 1), function(j) {
        st <- ifelse(j==1,0,steps[j-1])
        if (st + (i-1)*k[j] == 0) return(NULL)
        ytest[[j]][1:(st + (i-1)*k[j]), , drop = FALSE]

      }, simplify = FALSE))
    ytest_now_reg <- c(
      sapply(1:(length(steps) + 1), function(j) {
        st <- ifelse(j==1,0,steps[j-1])
        # if (st + (i-1)*k[j] == 0) return(NULL)
        ytest[[j]][1:(i*k[j]), , drop = FALSE]$reg

      }, simplify = FALSE))
    predict(update(fits, newdata=ytest_now_tmp), h=1, partial=TRUE, newxreg=ytest_now_reg)$pred # check xreg in update and predict step
  }, simplify=FALSE)

  test_fc_now <- do.call(rbind,test_fc_now)

  now_fvs_m <- FTATS:::hierarchy_list_to_matrix(now_fvs, k)
  fvs_m <- FTATS:::hierarchy_list_to_matrix(fvs, k)
  ym <- FTATS:::hierarchy_list_to_matrix(x, k)

  ytestm <- FTATS:::hierarchy_list_to_matrix(ytest, k)

  #stopifnot(nrow(ytestm) == 1)

  types <- expand.grid(var_type, cov_type)

  types <- types[ rowSums(types == "bu") != 1, , drop=FALSE ]

  colnames(types) <- c("var_type","cov_type")

  ym_now <- FTATS:::hierarchy_matrix_reduce_obs(ym, k, steps)
  now_fvs_m2 <- FTATS:::hierarchy_matrix_reduce_obs(now_fvs_m, k, steps)
  test_fc_now2 <- FTATS:::hierarchy_matrix_reduce_obs(test_fc_now, k, steps) #[, colnames(ym_now), drop=FALSE]

  G_bu <- cbind(matrix(0, ncol = sum(k) - k[length(k)], nrow = k[length(k)]), diag(1, k[length(k)], k[length(k)]))
  G_bu <- G_bu[-c(1:bottom_steps), which(colnames(ym) %in% colnames(ym_now)), drop=FALSE]

  S_reduced <- FTATS:::construct_S_matrix(k)
  S_reduced <- S_reduced[ ,-c(1:bottom_steps), drop=FALSE]
  S_reduced <- S_reduced[rowSums(S_reduced)>0, , drop=FALSE]

  res <- apply(types, 1, function(r){
    recon_result_0 <- FTATS:::reconcile_forecasts(ym,
                                                  fvs_m,
                                                  k = k,
                                                  var_type = r[1],
                                                  cov_type = r[2],
                                                  # S = S_reduced,
                                                  # G = if(r[1] == "bu") G_bu else NULL,
                                                  # lambda = seq(0, 1, by = 0.1),
                                                  # neig = 1:ncol(ym),
                                                  err_measure = err_measure,
                                                  sntz = sntz,
                                                  ...)

    if(r[1] == "bu") {
      G_now <- G_bu
    } else {
      G_now <- NULL
    }
    recon_result <- FTATS:::reconcile_forecasts(ym_now,
                                                now_fvs_m2,
                                                k = k - c(0,steps),
                                                var_type = r[1],
                                                cov_type = r[2],
                                                S = S_reduced,
                                                G = G_now,
                                                # lambda = seq(0, 1, by = 0.1),
                                                # neig = 1:ncol(ym_now),
                                                # full_cov=full_cov,
                                                err_measure = err_measure,
                                                sntz = FALSE,
                                                ...)


    # training errors
    if (!is.null(recon_result$G)){
      back_tilde_now <- FTATS:::hierarchy_matrix_augment_obs(recon_result$yhat_tilde, k , ym[, !(colnames(ym) %in% colnames(ym_now)), drop=FALSE])

      if (sntz == TRUE){
        back_tilde_now <- FTATS:::set_negative_to_zero(back_tilde_now, k)
      }

      back_tilde_now2 <- back_tilde_now
      if (!is.null(recon_result_0$G)){
        back_tilde_now2[, !(colnames(ym) %in% colnames(ym_now))] <- recon_result_0$yhat_tilde[, !(colnames(ym) %in% colnames(ym_now))]
      }
    } else {
      back_tilde_now <- NA
      back_tilde_now2 <- NA
    }

    now_fvs_m3 <- now_fvs_m
    now_fvs_m3[, !(colnames(ym) %in% colnames(ym_now))] <- fvs_m[, !(colnames(ym) %in% colnames(ym_now))]

    errors_train <-
      sapply(c(
        list(fvs_m),
        list(recon_result_0$yhat_tilde),
        list(now_fvs_m3),
        list(back_tilde_now2)
      ),
      function(x) {
        if (all(is.na(x))) return(rep(NA, length(k)+1))
        yml <- FTATS:::hierarchy_matrix_to_list(ym[rowSums(is.na(x))==0, , drop=FALSE], k)
        x <- x[rowSums(is.na(x))==0, , drop=FALSE]
        #x[, !(colnames(x) %in% colnames(ym_now))] <- (ym - fvs_m)[rowSums(is.na(ym - fvs_m))==0, , drop=FALSE][, !(colnames(x) %in% colnames(ym_now))]
        xl <- FTATS:::hierarchy_matrix_to_list(x, k)
        e <- sapply(seq_along(yml), function(i) err_measure(yml[[i]]$x, xl[[i]]$x))
        c(sum(e), e)
      })

    colnames(errors_train) <- c("train_hat_err","train_tilde_err","train_hat_now_err", "train_tilde_now_err")

    # test errors
    if (!is.null(recon_result_0$G)){
      ytest_tilde <- with(recon_result_0, t(S%*%G%*%t(test_fc)))
      if (sntz == TRUE){
        ytest_tilde <- FTATS:::set_negative_to_zero(ytest_tilde, k)
      }
    } else {
      ytest_tilde <- NA
    }

    if (!is.null(recon_result$G)){
      test_tilde_now <- with(recon_result, t(S%*%G%*%t(test_fc_now2)))
      back_test_tilde_now <- FTATS:::hierarchy_matrix_augment_obs(test_tilde_now, k , ytestm[, !(colnames(ytestm) %in% colnames(ym_now)), drop=FALSE])

      if (sntz == TRUE) {
        back_test_tilde_now <- FTATS:::set_negative_to_zero(back_test_tilde_now, k)
      }

      back_test_tilde_now2 <- back_test_tilde_now
      if (!is.null(recon_result_0$G)){
        back_test_tilde_now2[, !(colnames(ym) %in% colnames(ym_now))] <- ytest_tilde[, !(colnames(ym) %in% colnames(ym_now)), drop=FALSE]
      }

    } else {
      test_tilde_now <- NA
      back_test_tilde_now <- NA
      back_test_tilde_now2 <- NA
    }

    test_fc_now3 <- test_fc_now
    test_fc_now3[, !(colnames(ym) %in% colnames(ym_now))] <- test_fc[, !(colnames(ym) %in% colnames(ym_now)), drop=FALSE]

    errors_test <-
      sapply(c(
        list(test_fc),
        list(ytest_tilde),
        list(test_fc_now3),
        list(back_test_tilde_now2)
      ),
      function(x) {
        if (all(is.na(x))) return(rep(NA, length(k)+1))
        yml <- FTATS:::hierarchy_matrix_to_list(ytestm[rowSums(is.na(x))==0, , drop=FALSE], k)
        x <- x[rowSums(is.na(x))==0, , drop=FALSE]
        #x[, !(colnames(x) %in% colnames(ym_now))] <- (ym - fvs_m)[rowSums(is.na(ym - fvs_m))==0, , drop=FALSE][, !(colnames(x) %in% colnames(ym_now))]
        xl <- FTATS:::hierarchy_matrix_to_list(x, k)
        e <- sapply(seq_along(yml), function(i) err_measure(yml[[i]]$x, xl[[i]]$x))
        c(sum(e), e)
      })

    colnames(errors_test) <-c("test_hat_err", "test_tilde_err", "test_hat_now_err", "test_tilde_now_err")

    errors <- cbind(errors_train, errors_test)

    errors <- data.frame(level=factor(c("overall",seq_along(k))),
                         errors)

    # errors$train_now_rel_err <- errors[,3]/errors[,2] - 1
    # errors$train_now_tilde_rel_err <- errors[,4]/errors[,3] - 1
    # errors$train_tilde_rel_err <- errors[,4]/errors[,2] - 1
    #
    # errors$test_now_rel_err <- errors[,6]/errors[,5] - 1
    # errors$test_now_tilde_rel_err <- errors[,7]/errors[,6] - 1
    # errors$test_tilde_rel_err <- errors[,7]/errors[,5] - 1

    errors$train_tilde_now_base_rel_err <- errors$train_tilde_now_err/errors$train_hat_err - 1
    errors$train_tilde_base_rel_err <- errors$train_tilde_err/errors$train_hat_err - 1
    errors$train_tilde_now_rel_err <- errors$train_tilde_now_err/errors$train_tilde_err - 1

    errors$test_tilde_now_base_rel_err <- errors$test_tilde_now_err/errors$test_hat_err - 1
    errors$test_tilde_rel_err <- errors$test_tilde_err/errors$test_hat_err - 1
    errors$test_tilde_now_rel_err <- errors$test_tilde_now_err/errors$test_tilde_err - 1

    if (verbose) {
        resi <- list(errors=errors,
           hat = fvs_m,
           test_hat = test_fc,
           tilde = recon_result_0$yhat_tilde,
           test_tilde = ytest_tilde,
           hat_now = now_fvs_m,
           test_hat_now = test_fc_now,
           tilde_now = back_tilde_now,
           test_tilde_now = back_test_tilde_now,
           obs_dat = list(colnames(ym)[!(colnames(ym) %in% colnames(ym_now))]),
           SG = tryCatch(recon_result_0$S %*% recon_result_0$G, error=function(e)return(NULL)),
           SG_now = tryCatch(recon_result$S %*% recon_result$G, error=function(e)return(NULL))
         )
    } else {
      resi <- list(errors=errors)
    }

    return(resi)
  })

  fres <- types %>%
    mutate(res = res) %>%
    unnest_wider(res, simplify=FALSE)

  return(fres)
}
