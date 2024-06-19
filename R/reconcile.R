#' reconcile_forecasts
#'
#' Reconcile base forecasts using selected reconciliation method
#' @details
#' cov_shrink and corr_shrink types need additional lambda parameter between 0 and 1.
#' inv_glasso requires rho>0 parameter.
#' spectral scaling requires lambda and neig parameters.
#'
#' if multiple are given, then cv is performed to find optimal one.
#'
#'
#' @param y data matrix representing the hierarchy of time series
#' @param yhat data matrix of base forecasts with same dimension as y
#' @param k vector of aggregation levels
#' @param var_type type of variance to be used on each aggregation level
#' @param cov_type type of covariance to be used between various aggregation levels
#' @param yhat_test optional base test forecasts with same number of columns as y and yhat
#' @param test_recon_cum should the reonciliation method be re-estimated cumulatively on the test data (T/F)?
#' @param ytest data matrix of test data in order to be able to use test_recon_cum=TRUE
#' @param S summing matrix. if not given, it is computed according to k
#' @param G mapping matrix. if not given, it is estimated by selected reconciliation method
#' @param SG transformation matrix
#' @param full_cov if not given, it is estimated using data y and yhat
#' @param ... additional parameters for selected reconciliation methods
#' @return list with
#' \itemize{
#'  \item yhat_tilde matrix of reconciled forecasts with same dimension as yhat
#'  \item cov_est
#'  \item S
#'  \item G
#'  \item yhat_test_tilde matrix of reconciled test forecasts with same dimension as yhat_test (if given)
#' }
#' @import glasso
#' @import Matrix
#' @seealso [reconcile_forecasts.cv()]
#' @export
reconcile_forecasts <- function(y,
                                yhat,
                                k,
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
                                yhat_test=NULL,
                                test_recon_cum=FALSE,
                                ytest=NULL,
                                S=NULL,
                                G=NULL,
                                SG=NULL,
                                full_cov=NULL,
                                ...){

  #stopifnot(ncol(y) == ncol(yhat) & ncol(y) == nrow(S) & nrow(S) == sum(k))
  params <- list(...)

  if (!is.null(params$err_measure)) {
    err_measure <- params$err_measure
  }
  else {
    err_measure <- formals(reconcile_forecasts.cv)$err_measure
  }

  cov_est <- NULL

  if (is.null(S)){
    S <- construct_S_matrix(k)
  }

  if (is.null(G) & is.null(SG)) {

    cov_type <- match.arg(cov_type)
    var_type <- match.arg(var_type)

    if (is.null(full_cov)) {
      full_cov_est <- cov(y-yhat, use="na.or.complete")
    }
    else {
      full_cov_est <- full_cov
    }

    if (var_type == "bu" |
        cov_type == "bu") {
      G <-
        cbind(matrix(0, ncol = sum(k) - k[length(k)], nrow = k[length(k)]), diag(1, k[length(k)], k[length(k)]))
    }
    else {
      if (var_type == "var") {
        var_diag <- diag(diag(full_cov_est))
      }
      else if (var_type == "struc") {
        var_diag <- diag(as.numeric(S %*% rep(1, ncol(S))))
      }
      else if (var_type == "series") {
        variances <- diag(full_cov_est)
        pool_mat <-
          Matrix::bdiag(Map(`/`, sapply(
            k, matrix, data = 1, nrow = 1
          ), k))
        var_diag <-
          diag(rep(as.numeric(pool_mat %*% variances), times = k))
      }

      if (cov_type == "ols") {
        cov_est <- diag(1, nrow(full_cov_est), ncol(full_cov_est))
      }
      else if (cov_type == "scaling") {
        cov_est <- var_diag
      }
      else if (cov_type == "full") {
        cov_est <- full_cov_est
      }
      else if (cov_type == "cov_shrink") {
        if (!is.null(params$lambda)) {
          if (length(params$lambda) == 1) {
            cov_est <-
              (1 - params$lambda) * full_cov_est + params$lambda * diag(1, ncol(full_cov_est), ncol(full_cov_est))
          } else {
            #do cv
            tune_grid <-
              matrix(
                params$lambda,
                nrow = length(params$lambda),
                dimnames = c(list(NULL), "lambda")
              )
            return(reconcile_forecasts.cv(y = y,
                                          yhat = yhat,
                                          k = k,
                                          var_type = var_type,
                                          cov_type = cov_type,
                                          tune_grid = tune_grid,
                                          yhat_test = yhat_test,
                                          ytest = ytest,
                                          test_recon_cum = test_recon_cum,
                                          full_cov = full_cov,
                                          S=S,
                                          err_measure = err_measure))
          }
        }
      }
      else if (cov_type == "autocov_scaling") {
        ind_mat <-
          as.matrix(Matrix::bdiag(sapply(k, function(ki)
            matrix(TRUE, nrow = ki, ncol = ki))))
        cov_est <- full_cov_est
        cov_est[!ind_mat] <- 0
      }
      else if (cov_type == "markov_scaling") {
        level_inds <- cbind(1 + c(0, k[-length(k)]), cumsum(k))

        if (is.null(params$level_lag1_corrs)) {
          level_lag1_corrs <- apply(level_inds, 1, function(i) {
            if (i[1] == 1)
              return(1)
            cor(do.call(rbind, sapply(1:(i[1] - 1), function(j)
              (y - yhat)[, i[1]:i[2]][, j:(j + 1)], simplify = FALSE)))[1, 2]
          })
        }
        else {
          level_lag1_corrs <- params$level_lag1_corrs
        }

        level_mats <- sapply(seq_along(k), function(i) {
          if (k[i] == 1)
            return(1)
          tmp <- diag(1, k[i], k[i])
          ind <- expand.grid(1:k[i], 1:k[i])
          matrix(level_lag1_corrs[i] ^ abs(apply(ind, 1, diff)), ncol = k[i])
        })

        cov_est <-
          sqrt(var_diag) %*% as.matrix(Matrix::bdiag(level_mats)) %*% sqrt(var_diag)
      }
      else if (cov_type == "inv_glasso") {
        if (!is.null(params$rho)) {
          if (length(params$rho) == 1) {
            glasso_mat <- glasso::glasso(
              cov2cor(full_cov_est),
              rho = params$rho,
              penalize.diagonal = FALSE,
              trace = FALSE,
              approx = FALSE
            )$w

            cov_est <-
              sqrt(var_diag) %*% glasso_mat %*% sqrt(var_diag)
          } else {
            tune_grid <-
              matrix(
                params$rho,
                nrow = length(params$rho),
                dimnames = c(list(NULL), "rho")
              )
            return(reconcile_forecasts.cv(y = y,
                                          yhat = yhat,
                                          k = k,
                                          var_type = var_type,
                                          cov_type = cov_type,
                                          tune_grid = tune_grid,
                                          yhat_test = yhat_test,
                                          ytest = ytest,
                                          test_recon_cum = test_recon_cum,
                                          full_cov = full_cov,
                                          S=S,
                                          err_measure = err_measure))
          }
        }
      }
      else if (cov_type == "corr_shrink") {
        if (!is.null(params$lambda)) {
          if (length(params$lambda) == 1) {
            full_corr_shrink <-
              (1 - params$lambda) * cov2cor(full_cov_est) + params$lambda * diag(1, ncol(full_cov_est))
            cov_est <-
              sqrt(var_diag) %*% full_corr_shrink %*% sqrt(var_diag)
          } else {
            #do cv
            tune_grid <-
              matrix(
                params$lambda,
                nrow = length(params$lambda),
                dimnames = c(list(NULL), "lambda")
              )
            return(reconcile_forecasts.cv(y = y,
                                          yhat = yhat,
                                          k = k,
                                          var_type = var_type,
                                          cov_type = cov_type,
                                          tune_grid = tune_grid,
                                          yhat_test = yhat_test,
                                          ytest = ytest,
                                          test_recon_cum = test_recon_cum,
                                          full_cov = full_cov,
                                          S=S,
                                          err_measure = err_measure))
          }
        }
      }
      else if (cov_type == "spectral") {
        if (!is.null(params$lambda) &
            !is.null(params$neig)) {
          if (length(params$lambda) == 1 &
              length(params$neig) == 1) {
            full_corr_shrink <-
              (1  -  params$lambda) * cov2cor(full_cov_est) + params$lambda * diag(1, ncol(full_cov_est))
            eig_full_corr_shrink <- eigen(full_corr_shrink)
            W <-
              eig_full_corr_shrink$vectors[, 1:params$neig, drop = FALSE]
            sigma2 <-
              mean(eig_full_corr_shrink$values[-(1:params$neig)])
            Q <-
              W %*% diag(
                eig_full_corr_shrink$values[1:params$neig]  -  sigma2,
                ncol = params$neig,
                nrow = params$neig
              ) %*% t(W) + params$lambda *
              diag(sigma2, ncol  =  nrow(W), nrow  =  nrow(W))

            cov_est <- sqrt(var_diag) %*% Q %*% sqrt(var_diag)
          } else {
            tune_grid <-
              as.matrix(expand.grid(lambda = params$lambda, neig = params$neig))
            return(reconcile_forecasts.cv(y = y,
                                          yhat = yhat,
                                          k = k,
                                          var_type = var_type,
                                          cov_type = cov_type,
                                          tune_grid = tune_grid,
                                          yhat_test = yhat_test,
                                          ytest = ytest,
                                          test_recon_cum = test_recon_cum,
                                          full_cov = full_cov,
                                          S=S,
                                          err_measure = err_measure))
          }
        }
      }

      G <-
        tryCatch(
          solve(t(S) %*% solve(cov_est) %*% S) %*% t(S) %*% solve(cov_est),
          error = function(e) {
            #print(e)
            return(NULL)
          }
        )
    }
  }

  if (is.null(G) & is.null(SG)){
    yhat_tilde <- NA
    yhat_test_tilde <- NA
  } else {
    if (is.null(SG)) SG <- S %*% G
    yhat_tilde <- t(SG %*% t(yhat))
    colnames(yhat_tilde) <- colnames(yhat)
    if (!is.null(yhat_test)) {

      if (test_recon_cum == TRUE & nrow(yhat_test) > 1){
        yhat_test_tilde <- t(SG %*% t(yhat_test[1, , drop=FALSE]))
        if (is.null(ytest)) { # ex ante
          y2 <- rbind(y,yhat_test[1, , drop=FALSE])
          ytest2 <- NULL
        } else { # ex post
          y2 <- rbind(y, ytest[1, , drop=FALSE])
          ytest2 <- ytest[-1, , drop=FALSE]
        }
        yhat2 <- rbind(yhat, yhat_test[1, , drop=FALSE])
        yhat_test2 <- yhat_test[-1, , drop=FALSE]

        yhat_test_tilde <- rbind(yhat_test_tilde,
                                 reconcile_forecasts(y = y2,
                                                     yhat = yhat2,
                                                     k = k,
                                                     var_type=var_type,
                                                     cov_type=cov_type,
                                                     yhat_test=yhat_test2,
                                                     test_recon_cum=TRUE,
                                                     ytest=ytest2,
                                                     S=S,
                                                     G=NULL,
                                                     full_cov=NULL,
                                                     ...
                                 )$yhat_test_tilde)
      } else {
        yhat_test_tilde <- t(SG %*% t(yhat_test))
      }

      colnames(yhat_test_tilde) <- colnames(yhat_test)
    }
    dimnames(G) <- if (!is.null(G)) list(colnames(y)[(ncol(y)-k[length(k)]+1):ncol(y)], colnames(y))

    dimnames(cov_est) <- if (!is.null(cov_est)) list(colnames(y), colnames(y))
  }

  dimnames(S) <- if (!is.null(S)) list(colnames(y), colnames(y)[(ncol(y)-k[length(k)]+1):ncol(y)])

  result <- if (is.null(yhat_test)) {
    list(
      yhat_tilde = yhat_tilde,
      cov_est = cov_est,
      S = S,
      G = G)
  } else {
    list(
      yhat_tilde = yhat_tilde,
      yhat_test_tilde = yhat_test_tilde,
      cov_est = cov_est,
      S = S,
      G = G)
  }

  return(result)
}

#' reconcile_forecasts.cv
#'
#' Cross-Validation for selected hyperparameters for recon. methods where necessary
#' @details
#' cov_shrink and corr_shrink types need additional lambda parameter between 0 and 1.
#' inv_glasso requires rho>0 parameter.
#' spectral scaling requires lambda and neig parameters.
#'
#'
#' @inheritParams reconcile_forecasts
#' @param tune_grid matrix with hyperparameters as variables, and combinations in rows
#' @param ... additional parameters for reconcile_forecasts()
#' @param err_measure error function f(x,xtilde) to be minimized and choose optimal hyperparameters. default MSE.
#' @return list of optimal recon result, and cv results
#' @import glasso
#' @import Matrix
#' @seealso [reconcile_forecasts()]
#' @export
reconcile_forecasts.cv <- function(y,
                                   yhat,
                                   k,
                                   var_type=c("var","struc","series"),
                                   cov_type=c("cov_shrink",
                                              "inv_glasso",
                                              "corr_shrink",
                                              "spectral"),
                                   tune_grid,
                                   yhat_test=NULL,
                                   ytest=NULL,
                                   full_cov=NULL,
                                   S=NULL,
                                   err_measure = function(x, xtilde) mean( (x-xtilde)^2),
                                   ...) {

  nfold <- min(nrow(y)/2, 5)

  # rand <- sample(nrow(y))
  # folds <- sapply(1:nfold, function(i) rand[rand %% nfold + 1 == i], simplify=FALSE)

  tfold_size <- ceiling(nrow(y)/2/(nfold+1))
  folds <- sapply(head(seq(ceiling(nrow(y)/2), nrow(y), length.out=nfold+1), n=-1), function(j){
    list(train=1:(j-1), test=j:(j+tfold_size))
  }, simplify=FALSE)

  # tune_grid <- if(!is.matrix(tune_grid)) {
  #   matrix(tune_grid, nrow=length(tune_grid))
  # }

  stopifnot(is.matrix(tune_grid))

  cv_result <- apply(tune_grid, 1, function(par) {
    cv_rel_resi <- sapply(seq_along(folds), function(i) {
      trfi <- folds[[i]]$train
      tefi <- folds[[i]]$test
      tmp_res <- do.call(reconcile_forecasts,
                         c(
                           list(
                             y = y[trfi, , drop=FALSE],
                             yhat = yhat[trfi, , drop=FALSE],
                             k = k,
                             var_type = var_type,
                             cov_type = cov_type,
                             yhat_test = yhat[tefi, , drop=FALSE],
                             test_recon_cum = FALSE,
                             full_cov = full_cov,
                             S=S
                           ),
                           par
                         ))


      if (is.null(tmp_res$G)) return(NA)

      y_fi_l <- hierarchy_matrix_to_list(y[tefi, , drop=FALSE], k = k)
      ytest_fi_l <- hierarchy_matrix_to_list(tmp_res$yhat_test_tilde, k = k)
      # rel_resi_tilde <-
      #   hierarchy_matrix_to_list(y[fi, , drop=FALSE] - tmp_res$yhat_test_tilde, k = k)
      # rel_resi_hat <-
      #   hierarchy_matrix_to_list(y[fi, , drop=FALSE] - yhat[fi, , drop=FALSE], k = k)
      # mean_rel_rmse <- unlist(lapply(rel_resi_tilde, function(r)sqrt(mean(r$x^2)))) /
      # unlist(lapply(rel_resi_hat, function(r)sqrt(mean(r$x^2))))

      mean_errors <- sapply(seq_along(y_fi_l), function(j) err_measure(y_fi_l[[j]]$x, ytest_fi_l[[j]]$x))

      return(mean(mean_errors, na.rm=TRUE))
    })
    return(cv_rel_resi)
  })
  opt_par <- tune_grid[which.min(colMeans(cv_result, na.rm=TRUE)),]
  best_recon_result <- do.call(reconcile_forecasts,
                               c(
                                 list(
                                   y = y,
                                   yhat = yhat,
                                   k = k,
                                   var_type = var_type,
                                   cov_type = cov_type,
                                   yhat_test = yhat_test,
                                   ytest = ytest,
                                   full_cov=full_cov,
                                   S=S,
                                   ...
                                 ),
                                 opt_par
                               ))

  return(c(best_recon_result,
           cv = list(c(
             opt_par, "cv_scores" = list(colMeans(cv_result, na.rm=TRUE))
           ))))
}
