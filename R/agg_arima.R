#' fit_agg_arima
#'
#' @details
#' Each entry of x must be a data.frame with columns index, corresponding to the observed timestep. This should be aligned across all levels.
#' The formula must align with the columns available in the data.frame, i.e. if response is x with regressor a, then formula = x~a.
#'
#' Fit ARIMA models to multiple aggregation levels of a time series
#' @param formula formula object for modeling the data. if no regressor, use intercept only.
#' @param x list of data.frames, containing variables according to formula
#' @param k vector of aggregation levels
#' @param orders list of ARIMA orders per level of aggregation
#' @param seasonal list of seasonal ARIMA orders per level of aggregation
#' @param m vector of frequencies. required if seasonal order is given
#' @param mean include mean (T/F)
#' @param fixed list of vectors to fix parameters in the models
#' @param ... additional parameters for model function
#' @return object of class agg_arima containing
#' \itemize{
#'  \item models
#'  \item formula
#'  \item k
#'  \item orders
#'  \item fixed
#'  \item xtest_idx
#'  \item data x
#' }
#' @import forecast
#' @export

fit_agg_arima <- function(formula, x, k, orders=rep(NA,length(x)), seasonal=rep(NA,length(x)), m=rep(NA,length(x)), mean=TRUE, fixed=NULL, ...) {
  models <- sapply(seq_along(x), function(i){

    fmla_x <- extract_formula(formula, x[[i]])

    if (all(is.na(orders[[i]])) & all(is.na(seasonal[[i]]))) { # "auto.arima"
      if (!is.null(fixed)) message(sprintf("fixed values are not used for level %d", i))
      if (is.null(fmla_x$regs)){
        regs <- NULL
      } else {
        regs <- as.matrix(fmla_x$regs)
      }
      m <- forecast::auto.arima(fmla_x$response,
                                xreg = regs,
                                allowdrift = FALSE,
                                allowmean = mean,
                                method = "ML",
                                ...)

      m$series <- sprintf("Level %d with %d obsn.", i, nrow(x[[i]]))
      m$call$xreg <- regs
      return(m)

    } else { # regular arima with possible model selection

      stopifnot(length(x) == length(orders))

      stopifnot((length(orders[[i]]) == 3 |
                   (length(orders[[i]]) == 1 &
                      any(is.na(orders[[i]])))) &
                  ((length(seasonal[[i]]) == 3 &
                      !is.na(m[i])) |
                     (length(seasonal[[i]]) == 1 &
                        any(is.na(seasonal[[i]])))))
      fixed2 <- list()
      orders2 <- list()
      seasonal2 <- list()

      fmla_x <- extract_formula(formula, x[[i]])
      y <- fmla_x$response

      if (is.na(orders[[i]][1])){ #p
        orders2 <- list(0:5)
        fixed2 <- list(NA)
      } else {
        orders2 <- list(orders[[i]][1])
        if (orders[[i]][1] >= 1){
          fixed2 <- list(fixed[[i]][1:orders[[i]][1]])
        } else {
          fixed2 <- list(NA)
        }
      }

      if (is.na(orders[[i]][2])) { #d
        #orders2 <- c(orders2,list(0:2) )
        d <- forecast::ndiffs(y)
        orders2 <- c(orders2, list(d))
      }
      else {
        orders2 <- c(orders2,list(orders[[i]][2]))
      }
      fixed2 <- c(fixed2, list(NA))

      if (is.na(orders[[i]][3])) { #q
        orders2 <- c(orders2, list(0:5))
        fixed2 <- c(fixed2, list(NA))
      } else {
        orders2 <- c(orders2,list(orders[[i]][3]))
        if (orders[[i]][3] >= 1){
          fixed2 <- c(fixed2, list(fixed[[i]][orders[[i]][1] +(1:orders[[i]][3])]))
        } else {
          fixed2 <- c(fixed2, list(NA))
        }
      }

      if (is.na(seasonal[[i]][1])){ #P
        seasonal2 <- list(0:5)
        fixed2 <- c(fixed2, list(NA))
      } else {
        seasonal2 <- list(seasonal[[i]][1])
        if (seasonal[[i]][1] >= 1){
          fixed2 <- c(fixed2, list(fixed[[i]][sum(orders[[i]],na.rm=TRUE) + (1:seasonal[[i]][1])]))
        } else {
          fixed2 <- c(fixed2,list(NA))
        }
      }

      if (is.na(seasonal[[i]][2])) { #D
        #seasonal2 <- c(seasonal2,list(0:2) )
        if (!is.na(m[i]) &
            length(y) > 2*m[i]) {
          D <- forecast::nsdiffs(y)
        } else {
          D <- 0
        }
        seasonal2 <- c(seasonal2, list(D))
      }
      else {
        seasonal2 <- c(seasonal2,list(seasonal[[i]][2]))
      }
      fixed2 <- c(fixed2, list(NA))

      if (is.na(seasonal[[i]][3])) { #Q
        seasonal2 <- c(seasonal2, list(0:5))
        fixed2 <- c(fixed2, list(NA))
      } else {
        seasonal2 <- c(seasonal2,list(seasonal[[i]][3]))
        if (seasonal[[i]][3] >= 1){
          fixed2 <- c(fixed2, list(fixed[[i]][sum(orders[[i]],na.rm=TRUE) + seasonal[[i]][1] + (1:seasonal[[i]][3])]))
        } else {
          fixed2 <- c(fixed2, list(NA))
        }
      }

      # if (length(list(...))>0){
      #   if (any(grepl("include.mean",names(list(...)))) & list(...)[[grep("include.mean",names(list(...)))]] ){
      #     fixed2 <- c(fixed2, tail(fixed[[i]],n=1))
      #   }
      # }

      if (mean){
        fixed2 <- c(fixed2, tail(fixed[[i]],n=1))
      }

      if (all(sapply(orders2, length) == 1) &
          all(sapply(seasonal2, length) == 1)) {
        possible_orders <- matrix(unlist(c(orders2, seasonal2)), nrow = 1)
      } else {
        possible_orders <- expand.grid(c(orders2, seasonal2))
        possible_orders <-
          possible_orders[rowSums(possible_orders[, sapply(c(orders2, seasonal2),length) > 1]) <= 5,]


      }

      if (is.na(m[[i]]) |
          m[[i]] == 1)
        possible_orders <-
        possible_orders[rowSums(possible_orders[, 4:6, drop=FALSE]) == 0, , drop = FALSE]

      stopifnot(nrow(possible_orders)>0)

      possible_models <- apply(possible_orders, 1, function(o) {
        f <- NULL
        if (!is.null(fixed)) {
          if (o[1] > 0) {
            if (all(is.na(fixed2[[1]]))) {
              f <- rep(NA, o[1])
            }
            else if (o[1] >= 1) {
              f <- fixed2[[1]]
            }
          }
          if (o[3] > 0) {
            if (all(is.na(fixed2[[3]]))) {
              f <- c(f, rep(NA, o[3]))
            }
            else if (o[3] >= 1) {
              f <- c(f, fixed2[[3]])
            }
          }
          if (o[4] > 0) {
            if (all(is.na(fixed2[[4]]))) {
              f <- c(f, rep(NA, o[4]))
            }
            else if (o[4] >= 1) {
              f <- c(f, fixed2[[4]])
            }
          }
          if (o[6] > 0) {
            if (all(is.na(fixed2[[6]]))) {
              f <- c(f, rep(NA, o[6]))
            }
            else if (o[6] >= 1) {
              f <- c(f, fixed2[[6]])
            }
          }
          if (length(fixed2) > 6 ){
            f <- c(f, fixed2[[7]])
          }
        }


        if (is.null(fmla_x$regs)){
          regs <- NULL
        } else {
          regs <- as.matrix(fmla_x$regs)
        }

        m <- tryCatch(
          forecast::Arima(
            y = y,
            xreg = regs,
            order = o[1:3],
            fixed = f,
            seasonal = list(order = o[4:6], period = m[i]),
            include.mean=mean,
            include.drift=FALSE,
            method = "ML",
            ...
          ),
          error = function(e) {
            print(e)
            return(NULL)
          }
        )
        m$series <- sprintf("Level %d with %d obsn.", i, nrow(x[[i]]))
        m$call$xreg <- regs
        return(m)
      })

      aiccs <- lapply(possible_models, function(m) m$aic + (2*length(m$coef)^2 + 2*length(m$coef))/(m$nobs-length(m$coef)-1))

      return(possible_models[[which.min(aiccs)]])

    }
  }, simplify=FALSE)

  # return list of models with class agg_arima
  res <- list(models=models,
              formula=formula,
              k=k,
              orders=lapply(models,forecast::arimaorder),
              fixed=fixed,
              x=x)
  attr(res, "class") <- "agg_arima"
  return(res)
}

#' Print aggregated ARIMA models
#'
#' @param x result of fit_agg_arima function
#' @param ... additional parameters to print of individual ARIMA models
#' @method print agg_arima
#' @export
print.agg_arima <- function(x, ...){
  lapply(x$models, function(m){
    print(m, ...)
    cat("---------------------------------------------\n")
  })
  invisible(x)
}

#' Update aggregated ARIMA models with new data

#' @param object result of fit_agg_arima
#' @param newdata list of newdata, with same format used for initially fitting the models. It is also possible to only update parts of the aggregation hierarchy (in an on-line setting)
#' @param fix what to fix in the models (models or only the orders)
#' @return updated agg_arima object
#' @import forecast
#' @method update agg_arima
#' @export

update.agg_arima <-
  function(object, newdata, fix = c("model", "order")) {
    fix <- match.arg(fix)

    object$x <- Map(rbind, object$x, newdata)
    if (fix == "model") {
      object$models <- sapply(seq_along(object$models), function(i) {
        fmla_x <- extract_formula(object$formula, object$x[[i]])
        if (is.null(fmla_x$regs)){
          regs <- NULL
        } else {
          regs <- as.matrix(fmla_x$regs)
        }
        m <- forecast::Arima(fmla_x$response, model = object$models[[i]], xreg = regs)
        m$series <- sprintf("Level %d with %d obsn.", i, nrow(object$x[[i]]))
        m$call$xreg <- regs
        return(m)
      }, simplify = FALSE)
    } else if (fix == "order") {
      object <-
        fit_agg_arima(
          formula = object$formula,
          x = object$x,
          k = object$k,
          orders = object$orders,
          fixed = object$fixed,
          lambda = object$models[[1]]$lambda
        )
    }

    return(object)
  }

#' Predict/Forecast agg. ARIMA models
#'
#' @param object of class agg_arima
#' @param h forecast horizon based on top level of the aggregation hierarchy
#' @param partial allow partially observed hierarchies (TRUE some lower level forecasts may be actual observations, FALSE replace partial observations with fitted values)
#' @param newxreg list of regressors if initially used to fit the models
#' @return list of data.frames with variables pred and obs, indicating if the forecasts is an actual forecast or not
#' @import forecast
#' @method predict agg_arima
#' @export
predict.agg_arima <- function(object, h=1, partial=TRUE, newxreg=NULL){
  k <- object$k
  max_index <- max(object$x[[1]]$index)
  # diffs <- c(1, diff(max_indices))
  # diffs_mod <- diffs %% k

  # diffs_mod <- c(1,max_indices[-1]-max_indices[1]) %/% c(k[-1]/k[-length(k)],1)
  diffs_mod <- unlist(lapply(object$x, function(xi) nrow(xi[xi$index>max_index, ])))

  # check what the horizon is
  if (length(h) > 1){
    stopifnot(!is.unsorted(h) & all(diff(h)==1))
    h2 <- h[length(h)]
  } else {
    h2 <- h
  }

  fcs <- sapply(seq_along(k), function(i){
    if (partial | i == 1){
      non_fc_values <- tail(object$x[[i]]$x, n = diffs_mod[i])
    } else {
      non_fc_values <- sapply(1:diffs_mod[i], function(j){
        fvals <- fitted(object$models[[i]],h=j)
        fvals[length(fvals)-(diffs_mod[i]-j)]
      })
    }
    newxregg <- if (diffs_mod[i] > 0) {
      newxreg[[i]][-(1:diffs_mod[i]), , drop=FALSE]
    } else {
      newxreg[[i]]
    }
    fc_values <- predict(object$models[[i]],
                         n.ahead = k[i] - diffs_mod[i] + k[i] * (h2 - 1),
                         newxreg=newxregg)$pred
    mat <- matrix(c(non_fc_values,
                    as.numeric(fc_values)),
                  nrow = h2,
                  byrow = TRUE)
    colnames(mat) <- paste(i,1:k[i],sep="-")
    return(mat)
  })
  fcs_mat <- do.call(cbind, fcs)
  rownames(fcs_mat) <- max(object$x[[1]]$index) + (1:h2)*k[length(k)]

  obs <- sapply(seq_along(k), function(i){
    idx <- 1:(k[i]*h2)
    mat <- matrix(ifelse(idx > diffs_mod[i],
                         "fc",
                         ifelse(partial, "obs","fitted")),
                  nrow = h2,
                  byrow = TRUE)
    colnames(mat) <- paste(i,1:k[i],sep="-")
    return(mat)
  })

  obs_mat <- do.call(cbind, obs)
  rownames(obs_mat) <- rownames(fcs_mat)

  if (length(h) == 1){
    fcs_mat <- fcs_mat[nrow(fcs_mat), , drop=FALSE]
    obs_mat <- obs_mat[nrow(obs_mat), , drop=FALSE]
  }
  return(list(pred=fcs_mat,
              type=obs_mat))
}

#' Obtain fitted values for agg. ARIMA models
#'
#' @param object result of fit_agg_arima
#' @param h horizon used to internally compute fitted values
#' @param simplify simplify the computation of fitted values with computing time in mind (if TRUE then only h-step forecasts are computed on each level. if FALSE, then 1:(kh) step forecasts are computed for aggregation k)
#' @return list of data.frames with variables index and x
#' @import forecast
#' @method fitted agg_arima
#' @export
fitted.agg_arima <- function(object, h=1, simplify=FALSE){
  stopifnot(h>=1)
  lapply(seq_along(object$models), function(i){
    if (simplify == TRUE) return(data.frame(index = object$x[[i]]$index, x = as.numeric(fitted(object$models[[i]], h))))
    if (object$k[i] > 1){
      s <- object$k[i]*(h-1) + 1
      fi <- sapply(s:(object$k[i]*h),function(hi){
        f <- fitted(object=object$models[[i]],h=hi)
        f[seq(hi,by=object$k[i],length.out=length(f)/object$k[i])]
      })
      fi <- as.numeric(t(fi))
    } else {
      fi <- as.numeric(fitted(object$models[[i]], object$k[i]*h))
    }
    data.frame(index = object$x[[i]]$index, x = fi)
  })
}
#' Obtain fitted values for agg. ARIMA models in a more efficient way
#'
#' @param object result of fit_agg_arima
#' @param steps number of new data steps steps for each level of the hierarchy. Useful to incorporate on-line data.
#' @param h horizon used to internally compute fitted values (only h=1 allowed so far)
#' @return list of data.frames with variables index and x
#' @import forecast
#' @method fitted agg_arima
#' @export
fitted2 <- function(object, steps, h=1){
  UseMethod("fitted2",object)
}

fitted2.agg_arima <- function(object, steps, h=1){
  stopifnot(h==1)
  fits <- object
  x <- fits$x
  k <- fits$k
  now_fvs <- lapply(x, function(xi) xi[, c("index","x")])
  for (i in seq_along(steps)){
    #if (steps[length(steps)-i+1] == 0) next
    i2 <- length(x)-i+1
    # now_ind <- seq(sind + timestep + steps[length(steps)-i+1], length.out=nrow(x[[1]])-1, by=k[length(k)])
    now_ind <- seq(k[i2] + steps[length(steps)-i+1], length.out=nrow(x[[1]])-1, by=k[i2])
    new_fcs <- sapply(now_ind,
                      function(j){

                        if (steps[length(steps)-i+1] > 0) {
                          obs <- x[[i2]]$x[(j-steps[length(steps)-i+1]+1):j]
                        } else {
                          obs <- numeric(0)
                        }

                        fcs <- predict(forecast::Arima(x[[i2]]$x[1:j],
                                                       model = fits$models[[i2]],
                                                       xreg = x[[i2]]$reg[1:j, , drop=FALSE]),
                                       n.ahead = k[i2] - steps[length(steps) - i + 1],
                                       newxreg = x[[i2]]$reg[(j+1):(j+k[i2] - steps[length(steps) - i + 1]), , drop=FALSE])$pred

                        c(obs,fcs)
                      })

    now_fvs[[i2]]$x[(1:k[i2])] <- NA
    now_fvs[[i2]]$x[-(1:k[i2])] <- as.numeric(new_fcs)
  }

  return(now_fvs)
}

#' update_and_predict
#'
#' Update aggregated ARIMA models with new data and forecast after
#' @param object of class agg_arima
#' @param newdata_mat matrix of new data with number of columns = size of hierarchy
#' @param h forecast horizon based on top level of the hierarchy
#' @param cumulative perform cumulative one-step ahead forecasts or multi-step forecast (T/F)
#' @param fix what to fix in the models (models or only the orders)
#' @param newxreg regressors if initially used in fitting the models
#' @param ... not used
#' @return list of predict data.frames
#' @import forecast
#' @seealso [update.agg_arima()] [predict.agg_arima()]
#' @export
update_and_predict <- function(object, newdata_mat, h=1, cumulative=TRUE, fix=c("model", "order"), newxreg=NULL, ...){
  UseMethod("update_and_predict",object)
}

#' @method update_and_predict agg_arima
#' @export
update_and_predict.agg_arima <- function(object, newdata_mat, h=1, cumulative=TRUE, fix=c("model", "order"), newxreg=NULL, ...){
  fix <- match.arg(fix)
  k <- object$k
  if (cumulative) {
    fcs <- do.call(rbind,sapply(1:nrow(newdata_mat), function(i) {
      new_data <- hierarchy_matrix_to_list(newdata_mat[1:i, , drop=FALSE], object$k, value_name = colnames(object$x[[1]])[2])
      new_data <- sapply(seq_along(new_data), function(j){
        if (!is.null(newxreg[[j]])){
          tibble(new_data[[j]],newxreg[[j]][1:(i*k[j]),])
        } else {
          new_data[[j]]
        }
      }, simplify=FALSE)
      new_fit <-
        update(object, new_data, fix = fix)
      newxregg <- if (!is.null(newxreg)) {
        sapply(seq_along(newxreg), function(j) {
          if (!is.null(newxreg[[j]])) {
            as.matrix(newxreg[[j]][(i * k[j] + 1):((i + h) * k[j]), , drop = FALSE])
          } else {
            NULL
          }
        }, simplify = FALSE)
      } else {
        NULL
      }
      predict(new_fit, h = h, newxreg=newxregg)$pred
    }, simplify=FALSE))
  } else {
    new_data <- hierarchy_matrix_to_list(newdata_mat, object$k, value_name = colnames(object$x[[1]])[2])
    new_data <- sapply(1:length(new_data), function(j){
      if (!is.null(newxreg[[j]])){
        tibble(new_data[[j]],xreg=newxreg[[j]][1:nrow(new_data[[j]]), ])
      } else {
        new_data[[j]]
      }
    })
    new_fit <-
      update(object, new_data, fix = fix)
    newxregg <- sapply(1:length(newxreg), function(i){
      newxreg[[i]][-(1:nrow(new_data[[i]])), ]
    }, simplify=FALSE)
    fcs <- predict(new_fit, h = h, newxreg=newxregg)$pred
  }

  return(fcs)
}
