fit_agg_arima2 <- function(formula,
  x,
  k,
  orders = rep(NA, length(x)),
  seasonal = rep(NA, length(x)),
  m = rep(NA, length(x)),
  mean = TRUE,
  fixed = NULL,
  trans = rep(list(NULL), length(x)),
  trans_back = rep(list(NULL), length(x)),
  ...
){
  
  x2 <- sapply(seq_along(x), function(i) {
    if (!is.null(trans[[i]])) {
      as_tibble(sapply(names(x[[i]]), function(v) {
        if (v %in% names(trans[[i]])) {
          trans[[i]][[v]](x[[i]][[v]])
        } else {
          x[[i]][[v]]
        }
      }, simplify = FALSE, USE.NAMES = TRUE))
    } else {
      x[[i]]
    }
  }, simplify=FALSE)
  
  attributes(x2) <- attributes(x)
  
  res <- fit_agg_arima(
    formula,
    x2,
    k,
    orders,
    seasonal,
    m,
    mean,
    fixed,
    ...
  )
  
  res$trans <- trans
  res$trans_back <- trans_back
  
  class(res) <- "agg_arima2"
  
  return(res)
}

print.agg_arima2 <- function(object) {
  class(object) <- "agg_arima"
  FTATS:::print.agg_arima(object)
}

predict.agg_arima2 <- function(object, h=1, partial=TRUE, newxreg=NULL) {
  if (!is.null(newxreg)) {
    newxreg2 <- sapply(seq_along(newxreg), function(i) {
      if (!is.null(object$trans[[i]]) & !is.null(newxreg[[i]]) & !is.null(names(newxreg[[i]]))) {
        data.frame(sapply(names(newxreg[[i]]), function(v) {
          if (v %in% names(object$trans[[i]])) {
            object$trans[[i]][[v]](newxreg[[i]][, v])
          } else {
            newxreg[[i]][, v]
          }
        }, simplify = TRUE),
        row.names = rownames(newxreg[[i]]))
      } else {
        newxreg[[i]]
      }
    }, simplify = FALSE)
  } else {
    newxreg2 <- NULL
  }
  
  class(object) <- "agg_arima"
  res <- FTATS:::predict.agg_arima(object, h=h, partial=partial, newxreg=newxreg2)
  preds <- FTATS:::hierarchy_matrix_to_list(res$pred, k=object$k)
  
  preds <- sapply(seq_along(preds), function(i) {
    if (!is.null(object$trans_back[[i]])) {
      data.frame(sapply(names(preds[[i]]), function(v) {
        if (v %in% names(object$trans_back[[i]])) {
          object$trans_back[[i]][[v]](preds[[i]][, v])
        } else {
          preds[[i]][, v]
        }
      }, simplify = FALSE))
    } else {
      preds[[i]]
    }
  }, simplify=FALSE)
  
  res$pred <- FTATS:::hierarchy_list_to_matrix(preds, k=object$k)
  
  return(res)
}

fitted.agg_arima2 <- function(object, h=1, simplify=FALSE) {
  class(object) <- "agg_arima"
  res <- FTATS:::fitted.agg_arima(object, h=h, simplify=simplify)
  res2 <- res
  
  # res2$x <- FTATS:::hierarchy_matrix_to_list(res2$x, k=object$k)
  res2 <- sapply(seq_along(res2), function(i) {
    if (!is.null(object$trans_back[[i]])) {
      data.frame(sapply(names(res2[[i]]), function(v) {
        if (v %in% names(object$trans_back[[i]])) {
          object$trans_back[[i]][[v]](res2[[i]][, v])
        } else {
          res2[[i]][, v]
        }
      }, simplify = FALSE),
      row.names = rownames(res2[[i]]))
    } else {
      res2[[i]]
    }
  }, simplify=FALSE)
  
  # res2$x <- FTATS:::hierarchy_list_to_matrix(res2$x, k=object$k)
  
  return(res2)
}


fitted2 <- function(object, steps, h=1) {
  UseMethod("fitted2",object)
}
fitted2.agg_arima2 <- function(object, steps, h=1) {
  class(object) <- "agg_arima"
  res <- FTATS:::fitted2.agg_arima(object, steps=steps, h=h)
  res2 <- res
  
  # res2$x <- FTATS:::hierarchy_matrix_to_list(res2$x, k=object$k)
  res2 <- sapply(seq_along(res2), function(i) {
    if (!is.null(object$trans_back[[i]])) {
      data.frame(sapply(names(res2[[i]]), function(v) {
        if (v %in% names(object$trans_back[[i]])) {
          object$trans_back[[i]][[v]](res2[[i]][, v])
        } else {
          res2[[i]][, v]
        }
      }, simplify = FALSE),
      row.names = rownames(res2[[i]]))
    } else {
      res2[[i]]
    }
  }, simplify=FALSE)
  
  # res2$x <- FTATS:::hierarchy_list_to_matrix(res2$x, k=object$k)
  
  return(res2)
}



update.agg_arima2 <-
  function(object, newdata, fix = c("model", "order")) {
    class(object) <- "agg_arima"
    trans <- object$trans
    newdata <- sapply(seq_along(newdata), function(i) {
      if (!is.null(trans[[i]])) {
        as_tibble(sapply(names(newdata[[i]]), function(v) {
          if (v %in% names(trans[[i]])) {
            trans[[i]][[v]](newdata[[i]][[v]])
          } else {
            newdata[[i]][[v]]
          }
        }, simplify = FALSE, USE.NAMES = TRUE))
      } else {
        newdata[[i]]
      }
    }, simplify=FALSE)
    
    res <- FTATS:::update.agg_arima(object, newdata, fix)
    
    class(res) <- "agg_arima2"
    
    return(res)
}

update_and_predict <- function(object, newdata_mat, h=1, cumulative=TRUE, fix=c("model", "order"), newxreg=NULL, ...){
  UseMethod("update_and_predict",object)
}

update_and_predict.agg_arima2 <- function(object, newdata_mat, h=1, cumulative=TRUE, fix=c("model", "order"), newxreg=NULL, ...){
  class(object) <- "agg_arima"
  
  res <- FTATS:::update_and_predict.agg_arima(object, newdata_mat, h, cumulative, fix, newxreg, ...)
  
  preds <- FTATS:::hierarchy_matrix_to_list(res, k=object$k)
  
  preds <- sapply(seq_along(preds), function(i) {
    if (!is.null(object$trans_back[[i]])) {
      data.frame(sapply(names(preds[[i]]), function(v) {
        if (v %in% names(object$trans_back[[i]])) {
          object$trans_back[[i]][[v]](preds[[i]][, v])
        } else {
          preds[[i]][, v]
        }
      }, simplify = FALSE))
    } else {
      preds[[i]]
    }
  }, simplify=FALSE)
  
  res <- FTATS:::hierarchy_list_to_matrix(preds, k=object$k)
  
  return(res)
}

# test <- fit_agg_arima2(x~1, 
#                        x=x,
#                        k = c(1,4), 
#                        trans=c(list(NULL), list(list(x=function(y)log(y+10)))),
#                        trans_back = c(list(NULL), list(list(x=function(y)exp(y)-10))))
# print(predict(test))
# print(fitted(test))