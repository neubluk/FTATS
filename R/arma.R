sim_arma <- function(n, phi=NULL, theta=NULL, Phi=NULL, Theta=NULL, m=1, sigma.sq=1){
  p <- length(phi)
  q <- length(theta)
  P <- length(Phi)
  Q <- length(Theta)

  if (P == 0 & Q == 0) m<-1

  model <- forecast::Arima(ts(rnorm(n,sd=sqrt(sigma.sq)),freq=m), order=c(p,0,q), seasonal=c(P,0,Q),
                           fixed=c(phi=phi, theta=theta, Phi=Phi, Theta=Theta, mu=0))
  x <- as.numeric(simulate(model, nsim=n))

  return(x)
}

sim_agg_arma <- function(n, k, phi=NULL, theta=NULL, Phi=NULL, Theta=NULL, m=1, bottom_mean = 0, effects=NULL, sigma.sq=1){
  nbottom <- n*k[length(k)]

  x <- sim_arma(nbottom, phi=phi, theta=theta, Phi=Phi, Theta=Theta, m=m, sigma.sq=sigma.sq)

  if (!is.null(effects)){
    #stopifnot(length(effects) == m)
    x <- x + rep(effects, nbottom/length(effects))
  }

  x <- x + bottom_mean

  aggs <- sapply(k, function(ki) {
    w <- k[length(k)] / ki
    y <- sapply(1:(n * ki), function(i) {
      sum(x[((i - 1) * w + 1):(i * w)])
    })
    data.frame(index=seq(w,length.out=length(y),by=w),x=y)
  }, simplify=FALSE)

  if (!is.null(effects)){
    aggs[[length(aggs)]] <- tibble(aggs[[length(aggs)]])
    aggs[[length(aggs)]]$xreg <- model.matrix(~factor(rep(effects, nbottom/length(effects)),
                                                      levels=unique(effects)))[,-1]
    colnames(aggs[[length(aggs)]]$xreg) <- unique(effects)[-1]
  }

  attr(aggs, "class") <- "agg_ts"
  attr(aggs, "k") <- k
  return(aggs)
}

#' Plot temporally aggregated ts
#'
#' @param object of class agg_ts
#' @param type of plot. ts to plot trajectories, acf/pacf to plot (partial-) autocorrelation functions
#' @param ... not used
#' @return ggplot2 object
#' @import ggplot2
#' @method plot agg_ts
#' @export
plot.agg_ts <- function(object, type=c("ts","acf","pacf"), ...){

  type <- match.arg(type)

  k <- attr(object, "k")
  k2 <- rev(cumprod(rev(k[-1]/k[-length(k)])))

  if (type == "ts") {
    data <- do.call(rbind,sapply(seq_along(object), function(i) {
      object[[i]]$level <- i
      if (i < length(object)) {
        object[[i]]$index <- object[[i]]$index - k2[i] + 1
        last_obs_complete <- tail(object[[i]],n=1)
        last_obs_complete$index <- max(object[[i]]$index) + k2[i] -1
        last_obs_complete$level <- i
        object[[i]] <- rbind(object[[i]],last_obs_complete)
      }
      return(data.frame(object[[i]]))
    }, simplify = FALSE))


    params <- list(...)
    if (!is.null(params$xvar)){
      xvar <- params$xvar
    } else {
      xvar <- "index"
    }
    if (is.null(data$fc)) {
      p <- ggplot(mapping = aes(
        x = .data[[xvar]],
        y = x,
        linetype = factor(level)
      ))
    } else {
      # if (length(unique(data$fc))>3){
      #   fc_split <- split(unique(data$fc), ceiling(seq_along(unique(data$fc))/3))
      #   group <- as.numeric(sapply(seq_along(fc_split), function(i)ifelse(data$fc %in% fc_split[[i]], i, NA)))
      #   data$group <- group[!is.na(group)]
      # }

      p <-
        ggplot(
          mapping = aes(
            x = .data[[xvar]],
            y = x,
            linetype = factor(level),
            color = fc,
          ),
          size = ifelse(data$fc == "obs", 1.5, 0.7)
        )

    }
    p <- p +
      geom_line(data = data[data$level == max(data$level),]) +
      geom_step(data = data[data$level != max(data$level),], direction="hv") +
      #geom_vline(data=data[data$level!=max(data$level), ], aes(xintercept=index), linetype="dotted", alpha=0.2) +
      scale_linetype_discrete("Level")

  } else if (type == "acf"){
    acfs <- sapply(seq_along(object), function(i) {
      a <- acf(object[[i]]$x, plot=FALSE)

      data.frame(lag = a$lag[,,1] + (i-1)/length(object)/3,
                 acf = a$acf[,,1],
                 ci = qnorm((1 + 0.95)/2)/sqrt(nrow(object[[i]])),
                 level = i)

    }, simplify=FALSE)

    acfs <- do.call(rbind, acfs)

    p <- ggplot(data = acfs,
                mapping = aes(x=lag, y=acf, color=factor(level)))+
      geom_point()+
      geom_hline(aes(yintercept = -ci, color=factor(level)), linetype="dotted")+
      geom_hline(aes(yintercept = ci, color=factor(level)), linetype="dotted")+
      geom_segment(aes(xend=lag, yend=0))+
      scale_x_continuous(breaks = 0:max(acfs$lag))

  } else if (type == "pacf"){
    acfs <- sapply(seq_along(object), function(i) {
      a <- pacf(object[[i]]$x, plot=FALSE)

      data.frame(lag = a$lag[,,1] + (i-1)/length(object)/3,
                 pacf = a$acf[,,1],
                 ci = qnorm((1 + 0.95)/2)/sqrt(nrow(object[[i]])),
                 level = i)

    }, simplify=FALSE)

    acfs <- do.call(rbind, acfs)

    p <- ggplot(data = acfs,
                mapping = aes(x=lag, y=pacf, color=factor(level)))+
      geom_point()+
      geom_hline(aes(yintercept = -ci, color=factor(level)), linetype="dotted")+
      geom_hline(aes(yintercept = ci, color=factor(level)), linetype="dotted")+
      geom_segment(aes(xend=lag, yend=0))+
      scale_x_continuous(breaks = 1:max(acfs$lag))

  }
  print(p)
}
