generate_ar_params <- function(p) {

  if (p==0) return(NULL)

  if (p==1) return(runif(1,-1,1))

  phi_helper <- matrix(0, p, p)
  pacorr <- sapply(1:p, function(k) {
    #rbeta(1, floor((k + 1) / 2), floor(k / 2) + 1)
    ExtDist::rBeta_ab(1, floor((k + 1) / 2), floor(k / 2) + 1, -1, 1)
  })

  diag(phi_helper) <- pacorr

  for (k in 2:p) {
    for (i in 1:(k - 1)) {
      phi_helper[k, i] <-
        phi_helper[k - 1, i] - pacorr[k] * phi_helper[k - 1, k - i]
    }
  }

  return(phi_helper[p, ])
}

compute_arma11_params <- function(phi, n, sigma.sq=1){

  if (is.na(phi)) return(
    c(beta=NA,
      eta=NA,
      sigma.sq_agg=NA)
  )

  phi_vec <- phi^(0:(n-1))

  A <- cbind(upper.tri(matrix(1,nrow=n,ncol=n),diag=TRUE), rbind(0,lower.tri(matrix(1,nrow=n-1,ncol=n-1),diag=TRUE)))
  cov_lagn <- cbind(
    rbind(matrix(0,n,n-1), diag(1,n-1)),
    rbind(matrix(0,2*n-1,n)))

  rho <- (t(phi_vec)%*%A%*%cov_lagn%*%t(A)%*%phi_vec) /
    (t(phi_vec)%*%A%*%t(A)%*%phi_vec)
  res_eta <- uniroot(function(x)x^2-x/rho+1,interval=c(-1,1))$root

  return(
    c(beta=phi^n,
      eta=res_eta,
      sigma.sq_agg = sigma.sq*(t(phi_vec)%*%A%*%t(A)%*%phi_vec)/(1+res_eta^2))
  )
}

compute_cov_agg <- function(phi, n, h, sigma.sq=1){
  return(
    sigma.sq/(1-phi)*((1-phi^h)/(1-phi)-phi^(n-h+1)*(1-phi^(2*h))/(1-phi^2))
  )
}

agg_ar_cov <- function(phi, n, sigma.sq=1){
  stopifnot(length(phi)==1)

  tmp <- matrix(0, nrow = n, ncol = n)
  phi_vec <- phi ^ (0:(n - 1))
  tmp[, 1] <- phi_vec
  for (i in 2:n) {
    tmp[-(1:(i - 1)), i] <- phi_vec[1:(n - i + 1)]
  }
  cov_bottom <- sigma.sq* tcrossprod(tmp)

  agg_params <- compute_arma11_params(phi,n,sigma.sq)

  covs_agg <- sapply(1:n, compute_cov_agg, phi=phi, n=n, sigma.sq=sigma.sq)
  first_col <- c(agg_params["sigma.sq_agg"], covs_agg, use.names=FALSE)
  res <- cbind(first_col, rbind(covs_agg, cov_bottom, deparse.level = 0) ,deparse.level = 0)

  return(res)
}

get_agg_order <- function(p,
                          k,
                          d = 0,
                          q = 0,
                          P = 0,
                          D = 0,
                          Q = 0,
                          m = c(0,0)
) {

  agg_order <- c("p" = p,
                 "d" = d,
                 "q" = floor(((p+1)*(k-1)+d*(k-1)+q)/k),
                 "P" = P,
                 "D" = D,
                 "Q" = floor( ( (P+D)*m[1]*k + (Q-P-D)*m[2])/k),
                 use.names=FALSE)

  return(agg_order)
}
