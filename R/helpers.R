construct_S_matrix <- function(k){
  do.call(rbind, sapply(k, function(ki) kronecker(diag(1,ki,ki), t(rep(1,k[length(k)]/ki)))))
}

hierarchy_list_to_matrix <- function(xlist, k, index_pos=1, value_pos=2){
  max_nrow <- max(sapply(seq_along(xlist), function(i) ceiling(nrow(xlist[[i]])/k[i])))
  mat <- do.call(cbind,
                 sapply(seq_along(xlist), function(i) {
                   x <- xlist[[i]][[value_pos]]
                   length(x) <- max_nrow*k[i]
                   matrix(x, ncol = k[i], byrow = TRUE, dimnames = list(NULL, paste(i,1:k[i],sep="-")))
                 }))
  if (nrow(xlist[[1]]) > 0){
    rownames(mat) <- seq(xlist[[1]][[index_pos]][1], length.out=nrow(mat), by=k[length(k)])
  } else {
    rownames(mat) <- seq(ceiling(xlist[[length(xlist)]][[index_pos]][1]/k[length(k)])*k[length(k)], length.out=nrow(mat), by=k[length(k)])
  }
  return(mat)
}

hierarchy_matrix_to_list <- function(xmat, k, drop_na=TRUE, value_name="x"){
  inds <- cbind(1+c(0,cumsum(k[-length(k)])),
                cumsum(k),
                k,
                deparse.level = 0)

  l <- apply(inds, 1, function(i){
    df <- data.frame(
      index = as.numeric(sapply(rownames(xmat), function(ind) rev(as.numeric(ind) - seq(from = 0, length.out = diff(i[1:2])+1, by = k[length(k)] / i[3])),USE.NAMES = FALSE)),
      x = as.numeric(t(xmat[, i[1]:i[2]]))
    )
    rownames(df) <- NULL
    colnames(df)[2] <- value_name
    if (drop_na) df <- df[rowSums(is.na(df))==0, ]
    return(df)
  }, simplify=FALSE)

  class(l) <- "agg_ts"
  attr(l, "k") <- k

  return(l)
}

extract_formula <- function(fmla, data){
  terms_fmla <- terms(fmla, data=data)
  attr_terms_fmla <- attributes(terms_fmla)
  vars_fmla <- as.character(attr_terms_fmla$variables)[-1]
  response <- as.numeric(data[[vars_fmla[attr_terms_fmla$response]]])

  if (length(labels(terms_fmla)) > 0 & any(colnames(data) %in% labels(terms_fmla))){
    regs <- (data[ ,labels(terms_fmla)])
  } else {
    regs <- NULL
  }

  return(list(response=response,
              regs=regs))
}

hierarchy_list_to_plot_data <- function(xlist, k=NULL){
  if (!is.null(attr(xlist, "k"))){
    k <- attr(xlist, "k")
  }
  stopifnot(!is.null(k))

  k2 <- rev(cumprod(rev(k[-1]/k[-length(k)])))
  bind_rows(sapply(seq_along(xlist), function(i) {
    xlist[[i]]$level <- i
    if (i < length(xlist)) {
      xlist[[i]]$index <- xlist[[i]]$index - k2[i] + 1
      xlist[[i]] <- rbind(xlist[[i]],
                          data.frame(index=max(xlist[[i]]$index) + k2[i] -1,
                                     x=tail(xlist[[i]]$x,n=1),
                                     level=i))
    }
    return(data.frame(xlist[[i]]))
  }, simplify = FALSE))
}

hierarchy_matrix_reduce_obs <- function(xmat, k, steps){
  obs_hid <- tail(cbind(1+c(0,cumsum(k[-length(k)])),
                        c(0,cumsum(k[-length(k)]))+c(0,steps)),
                  n=-1)
  k2 <- k[-1]/k[-length(k)]

  obs_hid <- unlist(apply(obs_hid, 1, function(x) {
    if (x[1] <= x[2]) x[1]:x[2]
  }))
  obs_hid <- colnames(xmat)[obs_hid]

  xmat_long <-
    reshape2::melt(xmat, varnames = c("index", "hid"), as.is = TRUE)

  xmat_long_obs <- xmat_long[xmat_long$hid %in% obs_hid, ]

  obs_df_list <- list(xmat_long_obs)
  i <- 1
  for (l in (length(k)-1):1){
    p_mat <- t(sapply(str_split(unlist(obs_df_list[[i]]$hid), "-"), function(x) {
      x <- as.numeric(x)
      if (x[1] == 1)
        return(c(0,0))
      c(x[1] - 1, paste(x[1] - 1, ((x[2] - 1) %/% k2[x[1] - 1]) + 1, sep = "-"))
    }))
    tmp <- obs_df_list[[i]] %>%
      mutate(parent_hid = p_mat[,2],
             parent_level = as.numeric(p_mat[,1])) %>%
      filter(parent_level == l) %>%
      group_by(index, parent_hid) %>%
      summarise(value=sum(value)) %>%
      rename(hid = parent_hid) %>%
      ungroup()

    obs_df_list <- c(obs_df_list, list(tmp))
    i <- i + 1
  }

  obs_df <- do.call(rbind, obs_df_list[-1])

  xmat_long <- left_join(xmat_long, obs_df, by=c("index","hid"), suffix=c(".x",".o"))
  xmat_long$value2 <- xmat_long$value.x - replace_na(xmat_long$value.o,0)

  xmat2 <- reshape2::acast(xmat_long, index~hid, value.var="value2")[rownames(xmat),colnames(xmat),drop=FALSE]

  return(xmat2[, !(colnames(xmat2) %in% obs_hid), drop=FALSE])
}

hierarchy_matrix_augment_obs <- function(xmat, k, obs_mat){
  k2 <- k[-1]/k[-length(k)]
  xmat_long <-
    reshape2::melt(xmat, varnames = c("index", "hid"), as.is = TRUE)
  obs_mat_long <-
    reshape2::melt(obs_mat, varnames = c("index", "hid"), as.is = TRUE)

  obs_df_list <- list(obs_mat_long)
  i <- 1
  for (l in (length(k)-1):1){
    p_mat <- t(sapply(str_split(unlist(obs_df_list[[i]]$hid), "-"), function(x) {
      x <- as.numeric(x)
      if (x[1] == 1)
        return(c(0,0))
      c(x[1] - 1, paste(x[1] - 1, ((x[2] - 1) %/% k2[x[1] - 1]) + 1, sep = "-"))
    }))
    tmp <- obs_df_list[[i]] %>%
      mutate(parent_hid = p_mat[,2],
             parent_level = as.numeric(p_mat[,1])) %>%
      filter(parent_level == l) %>%
      group_by(index, parent_hid) %>%
      summarise(value=sum(value)) %>%
      rename(hid = parent_hid) %>%
      ungroup()

    obs_df_list <- c(obs_df_list, list(tmp))
    i <- i + 1
  }

  obs_df <- do.call(rbind, obs_df_list[-1])

  xmat_long <- left_join(xmat_long, obs_df, by=c("index","hid"), suffix = c(".x",".o"))
  xmat_long$value2 <- xmat_long$value.x + replace_na(xmat_long$value.o, 0)

  xmat2 <- reshape2::acast(xmat_long, index~hid, value.var="value2")[rownames(xmat),colnames(xmat),drop=FALSE]

  resmat <- matrix(nrow=nrow(xmat),ncol=sum(k), dimnames = list(rownames(xmat),unlist(sapply(seq_along(k), function(i)paste(i,1:k[i],sep="-")))))
  resmat[, colnames(resmat) %in% colnames(xmat2)] <- xmat2
  resmat[, !(colnames(resmat) %in% colnames(xmat2))] <- obs_mat

  resmat[rowSums(is.na(resmat))>0, ] <- NA

  return(resmat)
}
