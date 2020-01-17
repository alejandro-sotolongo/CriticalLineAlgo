
init_cla <- function(mu_vec, up_bound = 1, low_bound = 0) {

  n_assets <- nrow(mu_vec)
  up_bound_vec <- matrix(up_bound, ncol = 1, nrow = n_assets)
  low_bound_vec <- matrix(low_bound, ncol = 1, nrow = n_assets)
  if (sum(low_bound_vec) > 1) {
    stop('sum of lower bound weights is greater than 1')
  }
  mu_vec_rank <- n_assets - rank(mu_vec) + 1
  wgt_vec <- low_bound_vec
  remain_wgt <- 1 - sum(wgt_vec)
  for (i in 1:length(mu_vec)) {
    if (remain_wgt == 0) {
      break
    }
    max_mu_idx <- which(mu_vec_rank == i)
    up_bound_i <- up_bound_vec[max_mu_idx, 1]
    add_wgt <- up_bound_i - wgt_vec[max_mu_idx]
    if (add_wgt <= remain_wgt) {
      remain_wgt <- remain_wgt - (up_bound_i - wgt_vec[max_mu_idx])
      wgt_vec[max_mu_idx] <- up_bound_i
    } else {
      wgt_vec[max_mu_idx] <- wgt_vec[max_mu_idx] + remain_wgt
      remain_wgt <- 0
    }
  }
  res <- list()
  res$wgt_vec <- wgt_vec
  res$i_free <- max_mu_idx
  res$wgt_vec <- wgt_vec
  res$up_bound <- up_bound_vec
  res$low_bound <- low_bound_vec
  return(res)
}

get_b <- function(f, n_assets) {
  b <- 1:n_assets
  b[!b %in% f]
}

sub_mat <- function(mu_vec, cov_mat, wgt_vec, f) {

  res <- list()
  b <- get_b(f, nrow(mu_vec))
  res$cov_f <- cov_mat[f, f]
  res$cov_f_inv <- MASS::ginv(cov_mat[f, f])
  res$cov_fb <- cov_mat[f, b]
  res$mu_f <- mu_vec[f, 1]
  res$wgt_b <- wgt_vec[b, 1]
  return(res)
}

calc_lambda <- function(cov_f_inv, cov_fb, mu_f, wgt_b, i, b_i) {

  one_f <- matrix(1, nrow = length(mu_f), ncol = 1)
  one_b <- matrix(1, nrow = length(wgt_b), ncol = 1)
  c1 <- t(one_f) %*% cov_f_inv %*% one_f
  c2 <- cov_f_inv %*% mu_f
  c3 <- t(one_f) %*% cov_f_inv %*% mu_f
  c4 <- cov_f_inv %*% one_f
  c_i <- -c1[1] * c2[i] + c3[1] * c4[i]
  lam1 <- t(one_b) %*% wgt_b
  lam2 <- cov_f_inv %*% cov_fb %*% wgt_b
  lambda_i <- ((1 - lam1 + t(one_f) %*% lam2) * c4[i] - c1 * (b_i + lam2[i])) / c_i
  res <- list()
  res$l_i <- lambda_i[1]
  res$b_i <- b_i
  return(res)
}

calc_wgt_f <- function(cov_f_inv, cov_fb, mu_f, wgt_b, l) {

  one_f <- matrix(1, nrow = length(mu_f), ncol = 1)
  one_b <- matrix(1, nrow = length(wgt_b), ncol = 1)
  g1 <- t(one_f) %*% cov_f_inv %*% mu_f
  g2 <- t(one_f) %*% cov_f_inv %*% one_f
  g3 <- t(one_b) %*% wgt_b
  g4 <- t(one_f) %*% cov_f_inv %*% cov_fb %*% wgt_b
  gam <- -l * g1 / g2 + (1 - g3 + g4) / g2
  w1 <- cov_f_inv %*% cov_fb %*% wgt_b
  w2 <- cov_f_inv %*% one_f
  w3 <- cov_f_inv %*% mu_f
  wgt_f <- -w1 + gam[1] * w2 + l * w3
  return(wgt_f)
}

run_cla <- function() {

  # delete later
  data(ETF)
  mu_vec <- colMeans(ret[, 2:ncol(ret)])
  mu_vec <- matrix(mu_vec, ncol = 1)
  cov_mat <- cov(ret[, 2:ncol(ret)])
  # ----------
  store <- init_cla(mu_vec)
  n_assets <- nrow(mu_vec)
  store$l <- NA
  store$wgt_list <- list(store$wgt_vec)

  # Case a)
  l_in <- NA


  # Case b)
  l_out <- NA
  if (length(store$i_free) < n_assets) {
    b <- get_b(store$i_free, n_assets)
    l_vec <- rep(0, length(b))
    j <- 1
    for (i in b) {
      s <- sub_mat(mu_vec, cov_mat, store$wgt_vec, c(store$i_free, i))
      l <- calc_lambda(s$cov_f_inv, s$cov_fb, s$mu_f, s$wgt_b, length(s$mu_f),
                       store$wgt_vec[i])
      j <- j + 1
      if ((is.na(store$l) | l$l_i < store$l) & (l$l_i > 0 | l$l_i > l_out)) {
        l_out <- l$l_i
        i_out <- i
      }
    }
  }
  # descide lambda
  if (max(l_in, 0, na.rm) > l_out) {

  } else {
    store$l <- c(store$l, l_out)
    store$i_free <- c(store$i_free, i_out)
  }
  s <- sub_mat(mu_vec, cov_mat, wgt_vec, store$i_free)
  wgt_f <- calc_wgt_f(s$cov_f_inv, s$cov_fb, s$mu_f, s$wgt_b,
                      store$l[length(store$l)])
  store$wgt_vec[store$i_free] <- wgt_f
  store$wgt_list[[length(store$wgt_list) + 1]] <- store$wgt_vec
}
