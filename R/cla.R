#' @title Run the critical line algo
#' @param mu_vec vector or column matrix of expected returns
#' @param cov_mat expected covariance matrix
#' @param low_bound lower bound weights, see details
#' @param up_bound upper bound weights, see details
#' @param max_iter maximum number of iterations to run while searching for corner
#' @param low_vol_break optional volatility level to exit the solver if reached
#' portfolios
#' @return a list containing the weights for each corner portfolio \code{wgt_list}
#' @details The \code{low_bound} and \code{up_bound} weights are the corresponding
#' upper and lower bound weights for each asset in the \code{mu_vec}. They can
#' be entered as a single numeric value if all assets have the same upper or lower
#' bound weight or a vector of upper and lower bound weights can be entered.
#' The \code{low_vol_break} is an option to truncate the frontier at a certain
#' volatility level. For example if for a given problem you don't care about
#' portfolios with less than 2% volatility you can enter \code{0.02} to stop
#' the solver once it finds a corner portfolio with 2% vol.
#' @export
run_cla <- function(mu_vec, cov_mat, low_bound = 0, up_bound = 1,
                    max_iter = 1000, low_vol_break = 0, clean = TRUE) {

  mu_vec <- matrix(mu_vec, ncol = 1)
  store <- init_cla(mu_vec, low_bound, up_bound)
  n_assets <- nrow(mu_vec)
  store$l <- NA
  store$l_in <- NA
  store$l_out <- NA
  store$wgt_list <- list(store$wgt_vec)
  store$s <- list(NA)
  for (iter in 1:max_iter) {
    # Case a) bound a free weight
    l_in <- -100
    if (length(store$i_free) > 1) {
      s <- sub_mat(mu_vec, cov_mat, store$wgt_vec, store$i_free)
      j <- 1
      for (i in store$i_free) {
        l <- calc_lambda(s$cov_f_inv, s$cov_fb, s$mu_f, s$wgt_b, j,
                         c(store$low_bound[i], store$up_bound[i]))
        if (l$l_i > l_in) {
          l_in <- l$l_i
          i_in <- i
          b_i_in <- l$b_i
        }
        j <- j + 1
      }
    }
    # Case b) free an unbounded weight
    l_out <- -100
    if (length(store$i_free) < n_assets) {
      b <- get_b(store$i_free, n_assets)
      l_vec <- rep(0, length(b))
      for (i in b) {
        s <- sub_mat(mu_vec, cov_mat, store$wgt_vec, c(store$i_free, i))
        l <- calc_lambda(s$cov_f_inv, s$cov_fb, s$mu_f, s$wgt_b, length(s$mu_f),
                         store$wgt_vec[i])
        old_l <- store$l[length(store$l)]
        if ((is.na(old_l) | l$l_i < old_l) & l$l_i > l_out) {
          l_out <- l$l_i
          i_out <- i
        }
      }
    }
    # if lambda hits zero we are at the last corner portfolio and can exit
    if (l_out <= 0 & l_in <= 0) {
      break
    }
    # choose lambda
    if (max(c(l_in, -Inf), na.rm = TRUE) > l_out) {
      store$l <- c(store$l, l_in)
      store$i_free <- store$i_free[which(store$i_free != i_in)]
      store$wgt_vec[i_in] <- b_i_in
    } else {
      store$l <- c(store$l, l_out)
      store$i_free <- c(store$i_free, i_out)
    }
    store$l_in <- c(store$l_in, l_in)
    store$l_out <- c(store$l_out, l_out)
    # if all weights are bound we have hit the lowest corner portfolio and can
    # exit
    if (length(store$i_free) == 0) {
      break
    }
    # solve free weights and place into the stored weight list
    s <- sub_mat(mu_vec, cov_mat, store$wgt_vec, store$i_free)
    store$s[[iter + 1]] <- s
    wgt_f <- calc_wgt_f(s$cov_f_inv, s$cov_fb, s$mu_f, s$wgt_b,
                        store$l[length(store$l)])
    store$wgt_vec[store$i_free] <- wgt_f
    store$wgt_list[[length(store$wgt_list) + 1]] <- store$wgt_vec
    port_vol <- sqrt(t(store$wgt_vec) %*% cov_mat %*% store$wgt_vec)
    if (is.na(port_vol)) {
      port_vol <- Inf
    }
    if (port_vol <= low_vol_break) {
      break
    }
  }
  if (clean) {
    return(clean_store(store))
  } else {
    return(store)
  }
}


#' @title Solve weights for a target volatility
#' @param store output of \code{run_cla}
#' @param cov_mat covariance matrix
#' @param target_vol volatiltiy level to solve for
#' @param tol tolerance for volatility solution difference from \code{target_vol}
#' @param max_iter maximum iterations for attempted solution
calc_target_vol <- function(store, cov_mat, target_vol, tol = 0.0001,
                            max_iter = 10000) {

  port_vol <- sapply(store$wgt_list, calc_port_vol, cov_mat = cov_mat)
  vol_idx <- min(which(port_vol <= target_vol))
  if (vol_idx == Inf) {
    stop('target_vol is below possible vol range')
  }
  if (target_vol > port_vol[1]) {
    stop('target_vol is above possible vol range')
  }
  lambda <- store$l[c(vol_idx - 1, vol_idx)]
  wgt <- store$wgt_list[[vol_idx - 1]]
  for (i in 1:max_iter) {
    mean_lambda <- mean(lambda)
    s <- store$s[[vol_idx - 1]]
    wgt_f <- calc_wgt_f(s$cov_f_inv, s$cov_fb, s$mu_f, s$wgt_b, mean_lambda)
    wgt[s$f] <- wgt_f
    vol_i <- calc_port_vol(cov_mat, wgt)
    if (abs(lambda[1] - lambda[2]) <= tol) {
      break
    }
    if (vol_i < target_vol) {
      lambda[2] <- mean_lambda
    } else {
      lambda[1] <- mean_lambda
    }
  }
  return(wgt)
}


#' @title Initialize Critical Line Algo
#' @param mu_vec vector or column matrix of expected returns
#' @param low_bound lower bound weights
#' @param up_bound upper bound weights
#' @return list containing initial weight vector and free weight position
#' @export
init_cla <- function(mu_vec, low_bound, up_bound) {

  # param ----
  # mu_vec = vector of expected returns
  # low_bound = lower bound weights
  # up_bound = upper bound weights
  # ----
  # To initialize the algo, we first find the portfolio highest return given the
  # mu_vec and lower and upper bound constraints. The mu_vec and bounds are
  # structured as column matrixes and returned in the res list. The first free
  # weight index is returned in the res list as i_free.
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


#' @title Remove weight solutions that violate the lowerbound constraint
#' @param store storage list from \code{run_cla}
#' @note Poorly conditioned covariance matrices can lead to solutions with negative
#' weights.
#' @export
clean_store <- function(store) {
  wgt_err <- sapply(store$wgt_list, function(x) {sum(abs(x))}) > 1.01
  store$wgt_list <- store$wgt_list[!wgt_err]
  store$l <- store$l[!wgt_err]
  return(store)
}


#' @title Get bound weights from free weights
#' @param f free weight positions
#' @param n_assets total number of assets
#' @export
get_b <- function(f, n_assets) {

  b <- 1:n_assets
  b[!b %in% f]
}


#' @title Subset matrices into free and bound sections
#' @param mu_vec expected returns
#' @param cov_mat expected covariance matrix
#' @param wgt_vec asset weights
#' @param f free weight positions
#' @export
#' @importFrom MASS ginv
sub_mat <- function(mu_vec, cov_mat, wgt_vec, f) {

  res <- list()
  b <- get_b(f, nrow(mu_vec))
  res$cov_f <- cov_mat[f, f]
  res$cov_f_inv <- MASS::ginv(cov_mat[f, f])
  res$cov_fb <- cov_mat[f, b]
  res$mu_f <- mu_vec[f, 1]
  res$wgt_b <- wgt_vec[b, 1]
  res$f <- f
  return(res)
}


#' @title Calculate lambda
#' @param cov_f_inv inverse of covaraince of free weights
#' @param cov_fb covariance subset of free weight row and bound weight column
#' @param mu_f expected returns of free weights
#' @param wgt_b weights of bound assets
#' @param i position in lambda_i
#' @param b_i upper or lowerbound weight depending on intermediate calc
#' @export
calc_lambda <- function(cov_f_inv, cov_fb, mu_f, wgt_b, i, b_i) {

  one_f <- matrix(1, nrow = length(mu_f), ncol = 1)
  one_b <- matrix(1, nrow = length(wgt_b), ncol = 1)
  c1 <- t(one_f) %*% cov_f_inv %*% one_f
  c2 <- cov_f_inv %*% mu_f
  c3 <- t(one_f) %*% cov_f_inv %*% mu_f
  c4 <- cov_f_inv %*% one_f
  c_i <- -c1[1] * c2[i] + c3[1] * c4[i]
  if (length(b_i) > 1) {
    if (c_i > 0) {
      b_i <- b_i[2]
    } else {
      b_i <- b_i[1]
    }
  }
  lam1 <- t(one_b) %*% wgt_b
  lam2 <- cov_f_inv %*% cov_fb %*% wgt_b
  lambda_i <- ((1 - lam1 + t(one_f) %*% lam2) * c4[i] - c1 * (b_i + lam2[i])) / c_i
  res <- list()
  res$l_i <- lambda_i[1]
  res$b_i <- b_i
  return(res)
}


#' @title Calculate free weights
#' @param cov_f_inv inverse of covaraince of free weights
#' @param cov_fb covariance subset of free weight row and bound weight column
#' @param mu_f expected returns of free weights
#' @param wgt_b weights of bound assets
#' @param l lambda
#' @export
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


#' @title Portfolio return
#' @param mu expected return column matrix
#' @param w portfolio weight column
#' @export
calc_port_mu <- function(mu, w) {
  t(w) %*% mu
}


#' @title Portfolio volatility
#' @param cov_mat expected covariance matrix
#' @param w portfolio weight column
#' @export
calc_port_vol <- function(cov_mat, w) {
  sqrt(t(w) %*% cov_mat %*% w)[1]
}
