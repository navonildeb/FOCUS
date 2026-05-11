library(vars)
library(glmnet)
library(tseries)
library(astsa)
library(zoo)
library(knitr)
library(mgcv)
library(irlba)
library(foreach)
library(doParallel)
library(doMC)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(forecast)
library(reticulate)
Sys.setenv(RETICULATE_PYTHON='/usr/bin/python3')
# Sys.unsetenv("RETICULATE_PYTHON") # for unsetting the environment
# use_virtualenv("r-reticulate")
# np <- import("numpy") #test function to check numpy status
pd <- import("pandas")
os <- import("os")
ms <- import("mssa.mssa")
wg <- import("warnings")
wg$simplefilter("ignore")

pb <- txtProgressBar(min = 0, max = 100, 
                     style = 3, width = 50, char = "|")

colSD <- function(x, na.rm = FALSE) {
  apply(x, 2, sd, na.rm = na.rm)
}

heat_v0 <- function(X, main = "Heatmap") {
  p <- nrow(X)
  x_coords <- seq(0, 1, length.out = p + 1)
  y_coords <- seq(0, 1, length.out = p + 1)
  
  blue_palette <- colorRampPalette(c("white", "blue"))
  X <- X[, nrow(X):1]
  image.plot(x = x_coords, y = y_coords,
             z = X,
             col = blue_palette(100),
             zlim = c(0,1),
             xaxt = "n",
             yaxt = "n",
             xlab = "",
             ylab = "",
             legend.shrink = 0.6,
             main = main,
             cex.main = 2)
}


block_cross_validation <- function(time_series, block_size, spar_values) {
  n <- length(time_series)
  num_blocks <- floor(n / block_size)
  blocks <- split(time_series, rep(1:num_blocks, each = block_size, length.out = n))
  
  cv_errors <- sapply(spar_values, function(spar) {
    errors <- foreach(block_idx = 2:num_blocks, .combine = c) %do% {
      train_blocks <- blocks[c(1:(block_idx-1))]
      test_block <- blocks[[block_idx]]
      
      train_data <- unlist(train_blocks)
      fit <- smooth.spline(seq_along(train_data), train_data, spar = spar)
      
      test_indices <- seq_along(test_block) + (block_idx - 1) * block_size
      predictions <- predict(fit, test_indices)$y
      mean((test_block - predictions)^2)
    }
    mean(errors)
  })
  
  spar_values[which.min(cv_errors)]
}

# Main function to fit smoothing spline with optimal spar
fit_smoothing_spline <- function(time_series, block_size, spar_values = seq(0.1, 1, by = 0.1)) {
  best_spar <- block_cross_validation(time_series, block_size, spar_values)
  fit <- smooth.spline(seq_along(time_series), time_series, spar = best_spar)
  list(fit = fit, best_spar = best_spar)
}


# staggered adoption treatment pattern
staggered_adoption_matrix_xp1 <- function(S, TT, 
                                         prop_trt_row1,
                                         prop_available1,
                                         prop_trt_row0,
                                         prop_available0) {
  N = length(S)
  mat <- matrix(1, nrow = N, ncol = TT)
  propensity_scores <- matrix(1, nrow = N, ncol = TT)
  
  S1 <- which(S)
  S0 <- which(!S)
  
  # Assignments in rows with S = 1
  nrows_modify1 <- ceiling(prop_trt_row1 * sum(S==TRUE))
  rows_modify1 <- sample(S1, nrows_modify1)
  col_index1 <- ceiling(prop_available1 * TT)
  propensity_scores[S1, col_index1] <- 1
  propensity_scores[S1, (col_index1+1):TT] <- prop_available1
  for (i in rows_modify1) {
    mat[i, (col_index1 + 1):TT] <- 0
  }
  
  # Assignments in rows with S = 0
  nrows_modify0 <- ceiling(prop_trt_row0 * sum(S==FALSE))
  rows_modify0 <- sample(S0, nrows_modify0)
  col_index0 <- ceiling(prop_available0 * TT)
  propensity_scores[S0, col_index0] <- 1
  propensity_scores[S0, (col_index0+1):TT] <- prop_available0
  for (i in rows_modify0) {
    mat[i, (col_index0 + 1):TT] <- 0
  }
  
  return(list(mat = mat, propensity_scores = propensity_scores))
}


gen.var1 <- function(A, n.time, se.var, burn.in = 1000){
  if(is.null(dim(A)))
    r = 1
  else
    r = nrow(A)
  
  F.init <- array(rnorm(r, 0, se.var))
  for(k in 1:burn.in){
    F.init <- A %*% F.init + rnorm(r, mean = 0, sd = se.var)
  }
  
  F0 <- matrix(nrow = n.time, ncol = r)
  F0[1, ] <- A %*% F.init + 
    rnorm(r, mean = 0, sd = se.var)
  for(t in 2:n.time){
    F0[t, ] <- A %*% F0[(t-1), ] + 
      rnorm(r, mean = 0, sd = se.var)
  }
  return(F0)
}


gen.var1.general <- function(A, n.time, Sigma.var, burn.in = 1000){
  if(is.null(dim(A))){
    var_dim = 1
  } else {
    var_dim = nrow(A)
  }
  
  X.init <- matrix(mvtnorm::rmvnorm(1, mean = rep(0, var_dim), sigma = Sigma.var))
  for(k in 1:burn.in){
    X.init <- A %*% X.init + 
      matrix(mvtnorm::rmvnorm(1, mean = rep(0, var_dim), sigma = Sigma.var))
  }
  
  X0 <- matrix(nrow = n.time, ncol = var_dim)
  X0[1, ] <- A %*% X.init + 
    matrix(mvtnorm::rmvnorm(1, mean = rep(0, var_dim), sigma = Sigma.var))
  for(t in 2:n.time){
    X0[t, ] <- A %*% X0[(t-1), ] + 
      matrix(mvtnorm::rmvnorm(1, mean = rep(0, var_dim), sigma = Sigma.var))
  }
  return(X0)
}


gen.var1.sw2002 <- function(d, n.user, n.time, Sigma.var, burn.in = 1000){
  var_dim <- n.user
  
  X.init <- matrix(mvtnorm::rmvnorm(1, mean = rep(0, var_dim), sigma = Sigma.var))
  for(k in 1:burn.in){
    X.init <- d * X.init + 
      matrix(mvtnorm::rmvnorm(1, mean = rep(0, var_dim), sigma = Sigma.var))
  }
  
  X0 <- matrix(nrow = n.time, ncol = var_dim)
  X0[1, ] <- d * X.init + 
    matrix(mvtnorm::rmvnorm(1, mean = rep(0, var_dim), sigma = Sigma.var))
  for(t in 2:n.time){
    X0[t, ] <- d * X0[(t-1), ] + 
      matrix(mvtnorm::rmvnorm(1, mean = rep(0, var_dim), sigma = Sigma.var))
  }
  return(X0)
}

#### VARMA Simulator functions

sim_VARMA11 <- function(TT, A1, B1, Sigma, burn_in = 100) {
  r <- nrow(A1)
  TT_full <- TT + burn_in
  U <- matrix(0, TT_full, r)
  err <- MASS::mvrnorm(TT_full, mu = rep(0, r), Sigma = Sigma)
  U[1, ] <- err[1, ]
  for (t in 2:TT_full) {
    U[t, ] <- A1 %*% U[t - 1, ] + err[t, ] + B1 %*% err[t - 1, ]
  }
  U_out <- U[(burn_in + 1):TT_full, , drop = TRUE]
  return(U_out)
}

### A function for quantifying the population forecast of VARMA(1,1)
### given the past observations

VARMA11_forecast <- function(X, A1, B1, h) {
  TT <- nrow(X)
  r  <- ncol(X)
  eps <- matrix(0, nrow = TT, ncol = r)
  
  for (t in 1:TT) {
    X_lag   <- if (t - 1 >= 1) X[t-1, ] else rep(0, r)
    eps_lag <- if (t - 1 >= 1) eps[t-1, ] else rep(0, r)
    
    eps[t, ] <- X[t, ] -
      as.vector(A1 %*% X_lag) -
      as.vector(B1 %*% eps_lag)
  }
  # forecasts
  fcast <- matrix(0, nrow = h, ncol = r)
  
  # h = 1 MA part appears
  fcast[1, ] <- as.vector(A1 %*% X[TT, ]) +
    as.vector(B1 %*% eps[TT, ])
  
  # h >= 2 pure VAR recursion
  if (h >= 2) {
    for (k in 2:h) {
      fcast[k, ] <- as.vector(A1 %*% fcast[k-1, ])
    }
  }
  return(fcast)
}

v_VARMA11_cov <- function(A1, B1, se.varma){
  r <- nrow(A1)
  stopifnot(ncol(A1) == r, nrow(B1) == r, ncol(B1) == r)
  Sigma <- (se.varma^2) * diag(r)   # Var(eps_t)
  C <- (A1 + B1) %*% Sigma %*% t(A1 + B1)   # C = (A+B) Sigma (A+B)'
  # Solve S = C + A S A'  via vec(S) = (I - A \otimes A)^{-1} vec(C)
  K <- diag(r * r) - kronecker(A1, A1)
  vecS <- solve(K, as.vector(C))
  S <- matrix(vecS, nrow = r, ncol = r)
  Gamma0 <- Sigma + S  # Var(X_t) = Sigma + sum_{j>=0} A^j C (A^j)'
  v_var <- sum(diag(Gamma0)) / r
  return(v_var)
}

# v_target = v_VARMA11_cov(A1, B1, se.varma) in the data generation code
periodic.trend.DGP4_init <- function(TT, r, v_target, omega_range = c(1, 50), coef_range  = c(-1, 10), nsr = 1){
  P <- matrix(0, nrow = TT, ncol = r)
  for (k in 1:r) {
    a1 <- runif(1, min = coef_range[1], max = coef_range[2])
    a2 <- runif(1, min = coef_range[1], max = coef_range[2])
    w1 <- runif(1, min = omega_range[1], max = omega_range[2])
    w2 <- runif(1, min = omega_range[1], max = omega_range[2])
    pk <- a1 * cos(w1*(1:TT)/TT) + a2 * cos(w2*(1:TT)/TT)
    s2 <- var(pk)  # sample variance over t=1..TT
    scale <- if (s2 > 0) sqrt(v_target / (nsr * s2)) else 0 # change 3 to something else to control the SNR
    P[, k] <- scale * pk
  }
  return(P)
}

periodic.trend.DGP4 <- function(TT, r, v_target, nsr = 1){
  P <- matrix(0, nrow = TT, ncol = r)
  for (k in 1:r) {
    a1 <- 3
    a2 <- 6
    w1 <- 10
    w2 <- 20
    pk <- a1 * cos(w1*(1:TT)/TT) + a2 * cos(w2*(1:TT)/TT)
    s2 <- var(pk)  # sample variance over t=1..TT
    scale <- if (s2 > 0) sqrt(v_target / (nsr * s2)) else 0 # change 3 to something else to control the SNR
    P[, k] <- scale * pk
  }
  return(P)
}


#Estimate the factors and the loadings if the dimension of the latent space is given
est_pca <- function(Y, r){
  N_Y <- nrow(Y)
  TT_Y <- ncol(Y)
  Sigma_hat <- t(Y) %*% Y
  
  pc_Sigma <- prcomp_irlba(Sigma_hat/TT_Y, n = r, center = TRUE, scale. = FALSE)
  F_est <- sqrt(TT_Y) * pc_Sigma$rotation
  
  L_est <- Y %*% F_est / TT_Y
  return(list(factor_est = F_est, loading_est = L_est))
}


est_pca_missing_v0 <- function(Y, W, r){
  N <- nrow(Y)
  TT <- ncol(Y)
  Q <- matrix(list(), nrow = N, ncol = N)
  for(i in 1:N){
    for(j in 1:N){
      Q[[i, j]] <- which(W[i,] == 1 
                         & W[j,] == 1)
    }
  }
  Y.train <- ifelse(control_mat, Y, NA)

  Sigma_tilde <- matrix(nrow = N, ncol = N)
  for(i in 1:N){
    for(j in 1:N){
      run_ind <- Q[[i, j]]
      Sigma_tilde[i, j] <- sum(Y[i, run_ind] * Y[j, run_ind])/length(run_ind)
    }
  }
  
  Sigma.eigen <- eigen(Sigma_tilde/N)
  L.est <- sqrt(N) * as.matrix(Sigma.eigen$vectors[, (1:r)])
  F.est <- array(0, dim = c(TT, r))
  tmp1 <- control_mat
  tmp2 <- tmp1 * X
  
  for(j in 1:(TT-1)){
    tmp3 <- matrix(0, nrow = r, ncol = r)
    tmp4 <- numeric(r)
    for(i in 1:N){
      Li <- L.est[i,,drop = F]
      tmp3 <- tmp3 + tmp1[i,j] * (t(Li) %*% Li)
      tmp4 <- tmp4 + tmp2[i,j] * Li
    }
    F.est[j,] <- solve(tmp3) %*% tmp4
  }
  return(list(factor.est = F.est, loading.est = L_est))
}

est_pca_missing <- function(Y, W, r){
  N_Y <- nrow(Y)
  TT_Y <- ncol(Y)
  Q <- outer(1:TT_Y, 1:TT_Y, Vectorize(function(s, t) {
    which(W[, s] == 1 & W[, t] == 1)
  }))
  
  # Compute the covariance matrix estimator
  Sigma_hat <- outer(1:TT_Y, 1:TT_Y, Vectorize(function(s, t) {
    run_ind <- Q[[s, t]]
    sum(Y[run_ind, s] * Y[run_ind, t]) / length(run_ind)
  }))
  
  # PCA & estimate the factors
  # Sigma_eigen <- prcomp(Sigma_tilde/TT_Y)
  # F_est <- sqrt(TT_Y) * Sigma_eigen$rotation[,c(1:r)]

  Sigma_eigen <- prcomp_irlba(Sigma_hat/TT_Y, n = r, center = TRUE, scale. = FALSE)
  F_est <- sqrt(TT_Y) * Sigma_eigen$rotation
  if(is.null(dim(F_est))){
    F_est <- as.matrix(F_est)
  }
  
  L_est <- array(0, dim = c(N_Y, r))
  
  for(i in 1:N_Y){
    WFF <- matrix(0, nrow = r, ncol = r)
    WFY <- numeric(r)
    for(t in 1:TT_Y){
      
      Ft <- F_est[t,,drop = F]
      WFF <- WFF + W[i,t] * (t(Ft) %*% Ft)
      WFY <- WFY + (W[i,t] * Y[i,t]) * Ft
    }
    # if(any(eigen(WFF)$values < 10^-10)){
    #   ### regularize the covariance matrix with a small penalty
    #   WFF <- WFF + 10^-6 * diag(ncol(WFF))
    # }
    L_est[i,] <- WFY %*% solve(WFF)
  }
  return(list(cov_mat = Sigma_eigen, factor_est = F_est, loading_est = L_est))
}

# The preliminary var1 estimation function for the factors

est_var1 <- function(X, var.order = 1){
  if(!is.data.frame(X))
    X <- data.frame(X)
  colnames(X) <- as.vector(sapply('F', paste0, 1:r))
  x.var.list <- VAR(X, p = 1, type = "none", season = NULL)
  
  A_est <- matrix(nrow = r, ncol = r)
  for(i in 1:r){
    A_est[i, ] <- summary(x.var.list$varresult[[i]])$coefficients[,1]
  }
  return(list(x.var.list = x.var.list, A_est = A_est))
}

fit_VAR <- function(X, p = NULL, 
                    criterion = c("AIC", "BIC", "HQ", "FPE"), lag.max = 10,...) {
  if (!is.data.frame(X)) {
    X <- data.frame(X)
  }
  r <- ncol(X)
  colnames(X) <- paste0('F', 1:r)
  if (!requireNamespace("vars", quietly = TRUE)) {
    stop("Package 'vars' is required for this function.")
  }
  library(vars)
  
  # Determine the lag order if p is NULL
  if (is.null(p)) {
    criterion <- match.arg(criterion)
    selection <- VARselect(X, lag.max = lag.max, type = "none",...)
    # Check up to 10 lags
    p <- switch(criterion,
                AIC = selection$selection["AIC(n)"],
                BIC = selection$selection["SC(n)"],
                HQ = selection$selection["HQ(n)"],
                FPE = selection$selection["FPE(n)"])
    
    # message("Optimal lag order selected: p = ", p)
  }
  x.var.list <- VAR(X, p = p, type = "none",...)
  
  A_est <- vector("list", p)
  for (i in 1:p) {
    A_est[[i]] <- matrix(nrow = r, ncol = r)
  }
  
  for (i in 1:r) {
    coefficients <- summary(x.var.list$varresult[[i]])$coefficients[, 1]
    for (j in 1:p) {
      A_est[[j]][i, ] <- coefficients[((j - 1) * r + 1):(j * r)]
    }
  }
  
  return(list(p = p, x.var.list = x.var.list, A_est = A_est))
}

ar1_cov <- function(a, k){
  Sigma <- array(dim = c(k, k))
  for(i in 1:k)
    for(j in 1:k)
      Sigma[i, j] <- a^(abs(i-j))
  return(Sigma)
}

### Kalman filter for AR(1)
k_filter_ar1 <- function(x, a){
  n <- length(x)
  y <- numeric(n) 
  
  y[1] <- x[1]
  for(i in 2:n){
    y[i] <- a * x[i-1]
  }
  return(y)
}

### Kalman smoother for AR(1): variant 1
ks_smooth1_ar1 <- function(x, a) {
  n <- length(x)
  y <- numeric(n)  # initialize result vector
  
  # if (n < 2) {
  #   stop("Input vector 'x' must have at least 2 elements.")
  # }
  
  # Boundary cases
  y[1] <- a * x[2]
  y[n] <- a * x[n - 1]
  
  # Interior elements
  if (n > 2) {
    for (i in 2:(n - 1)) {
      y[i] <- (x[i - 1] + x[i + 1]) * a / (a^2 + 1)
    }
  }
  
  return(y)
}

### Kalman smoother for AR(1): variant 2
ks_smooth2_ar1 <- function(x, a){
  n <- length(x)
  y <- numeric(n)
  
  Sigma <- ar1_cov(a, n)
  
  for(i in 1:n){
    y[i] <- t(Sigma[i, -i]) %*% solve(Sigma[-i, -i]) %*% x[-i]
  }
  
  return(y)
}


# The preliminary var1 forecast function for the factors
forecast.var1 <- function(x.in, A, h){
  r <- nrow(A)
  x.for <- array(dim = c(h, r))
  x.iter <- x.in
  for(k in 1:h){
    x.iter <- A %*% x.iter
    x.for[k,] <- x.iter
  }
  return(x.for)
}

# A general function for forecasting a VAR(p) process up to a horizon h
forecast_VAR <- function(X, A_list, h = 1) { # h = forecast horizon
  r <- ncol(X)
  p <- length(A_list)
  n <- nrow(X)
  if(n < p){
    stop("Insufficient number of time points w.r.t the VAR order.")
  }

  history <- as.matrix(X[(n - p + 1):n, ])
  if(p == 1){
    history <- t(history)
  }
  forecast <- matrix(0, nrow = h, ncol = r)
  for (step in 1:h) {
    next_forecast <- rep(0, r)
    for (j in 1:p) {
      if (step - j <= 0) {  # Use history for initial steps
        next_forecast <- next_forecast + as.vector(A_list[[j]] %*% history[p + step - j, ])
      } else {  # Use previously forecasted values for further steps
        next_forecast <- next_forecast + as.vector(forecast[step - j, ] %*% A_list[[j]])
      }
    }
    forecast[step, ] <- next_forecast
  }
  return(forecast)
}


ridge_solve <- function(A, B, lambda = 1e-6) {
  solve(crossprod(A) + lambda * diag(ncol(A)), crossprod(A, B))
}

# map FFT frequencies -> your omega scale (w in cos(w*t/TT))
top_omega_peaks <- function(y, TT, omega_range = c(1, 50), M = 6) {
  y <- as.numeric(y)
  n <- length(y)
  Y <- fft(y - mean(y))
  spec <- Mod(Y)^2
  k <- 2:floor(n/2)                    # positive freqs
  angle <- 2*pi*(k-1)/n                # radians per step
  w_vals <- angle * TT                 # your omega scale
  
  keep <- which(w_vals >= omega_range[1] & w_vals <= omega_range[2])
  if (length(keep) == 0) return(integer(0))
  
  ord <- order(spec[k][keep], decreasing = TRUE)
  w_pick <- w_vals[keep][ord][1:min(M, length(ord))]
  unique(as.integer(round(w_pick)))
}

make_detX_omega <- function(tt, TT_full, omega = integer(0), include_trend = TRUE) {
  Xd <- matrix(1, length(tt), 1)
  colnames(Xd) <- "Intercept"
  if (include_trend) {
    Xd <- cbind(Xd, tt)
    colnames(Xd)[ncol(Xd)] <- "Trend"
  }
  for (w in omega) {
    Xd <- cbind(Xd, cos(w * tt / TT_full), sin(w * tt / TT_full))
    colnames(Xd)[(ncol(Xd)-1):ncol(Xd)] <- c(paste0("w", w, "_cos"), paste0("w", w, "_sin"))
  }
  Xd
}

# fit det with omega pair by BIC, then VAR(1) on residuals, forecast h
forecast_detpair_var1 <- function(X, h,
                                  omega_range = c(1, 50),
                                  M_peaks = 6,
                                  include_trend = TRUE,
                                  ridge_lambda = 1e-6) {
  if (!requireNamespace("vars", quietly = TRUE)) {
    stop("Install 'vars': install.packages('vars')")
  }
  
  X <- as.matrix(X)
  TT <- nrow(X); r <- ncol(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("V", 1:r)
  stopifnot(r == 2)
  
  t <- 1:TT
  t_future <- (TT+1):(TT+h)
  
  # --- Step 0: light prewhiten for detection ---
  var0 <- vars::VAR(X, p = 1, type = "const")
  U0 <- residuals(var0)          # (TT-1) x r
  # align indices (drop first time point)
  t_u <- 2:TT
  
  # --- Step 1: candidate omegas from periodogram peaks (union across cols) ---
  cand1 <- top_omega_peaks(U0[,1], TT = TT, omega_range = omega_range, M = M_peaks)
  cand2 <- top_omega_peaks(U0[,2], TT = TT, omega_range = omega_range, M = M_peaks)
  cand <- sort(unique(c(cand1, cand2)))
  
  # always allow "none" and "single"; and pairs among candidates
  models <- list(integer(0))
  models <- c(models, lapply(cand, function(w) c(w)))
  if (length(cand) >= 2) {
    pairs <- combn(cand, 2, simplify = FALSE)
    models <- c(models, pairs)
  }
  
  # --- Step 2: pick best deterministic model by BIC (fit on full data) ---
  best_bic <- Inf
  best_omega <- integer(0)
  best_Beta <- NULL
  best_Dhat <- NULL
  
  for (omega in models) {
    Xd <- make_detX_omega(t, TT_full = TT, omega = omega, include_trend = include_trend)
    Beta <- ridge_solve(Xd, X, lambda = ridge_lambda)
    Dhat <- Xd %*% Beta
    Uhat <- X - Dhat
    
    # Gaussian BIC on SSE (multivariate)
    sse <- sum(Uhat^2)
    kpar <- ncol(Xd) * r
    bic <- TT * log(sse / (TT*r)) + kpar * log(TT)
    
    if (bic < best_bic) {
      best_bic <- bic
      best_omega <- omega
      best_Beta <- Beta
      best_Dhat <- Dhat
    }
  }
  
  Uhat <- X - best_Dhat
  
  # VAR(1) on residuals + forecast
  var_fit <- vars::VAR(Uhat, p = 1, type = "const")
  pred <- predict(var_fit, n.ahead = h)
  
  U_fc <- cbind(pred$fcst[[1]][,"fcst"], pred$fcst[[2]][,"fcst"])
  colnames(U_fc) <- colnames(X)
  
  # deterministic forecast
  Xd_f <- make_detX_omega(t_future, TT_full = TT, omega = best_omega, include_trend = include_trend)
  D_fc <- Xd_f %*% best_Beta
  
  list(
    omega_used = best_omega,
    det_fitted = best_Dhat,
    resid_fitted = Uhat,
    forecast = D_fc + U_fc
  )
}



forecast.mssa <- function(data_in, for_ind, t_start, t_end = NA){
  if(!is.data.frame(data_in))
    data_in <- as.data.frame(data_in)
  colnames(data_in) <- paste0("X", 1:ncol(data_in))
  
  data_in <- r_to_py(data_in)
  ms_model = ms$mSSA()
  ms_model$update_model(data_in)
  
  if(is.na(t_end)) t_end = t_start
  if(t_start <= t_end){
    ms_predict <- ms_model$predict(paste0("X", for_ind), as.integer(t_start - 1), as.integer(t_end - 1)) 
  }
  return(ms_predict$`Mean Predictions`)
}

########################################################################################################
### Function for calculating coverage for DGP-1 in the paper

fit_ar1_scalar <- function(f) {
  f <- as.numeric(f)
  x <- f[-length(f)]
  y <- f[-1]
  phi_hat <- sum(x * y) / sum(x^2)
  eta_hat <- y - phi_hat * x
  
  list(
    phi = phi_hat,
    sigma_eta2 = mean(eta_hat^2),
    sigma_F2 = mean(f^2)   # sample variance proxy for stationary factor variance
  )
}

focus_ci_mcar_ar1 <- function(Y, W, h_max = 3, r = 1, alpha = 0.05,
                              obs_pattern = "mcar") {
  N <- nrow(Y)
  TT <- ncol(Y)
  deltaNT <- min(sqrt(N), sqrt(TT))
  
  Sigma_eigen <- est_pca_missing(Y, W, r)
  F_est <- as.numeric(Sigma_eigen$factor_est[, 1])
  L_est <- as.numeric(Sigma_eigen$loading_est[, 1])
  
  ar1 <- fit_ar1_scalar(F_est)
  phi_hat <- ar1$phi
  sigma_F2_hat <- mean(F_est^2)
  
  one_minus_phi2_hat <- max(0, 1 - phi_hat^2)
  
  W_mat <- as.matrix(W)
  Y_mat <- as.matrix(Y)
  
  p_hat <- mean(W_mat)
  
  Y_hat_in <- outer(L_est, F_est)
  resid_obs <- Y_mat - Y_hat_in
  sigma_eps2_hat <- mean(resid_obs[W_mat == 1]^2)
  
  sigma_L2_hat <- mean(L_est^2)
  
  kappa_L_hat <- mean((L_est^2 / sigma_L2_hat - 1)^2)
  
  F_T_hat <- F_est[TT]
  
  theta_hat <- sapply(
    1:h_max,
    function(h) L_est * (phi_hat^h) * F_T_hat
  )
  
  sigma_hat2 <- matrix(NA_real_, nrow = N, ncol = h_max)
  xi2_hat <- matrix(NA_real_, nrow = N, ncol = h_max)
  tau2_hat <- matrix(NA_real_, nrow = N, ncol = h_max)
  
  for (h in 1:h_max) {
    xi2_common <- deltaNT^2 * (phi_hat^(2 * h)) * sigma_eps2_hat *
      (L_est^2 / (N * sigma_L2_hat) +
         F_T_hat^2 / (TT * sigma_F2_hat))
    
    if (obs_pattern == "mcar") {
      xi2_base <- xi2_common / p_hat
      
      xi2_mcar_extra <- deltaNT^2 * (phi_hat^(2 * h)) *
        (L_est^2 * F_T_hat^2 / (N * sigma_L2_hat)) *
        (1 / p_hat - 1) *
        kappa_L_hat
      
      xi2_h <- xi2_base + xi2_mcar_extra
      
    } else if (obs_pattern == "staggered" || obs_pattern == "simultaneous") {
      xi2_h <- xi2_common
      
    } else {
      stop("obs_pattern must be one of: 'mcar', 'staggered', 'simultaneous'")
    }
    
    tau2_h <- deltaNT^2 * (h^2 / TT) *
      (phi_hat^(2 * h - 2)) *
      one_minus_phi2_hat *
      L_est^2 *
      F_T_hat^2
    
    xi2_hat[, h] <- xi2_h
    tau2_hat[, h] <- tau2_h
    sigma_hat2[, h] <- xi2_h + tau2_h
  }
  
  z <- qnorm(1 - alpha / 2)
  half_width <- z * sqrt(sigma_hat2) / deltaNT
  
  lower <- theta_hat - half_width
  upper <- theta_hat + half_width
  
  colnames(theta_hat) <- paste0("h", 1:h_max)
  colnames(lower) <- paste0("h", 1:h_max)
  colnames(upper) <- paste0("h", 1:h_max)
  colnames(sigma_hat2) <- paste0("h", 1:h_max)
  colnames(xi2_hat) <- paste0("h", 1:h_max)
  colnames(tau2_hat) <- paste0("h", 1:h_max)
  
  list(
    theta_hat = theta_hat,
    lower = lower,
    upper = upper,
    sigma_hat2 = sigma_hat2,
    xi2_hat = xi2_hat,
    tau2_hat = tau2_hat,
    deltaNT = deltaNT,
    
    estimates = list(
      F_est = F_est,
      L_est = L_est,
      phi_hat = phi_hat,
      sigma_F2_hat = sigma_F2_hat,
      sigma_eps2_hat = sigma_eps2_hat,
      sigma_L2_hat = sigma_L2_hat,
      kappa_L_hat = kappa_L_hat,
      p_hat = p_hat,
      F_T_hat = F_T_hat,
      one_minus_phi2_hat = one_minus_phi2_hat
    )
  )
}

coverage_pipeline_ar1 <- function(
    DGP = "DGP1",
    N.arr = 2^c(5:9),
    TT.arr = 2^c(5:9),
    R.max = 30,
    alpha = 0.05,
    h_max = 3,
    r = 1,
    obs_pattern = "mcar",
    ncores = 30) {
  
  cp_array <- array(NA_real_, dim = c(length(N.arr), length(TT.arr), h_max))
  avg_len  <- array(NA_real_, dim = c(length(N.arr), length(TT.arr), h_max))
  
  hist_out <- vector("list", length(N.arr) * length(TT.arr))
  names(hist_out) <- as.vector(outer(
    paste0("N", N.arr),
    paste0("T", TT.arr),
    paste,
    sep = "_"
  ))
  
  base_dir <- file.path("data_files", paste0("Covg_", DGP))
  
  counter <- 1
  
  for (n_ind in seq_along(N.arr)) {
    N <- N.arr[n_ind]
    
    for (t_ind in seq_along(TT.arr)) {
      TT <- TT.arr[t_ind]
      
      registerDoMC(cores = ncores)
      
      out <- foreach(mc_iter = 1:R.max) %dopar% {
        
        file_name_C <- sprintf(
          "%s/C0_files/%s_C0_N%d_T%d_iter%d.csv",
          base_dir, DGP, N, TT, mc_iter
        )
        C0 <- as.matrix(read.csv(file_name_C))
        C0 <- C0[, 1:h_max, drop = FALSE]
        
        file_name_Y <- sprintf(
          "%s/Y_files/%s_Y_N%d_T%d_iter%d.csv",
          base_dir, DGP, N, TT, mc_iter
        )
        Y_full <- as.matrix(read.csv(file_name_Y))
        Y <- Y_full[, 1:TT, drop = FALSE]
        
        if (obs_pattern == "mcar") {
          file_name_W <- sprintf(
            "%s/W_files/%s_mcarW_N%d_T%d_iter%d.csv",
            base_dir, DGP, N, TT, mc_iter
          )
        } else {
          file_name_W <- sprintf(
            "%s/W_files/%s_simultW_N%d_T%d_iter%d.csv",
            base_dir, DGP, N, TT, mc_iter
          )
        }
        
        W <- as.matrix(read.csv(file_name_W))
        
        fit <- focus_ci_mcar_ar1(
          Y = Y,
          W = W,
          h_max = h_max,
          r = r,
          alpha = alpha,
          obs_pattern = obs_pattern
        )
        
        covered <- C0 >= fit$lower & C0 <= fit$upper
        ci_len <- fit$upper - fit$lower
        
        cover_h <- colMeans(covered)
        len_h <- colMeans(ci_len)
        z_forecast <- fit$deltaNT * (fit$theta_hat - C0) / sqrt(fit$sigma_hat2)
        
        list(
          cover = cover_h,
          len = len_h,
          z_forecast = z_forecast
        )
      }
      
      registerDoSEQ()
      
      cover_mat <- do.call(rbind, lapply(out, `[[`, "cover"))
      len_mat <- do.call(rbind, lapply(out, `[[`, "len"))
      
      cp_array[n_ind, t_ind, ] <- colMeans(cover_mat)
      avg_len[n_ind, t_ind, ] <- colMeans(len_mat)
      
      z_array <- array(
        NA_real_,
        dim = c(N, h_max, R.max),
        dimnames = list(
          unit = paste0("i", 1:N),
          horizon = paste0("h", 1:h_max),
          iter = paste0("iter", 1:R.max)
        )
      )
      
      for (mc_iter in 1:R.max) {
        z_array[, , mc_iter] <- out[[mc_iter]]$z_forecast
      }
      
      hist_out[[counter]] <- list(
        N = N,
        TT = TT,
        z_forecast = z_array
      )
      
      cat(sprintf("Done: %s, N=%d, T=%d\n", DGP, N, TT))
      
      counter <- counter + 1
    }
  }
  
  list(
    DGP = DGP,
    obs_pattern = obs_pattern,
    cp_array = cp_array,
    avg_len = avg_len,
    hist_out = hist_out
  )
}
