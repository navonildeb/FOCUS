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

# Main function to fit smoothing spline with optimal spar
fit_smoothing_spline <- function(time_series, block_size, spar_values = seq(0.1, 1, by = 0.1)) {
  best_spar <- block_cross_validation(time_series, block_size, spar_values)
  fit <- smooth.spline(seq_along(time_series), time_series, spar = best_spar)
  list(fit = fit, best_spar = best_spar)
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

### A function for quantifying the population forecast of ARMA(3,1)
### given the past observations

# F0: vector of observed series (F0_1,...,F0_{TT}), with dimension r = 1
# phi: AR(3) coefficients (e.g. c(phi1, phi2, phi3))
# theta: MA(1) coefficient
# h: forecast horizon (integer >= 1)
arma31_forecast <- function(F0, phi, theta, h) {
  TT <- length(F0)
  phi1 <- phi[1]; phi2 <- phi[2]; phi3 <- phi[3]
  theta1 <- theta
  
  # Compute innovations recursively
  eps <- numeric(TT)
  for (t in 1:TT) {
    F0_1 <- ifelse(t-1 >= 1, F0[t-1], 0)
    F0_2 <- ifelse(t-2 >= 1, F0[t-2], 0)
    F0_3 <- ifelse(t-3 >= 1, F0[t-3], 0)
    eps1 <- ifelse(t-1 >= 1, eps[t-1], 0)
    eps[t] <- F0[t] - (phi1 * F0_1 + phi2 * F0_2 + phi3 * F0_3 + theta1 * eps1)
  }
  
  # Storage for forecasts
  fcast <- numeric(h)
  
  # h = 1 forecast (includes MA part)
  fcast[1] <- phi1 * F0[TT] + phi2 * F0[TT-1] + phi3 * F0[TT-2] + theta1 * eps[TT]
  
  # h >= 2 forecasts: pure AR recursion
  if (h >= 2) {
    for (k in 2:h) {
      f1 <- fcast[k-1]
      f2 <- if (k > 2) fcast[k-2] else F0[TT]      # for h=2
      f3 <- if (k > 3) fcast[k-3] else F0[TT-1]    # for h=2,3
      fcast[k] <- phi1 * f1 + phi2 * f2 + phi3 * f3
    }
  }
  return(fcast)
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


