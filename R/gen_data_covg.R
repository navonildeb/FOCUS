source("library_causalTS.R")

arma_forecast_known <- function(x, ar = numeric(0), ma = numeric(0), h = 3) {
  x <- as.numeric(x)
  TT <- length(x)
  p <- length(ar)
  q <- length(ma)
  
  # Recover innovations recursively under the known ARMA parameters.
  e <- numeric(TT)
  for (t in seq_len(TT)) {
    ar_part <- if (p > 0) {
      sum(sapply(seq_len(p), function(j) {
        if (t - j >= 1) ar[j] * x[t - j] else 0
      }))
    } else 0
    
    ma_part <- if (q > 0) {
      sum(sapply(seq_len(q), function(j) {
        if (t - j >= 1) ma[j] * e[t - j] else 0
      }))
    } else 0
    
    e[t] <- x[t] - ar_part - ma_part
  }
  
  x_ext <- c(x, rep(NA_real_, h))
  
  for (ell in seq_len(h)) {
    t <- TT + ell
    
    ar_part <- if (p > 0) {
      sum(sapply(seq_len(p), function(j) {
        ar[j] * x_ext[t - j]
      }))
    } else 0
    
    ma_part <- if (q > 0) {
      sum(sapply(seq_len(q), function(j) {
        # Future innovations are zero in conditional mean forecasts.
        if (t - j <= TT) ma[j] * e[t - j] else 0
      }))
    } else 0
    
    x_ext[t] <- ar_part + ma_part
  }
  
  matrix(x_ext[(TT + 1):(TT + h)], ncol = h)
}

get_dgp_spec <- function(DGP, A1, A2, A3, B1) {
  if (DGP == "DGP1") {
    return(list(ar = c(A1), ma = numeric(0)))
  }
  
  if (DGP == "DGP2") {
    return(list(ar = c(A1, A2, A3), ma = numeric(0)))
  }
  
  if (DGP == "DGP3") {
    return(list(ar = c(A1, A2, A3), ma = c(B1)))
  }
  
  stop("DGP must be one of: DGP1, DGP2, DGP3")
}

generate_data <- function(N.arr, TT.arr, R.max, DGP,
                          r = 1, A1, A2, A3, B1,
                          se.idio, se.loading, se.arma,
                          p_mcar, obs_pattern = "mcar") {
  
  dgp_spec <- get_dgp_spec(DGP, A1, A2, A3, B1)
  
  base_dir <- file.path("data_files", paste0("Covg_", DGP))
  C_dir <- file.path(base_dir, "C0_files")
  Y_dir <- file.path(base_dir, "Y_files")
  W_dir <- file.path(base_dir, "W_files")
  
  dir.create(C_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(Y_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(W_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (n_ind in seq_along(N.arr)) {
    N <- N.arr[n_ind]
    set.seed(n_ind)
    
    L0 <- matrix(rnorm(N * r, mean = 0, sd = se.loading), nrow = N, ncol = r)
    
    S <- rowSums(L0) >= 0
    
    for (t_ind in seq_along(TT.arr)) {
      TT <- TT.arr[t_ind]
      
      registerDoMC()
      options(cores = 30)
      
      total_arr <- foreach(mc_iter = 1:R.max) %dopar% {
        set.seed(mc_iter)
        
        arma_model <- list(ar = dgp_spec$ar)
        if (length(dgp_spec$ma) > 0) {
          arma_model$ma <- dgp_spec$ma
        }
        
        F0 <- arima.sim(
          model = arma_model,
          n = 1.5 * TT,
          sd = se.arma
        )
        
        F_forecast <- arma_forecast_known(
          F0[1:TT],
          ar = dgp_spec$ar,
          ma = dgp_spec$ma,
          h = 3
        )
        
        C0 <- L0 %*% F_forecast
        
        file_name_C <- sprintf(
          "%s/%s_C0_N%d_T%d_iter%d.csv",
          C_dir, DGP, N, TT, mc_iter
        )
        write.csv(C0, file = file_name_C, row.names = FALSE)
        
        Y_full <- L0 %*% t(F0) +
          array(rnorm(N * (1.5 * TT), mean = 0, sd = se.idio),
                dim = c(N, 1.5 * TT))
        
        file_name_Y <- sprintf(
          "%s/%s_Y_N%d_T%d_iter%d.csv",
          Y_dir, DGP, N, TT, mc_iter
        )
        write.csv(Y_full, file = file_name_Y, row.names = FALSE)
        
        if (obs_pattern == "mcar") {
          W <- matrix(rbinom(N * TT, size = 1, prob = p_mcar),
                      nrow = N, ncol = TT)
          
          file_name_W <- sprintf(
            "%s/%s_mcarW_N%d_T%d_iter%d.csv",
            W_dir, DGP, N, TT, mc_iter
          )
          write.csv(W, file = file_name_W, row.names = FALSE)
          
        } else if (obs_pattern == "simultaneous") {
          W <- staggered_adoption_matrix_xp1(S, TT, 0.75, 0.75, 0.625, 0.375)$mat
          
          file_name_W <- sprintf(
            "%s/%s_simultW_N%d_T%d_iter%d.csv",
            W_dir, DGP, N, TT, mc_iter
          )
          write.csv(W, file = file_name_W, row.names = FALSE)
        }
      }
      
      registerDoSEQ()
      
      iter.cur <<- iter.cur + 1
      setTxtProgressBar(pb, 100 * iter.cur / (length(TT.arr) * length(N.arr)))
    }
  }
}




R.max <- 30
r <- 1

DGP <- "DGP1"   # options: "DGP1", "DGP2", "DGP3"

A1 <- 0.5
A2 <- -0.4
A3 <- 0.2
B1 <- 0.5

se.idio <- 0.1
se.loading <- 0.5
se.arma <- 0.7
p_mcar <- 0.7

TT.arr <- 2^c(5:9)
N.arr <- 2^c(5:9)
iter.cur <- 0

obs_pattern <- "mcar"

generate_data(
  N.arr, TT.arr, R.max, DGP,
  r = r,
  A1 = A1, A2 = A2, A3 = A3, B1 = B1,
  se.idio = se.idio,
  se.loading = se.loading,
  se.arma = se.arma,
  p_mcar = p_mcar,
  obs_pattern = obs_pattern
)
