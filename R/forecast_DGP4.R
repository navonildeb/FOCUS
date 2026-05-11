source('library_causalTS.R')

r = 2 # rank is assumed to be known from the DGP

# range of T and N is 2^5 to 2^9
TT.arr = 2^c(5:9)
N.arr = 2^c(5:9)
R.max = 30
iter.cur = 0

rep_focus_h1 <- array(dim = c(length(N.arr), length(TT.arr), R.max))
rep_focus_h2 <- array(dim = c(length(N.arr), length(TT.arr), R.max))
rep_focus_h3 <- array(dim = c(length(N.arr), length(TT.arr), R.max))

rep_ms_h1 <- array(dim = c(length(N.arr), length(TT.arr), R.max))
rep_ms_h2 <- array(dim = c(length(N.arr), length(TT.arr), R.max))
rep_ms_h3 <- array(dim = c(length(N.arr), length(TT.arr), R.max))

timerep_focus <- array(dim = c(length(N.arr), length(TT.arr), R.max))
timerep_ms <- array(dim = c(length(N.arr), length(TT.arr), R.max))

for(n_ind in seq_along(N.arr)){
  N = N.arr[n_ind]
  for(t_ind in seq_along(TT.arr)){
    TT = TT.arr[t_ind]
    
    registerDoMC()
    options(cores = 15)
    total_mat = foreach(mc_iter = 1:R.max, .combine = rbind) %dopar% {
      file_name_C <- sprintf("data_files/DGP4/C0_files/DGP4_C0_N%d_T%d_iter%d.csv", N, TT, mc_iter)
      C0 <- read.csv(file_name_C)
      
      file_name_Y <- sprintf("data_files/DGP4/Y_files/DGP4_Y_N%d_T%d_iter%d.csv", N, TT, mc_iter)
      Y_full <- read.csv(file_name_Y)
      
      file_name_W <- sprintf("data_files/DGP4/W_files/DGP4_simultW_N%d_T%d_iter%d.csv", N, TT, mc_iter)
      W <- read.csv(file_name_W)
      
      Y <- Y_full[, c(1:TT)]
      
      t1 <- proc.time()
      # Estimate the factors and the loadings
      Sigma_eigen <- est_pca_missing(Y, W, r)
      F_est <- Sigma_eigen$factor_est
      L_est <- Sigma_eigen$loading_est
      
      # Fit VAR model + periodicity extraction
      res <- forecast_detpair_var1(F_est, h = 3, omega_range = c(1,50), M_peaks = 6, include_trend = TRUE)
      F_for <- res$forecast
      C_hat1 <- L_est[c(1:32), ] %*% t(F_for)
      focus_time <- (proc.time() - t1)[3]
      
      t1 <- proc.time()
      Y_ms <- Y
      Y_ms[W == 0] <- NA
      rownames(Y_ms) <- NULL
      colnames(Y_ms) <- NULL
      
      C_hat2 <- t(
        vapply(
          1:32,
          function(i) forecast.mssa(t(Y_ms), for_ind = i, t_start = TT+1, t_end = TT+3),
          numeric(3)
        )
      )
      ms_time <- (proc.time() - t1)[3]
      
      # Forecast error
      av_err1 = colMeans((C_hat1 - as.matrix(C0[1:32, ]))^2) #Focus
      av_err2 = colMeans((C_hat2 - as.matrix(C0[1:32, ]))^2)  #mSSA
      
      cat(paste("log(N) = ", log2(N), ", log(T) = ", log2(TT), ", mc_iter = ", mc_iter, "done \n"))
      c(av_err1, av_err2, focus_time, ms_time)
    }
    registerDoSEQ()
    
    total_err_mat_h1 <- total_mat[,c(1,4)]
    total_err_mat_h2 <- total_mat[,c(2,5)]
    total_err_mat_h3 <- total_mat[,c(3,6)]
    total_time_mat <- total_mat[,c(7,8)]
    
    rep_focus_h1[n_ind, t_ind, ] <- total_err_mat_h1[,1]
    rep_focus_h2[n_ind, t_ind, ] <- total_err_mat_h2[,1]
    rep_focus_h3[n_ind, t_ind, ] <- total_err_mat_h3[,1]
    
    rep_ms_h1[n_ind, t_ind, ] <- total_err_mat_h1[,2]
    rep_ms_h2[n_ind, t_ind, ] <- total_err_mat_h2[,2]
    rep_ms_h3[n_ind, t_ind, ] <- total_err_mat_h3[,2]
    
    timerep_focus[n_ind, t_ind, ] <- total_time_mat[,1]
    timerep_ms[n_ind, t_ind, ] <- total_time_mat[,2]
    
    iter.cur <- iter.cur + 1
    setTxtProgressBar(pb, 100 * iter.cur/(length(TT.arr) * length(N.arr)))
  }
}

time2 <- Sys.time()
cat(paste("Total simulation time:", time2 - time1))

out_dgp4 <- list(rep_focus_h1=rep_focus_h1,
                 rep_focus_h2=rep_focus_h2,
                 rep_focus_h3=rep_focus_h3,
                 rep_ms_h1=rep_ms_h1,
                 rep_ms_h2=rep_ms_h2,
                 rep_ms_h3=rep_ms_h3,
                 timerep_focus = timerep_focus,
                 timerep_ms = timerep_ms)

save(out_dgp4, file = "data_files/DGP4/DGP4_out.RData")

