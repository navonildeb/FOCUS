source('library_causalTS.R')
r = 1 # rank is assumed to be known from the DGP

# range of T and N is 2^5 to 2^9
TT.arr = 2^c(5:9)
N.arr = 2^c(5:6)
R.max = 30
iter.cur = 0

rep_focus_h1 <- array(dim = c(length(N.arr), length(TT.arr), R.max))
rep_focus_h2 <- array(dim = c(length(N.arr), length(TT.arr), R.max))
rep_focus_h3 <- array(dim = c(length(N.arr), length(TT.arr), R.max))

rep_ms_h1 <- array(dim = c(length(N.arr), length(TT.arr), R.max))
rep_ms_h2 <- array(dim = c(length(N.arr), length(TT.arr), R.max))
rep_ms_h3 <- array(dim = c(length(N.arr), length(TT.arr), R.max))

time1 <- Sys.time()

for(n_ind in seq_along(N.arr)){
  N = N.arr[n_ind]
  for(t_ind in seq_along(TT.arr)){
    TT = TT.arr[t_ind]
    
    registerDoMC()
    options(cores = 30)
    total_err_mat = foreach(mc_iter = 1:R.max, .combine = rbind) %dopar% {
      file_name_C <- sprintf("data_files/DGP0/C0_files/DGP0_C0_N%d_T%d_iter%d.csv", N, TT, mc_iter)
      C0 <- read.csv(file_name_C)
      
      file_name_Y <- sprintf("data_files/DGP0/Y_files/DGP0_Y_N%d_T%d_iter%d.csv", N, TT, mc_iter)
      Y_full <- read.csv(file_name_Y)
      Y <- Y_full[, c(1:TT)]
      
      file_name_W <- sprintf("data_files/DGP0/W_files/DGP0_mcarW_N%d_T%d_iter%d.csv", N, TT, mc_iter)
      W <- read.csv(file_name_W)
      
      Sigma_eigen <- est_pca_missing(Y, W, r)
      F_est <- Sigma_eigen$factor_est
      L_est <- Sigma_eigen$loading_est

      ### Task 1: Fit a VAR model
      ar_fit <- ar(F_est, aic = TRUE, order.max = 5, demean = FALSE)
      F_for_var1 <- predict(ar_fit, n.ahead = 3)$pred

      ### Forecast the common components

      C_hat1 <- L_est %*% t(F_for_var1)
      
      Y_ms <- Y
      Y_ms[W == 0] <- NA
      rownames(Y_ms) <- NULL
      colnames(Y_ms) <- NULL
      
      C_hat2 <- t(
        vapply(
          1:N,
          function(i) forecast.mssa(t(Y_ms), for_ind = i, t_start = TT+1, t_end = TT+3),
          numeric(3)   # expected length of forecast
        )
      )
      
      # Forecast error
      av_err1 = colMeans((C_hat1 - C0)^2) #Focus
      av_err2 = colMeans((C_hat2 - C0)^2)  #mSSA
      
      
      cat(paste("log(N) = ", log2(N), ", log(T) = ", log2(TT), ", mc_iter = ", mc_iter, "done \n"))
      
      c(av_err1, av_err2)
    }
    
    registerDoSEQ()
    
    total_err_mat_h1 <- total_err_mat[,c(1,4)]
    total_err_mat_h2 <- total_err_mat[,c(2,5)]
    total_err_mat_h3 <- total_err_mat[,c(3,6)]
    
    rep_focus_h1[n_ind, t_ind, ] <- total_err_mat_h1[,1]
    rep_focus_h2[n_ind, t_ind, ] <- total_err_mat_h2[,1]
    rep_focus_h3[n_ind, t_ind, ] <- total_err_mat_h3[,1]
    
    rep_ms_h1[n_ind, t_ind, ] <- total_err_mat_h1[,2]
    rep_ms_h2[n_ind, t_ind, ] <- total_err_mat_h2[,2]
    rep_ms_h3[n_ind, t_ind, ] <- total_err_mat_h3[,2]
    
    iter.cur <- iter.cur + 1
    setTxtProgressBar(pb, 100 * iter.cur/(length(TT.arr) * length(N.arr)))
  }
}

time2 <- Sys.time()
cat(paste("Total simulation time:", time2 - time1))

out_dgp0 <- list(rep_focus_h1=rep_focus_h1,
                 rep_focus_h2=rep_focus_h2,
                 rep_focus_h3=rep_focus_h3,
                 rep_ms_h1=rep_ms_h1,
                 rep_ms_h2=rep_ms_h2,
                 rep_ms_h3=rep_ms_h3)


save(out_dgp0, file = "data_files/DGP0/DGP0_out.RData")

