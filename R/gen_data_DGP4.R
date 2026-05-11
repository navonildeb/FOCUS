source("library_causalTS.R")
generate_data <- function(N.arr, TT.arr, R.max, r = 2, A1, B1,
                          se.idio, se.loading, se.varma, p_mcar, obs_pattern = "mcar"){
  for(n_ind in seq_along(N.arr)){
    # Fix the number of users N
    N = N.arr[n_ind]
    set.seed(n_ind)
    
    ### Generate the loadings
    L0 <- matrix(rnorm(N*r, mean = 0, sd = se.loading), nrow = N, ncol = r)
    
    ### Generate the unit-specific unobserved features
    S <- (rowSums(L0) >= 0)
    
    for(t_ind in seq_along(TT.arr)){
      # Fix the number of time points T
      TT = TT.arr[t_ind]
      registerDoMC()
      options(cores = 30)
      
      total_arr = foreach(mc_iter = 1:R.max) %dopar% {
        set.seed(mc_iter)
        # Generate factors up to 3T/2 time points to include T/2 post treatment time points for SyN-BEATS
        F0_noise <- sim_VARMA11(1.5*TT, A1, B1, Sigma = se.varma^2 *diag(r), burn_in = 1000)
        
        v_target = v_VARMA11_cov(A1, B1, se.varma)
        F0_per <- periodic.trend.DGP4(1.5*TT, r, v_target, nsr = 1)
        
        F0 <- F0_per + F0_noise
        
        # Target of forecast 1,2,3-step ahead
        F0_for <- VARMA11_forecast(F0_noise[1:TT, ], A1, B1, 3) + F0_per[TT + c(1:3), ]
        C0 <- L0 %*% t(F0_for)
        
        ### Save the mean future outcomes to CSV
        file_name_C <- sprintf("data_files/DGP4/C0_files/DGP4_C0_N%d_T%d_iter%d.csv", N, TT, mc_iter)
        write.csv(C0, file = file_name_C, row.names = FALSE)
        
        Y_full <- L0 %*% t(F0) + array(rnorm(N*(1.5*TT), mean = 0, sd = se.idio), dim = c(N, 1.5*TT))
        
        ### Save the generated matrices to CSV
        file_name_Y <- sprintf("data_files/DGP4/Y_files/DGP4_Y_N%d_T%d_iter%d.csv", N, TT, mc_iter)
        write.csv(Y_full, file = file_name_Y, row.names = FALSE)
        
        if(obs_pattern == "mcar"){
          W <- matrix(rbinom(N * TT, size = 1, prob = p_mcar), nrow = N, ncol = TT)
          file_name_W <- sprintf("data_files/DGP4/W_files/DGP4_mcarW_N%d_T%d_iter%d.csv", N, TT, mc_iter)
          write.csv(W, file = file_name_W, row.names = FALSE)
          
        } else if (obs_pattern == "simultaneous"){
          W <- staggered_adoption_matrix_xp1(S, TT, 0.75, 0.75, 0.625, 0.375)$mat
          file_name_W <- sprintf("data_files/DGP4/W_files/DGP4_simultW_N%d_T%d_iter%d.csv", N, TT, mc_iter)
          write.csv(W, file = file_name_W, row.names = FALSE)
        }
      }
      registerDoSEQ()
      
      iter.cur <- iter.cur + 1
      setTxtProgressBar(pb, 100 * iter.cur/(length(TT.arr) * length(N.arr)))
    }
  }
}


### Set the options
R.max <- 30 # Max-number of Monte Carlo replicates
r = 2 # dimension of the factors/loadings
### Create the coefficient scalars

A1 <- matrix(c(0.5, 0.3, -0.2, 0.5), nrow = 2, byrow = 2)
B1 <- matrix(c(0.3, 0, 0, 0.3), nrow = 2, byrow = 2)

se.idio = 0.1 # Idiosyncratic error sd
se.loading = 0.5 #Loading matrix sd
se.varma = 0.5 # Factor-VARMA(1,1) process error sd
p_mcar = 0.7 # observation probability in MCAR
TT.arr = 2^c(5:9) # Array of T's 
N.arr = 2^c(5:9) # Array of N's
iter.cur = 0 # A counter variable for the Monte Carlo loops

obs_pattern = "simultaneous" # options: mcar, simultaneous

### Now generate the date with the set options
generate_data(N.arr, TT.arr, R.max, r, A1, B1, se.idio, se.loading, se.varma, p_mcar, obs_pattern)
