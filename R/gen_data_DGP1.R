source("library_causalTS.R")

generate_data <- function(N.arr, TT.arr, R.max, r = 1, A0,
                          se.idio, se.loading, se.var, p_mcar, obs_pattern = "mcar"){
  for(n_ind in seq_along(N.arr)){
    # Fix the number of users N
    N = N.arr[n_ind]
    set.seed(n_ind)
    
    ### Generate the loadings
    L0 <- matrix(rnorm(N * r, mean = 0, sd = se.loading), nrow = N, ncol = r)
    
    ### Generate the unit-specific unobserved features
    S <- (rowSums(L0) >= 0)
    
    for(t_ind in seq_along(TT.arr)){
      ### Step 1: Generate the outcome matrix Y
      
      # Fix the number of time points T
      TT = TT.arr[t_ind]
      registerDoMC()
      options(cores = 30)
      
      total_arr = foreach(mc_iter = 1:R.max) %dopar% {
        set.seed(mc_iter)
        
        ## DGP1 (later convert to a if-else loop with DGP name)
        
        # Generate factors up to 3T/2 time points to include T/2 post treatment time points for SyN-BEATS
        F0_noise <- gen.var1(A0, n.time = 1.5*TT, se.var)
        F0_trend <- replicate(r,2*rep((1:(1.5*TT))/TT)^2)
        F0 <- F0_trend + F0_noise
        
        # Target of forecast 1,2,3-step ahead
        C0 <- cbind(L0 %*% (A0 %*%F0_noise[TT,] + F0_trend[TT+1, ]), 
                    L0 %*% (A0^2 %*%F0_noise[TT,] + F0_trend[TT+2, ]), 
                    L0 %*% (A0^3 %*%F0_noise[TT,] + F0_trend[TT+3, ]))
        
        ### Save the mean future outcomes to CSV
        file_name_C <- sprintf("data_files/DGP1/C0_files/DGP1_C0_N%d_T%d_iter%d.csv", N, TT, mc_iter)
        write.csv(C0, file = file_name_C, row.names = FALSE)
        
        Y_full <- L0 %*% t(F0) + array(rnorm(N*(1.5*TT), mean = 0, sd = se.idio), dim = c(N, 1.5*TT))
        
        ### Save the generated matrices to CSV
        ### Step 1: Save Y
        file_name_Y <- sprintf("data_files/DGP1/Y_files/DGP1_Y_N%d_T%d_iter%d.csv", N, TT, mc_iter)
        write.csv(Y_full, file = file_name_Y, row.names = FALSE)
        
        ### Step 2: Generate and save W
        if(obs_pattern == "mcar"){
          W <- matrix(rbinom(N * TT, size = 1, prob = p_mcar), nrow = N, ncol = TT)
          file_name_W <- sprintf("data_files/DGP1/W_files/DGP1_mcarW_N%d_T%d_iter%d.csv", N, TT, mc_iter)
          write.csv(W, file = file_name_W, row.names = FALSE)
          
        } else if (obs_pattern == "simultaneous"){
          W <- staggered_adoption_matrix_xp1(S, TT, 0.75, 0.75, 0.625, 0.375)$mat
          file_name_W <- sprintf("data_files/DGP1/W_files/DGP1_simultW_N%d_T%d_iter%d.csv", N, TT, mc_iter)
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
r = 1 # dimension of the factors/loadings
### Create the coefficient scalar if r = 1 (DGP1) /matrix if r > 1 (DGP2)
# A0 = matrix(c(0.5, -0.1, 0, 0.5), nrow = 2, ncol = 2, byrow = TRUE)
A0 = 0.5
se.idio = 0.1 # Idiosyncratic error sd
se.loading = 0.5 #Loading matrix sd
se.var = 0.5 # Factor-VAR process error sd
p_mcar = 0.7 # observation probability in MCAR
TT.arr = 2^c(5:9) # Array of T's 
N.arr = 2^c(5:9) # Array of N's
iter.cur = 0 # A counter variable for the Monte Carlo loops

obs_pattern = "simultaneous" # options: mcar, simultaneous

### Now generate the date with the set options
generate_data(N.arr, TT.arr, R.max, r = 1, A0, se.idio, se.loading, se.var, p_mcar, obs_pattern)
