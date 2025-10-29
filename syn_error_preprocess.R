source("library_causalTS.R")

dgp_ind = 0 # change dgp_ind = 0, 1, 2 accordingly
TT_arr <- 2^c(5:8)
# N <- ifelse(dgp_ind == 0, 32, 64)
N <- 64
R.max = 30

load(sprintf("data_files/DGP%d/DGP%d_out.RData", dgp_ind, dgp_ind))

rep_syn_h1 <- array(dim = c(length(TT_arr), R.max))
rep_syn_h2 <- array(dim = c(length(TT_arr), R.max))
rep_syn_h3 <- array(dim = c(length(TT_arr), R.max))

for(t_ind in seq_along(TT_arr)){
  TT <- TT_arr[t_ind]
  
  syn_out <- foreach(mc_iter = 1:R.max, .combine = rbind) %dopar% {
    registerDoMC()
    options(cores = 30)
    file_name_syn <- sprintf("data_files/DGP%d/synout_files/N%d_T%d/synout_DGP%d_Y_N%d_T%d_iter%d.csv",
                             dgp_ind, N, TT, dgp_ind, N, TT, mc_iter)
    file_name_C <- sprintf("data_files/DGP%d/C0_files/DGP%d_C0_N%d_T%d_iter%d.csv",
                           dgp_ind, dgp_ind, N, TT, mc_iter)
    Chat_3 <- read.csv(file_name_syn) %>% dplyr::select(h1, h2, h3)
    C0 <- read.csv(file_name_C)
    
    mean_syn <- colMeans((Chat_3 - C0[1:32, ])^2)
    
    # return as numeric vector
    c(mean_syn)
  }
  registerDoSEQ()
  # syn_out is now an R.max × 3 matrix
  rep_syn_h1[t_ind,] <- syn_out[, "h1"]
  rep_syn_h2[t_ind,] <- syn_out[, "h2"]
  rep_syn_h3[t_ind,] <- syn_out[, "h3"]
}

syn_errors_final <- list(rep_syn_h1 = rep_syn_h1,
                         rep_syn_h2 = rep_syn_h2,
                         rep_syn_h3 = rep_syn_h3)
  
save(syn_errors_final, file = sprintf("data_files/DGP%d/synout_files/synout_DGP%d.RData", dgp_ind, dgp_ind))

