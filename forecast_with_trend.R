# Estimate factors by PCA, then do forecast on var1, TREND INCLUDED
source('library_causalTS.R')
library(foreach)
library(doMC)
registerDoMC()
options(cores = 50)

r = 2
A <- matrix(0, nrow = r, ncol = r)
diag(A) <- 0.8
for(i in 2:r){
  A[(i-1), i] = -0.1
}

# The noise standard errors
se.idio = 0.1 # Idiosyncratic error
se.loading = 0.5 #Loading matrix
se.var = 0.5 # Factor-VAR process


TT.arr = 2^c(5:10)
N.arr = 2^c(5:10)
R.max = 50

iter.cur = 0
total.arr = array(dim = c(length(N.arr), length(TT.arr), 3))

time1 <- proc.time()
for(i in seq_along(N.arr)){
  N = N.arr[i]
  L0 <- matrix(rnorm(N*r, mean = 1, sd = se.loading), nrow = N, ncol = r)
  for(j in seq_along(TT.arr)){
    TT = TT.arr[j]
    av.err.arr = foreach(rep = 1:R.max, .combine = rbind) %dopar% {
      source('library_causalTS.R')
      F0_noise <- gen.var1(A, TT, se.var)
      # v = ((1:TT)/TT)^2
      v = rep(0, TT)
      F0_trend <- replicate(r, 10 * v)
      F0 = F0_trend + F0_noise
      # Target coordinate of C(i,t) for imputation/forecast
      # C0 <- sum(L0[1,] * F0[1,]) # Target of estimation
      C0 <- L0[1,] %*% (A %*% F0_noise[TT-1,] + 
                          F0_trend[TT, ]) # Target of forecast
      
      X <- L0 %*% t(F0) + array(rnorm(N*TT, mean = 0, sd = se.idio), dim = c(N, TT))
      
      X.train <- X[,-TT]
      
      # X.test <- X[,TT]
      
      # Step 1. PCA
      XX.eigen <- eigen(t(X.train) %*% X.train)
      F.est <-  sqrt(nrow(XX.eigen$vectors)) * XX.eigen$vectors[, (1:r)]
      L.est <- X.train %*% F.est / nrow(XX.eigen$vectors)
      
      # Step 2: Detrend with fitting spline
      F.est.noise <- array(dim = dim(F.est))
      F.for.trend <- numeric(r)
      
      for(j in 1:r){
        tmp <- F.est[,j]
        result <- fit_smoothing_spline(tmp, block_size = ceiling(length(tmp)/10))
        F.est.noise[,j] <- tmp - result$fit$y
        F.for.trend[j] <- predict(result$fit,TT)$y
      }
      
      A.est <- est.var1(F.est.noise)$A.est
      
      F.for <- A.est %*% F.est.noise[nrow(F.est.noise), ] + F.for.trend
      
      
      av.err1 = abs(sum(L.est[1,] * F.for) - C0) # Common component forecast
      # av.err1 = (sum(L.est[1,] * F.est[1,]) - C0)^2 # Common component estimation
      av.err2 = abs(forecast.mssa(X.train[1,], for_ind = 1, t_start = TT) - C0)
      
      av.err3 = abs(forecast.mssa(t(X.train), for_ind = 1, t_start = TT) - C0)
      
      c(av.err1, av.err2, av.err3)
    }
    av.err.arr = colMeans(av.err.arr)
    total.arr[i,j,] = av.err.arr
    
    iter.cur <- iter.cur + 1
    setTxtProgressBar(pb, 100 * iter.cur/(length(TT.arr) * length(N.arr)))
    # print(c(i,j))
  }
}
time2 <- proc.time()
print(time2 - time1)

save(total.arr, file = "N5_10_T5_10_R50_notrend.RData")
