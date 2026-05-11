source("library_causalTS.R")
R.max <- 30
r <- 1
h_max <- 3
alpha <- 0.05
ncores <- 30

N.val <- 2^9
T.val <- 2^7
res_dgp1 <- coverage_pipeline_ar1(
  DGP = "DGP1",
  N.arr = N.val,
  TT.arr = T.val,
  R.max = R.max,
  alpha = alpha,
  h_max = h_max,
  r = r,
  obs_pattern = obs_pattern,
  ncores = ncores
)
print(sprintf("Coverage and length: (%.2f, %.2f)",round(100 * res_dgp1$cp_array[,,1], 2),res_dgp1$avg_len[,,1]))

key <- paste0("N", N.val, "_T", T.val)
x <- as.numeric(res_dgp1$hist_out[[key]]$z_forecast)
df <- data.frame(z_forecast = x)
p_hist <- ggplot(df, aes(x = z_forecast)) +
  geom_histogram(aes(y = after_stat(density)),bins = 30,fill = "#0072B2",
                 color = "white",alpha = 0.8) +
  stat_function(fun = dnorm,args = list(mean = 0, sd = 1),color = "lightcoral",linewidth = 1.1) +
  labs(title = paste0("DGP-1, N = ", N.val, ", T = ", T.val),
       x = expression(
         sqrt(delta[NT]) *
           (hat(theta)[paste(i, ",", T, ":", T+h)] -
              theta[paste(i, ",", T, ":", T+h)]) /
           hat(sigma)[paste(i, ",", T, ",", h)],
       )) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 18),
    panel.grid.minor = element_blank()
  )
p_hist



