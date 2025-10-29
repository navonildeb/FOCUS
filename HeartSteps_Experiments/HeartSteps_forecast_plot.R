source('library_causalTS.R')
load("HeartSteps_Experiments/HeartStepsV1-main/data_files/hs_out.RData")

### Plot 1: Comparison of predicted future outcomes across users
### with non-zero steps at fixed time T

TT_arr <- seq(100, 200, 10)


for(t_ind in seq_along(TT_arr)){

TT <- TT_arr[t_ind]
x0 = hs_out[[t_ind]]$Y_nz # Test data point
x1 = hs_out[[t_ind]]$Y_nz_focus # Prediction of FOCUS
x2 = hs_out[[t_ind]]$Y_nz_ms # Prediction of 

n <- length(x0)
df <- data.frame(
  User = 1:n,
  diff_x1 = x1 - x0,
  diff_x2 = x2 - x0
)

df_long <- df %>%
  pivot_longer(cols = c(diff_x1, diff_x2),
               names_to = "Method", values_to = "Error")

# Example: custom names for the x-axis
custom_labels <- hs_out[[t_ind]]$nz  

plt1 <- ggplot(df_long, aes(x = User, y = Error, color = Method, shape = Method)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 1.0) +
  scale_color_manual(
    values = c("diff_x2" = "#D55E00", "diff_x1" = "#0072B2"),  # mSSA red, FOCUS blue
    labels = c("diff_x2" = "mSSA", "diff_x1" = "FOCUS"),
    breaks = c("diff_x2", "diff_x1")
  ) +
  scale_shape_manual(
    values = c("diff_x2" = 16, "diff_x1" = 17),  # square for mSSA, triangle for FOCUS
    labels = c("diff_x2" = "mSSA", "diff_x1" = "FOCUS"),
    breaks = c("diff_x2", "diff_x1")
  ) +
  scale_x_continuous(
    breaks = 1:n,
    labels = custom_labels
  ) + 
  # scale_y_continuous(
  #   limits = c(-6, 2)
  # ) +
  labs(
    x = "Users with Positive Steps",
    y = expression("Prediction Error (" * hat(Y) - Y * ")"),
    title = sprintf("T = %d", TT_arr[t_ind]),
    color = "Series",
    shape = "Series"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "plain"),
    axis.title.x = element_text(size = 20, face = "plain"),
    axis.title.y = element_text(size = 20, face = "plain"),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    legend.title = element_blank(),
    legend.position = c(0.70, 0.20),
    legend.justification = c(0, 1),
    legend.text = element_text(size = 16),
    legend.background = element_rect(fill = scales::alpha("white", 0.6), color = NA), # semi-transparent
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray85", linewidth = 0.3),
    axis.line = element_line(color = "black"),
    axis.line.y.right = element_blank(),
    axis.line.x.top   = element_blank()
  )

plt1_name <- sprintf("HeartSteps_Experiments/HeartSteps_plots/Steps_prediction_%d.pdf", TT)

ggsave(plt1_name, plt1, width = 6, height = 6)

}


### Plot 2: Comparison of MSPE across T

MS_focus_arr = MS_ms_arr = numeric(length(TT_arr))
for(t_ind in seq_along(TT_arr)){
  MS_focus_arr[t_ind] <- sum(hs_out[[t_ind]]$MS_nz_focus)
  MS_ms_arr[t_ind] <- sum(hs_out[[t_ind]]$MS_nz_ms)
}

###  Plot the relative errors with the benchmark mSSA
df <- data.frame(
  TT    = TT_arr,
  Diff = MS_focus_arr - MS_ms_arr
)

plt2 <- ggplot(df, aes(x = TT, y = Diff)) +
  geom_hline(yintercept = 0, color = "red3", linetype = "dashed", linewidth = 1.0) +
  geom_line(color = "#0072B2", linetype = "solid", linewidth = 1, alpha = 0.5) +
  geom_point(color = "#0072B2", shape = 15, size = 4) +
  scale_x_continuous(
    breaks = seq(100, 200, by = 20),
    minor_breaks = seq(100, 200, by = 10)
  ) +
  scale_y_continuous(
    limits = c(-20, 5)
  ) +
  labs(
    x = "T",
    y = expression("Diff of MSPE (Positive Steps)"),
    title = "FOCUS - mSSA"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 24, hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    axis.line = element_line(color = "black"),
    panel.grid.minor.x = element_line(color = "grey85", linewidth = 0.5),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.5),
    panel.grid.minor.y = element_blank()
  )

ggsave("HeartSteps_Experiments/HeartSteps_plots/hs_predictions_vs_T.pdf", plot = plt2, width = 6, height = 6)



