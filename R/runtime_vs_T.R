source("library_causalTS.R")
library(ggforce)

dgp_ind <- 2 # dgp_ind = 0,1,2,4
N.arr <- 2^c(5:9)
TT.arr <- 2^(5:9)

load(sprintf("data_files/DGP%d/DGP%d_out.RData", dgp_ind, dgp_ind))
load(sprintf("data_files/DGP%d/synout_files/synout_DGP%d.RData", dgp_ind, dgp_ind))

out_dgp <- switch(
  as.character(dgp_ind),
  "0" = out_dgp0,
  "1" = out_dgp1,
  "2" = out_dgp2,
  "4" = out_dgp4
)

n_ind <- 2
N <- N.arr[n_ind]

## Reshape runtime matrix
reshape_time <- function(time_mat, method_name, TT_vec = TT.arr) {
  
  time_mat <- drop(time_mat)
  
  ## If rows are reps and columns are T, transpose it
  if (length(dim(time_mat)) == 2 &&
      nrow(time_mat) != length(TT_vec) &&
      ncol(time_mat) == length(TT_vec)) {
    time_mat <- t(time_mat)
  }
  
  stopifnot(nrow(time_mat) == length(TT_vec))
  
  as.data.frame(time_mat) %>%
    mutate(TT = TT_vec) %>%
    pivot_longer(-TT, names_to = "rep", values_to = "time") %>%
    mutate(
      rep     = as.integer(gsub("V", "", rep)),
      method  = method_name,
      logT    = log2(TT),
      logtime = log2(time)
    ) %>%
    select(method, TT, logT, rep, time, logtime)
}

## FOCUS and mSSA use all T values
df_focus <- reshape_time(drop(out_dgp$timerep_focus[n_ind,,]), "FOCUS", TT.arr)
df_ms    <- reshape_time(drop(out_dgp$timerep_ms[n_ind,,]),    "mSSA",  TT.arr)

## SyNBEATS currently has only 4 T values: dim = 4 x 30
syn_time <- drop(syn_errors_final$timerep_syn)
TT_syn <- TT.arr[seq_len(nrow(syn_time))]

df_syn <- reshape_time(syn_time, "SyNBEATS", TT_syn)

df_all <- bind_rows(df_ms, df_syn, df_focus)

## enforce method order
df_all$method <- factor(
  df_all$method,
  levels = c("mSSA", "SyNBEATS", "FOCUS")
)

df_summary <- df_all %>%
  group_by(method, TT, logT) %>%
  summarise(
    center = median(time),
    spread = mad(time, constant = 1),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = pmax(center - spread, 0),
    ymax = center + spread
  )

df_summary$method <- factor(
  df_summary$method,
  levels = c("mSSA", "SyNBEATS", "FOCUS")
)

method_colors <- c(
  "mSSA"     = "#D55E00",
  "FOCUS"    = "#0072B2",
  "SyNBEATS" = "#009E73"
)

method_shapes <- c(
  "mSSA"     = 16,
  "FOCUS"    = 17,
  "SyNBEATS" = 18
)

method_lty <- c(
  "mSSA"     = 2,
  "FOCUS"    = 1,
  "SyNBEATS" = 3
)

plt_time <- ggplot(
  df_summary,
  aes(
    x = logT,
    y = center,
    color = method,
    shape = method
  )
) +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    width = 0.2,
    linewidth = 0.9,
    alpha = 0.5
  ) +
  geom_line(
    aes(group = method, linetype = method),
    linewidth = 1
  ) +
  geom_point(size = 5, alpha = 0.8) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  scale_linetype_manual(values = method_lty) +
  scale_x_continuous(
    breaks = log2(TT.arr),
    labels = as.expression(
      lapply(log2(TT.arr), function(x) bquote(2^{.(x)}))
    )
  ) +
  labs(
    x = expression(T),
    y = "Runtime (seconds)",
    title = sprintf("DGP-%d", dgp_ind)
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 18),
    plot.title = element_text(size = 22, hjust = 0.5),
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.55, 0.85),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.3),
      color = NA
    )
  )

plt_name <- sprintf("figures/runtime_plots_vs_T/runtime_DGP%d.pdf", dgp_ind)

ggsave(
  plt_name,
  plot = plt_time,
  width = 7,
  height = 4
)