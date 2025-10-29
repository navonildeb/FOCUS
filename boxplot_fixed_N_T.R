source('library_causalTS.R')

methods <- c("FOCUS", "mSSA", "SyN-BEATS")
TT_arr <- 2^c(5:8)
logT_arr <- 5:8
R.max <- 30
n_ind = 2
dgp_ind = 2
h = 3

load(sprintf("data_files/DGP%d/DGP%d_out.RData", dgp_ind, dgp_ind))

out_dgp <- switch (as.character(dgp_ind),
  "0" = out_dgp0,
  "1" = out_dgp1,
  "2" = out_dgp2
)

err_focus <- t(switch(as.character(h),
                    "1" = out_dgp$rep_focus_h1[n_ind, seq_along(TT_arr), ],
                    "2" = out_dgp$rep_focus_h2[n_ind, seq_along(TT_arr), ],
                    "3" = out_dgp$rep_focus_h3[n_ind, seq_along(TT_arr), ]
))

err_ms <- t(switch(as.character(h),
                 "1" = out_dgp$rep_ms_h1[n_ind, seq_along(TT_arr), ],
                 "2" = out_dgp$rep_ms_h2[n_ind, seq_along(TT_arr), ],
                 "3" = out_dgp$rep_ms_h3[n_ind, seq_along(TT_arr), ]
))

load(sprintf("data_files/DGP%d/synout_files/synout_DGP%d.RData", dgp_ind, dgp_ind))
err_syn <- t(switch(as.character(h),
                  "1" = syn_errors_final$rep_syn_h1,
                  "2" = syn_errors_final$rep_syn_h2,
                  "3" = syn_errors_final$rep_syn_h3
))



# Reshape into long format (raw differences only)
df_long <- lapply(seq_along(logT_arr), function(j) {
  tibble(
    rep = 1:nrow(err_focus),
    logT = logT_arr[j],
    diff_ms  = err_focus[, j] - err_ms[, j],
    diff_syn = err_focus[, j] - err_syn[, j]
  )
}) %>% bind_rows()

df_diff <- df_long %>%
  pivot_longer(c(diff_ms, diff_syn),
               names_to = "Comparison", values_to = "diff") %>%
  mutate(Comparison = recode(Comparison,
                             "diff_ms" = "FOCUS - mSSA",
                             "diff_syn" = "FOCUS - SyNBEATS"))

q10 <- df_diff %>%
  group_by(logT, Comparison) %>%
  summarise(p10 = quantile(diff, 0.1, na.rm = TRUE), .groups = "drop")

ymin <- min(q10$p10)

# --- Wilcoxon tests (same as your code, but build plotmath labels "ab.cd %*% 10^{x}")
pvals <- df_diff %>%
  group_by(logT, Comparison) %>%
  summarise(
    pval = wilcox.test(diff, rep(0, R.max), alternative = "less", paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    expo  = ifelse(pval > 0, floor(log10(pval)), NA_integer_),
    mant  = ifelse(pval > 0, pval / (10^expo), NA_real_),
    # plotmath string: e.g. "1.23 %*% 10^{-4}"
    label = ifelse(pval > 0,
                   sprintf("%.2f %%*%% 10^{%d}", mant, expo),
                   "0")
  )

pvals <- pvals %>%
  mutate(ypos = ifelse(Comparison == "FOCUS - mSSA",
                       -0.22 * ymin,   # higher
                       -0.15 * ymin))  # lower

# --- Plot (only change: parse = TRUE in geom_text)
box_plt <- ggplot(df_diff, aes(x = factor(logT), y = diff, fill = Comparison)) +
  geom_boxplot(position = position_dodge(width = 0.5), alpha = 0.6) +
  scale_fill_manual(values = c("FOCUS - mSSA" = "#D55E00",
                               "FOCUS - SyNBEATS" = "#009E73")) +
  scale_y_continuous(limits = c(ymin, -0.25*ymin)) +
  scale_x_discrete(
    labels = c("5" = expression(2^5),
               "6" = expression(2^6),
               "7" = expression(2^7),
               "8" = expression(2^8))
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_text(
    data = pvals,
    aes(x = factor(logT), y = ypos, group = Comparison, label = label),
    position = position_dodge(width = 0.5),
    vjust = 0,
    parse = TRUE     # <-- this makes "1.23 %*% 10^{-4}" render as 1.23 × 10^{-4}
  ) +
  labs(
    x = expression(T),
    y = "Diff. of Mean sq. Forecast Error",
    title = sprintf("DGP-%d, h = %d", dgp_ind+1, h)
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.10, 0.23),
    legend.justification = c(0, 1),
    legend.text = element_text(size = 13),
    legend.background = element_rect(fill = scales::alpha("white", 0.3), color = NA),
    legend.box.background = element_blank(),
    axis.title.x = element_text(size = 20, face = "plain"),
    axis.title.y = element_text(size = 20, face = "plain"),
    axis.text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5, size = 24, face = "plain"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray85", linewidth = 0.3),
    axis.line = element_line(color = "black"),
    axis.line.y.right = element_blank(),
    axis.line.x.top   = element_blank()
  )


box_plt

plt_name <- sprintf("figures_out/boxplots/boxplot_DGP%d_h%d.pdf", dgp_ind, h)
ggsave(plt_name, plot = box_plt, width = 6, height = 6)
