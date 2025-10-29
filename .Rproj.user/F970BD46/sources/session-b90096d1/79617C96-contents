source('library_causalTS.R')
library(dplyr)
library(tidyverse)
library(hms)
library(ggplot2)
library(superheat)
library(gridExtra)
library(tseries)
library(forecast)
library(astsa)
library(lubridate)

# Reading the cleaned HeartSteps data set

hs <- read.csv("./HeartSteps_Experiments/HeartStepsV1-main/data_files/cleaned_data/suggestions_select.csv")
hs <- hs %>%  group_by(user.index)
hs_list <- hs %>% group_split()
hs_final_list <- list()
n_users <- hs %>% n_groups()
## max no of days = 300
# total_step_mat <- array(dim = c(n_users, 100))
## Set the user index and perform user-specific analysis
# variation_diff = numeric(user_ind)

for(user_ind in 1:n_users){
  # user_ind <- 10
  hs_user <- hs_list[[user_ind]] %>% 
    filter(sugg.select.slot %in% c(1:5))
  s <- hs_user$sugg.select.slot
  d <- as.Date(hs_user$sugg.select.utime)
  k <- nrow(hs_user)
  td <- difftime(hs_user$sugg.select.utime[2:k], hs_user$sugg.select.utime[1:(k-1)])
  z <- numeric(k) # indicator for keeping the date change
  
  # Create the date change indicator according to the 12 hour gap criterion: see notes
  z[1] = 0
  for(i in 1:(k-1)){
    if(s[i+1] > s[i]){
      z[i+1] <- ifelse(td[i] > 12, 1, 0)
    } else if(s[i+1] == s[i]){
      z[i+1] <- ifelse(td[i] > 12, 1, -1)
    } else{
      z[i+1] <- 1
    }
  }
  
  # Create the final data set tailored for pairwise slot analysis
  hs_final <- hs_user %>%
    dplyr::mutate(
      trt = send,
      date_grp = 1+ cumsum(z),
      jb10 = log(1 + jbsteps10),
      jb30 = log(1 + jbsteps30),
      jb40 = log(1 + jbsteps40),
      jb60 = log(1 + jbsteps60),
      jb90 = log(1 + jbsteps90),
      jb120 = log(1 + jbsteps120),
      s0 = sugg.select.slot,
      dt0 = sugg.select.utime,
      dt = sugg.decision.utime) %>%
    dplyr::select(dt0, 
                  s0, 
                  dt, 
                  date_grp, 
                  trt, 
                  jb10,
                  jb30, 
                  jb40,
                  jb60,
                  jb90,
                  jb120) %>% pivot_wider(
      id_cols = date_grp,
      names_from = s0,
      values_from = c(dt, 
                      trt, 
                      jb10,
                      jb30,
                      jb40,
                      jb60,
                      jb90,
                      jb120),
      names_sep = "_",
      values_fill = NA
    ) 
  
  # order the columns according to the names
  hs_names <- c("date_grp",
                paste0(rep(c("dt_", 
                             "trt_", 
                             "jb10_",
                             "jb30_",
                             "jb40_",
                             "jb60_",
                             "jb90_",
                             "jb120_"), 
                           each = 5), 1:5))
  hs_final <-hs_final[, hs_names]
  # # include the cleaned data in a list
  hs_final_list[[user_ind]] <- hs_final
}

### Create the full W and Y matrices

W_list <- list()
Y_list <- list()

for(user in 1:n_users){
  hs_user <- hs_final_list[[user]]
  trt_names <- paste0("trt_",1:5)
  step_names <- paste0("jb30_",1:5)
  W_list[[user]] <- c(t(hs_user[,trt_names]))
  Y_list[[user]] <- c(t(hs_user[,step_names]))
}

max_length <- max(sapply(W_list, length))
W_list_padded <- lapply(W_list, function(x) c(x, rep(NA, max_length - length(x))))
Y_list_padded <- lapply(Y_list, function(x) c(x, rep(NA, max_length - length(x))))

W_full <- do.call(rbind, W_list_padded)
Y_full <- do.call(rbind, Y_list_padded)


trt_prop <- numeric(nrow(W_full))

for(i in seq_along(trt_prop)){
  w <- W_full[i,] + 0
  sum(1 - w, na.rm = TRUE) / (length(w) - sum(is.na(w)))
}




### Pair-wise comparison of consecutive slots for a fixed user

diff_vec <- numeric(n_users)
diff_vec_trt <- numeric(n_users)
plt_list <- list()
for(user in 1:n_users){
  hs_user <- hs_final_list[[user]]
  x = hs_user$jb30_1
  y = hs_user$jb30_2
  z = hs_user$trt_2
  
  t1 <- prop.table(table(x[!z]!=0, y[!z]!=0), 1)
  if(ncol(t1) > 1){
    p <- t1[,2]
    diff_vec[user] <- ifelse(length(p)==2, diff(p), NA)
  } else {
    diff_vec[user] <- NA
  }
  
  t2 <- prop.table(table(x[z]!=0, y[z]!=0), 1)
  if(ncol(t2) > 1){
    p <- t2[,2]
    diff_vec_trt[user] <- ifelse(length(p)==2, diff(p), NA)
  } else {
    diff_vec_trt[user] <- NA
  }
  
  # Plot the pairwise outcomes for a fixed user
  # x_jitter <- ifelse(x == 0 & y == 0, jitter(x, amount = 0.1), x)
  # y_jitter <- ifelse(x == 0 & y == 0, jitter(y, amount = 0.1), y)
  # plot(x_jitter, y_jitter, col = z+1)
  
  # plt_data <- data.frame(x = x, y = y, z = z)
  # 
  # 
  # # Plot
  # 
  # plt <- ggplot(data = plt_data, aes(x = x, y = y, color = z)) +
  #   
  #   scale_color_discrete(na.translate = FALSE,
  #                        labels = c("Nudge = 0", "Nudge = 1")) +
  #   
  #   geom_point(size = 4, shape = 16, alpha = 0.6, stroke = 1.5) +
  #   
  #   labs(title = paste0("User ", user_ind),
  #        x = "Slot 2",
  #        y = "Slot 3") + labs(color = " ") +
  #   
  #   theme_minimal() +
  #   
  #   theme(# text=element_text(family="CM Roman", face="bold", size=12),
  #     plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),      # Plot title font size
  #     axis.title.x = element_text(size = 16),   # X-axis label font size
  #     axis.title.y = element_text(size = 16),   # Y-axis label font size
  #     axis.text = element_text(size = 12),      # Axis tick label font size
  #     axis.line = element_line(color = "black"),# Axis color
  #     legend.title = element_text(size = 14),   # Legend title font size
  #     legend.text = element_text(size = 12),
  #     legend.position = c(0.5, 1.1),  # Position the legend above the plot (relative coordinates)
  #     legend.justification = c(1.8, 1.1),  # Center the legend horizontally
  #     plot.margin = unit(c(1, 1, 1, 1), "cm")  # Set margins to 2 cm on all sides
  #   )
  # plt
  #
  # ggsave(filename = paste0("user", user_ind, "_3_4.pdf"),
  #        plot = plt,
  #        path = paste0("./HeartSteps_Experiments/Figures/HeartSteps_TimePlot/slot_3_4"),
  #        width=7, height=7)
}
density1 <- density(diff_vec, na.rm = T)
density2 <- density(diff_vec_trt, na.rm = T)

### Base plot version of the histogram
# hist(diff_vec, breaks = 10, probability = TRUE, 
#      col = rgb(0, 0, 1, 0.3), border = "white",
#      xlim = c(-1,1),
#      xlab = "Value", main = "Histogram with Density Overlay")
# 
# lines(density1, col = "blue", lwd = 2)
# 
# hist(diff_vec_trt, breaks = 10, probability = TRUE, 
#      col = rgb(1, 0, 0, 0.3), border = "white", add = T)
# lines(density2, col = "red", lwd = 2)

data_control <- data.frame(Value = diff_vec, Group = "Control")
data_treatment <- data.frame(Value = diff_vec_trt, Group = "Treatment")

# Combine data
data_combined <- rbind(data_control, data_treatment)

# Create the plot
ggplot(data_combined, aes(x = Value)) +
  # Histogram for both groups
  geom_histogram(data = subset(data_combined, Group == "Control"),
                 aes(y = ..density.., fill = "Control"), 
                 bins = 15, color = NA, alpha = 0.3) +
  geom_histogram(data = subset(data_combined, Group == "Treatment"),
                 aes(y = ..density.., fill = "Treatment"), 
                 bins = 15, color = NA, alpha = 0.3) +
  # Density lines
  geom_density(data = subset(data_combined, Group == "Control"),
               aes(color = "Control"), size = 1.2) +
  geom_density(data = subset(data_combined, Group == "Treatment"),
               aes(color = "Treatment"), size = 1.2) +
  # Custom axis lines
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"), # Add axis lines, 
    panel.spacing = unit(1, "lines"),  # Adjust spacing between panels
    plot.margin = margin(30, 30, 30, 30),  # Increase plot margins
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 12),     # Increase legend text size
    legend.title = element_text(size = 14),
    legend.position = c(0.85, 0.85),           # Place legend in top-right corner
    legend.spacing = unit(0.2, "cm"),         # Reduce spacing between items
    legend.key.size = unit(1, "lines")
  ) +
  # Customize scales and labels
  scale_fill_manual(values = c("Control" = rgb(0,0,1,0.5), "Treatment" = rgb(1,0,0,0.5))) +
  scale_color_manual(values = c("Control" = rgb(0,0,1,0.5), "Treatment" = rgb(1,0,0,0.5))) +
  scale_x_continuous(limits = c(-1.0, 1.0)) + # Set x-axis limits +
  labs(
       x = "Value",
       y = "Density",
       fill = "Legend",
       color = "Legend") +
  # ggtitle(expression(k==3))
  ggtitle(expression("Histogram of " ~ Delta[i]^{(2)}))




# Only plot slot pairs (4,5) for all users, points are colour-coded by the users

out_45 <- c()

for(user in 1:n_users){
  hs_user <- hs_final_list[[user]]
  user_45 <- hs_user[, c('jb30_3','jb30_4', 'jb30_5', 'trt_5')]
  user_vec <- rep(user, nrow(user_45))
  out_45 <- rbind(out_45, cbind(user_vec, user_45))
}

x = out_45$jb30_3 + out_45$jb30_4
y = out_45$jb30_5
z = out_45$trt_5
x_jitter <- ifelse(x == 0 & y == 0, jitter(x, amount = 0.1), x)
y_jitter <- ifelse(x == 0 & y == 0, jitter(y, amount = 0.1), y)
plot(x_jitter, y_jitter, col = z)

prop.table(table(x[z]!=0, y[z]!=0), margin = 1)



# box plot of the outcomes across the slots of the day
par(mfrow = c(10, 4))
for(user_ind in 1:n_users){
  astsa::tsplot(total_step_mat[user_ind,1:50], 
                col = "pink3", lwd = 2, ylim = c(0, 10),
                ylab = "jbsteps30", main = paste0("User ", user_ind))
  }

# Effect of slots on the outcomes


boxplot(log(1+jbsteps30) ~ sugg.decision.slot, 
        data = hs_subset, 
        subset = which(sugg.decision.slot %in% c(1:5)),
        xlab = "Decision Slot",
        ylab = "log(1 + jbsteps30)",
        names = c("10:00 am", "3:30 pm", "6:30 pm", "10:00 pm", "11:30 pm"),
        col = c("red", "skyblue3", "darkgreen", "purple", "darkorange"),
        cex.lab = 2, cex.axis = 1.5,
        main = "User 7", cex.main = 2)

boxplot(log(1+jbsteps30) ~ sugg.decision.slot, 
        data = hs_subset, 
        subset = which(sugg.decision.slot %in% c(1:5)),
        xlab = "Decision Slot",
        ylab = "log(1 + jbsteps30)",
        names = c("11:30 am", "4:00 pm", "6:30 pm", "9:50 pm", "11:30 pm"),
        col = c("red", "skyblue3", "darkgreen", "purple", "darkorange"),
        cex.lab = 2, cex.axis = 1.5, 
        main = "User 11", cex.main = 2)



# State space transition plots

prop_mat <- array(dim = c(n_users, 5)) # number of slots/day = 5
plot_list <- list()

for(user_ind in 1:9){
  hs_user <- hs_list[[user_ind]]
  x <- log(1+hs_user$jbsteps30)
  s <- hs_user$sugg.select.slot
  
  prop_mat[user_ind, ] <- prop.table(table(x!=0,s), 2)[2,]
  
  data_plot <- data.frame(s = as.factor(s), x = x)
  plot_list[[user_ind]] <- ggplot(data_plot, aes(x = s, y = x)) +
      geom_boxplot(aes(fill = s), na.rm = FALSE, orientation = "x") +
      
      theme(legend.position = "none",
            plot.margin = margin(10, 10, 10, 10),  # Increase plot margins
            plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
            axis.text.x = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            strip.text = element_text(size = 14)) +
      labs(x = "Slot", y = "log(1 + jbsteps30)", title = paste0("User ", user_ind))
}
combined_plot <- patchwork::wrap_plots(plot_list, ncol = 3)  # Adjust ncol to control the layout
print(combined_plot)


# Convert matrix to a data frame for ggplot2
df <- as.data.frame(prop_mat)
df$row <- 1:nrow(df)  # Add row indices
df_long <- pivot_longer(df, cols = -row, names_to = "Column", values_to = "Value")
df_long$Column <- as.factor(gsub("V", "", df_long$Column))  # Convert column names to factor
levels(df_long$Column) <- paste("k = ", 1:5)
math_labels <- setNames(lapply(1:5, function(i) bquote(s[t] == .(i))), as.character(1:5))

# Plot the histograms
ggplot(df_long, aes(x = Value, fill = Column)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.05, color = "black", alpha = 0.7) +  # Histogram with fixed binwidth
  facet_grid(. ~ Column, scales = "fixed", labeller = labeller(Column = math_labels)) +
  # facet_grid(. ~ Column, scales = "fixed") +  # Facet columns vertically and fix scales
  scale_fill_brewer(palette = "Set2") +       # Add color palette
  # theme_minimal() +
  theme(
    legend.position = "none",
    panel.spacing = unit(1, "lines"),  # Adjust spacing between panels
    plot.margin = margin(30, 30, 30, 30),  # Increase plot margins
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text.x = element_blank(),  # Center the title,
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text = element_text(size = 14)
  ) +
  labs(
    title = expression("Histograms of estimated " ~ P(Y[t] > 0 ~ "|" ~ s[t] == k)),
    x = "Values",
    y = " Density"
  ) +
  coord_flip()


# Peparation of boxplots
plot_list <- lapply(hs_list, function(data) {
  ggplot(data, aes(x = sugg.select.slot, y = log(1+jbsteps30))) +
    geom_boxplot(aes(fill = sugg.select.slot)) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(x = "Category", y = "Value", title = "Boxplot of x ~ s")
})


# Combine the plots into an array/grid
combined_plot <- wrap_plots(plot_list, ncol = 2)  # Adjust ncol to control the layout
print(combined_plot)


# Slot-slot scatterplots


w <- numeric(n_users)
plot_list <- list()
for(user_ind in 1:n_users){
  hs_user <- hs_list[[user_ind]]
  
  TT <- nrow(hs_user)
  k <- 1
  
  x <- log(1+hs_user$jbsteps30[1:(TT-k)])
  y <- log(1+hs_user$jbsteps30[(k+1):TT])
  z <- hs_user$sugg.select.slot[(k+1):TT]
  w[user_ind] = diff(prop.table(table(x!=0, y!=0), 1)[,2])
  # set.seed(123)
  # plt_data <- data.frame(x = x, y = y, z = z)
  # plt <- ggplot(data = plt_data, aes(x = x, y = y, color = z)) +
  #   
  #   scale_color_discrete(na.translate = FALSE,
  #                        labels = c("Nudge = 0", "Nudge = 1")) +
  #   
  #   geom_point(size = 4, shape = 16, alpha = 0.6, stroke = 1.5) +
  #   geom_jitter(
  #     data = plt_data[plt_data$x == 0 | plt_data$y == 0, ], # Subset for (0, 0) points
  #     aes(x = x, y = y),
  #     width = 0.1, height = 0.1, # Jitter amount
  #     size = 4, shape = 16, alpha = 0.6, stroke = 1.5
  #   ) +
  #   labs(title = paste0("User ", user_ind, ", Nudge proportion = ", round(mean(z), 2)),
  #        x = expression(Y[t-1]),
  #        y = expression(Y[t])) + labs(color = " ") +
  #   
  #   theme_minimal() +
  #   
  #   theme(# text=element_text(family="CM Roman", face="bold", size=12),
  #     plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),      # Plot title font size
  #     axis.title.x = element_text(size = 20),   # X-axis label font size
  #     axis.title.y = element_text(size = 20),   # Y-axis label font size
  #     axis.text = element_text(size = 12),      # Axis tick label font size
  #     axis.line = element_line(color = "black"),# Axis color
  #     legend.title = element_text(size = 14),   # Legend title font size
  #     legend.text = element_text(size = 12),
  #     legend.position = c(0.5, 1.1),  # Position the legend above the plot (relative coordinates)
  #     legend.justification = c(1.8, 1.1),  # Center the legend horizontally
  #     plot.margin = unit(c(1, 1, 1, 1), "cm")  # Set margins to 2 cm on all sides
  #   )
  # 
  # plot_list[[user_ind]] <- plt
  # ggsave(filename = paste0("user", user_ind, "_lag1.pdf"),
  #        plot = plt,
  #        path = paste0("./HeartSteps_Experiments/Figures/HeartSteps_TimePlot/lag1_plots"),
  #        width=7, height=7)
}


w_data = as.data.frame(w)
ggplot(w_data, aes(x=w)) + 
  geom_histogram(binwidth = 0.05, aes(y=..density..), color = "gray90",
                 fill="steelblue1", alpha = 0.6)+
  labs(title = " ",
       x = "Proportion difference",
       y = "Density")+
  theme(# text=element_text(family="CM Roman", face="bold", size=12),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),      # Plot title font size
    axis.title.x = element_text(size = 20),   # X-axis label font size
    axis.title.y = element_text(size = 20),   # Y-axis label font size
    axis.text = element_text(size = 12),      # Axis tick label font size
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # Set margins to 2 cm on all sides
  )


ggplot(w_data, aes(x = 1:37, y = w)) + 
  geom_hline(yintercept = 0, color = "pink", linetype = "dashed", size = 1) +
  geom_point(size = 3, col = 'red3', alpha = 0.5) +
  theme_minimal() +
  labs(x = "User", y = "Proportion difference") +
  theme(# text=element_text(family="CM Roman", face="bold", size=12),
  plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),      # Plot title font size
  axis.title.x = element_text(size = 20),   # X-axis label font size
  axis.title.y = element_text(size = 20),   # Y-axis label font size
  axis.text = element_text(size = 14),      # Axis tick label font size
  plot.margin = unit(c(1, 1, 1, 1), "cm")  # Set margins to 2 cm on all sides
)

plot_list2 = plot_list[1:9]

combined_plot <- patchwork::wrap_plots(plot_list2, ncol = 3)  # Adjust ncol to control the layout
print(combined_plot)











plt_data <- data.frame(x = x, y = y, z = z)
# Plot


prop_nudge <- array(dim = c(n_users, 2))
for(user_ind in 1:n_users){

hs_user <- hs_list[[user_ind]]
hs_final <- hs_final_list[[user_ind]]

x = log(1+hs_user$jbsteps30)[1:(TT-1)]
y = log(1+hs_user$jbsteps30)[2:TT]
z = hs_user$send[2:TT]
p_tab = prop.table(table(interaction(z, x!=0), y!=0), 1)[,2]
prop_nudge[user_ind, ] = c(p_tab[2] - p_tab[1], p_tab[4] - p_tab[3])
}

df <- data.frame(
  value = c(prop_nudge[, 1], prop_nudge[, 2]),  # Combine both columns
  column = rep(c("Column 1", "Column 2"), each = 37)  # Add a column identifier
)
ggplot(df, aes(x = value, fill = column)) +
  geom_histogram(position = "identity", alpha = 0.4, binwidth = 0.05, color = "black") +
  labs(title = "", 
       x = "Value",
       y = "Frequency",
       fill = "") +
  scale_fill_manual(
    values = c("blue", "red"),          # Custom colors
    labels = c(expression(Y[t-1]~"="~0), expression(Y[t-1]~">"~0))  # Custom legend labels
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 14),
    panel.spacing = unit(1, "lines"),  # Adjust spacing between panels
    plot.margin = margin(30, 30, 30, 30),  # Increase plot margins
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text = element_text(size = 14)
  )



TT <- nrow(hs_user)
steps <- log(1+hs_user$jbsteps30)
s <- hs_user$sugg.select.slot

x <- steps[s==5]
y <- steps[s==2]
s <- s[2:TT]
# View(cbind(x, y, s1))

plot(x[s!=1], y[s!=1], pch = 20)
prop.table(table(x[s!=1]!=0, y[s!=1]!=0) , 1)
prop.table(table(x!=0, y!=0) , 1)


plt <- ggplot(data = plt_data, aes(x = x, y = y, color = z)) +
  
  scale_color_discrete(na.translate = FALSE,
                       labels = c("Nudge = 0", "Nudge = 1")) +
  
  geom_point(size = 4, shape = 16, alpha = 0.6, stroke = 1.5) +
  geom_jitter(
    data = plt_data[plt_data$x == 0 | plt_data$y == 0, ], # Subset for (0, 0) points
    aes(x = x, y = y),
    width = 0.3, height = 0.3, # Jitter amount
    size = 4, shape = 16, alpha = 0.6, stroke = 1.5
  ) +
  labs(title = paste0("User ", user_ind, ", lag = ", k, ", Nudge proportion = ", round(mean(z), 2)),
       x = expression(y[t-lag]),
       y = expression(y[t])) + labs(color = " ") +
  
  theme_minimal() +
  
  theme(# text=element_text(family="CM Roman", face="bold", size=12),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),      # Plot title font size
    axis.title.x = element_text(size = 16),   # X-axis label font size
    axis.title.y = element_text(size = 16),   # Y-axis label font size
    axis.text = element_text(size = 12),      # Axis tick label font size
    axis.line = element_line(color = "black"),# Axis color
    legend.title = element_text(size = 14),   # Legend title font size
    legend.text = element_text(size = 12),
    legend.position = c(0.5, 1.1),  # Position the legend above the plot (relative coordinates)
    legend.justification = c(1.8, 1.1),  # Center the legend horizontally
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # Set margins to 2 cm on all sides
  )
plt
