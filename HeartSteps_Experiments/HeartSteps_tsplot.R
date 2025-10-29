# Continued from HeartSteps_Clean.R
# library(tikzDevice)
# Check: https://cran.r-project.org/web/packages/tikzDevice/index.html
hs_list <- hs %>% group_split()
n_users <- hs %>% n_groups()

parent_dir <- "./HeartSteps_Experiments/Figures/HeartSteps_TimePlot"  # replace with your path

# Create folders "User_1" to "User_20"
for (user_ind in 1:n_users) {
  dir_name <- sprintf("User%d", user_ind)  # Format the folder name
  dir.create(file.path(parent_dir, dir_name))  # Create the folder
}


user_ind = 10
hs_user = hs_list[[user_ind]]
jb30_user <- log(1+hs_user$jbsteps30)

date_time_user <- hs_user$sugg.decision.utime
date_user <- as.Date(date_time_user)
time_user <- strftime(date_time_user, format = "%H:%M:%S")
trt_user <- hs_user$send

jb30_trt <- jb30_user[which(trt_user)]
jb30_ctrl <- jb30_user[which(!trt_user)]

datetime_trt <- date_time_user[trt_user]
datetime_ctrl <- date_time_user[!trt_user]

plot_data <- data.frame(
  date_time = date_time_user,
  jb30 = jb30_user,
  trt = trt_user
)

plot_data <- plot_data[5:34,]%>%
  mutate(date_time = as.POSIXct(date_time))

# Extract 12am of each day
midnights <- as.POSIXct(unique(as.Date(plot_data$date_time)))

# Plot
ggplot(plot_data, aes(x = date_time, y = jb30)) +
  # Solid vertical lines for midnight
  geom_vline(xintercept = as.numeric(midnights),
             linetype = "solid", color = "blue", linewidth = 0.5) +
  
  # Dashed lines for each time point
  geom_vline(xintercept = as.numeric(plot_data$date_time),
             linetype = "dashed", color = "skyblue1", linewidth = 0.5) +
  
  # Scatter points colored by treatment
  geom_point(aes(color = factor(trt)), size = 3, alpha = 0.8) +
  
  # Connecting lines
  geom_line(alpha = 0.4, color = "purple", linewidth = 0.8) +
  
  # Set x-axis ticks at midday
  scale_x_datetime(breaks = midday_labels, date_labels = "%b %d") +
  
  # Remove legend
  scale_color_manual(values = c("green3", "red3")) +
  
  labs(
    title = "",
    x = " ",
    y = "log (Steps)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 2, hjust = 1, size = 11),
    axis.title = element_text(size = 18),
    legend.position = "none"
  )




  
check = 0
for(user_ind in c(1:n_users)){
  hs_user = hs_list[[user_ind]]
  jb30_user <- log(1+hs_user$jbsteps30)
  
  date_time_user <- hs_user$sugg.decision.utime
  date_user <- as.Date(date_time_user)
  time_user <- strftime(date_time_user, format = "%H:%M:%S")
  trt_user <- hs_user$send
  
  jb30_trt <- jb30_user[which(trt_user)]
  jb30_ctrl <- jb30_user[which(!trt_user)]
  
  datetime_trt <- date_time_user[trt_user]
  datetime_ctrl <- date_time_user[!trt_user]
  
  
  # The codes for plotting the stepcounts
  
  # tikz("demo_tsplot.tex", standAlone = TRUE)
  
  
  #### EXPERIMENT
  # TT = length(jb30_user)
  # plot(1, type = 'n', xlim = c(0, 300), ylim = c(0, 9),
  #      ylab = 'log steps')
  # points(jb30_user, col = trt_user + 1, pch = 16)
  # lines(jb30_user, lty = 1, lwd = 1, col = 'skyblue')
  # unique_dates <- c(1, diff(date_user))
  # abline(v = which(unique_dates == 1), 
  #        lty = 1, 
  #        col = 'blue')
  # 
  # 
  # x = date_time_user
  # y = jb30_user
  # z = trt_user
  # unique_dates <- unique(format(as.Date(x), "%Y-%m-%d"))
  # unique_months <- unique(format(as.Date(x), "%m-%d"))
  # 
  # plot(1, type = "n",
  #      xaxt = "n",
  #      xlab = "",
  #      ylab = "",
  #      xlim = range(x, na.rm = TRUE),
  #      ylim = range(y, na.rm = TRUE),
  #      main = paste0("User ", user_ind))
  # 
  # abline(v = as.POSIXct(unique_dates), 
  #        lty = 1, 
  #        col = rgb(0,0,1,0.6))
  # axis(1, at = x, labels = FALSE)
  # abline(v = x, 
  #        lty = 3,
  #        col = 'skyblue1')
  # points(x, y,
  #        col = case_when(z == TRUE ~ "green3", 
  #                        z == FALSE ~ 'red3'),
  #        pch = case_when(z == TRUE ~ 15, 
  #                        z == FALSE ~ 16),
  #        cex = 1)
  # 
  # lines(x, y,
  #       col = rgb(0.5, 0, 0.5, alpha = 0.4),
  #       lty = 1,
  #       lwd = 2)
  # 
  # text(x = as.POSIXct(paste(unique_dates, "12:00:00")),
  #      y = min(y, na.rm = TRUE) - 0.7,
  #      labels = unique_months,
  #      srt = 90, 
  #      adj = 1, 
  #      xpd = TRUE,
  #      cex = 0.6)
  # legend(max(x) + as.difftime(0.1, units = "hours"),
  #        5,
  #        legend = c("Trt", "Ctrl"),
  #        col = c("green3", "red3"),
  #        pch = c(15, 16, 17),
  #        bg = "transparent",
  #        bty = "n",
  #        cex = 1,
  #        xpd = TRUE)
  # 
  # 
  # 
  # 
  
  
  plot_data <- data.frame(
    date_time = date_time_user,
    jb30 = jb30_user,
    trt = trt_user
  )
  
  plot_data <- plot_data[1:20,]
  
  # Convert dates for unique dates and labels
  unique_dates <- as.Date(unique(format(as.Date(plot_data$date_time), "%Y-%m-%d")))
  unique_months <- unique(format(as.Date(plot_data$date_time), "%m-%d"))
  
  # Base plot
  p <- ggplot(plot_data, aes(x = date_time, y = jb30)) +
    geom_vline(xintercept = as.numeric(as.POSIXct(unique_dates)), 
               linetype = "solid", color = rgb(0, 0, 1, 0.6), linewidth = 0.5) +
    # Add vertical dashed lines for all x points
    geom_vline(xintercept = as.numeric(plot_data$date_time), 
               linetype = "dotted", color = 'skyblue1', linewidth = 0.5) +
    
    # Add points with color and shape based on treatment variable
    geom_point(aes(color = factor(trt), shape = factor(trt)), size = 2) +
    # Add lines
    geom_line(color = rgb(0.5, 0, 0.5, alpha = 0.4), size = 1) +
    
    # Customize scales and labels
    scale_color_manual(values = c("green3", "red3"), labels = c("Trt", "Ctrl")) +
    scale_shape_manual(values = c(15, 16), labels = c("Trt", "Ctrl")) +
    labs(
      title = paste0("User ", user_ind),
      x = NULL,
      y = "log steps",
      color = " ",
      shape = " "
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      axis.title.x = element_text(size = 20),   # X-axis label font size
      axis.title.y = element_text(size = 20),   # Y-axis label font size
      axis.text.x = element_blank(),
      legend.text = element_text(size = 12),  # Increase legend text size
      axis.ticks.x = element_blank(),
      plot.margin = margin(10, 10, 10, 10),  # Further increase bottom margin
      plot.caption = element_text(hjust = 0, vjust = 1)
    ) +
    coord_cartesian(clip = "off") + 
    geom_text(data = data.frame(x = as.POSIXct(paste(unique_dates, "12:00:00")), 
                                y = min(plot_data$jb30, na.rm = TRUE) - 1, 
                                label = unique_months),
              aes(x = x, y = y, label = label),
              angle = 90, hjust = 1, size = 3) +
    theme(
      plot.margin = margin(10, 10, 50, 10)  # Ensure enough space for axis labels
    )
  
  # Add legend outside the plot
  p <- p +
    theme(legend.position = "right")
  
  # Print plot
  # print(p)
  
  
  ggsave(filename = paste0("User", user_ind, "_jbsteps30_full.pdf"),
         plot = p,
         path = paste0("./HeartSteps_Experiments/Figures/HeartSteps_TimePlot/full_plots"),
         # path = paste0("./HeartSteps_Experiments/Figures/HeartSteps_TimePlot/User", user_ind),
         width=25, height=3)
  ### EXPERIMENT
}

cat(paste("Plot complete for", check, "many users."))












