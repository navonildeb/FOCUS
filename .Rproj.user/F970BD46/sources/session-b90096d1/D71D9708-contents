library(tidyverse)
library(hms)
library(superheat)
library(gridExtra)

### Cleaning HeartSteps

hs_full <- read_csv("./HeartSteps_Experiments/HeartStepsV1-main/data_files/suggestions.csv")
hs <- hs_full %>%
  group_by(user.index) %>%
  arrange(sugg.decision.utime, .by_group = TRUE) %>%
  dplyr::select(user.index,
                sugg.select.utime,
                sugg.select.slot,
                sugg.select.update,
                sugg.decision.utime, 
                sugg.decision.slot,
                avail, 
                send,
                jbsteps10,
                jbsteps30,
                jbsteps40,
                jbsteps60,
                jbsteps90,
                jbsteps120)

hs <- hs %>% 
  drop_na(user.index,
          sugg.select.utime,
          sugg.select.slot,
          sugg.select.update,
          sugg.decision.utime, 
          sugg.decision.slot) %>%
  filter(avail == TRUE) %>% dplyr::select(-avail)
hs_list <- hs %>% group_split()
n_users <- hs %>% n_groups()

write.csv(hs_full, 
          "./HeartSteps_Experiments/HeartStepsV1-main/data_files/cleaned_data/suggestions_full.csv",
          row.names = FALSE)
write.csv(hs, 
          "./HeartSteps_Experiments/HeartStepsV1-main/data_files/cleaned_data/suggestions_select.csv",
          row.names = FALSE)


# Create the full matrices of outcome and treatment variables arranged with the slots

hs_list <- hs %>% group_split()
hs_final_list <- list()
n_users <- hs %>% n_groups()

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
    dplyr::select(dt0, s0, dt, date_grp, trt, jb10, jb30, jb40,jb60,
                  jb90, jb120) %>% pivot_wider(
                    id_cols = date_grp,
                    names_from = s0,
                    values_from = c(dt, trt, jb10, jb30, jb40, jb60,
                                    jb90, jb120),
                    names_sep = "_",
                    values_fill = NA
                  ) 
  
  # order the columns according to the names
  hs_names <- c("date_grp",
                paste0(rep(c("dt_", "trt_", "jb10_", "jb30_", "jb40_",
                             "jb60_", "jb90_", "jb120_"), 
                           each = 5), 1:5))
  hs_final <-hs_final[, hs_names]
  ### include the cleaned data in a list
  hs_final_list[[user_ind]] <- hs_final
}

save(hs_final_list, file = "HeartSteps_Experiments/HeartStepsV1-main/data_files/cleaned_data/suggestions_select_list.RData")

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

write.csv(W_full, file = "HeartSteps_Experiments/HeartStepsV1-main/data_files/cleaned_data/W_full.csv")
write.csv(Y_full, file = "HeartSteps_Experiments/HeartStepsV1-main/data_files/cleaned_data/log_jb30_full.csv")
