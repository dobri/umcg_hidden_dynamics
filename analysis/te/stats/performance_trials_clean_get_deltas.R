library(janitor)
library(tidyverse)

save_files = TRUE

#----------------------------------------------
# Part A. Import, re-label, clean, conditions, % improvement
#----------------------------------------------

# setwd('~/logos/umcg_hidden_dynamics/analysis/te/stats')
# setwd('C:\\Users\\ddotov\\Downloads')

filename_scores_in <- 'Scores_TE_2024-05-26.csv'
x<-read.csv(filename_scores_in) %>%
  clean_names()
x$training_phase_label <- sub("_", "", x$training_phase_label)
x$training_phase_label <- sub("est_", "est", x$training_phase_label)
summary(x)


# Create a 'training condition' label in the long table
xcond <- x[x$training_phase_label=='Training',] %>%
  select(pp, task_code) %>%
  group_by(pp) %>%
  summarise(training_condition = mean(task_code)) %>%
  ungroup()
x <- x %>% left_join(., xcond, by = "pp")
x$training_condition <- as.factor(x$training_condition)
rm(xcond)

x$Training <- 'Periodic Fixed'
x$Training[x$training_condition==25] <- 'Chaotic Interactive'
x$Training[x$training_condition==30] <- 'Chaotic Non-interactive'
x$Training <- as.factor(x$Training)

x$task_code <- as.factor(x$task_code)
x$task_label <- 'Periodic Fixed'
x$task_label[x$task_code==12] <- 'Periodic Fadeout (SCT)'
x$task_label[x$task_code==20] <- 'Chaotic'
x$task_label[x$task_code==21] <- 'Chaotic Visual'
x$task_label[x$task_code==10] <- 'Periodic Interactive Training'
x$task_label[x$task_code==25] <- 'Chaotic Interactive Training'
x$task_label[x$task_code==30] <- 'Chaotic Non-Interactive Training'
x$task_label <- as.factor(x$task_label)

x$training_phase_label <- as.factor(x$training_phase_label)
x$training_phase_label <- relevel(x$training_phase_label, ref='PreTest')


# Pseudo-trial number vector
trial_counter = 1
x$trial = 0
x$trial[1] = trial_counter
for (n in seq(2,dim(x)[1])) {
  if ((x$pp[n]==x$pp[n-1]) && 
      (x$task_code[n]==x$task_code[n-1]) &&
      (x$training_phase[n]==x$training_phase[n-1])) {
    trial_counter = trial_counter + 1
  }
  else {
    trial_counter = 1
  }
  x$trial[n] <- trial_counter
}
x$pp <- factor(x$pp)


# Remove one single weird trial?
x <- x[!((x$training_phase == 'Retention') & (x$trial==6)),]


# Mean pre-test
xd <- x[x$training_phase_label=='PreTest',] %>%
  group_by(pp,task_code) %>%
  summarise(pre_score = mean(score),
            pre_c = mean(c),
            pre_tau = mean(tau),
            pre_pitch_error = mean(pitch_error),
            pre_testimtouser = mean(te12rescaled),
            pre_teusertostim = mean(te21rescaled)) %>%
  ungroup()
x <- x %>% left_join(., xd, by = c("pp","task_code"))
rm(xd)


# Rescale as %improvement
x$score_delta <- (x$score - x$pre_score)/x$pre_score
x$c_delta <- 100*(x$c - x$pre_c)/x$pre_c
x$tau_delta <- 100*(x$tau - x$pre_tau)/x$pre_tau
x$pitch_error_delta <- 100*(x$pitch_error - x$pre_pitch_error)/x$pre_pitch_error
x$testimtouser_delta <- 100*(x$te12rescaled - x$pre_testimtouser)/x$pre_testimtouser
x$teusertostim_delta <- 100*(x$te21rescaled - x$pre_teusertostim)/x$pre_teusertostim


# Save the cleaned and nicely labeled scores table.
if(save_files){
  filename = paste(filename_scores_in,'cleaned','.csv',sep='_')
  write.csv(x, file=filename, row.names=FALSE)
}

# End Part A


#----------------------------------------------
# Part B. Test data only, averaged per task
#----------------------------------------------

# rm(list=ls())
# setwd('~/logos/umcg_hidden_dynamics/analysis/te/stats')
# setwd('C:\\Users\\ddotov\\Downloads')
# filename = 'Scores_TE_2024-05-26.csv_cleaned_.csv'
# w/out the training conditions
x_test <- filter(x, training_phase_label != "Training")
# x_test <- filter(read.csv(filename,sep=','), training_phase != "Training")

x_test_ave <- x_test %>%
  group_by(pp,task_label,training_phase_label) %>%
  summarise(score = mean(score),
            c = mean(c),
            pitch_error = mean(pitch_error),
            te12rescaled = mean(te12rescaled),
            te21rescaled = mean(te21rescaled)) %>%
  ungroup() %>%
  # Add back columns of interest
  inner_join(select(x_test, pp, task_label, Training), relationship = "many-to-many") %>%
  # Get rid of duplicated rows
  distinct(pp, task_label, training_phase_label, .keep_all = TRUE)

if(save_files){
  filename_ave = paste(filename,'_ave','.csv',sep='_')
  write.csv(x_test_ave, file=filename_ave, row.names=FALSE)
}

x_test_ave_wide <- x_test_ave %>%
  pivot_wider(names_from = training_phase_label, values_from = c('score','c','pitch_error','te12rescaled','te21rescaled'))

if(save_files){
  filename_ave = paste(filename,'_ave_wide','.csv',sep='_')
  write.csv(x_test_ave_wide, file=filename_ave, row.names=FALSE)
}

# End Part B


#----------------------------------------------
# Part C. Test data only, delta scores averaged per task
#----------------------------------------------

# rm(list=ls())
# setwd('~/logos/umcg_hidden_dynamics/analysis/performance')
# setwd('C:\\Users\\ddotov\\Downloads')
# filename = 'Scores_TE_2024-05-26.csv_cleaned_.csv'
# w/out the training conditions
# x_test <- filter(read.csv(filename,sep=','), training_phase != "Training")
x_test <- filter(x, training_phase_label != "Training")


#----------------------------------------------
# Delta scores as percentage (per pp) of the average PreTest score 
# Each task_label, pp, and training_phase
# PreTest vs PostTest; PreTest vs Retention
# Widen the table such that Pre/Post/Retention are in their own columns
# and calculate difference scores between target columns.
#----------------------------------------------
fun_ave_pivot_diff <- function(x, y){
  x$dv <- y
  diff_data <- x %>%
    group_by(pp, training_phase_label, task_label) %>%
    summarise(mean_PerTask = mean(dv)) %>%
    ungroup() %>%
    pivot_wider(names_from = training_phase_label, values_from = mean_PerTask) %>%
    group_by(pp) %>%
    mutate(PostDelta = 100*(PostTest - PreTest)/PreTest, RetentionDelta = 100*(Retention - PreTest)/PreTest) %>% # Calc diff scores
    inner_join(select(x, pp, Training), relationship = "many-to-many") %>% # Add back columns of interest
    distinct(pp, task_label, .keep_all = TRUE) %>% # Get rid of duplicated rows
    pivot_longer(cols=c("PostDelta","RetentionDelta"), 
                 names_to = "TrainingPhase", values_to = "Delta")
  return (diff_data)
}


x_delta <- fun_ave_pivot_diff(x_test, x_test$score)
diff_data <- fun_ave_pivot_diff(x_test, x_test$c)
x_delta <- x_delta %>% left_join(., diff_data, by = c("pp","task_label","Training","TrainingPhase"))

diff_data <- fun_ave_pivot_diff(x_test, x_test$pitch_error)
x_delta <- x_delta %>% left_join(., diff_data, by = c("pp","task_label","Training","TrainingPhase"))

diff_data <- fun_ave_pivot_diff(x_test, x_test$te12rescaled)
x_delta <- x_delta %>% left_join(., diff_data, by = c("pp","task_label","Training","TrainingPhase"))

diff_data <- fun_ave_pivot_diff(x_test, x_test$te21rescaled)
x_delta <- x_delta %>% left_join(., diff_data, by = c("pp","task_label","Training","TrainingPhase"))
rm(diff_data)

# Sanity check
# head(x_delta, 20)

names(x_delta) <- c("pp","task_label",
                  "score_PreTest","score_PostTest","score_Retention",
                  "Training","TrainingPhase",
                  "score_Delta",
                  "c_PreTest","c_PostTest","c_Retention","c_Delta",
                  "pitch_error_PreTest","pitch_error_PostTest","pitch_error_Retention","pitch_error_Delta",
                  "te12_PreTest","te12_PostTest","te12_Retention","te12_Delta",
                  "te21_PreTest","te21_PostTest","te21_Retention","te21_Delta")

if(save_files){
  filename_delta = paste(filename,'_delta_ave','.csv',sep='_')
  write.csv(x_delta, file=filename_delta, row.names=FALSE)
}


x_delta_wider <- x_delta %>%
  select(pp, task_label, Training, TrainingPhase, c_Delta, pitch_error_Delta, score_Delta,
         te12_Delta, te21_Delta) %>%
  pivot_wider(names_from = TrainingPhase, 
              values_from = c('c_Delta','pitch_error_Delta','score_Delta','te12_Delta','te21_Delta'))
names(x_delta_wider) <- c("pp","Task Label","Training","c Delta% PostTest","c Delta% Retention",
                        "Pitch Error Delta% PostTest","Pitch Error Delta% Retention",
                        "Score Delta% PostTest","Score Delta% Retention",
                        "TE S-U Delta% PostTest","TE S-U Delta% Retention",
                        "TE U-S Delta% PostTest","TE U-S Delta% Retention")

if(save_files){
  filename_delta = paste(filename,'_delta_ave_wide','.csv',sep='_')
  write.csv(x_delta_wider, file=filename_delta, row.names=FALSE)
}
