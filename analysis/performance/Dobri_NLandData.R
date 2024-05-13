#----------------------------------------------
library(tidyverse)
#----------------------------------------------

X <- read.csv('Scores_2024-05-09.csv_cleaned_2024-05-12_.csv',sep=',')

# Important columns:
# pp (participant)
# training_phase
# score (DV)
# task_label
# Training


# Filter out the training conditions from test conditions
# test_data <- filter(X, training_phase != "Training")
test_data <- filter(X,
                    task_label %in% c("Periodic Fixed", "Periodic Fadeout (SCT)",
                                      "Chaotic", "Chaotic Visual"))

#----------------------------------------------
# Average the score for each task_label per pp and training_phase
# Delta scores as percentage (per pp)
#    PreTest vs PostTest; PreTest vs Retention
# Widen the table such that Pre/Post/Retention are in their own columns
# and calculate difference scores between target columns
#----------------------------------------------
fun_ave_pivot_diff <- function(test_data, y){
  test_data$dv <- y
  diff_data <- test_data %>%
    group_by(pp, training_phase, task_label) %>%
    summarise(mean_PerTask = mean(dv)) %>%
    ungroup() %>%
    pivot_wider(names_from = training_phase, values_from = mean_PerTask) %>%
    group_by(pp) %>%
    mutate(PostDelta = (PostTest - PreTest)/PreTest, RetentionDelta = (Retention - PreTest)/PreTest) %>% # Calc diff scores
    inner_join(select(X, pp, Training), relationship = "many-to-many") %>% # Add back columns of interest
    distinct(pp, task_label, .keep_all = TRUE) %>% # Get rid of duplicated rows
    pivot_longer(cols=c("PostDelta","RetentionDelta"), 
                 names_to = "TrainingPhase", values_to = "Delta")
    return (diff_data)
}


diff_data <- fun_ave_pivot_diff(test_data, test_data$score)
Delta <- diff_data
diff_data <- fun_ave_pivot_diff(test_data, test_data$c)
Delta <- Delta %>% left_join(., diff_data, by = c("pp","task_label","Training","TrainingPhase"))
diff_data <- fun_ave_pivot_diff(test_data, test_data$pitch_error)
Delta <- Delta %>% left_join(., diff_data, by = c("pp","task_label","Training","TrainingPhase"))

# Sanity check
# head(diff_data, 20)

names(Delta) <- c("pp","task_label",
                  "score_PostTest","score_PreTest","score_Retention",
                  "Training","TrainingPhase",
                  "score_Delta","c_PostTest","c_PreTest","c_Retention","c_Delta",
                  "pitch_error_PostTest","pitch_error_PreTest","pitch_error_Retention","pitch_error_Delta")

# sum_data <- test_data %>%
#   group_by(pp, training_phase, task_label) %>%
#   summarise(mean_PerTask = mean(dv)) %>%
#   ungroup()
# 
# pivot_data <- sum_data %>%
#   pivot_wider(names_from = training_phase, values_from = mean_PerTask)
# 
# diff_data <- pivot_data %>%
#   group_by(pp) %>%
#   mutate(Post_Pre = PostTest - PreTest, Ret_Pre = Retention - PreTest) %>% # Calc diff scores
#   inner_join(select(X, pp, Training), relationship = "many-to-many") %>% # Add back columns of interest
#   distinct(pp, task_label, .keep_all = TRUE) # Get rid of duplicated rows

#----------------------------------------------
# Stats
library(ggplot2)
library(ghibli)
#----------------------------------------------

colors<-ghibli_palette("PonyoMedium",7,type=("continuous"))[c(3,5,6)]
colors[1]<-ghibli_palette("MarnieMedium2",7,type=("continuous"))[6]
colors[2]<-ghibli_palette("MononokeMedium",7,type=("continuous"))[5]
colors[3]<-ghibli_palette("YesterdayMedium",7,type=("continuous"))[6]

Delta$TrainingPhase <- factor(Delta$TrainingPhase)
Delta$Training <- factor(Delta$Training)
Delta$task_label <- factor(Delta$task_label)
for (dv in c('score_Delta','c_Delta','pitch_error_Delta')) {
  Delta$dv <- as.numeric(unlist(as.data.frame(Delta[,dv])))
  if (dv=='score_Delta') {dv_lab = 'Score Δ%'}
  if (dv=='c_Delta') {dv_lab = 'Sync Δ%'}
  if (dv=='pitch_error_Delta') {dv_lab = 'Pitch Error Δ%'}
  g<-ggplot(data=Delta, aes(x=Training, y=dv, colour=Training)) +
    facet_grid(task_label~TrainingPhase, scales="free_x") +
    geom_jitter(size=1, alpha=1, width=.1, height=0) +
    geom_violin(linewidth=1, alpha=.5, fill='grey') +
    geom_hline(aes(yintercept=0),
               color="blue", linetype="dashed", linewidth=1) +
    theme_classic() +
    theme(panel.background = element_rect(fill = "#ffffff",
                                          colour = "#888888", linewidth = 1, linetype = "solid")) +
    theme(legend.position="none", legend.title=element_blank()) +
    labs(y = dv_lab) +
    coord_cartesian(ylim=c(-.5,1.)) +
    scale_x_discrete(guide = guide_axis(angle = 10)) +
    scale_colour_manual(values=colors)
  print(g)
  if (TRUE){
    filename = paste("performance",dv,Sys.Date(),'.png',sep='_')
    filename <- sub(" ", "_", filename)
    ggsave(filename, width=6, height=10, units='in', dpi=600)
  }
}
