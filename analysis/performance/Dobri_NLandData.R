#----------------------------------------------
library(tidyverse)
#----------------------------------------------

filename = 'Scores_2024-05-09.csv_cleaned_2024-05-12_.csv'
X <- read.csv(filename,sep=',')

# Important columns:
# pp (participant)
# training_phase
# score (DV)
# task_label
# Training


# Filter out the training conditions from test conditions
test_data <- filter(X, training_phase != "Training")
# test_data <- filter(X, task_label %in% c("Periodic Fixed", "Periodic Fadeout (SCT)",
#                                       "Chaotic", "Chaotic Visual"))

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
    mutate(PostDelta = 100*(PostTest - PreTest)/PreTest, RetentionDelta = 100*(Retention - PreTest)/PreTest) %>% # Calc diff scores
    inner_join(select(test_data, pp, Training), relationship = "many-to-many") %>% # Add back columns of interest
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

Delta_wider <- Delta %>%
  select(pp, task_label, Training, TrainingPhase, c_Delta, pitch_error_Delta, score_Delta) %>%
  pivot_wider(names_from = TrainingPhase, values_from = c('c_Delta','pitch_error_Delta','score_Delta'))
names(Delta_wider) <- c("pp","Task Label","Training","c Delta% PostTest","c Delta% Retention",
                        "Pitch Error Delta% PostTest","Pitch Error Delta% Retention",
                        "Score Delta% PostTest","Score Delta% Retention")
# filename_delta = paste(filename,'_delta_ave_wide','.csv',sep='_')
# write.csv(Delta_wider, file=filename_delta, row.names=FALSE)


test_ave <- test_data %>%
  group_by(pp,task_label,training_phase) %>%
  summarise(score = mean(score),
            c = mean(c),
            pitch_error = mean(pitch_error)) %>%
  ungroup() %>%
  # Add back columns of interest
  inner_join(select(test_data, pp, task_label, Training), relationship = "many-to-many") %>%
  # Get rid of duplicated rows
  distinct(pp, task_label, training_phase, .keep_all = TRUE)

test_ave_wide <- test_ave %>%
  pivot_wider(names_from = training_phase, values_from = c('score','c','pitch_error'))
# filename_ave = paste(filename,'_ave_wide','.csv',sep='_')
# write.csv(test_ave_wide, file=filename_ave, row.names=FALSE)


#----------------------------------------------
# Figures
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


test_ave$training_phase <- factor(test_ave$training_phase)
test_ave$training_phase <- relevel(test_ave$training_phase, ref='PreTest')
test_ave$Training <- factor(test_ave$Training)
test_ave$Training <- relevel(test_ave$Training, ref='Periodic Fixed')
test_ave$pp <- factor(test_ave$pp)
for (dv in c('score','c','pitch_error')) {
  test_ave$dv <- as.numeric(unlist(as.data.frame(test_ave[,dv])))
  if (dv=='score') {dv_lab = 'Score'}
  if (dv=='c') {dv_lab = 'Sync'}
  if (dv=='pitch_error') {dv_lab = 'Pitch Error'}
  g<-ggplot(data=test_ave, aes(x=training_phase, y=dv)) +
    facet_grid(task_label~Training, scales="free_y") +
    geom_violin(linewidth=1, alpha=.5, fill='grey') +
    geom_line(aes(x=training_phase, y=dv, group=pp), linewidth=.5, alpha=.4, colour = 'black') +
    geom_jitter(size=2, alpha=1, width=0, height=0, colour = 'black') +
    # linetype='solid',
    # stat_summary(geom = "line", fun.y = mean, size = 3, colour = 'red') +
    theme_classic() +
    theme(panel.background = element_rect(fill = "#ffffff",
                                          colour = "#888888", linewidth = 1, linetype = "solid")) +
    theme(legend.position="top", legend.title=element_blank()) +
    labs(y=dv_lab, x='') +
    # coord_cartesian(ylim=c(-.5,1.)) +
    scale_colour_manual(values='black')
  print(g)
  if (FALSE){
    filename = paste("performance_lines",dv,Sys.Date(),'.png',sep='_')
    filename <- sub(" ", "_", filename)
    ggsave(filename, width=7, height=6, units='in', dpi=600)
  }
}


#----------------------------------------------
# Stats
#----------------------------------------------

# Check assumptions: normalcy; variance is tested in Sphericity
shapiro.test(Delta$score_Delta)

ggplot(data = Delta, aes(sample = score_Delta)) + stat_qq() +
  stat_qq_line()



# ANOVA
Delta_aov <- aov(score_Delta ~ task_label * Training, data = Delta)
summary(Delta_aov)

Delta_lm <- lm(score_Delta ~ task_label * Training, data = Delta)
anova(Delta_lm)
summary(Delat_lm) # this one has the intercept
