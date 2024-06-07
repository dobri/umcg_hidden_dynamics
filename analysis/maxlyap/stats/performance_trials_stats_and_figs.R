library(tidyverse)
library(ggplot2)
library(ghibli)
library(lme4)
library(lmerTest)
library(texreg)
library(rstatix)


#----------------------------------------------
# Part 0. Training scores ~ Trial
# Step 1. Plot
#----------------------------------------------

rm(list=ls())
# setwd('~/logos/umcg_hidden_dynamics/analysis/maxlyap/stats')
# setwd('C:\\Users\\ddotov\\Downloads')
filename_scores_in <- 'Scores_MaxLyap_2024-06-07.csv_cleaned_.csv'
x <- read.csv(filename_scores_in,sep=',')
x$training_phase_label <- factor(x$training_phase_label)
x$training_phase_label <- relevel(x$training_phase_label, ref='PreTest')
x$training_condition <- factor(x$training_condition)
x$Training <- factor(x$Training)
x$task_label <- factor(x$task_label)
x$task_code <- factor(x$task_code)
x$pp <- factor(x$pp)


# Visualize
colors<-ghibli_palette("PonyoMedium",7,type=("continuous"))[c(3,5,6)]
colors[1]<-ghibli_palette("MarnieMedium2",7,type=("continuous"))[6]
colors[2]<-ghibli_palette("MononokeMedium",7,type=("continuous"))[5]
colors[3]<-ghibli_palette("YesterdayMedium",7,type=("continuous"))[6]

for (dv in c('max_lyap_stim_short','max_lyap_stim_long','max_lyap_user_short','max_lyap_user_long')) {
  xs <- x[(x$training_phase_label=='Training'),]
  xs$dv <- xs[,dv]
  
  if (dv=='max_lyap_stim_short') {dv_lab = 'Max Lyapunov Short, Stimulus'}
  if (dv=='max_lyap_stim_long') {dv_lab = 'Max Lyapunov Long, Stimulus'}
  if (dv=='max_lyap_user_short') {dv_lab = 'Max Lyapunov Short, User'}
  if (dv=='max_lyap_user_long') {dv_lab = 'Max Lyapunov Long, User'}
  
  g<-ggplot(data=xs, aes(x=trial, y=dv, colour=Training)) +
    facet_grid(~Training) +
    geom_jitter(size=1, alpha=.2, width=.1, height=0) +
    geom_line(aes(x=trial, y=0), col='white', linewidth=1.2, alpha=.7) +
    stat_summary(aes(x=trial, y=dv, colour=Training), geom="ribbon", fun.data=mean_cl_boot, alpha=.7) +
    stat_summary(aes(x=trial, y=dv, colour=Training), fun='mean', geom='line', linewidth=2.2, alpha=1) +
    theme_classic() +
    theme(panel.background = element_rect(fill = "#111111",
                                          colour = "#000000", linewidth = 1, linetype = "solid")) +
    theme(legend.position="none") +
    labs(y = dv_lab) +
    labs(x = "Trial") + 
    #coord_cartesian(ylim=c(0,2)) +
    scale_colour_manual(values=colors)
  print(g)
  
  if (FALSE){
    filename = paste("ml_training_trials_all",dv,Sys.Date(),'.png',sep='_')
    filename <- sub(" ", "_", filename)
    ggsave(filename, width=10, height=4, units='in', dpi=600)
  }
}


#----------------------------------------------
# Part A. Test scores ~ Trial
# Step 1. Plot
#----------------------------------------------


#----------------------------------------------
# Step 2. Stats training TE ~ trial
#----------------------------------------------
sink("ml_training_stats.txt")
print(Sys.Date())
sink()
for (dv in c('max_lyap_stim_short','max_lyap_stim_long','max_lyap_user_short','max_lyap_user_long')) {
  for (group in c('Periodic Fixed','Chaotic Interactive','Chaotic Non-interactive')) {
    xs <- x[(x$Training==group) & (x$training_phase_label=='Training'),]
    xs$dv <- xs[,dv]
    
    if (dv=='max_lyap_stim_short') {dv_lab = 'Max Lyapunov Short, Stimulus'}
    if (dv=='max_lyap_stim_long') {dv_lab = 'Max Lyapunov Long, Stimulus'}
    if (dv=='max_lyap_user_short') {dv_lab = 'Max Lyapunov Short, User'}
    if (dv=='max_lyap_user_long') {dv_lab = 'Max Lyapunov Long, User'}

    m1=lmer(dv ~ 1 + (1|pp),data=xs,REML=0)
    m2=lmer(dv ~ 1 + trial + (1|pp),data=xs,REML=0)
    #m3=lmer(dv ~ 1 + trial + Training + (1|pp),data=xs,REML=0)
    #m4=lmer(dv ~ 1 + trial * Training + (1|pp),data=xs,REML=0)
    sink("ml_training_stats.txt", append = T)
    cat("\n###########\n\n")
    cat(paste('Training'))
    cat("\n")
    cat('DV: ',dv,sep='')
    cat("\n")
    cat('Group: ',group,sep='')
    cat("\n")
    # print(anova(m1,m2,m3,m4))
    print(anova(m1,m2))
    print(summary(m1))
    print(summary(m2))
    # print(summary(m3))
    # print(summary(m4))
    # print(screenreg(list(m1,m2,m3,m4)))
    print(screenreg(list(m1,m2)))
    sink()
  }
}


#----------------------------------------------
# Part B. Plot test data, averaged per task. Line plots
#----------------------------------------------
rm(list=ls())
filename = 'Scores_MaxLyap_2024-06-07.csv_cleaned_.csv__ave_.csv'
x_test_ave <- read.csv(filename,sep=',')

x_test_ave$training_phase_label <- factor(x_test_ave$training_phase_label)
x_test_ave$training_phase_label <- relevel(x_test_ave$training_phase_label, ref='PreTest')
x_test_ave$Training <- factor(x_test_ave$Training)
x_test_ave$Training <- relevel(x_test_ave$Training, ref='Periodic Fixed')
x_test_ave$pp <- factor(x_test_ave$pp)
for (dv in c('max_lyap_stim_short','max_lyap_stim_long','max_lyap_user_short','max_lyap_user_long')) {
  x_test_ave$dv <- as.numeric(unlist(as.data.frame(x_test_ave[,dv])))
  if (dv=='max_lyap_stim_short') {dv_lab = 'Max Lyapunov Short, Stimulus'}
  if (dv=='max_lyap_stim_long') {dv_lab = 'Max Lyapunov Long, Stimulus'}
  if (dv=='max_lyap_user_short') {dv_lab = 'Max Lyapunov Short, User'}
  if (dv=='max_lyap_user_long') {dv_lab = 'Max Lyapunov Long, User'}
  
  g<-ggplot(data=x_test_ave, aes(x=training_phase_label, y=dv)) +
    facet_grid(task_label~Training, scales="free_y") +
    # geom_line(aes(x=training_phase_label, y=dv*0, group=pp), linetype = 'dashed', linewidth=.5, alpha=.5, colour = 'black') +
    geom_violin(linewidth=1, alpha=.5, fill='grey') +
    geom_line(aes(x=training_phase_label, y=dv, group=pp), linewidth=.5, alpha=.4, colour = 'black') +
    geom_jitter(size=2, alpha=1, width=0, height=0, colour = 'black') +
    theme_classic() +
    theme(panel.background = element_rect(fill = "#ffffff",
                                          colour = "#888888", linewidth = 1, linetype = "solid")) +
    theme(legend.position="top", legend.title=element_blank()) +
    labs(y=dv_lab, x='') +
    scale_colour_manual(values='black')
  print(g)
  if (TRUE){
    filename = paste("ml_lines_ave",dv,Sys.Date(),'.png',sep='_')
    filename <- sub(" ", "_", filename)
    ggsave(filename, width=7, height=7, units='in', dpi=600)
  }
}


#----------------------------------------------
# Part C. Plot test data, deltas averaged per task. Violins.
#----------------------------------------------

rm(list=ls())
filename = 'Scores_MaxLyap_2024-06-07.csv_cleaned_.csv__delta_ave_.csv'
x_delta <- read.csv(filename,sep=',')
x_delta$TrainingPhase <- factor(x_delta$TrainingPhase)
x_delta$Training <- factor(x_delta$Training)
x_delta$task_label <- factor(x_delta$task_label)

colors<-ghibli_palette("PonyoMedium",7,type=("continuous"))[c(3,5,6)]
colors[1]<-ghibli_palette("MarnieMedium2",7,type=("continuous"))[6]
colors[2]<-ghibli_palette("MononokeMedium",7,type=("continuous"))[5]
colors[3]<-ghibli_palette("YesterdayMedium",7,type=("continuous"))[6]
for (dv in c('mlss_Delta','mlsl_Delta','mlus_Delta','mlul_Delta')) {
  x_delta$dv <- as.numeric(unlist(as.data.frame(x_delta[,dv])))
  if (dv=='mlss_Delta') {dv_lab = 'Maximum Lyapunov Short, Stimulus Δ%'}
  if (dv=='mlsl_Delta') {dv_lab = 'Maximum Lyapunov Long, Stimulus Δ%'}
  if (dv=='mlus_Delta') {dv_lab = 'Maximum Lyapunov Short, User Δ%'}
  if (dv=='mlul_Delta') {dv_lab = 'Maximum Lyapunov Long, User Δ%'}
  g<-ggplot(data=x_delta, aes(x=Training, y=dv, colour=Training)) +
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
    # coord_cartesian(ylim=c(-50,100.)) +
    scale_x_discrete(guide = guide_axis(angle = 10)) +
    scale_colour_manual(values=colors)
  print(g)
  if (TRUE){
    filename = paste("ml_tests_ave",dv,Sys.Date(),'.png',sep='_')
    filename <- sub(" ", "_", filename)
    ggsave(filename, width=6, height=10, units='in', dpi=600)
  }
}


#----------------------------------------------
# Part D. Stats on test data.
#----------------------------------------------

rm(list=ls())
filename = 'Scores_MaxLyap_2024-06-07.csv_cleaned_.csv__ave_.csv'
x_test_ave <- read.csv(filename,sep=',')

x_test_ave$training_phase_label <- factor(x_test_ave$training_phase_label)
x_test_ave$training_phase_label <- relevel(x_test_ave$training_phase_label, ref='PreTest')
x_test_ave$Training <- factor(x_test_ave$Training)
x_test_ave$Training <- relevel(x_test_ave$Training, ref='Periodic Fixed')
x_test_ave$pp <- factor(x_test_ave$pp)
sink("ml_tests_stats.txt")
sink()
for (dv in c('max_lyap_stim_short','max_lyap_stim_long','max_lyap_user_short','max_lyap_user_long')) {
  x_test_ave$dv <- as.numeric(unlist(as.data.frame(x_test_ave[,dv])))
  if (dv=='max_lyap_stim_short') {dv_lab = 'Max Lyapunov Short, Stimulus'}
  if (dv=='max_lyap_stim_long') {dv_lab = 'Max Lyapunov Long, Stimulus'}
  if (dv=='max_lyap_user_short') {dv_lab = 'Max Lyapunov Short, User'}
  if (dv=='max_lyap_user_long') {dv_lab = 'Max Lyapunov Long, User'}
  sink("ml_tests_stats.txt", append = T)
  for (task in c('Periodic Fixed','Periodic Fadeout (SCT)','Chaotic','Chaotic Visual')) {
    cat("\n###########\n\n")
    cat(paste('DV: ',dv_lab,sep=''))
    cat("\n")
    cat(paste('Test: ',task,sep=''))
    cat("\n")
    xs <- x_test_ave[x_test_ave$task_label==task,]
    m1=lmer(dv ~ 1 + (1|pp),data=xs,REML=0)
    m2=lmer(dv ~ 1 + training_phase_label + (1|pp),data=xs,REML=0)
    m3=lmer(dv ~ 1 + training_phase_label + Training + (1|pp),data=xs,REML=0)
    m4=lmer(dv ~ 1 + training_phase_label * Training + (1|pp),data=xs,REML=0)
    print(anova(m1,m2,m3,m4))
    print(summary(m1))
    print(summary(m2))
    print(summary(m3))
    print(summary(m4))
    print(screenreg(list(m1,m2,m3,m4)))
  }
  sink()
}
