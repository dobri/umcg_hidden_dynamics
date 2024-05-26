library(tidyverse)
library(ggplot2)
library(ghibli)
library(lme4)
library(lmerTest)
library(texreg)
library(rstatix)


#----------------------------------------------
# Part A. Scores ~ Trial
# Step 1. Plot performance scores ~ trial
#----------------------------------------------

rm(list=ls())
# setwd('~/logos/umcg_hidden_dynamics/analysis/performance')
# setwd('C:\\Users\\ddotov\\Downloads')
filename_scores_in <- 'Scores_2024-05-09.csv_cleaned_.csv'
x <- read.csv(filename_scores_in,sep=',')
x$training_phase <- factor(x$training_phase)
x$training_phase <- relevel(x$training_phase, ref='PreTest')
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

# for (dv in c('score','c','pitch_error','tau')) {
for (dv in c('score_delta','c_delta','pitch_error_delta','tau_delta')) {
  for (task in c('Periodic Fixed','Periodic Fadeout (SCT)','Chaotic','Chaotic Visual')) {
    xs <- x[(x$task_label==task) & (x$training_phase!='Training'),]
    xs$dv <- xs[,dv]
    
    if (dv=='score') {dv_lab = 'Score'}
    if (dv=='c') {dv_lab = 'Sync'}
    if (dv=='pitch_error') {dv_lab = 'Pitch Error'}
    if (dv=='tau') {dv_lab = 'Anticipation'}
    if (dv=='score_delta') {dv_lab = 'Score Δ%'}
    if (dv=='c_delta') {dv_lab = 'Sync Δ%'}
    if (dv=='pitch_error_delta') {dv_lab = 'Pitch Error Δ%'}
    if (dv=='tau_delta') {dv_lab = 'Anticipation Δ%'}
    
    g<-ggplot(data=xs, aes(x=trial, y=dv, colour=Training)) +
      facet_grid(~training_phase, scales="free_x") +
      geom_jitter(size=1, alpha=.2, width=.1, height=0) +
      geom_line(aes(x=trial, y=0), col='black', linewidth=1.2, alpha=.7) +
      stat_summary(aes(x=trial, y=dv, colour=Training), fun='mean', geom='line', linewidth=2.2, alpha=1) +
      stat_summary(aes(x=trial, y=dv, colour=Training), geom="ribbon", fun.data=mean_cl_boot, alpha=.5) +
      # geom_line(aes(x=trial, y=fitted, colour=Training), linewidth=1.2, alpha=.8) +
      theme_classic() +
      theme(panel.background = element_rect(fill = "#111111",
                                            colour = "#000000", linewidth = 1, linetype = "solid")) +
      theme(legend.position="top") +
      labs(y = dv_lab) +
      labs(x = "Trial") + 
      scale_colour_manual(values=colors) + 
      ggtitle(paste('Test task: ', task, sep=''))
    print(g)

    if (FALSE){
      filename = paste("performance_trials_all",dv,task,Sys.Date(),'.png',sep='_')
      filename <- sub(" ", "_", filename)
      ggsave(filename, width=8, height=4, units='in', dpi=600)
    }
  }
}


#----------------------------------------------
# Step 2. Stats delta scores ~ trial
#----------------------------------------------
sink("performance_scores_stats.txt")
print(Sys.Date())
sink()
# Plot performance scores and stats ~ trial
# for (dv in c('score','c','pitch_error','tau')) {
for (dv in c('score_delta','c_delta','pitch_error_delta','tau_delta')) {
  for (task in c('Periodic Fixed','Periodic Fadeout (SCT)','Chaotic','Chaotic Visual')) {
    xs <- x[(x$task_label==task) & (x$training_phase!='Training'),]
    xs$dv <- xs[,dv]
    
    if (dv=='score') {dv_lab = 'Score'}
    if (dv=='c') {dv_lab = 'Sync'}
    if (dv=='pitch_error') {dv_lab = 'Pitch Error'}
    if (dv=='tau') {dv_lab = 'Anticipation'}
    if (dv=='score_delta') {dv_lab = 'Score Δ%'}
    if (dv=='c_delta') {dv_lab = 'Sync Δ%'}
    if (dv=='pitch_error_delta') {dv_lab = 'Pitch Error Δ%'}
    if (dv=='tau_delta') {dv_lab = 'Anticipation Δ%'}
    
    m1=lmer(dv ~ 1 + (1|pp),data=xs,REML=0)
    m2=lmer(dv ~ 1 + trial + (1|pp),data=xs,REML=0)
    m3=lmer(dv ~ 1 + trial + training_phase + (1|pp),data=xs,REML=0)
    m4=lmer(dv ~ 1 + trial + training_phase + training_phase:Training + (1|pp),data=xs,REML=0)
    sink("performance_scores_stats.txt", append = T)
    cat("\n###########\n\n")
    cat(paste('Test: ',task,'; DV: ',dv,sep=''))
    print(anova(m1,m2,m3,m4))
    print(summary(m1))
    print(summary(m2))
    print(summary(m3))
    print(summary(m4))
    print(screenreg(list(m1,m2,m3,m4)))
    sink()
  }
}


#----------------------------------------------
# Part B. Plot test data only, averaged per task. Line plots
#----------------------------------------------
rm(list=ls())
filename = 'Scores_2024-05-09.csv_cleaned_.csv__ave_.csv'
x_test_ave <- read.csv(filename,sep=',')

x_test_ave$training_phase <- factor(x_test_ave$training_phase)
x_test_ave$training_phase <- relevel(x_test_ave$training_phase, ref='PreTest')
x_test_ave$Training <- factor(x_test_ave$Training)
x_test_ave$Training <- relevel(x_test_ave$Training, ref='Periodic Fixed')
x_test_ave$pp <- factor(x_test_ave$pp)
for (dv in c('score','c','pitch_error')) {
  x_test_ave$dv <- as.numeric(unlist(as.data.frame(x_test_ave[,dv])))
  if (dv=='score') {dv_lab = 'Score'}
  if (dv=='c') {dv_lab = 'Sync'}
  if (dv=='pitch_error') {dv_lab = 'Pitch Error'}
  g<-ggplot(data=x_test_ave, aes(x=training_phase, y=dv)) +
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
    filename = paste("performance_lines_ave",dv,Sys.Date(),'.png',sep='_')
    filename <- sub(" ", "_", filename)
    ggsave(filename, width=7, height=6, units='in', dpi=600)
  }
}


#----------------------------------------------
# Part C. Plot test data only, deltas averaged per task. Violins.
#----------------------------------------------

rm(list=ls())
filename = 'Scores_2024-05-09.csv_cleaned_.csv__delta_ave_.csv'
x_delta <- read.csv(filename,sep=',')
x_delta$TrainingPhase <- factor(x_delta$TrainingPhase)
x_delta$Training <- factor(x_delta$Training)
x_delta$task_label <- factor(x_delta$task_label)


colors<-ghibli_palette("PonyoMedium",7,type=("continuous"))[c(3,5,6)]
colors[1]<-ghibli_palette("MarnieMedium2",7,type=("continuous"))[6]
colors[2]<-ghibli_palette("MononokeMedium",7,type=("continuous"))[5]
colors[3]<-ghibli_palette("YesterdayMedium",7,type=("continuous"))[6]
for (dv in c('score_Delta','c_Delta','pitch_error_Delta')) {
  x_delta$dv <- as.numeric(unlist(as.data.frame(x_delta[,dv])))
  if (dv=='score_Delta') {dv_lab = 'Score Δ%'}
  if (dv=='c_Delta') {dv_lab = 'Sync Δ%'}
  if (dv=='pitch_error_Delta') {dv_lab = 'Pitch Error Δ%'}
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
    coord_cartesian(ylim=c(-50,100.)) +
    scale_x_discrete(guide = guide_axis(angle = 10)) +
    scale_colour_manual(values=colors)
  print(g)
  if (FALSE){
    filename = paste("performance_ave",dv,Sys.Date(),'.png',sep='_')
    filename <- sub(" ", "_", filename)
    ggsave(filename, width=6, height=10, units='in', dpi=600)
  }
}


#----------------------------------------------
# Part D. Stats on test data.
#----------------------------------------------

sink("performance_deltas_stats.txt")
sink()
for (dv in c('score_Delta','c_Delta','pitch_error_Delta')) {
  x_delta$dv <- as.numeric(unlist(as.data.frame(x_delta[,dv])))
  if (dv=='score_Delta') {dv_lab = 'Score Δ%'}
  if (dv=='c_Delta') {dv_lab = 'Sync Δ%'}
  if (dv=='pitch_error_Delta') {dv_lab = 'Pitch Error Δ%'}
  sink("performance_deltas_stats.txt", append = T)
  for (task in c('Periodic Fixed','Periodic Fadeout (SCT)','Chaotic','Chaotic Visual')) {
    cat("\n###########\n\n")
    cat(paste('DV: ',dv_lab,sep=''))
    cat(paste('Test: ',task,sep=''))
    xs <- x_delta[x_delta$task_label==task,]
    m1=lmer(dv ~ 1 + (1|pp),data=xs,REML=0)
    m2=lmer(dv ~ 1 + TrainingPhase + (1|pp),data=xs,REML=0)
    m3=lmer(dv ~ 1 + TrainingPhase + Training + (1|pp),data=xs,REML=0)
    m4=lmer(dv ~ 1 + TrainingPhase * Training + (1|pp),data=xs,REML=0)
    print(anova(m1,m2,m3,m4))
    print(summary(m1))
    print(summary(m2))
    print(summary(m3))
    print(summary(m4))
    print(screenreg(list(m1,m2,m3,m4)))
  }
  sink()
}


sink("performance_deltas_stats_ttests.txt")
sink()
for (dv in c('score_Delta','c_Delta','pitch_error_Delta')) {
  x_delta$dv <- as.numeric(unlist(as.data.frame(x_delta[,dv])))
  if (dv=='score_Delta') {dv_lab = 'Score Δ%'}
  if (dv=='c_Delta') {dv_lab = 'Sync Δ%'}
  if (dv=='pitch_error_Delta') {dv_lab = 'Pitch Error Δ%'}
  sink("performance_deltas_stats_ttests.txt", append = T)
  for (task in c('Periodic Fixed','Periodic Fadeout (SCT)','Chaotic','Chaotic Visual')) {
    cat("\n########### t-tests ~ 0 \n\n")
    cat(paste('DV: ',dv_lab,sep=''))
    cat(paste('Test: ',task,sep=''))

    xs <- x_delta[x_delta$task_label==task,]
    stat.test <- xs %>%
      group_by(Training,TrainingPhase) %>%
      t_test(dv ~ 1, mu = 0) %>%
      adjust_pvalue(method = "fdr") %>%
      add_significance()
    print(stat.test)
  }
  sink()
}
