library(janitor)
library(tidyverse)
library(ggplot2)
library(ghibli)
library(lme4)
library(texreg)

setwd('~/logos/umcg_hidden_dynamics/analysis/performance')
# setwd('C:\\Users\\ddotov\\Downloads')
filename_scores_in <- 'Scores_2024-05-09.csv'
x<-read.csv(filename_scores_in) %>%
  clean_names()
x$training_phase <- sub("_", "", x$training_phase)
x$training_phase <- sub("est_", "est", x$training_phase)
summary(x)
x<-x[x$trial_duration>50,]
summary(x)


# Create a 'training condition' label in the long table
xcond <- x[x$training_phase=='Training',] %>%
  select(pp, task_code) %>%
  group_by(pp) %>%
  summarise(training_condition = mean(task_code)) %>%
  ungroup()
x <- x %>% left_join(., xcond, by = "pp")
as.factor(x$training_condition)

x$Training <- 'Periodic Fixed'
x$Training[x$training_condition==25] <- 'Chaotic Interactive'
x$Training[x$training_condition==30] <- 'Chaotic Non-interactive'
x$Training <- as.factor(x$Training)

x$task_code <- x$task_code + x$visual/2
x$task_code <- as.factor(x$task_code)
x$task_label <- 'Periodic Fixed'
x$task_label[x$task_code==12] <- 'Periodic Fadeout (SCT)'
x$task_label[x$task_code==20] <- 'Chaotic'
x$task_label[x$task_code==20.5] <- 'Chaotic Visual'
x$task_label[x$task_code==10] <- 'Chaotic Interactive Training'
x$task_label[x$task_code==25] <- 'Chaotic Interactive Training'
x$task_label[x$task_code==30] <- 'Chaotic Non-Interactive Training'
x$task_label <- as.factor(x$task_label)

x$training_phase <- as.factor(x$training_phase)
x$training_phase <- relevel(x$training_phase, ref='PreTest')


# Pseudo-trial number vector
trial_counter = 1
x$trial = 0
x$trial[1] = trial_counter
for (n in seq(2,dim(x)[1])) {
  if ((x$pp[n]==x$pp[n-1]) && (x$task_code[n]==x$task_code[n-1]) ) {
    trial_counter = trial_counter + 1
  }
  else {
    trial_counter = 1
  }
  x$trial[n] <- trial_counter
}
x$pp <- factor(x$pp)

# One single weird trial.
x <- x[!((x$training_phase == 'Retention') & (x$trial==6)),]

# filename = paste(filename_scores_in,'cleaned',Sys.Date(),'.csv',sep='_')
# write.csv(x, file=filename, row.names=FALSE)

# Visualize
colors<-ghibli_palette("PonyoMedium",7,type=("continuous"))[c(3,5,6)]
colors[1]<-ghibli_palette("MarnieMedium2",7,type=("continuous"))[6]
colors[2]<-ghibli_palette("MononokeMedium",7,type=("continuous"))[5]
colors[3]<-ghibli_palette("YesterdayMedium",7,type=("continuous"))[6]

sink("performance_scores_stats.txt")
print(Sys.Date())
sink()
# Plot performance scores and stats ~ trial
for (dv in c('score','c','pitch_error','tau')) {
# for (dv in c('tau')) {
  for (task in c('Periodic Fixed','Periodic Fadeout (SCT)','Chaotic','Chaotic Visual')) {
    xs <- x[(x$task_label==task) & (x$training_phase!='Training'),]
    xs$dv <- xs[,dv]

    if (dv=='score') {dv_lab = 'Score'}
    if (dv=='c') {dv_lab = 'Sync'}
    if (dv=='pitch_error') {dv_lab = 'Pitch Error'}
    if (dv=='tau') {dv_lab = 'Anticipation'}
    
    m1=lmer(dv ~ 1 + (1|pp),data=xs,REML=0)
    m2=lmer(dv ~ 1 + trial + (1|pp),data=xs,REML=0)
    m3=lmer(dv ~ 1 + trial + training_phase + (1|pp),data=xs,REML=0)
    m4=lmer(dv ~ 1 + trial + training_phase + training_phase:Training + (1|pp),data=xs,REML=0)
    sink("performance_scores_stats.txt", append = T)
    cat("\n###########\n\n")
    print(paste('Test: ',task,'; DV: ',dv,sep=''))
    print(anova(m1,m2,m3,m4))
    print(summary(m1))
    print(summary(m2))
    print(summary(m3))
    print(summary(m4))
    print(screenreg(list(m1,m2,m3,m4)))
    sink()
  
    if (dv=='score') {xs$fitted <- getME(m4,'X') %*% fixef(m4)}
    if (dv=='c') {xs$fitted <- getME(m4,'X') %*% fixef(m4)}
    if (dv=='pitch_error') {xs$fitted <- getME(m4,'X') %*% fixef(m4)}
    if (dv=='tau') {xs$fitted <- getME(m4,'X') %*% fixef(m4)}

    if (TRUE){
      g<-ggplot(data=xs, aes(x=trial, y=dv, colour=Training)) +
        facet_grid(~training_phase, scales="free_x") +
        geom_jitter(size=1, alpha=.2, width=.1, height=0) +
        # geom_line(aes(x=trial, y=0), col='black', size=1.2, alpha=.7) +
        # stat_summary(aes(x=trial, y=dv, colour=Training), fun='mean', geom='line', linewidth=2.2, alpha=1) +
        # stat_summary(aes(x=trial, y=dv, colour=Training), geom="ribbon", fun.data=mean_cl_boot, alpha=.5) +
        geom_line(aes(x=trial, y=fitted, colour=Training), linewidth=1.2, alpha=.8) +
        theme_classic() +
        theme(panel.background = element_rect(fill = "#111111",
              colour = "#000000", linewidth = 1, linetype = "solid")) +
        theme(legend.position="top") +
              #legend.position="top", legend.title=element_blank()) +
        #scale_x_continuous(breaks=seq(0,limit_lags,1), limits=c(0,limit_lags)) +
        #scale_y_continuous(limits=c(-.6,.4)) +
        labs(y = dv_lab) +
        labs(x = "Trial") + 
        scale_colour_manual(values=colors) + 
        ggtitle(paste('Test task: ', task, sep=''))
      print(g)
    }
    
    if (TRUE){
      filename = paste("performance",dv,task,Sys.Date(),'.png',sep='_')
      filename <- sub(" ", "_", filename)
      ggsave(filename, width=8, height=4, units='in', dpi=600)
    }
  }
}
