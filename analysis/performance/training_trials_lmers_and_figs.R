library(janitor)
library(tidyverse)
library(ggplot2)
library(ghibli)


setwd('~/logos/umcg_hidden_dynamics/analysis/performance')
# setwd('C:\\Users\\ddotov\\Downloads')
x<-read.csv('Scores_2024-05-09.csv') %>%
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

# Visualize
colors<-ghibli_palette("PonyoMedium",7,type=("continuous"))[c(3,5,6)]
colors[1]<-ghibli_palette("MarnieMedium2",7,type=("continuous"))[6]
colors[2]<-ghibli_palette("MononokeMedium",7,type=("continuous"))[5]
colors[3]<-ghibli_palette("YesterdayMedium",7,type=("continuous"))[6]

# Plot performance scores and stats ~ trial
for (dv in c('score','c','pitch_error','tau')) {
# for (dv in c('tau')) {
  for (task in c('Periodic Fixed','Periodic Fadeout (SCT)','Chaotic','Chaotic Visual')) {
    print(paste('Test:',task,', DV:',dv))

    xs <- x[(x$task_label==task) & (x$training_phase!='Training'),]
    xs$dv <- xs[,dv]

    if (dv=='score') {dv_lab = 'Score'}
    if (dv=='c') {dv_lab = 'Sync'}
    if (dv=='pitch_error') {dv_lab = 'Pitch Error'}
    if (dv=='tau') {dv_lab = 'Anticipation'}
    
    g<-ggplot(data=xs, aes(x=trial, y=dv, colour=Training)) +
      facet_grid(~training_phase, scales="free_x") +
      geom_jitter(size=1, alpha=.2, width=.1, height=0) +
      #geom_line(aes(x=trial, y=0), col='black', size=1.2, alpha=.7) +
      stat_summary(aes(x=trial, y=dv, colour=Training), fun='mean', geom='line', linewidth=2.2, alpha=1) +
      stat_summary(aes(x=trial, y=dv, colour=Training), geom="ribbon", fun.data=mean_cl_boot, alpha=.5) +
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
    
    if (TRUE){
      filename = paste("performance",dv,task,Sys.Date(),'.png',sep='_')
      filename <- sub(" ", "_", filename)
      ggsave(filename, width=8, height=4, units='in', dpi=600)
    }
  }
}


### END
g<-list('vector',4)
counter = 0
for (dv in c(2,3,5,4)){
  
  x$dv<-x[,dv]
  
  if (names(x)[dv]=='score') {dv_lab = 'Sync & Match Score [C/RMSE]'}
  if (names(x)[dv]=='cmax') {dv_lab = 'C'}
  if (names(x)[dv]=='tau') {dv_lab = 'τ'}
  if (names(x)[dv]=='rmse') {dv_lab = 'RMSE'}
  
  m1=lmer(dv ~ 1 + (-1+trial|pp),data=x,REML=0)
  m2=lmer(dv ~ 1 + trial + (-1+trial|pp),data=x,REML=0)
  m3=lmer(dv ~ 1 + trial+task + (-1+trial|pp),data=x,REML=0)
  m4=lmer(dv ~ 1 + trial*task + (-1+trial|pp),data=x,REML=0)

  # if (dv==2) {x$fitted <- getME(m4,'X') %*% fixef(m4)} 
  # if (dv==3) {x$fitted <- getME(m4,'X') %*% fixef(m4)}
  # if (dv==4) {x$fitted <- getME(m3,'X') %*% fixef(m3)}
  # if (dv==5) {x$fitted <- getME(m2,'X') %*% fixef(m2)}
  if (dv==2) {x$fitted <- getME(m3,'X') %*% fixef(m3)}
  if (dv==3) {x$fitted <- getME(m3,'X') %*% fixef(m3)}
  if (dv==4) {x$fitted <- getME(m3,'X') %*% fixef(m3)}
  if (dv==5) {x$fitted <- getME(m3,'X') %*% fixef(m3)}
  
  counter = counter + 1
  g[[counter]] <- ggplot(data=x, aes(x=trial, y=dv, colour=as.factor(task))) +
    geom_jitter(size=1, alpha=.5, width=.5, height=.1) +
    geom_line(aes(x=trial, y=fitted, colour=task), size=1.2, alpha=.7) +
    #geom_line(aes(x=trial, y=0), col='black', size=1.2, alpha=.7) +
    theme_classic() +
    theme(panel.background = element_rect(fill = "#111111",
                                          colour = "#000000",size = 1, linetype = "solid")) +
    theme(legend.position="top",legend.title=element_blank()) +
    labs(y = dv_lab) +
    labs(x = "Trial")
  g[[counter]] <- g[[counter]] + scale_colour_manual(values=colors)
}
multiplot(plotlist=g,cols=4)

if (FALSE){
  filename=paste("perf_scores_training_with_lin_model",'_',Sys.Date(),'.png',sep='')
  png(filename=filename,width=16,height=6,units="in",res=300)
  multiplot(plotlist=g,cols=4)
  dev.off()
}



