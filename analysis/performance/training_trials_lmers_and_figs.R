
setwd('~/logos/c3/umcg_hidden_dynamics/extract_measures_notes/')
x<-read.csv('all_scores.csv')

summary(x)

trial_counter = 0
ind1 <- regexpr('_pp',as.character(x$raw.data.file[1]))
ind2 <- regexpr('_d',as.character(x$raw.data.file[1]))
pp0 <- as.numeric(substr(as.character(x$raw.data.file[1]),ind1[1]+3,ind2[1]-1))
for (n in seq(1,dim(x)[1])) {
  ind1 <- regexpr('_pp',as.character(x$raw.data.file[n]))
  ind2 <- regexpr('_d',as.character(x$raw.data.file[n]))
  pp <- as.numeric(substr(as.character(x$raw.data.file[n]),ind1[1]+3,ind2[1]-1))
  if (pp!=pp0) {trial_counter = 1; pp0=pp} else {trial_counter = trial_counter + 1}
  x$trial[n] <- trial_counter
  x$pp[n] <- pp
  
  ind<-regexpr('task',as.character(x$raw.data.file[n]))
  x$task_num[n] <- as.numeric(substr(as.character(x$raw.data.file[n]),ind[1]+4,ind[1]+6))
  if (x$task_num[n]==25) x$task[n]='coupled_unstable'
  if (x$task_num[n]==30) x$task[n]='uncoupl_unstable'
  if (x$task_num[n]==10) x$task[n]='uncoupl_periodic'
  
  x$trial_id[n] <- as.numeric(substr(as.character(x$raw.data.file[n]),ind[1]-7,ind[1]-2))
  
  ind<-regexpr('aud',as.character(x$raw.data.file[n]))
  x$aud[n] <- as.numeric(substr(as.character(x$raw.data.file[n]),ind[1]+3,ind[1]+3))
  ind<-regexpr('vis',as.character(x$raw.data.file[n]))
  x$vis[n] <- as.numeric(substr(as.character(x$raw.data.file[n]),ind[1]+3,ind[1]+3))
  ind<-regexpr('eps',as.character(x$raw.data.file[n]))
  x$eps[n] <- as.numeric(substr(as.character(x$raw.data.file[n]),ind[1]+3,ind[1]+5))
  ind<-regexpr('_d',as.character(x$raw.data.file[n]))
  x$dn[n] <- as.numeric(substr(as.character(x$raw.data.file[n]),ind[1]+2,ind[1]+2))
  
}
x$pp <- factor(x$pp)


# Linear mixed-effects models
library(lme4)
library(texreg)

for (dv in c(2,3,4,5)){
  x$dv<-x[,dv]
  
  if (names(x)[dv]=='score') {dv_lab = 'Sync & Match Score [C/RMSE]'}
  if (names(x)[dv]=='cmax') {dv_lab = 'C'}
  if (names(x)[dv]=='tau') {dv_lab = 'τ'}
  if (names(x)[dv]=='rmse') {dv_lab = 'RMSE'}
  
  m1=lmer(dv ~ 1 + (-1+trial|pp),data=x,REML=0)
  m2=lmer(dv ~ 1 + trial + (-1+trial|pp),data=x,REML=0)
  m3=lmer(dv ~ 1 + trial+task + (-1+trial|pp),data=x,REML=0)
  m4=lmer(dv ~ 1 + trial*task + (-1+trial|pp),data=x,REML=0)
  # These don't seem to work.
  m5a=lmer(dv ~ 1 + trial*task + (1 + trial|pp),data=x,REML=0)
  m5b=lmer(dv ~ 1 + trial*task + (1|pp) + (trial|pp),data=x,REML=0)
  
  sink(paste("diary_lmems_trial_task_",'_',Sys.Date(),sep=''),append=TRUE)
  print(dv_lab)
  print(dv_lab)
  print(dv_lab)
  print(anova(m1,m2,m3,m4))
  print(summary(m1))
  print(summary(m2))
  print(summary(m3))
  print(summary(m4))
  print(screenreg(list(m1,m2,m3,m4)))
  sink()
  htmlreg(list(m1,m2,m3,m4), file = paste("texreg_lmems_",names(x)[dv],'_',Sys.Date(),'.doc',sep=''), single.row = FALSE, digits=3, inline.css=FALSE, doctype=TRUE, html.tag=TRUE, head.tag=TRUE, body.tag=TRUE)
}


# Plot performance scores and stats ~ trial
source('~/logos/c3/umcg_hidden_dynamics/extract_measures_notes/multiplot.R')
library(ggplot2)
library(ghibli)

colors<-ghibli_palette("PonyoMedium",7,type=("continuous"))[c(3,5,6)]
# colors[1]<-ghibli_palette("MarnieDark2",7,type=("continuous"))[6]
# colors[2]<-ghibli_palette("MononokeDark",7,type=("continuous"))[5]
# colors[3]<-ghibli_palette("YesterdayDark",7,type=("continuous"))[6]
# Autumn color palette!
colors[1]<-ghibli_palette("MarnieMedium2",7,type=("continuous"))[6]
colors[2]<-ghibli_palette("MononokeMedium",7,type=("continuous"))[5]
colors[3]<-ghibli_palette("YesterdayMedium",7,type=("continuous"))[6]

g<-list('vector',4)
counter = 0
for (dv in c(2,3,5,4)){

  x$dv<-x[,dv]

  if (names(x)[dv]=='score') {dv_lab = 'Sync & Match Score [C/RMSE]'}
  if (names(x)[dv]=='cmax') {dv_lab = 'C'}
  if (names(x)[dv]=='tau') {dv_lab = 'τ'}
  if (names(x)[dv]=='rmse') {dv_lab = 'RMSE'}
  
  counter = counter + 1
  g[[counter]] <- ggplot(data=x, aes(x=trial, y=dv, colour=as.factor(task))) +
    geom_jitter(size=1, alpha=.2, width=.1, height=0) +
    #geom_line(aes(x=trial, y=0), col='black', size=1.2, alpha=.7) +
    stat_summary(aes(x=trial, y=dv, colour=as.factor(task)), fun='mean', geom='line', size=2.2, alpha=1) +
    stat_summary(aes(x=trial, y=dv, colour=as.factor(task)), geom="ribbon", fun.data=mean_cl_boot, alpha=.5) +
    theme_classic() +
    theme(panel.background = element_rect(fill = "#111111",
          colour = "#000000",size = 1, linetype = "solid")) +
    theme(legend.position="top",legend.title=element_blank()) +
    #scale_x_continuous(breaks=seq(0,limit_lags,1), limits=c(0,limit_lags)) +
    #scale_y_continuous(limits=c(-.6,.4)) +
    labs(y = dv_lab) +
    labs(x = "Trial")
    g[[counter]] <- g[[counter]] + scale_colour_manual(values=colors)
}
multiplot(plotlist=g,cols=4)

if (FALSE){
  filename=paste("perf_scores_training",'_',Sys.Date(),'.png',sep='')
  png(filename=filename,width=16,height=4,units="in",res=300)
  multiplot(plotlist=g,cols=4)
  dev.off()
}



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
