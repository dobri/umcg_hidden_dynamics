setwd('~/logos/c3/umcg_hidden_dynamics/matlyap/')
x<-read.csv('scores_lyapunov-14-Oct-2020.csv',sep=',')
summary(x)
x$pp<-factor(x$pp)

# Name the conditions
trial_counter = 0
for (n in seq(1,dim(x)[1])) {
  if (x$condition[n]==25) x$task[n]='coupled_unstable'
  if (x$condition[n]==30) x$task[n]='uncoupl_unstable'
  if (x$condition[n]==10) x$task[n]='uncoupl_periodic'
}


# Linear mixed-effects models
library(lme4)
library(texreg)

for (dv in c(5,6,7,8)){
  x$dv<-x[,dv]
  
  if (names(x)[dv]=='LyapsTutor') {dv_lab = 'Max-Lyapunov-short-Tutor'}
  if (names(x)[dv]=='LyaplTutor') {dv_lab = 'Max-Lyapunov-long-Tutor'}
  if (names(x)[dv]=='LyapsTrainee') {dv_lab = 'Max-Lyapunov-short-Trainee'}
  if (names(x)[dv]=='LyaplTrainee') {dv_lab = 'Max-Lyapunov-long-Trainee'}
  
  # m1=lmer(dv ~ 1 + (-1+trial|pp),data=x,REML=0)
  # m2=lmer(dv ~ 1 + trial + (-1+trial|pp),data=x,REML=0)
  # m3=lmer(dv ~ 1 + trial+task + (-1+trial|pp),data=x,REML=0)
  # m4=lmer(dv ~ 1 + trial*task + (-1+trial|pp),data=x,REML=0)
  m1=lmer(dv ~ 1 + (1|pp),data=x,REML=0)
  m2=lmer(dv ~ 1 + trial + (1|pp),data=x,REML=0)
  m3=lmer(dv ~ 1 + trial+task + (1|pp),data=x,REML=0)
  m4=lmer(dv ~ 1 + trial*task + (1|pp),data=x,REML=0)
  # These don't seem to work.
  #m5a=lmer(dv ~ 1 + trial*task + (1 + trial|pp),data=x,REML=0)
  #m5b=lmer(dv ~ 1 + trial*task + (1|pp) + (trial|pp),data=x,REML=0)
  
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
source('~/logos/c3/umcg_hidden_dynamics/matlyap/multiplot.R')
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
for (dv in c(5,6,7,8)){

  x$dv<-x[,dv]

  if (names(x)[dv]=='LyapsTutor') {dv_lab = 'Max-Lyapunov-short-Tutor'}
  if (names(x)[dv]=='LyaplTutor') {dv_lab = 'Max-Lyapunov-long-Tutor'}
  if (names(x)[dv]=='LyapsTrainee') {dv_lab = 'Max-Lyapunov-short-Trainee'}
  if (names(x)[dv]=='LyaplTrainee') {dv_lab = 'Max-Lyapunov-long-Trainee'}
  
  counter = counter + 1
  g[[counter]] <- ggplot(data=x, aes(x=trial, y=dv, colour=as.factor(task))) +
    geom_jitter(size=1, alpha=.2, width=0, height=0) +
    #geom_line(aes(x=trial, y=0), col='black', size=1.2, alpha=.7) +
    stat_summary(aes(x=trial, y=dv, colour=as.factor(task)), fun='mean', geom='line', size=2.2, alpha=.7) +
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
multiplot(plotlist=g,cols=2)

if (FALSE){
  filename=paste("lyap_in_training",'_',Sys.Date(),'.png',sep='')
  png(filename=filename,width=8,height=6,units="in",res=300)
  multiplot(plotlist=g,cols=2)
  dev.off()
}



g<-list('vector',4)
counter = 0
for (dv in c(5,6,7,8)){
  
  x$dv<-x[,dv]
  
  if (names(x)[dv]=='LyapsTutor') {dv_lab = 'Max-Lyapunov-short-Tutor'}
  if (names(x)[dv]=='LyaplTutor') {dv_lab = 'Max-Lyapunov-long-Tutor'}
  if (names(x)[dv]=='LyapsTrainee') {dv_lab = 'Max-Lyapunov-short-Trainee'}
  if (names(x)[dv]=='LyaplTrainee') {dv_lab = 'Max-Lyapunov-long-Trainee'}
  
  m1=lmer(dv ~ 1 + (1|pp),data=x,REML=0)
  m2=lmer(dv ~ 1 + trial + (1|pp),data=x,REML=0)
  m3=lmer(dv ~ 1 + trial+task + (1|pp),data=x,REML=0)
  m4=lmer(dv ~ 1 + trial*task + (1|pp),data=x,REML=0)

  # if (dv==5) {x$fitted <- getME(m4,'X') %*% fixef(m4)} 
  # if (dv==6) {x$fitted <- getME(m4,'X') %*% fixef(m4)}
  # if (dv==7) {x$fitted <- getME(m3,'X') %*% fixef(m3)}
  # if (dv==8) {x$fitted <- getME(m2,'X') %*% fixef(m2)}
  if (dv==5) {x$fitted <- getME(m4,'X') %*% fixef(m4)}
  if (dv==6) {x$fitted <- getME(m3,'X') %*% fixef(m3)}
  if (dv==7) {x$fitted <- getME(m3,'X') %*% fixef(m3)}
  if (dv==8) {x$fitted <- getME(m4,'X') %*% fixef(m4)}
  
  counter = counter + 1
  g[[counter]] <- ggplot(data=x, aes(x=trial, y=dv, colour=as.factor(task))) +
    geom_jitter(size=1, alpha=.5, width=.0, height=.0) +
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
multiplot(plotlist=g,cols=2)

if (FALSE){
  filename=paste("lyap_scores_training_with_lin_model",'_',Sys.Date(),'.png',sep='')
  png(filename=filename,width=8,height=6,units="in",res=300)
  multiplot(plotlist=g,cols=2)
  dev.off()
}


# λ ~ λ
#colors<-ghibli_palette("MononokeLight",7,type=("continuous"))[2:4]
colors<-ghibli_palette("LaputaMedium",7,type=("continuous"))[c(6,3,7)]
colors[3]<-ghibli_palette("MononokeMedium",7,type=("continuous"))[6]
g<-list('vector',2)
counter = 0
for (dv in c(5,6)){
  
  x$iv<-x[,dv]
  x$dv<-x[,dv+2]
  
  if (names(x)[dv]=='LyapsTutor') {
    dv_lab = 'Max-Lyapunov-short-Trainee'
    iv_lab = 'Max-Lyapunov-short-Tutor'
  }
  if (names(x)[dv]=='LyaplTutor') {
    dv_lab = 'Max-Lyapunov-long-Trainee'
    iv_lab = 'Max-Lyapunov-long-Tutor'
  }
  
  m1=lmer(dv ~ 1 + (1|pp),data=x,REML=0)
  m2=lmer(dv ~ 1 + iv + (1|pp),data=x,REML=0)
  m3=lmer(dv ~ 1 + iv+task + (1|pp),data=x,REML=0)
  m4=lmer(dv ~ 1 + iv*task + (1|pp),data=x,REML=0)
  x$fitted <- getME(m4,'X') %*% fixef(m4)
  
  counter = counter + 1
  g[[counter]] <- ggplot(data=x, aes(x=iv, y=dv, colour=as.factor(task))) +
    geom_jitter(size=1., alpha=.9, width=.0, height=.0) +
    geom_line(aes(x=iv, y=iv), col='grey', size=1., alpha=.7, linetype = "longdash") +
    geom_line(aes(x=iv, y=fitted, colour=task), size=2, alpha=1.) +
    theme_classic() +
    theme(panel.background = element_rect(fill = "#000000",
                                          colour = "#111111",size = 1, linetype = "solid")) +
    theme(legend.position="top",legend.title=element_blank()) +
    labs(y = dv_lab) +
    labs(x = iv_lab)
  g[[counter]] <- g[[counter]] + scale_colour_manual(values=colors)
}
multiplot(plotlist=g,cols=2)

if (FALSE){
  filename=paste("lyap_vs_lyap_training_with_lin_model",'_',Sys.Date(),'.png',sep='')
  png(filename=filename,width=8,height=4,units="in",res=300)
  multiplot(plotlist=g,cols=2)
  dev.off()
}



# λ ~ C
xl = x
names(xl)<-c('ppl','triall','condition','period','LyapsTutor','LyaplTutor','LyapsTrainee','LyaplTrainee','taskl')
# go training_trials_lmers_and_figs.R and run the first block of code.

plot(x$task_num-xl$condition+as.numeric(x$pp)-as.numeric(xl$pp))

x <- cbind(xl,x)

#colors<-ghibli_palette("MononokeLight",7,type=("continuous"))[2:4]
colors<-ghibli_palette("LaputaMedium",7,type=("continuous"))[c(6,3,7)]
colors[3]<-ghibli_palette("MononokeMedium",7,type=("continuous"))[6]
g<-list('vector',4)
counter = 0
for (dv in c(5,6,7,8)){
  
  x$dv<-x[,dv]
  # x$iv<-x$cmax
  x$iv<-x$score
  
  if (names(x)[dv]=='LyapsTutor') {dv_lab = 'Max-Lyapunov-short-Tutor'}
  if (names(x)[dv]=='LyaplTutor') {dv_lab = 'Max-Lyapunov-long-Tutor'}
  if (names(x)[dv]=='LyapsTrainee') {dv_lab = 'Max-Lyapunov-short-Trainee'}
  if (names(x)[dv]=='LyaplTrainee') {dv_lab = 'Max-Lyapunov-long-Trainee'}
  
  m1=lmer(dv ~ 1 + (1|pp),data=x,REML=0)
  m2=lmer(dv ~ 1 + iv + (1|pp),data=x,REML=0)
  m3=lmer(dv ~ 1 + iv+task + (1|pp),data=x,REML=0)
  m4=lmer(dv ~ 1 + iv*task + (1|pp),data=x,REML=0)
  
  if (dv==5) {x$fitted <- getME(m4,'X') %*% fixef(m4)}
  if (dv==6) {x$fitted <- getME(m4,'X') %*% fixef(m4)}
  if (dv==7) {x$fitted <- getME(m4,'X') %*% fixef(m4)}
  if (dv==8) {x$fitted <- getME(m4,'X') %*% fixef(m4)}
  
  counter = counter + 1
  g[[counter]] <- ggplot(data=x, aes(x=iv, y=dv, colour=as.factor(task))) +
    geom_jitter(size=1., alpha=.8, width=.0, height=.0) +
    geom_line(aes(x=iv, y=fitted, colour=task), size=2, alpha=1.) +
    theme_classic() +
    theme(panel.background = element_rect(fill = "#000000",
            colour = "#000000",size = 1, linetype = "solid")) +
    theme(legend.position="top",legend.title=element_blank()) +
    labs(y = dv_lab) +
    labs(x = "Sync")
  # if ((dv %% 2) == 1){
  #   g[[counter]] <- g[[counter]] + scale_y_continuous(limits=c(.3, 2.6))
  # } else {
  #   g[[counter]] <- g[[counter]] + scale_y_continuous(limits=c(-.0, .9))
  # }
  g[[counter]] <- g[[counter]] + scale_colour_manual(values=colors)
}
multiplot(plotlist=g,cols=2)

if (FALSE){
  filename=paste("lyap_vs_scores_training_with_lin_model",'_',Sys.Date(),'.png',sep='')
  png(filename=filename,width=8,height=6,units="in",res=300)
  multiplot(plotlist=g,cols=2)
  dev.off()
}
