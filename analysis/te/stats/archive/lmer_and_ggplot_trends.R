library(lme4)
library(texreg)
library(lattice)

setwd('/home/dobri/logos/c3/umcg_hidden_dynamics/analysis')
X<-read.csv('te_trial_shuffle_surrs_14pps_20210214-202026.csv',sep=',')

X$PP <- factor(X$PP)
X$Task <- factor(X$Task)

boxplot(TEstim~Task,data=X)
boxplot(TEpp~Task,data=X)
boxplot(taustim~Task,data=X)
boxplot(taupp~Task,data=X)
boxplot(interactdelay12~Task,data=X)
boxplot(interactdelay21~Task,data=X)

m1=lmer(TEpp ~ 1 + (1|PP),data=X,REML=0)
m2=lmer(TEpp ~ 1 + trial + (1|PP),data=X,REML=0)
m3=lmer(TEpp ~ 1 + trial + Task + (1|PP),data=X,REML=0)
m4=lmer(TEpp ~ 1 + trial*Task + (1|PP),data=X,REML=0)

anova(m1,m2,m3,m4)
summary(m3)


library(wesanderson)
library(ggplot2)
library(ghibli)
source('~/mcmc/collective_intelligence/square/analysis_8/ranalysis/multiplot.R')

colors<-ghibli_palette("PonyoMedium",7,type=("continuous"))[c(3,5,6)]
# colors[1]<-ghibli_palette("MarnieDark2",7,type=("continuous"))[6]
# colors[2]<-ghibli_palette("MononokeDark",7,type=("continuous"))[5]
# colors[3]<-ghibli_palette("YesterdayDark",7,type=("continuous"))[6]
# Autumn color palette!
colors[1]<-ghibli_palette("MarnieMedium2",7,type=("continuous"))[6]
colors[2]<-ghibli_palette("MononokeMedium",7,type=("continuous"))[5]
colors[3]<-ghibli_palette("YesterdayMedium",7,type=("continuous"))[6]

# No. The next one.
g<-list('vector',4)
counter = 0
for (dv in c(9,10,27,28)){
  
  X$dv<-X[,dv]
  
  if (names(X)[dv]=='TEstim') {dv_lab = 'TE stim->pp'}
  if (names(X)[dv]=='TEpp') {dv_lab = 'TE pp->stim'}
  if (names(X)[dv]=='surrprop12') {dv_lab = 'Surr Test 1->2 %'}
  if (names(X)[dv]=='surrprop21') {dv_lab = 'Surr Test 2->1 %'}
  
  counter = counter + 1
  g[[counter]] <- ggplot(data=X, aes(x=trial, y=dv, colour=as.factor(Task))) +
    geom_jitter(size=1, alpha=.2, width=.0, height=.0) +
    #geom_line(aes(x=trial, y=0), col='black', size=1.2, alpha=.7) +
    stat_summary(aes(x=trial, y=dv, colour=as.factor(Task)), fun='mean', geom='line', size=2.2, alpha=1) +
    stat_summary(aes(x=trial, y=dv, colour=as.factor(Task)), geom="ribbon", fun.data=mean_cl_boot, alpha=.5) +
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
  filename=paste("te1_scores_training",'_',Sys.Date(),'.png',sep='')
  png(filename=filename,width=16,height=4,units="in",res=300)
  multiplot(plotlist=g,cols=4)
  dev.off()
}


g<-list('vector',4)
counter = 0
for (dv in c(9,10,27,28)){
  
  X$dv<-X[,dv]
  
  if (names(X)[dv]=='TEstim') {dv_lab = 'TE stim->pp'}
  if (names(X)[dv]=='TEpp') {dv_lab = 'TE pp->stim'}
  if (names(X)[dv]=='surrprop12') {dv_lab = 'Surr Test 1->2 %'}
  if (names(X)[dv]=='surrprop21') {dv_lab = 'Surr Test 2->1 %'}
  
  m1=lmer(dv ~ 1 + (1|PP),data=X,REML=0)
  m2=lmer(dv ~ 1 + trial + (1|PP),data=X,REML=0)
  m3=lmer(dv ~ 1 + trial+Task + (1|PP),data=X,REML=0)
  m4=lmer(dv ~ 1 + trial*Task + (1|PP),data=X,REML=0)
  
  if (dv==9) {X$fitted <- getME(m4,'X') %*% fixef(m4)}
  if (dv==10) {X$fitted <- getME(m4,'X') %*% fixef(m4)}
  if (dv==27) {X$fitted <- getME(m4,'X') %*% fixef(m4)}
  if (dv==28) {X$fitted <- getME(m4,'X') %*% fixef(m4)}

  counter = counter + 1
  g[[counter]] <- ggplot(data=X, aes(x=trial, y=dv, colour=as.factor(Task))) +
    geom_jitter(size=1, alpha=.5, width=.0, height=.0) +
    geom_line(aes(x=trial, y=fitted, colour=Task), size=1.2, alpha=.7) +
    #geom_line(aes(x=trial, y=0), col='black', size=1.2, alpha=.7) +
    theme_classic() +
    theme(panel.background = element_rect(fill = "#111111",
                colour = "#000000", size = 1, linetype = "solid")) +
    theme(legend.position="top", legend.title=element_blank()) +
    labs(y = dv_lab) +
    labs(x = "Trial")
  g[[counter]] <- g[[counter]] + scale_colour_manual(values=colors)
}
multiplot(plotlist=g,cols=4)

if (FALSE){
  filename=paste("te_scores_training",'_',Sys.Date(),'.png',sep='')
  png(filename=filename,width=16,height=4,units="in",res=300)
  multiplot(plotlist=g,cols=4)
  dev.off()
}
