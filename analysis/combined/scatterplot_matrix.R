# Libraries for the graphing part
library(ggplot2)
library(wesanderson)
library(plyr)
library(tidyverse)


panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19, col = colour_palette[as.factor(x_task$Training)], cex = 2)
}

x_delta_ave_wide <- read.csv("Scores_TE_2024-05-26.csv_cleaned_.csv__delta_ave_wide_.csv")
names(x_delta_ave_wide) <- c("pp","Task Label","Training","c Delta% PostTest","c Delta% Retention",
                          "Pitch Error Delta% PostTest","Pitch Error Delta% Retention",
                          "Score Delta% PostTest","Score Delta% Retention",
                          "TE S-U Delta% PostTest","TE S-U Delta% Retention",
                          "TE U-S Delta% PostTest","TE U-S Delta% Retention",
                          "TE S-U Training","TE U-S Training")

x_delta_ave_wide2 <- read.csv("Scores_MaxLyap_2024-06-07.csv_cleaned_.csv__delta_ave_wide_.csv")
names(x_delta_ave_wide2) <- c("pp","Task Label","Training","c Delta% PostTest","c Delta% Retention",
                          "Pitch Error Delta% PostTest","Pitch Error Delta% Retention",
                          "Score Delta% PostTest","Score Delta% Retention",
                          "MLS Stimulus Delta% PostTest","MLS Stimulus Delta% Retention",
                          "MLL Stimulus Delta% PostTest","MLL Stimulus Delta% Retention",
                          "MLS User Delta% PostTest","MLS User Delta% Retention",
                          "MLL User Delta% PostTest","MLL User Delta% Retention",
                          "MLS Stimulus Training","MLL Stimulus Training",
                          "MLS User Training","MLL User Training")

X <- x_delta_ave_wide[,c(1:2,10:15)] %>% left_join(., x_delta_ave_wide2, by = c("pp","Task Label"))
rm(x_delta_ave_wide,x_delta_ave_wide2)

# The dependent variables that we're potentially interested in
dvs<-c(10:15,7:8,24:25,26:27)
# The colors of the bars
colour_palette <- wes_palette("FantasticFox1",10,type=("continuous"))[c(5,8,9)]


save_figs <- TRUE
task_list <- c("Chaotic", "Chaotic Visual", "Periodic Fadeout (SCT)", "Periodic Fixed")
# Plot the data
for (t in task_list) {
  x_task <- X[which(X$'Task Label'==t), ] #segment only the data for one task
  if (save_figs){ #save_figs above needs to be true
    filename=paste("DVs_corr_mat_",gsub(" ", "_", t),'_',Sys.Date(),'.png',sep='')
    png(filename=filename,width=20,height=20,units="in",res=400)
  }
  # Create the plots
  pairs(x_task[,dvs],
        lower.panel = panel.cor,
        # lower.panel = NULL,
        upper.panel = upper.panel)
  if (save_figs){
    dev.off()
  }
}