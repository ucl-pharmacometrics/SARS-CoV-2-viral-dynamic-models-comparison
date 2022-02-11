rm(list=ls())
###### set working directory  ######
mod.dir <- "/home/akosua/Documents/Model_evaluation/SARS_CoV-2_viral_dynamic_models_comparison/Gastine_et_al_2021/Gastine_TCLE/"
setwd(mod.dir)
library(ggplot2)
library(scales)
library(deSolve)
library(nlmixr)
#
##### read data ######################################################################## 
vl_dat <- read.csv("vl_dat_14days_post_symptom.csv", h = TRUE) 
#
################ GOF plots function ######################
basic.gof <- function(run.dat){
  #
  library(ggplot2)
  library(gridExtra) 
  dv.pred <- ggplot() + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "PRED", y = "DV",
         title = paste("DV vs. PRED")) +
    geom_point(aes(x = PRED, y = DV, group = ID), 
               data = run.dat, col = grey(0.8), cex = 2) +
    geom_line(aes(x = PRED, y = DV, group = ID), 
              data = run.dat, col = grey(0.8)) +
    geom_point(aes(x = PRED, y = DV), 
               data = run.dat[run.dat$CENS == 1, ], col = "black", cex = 2) +
    geom_abline(slope = 1, intercept = 0, col = "black", lwd = 1) + 
    geom_smooth(aes(x = PRED, y = DV), 
                data = run.dat, span = 5/3, se = FALSE, col = "red", lwd = 1) +
    scale_x_continuous(limits = c(7, 18), breaks = seq(0,18,2))
#
  dv.ipred <- ggplot() + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "IPRED", y = "DV",
         title = paste("DV vs. IPRED")) +
    geom_point(aes(x = IPRED, y = DV, group = ID), 
               data = run.dat, col = grey(0.8), cex = 2) +
    geom_line(aes(x = IPRED, y = DV, group = ID), 
              data = run.dat, col = grey(0.8)) +
    geom_point(aes(x = IPRED, y = DV), 
               data = run.dat[run.dat$CENS == 1, ], col = "black", cex = 2) +
    geom_abline(slope = 1, intercept = 0, col = "black", lwd = 1) + 
    geom_smooth(aes(x = IPRED, y = DV), 
                data = run.dat, span = 5/3, se = FALSE, col = "red", lwd = 1) 
  #
  npde.tad <- ggplot() + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "Time since infection (d)", y = "NPDE",
         title = paste("NPDE vs. TIME")) +
    geom_point(aes(x = TIME, y = NPDE, group = ID), 
               data = run.dat, col = grey(0.8), cex = 2) +
    geom_line(aes(x = TIME, y = NPDE, group = ID), 
              data = run.dat, col = grey(0.8)) +
    geom_point(aes(x = TIME, y = NPDE), 
               data = run.dat[run.dat$CENS == 1, ], col = "black", cex = 2) +
    geom_hline(yintercept = c(-2, 0, 2), lty = c(2, 1, 2)) + 
    geom_smooth(aes(x = TIME, y = NPDE), 
                data = run.dat, span = 5/3, se = FALSE, col = "red", lwd = 1) 
  #
  npde.pred <- ggplot() + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "PRED", y = "NPDE",
         title = paste("NPDE vs. PRED")) +
    geom_point(aes(x = PRED, y = NPDE, group = ID), 
               data = run.dat, col = grey(0.8), cex = 2) +
    geom_line(aes(x = PRED, y = NPDE, group = ID), 
              data = run.dat, col = grey(0.8)) +
    geom_point(aes(x = PRED, y = NPDE), 
               data = run.dat[run.dat$CENS == 1, ], col = "black", cex = 2) +
    geom_hline(yintercept = c(-2, 0, 2), lty = c(2, 1, 2)) + 
    geom_smooth(aes(x = PRED, y = NPDE), 
                data = run.dat, span = 5/3, se = FALSE, col = "red", lwd = 1) +
    scale_x_continuous(limits = c(7, 18), breaks = seq(0,18,2))
  grid.arrange(dv.pred, dv.ipred, npde.tad, npde.pred, ncol = 2) 
  
}
###################### end of GOF function ####################################
#
### TCLE model ######
## start the run record:
run.rec <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(run.rec) <- c("Time since infection (d)", "dataset", "description", "OFV", "AIC", "BIC")
#
TCLE <- function() {
  ini({
    tbeta  <- -10.72       # Rate constant for virus infection
    tdelta <- -0.51         # Death rate of infected cells
    trho   <- 3.12          # Viral replication
    tc     <-  2.30        # Viral clearance  
    tk     <- 1.10        # productivity rate of infected cells
    
    
    eta.beta ~ 0.1
    eta.delta ~ 0.1
    eta.rho ~ 0.1
    eta.c ~ 0.1
    eta.k ~ 0.1
    
    add.err <- 1     # residual variability
  })
  model({
    beta  <- exp(tbeta + eta.beta)    # individual value of beta
    delta <- exp(tdelta + eta.delta)  # individual value of delta
    rho  <- exp(trho + eta.rho)      # individual value of rho
    c    <- exp(tc + eta.c)          # individual value of c
    k    <- exp(tk + eta.k)          # individual value of k
    
    #### delta= 0.6, rho= 22.7, c= 10, k=3, beta = 2.21e-05, T0 = 1.3*(10^5), V0 = 0.1), initial estimates from Dodds et al doi: 10.1111/bcp.14486  
    
    A_T(0) = 1.3*(10^5)
    A_I1(0) = 1/30
    A_I2(0) = 0
    A_v(0) = 0.1
    
    ##### A_T = 1.3*(10^5), A_I1 = 1/30 , A_I2 = 0 , A_v = 0, initial estimates from Neant et al doi:10.1073/pnas.2017962118  
    
    d/dt(A_T) <- -beta * A_v * A_T
    d/dt(A_I1) <- beta * A_v * A_T - k * A_I1
    d/dt(A_I2) <- k * A_I1 - delta *  A_I2
    d/dt(A_v) <- rho * A_I2 - c * A_v - beta*A_v*A_T
    
    
    vt = log(A_v)
    vt ~ add(add.err)       # define error model
  })
}
# Check the model
nlmixr(TCLE)
#
##### read data ######################################################################## 0.5 Day
vl_dat <- read.csv("vl_dat_14days_post_symptom.csv", h = TRUE) 
##### Create a column for time of infection (i.e., time + 0.5) 
vl_dat$time_of_infection <- vl_dat$time + 0.5
###### rename time as time_since_symptom_onset 
names(vl_dat)[names(vl_dat) == "time"] <- "time_since_symptom_onset"
###### rename time_of_infection as time
names(vl_dat)[names(vl_dat) == "time_of_infection"] <- "time"
#
## Fit the model 
run0.5.lst <- nlmixr(TCLE, vl_dat, est = "saem",
                     control = saemControl(
                       seed = 99,
                       nu = c(2,2,2)),
                     table = tableControl(npde = TRUE, censMethod = c("cdf")))
## save output to use later
saveRDS(run0.5.lst, "run0.5.rds")
run0.5.lst <- readRDS("run0.5.rds")
######## plot GOF#################################
pdf("0.5_day_gof",width = 6, height = 4)
basic.gof(run.dat = run0.5.lst)
dev.off()
#
############ VPC + GGPLOT #########################################
#### Make a VPC with nlmixr 
library(vpc)
library(ggpubr)
run_vpc <- nlmixr::vpc(run0.5.lst, show = list(obs_dv = TRUE))
# extract the simulated data
sim_data <- run_vpc$sim
# specify that time is the idv column or it gives an error message asking for idv
sim_data$idv <- sim_data$time
# then make the vpc with the vpc package:
vpc1 <- vpc::vpc(sim = sim_data,
                 obs = run0.5.lst,
                 show = list(obs_dv = TRUE),
                 lloq = 6.793,
                 facet = "rows") 
vpc2 <- vpc::vpc_cens(sim = sim_data,
                      obs = run0.5.lst,
                      lloq = 6.793,
                      facet = "rows") 
plot1 <- vpc1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "", y = "DV",
       title = paste("VPC")) + 
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(-10, 30), breaks = seq(0,30,10)) 
plot2 <- vpc2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + theme (axis.title.y = element_text(size = 10))+
  labs(x = "Time since infection (d)", y = "Probability of < LOD",
       title = paste("")) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.5)) 
###################################################### end VPC + GGPLOT
pdf("0.5_days_vpc",width = 6, height = 4)
ggarrange(plot1, plot2, nrow=2, ncol=1, heights = c(0.7, 0.3))
dev.off()
#
run.rec[1, ] <- c("0.5", "vl_dat_14days_post_symptom", "TCLE model", round(run0.5.lst$OBJF),round(run0.5.lst$AIC), round(run0.5.lst$BIC) )
#
#################
#
##### read data ######################################################################## 1 Day
vl_dat <- read.csv("vl_dat_14days_post_symptom.csv", h = TRUE) 
##### Create a column for time of infection (i.e., time + 1) 
vl_dat$time_of_infection <- vl_dat$time + 1
###### rename time as time_since_symptom_onset 
names(vl_dat)[names(vl_dat) == "time"] <- "time_since_symptom_onset"
###### rename time_of_infection as time
names(vl_dat)[names(vl_dat) == "time_of_infection"] <- "time"
#
## Fit the model 
run1.lst <- nlmixr(TCLE, vl_dat, est = "saem",
                   control = saemControl(
                     seed = 99,
                     nu = c(2,2,2)),
                   table = tableControl(npde = TRUE, censMethod = c("cdf")))
## save output to use later
saveRDS(run1.lst, "run1.rds")
run1.lst <- readRDS("run1.rds")
######## plot GOF#################################
pdf("1_day_gof",width = 6, height = 4)
basic.gof(run.dat = run1.lst)
dev.off()
#
############ VPC + GGPLOT #########################################
#### Make a VPC with nlmixr 
library(vpc)
library(ggpubr)
run_vpc <- nlmixr::vpc(run1.lst, show = list(obs_dv = TRUE))
# extract the simulated data
sim_data <- run_vpc$sim
# specify that time is the idv column or it gives an error message asking for idv
sim_data$idv <- sim_data$time
# then make the vpc with the vpc package:
vpc1 <- vpc::vpc(sim = sim_data,
                 obs = run1.lst,
                 show = list(obs_dv = TRUE),
                 lloq = 6.793,
                 facet = "rows") 
vpc2 <- vpc::vpc_cens(sim = sim_data,
                      obs = run1.lst,
                      lloq = 6.793,
                      facet = "rows") 
plot1 <- vpc1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "", y = "DV",
       title = paste("VPC")) + 
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(-10, 30), breaks = seq(0,30,10)) 
plot2 <- vpc2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + theme (axis.title.y = element_text(size = 10))+
  labs(x = "Time since infection (d)", y = "Probability of < LOD",
       title = paste("")) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0,0.6,0.2)) 
###################################################### end VPC + GGPLOT
pdf("1_day_vpc",width = 6, height = 4)
ggarrange(plot1, plot2, nrow=2, ncol=1, heights = c(0.7, 0.3))
dev.off()
#
run.rec[2, ] <- c("1", "vl_dat_14days_post_symptom", "TCLE model", round(run1.lst$OBJF),round(run1.lst$AIC), round(run1.lst$BIC) )
#
#################
#
##### read data ######################################################################## 1.5 Days
vl_dat <- read.csv("vl_dat_14days_post_symptom.csv", h = TRUE) 
##### Create a column for time of infection (i.e., time + 1.5) 
vl_dat$time_of_infection <- vl_dat$time + 1.5
###### rename time as time_since_symptom_onset 
names(vl_dat)[names(vl_dat) == "time"] <- "time_since_symptom_onset"
###### rename time_of_infection as time
names(vl_dat)[names(vl_dat) == "time_of_infection"] <- "time"
#
## Fit the model 
run1.5.lst <- nlmixr(TCLE, vl_dat, est = "saem",
                     control = saemControl(
                       seed = 99,
                       nu = c(2,2,2)),
                     table = tableControl(npde = TRUE, censMethod = c("cdf")))
## save output to use later
saveRDS(run1.5.lst, "run1.5.rds")
run1.5.lst <- readRDS("run1.5.rds")
######## plot GOF#################################
pdf("1.5_day_gof",width = 6, height = 4)
basic.gof(run.dat = run1.5.lst)
dev.off()
#
############ VPC + GGPLOT #########################################
#### Make a VPC with nlmixr 
library(vpc)
library(ggpubr)
run_vpc <- nlmixr::vpc(run1.5.lst, show = list(obs_dv = TRUE))
# extract the simulated data
sim_data <- run_vpc$sim
# specify that time is the idv column or it gives an error message asking for idv
sim_data$idv <- sim_data$time
# then make the vpc with the vpc package:
vpc1 <- vpc::vpc(sim = sim_data,
                 obs = run1.5.lst,
                 show = list(obs_dv = TRUE),
                 lloq = 6.793,
                 facet = "rows") 
vpc2 <- vpc::vpc_cens(sim = sim_data,
                      obs = run1.5.lst,
                      lloq = 6.793,
                      facet = "rows") 
plot1 <- vpc1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "", y = "DV",
       title = paste("VPC")) + 
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(-10, 30), breaks = seq(0,30,10)) 
plot2 <- vpc2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + theme (axis.title.y = element_text(size = 10))+
  labs(x = "Time since infection (d)", y = "Probability of < LOD",
       title = paste("")) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0,0.6,0.2)) 
###################################################### end VPC + GGPLOT
pdf("1.5_days_vpc",width = 6, height = 4)
ggarrange(plot1, plot2, nrow=2, ncol=1, heights = c(0.7, 0.3))
dev.off()
#
run.rec[3, ] <- c("1.5", "vl_dat_14days_post_symptom", "TCLE model", round(run1.5.lst$OBJF),round(run1.5.lst$AIC), round(run1.5.lst$BIC) )
#
#################
#
##### read data ######################################################################## 2 Days
vl_dat <- read.csv("vl_dat_14days_post_symptom.csv", h = TRUE) 
##### Create a column for time of infection (i.e., time + 2) 
vl_dat$time_of_infection <- vl_dat$time + 2
###### rename time as time_since_symptom_onset 
names(vl_dat)[names(vl_dat) == "time"] <- "time_since_symptom_onset"
###### rename time_of_infection as time
names(vl_dat)[names(vl_dat) == "time_of_infection"] <- "time"
#
## Fit the model 
run2.lst <- nlmixr(TCLE, vl_dat, est = "saem",
                   control = saemControl(
                     seed = 99,
                     nu = c(2,2,2)),
                   table = tableControl(npde = TRUE, censMethod = c("cdf")))
## save output to use later
saveRDS(run2.lst, "run2.rds")
run2.lst <- readRDS("run2.rds")
######## plot GOF#################################
pdf("2_day_gof",width = 6, height = 4)
basic.gof(run.dat = run2.lst)
dev.off()
#
############ VPC + GGPLOT #########################################
#### Make a VPC with nlmixr 
library(vpc)
library(ggpubr)
run_vpc <- nlmixr::vpc(run2.lst, show = list(obs_dv = TRUE))
# extract the simulated data
sim_data <- run_vpc$sim
# specify that time is the idv column or it gives an error message asking for idv
sim_data$idv <- sim_data$time
# then make the vpc with the vpc package:
vpc1 <- vpc::vpc(sim = sim_data,
                 obs = run2.lst,
                 show = list(obs_dv = TRUE),
                 lloq = 6.793,
                 facet = "rows") 
vpc2 <- vpc::vpc_cens(sim = sim_data,
                      obs = run2.lst,
                      lloq = 6.793,
                      facet = "rows") 
plot1 <- vpc1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "", y = "DV",
       title = paste("VPC")) + 
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(-10, 30), breaks = seq(0,30,10)) 
plot2 <- vpc2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + theme (axis.title.y = element_text(size = 10))+
  labs(x = "Time since infection (d)", y = "Probability of < LOD",
       title = paste("")) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0,0.6,0.2)) 
###################################################### end VPC + GGPLOT
pdf("2_days_vpc",width = 6, height = 4)
ggarrange(plot1, plot2, nrow=2, ncol=1, heights = c(0.7, 0.3))
dev.off()
#
run.rec[4, ] <- c("2", "vl_dat_14days_post_symptom", "TCLE model", round(run2.lst$OBJF),round(run2.lst$AIC), round(run2.lst$BIC) )
#
#################
#
##### read data ######################################################################## 3 Days
vl_dat <- read.csv("vl_dat_14days_post_symptom.csv", h = TRUE) 
##### Create a column for time of infection (i.e., time + 3) 
vl_dat$time_of_infection <- vl_dat$time + 3
###### rename time as time_since_symptom_onset 
names(vl_dat)[names(vl_dat) == "time"] <- "time_since_symptom_onset"
###### rename time_of_infection as time
names(vl_dat)[names(vl_dat) == "time_of_infection"] <- "time"
#
## Fit the model 
run3.lst <- nlmixr(TCLE, vl_dat, est = "saem",
                   control = saemControl(
                     seed = 99,
                     nu = c(2,2,2)),
                   table = tableControl(npde = TRUE, censMethod = c("cdf")))
## save output to use later
saveRDS(run3.lst, "run3.rds")
run3.lst <- readRDS("run3.rds")
######## plot GOF#################################
pdf("3_day_gof",width = 6, height = 4)
basic.gof(run.dat = run3.lst)
dev.off()
#
############ VPC + GGPLOT #########################################
#### Make a VPC with nlmixr 
library(vpc)
library(ggpubr)
run_vpc <- nlmixr::vpc(run3.lst, show = list(obs_dv = TRUE))
# extract the simulated data
sim_data <- run_vpc$sim
# specify that time is the idv column or it gives an error message asking for idv
sim_data$idv <- sim_data$time
# then make the vpc with the vpc package:
vpc1 <- vpc::vpc(sim = sim_data,
                 obs = run3.lst,
                 show = list(obs_dv = TRUE),
                 lloq = 6.793,
                 facet = "rows") 
vpc2 <- vpc::vpc_cens(sim = sim_data,
                      obs = run3.lst,
                      lloq = 6.793,
                      facet = "rows") 
plot1 <- vpc1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "", y = "DV",
       title = paste("VPC")) + 
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(-10, 30), breaks = seq(0,30,10)) 
plot2 <- vpc2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + theme (axis.title.y = element_text(size = 10))+
  labs(x = "Time since infection (d)", y = "Probability of < LOD",
       title = paste("")) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0,0.6,0.2)) 
###################################################### end VPC + GGPLOT
pdf("3_days_vpc",width = 6, height = 4)
ggarrange(plot1, plot2, nrow=2, ncol=1, heights = c(0.7, 0.3))
dev.off()
#
run.rec[5, ] <- c("3", "vl_dat_14days_post_symptom", "TCLE model", round(run3.lst$OBJF),round(run3.lst$AIC), round(run3.lst$BIC) )
#
#################
#
##### read data ######################################################################## 4 Days
vl_dat <- read.csv("vl_dat_14days_post_symptom.csv", h = TRUE) 
##### Create a column for time of infection (i.e., time + 4) 
vl_dat$time_of_infection <- vl_dat$time + 4
###### rename time as time_since_symptom_onset 
names(vl_dat)[names(vl_dat) == "time"] <- "time_since_symptom_onset"
###### rename time_of_infection as time
names(vl_dat)[names(vl_dat) == "time_of_infection"] <- "time"
#
## Fit the model 
run4.lst <- nlmixr(TCLE, vl_dat, est = "saem",
                   control = saemControl(
                     seed = 99,
                     nu = c(2,2,2)),
                   table = tableControl(npde = TRUE, censMethod = c("cdf")))
## save output to use later
saveRDS(run4.lst, "run4.rds")
run4.lst <- readRDS("run4.rds")
######## plot GOF#################################
pdf("4_day_gof",width = 6, height = 4)
basic.gof(run.dat = run4.lst)
dev.off()
#
############ VPC + GGPLOT #########################################
#### Make a VPC with nlmixr 
library(vpc)
library(ggpubr)
run_vpc <- nlmixr::vpc(run4.lst, show = list(obs_dv = TRUE))
# extract the simulated data
sim_data <- run_vpc$sim
# specify that time is the idv column or it gives an error message asking for idv
sim_data$idv <- sim_data$time
# then make the vpc with the vpc package:
vpc1 <- vpc::vpc(sim = sim_data,
                 obs = run4.lst,
                 show = list(obs_dv = TRUE),
                 lloq = 6.793,
                 facet = "rows") 
vpc2 <- vpc::vpc_cens(sim = sim_data,
                      obs = run4.lst,
                      lloq = 6.793,
                      facet = "rows") 
plot1 <- vpc1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "", y = "DV",
       title = paste("VPC")) + 
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(-10, 30), breaks = seq(0,30,10)) 
plot2 <- vpc2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + theme (axis.title.y = element_text(size = 10))+
  labs(x = "Time since infection (d)", y = "Probability of < LOD",
       title = paste("")) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0,0.6,0.2)) 
###################################################### end VPC + GGPLOT
pdf("4_days_vpc",width = 6, height = 4)
ggarrange(plot1, plot2, nrow=2, ncol=1, heights = c(0.7, 0.3))
dev.off()
#
run.rec[6, ] <- c("4", "vl_dat_14days_post_symptom", "TCLE model", round(run4.lst$OBJF),round(run4.lst$AIC), round(run4.lst$BIC) )
#
#################
#
##### read data ######################################################################## 5 Days
vl_dat <- read.csv("vl_dat_14days_post_symptom.csv", h = TRUE) 
##### Create a column for time of infection (i.e., time + 5) 
vl_dat$time_of_infection <- vl_dat$time + 5
###### rename time as time_since_symptom_onset 
names(vl_dat)[names(vl_dat) == "time"] <- "time_since_symptom_onset"
###### rename time_of_infection as time
names(vl_dat)[names(vl_dat) == "time_of_infection"] <- "time"
#
## Fit the model 
run5.lst <- nlmixr(TCLE, vl_dat, est = "saem",
                   control = saemControl(
                     seed = 99,
                     nu = c(2,2,2)),
                   table = tableControl(npde = TRUE, censMethod = c("cdf")))
## save output to use later
saveRDS(run5.lst, "run5.rds")
run5.lst <- readRDS("run5.rds")
######## plot GOF#################################
pdf("5_day_gof",width = 6, height = 4)
basic.gof(run.dat = run5.lst)
dev.off()
#
############ VPC + GGPLOT #########################################
#### Make a VPC with nlmixr 
library(vpc)
library(ggpubr)
run_vpc <- nlmixr::vpc(run5.lst, show = list(obs_dv = TRUE))
# extract the simulated data
sim_data <- run_vpc$sim
# specify that time is the idv column or it gives an error message asking for idv
sim_data$idv <- sim_data$time
# then make the vpc with the vpc package:
vpc1 <- vpc::vpc(sim = sim_data,
                 obs = run5.lst,
                 show = list(obs_dv = TRUE),
                 lloq = 6.793,
                 facet = "rows") 
vpc2 <- vpc::vpc_cens(sim = sim_data,
                      obs = run5.lst,
                      lloq = 6.793,
                      facet = "rows") 
plot1 <- vpc1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "", y = "DV",
       title = paste("VPC")) + 
  scale_x_continuous(limits = c(0, 20), breaks = seq(0,20,5)) +
  scale_y_continuous(limits = c(-20, 30), breaks = seq(0,30,10)) 
plot2 <- vpc2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + theme (axis.title.y = element_text(size = 10))+
  labs(x = "Time since infection (d)", y = "Probability of < LOD",
       title = paste("")) +
  scale_x_continuous(limits = c(0, 20), breaks = seq(0,20,5)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.5)) 
###################################################### end VPC + GGPLOT
pdf("5_days_vpc",width = 6, height = 4)
ggarrange(plot1, plot2, nrow=2, ncol=1, heights = c(0.7, 0.3))
dev.off()
#
run.rec[7, ] <- c("5", "vl_dat_14days_post_symptom", "TCLE model", round(run5.lst$OBJF),round(run5.lst$AIC), round(run5.lst$BIC) )
#
#################
#
##### read data ######################################################################## 6 Days
vl_dat <- read.csv("vl_dat_14days_post_symptom.csv", h = TRUE) 
##### Create a column for time of infection (i.e., time + 6) 
vl_dat$time_of_infection <- vl_dat$time + 6
###### rename time as time_since_symptom_onset 
names(vl_dat)[names(vl_dat) == "time"] <- "time_since_symptom_onset"
###### rename time_of_infection as time
names(vl_dat)[names(vl_dat) == "time_of_infection"] <- "time"
#
## Fit the model 
run6.lst <- nlmixr(TCLE, vl_dat, est = "saem",
                   control = saemControl(
                     seed = 99,
                     nu = c(2,2,2)),
                   table = tableControl(npde = TRUE, censMethod = c("cdf")))
## save output to use later
saveRDS(run6.lst, "run6.rds")
run6.lst <- readRDS("run6.rds")
######## plot GOF#################################
pdf("6_day_gof",width = 6, height = 4)
basic.gof(run.dat = run6.lst)
dev.off()
#
############ VPC + GGPLOT #########################################
#### Make a VPC with nlmixr 
library(vpc)
library(ggpubr)
run_vpc <- nlmixr::vpc(run6.lst, show = list(obs_dv = TRUE))
# extract the simulated data
sim_data <- run_vpc$sim
# specify that time is the idv column or it gives an error message asking for idv
sim_data$idv <- sim_data$time
# then make the vpc with the vpc package:
vpc1 <- vpc::vpc(sim = sim_data,
                 obs = run6.lst,
                 show = list(obs_dv = TRUE),
                 lloq = 6.793,
                 facet = "rows") 
vpc2 <- vpc::vpc_cens(sim = sim_data,
                      obs = run6.lst,
                      lloq = 6.793,
                      facet = "rows") 
plot1 <- vpc1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "", y = "DV",
       title = paste("VPC")) + 
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(-10, 30), breaks = seq(0,30,10)) 
plot2 <- vpc2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + theme (axis.title.y = element_text(size = 10))+
  labs(x = "Time since infection (d)", y = "Probability of < LOD",
       title = paste("")) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0,0.6,0.2)) 
###################################################### end VPC + GGPLOT
pdf("6_days_vpc",width = 6, height = 4)
ggarrange(plot1, plot2, nrow=2, ncol=1, heights = c(0.7, 0.3))
dev.off()
#
run.rec[8, ] <- c("6", "vl_dat_14days_post_symptom", "TCLE model", round(run6.lst$OBJF),round(run6.lst$AIC), round(run6.lst$BIC) )
#
#################
#
##### read data ######################################################################## 7 Days
vl_dat <- read.csv("vl_dat_14days_post_symptom.csv", h = TRUE) 
##### Create a column for time of infection (i.e., time + 7) 
vl_dat$time_of_infection <- vl_dat$time + 7
###### rename time as time_since_symptom_onset 
names(vl_dat)[names(vl_dat) == "time"] <- "time_since_symptom_onset"
###### rename time_of_infection as time
names(vl_dat)[names(vl_dat) == "time_of_infection"] <- "time"
#
## Fit the model 
run7.lst <- nlmixr(TCLE, vl_dat, est = "saem",
                   control = saemControl(
                     seed = 99,
                     nu = c(2,2,2)),
                   table = tableControl(npde = TRUE, censMethod = c("cdf")))
## save output to use later
saveRDS(run7.lst, "run7.rds")
run7.lst <- readRDS("run7.rds")
######## plot GOF#################################
pdf("7_day_gof",width = 6, height = 4)
basic.gof(run.dat = run7.lst)
dev.off()
#
############ VPC + GGPLOT #########################################
#### Make a VPC with nlmixr 
library(vpc)
library(ggpubr)
run_vpc <- nlmixr::vpc(run7.lst, show = list(obs_dv = TRUE))
# extract the simulated data
sim_data <- run_vpc$sim
# specify that time is the idv column or it gives an error message asking for idv
sim_data$idv <- sim_data$time
# then make the vpc with the vpc package:
vpc1 <- vpc::vpc(sim = sim_data,
                 obs = run7.lst,
                 show = list(obs_dv = TRUE),
                 lloq = 6.793,
                 facet = "rows") 
vpc2 <- vpc::vpc_cens(sim = sim_data,
                      obs = run7.lst,
                      lloq = 6.793,
                      facet = "rows") 
plot1 <- vpc1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "", y = "DV",
       title = paste("VPC")) + 
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(-10, 30), breaks = seq(0,30,10)) 
plot2 <- vpc2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + theme (axis.title.y = element_text(size = 10))+
  labs(x = "Time since infection (d)", y = "Probability of < LOD",
       title = paste("")) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0,0.6,0.2)) 
###################################################### end VPC + GGPLOT
pdf("7_days_vpc",width = 6, height = 4)
ggarrange(plot1, plot2, nrow=2, ncol=1, heights = c(0.7, 0.3))
dev.off()
#
run.rec[9, ] <- c("7", "vl_dat_14days_post_symptom", "TCLE model", round(run7.lst$OBJF),round(run7.lst$AIC), round(run7.lst$BIC) )
#
#################
#
##### read data ######################################################################## 8 Days
vl_dat <- read.csv("vl_dat_14days_post_symptom.csv", h = TRUE) 
##### Create a column for time of infection (i.e., time + 8) 
vl_dat$time_of_infection <- vl_dat$time + 8
###### rename time as time_since_symptom_onset 
names(vl_dat)[names(vl_dat) == "time"] <- "time_since_symptom_onset"
###### rename time_of_infection as time
names(vl_dat)[names(vl_dat) == "time_of_infection"] <- "time"
#
## Fit the model 
run8.lst <- nlmixr(TCLE, vl_dat, est = "saem",
                   control = saemControl(
                     seed = 99,
                     nu = c(2,2,2)),
                   table = tableControl(npde = TRUE, censMethod = c("cdf")))
## save output to use later
saveRDS(run8.lst, "run8.rds")
run8.lst <- readRDS("run8.rds")
######## plot GOF#################################
pdf("8_day_gof",width = 6, height = 4)
basic.gof(run.dat = run8.lst)
dev.off()
#
############ VPC + GGPLOT #########################################
#### Make a VPC with nlmixr 
library(vpc)
library(ggpubr)
run_vpc <- nlmixr::vpc(run8.lst, show = list(obs_dv = TRUE))
# extract the simulated data
sim_data <- run_vpc$sim
# specify that time is the idv column or it gives an error message asking for idv
sim_data$idv <- sim_data$time
# then make the vpc with the vpc package:
vpc1 <- vpc::vpc(sim = sim_data,
                 obs = run8.lst,
                 show = list(obs_dv = TRUE),
                 lloq = 6.793,
                 facet = "rows") 
vpc2 <- vpc::vpc_cens(sim = sim_data,
                      obs = run8.lst,
                      lloq = 6.793,
                      facet = "rows") 
plot1 <- vpc1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "", y = "DV",
       title = paste("VPC")) + 
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(-10, 30), breaks = seq(0,30,10)) 
plot2 <- vpc2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + theme (axis.title.y = element_text(size = 10))+
  labs(x = "Time since infection (d)", y = "Probability of < LOD",
       title = paste("")) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0,0.6,0.2)) 
###################################################### end VPC + GGPLOT
pdf("8_days_vpc",width = 6, height = 4)
ggarrange(plot1, plot2, nrow=2, ncol=1, heights = c(0.7, 0.3))
dev.off()
#
run.rec[10, ] <- c("8", "vl_dat_14days_post_symptom", "TCLE model", round(run8.lst$OBJF),round(run8.lst$AIC), round(run8.lst$BIC) )
#
#################
#
##### read data ######################################################################## 9 Days
vl_dat <- read.csv("vl_dat_14days_post_symptom.csv", h = TRUE) 
##### Create a column for time of infection (i.e., time + 9) 
vl_dat$time_of_infection <- vl_dat$time + 9
###### rename time as time_since_symptom_onset 
names(vl_dat)[names(vl_dat) == "time"] <- "time_since_symptom_onset"
###### rename time_of_infection as time
names(vl_dat)[names(vl_dat) == "time_of_infection"] <- "time"
#
## Fit the model 
run9.lst <- nlmixr(TCLE, vl_dat, est = "saem",
                   control = saemControl(
                     seed = 99,
                     nu = c(2,2,2)),
                   table = tableControl(npde = TRUE, censMethod = c("cdf")))
## save output to use later
saveRDS(run9.lst, "run9.rds")
run9.lst <- readRDS("run9.rds")
######## plot GOF#################################
pdf("9_day_gof",width = 6, height = 4)
basic.gof(run.dat = run9.lst)
dev.off()
#
############ VPC + GGPLOT #########################################
#### Make a VPC with nlmixr 
library(vpc)
library(ggpubr)
run_vpc <- nlmixr::vpc(run9.lst, show = list(obs_dv = TRUE))
# extract the simulated data
sim_data <- run_vpc$sim
# specify that time is the idv column or it gives an error message asking for idv
sim_data$idv <- sim_data$time
# then make the vpc with the vpc package:
vpc1 <- vpc::vpc(sim = sim_data,
                 obs = run9.lst,
                 show = list(obs_dv = TRUE),
                 lloq = 6.793,
                 facet = "rows") 
vpc2 <- vpc::vpc_cens(sim = sim_data,
                      obs = run9.lst,
                      lloq = 6.793,
                      facet = "rows") 
plot1 <- vpc1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "", y = "DV",
       title = paste("VPC")) + 
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(-10, 30), breaks = seq(0,30,10)) 
plot2 <- vpc2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + theme (axis.title.y = element_text(size = 10))+
  labs(x = "Time since infection (d)", y = "Probability of < LOD",
       title = paste("")) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0,0.6,0.2)) 
###################################################### end VPC + GGPLOT
pdf("9_days_vpc",width = 6, height = 4)
ggarrange(plot1, plot2, nrow=2, ncol=1, heights = c(0.7, 0.3))
dev.off()
#
run.rec[11, ] <- c("9", "vl_dat_14days_post_symptom", "TCLE model", round(run9.lst$OBJF),round(run9.lst$AIC), round(run9.lst$BIC) )
#
#################
#
##### read data ######################################################################## 10 Days
vl_dat <- read.csv("vl_dat_14days_post_symptom.csv", h = TRUE) 
##### Create a column for time of infection (i.e., time + 10) 
vl_dat$time_of_infection <- vl_dat$time + 10
###### rename time as time_since_symptom_onset 
names(vl_dat)[names(vl_dat) == "time"] <- "time_since_symptom_onset"
###### rename time_of_infection as time
names(vl_dat)[names(vl_dat) == "time_of_infection"] <- "time"
#
## Fit the model 
run10.lst <- nlmixr(TCLE, vl_dat, est = "saem",
                    control = saemControl(
                      seed = 99,
                      nu = c(2,2,2)),
                    table = tableControl(npde = TRUE, censMethod = c("cdf")))
## save output to use later
saveRDS(run10.lst, "run10.rds")
run10.lst <- readRDS("run10.rds")
######## plot GOF#################################
pdf("10_day_gof",width = 6, height = 4)
basic.gof(run.dat = run10.lst)
dev.off()
#
############ VPC + GGPLOT #########################################
#### Make a VPC with nlmixr 
library(vpc)
library(ggpubr)
run_vpc <- nlmixr::vpc(run10.lst, show = list(obs_dv = TRUE))
# extract the simulated data
sim_data <- run_vpc$sim
# specify that time is the idv column or it gives an error message asking for idv
sim_data$idv <- sim_data$time
# then make the vpc with the vpc package:
vpc1 <- vpc::vpc(sim = sim_data,
                 obs = run10.lst,
                 show = list(obs_dv = TRUE),
                 lloq = 6.793,
                 facet = "rows") 
vpc2 <- vpc::vpc_cens(sim = sim_data,
                      obs = run10.lst,
                      lloq = 6.793,
                      facet = "rows") 
plot1 <- vpc1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "", y = "DV",
       title = paste("VPC")) + 
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(-10, 30), breaks = seq(0,30,10)) 
plot2 <- vpc2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + theme (axis.title.y = element_text(size = 10))+
  labs(x = "Time since infection (d)", y = "Probability of < LOD",
       title = paste("")) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0,0.6,0.2)) 
###################################################### end VPC + GGPLOT
pdf("10_days_vpc",width = 6, height = 4)
ggarrange(plot1, plot2, nrow=2, ncol=1, heights = c(0.7, 0.3))
dev.off()
#
run.rec[12, ] <- c("10", "vl_dat_14days_post_symptom", "TCLE model", round(run10.lst$OBJF),round(run10.lst$AIC), round(run10.lst$BIC) )
#
###########################
#
##### read data ######################################################################## 11 Days
vl_dat <- read.csv("vl_dat_14days_post_symptom.csv", h = TRUE) 
##### Create a column for time of infection (i.e., time + 11) 
vl_dat$time_of_infection <- vl_dat$time + 11
###### rename time as time_since_symptom_onset 
names(vl_dat)[names(vl_dat) == "time"] <- "time_since_symptom_onset"
###### rename time_of_infection as time
names(vl_dat)[names(vl_dat) == "time_of_infection"] <- "time"
#
## Fit the model 
run11.lst <- nlmixr(TCLE, vl_dat, est = "saem",
                    control = saemControl(
                      seed = 99,
                      nu = c(2,2,2)),
                    table = tableControl(npde = TRUE, censMethod = c("cdf")))
## save output to use later
saveRDS(run11.lst, "run11.rds")
run11.lst <- readRDS("run11.rds")
######## plot GOF#################################
pdf("11_day_gof",width = 6, height = 4)
basic.gof(run.dat = run11.lst)
dev.off()
#
############ VPC + GGPLOT #########################################
#### Make a VPC with nlmixr 
library(vpc)
library(ggpubr)
run_vpc <- nlmixr::vpc(run11.lst, show = list(obs_dv = TRUE))
# extract the simulated data
sim_data <- run_vpc$sim
# specify that time is the idv column or it gives an error message asking for idv
sim_data$idv <- sim_data$time
# then make the vpc with the vpc package:
vpc1 <- vpc::vpc(sim = sim_data,
                 obs = run11.lst,
                 show = list(obs_dv = TRUE),
                 lloq = 6.793,
                 facet = "rows") 
vpc2 <- vpc::vpc_cens(sim = sim_data,
                      obs = run11.lst,
                      lloq = 6.793,
                      facet = "rows") 
plot1 <- vpc1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "", y = "DV",
       title = paste("VPC")) + 
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(-10, 30), breaks = seq(0,30,10)) 
plot2 <- vpc2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + theme (axis.title.y = element_text(size = 10))+
  labs(x = "Time since infection (d)", y = "Probability of < LOD",
       title = paste("")) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0,0.6,0.2)) 
###################################################### end VPC + GGPLOT
pdf("11_days_vpc",width = 6, height = 4)
ggarrange(plot1, plot2, nrow=2, ncol=1, heights = c(0.7, 0.3))
dev.off()
#
run.rec[13, ] <- c("11", "vl_dat_14days_post_symptom", "TCLE model", round(run11.lst$OBJF),round(run11.lst$AIC), round(run11.lst$BIC) )
#
###########################
#
##### read data ######################################################################## 12 Days
vl_dat <- read.csv("vl_dat_14days_post_symptom.csv", h = TRUE) 
##### Create a column for time of infection (i.e., time + 12) 
vl_dat$time_of_infection <- vl_dat$time + 12
###### rename time as time_since_symptom_onset 
names(vl_dat)[names(vl_dat) == "time"] <- "time_since_symptom_onset"
###### rename time_of_infection as time
names(vl_dat)[names(vl_dat) == "time_of_infection"] <- "time"
#
## Fit the model 
run12.lst <- nlmixr(TCLE, vl_dat, est = "saem",
                    control = saemControl(
                      seed = 99,
                      nu = c(2,2,2)),
                    table = tableControl(npde = TRUE, censMethod = c("cdf")))
## save output to use later
saveRDS(run12.lst, "run12.rds")
run12.lst <- readRDS("run12.rds")
######## plot GOF#################################
pdf("12_day_gof",width = 6, height = 4)
basic.gof(run.dat = run12.lst)
dev.off()
#
############ VPC + GGPLOT #########################################
#### Make a VPC with nlmixr 
library(vpc)
library(ggpubr)
run_vpc <- nlmixr::vpc(run12.lst, show = list(obs_dv = TRUE))
# extract the simulated data
sim_data <- run_vpc$sim
# specify that time is the idv column or it gives an error message asking for idv
sim_data$idv <- sim_data$time
# then make the vpc with the vpc package:
vpc1 <- vpc::vpc(sim = sim_data,
                 obs = run12.lst,
                 show = list(obs_dv = TRUE),
                 lloq = 6.793,
                 facet = "rows") 
vpc2 <- vpc::vpc_cens(sim = sim_data,
                      obs = run12.lst,
                      lloq = 6.793,
                      facet = "rows") 
plot1 <- vpc1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "", y = "DV",
       title = paste("VPC")) + 
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(-10, 30), breaks = seq(0,30,10)) 
plot2 <- vpc2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + theme (axis.title.y = element_text(size = 10))+
  labs(x = "Time since infection (d)", y = "Probability of < LOD",
       title = paste("")) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0,0.6,0.2)) 
###################################################### end VPC + GGPLOT
pdf("12_days_vpc",width = 6, height = 4)
ggarrange(plot1, plot2, nrow=2, ncol=1, heights = c(0.7, 0.3))
dev.off()
#
run.rec[14, ] <- c("12", "vl_dat_14days_post_symptom", "TCLE model", round(run12.lst$OBJF),round(run12.lst$AIC), round(run12.lst$BIC) )
#
###########################
#
##### read data ######################################################################## 13 Days
vl_dat <- read.csv("vl_dat_14days_post_symptom.csv", h = TRUE) 
##### Create a column for time of infection (i.e., time + 13) 
vl_dat$time_of_infection <- vl_dat$time + 13
###### rename time as time_since_symptom_onset 
names(vl_dat)[names(vl_dat) == "time"] <- "time_since_symptom_onset"
###### rename time_of_infection as time
names(vl_dat)[names(vl_dat) == "time_of_infection"] <- "time"
#
## Fit the model 
run13.lst <- nlmixr(TCLE, vl_dat, est = "saem",
                    control = saemControl(
                      seed = 99,
                      nu = c(2,2,2)),
                    table = tableControl(npde = TRUE, censMethod = c("cdf")))
## save output to use later
saveRDS(run13.lst, "run13.rds")
run13.lst <- readRDS("run13.rds")
######## plot GOF#################################
pdf("13_day_gof",width = 6, height = 4)
basic.gof(run.dat = run13.lst)
dev.off()
#
############ VPC + GGPLOT #########################################
#### Make a VPC with nlmixr 
library(vpc)
library(ggpubr)
run_vpc <- nlmixr::vpc(run13.lst, show = list(obs_dv = TRUE))
# extract the simulated data
sim_data <- run_vpc$sim
# specify that time is the idv column or it gives an error message asking for idv
sim_data$idv <- sim_data$time
# then make the vpc with the vpc package:
vpc1 <- vpc::vpc(sim = sim_data,
                 obs = run13.lst,
                 show = list(obs_dv = TRUE),
                 lloq = 6.793,
                 facet = "rows") 
vpc2 <- vpc::vpc_cens(sim = sim_data,
                      obs = run13.lst,
                      lloq = 6.793,
                      facet = "rows") 
plot1 <- vpc1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "", y = "DV",
       title = paste("VPC")) + 
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(-10, 30), breaks = seq(0,30,10)) 
plot2 <- vpc2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + theme (axis.title.y = element_text(size = 10))+
  labs(x = "Time since infection (d)", y = "Probability of < LOD",
       title = paste("")) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0,0.6,0.2)) 
###################################################### end VPC + GGPLOT
pdf("13_days_vpc",width = 6, height = 4)
ggarrange(plot1, plot2, nrow=2, ncol=1, heights = c(0.7, 0.3))
dev.off()
#
run.rec[15, ] <- c("13", "vl_dat_14days_post_symptom", "TCLE model", round(run13.lst$OBJF),round(run13.lst$AIC), round(run13.lst$BIC) )
#
##########################
#
##### read data ######################################################################## 14 Days
vl_dat <- read.csv("vl_dat_14days_post_symptom.csv", h = TRUE) 
##### Create a column for time of infection (i.e., time + 14) 
vl_dat$time_of_infection <- vl_dat$time + 14
###### rename time as time_since_symptom_onset 
names(vl_dat)[names(vl_dat) == "time"] <- "time_since_symptom_onset"
###### rename time_of_infection as time
names(vl_dat)[names(vl_dat) == "time_of_infection"] <- "time"
#
## Fit the model 
run14.lst <- nlmixr(TCLE, vl_dat, est = "saem",
                    control = saemControl(
                      seed = 99,
                      nu = c(2,2,2)),
                    table = tableControl(npde = TRUE, censMethod = c("cdf")))
## save output to use later
saveRDS(run14.lst, "run14.rds")
run14.lst <- readRDS("run14.rds")
######## plot GOF#################################
pdf("14_day_gof",width = 6, height = 4)
basic.gof(run.dat = run14.lst)
dev.off()
#
############ VPC + GGPLOT #########################################
#### Make a VPC with nlmixr 
library(vpc)
library(ggpubr)
run_vpc <- nlmixr::vpc(run14.lst, show = list(obs_dv = TRUE))
# extract the simulated data
sim_data <- run_vpc$sim
# specify that time is the idv column or it gives an error message asking for idv
sim_data$idv <- sim_data$time
# then make the vpc with the vpc package:
vpc1 <- vpc::vpc(sim = sim_data,
                 obs = run14.lst,
                 show = list(obs_dv = TRUE),
                 lloq = 6.793,
                 facet = "rows") 
vpc2 <- vpc::vpc_cens(sim = sim_data,
                      obs = run14.lst,
                      lloq = 6.793,
                      facet = "rows") 
plot1 <- vpc1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "", y = "DV",
       title = paste("VPC")) + 
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(-10, 30), breaks = seq(0,30,10)) 
plot2 <- vpc2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + theme (axis.title.y = element_text(size = 10))+
  labs(x = "Time since infection (d)", y = "Probability of < LOD",
       title = paste("")) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0,0.6,0.2)) 
###################################################### end VPC + GGPLOT
pdf("14_days_vpc",width = 6, height = 4)
ggarrange(plot1, plot2, nrow=2, ncol=1, heights = c(0.7, 0.3))
dev.off()
#
run.rec[16, ] <- c("14", "vl_dat_14days_post_symptom", "TCLE model", round(run14.lst$OBJF),round(run14.lst$AIC), round(run14.lst$BIC) )
#
## save run record to use later
saveRDS(run.rec, "runrec.rds")
run.rec <- readRDS("runrec.rds")
#
## Export run record output to csv ######################
write.csv(run.rec, "run_record.csv", row.names = F, quote = F)
## Export run record output to table format csv ######################
write.table(run.rec, "run_record_2.csv", sep = "\t",  row.names = F, quote = F)
#################
#
#### Calculations from final estimates in manuscript (Table S2). 
### Basic reproduction number (R0)
rho <- 2.96*10^3
beta <- 7.76*10^-5
A_T <- 1.3*10^5
delta <- 0.69
c <- 11.6
#
R0 <- (rho*beta*A_T)/(delta*(c + beta*A_T)) # 1995.39
#
### Duration of virus production (L)
L <- 1/delta # 1.45
####### end ###############
#
############### R0 with fixed parameters #######################
###### TCLERo ######
# Parameter estimates from Kamal et al doi:10.1128/AAC.00069-15
TCLERo <- function() {
  ini({
    tbeta  <- -10.72       # Rate constant for virus infection
    tdelta <- -0.51         # Death rate of infected cells
    #trho   <- 3.12          # Viral replication
    tc     <-  2.30        # Viral clearance  
    tk     <- 1.10        # productivity rate of infected cells
    
    
    eta.beta ~ 0.1
    eta.delta ~ 0.1
    #eta.rho ~ 0.1
    eta.c ~ 0.1
    eta.k ~ 0.1
    
    add.err <- 1     # residual variability
  })
  model({
    beta  <- exp(tbeta + eta.beta)    # individual value of beta
    delta <- exp(tdelta + eta.delta)  # individual value of delta
    #rho  <- exp(trho + eta.rho)      # individual value of rho
    c    <- exp(tc + eta.c)          # individual value of c
    k    <- exp(tk + eta.k)          # individual value of k
    rho <- 10
   
    
    #### delta= 0.6, rho= 22.7, c= 10, k=3, beta = 2.21e-05, T0 = 1.3*(10^5), V0 = 0.1), initial estimates from Dodds et al doi: 10.1111/bcp.14486  
    
    A_T(0) = 1.3*(10^5)
    A_I1(0) = 1/30
    A_I2(0) = 0
    A_v(0) = 0.1
    
    ##### A_T = 1.3*(10^5), A_I1 = 1/30 , A_I2 = 0 , A_v = 0, initial estimates from Neant et al doi:10.1073/pnas.2017962118  
    
    d/dt(A_T) <- -beta * A_v * A_T
    d/dt(A_I1) <- beta * A_v * A_T - k * A_I1
    d/dt(A_I2) <- k * A_I1 - delta *  A_I2
    d/dt(A_v) <- rho * A_I2 - c * A_v - beta*A_v*A_T
    
    
    vt = log(A_v)
    vt ~ add(add.err)       # define error model
  })
}
# Check the model
nlmixr(TCLERo)
#
##### read data ######################################################################## 0.5 Day Ro
vl_dat <- read.csv("vl_dat_14days_post_symptom.csv", h = TRUE) 
##### Create a column for time of infection (i.e., time + 0.5) 
vl_dat$time_of_infection <- vl_dat$time + 0.5
###### rename time as time_since_symptom_onset 
names(vl_dat)[names(vl_dat) == "time"] <- "time_since_symptom_onset"
###### rename time_of_infection as time
names(vl_dat)[names(vl_dat) == "time_of_infection"] <- "time"
#
## Fit the model 
run0.5.lstRo <- nlmixr(TCLERo, vl_dat, est = "saem",
                       control = saemControl(
                         seed = 99,
                         nu = c(2,2,2)),
                       table = tableControl(npde = TRUE, censMethod = c("cdf")))
#
### Basic reproduction number (R0)
rho <- 10
beta <- 0.000228
A_T <- 1.3*10^5
delta <- 0.43
c <- 1.4
#
R0 <- (rho*beta*A_T)/(delta*(c + beta*A_T)) # c_fixed (10) = 1727.54 , c/b_fixed (10/0.0005) = 2358.38 , rho_fixed (10) = 22.21
#
### Duration of virus production (L)
L <- 1/delta # c_fixed (10) = 1.43, c/b_fixed (10/0.0005) = 1.67 , rho_fixed (10) = 2.33
#
## save output to use later for rho_fixed = 10
saveRDS(run0.5.lstRo, "run0.5.lstRo")
run0.5.lstRo <- readRDS("run0.5.lstRo")
## OBF, AIC , BIC for rho_fixed = 10
run0.5.lstRo$OBJF # 3329.951
run0.5.lstRo$AIC # 5183.99
run0.5.lstRo$BIC # 5228.151
#
#### end ####



