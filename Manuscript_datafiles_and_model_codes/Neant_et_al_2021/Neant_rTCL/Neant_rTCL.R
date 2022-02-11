rm(list=ls())
###### set working directory  ######
mod.dir <- "/home/akosua/Documents/Model_evaluation/SARS_CoV-2_viral_dynamic_models_comparison/Neant_et_al_2021/Neant_rTCL/"
setwd(mod.dir)
library(ggplot2)
library(scales)
library(deSolve)
library(nlmixr)
##### read data ######
vl_dat <- read.csv("vl_dat_14days_post_symptom_under65.csv", h = TRUE)
#
########### GOF function plots ####################################################
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
    scale_x_continuous(limits = c(9, 17), breaks = seq(0,17,2)) +
    scale_y_continuous(limits = c(-10, 20), breaks = seq(0,20,10))
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
                data = run.dat, span = 5/3, se = FALSE, col = "red", lwd = 1) +
    scale_x_continuous(limits = c(0, 20), breaks = seq(0,20,5)) +
    scale_y_continuous(limits = c(-10, 20), breaks = seq(0,20,10))
  #
  npde.tad <- ggplot() + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "Time since symptom onset (d)", y = "NPDE",
         title = paste("NPDE vs. TIME")) +
    geom_point(aes(x = TIME, y = NPDE, group = ID), 
               data = run.dat, col = grey(0.8), cex = 2) +
    geom_line(aes(x = TIME, y = NPDE, group = ID), 
              data = run.dat, col = grey(0.8)) +
    geom_point(aes(x = TIME, y = NPDE), 
               data = run.dat[run.dat$CENS == 1, ], col = "black", cex = 2) +
    geom_hline(yintercept = c(-2, 0, 2), lty = c(2, 1, 2)) + 
    geom_smooth(aes(x = TIME, y = NPDE), 
                data = run.dat, span = 5/3, se = FALSE, col = "red", lwd = 1) +
    scale_x_continuous(limits = c(0, 16), breaks = seq(0,16,4))
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
    scale_x_continuous(limits = c(9, 17), breaks = seq(0,17,2))
  grid.arrange(dv.pred, dv.ipred, npde.tad, npde.pred, ncol = 2) 
}
################## end of function
#
###### rTCL model (Gamma fixed to 1)  ######
## start the run record:
run.rec <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(run.rec) <- c("Dataset", "Model", "OFV", "AIC", "BIC")
# Parameter estimates from Gastine et al: doi: 10.1002/cpt.2223
rTCL.gam.fix <- function() {
  ini({
    tbeta  <- -6.89 # Rate constant for virus infection
    tdelta <- -0.5  # Death rate of infected cells
    tv0    <- 13    # Viral load at symptom onset
    
    eta.beta +  eta.delta + eta.v0  ~ c(1,
                                        0.01, 1,
                                        0.01, 0.01, 1)

    add.err <- 1     # residual variability
    
  })
  model({
    beta  <- exp(tbeta + eta.beta)    # individual value of beta
    delta <- exp(tdelta + eta.delta)  # individual value of delta
    v0    <- exp(tv0 + eta.v0)           # individual value of v0
    gamma <- 1  # T(0)*p/c, fixed value   
    
    A_f(0) = 1
    A_v(0) = v0
    
    d/dt(A_f) = -beta * A_f * A_v
    d/dt(A_v) = gamma * A_f * A_v - delta * A_v 
    
    vt = log(A_v)
    vt ~ add(add.err)       # define error model
  })
}
# Check the model
nlmixr(rTCL.gam.fix)
## Fit the model 
run1.lst <- nlmixr(rTCL.gam.fix, vl_dat , est = "saem",
                   control = saemControl(
                     seed = 99,
                     nu = c(2,2,2)),
                   table = tableControl(npde = TRUE, censMethod = c("cdf")))

## save output to use later
saveRDS(run1.lst, "run1.rds")
run1.lst <- readRDS("run1.rds")
#
### plot GOF ###################################
pdf("gof_1",width = 6, height = 4)
basic.gof(run.dat = run1.lst)
dev.off()
#
#### VPC + GGPLOT #############################################
#### Make a VPC with nlmixr 
library(vpc)
library(ggpubr)
run1_vpc <- nlmixr::vpc(run1.lst, show = list(obs_dv = TRUE))
# extract the simulated data
sim_data <- run1_vpc$sim
# specify that time is the idv column or it gives an error message asking for idv
sim_data$idv <- sim_data$time
# then make the vpc with the vpc package:
vpc1 <- vpc::vpc(sim = sim_data,
                 obs = run1.lst,
                 show = list(obs_dv = TRUE),
                 lloq = 6.334,
                 facet = "rows") 
vpc2 <- vpc::vpc_cens(sim = sim_data,
                      obs = run1.lst,
                      lloq = 6.334,
                      facet = "rows") 
plot1 <- vpc1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "", y = "DV",
       title = paste("VPC")) + 
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(-10, 40), breaks = seq(0,40,10))
plot2 <- vpc2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + theme (axis.title.y = element_text(size = 10))+
  labs(x = "Time since symptom onset (d)", y = "Probability of < LOD",
       title = paste("")) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.5))
#
pdf("vpc_1",width = 6, height = 4)
ggarrange(plot1, plot2, nrow=2, ncol=1, heights = c(0.7, 0.3))
dev.off()
################################################################### end VPC + GGPLOT
#
run.rec[1, ] <- c("vl_dat_14days_post_symptom_under65", "rTCL model", round(run1.lst$OBJF),round(run1.lst$AIC), round(run1.lst$BIC) )
#
## save run record to use later
saveRDS(run.rec, "runrec.rds")
run.rec <- readRDS("runrec.rds")
#
## Export run record output to csv ######################
write.csv(run.rec, "run_record.csv",row.names = F, quote = F)
#
## Export run record output to table format csv ######################
write.table(run.rec, "run_record_2.csv", sep = "\t",  row.names = F, quote = F)
#
#### Calculations from final estimates in manuscript (Table S2). 
### Basic reproduction number (R0)
gamma <- 1
delta <- 0.80
R0 <- gamma/delta
#
### Duration of virus production (L)
L <- 1/delta
## end
