rm(list=ls())
###### set working directory  ######
mod.dir <- "/home/akosua/Documents/Model_evaluation/SARS_CoV-2_viral_dynamic_models_comparison/Neant_et_al_2021/Neant_SI/"
setwd(mod.dir)
library(ggplot2)
library(scales)
library(deSolve)
library(nlmixr)
##### read data ######
vl_dat <- read.csv("Neant_data.csv", h = TRUE)
#
str(vl_dat)
#### convert time from character class to numeric class
vl_dat$time_monolix = as.numeric(as.character(vl_dat$time_monolix))
str(vl_dat)
#####remove type = 2 and NA rows
vl_dat <- vl_dat[vl_dat$type!=2, ]
vl_dat <- vl_dat[!is.na(vl_dat$ID), ]
#
#### rename column y to dv_log10_cpml
names(vl_dat)[names(vl_dat) == "y"] <- "dv_log10_cpml"
##### create a column for time_since_symptom_onset and rename as time
vl_dat$time_since_symptom_onset <- vl_dat$time_monolix - 14 # time_monolix = time_since_symptom + 14 days
names(vl_dat)[names(vl_dat) == "time_since_symptom_onset"] <- "time"
#
######### plot data: Figure 1A from Neant et al., doi:  10.1073/pnas.2017962118 #######
pl.1 <- ggplot() + 
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Time since symptom onset", y = "log10 copies/mL",
       title = "SARS-CoV-2 viral loads") +
  geom_point(aes(x = time, y = dv_log10_cpml, group = ID) , data = vl_dat, cex = 2) + 
  geom_line(aes(x = time, y = dv_log10_cpml, group = ID), data = vl_dat) + 
  scale_x_continuous(limits = c(0, 63), breaks = c(seq(0, 63, 7))) + scale_y_continuous(limits = c(2, 12), breaks = c(seq(0, 12, 2)))
geom_smooth() 
facet_wrap(~ NULL)
pl.1
#
#### create a column for dv_vl_cpml (i.e., viral load in copies per ml)
vl_dat$dv_vl_cpml <- 10^(vl_dat$dv_log10_cpml)
##### create a column for logcpml (i.e., natural log (ln)  of viral load in copies per ml)
vl_dat$dv_vl_logcpml <- log(vl_dat$dv_vl_cpml)
###rename dv_vl_log_cpml as dv in a new column
vl_dat$dv <- vl_dat$dv_vl_logcpml
#
########## LIMIT data item:
#######  from the nlmixr guide: https://nlmixrdevelopment.github.io/RxODE/articles/RxODE-datasets.html
######## CENS = 1 the value is censored between (LIMIT, DV)
######## set LIMIT close to 0 for CENS = 1 
##### make sure dv is set to LOD for all censored items (Note: in Neant data set dv is already set to LOD for censores items)
########### Create nlmixr column called LIMIT and set to 0. 
vl_dat$LIMIT <- 0
vl_dat$LIMIT[vl_dat$CENS == 1] <- log(0.001)
#
##### Limit data set to <= 14 d time since symptom onset, under 65 years
vl_dat_days14_under65 <- vl_dat[vl_dat$time <= 14 & vl_dat$time >= 0 & vl_dat$age_cat_cov==1, ]
######## write output file for limited data set
write.csv(vl_dat_days14_under65 , "vl_dat_14days_post_symptom_under65.csv", row.names = F, quote = F) 
#
##### read in clean data  ######
vl_dat_clean <- read.csv("vl_dat_14days_post_symptom_under65.csv", h = TRUE)
#
########### GOF function plots ####################################################
#
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
################## end of function ####################
#
###### SI model ######
# Parameter estimates from Gastine et al: doi: 10.1002/cpt.2223
### start the run record:
run.rec <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(run.rec) <- c("Dataset", "Model", "OFV", "AIC", "BIC")
#
SI.model <- function() {
  ini({
    tdelta <- -0.5  # Death rate of infected cells
    tv0    <- 13    # Viral load at symptom onset)
    
    eta.delta + eta.v0  ~ c(1,
                            0.01, 1)
                                        

    add.err <- 1     # residual variability
    
  })
  model({
    delta <- exp(tdelta + eta.delta)  # individual value of delta
    v0    <- exp(tv0 + eta.v0)           # individual value of v0
   
    
    A_f(0) = 1
    A_v(0) = v0
    
  
    d/dt(A_v) = - delta * A_v 
    
    vt = log(A_v)
    vt ~ add(add.err)       # define error model
  })
}
# Check the model
nlmixr(SI.model)
## Fit the model 
run1.lst <- nlmixr(SI.model, vl_dat_clean , est = "saem",
                   control = saemControl(
                     seed = 99,
                     nu = c(2,2,2)),
                   table = tableControl(npde = TRUE, censMethod = c("cdf")))
#
## save output to use later
saveRDS(run1.lst, "run1.rds")
run1.lst <- readRDS("run1.rds")
#
### plot GOF ###################################
pdf("gof",width = 6, height = 4)
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

pdf("vpc",width = 6, height = 4)
ggarrange(plot1, plot2, nrow=2, ncol=1, heights = c(0.7, 0.3))
dev.off()
################################################################### end VPC + GGPLOT
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
#### Calculations from final estimates in manuscript (Table S3). 
### Duration of virus production (L)
delta <- 0.85 
L <- 1/delta
## end


