rm(list=ls())
###### set wd  ######
mod.dir <- "/home/akosua/Documents/Model_evaluation/SARS_CoV-2_viral_dynamic_models_comparison/Gastine_et_al_2021/Gastine_SI/"
setwd(mod.dir)
library(ggplot2)
library(scales)
library(deSolve)
library(nlmixr)
##### read data ######
vl_dat <- read.csv("Gastine_data.csv", h = TRUE)
#
# remove EVID=2 rows which were added because NONMEM carries backwards
# covariates whereas nlmixr carries forward.
vl_dat <- vl_dat[vl_dat$EVID!=2, ]
# standardise column names id, time and dv for nlmixr (we will wok with log(dv))
# CENSset to BLOD, LIMIT set to limit of detection
colnames(vl_dat) <- c("id",  "time", "sample_site_num", "sample_site_group_num",
                      "vl_quality", "assay", "paper", "dv_ct",                 
                      "dv_vl_log10_cpml", "dv", "bloq_1", "CENS",                
                      "loq_value_log10_cpml", "lod_value_log10_cpml", "loq_value_log_cpml", "LIMIT",    
                      "age_y", "male", "fever", "icu",                   
                      "days_on_icu", "vent", "days_ventilated", "died",                  
                      "time_to_death", "disease_status", "drug_quality", "drug_code",
                      "drug_lpvr", "drug_ifn", "drug_ifn_alpha", "drug_ifn_beta",
                      "drug_cqhcq", "drug_remd", "drug_azit", "drug_riba",
                      "drug_umif", "drug_thym", "sample_site_group_num2", "urt",
                      "lrt", "nose", "asymp", "mild", 
                      "moderate", "severe", "EVID", "L2",                    
                      "TYPE")
# LIMIT data item:
# from the nlmixr guide: https://nlmixrdevelopment.github.io/RxODE/articles/RxODE-datasets.html
# CENS = 1 the value is censored between (LIMIT, DV)
# make sure dv is set to LOD for all censored items
vl_dat$dv[vl_dat$CENS == 1] <-  vl_dat$LIMIT[vl_dat$CENS == 1]
# 2. set LIMIT close to 0 for CENS = 1 
vl_dat$LIMIT[vl_dat$CENS == 1] <- log(0.001)
#
##### plot data #######
pl.1 <- ggplot(aes(x = time, y = dv_vl_log10_cpml, colour = id), data = vl_dat) + 
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_gradientn(colours = rainbow(3)) +
  labs(x = "Time since symptom onset (d)", y = "log10 copies/mL",
       title = "SARS-CoV-2 viral loads") +
  geom_point(cex = 2) + 
  geom_line() +
  geom_smooth() +
  facet_wrap(~ sample_site_group_num)
pl.1
#
# make simplified datasets for modelling, limited to 14 days post symptom onset 
vl_dat_14days_post_symptom <- vl_dat[vl_dat$sample_site_group_num == 2 & vl_dat$time <= 14 &
                                      vl_dat$time >= 0 & vl_dat$drug_code == 0, ]
######### write output file for modelling
write.csv(vl_dat_14days_post_symptom , "vl_dat_14days_post_symptom.csv", row.names = F, quote = F)
#
##### read in clean data  ######
vl_dat_clean <- read.csv("vl_dat_14days_post_symptom.csv", h = TRUE)
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
################## end of GOF function
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
    tv0    <- 13    # Viral load at symptom onset
    
    eta.delta + eta.v0  ~ c(1,
                            0.01, 1)
    
    add.err <- 1     # residual variability
    
  })
  model({
    
    delta <- exp(tdelta + eta.delta)  # individual value of delta
    v0    <- exp(tv0 + eta.v0)           # individual value of v0
    
    
    A_v(0) = v0
    
    d/dt(A_v) =  - delta * A_v #copies/ml
    
    vt = log(A_v) #copies/ml to lncopies/ml
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
  scale_y_continuous(limits = c(-10, 30), breaks = seq(0,25,10))
#
plot2 <- vpc2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) + theme (axis.title.y = element_text(size = 10))+
  labs(x = "Time since symptom onset (d)", y = "Probability of <LOD",
       title = paste("")) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,5)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.5))
#
pdf("vpc",width = 6, height = 4)
ggarrange(plot1, plot2, nrow=2, ncol=1, heights = c(0.7, 0.3))
dev.off()
################################################################### end VPC + GGPLOT
run.rec[1, ] <- c("vl_dat_14days_post_symptom", "SI model", round(run1.lst$OBJF),round(run1.lst$AIC), round(run1.lst$BIC) )
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
### Duration of virus production (L)
delta <- 0.56
L <- 1/delta # 1.79
## end




