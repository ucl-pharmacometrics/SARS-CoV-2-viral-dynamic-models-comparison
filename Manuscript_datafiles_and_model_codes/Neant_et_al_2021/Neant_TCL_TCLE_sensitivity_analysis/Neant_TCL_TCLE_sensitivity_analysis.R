############set working directory
setwd("/home/akosua/Documents/Model_evaluation/SARS_CoV-2_viral_dynamic_models_comparison/Neant_et_al_2021/Neant_TCL_TCLE_sensitivity_analysis/")
library(ggplot2)
library(scales)
library(deSolve)
##### Read in data ########
dat <- read.csv("Neant_TCL_TCLE_sensitivity_analysis_BIC.csv", h = TRUE)

pl.1 <- ggplot() + 
  theme_bw() +
  labs(x = "Time since infection (d)", y = "BIC",
       title = "Sensitivity Analysis") +
  geom_point(aes(x = Day, y = BIC, group = Model, colour = Model) , data = dat, cex = 2) + 
  geom_line(aes(x = Day, y = BIC, group = Model, colour = Model), data = dat) + 
  scale_x_continuous(limits = c(0, 14), breaks = c(seq(0, 14, 1)))
geom_smooth() 
facet_wrap(~ NULL)
pdf("Neant_TCL_TCLE_Sensitivity_Analysis",width = 6, height = 4)
pl.1
dev.off() 




