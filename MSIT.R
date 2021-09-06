# Setup -------------------------------------------------------------------
rm(list=ls())
library(pmwg)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
setwd("~/RL-PMWG/")

source("utils/dataPrep.R")
source("utils/diagnostics.R")
source("utils/postPredict.R")
source("utils/run.R")

source("models/RW/RD.R")
source("models/RW/dists_ARD.R")
source("models/RW/dists.R")

library(dmcAdapt)

experiments <- list("MSIT" = list(name = "MSIT", rtCol = 'rt', subNCol = "subjectNumber", accCol = "correct", accCriterion = 0.5,
                                  "RL" = F, RTLimits = c(0.15, 2), factors = c("pos", "resp", "condit", "uniq", "flank", "nstims", "stimuli"), 
                                  transFunc = list(func = driftsMSIT3, earlyTransform = F),
                                  constants = list(s = 1, A = 0, Qs = 0, QMatch = 1), match = list(v = c("resp")),
                                  parNames = c("t0", "QFlank", "QSimon", "V0",  "wD", "B", "wS"),
                                  llFunc = likelihood.ARD, priorMean = c(-0.5, rep(0.1, 6)), startpoints = NULL,
                                  modelName = "MSITprocsAdv")
)
experiments <- prepDataJoint(experiments, splitSess = F, hardCutAcc = T, hardCutN = T) #If splitsess, we can run a joint model comparing session one with session two for the current task. 
pmwgRun(experiments, epsilon = 0.47)
# pmwg_post(experiments[2], PP = T) 
