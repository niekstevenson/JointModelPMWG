# Setup -------------------------------------------------------------------
rm(list=ls())
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
source("pmwg/sampling.R")


exp <- list("MSIT" = list(name = "MSIT", rtCol = 'rt', subNCol = "subjectNumber", accCol = "correct", accCriterion = 0.5,
                                  "RL" = F, RTLimits = c(0.15, 2), factors = c("pos", "resp", "condit", "uniq", "flank", "nstims", "stimuli"),
                                  transFunc = list(func = driftsMSIT, earlyTransform = F),
                                  constants = list(s = 1, A = 0, vPos.3 = 0), match = list(v = c("resp", "uniq")),
                                  parNames = c("t0", "vFlank", "vSimon", "vPos.1", "vPos.2", "vMatch", "v", "B"),
                                  llFunc = likelihood.RD, priorMean = c(-0.5, rep(0.1, 7)), startpoints = NULL,
                                  modelName = "MSIT_procs")
)
exp <- prepDataJoint(exp, splitSess = F, hardCutAcc = T, hardCutN = T) #If splitsess, we can run a joint model comparing session one with session two for the current task.
debug(runSampler)
pmwgRun(exp, n_cores = 12)