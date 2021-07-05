# Setup -------------------------------------------------------------------
rm(list=ls())
library(pmwg)
library(dplyr)
setwd("~/RL-PMWG/")

source("utils/utils-PMWG.R")
source("utils/DataPrep-PMWG.R")
source("jointUtils/DataPrepJoint-PMWG.R")
source("jointUtils/jointUtils-PMWG.R")
source("utils/pmwg_Run.R")
source("models/RW/RD.R")
source("models/RW/dists.R")
source("jointUtils/jointRun-PMWG.R")
#source("models/RW/dists_Timing.R")
source("jointUtils/plotSessionUtils.R")

library(dmcAdapt)

#Name: data name to load
#SubNCol: Where we'll look for the subject identifier
#accCol: Response correctness column in data, used to filter data
#accCriterion: Filter data based on accuracy per subject
#RL: Whether to employ RL learning rules in drift rate calculations (for now assumes an advantage framework)
#RTLimits: RTs above and below these will be filtered out
#Factors: Columns of the original data that should be included in the sampler, should also include the condition on which response correctness is determined if response-coded
#Match: a list with as name the parameter that should have matched entries in the sampler. If response coded first vector entry is the response, second vector entry is whether that response is correct. 
#Constants: Constants in the parameters that you want to add, e.g. to satisfy scaling constraints

#ParNames: The names of the parameters you want to include in the sampler. This has quite a specific format:
#Format explained: First parameter is the base parameter, then if included in match _1 means incorrect response, _2 correct (DDM convention).
#continued, then underscore for potential other factors and to which level of this factor it belongs to. Add additional factors with underscore. Separate factor and level with .
#If factor is the first entry of the match vector, it will assume that there should be separate values for the accumulators. e.g. b_response.1 means separate boundaries for the accumulators, this one for accumulator 1.

#llFunc, the likelihood function to use, for now only racing diffusion (Wald) included.
#priorMean: vector of lenght parnames with priors in order. Assumes normal priors and pars so transform in LL if necessary. 
#startPoints: Sensible startpoints to start sampling. (I was lazy)

modFunc <- function(a, b){(1 + a)*b}
plusFunc <- function(a,b){a+b}
simplexFunc <- function(b){1 -b}

experiments <- list("MSIT" = list(name = "MSIT", subNCol = "subjectNumber", accCol = "correct", accCriterion = 0.5, 
                                  "RL" = F, RTLimits = c(0.15, 2), factors = c("pos", "resp", "condit", "uniq"),
                                  constants = list(s = 1, A = 0), match = list(v = c("resp", "uniq"), B = c("resp", "pos")),
                                  parNames = c("B_1", "B_2", "t0", "v_1_condit.1", "v_1_condit.2", "v_1_condit.3", "v_1_condit.4",
                                               "v_2_condit.1", "v_2_condit.2", "v_2_condit.3", "v_2_condit.4"), 
                                  llFunc = likelihood.RD, priorMean = rep(0.1, 11), startpoints = NULL,
                                  modelName = "MSIT4v2B_t(t0)"),
                    
                    "RB" = list(name = "RB", subNCol = "subjectNumber", accCol = "correct", accCriterion = 0.55,
                                "RL" = F, RTLimits = c(0.15, 2.2), factors = c("type", "trans", "correctR"),
                                constants = list(s = 1, A= 0), match = list(v = c("resp", "correctR")),
                                parNames = c("v_1_type.comparison_trans.repeat", "v_1_type.reference_trans.repeat",  
                                             "v_1_type.comparison_trans.switch", "v_1_type.reference_trans.switch", 
                                             "v_2_type.comparison_trans.repeat", "v_2_type.reference_trans.repeat",  
                                             "v_2_type.comparison_trans.switch", "v_2_type.reference_trans.switch", 
                                             "t0", "B"), 
                                llFunc = likelihood.RD, priorMean = rep(0.1, 10), startpoints = NULL,
                                modelName = "RB8v_t(t0)"),
                    
                    "rlsat" = list(name = "rlsat", subNCol = "subjectNumber", accCol = "correct", accCriterion = 1.55, #the RL experiments have accuracy stored as 1 wrong, 2 correct. 
                                   "RL" = T, RTLimits = c(0.15, 2), factors = c("cue"), match = list(v = "correct"), 
                                   constants = list(s = 1, A= 0, SR = pnorm(-10)), 
                                   parNames = c("t0", "V0_cue.SPD", "V0_cue.ACC", "B_cue.SPD", "B_cue.ACC", "wS", "wD", "aV"), 
                                   llFunc = likelihood.RD, priorMean = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, -1.2), 
                                   startpoints = NULL, modelName = "rlsat2V02B_t(t0)"),
                    
                    "Rev" = list(name = "Rev", subNCol = "subjectNumber", accCol = "correctTotal", accCriterion = 1.55, 
                                 "RL" = T, RTLimits = c(0.15, 2), factors = NULL, match = list(v = "correct"),
                                 constants = list(s = 1, A= 0, SR = pnorm(-10)), 
                                 parNames = c("t0", "V0", "B", "wS", "wD", "aV"), 
                                 llFunc = likelihood.RD, priorMean = c(0.1, 0.1, 0.1, 0.1, 0.1, -1.2), 
                                 startpoints = NULL, modelName = "RevStd_t(t0)")
)
experiments <- prepDataJoint(experiments, splitSess = F, hardCutAcc = F, hardCutN = F) #If splitsess, we can run a joint model comparing session one with session two for the current task. 
# pmwg_jointRun(experiments[c(2,6)], epsilon = 0.5)
# pmwg_jointRun(experiments[c(3,7)], epsilon = 0.3)
# pmwg_jointRun(experiments[c(4,8)], epsilon = 0.35)
pmwg_run(experiments[1], epsilon = 0.4)
pmwg_run(experiments[2], epsilon = 0.5)
pmwg_run(experiments[3], epsilon = 0.45)
pmwg_run(experiments[4], epsilon = 0.45)
pmwg_jointRun(experiments, epsilon = 0.23)
# pmwg_post(experiments[c(1, 2)], PP = T) 

