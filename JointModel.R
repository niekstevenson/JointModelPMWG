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
source("jointUtils/jointRun.R")
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
experiments <- list("MSIT" = list(name = "MSIT",
                                  subNCol = "subjectNumber",
                                  accCol = "correct",
                                  rtCol = 'rt',
                                  accCriterion = 0.5,
                                  "RL" = F,
                                  RTLimits = c(0.15, 2),
                                  factors = c("pos", "resp", "condit", "uniq", "flank", "stimuli"),
                                  transFunc = list(func = driftsMSIT, earlyTransform = F),
                                  constants = list(s = 1, A = 0),
                                  match = list(vPos = c("resp", "uniq")),
                                  parNames = c("t0", "v", "vFlank", "vSimon","vPos.1", "vPos.2", "vPos.3", "B", "B2"),
                                  llFunc = likelihood.RD,
                                  priorMean =  c(-0.8, 0.5, 1.15, 1.15, 2.8, 2.56, 2.56, 1.6, 1.55), 
                                  startpoints = NULL,
                                  modelName = "MSITprocsB2"),
                    "RB" = list(name = "RB",
                                subNCol = "subjectNumber",
                                accCol = "correct",
                                rtCol = 'rt',
                                accCriterion = 0.55,
                                "RL" = F,
                                RTLimits = c(0.15, 2.2),
                                factors = c("type", "trans", "correctR", "correctRName", "typTrans"),
                                transFunc = list(func = driftsRePar, earlyTransform = T),
                                constants = list(s = 1, A= 0),
                                match = list(v = c("resp", "correctR")),
                                parNames = c("V0_type.comparison_typTrans.repeat", "V0_type.reference_typTrans.repeat",
                                             "V0_type.comparison_typTrans.switch", "V0_type.reference_typTrans.switch",
                                             "diff_correctRName.same", "diff_correctRName.diff", "B_trans.repeat", "B_trans.switch", "t0"),
                                llFunc = likelihood.RD,
                                priorMean = c(1.88, 1.55, 1.58, 1.47, 2.8, 2.5, 1.26, 1.48, -0.9),
                                startpoints = NULL,
                                modelName = "RB4V02Diff2B"),
                    "rlsat" = list(name = "rlsat",
                                   subNCol = "subjectNumber",
                                   accCol = "correct",
                                   accCriterion = 1.55, #the RL experiments have accuracy stored as 1 wrong, 2 correct.
                                   rtCol = 'rt',
                                   "RL" = T,
                                   RTLimits = c(0.15, 2),
                                   factors = c("cue"),
                                   match = list(v = "correct"),
                                   constants = list(s = 1, A= 0, Qs = pnorm(-10)),
                                   parNames = c("t0", "V0_cue.SPD", "V0_cue.ACC", "B_cue.SPD", "B_cue.ACC", "wS", "wD", "aV"),
                                   llFunc = likelihood.ARD,
                                   priorMean = c(-1.8, 3.6, 3.2, 2.14, 2.27, 0.5, 1.63, -1.25),
                                   startpoints = NULL,
                                   modelName = "rlsat2V02B"),
                    "Rev" = list(name = "Rev",
                                 subNCol = "subjectNumber",
                                 accCol = "correctTotal",
                                 rtCol = 'rt',
                                 accCriterion = 1.55,
                                 "RL" = T,
                                 RTLimits = c(0.15, 2),
                                 factors = NULL,
                                 match = list(v = "correct"),
                                 constants = list(s = 1, A= 0, Qs = pnorm(-10)),
                                 parNames = c("t0", "V0", "B", "wS", "wD", "aV"),
                                 llFunc = likelihood.ARD,
                                 priorMean = c(-1.32, 2, 1.56, 0.54, 1.47, -0.93),
                                 startpoints = NULL,
                                 modelName = "RevStd")
)
experiments <- prepDataJoint(experiments, splitSess = F, hardCutAcc = T, hardCutN = T) #If splitsess, we can run a joint model comparing session one with session two for the current task.
# pmwg_jointRun(experiments[c(2,6)], epsilon = 0.5)
# pmwg_jointRun(experiments[c(3,7)], epsilon = 0.3)
# pmwg_jointRun(experiments[c(4,8)], epsilon = 0.35)
# pmwgRun(experiments[1], epsilon = 0.45)
# pmwg_run(experiments[2], epsilon = 0.55)
# pmwg_run(experiments[3], epsilon = 0.45)
# pmwg_run(experiments[4], epsilon = 0.45)
pmwgJointRun(experiments, epsilon = 0.205)
# pmwg_post(experiments[2], PP = T)