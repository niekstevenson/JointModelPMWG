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

source("pmwg/sampling_factors.R")
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
#Match: if accuracy coded, only provide accuracy column. 
#Constants: Constants in the parameters that you want to add, e.g. to satisfy scaling constraints
#ParNames: The names of the parameters you want to include in the sampler. This has quite a specific format:
#Format explained: First parameter is the 'base' parameter, then if included in match _1 means incorrect response, _2 correct (DDM convention).
#continued, For all parameters add underscore for potential other factors and to which level of this factor it belongs to. Add additional factors with additional underscore. Separate factor and level with .
#If factor is the first entry of the match vector, it will assume that there should be separate values for the accumulators. e.g. b_response.1 means separate boundaries for the accumulators, this one for accumulator 1.
#llFunc, the likelihood function to use, for now only racing diffusion (Wald) included.
#priorMean: vector of lenght parnames with priors in order. Assumes normal priors and pars so transform in LL if necessary.
#startPoints: Sensible startpoints to start sampling. (I was lazy)

modFunc <- function(a, b){(1 + a)*b}
plusFunc <- function(a,b){a+b}
simplexFunc <- function(b){1 -b}

experiments <- list(
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
               parNames = c("V0", "B", "wD", "aV", "wS"),
               llFunc = likelihood.ARD,
               priorMean = rep(0.1, 5),
               startpoints = NULL,
               modelName = "RevStd"),
  "MSIT" = list(name = "MSIT",
                subNCol = "subjectNumber",
                accCol = "correct",
                rtCol = 'rt',
                accCriterion = 0.5,
                "RL" = F,
                RTLimits = c(0.15, 2),
                factors = c("pos", "resp", "condit", "uniq", "flank", "stimuli"),
                transFunc = list(func = driftsMSIT, earlyTransform = F),
                constants = list(s = 1, A = 0, vPos.3 = 0),
                match = list(vPos = c("resp", "uniq")),
                parNames = c("vFlank", "vSimon", "vPos.1", "vPos.2", "vMatch", "v", "B"),
                llFunc = likelihood.RD,
                priorMean =  rep(0.1, 7),
                startpoints = NULL,
                modelName = "MSITprocsBSimon"),
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
                           "diff_correctRName.same", "diff_correctRName.diff", "B_trans.repeat", "B_trans.switch"),
              llFunc = likelihood.RD,
              priorMean = rep(0.1, 8),
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
             parNames = c("V0_cue.SPD", "V0_cue.ACC", "B_cue.SPD", "B_cue.ACC", "wS", "wD", "aV"),
             llFunc = likelihood.ARD,
             priorMean = rep(0.1, 7),
             startpoints = NULL,
             modelName = "rlsat2V02B")
)
experiments <- prepDataJoint(experiments, splitSess = F, hardCutAcc = T, hardCutN = T, n_sessions = 2) #If splitsess, we can run a joint model comparing session one with session two for the current task.

sharedPars <- list(t0 = list(priors = -1, startpoints = NULL))
set.seed(1234)
pmwgJointRun(experiments, pstar = .7, n_factors = 5, n_cores = 20, sharedPars = sharedPars)
