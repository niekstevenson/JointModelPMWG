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
library(dmcAdapt)

exp <- list("rlsat" = list(name = "rlsat",
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
                             priorMean = rep(0.1, 8),
                             startpoints = NULL,
                             modelName = "rlsat2V02B")
)
exp <- prepDataJoint(exp, splitSess = F, hardCutAcc = T, hardCutN = T) #If splitsess, we can run a joint model comparing session one with session two for the current task.
n_cores <- 25


file <- "samples/RLSAT_try.RData"

priors <- list(
  theta_mu_mean = exp$rlsat$priorMean,
  theta_mu_var = diag(rep(1, length(exp$rlsat$parNames)))
)

# Create the Particle Metropolis within Gibbs sampler object
sampler <- pmwgs(
  data = exp$rlsat$preppedData,
  pars = exp$rlsat$parNames,
  ll_func = exp$rlsat$llFunc,
  prior = priors
)
sampler = init(sampler, start_mu = exp$rlsat$startPoints, n_cores = n_cores, useC = F)
save(sampler, file = file)
burned <- run_stage(sampler,
                    stage = "burn",
                    iter = 2500,
                    particles = 100,
                    n_cores = n_cores,
                    pstar = .6,
                    useC = F
)
save(burned, file = file)
adapted <- run_stage(burned,
                     stage = "adapt",
                     iter = 10000, #Set up quite high, should terminate early anyway if not, likely a problem
                     particles = 100,
                     n_cores = 5,
                     pstar = .6,
                     useC = F
)
save(adapted, file = file)
sampled <- run_stage(adapted, 
                     stage = "sample",
                     iter = 2000, 
                     particles = 100,
                     n_cores = 1,
                     pstar = .6,
                     useC = F
)
save(sampled, file = file)