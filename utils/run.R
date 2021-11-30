
runSampler <- function(sampler, file = NULL, epsilon = 0.5, pstar = .5, experiments = NULL, n_cores){
  burned <- run_stage(sampler,
                      stage = "burn",
                      iter = 1000,
                      useC = F,
                      particles = 100,
                      n_cores = n_cores,
                      pstar = pstar
  )
  save(burned, file = file)
  burned <- run_stage(burned,
                      stage = "burn",
                      iter = 1000,
                      useC = F,
                      particles = 100,
                      n_cores = n_cores,
                      pstar = pstar
  )
  save(burned, file = file)
  burned <- run_stage(burned,
                      stage = "burn",
                      iter = 1500,
                      useC = T,
                      particles = 100,
                      n_cores = n_cores,
                      pstar = pstar
  )
  save(burned, file = file)
  adapted <- run_stage(adapted,
                       stage = "adapt",
                       iter = 10000, #Set up quite high, should terminate early anyway if not, likely a problem
                       particles = 100,
                       useC = T,
                       n_cores = n_cores,
                       pstar = pstar
  )
  save(adapted, file = file)
  sampled <- run_stage(sampled, 
                       stage = "sample",
                       iter = 5000, 
                       particles = 100,
                       useC = F,
                       n_cores = n_cores,
                       pstar = pstar
  )
  save(sampled, file = file)
  sampled$epsilon = epsilon
  sampled$experiments = experiments
  save(sampled, file = file)
}

pmwgRun <- function(experiments, epsilon = NULL, n_cores){
  for(exp in experiments){
    
    priors <- list(
      theta_mu_mean = exp$priorMean,
      theta_mu_var = diag(rep(1, length(exp$parNames)))
    )
    
    # Create the Particle Metropolis within Gibbs sampler object
    sampler <- pmwgs(
      data = exp$preppedData,
      pars = exp$parNames,
      ll_func = exp$llFunc,
      prior = priors
    )
    sampler = init(sampler, start_mu = exp$startPoints, n_cores = n_cores)
    runSampler(sampler, file = paste0("samples/", exp$modelName, ".RData"), epsilon, pstar = .5, experiments = experiments, n_cores)
  }
}


