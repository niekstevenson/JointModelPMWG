
runSampler <- function(sampler, file = NULL, epsilon = 0.5, experiments = NULL){
  print(paste0("epsilon: ", epsilon))
  burned <- run_stage(sampler, 
                      stage = "burn",
                      iter = 5000,
                      particles = 100,
                      n_cores = 15,
                      epsilon = epsilon
  )
  save(burned, file = file)
  adapted <- run_stage(burned, 
                       stage = "adapt", 
                       iter = 10000, #Set up quite high, should terminate early anyway if not, likely a problem
                       particles = 100,
                       n_cores = 15
  )
  save(adapted, file = file)
  sampled <- run_stage(adapted, 
                       stage = "sample",
                       iter = 2000, 
                       particles = 100,
                       n_cores = 15,
                       epsilon = epsilon
  )
  save(sampled, file = file)
  sampled$epsilon = epsilon
  sampled$experiments = experiments
  save(sampled, file = file)
}

pmwgRun <- function(experiments, epsilon = NULL){
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
    sampler = init(sampler, start_mu = exp$startPoints)
    runSampler(sampler, file = paste0("samples/", exp$modelName, ".RData"), epsilon, experiments = experiments)
  }
}


