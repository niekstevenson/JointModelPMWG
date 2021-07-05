
pmwg_runSampler <- function(sampler, file = NULL, epsilon = 0.5){
  print(paste0("epsilon: ", epsilon))
  burned <- run_stage(sampler, 
                      stage = "burn",
                      iter = 5000,
                      particles = 100,
                      n_cores = 15,
                      epsilon = epsilon
  )
  adapted <- run_stage(burned, 
                       stage = "adapt", 
                       iter = 10000, #Set up quite high, should terminate early anyway if not, likely a problem
                       particles = 100,
                       n_cores = 15,
                       epsilon = epsilon
  )

  sampled <- run_stage(adapted, 
                       stage = "sample",
                       iter = 2000, 
                       particles = 100,
                       n_cores = 15,
                       epsilon = epsilon
  )
  
  save(sampled, file = file)
}

pmwg_run <- function(experiments, epsilon = NULL){
  for(exp in experiments){
    
    priors <- list(
      theta_mu_mean = exp$priorMean,
      theta_mu_var = diag(rep(9, length(exp$parNames)))
    )
    
    # Create the Particle Metropolis within Gibbs sampler object
    sampler <- pmwgs(
      data = exp$preppedData,
      pars = exp$parNames,
      ll_func = exp$llFunc,
      prior = priors
    )
    sampler = init(sampler, start_mu = exp$startPoints)
    pmwg_runSampler(sampler, file = paste0("samples/", exp$modelName, ".RData"), epsilon)
  }
}


pmwg_post<- function(experiments, PP = T){
  #Written so that automatically plots are made after running, doesn't work for predictive plots RL yet
  for (exp in experiments){
    path <- file.path("figures", exp$modelName)
    dir.create(path)
    load(paste0("samples/", exp$modelName, ".RData"))
    pmwg_IC(sampled, preppedData = exp$preppedData)
    pmwg_chainPlots(sampled, subjectParPlot = T, parameterPlot = T, subjectLLPlot = T, PDF = T, path = path)
    pmwg_parIntervals(sampled, PDF = T, path = path)
    pmwg_parHist(sampled, PDF = T, path = path)
    
    # pmwg_particlesPerSub(sampled)
    # pmwg_parValues(sampled)
    # pmwg_covariance(sampled)
    if(PP){
      dat <- exp$data
      factors = exp$factors
      tmp <- pmwg_Simulate(samples = sampled, dat = dat, ll_func = exp$llFunc, RL = exp$RL, combns = factors, n = 25, match = exp$match)
      ppInfo <- getDescriptives(tmp$pp, combns = factors)
      dataInfo <- getDescriptives(tmp$data, combns = factors)
      pmwg_ppPlot(ppInfo, dataInfo, factors, PDF = T, path = path)
    }

  }
}

