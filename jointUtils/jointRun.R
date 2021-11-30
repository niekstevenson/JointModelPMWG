source("jointUtils/jointPostPredict.R")
source("jointUtils/jointDiagnostics.R")

pmwgJointRun <- function(experiments, epsilon = NULL, pstar = NULL, n_factors = NULL, sharedPars = NULL, n_cores = 16){
  #Create a dataframe, with a row for each subject
  subjects <- sapply(experiments, FUN = function(x) return(x$preppedData[,'subject'][[1]])) %>% levels() %>% as.numeric()
  df <- data.frame(subject = subjects)
  pars <- character()
  priors <- numeric()
  llFuncs <- list()
  modelName <- "Joint"
  startpoints <- NULL
  for(exp in experiments){
    currentData <- data.frame(subject = as.numeric(levels(unique(exp$preppedData$subject))))
    currentData[exp$name] <- I(list(split(exp$preppedData, f = exp$preppedData$subject)))
    llFuncs[[exp$name]] <- exp$llFunc
    pars <- c(pars, paste(exp$name, "|", exp$parNames, sep = ""))
    priors <- c(priors, exp$priorMean)
    startpoints <- c(startpoints, exp$startpoints)
    modelName <- paste0(modelName, "_", exp$modelName)
    df <- base::merge(df, currentData, all = T)
  }
  
  attr(df, 'llFuncs') <- llFuncs
  modelName <- paste0(modelName, "_eps", epsilon)
  for (par in names(sharedPars)){
    pars <- c(par, pars)
    priors <- c(sharedPars[[par]]$priors, priors)
    startpoints <- c(sharedPars[[par]]$startpoints, startpoints)
  }
  priors <- list(
    theta_mu_mean = priors,
    theta_mu_var = diag(rep(1, length(pars)))
  )
  # Create the Particle Metropolis within Gibbs sampler object
  sampler <- pmwgs(
    data = df,
    pars = pars,
    ll_func = jointLL,
    n_factors = n_factors
  )
  sampler = init(sampler, start_mu = startpoints, useC = F, n_cores = n_cores)
  if(!is.null(sharedPars)) modelName <- paste0(modelName, "_sharedt0")
  path <- paste0("samples/Joint/factor/", n_factors, "/", modelName, ".RData")
  runSampler(sampler, path, epsilon, experiments = experiments, pstar, n_cores)
}

jointLL <- function(pars, data){
  #This makes sure that the likelihoods for every model are calculated separately based on prefixes and then added up.
  nonSharedParsIdx <- grepl("[|]", names(pars))
  sharedPars <- pars[!nonSharedParsIdx]
  pars <- pars[nonSharedParsIdx]
  parPreFixs <- gsub("[|].*", "", names(pars))
  llFuncs <- attr(data, 'llFuncs')
  totalSum <- 0
  for (model in unique(parPreFixs)){
    currentPars <- pars[which(parPreFixs == model)]
    names(currentPars) <- gsub(".*[|]", "", names(currentPars))
    currentPars <- c(currentPars, sharedPars)
    modelData <- data[[model]][[1]]
    if(!is.na(modelData)){
      totalSum <- totalSum + llFuncs[[model]](currentPars, modelData)
    }
  }
  return(totalSum)
}
