source("jointUtils/jointPostPredict.R")
source("jointUtils/jointDiagnostics.R")

pmwgJointRun <- function(experiments, epsilon = NULL, sharedPars = NULL){
  #Create a dataframe, with a row for each subject
  df <- data.frame(subject = as.numeric(levels(unique(experiments[[1]]$preppedData$subject))))
  pars <- character()
  priors <- numeric()
  llFuncs <- list()
  modelName <- "Joint"
  startpoints <- NULL
  #This makes sure that the pars, startpoints and priors for each model in the joint model are still identifiable
  for(exp in experiments){
    #Add columns to the df that hold 
    df[exp$modelName] <- I(list(split(exp$preppedData, f = exp$preppedData$subject)))
    llFuncs[[exp$modelName]] <- exp$llFunc
    pars <- c(pars, paste(exp$modelName, "|", exp$parNames, sep = ""))
    priors <- c(priors, exp$priorMean)
    startpoints <- c(startpoints, exp$startpoints)
    modelName <- paste0(modelName, "_", exp$modelName)
  }
  
  attr(df, 'llFuncs') <- llFuncs
  modelName <- paste0(modelName, "_eps", epsilon)
  for (par in names(sharedPars)){
    pars <- c(pars, par)
    priors <- c(priors, sharedPars[[par]]$priors)
    startpoints <- c(startpoints, sharedPars[[par]]$startpoints)
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
    prior = priors
  )
  sampler = init(sampler, start_mu = startpoints)
  path <- paste0("samples/Joint/", modelName, ".RData")
  if(!is.null(sharedPars)) path <- paste0("samples/Joint/sharedPars/", modelName, ".RData")
  runSampler(sampler, path, epsilon, experiments = experiments)
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
    if(nrow(modelData) > 0){
      totalSum <- totalSum + llFuncs[[model]](currentPars, modelData)
    }
  }
  return(totalSum)
}
