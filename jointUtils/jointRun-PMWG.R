
pmwg_jointRun <- function(experiments, epsilon = NULL, sharedPars = NULL){
  #This trick is a bit lame, but if I don't pass a df, the attributes will get lost in pmwg, so I must use a second column otherwise it 
  #is automatically converted to a vector by pmwg.
  #For now I've used an attribute approach, advantage is that there's no modifying of data structures, disadvantage is that complete datasets
  #are passed for every subject. I could check whether passing data specific to a subject is faster, or whether passing df as entries in a df could work. 
  df <- data.frame(subject = as.numeric(levels(unique(experiments[[1]]$preppedData$subject))), NAs = 0)
  pars <- character()
  priors <- numeric()
  modelName <- "Joint"
  startpoints <- NULL
  #This makes sure that the pars, startpoints and priors for each model in the joint model are still identifiable
  for(exp in experiments){
    # if(exp == "Rev"){
    #   exp$preppedData <- exp$preppedData[exp$preppedData$subject != 2,]
    # }
    attr(df, exp$modelName) <- exp$preppedData
    pars <- c(pars, paste(exp$modelName, "|", exp$parNames, sep = ""))
    priors <- c(priors, exp$priorMean)
    startpoints <- c(startpoints, exp$startpoints)
    modelName <- paste0(modelName, "_", exp$modelName)
  }
  for (par in names(sharedPars)){
    pars <- c(pars, par)
    priors <- c(priors, sharedPars[[par]][[1]])
    startpoints <- c(startpoints, sharedPars[[par]][[2]])
  }
  priors <- list(
    theta_mu_mean = priors,
    theta_mu_var = diag(rep(9, length(pars)))
  )
  
  # Create the Particle Metropolis within Gibbs sampler object
  sampler <- pmwgs(
    data = df,
    pars = pars,
    ll_func = jointLL,
    prior = priors
  )
  sampler = init(sampler, start_mu = startpoints)
  pmwg_runSampler(sampler, paste0("samples/Joint/", modelName, ".RData"), epsilon)
}



jointLL <- function(pars, data){
  #This makes sure that the likelihoods for every model are calculated separately based on prefixes and then added up. 
  nonSharedParsIdx <- grepl("[|]", names(pars))
  sharedPars <- pars[!nonSharedParsIdx]
  pars <- pars[nonSharedParsIdx]
  parPreFixs <- gsub("[|].*", "", names(pars))
  sum <- 0
  for (model in unique(parPreFixs)){
    currentPars <- pars[which(parPreFixs == model)]
    names(currentPars) <- gsub(".*[|]", "", names(currentPars))
    currentPars <- c(currentPars, sharedPars)
    modelData <- attr(data, model)[attr(data, model)[,"subject"] == data$subject,]
    if(nrow(modelData) > 0){
      sum <- sum + likelihood.RD(currentPars, modelData) 
    }
  }
  return(sum)
}