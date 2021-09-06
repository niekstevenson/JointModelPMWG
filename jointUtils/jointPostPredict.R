jointPostCheck <- function(samples, experiments, n = 25, sharedPars = NULL){
  parPreFixs <- gsub("[|].*", "", samples$par_names)
  i <- 0
  jointModelName <- paste(lapply(experiments, function(x){return(x$modelName)}), collapse = "_")
  dir.create(file.path("./figures", "joint", jointModelName), recursive = T)
  jointChainPlots(samples, PDF = T, path = file.path("./figures", "joint", jointModelName))
  chainCertainty(sampled, path = file.path("./figures", "joint", jointModelName))
  for (model in unique(parPreFixs)){
    Rev <- F
    i <- i + 1
    path <- file.path("./figures", "joint", jointModelName, experiments[[i]]$modelName)
    dir.create(path, recursive = T)
    if(i == 4){
      Rev <- T
    }

    origData <- experiments[[i]]$data
    preppedData <- experiments[[i]]$preppedData
    factors <- experiments[[i]]$factors
    llFunc <- experiments[[i]]$llFunc
    currentPars <- samples$par_names[which(parPreFixs == model)]
    load(paste0("samples/", experiments[[i]]$modelName, ".RData"))
    singleSamples <- sampled
    settings <- attr(preppedData, 'settings')
    tmp <- postSimulate(samples = singleSamples, dat = origData, ll_func = llFunc, preppedData = preppedData, RL = settings$RL, n = n,
                        match = settings$match, jointSamples = samples, inputPars = c(currentPars, sharedPars), Rev = Rev)
    if(i == 1 | i == 2){
      postCDF(tmp$pp, tmp$data, joint = tmp$pp_joint, factors, PDF = T, path = path)
      ppPlot(tmp$pp, tmp$data, joint = tmp$pp_joint, factors, PDF = T, path = path)
      postHist(tmp$pp, tmp$data, joint = tmp$pp_joint, factors, PDF = T, path = path)
    }
    if(i == 3){
      RLSATPlot(tmp$pp, tmp$data, joint = tmp$pp_joint, PDF = T, path = path)
    }
    if(i == 4){
      RevPlot(tmp$pp, tmp$data, joint = tmp$pp_joint, PDF = T, path = path)
    }
  }
}


