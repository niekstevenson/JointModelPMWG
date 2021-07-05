
loadRData <- function(name, dataRoot='./data', prefix = 'data_'){
  #loads an RData file, and returns it
  fileName <- file.path(dataRoot, paste0(prefix, name, '.RData'))
  load(fileName)
  get(ls()[ls() != fileName])
}


# Utility function for data loading
loadData <- function(origData, RL = FALSE, RTLimits = c(0.15, 2), factors = NULL, subNCol = "pp",
                     exclude = TRUE, constants = NULL, match = NULL, transFunc = NULL) {
  #Compat
  if(!is.data.frame(origData)) origData <- loadRData(origData)
  origData[,subNCol] <- as.factor(as.integer(origData[,subNCol]))
  if(exclude) origData <- droplevels(origData[!origData$excl,])
  # remove fast, slow, and null responses
  origData <- origData[origData$rt>RTLimits[1] & origData$rt < RTLimits[2] & !is.na(origData$rt),]
  print(paste0("using RT limits [", RTLimits[1], ",", RTLimits[2], "]"))

  # make PMWG-style
  preppedData <- origData
  preppedData$subject <- preppedData[,subNCol]
  if(min(preppedData[,match[[1]][1]]) == 0) warning("accumulator coding should start at 1, not 0")
  preppedData$R <- factor(preppedData[,match[[1]][1]])
  
  preppedData$RT <- origData$rt
  #Only keep the columns that we're interested in
  preppedData <- preppedData[,c('subject', 'R', 'RT', factors)]
  
  
  # Add covariates and other useful things as attributes
  attr(preppedData, 'constants') <- constants
  attr(preppedData, 'match') <- match
  attr(preppedData, 'n_v') <- length(unique(origData[,match[[1]][1]]))
  attr(preppedData, 'RL') <- RL
  attr(preppedData, 'transFunc') <- transFunc
  if (RL){
    # Prepare for fitting
    cvs <- list()
    choiceIdx <- list()
    for(sub in unique(origData[,subNCol])) {
      d <- prepareForFittingRL(origData[origData[,subNCol] == sub,])
      cvs[[sub]] <- d$outcomes
      choiceIdx[[sub]] <- d$VVchoiceIdx
    }
    attr(preppedData, 'VVchoiceIdx') <- choiceIdx
    attr(preppedData, 'cvs') <- cvs

  }

  
  return(list(preppedData=preppedData, origData=origData))
}


# Functions required for fitting RL-EAM taken from the Miletic et al. 2020
prepareForFittingRL <- function(data, n_bins=5) {
  # some checks first
  if(!'reward' %in% colnames(data)) stop('Column `reward` is missing')
  if(!'stimulus_set' %in% colnames(data)) stop('Column `stimulus_set` is missing')
  if(!'rt' %in% colnames(data)) stop('Column `rt` is missing (are you using `RT`?)')
  if(!'correct' %in% colnames(data)) stop('Column `correct` is missing')
  if('choice' %in% colnames(data)) warning('Warning! Column `choice` will be overwritten')
  
  # Fix stimulus set to be numerical and have no missing numbers
  data$stimulus_set <- data$stimulus_set-min(data$stimulus_set)
  data$stimulus_set <- match(data$stimulus_set, unique(data$stimulus_set))-1
  
  # Get response "correctness" - i.e., whether a choice corresponds to the optimal choice
  # code choice in terms of upper bound (2) or lower bound (1), using the DDM convention in `rtdists`
  data$choice <- data$correct
  
  # check for outcome column, for compatibility with older data. This will be removed later.
  if('outcome' %in% colnames(data)) {
    if(max(data$outcome) == 100) {
      data$reward <- data$outcome / 100
    }
  }
  
  # define bins per stimulus; useful for later checking model fit
  data$trialN_this_stim <- NA
  for(lvl in unique(data$stimulus_set)) {
    data$trialN_this_stim[data$stimulus_set==lvl] <- seq(1, sum(data$stimulus_set==lvl))
  }
  data$bin <- as.numeric(cut(data$trialN_this_stim, n_bins))
  
  # Set-up outcome matrix
  outcomes <- matrix(NA, nrow=nrow(data), ncol=length(unique(data$stimulus_set))*2)
  for(row in 1:nrow(data)) {
    cond = data$stimulus_set[row]
    outcomes[row,(cond)*2+ifelse(data$choice[row]==1, 2, 1)] <- data$reward[row]
  }
  
  # On the basis of which alternatives is chosen? 
  # Make a matrix of nTrials x 2; first column = trialN, second column = choice number
  VVchoiceIdx <- matrix(FALSE, nrow=nrow(data), ncol=ncol(outcomes))
  for(tr in 1:nrow(data)) {
    stimulus_set <- data[tr, 'stimulus_set']
    VVchoiceIdx[tr, ((stimulus_set)*2+1):((stimulus_set)*2+2)] <- TRUE
  }
  VVchoiceIdx <- which(t(VVchoiceIdx), arr.ind = TRUE)[,2:1]  # Gives for every trial each column that is chosen
  choice <- data$choice
  return(list(data=data, VVchoiceIdx=VVchoiceIdx, outcomes=outcomes))
}
