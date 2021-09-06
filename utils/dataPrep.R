prepDataJoint <- function(experiments, n_sessions = 2, slack = 0.1, hardCutAcc = TRUE, hardCutN = TRUE, splitSess = F, practiceBlock = 9){
  #The current function takes in an experiments object and makes (multiple) tasks ready for fitting
  for(name in names(experiments)){
    dat <- loadRData(name)
    dat <- dat[dat$block!= practiceBlock,] #Standard number for practice block = 9
    #Make sure that all columns have the proper name
    dat$subject <- dat[,experiments[[name]]$subNCol]
    dat$accuracy <- dat[,experiments[[name]]$accCol]
    dat$RT <- dat[,experiments[[name]]$rtCol]
    
    #Determine how many trials there are per session
    n_trialsPerSession <- trialsPerSession(dat, n_sessions)
    
    if(splitSess){
      #If we want to compare sessions, create two datasets and model + settings in the experiments list.
      for (i in 1:n_sessions){
        if(i != 1){
          experiments[[paste0(name, i)]] <- experiments[[name]]
          name <- paste0(name, i)
        } 
        tmp <- dat[dat$sessionNr == i,]
        #Filter based on acc criterion and necessary number of trials
        tmp <- filterSubs(tmp, n_trialsPerSession, experiments[[name]]$accCriterion, slack) 
        #Store the experiments and which subjects are excluded on what reasons
        experiments[[name]]$data <- tmp
        
        experiments[[name]]$completed <- unique(tmp$subject[!tmp$exclN])
        experiments[[name]]$accCritMet <- unique(tmp$subject[!tmp$exclAcc])
        
        if(i == 1) modelName <- experiments[[name]]$modelName
        experiments[[name]]$modelName <- paste0(modelName, "ses", i)
      }
    } else{
      #Filter based on acc criterion and necessary number of trials
      dat <- filterSubs(dat, n_trialsPerSession*n_sessions, experiments[[name]]$accCriterion, slack)
      
      #Store the experiments and which subjects are excluded on what reasons
      experiments[[name]]$data <- dat
      experiments[[name]]$completed <- unique(dat$subject[!dat$exclN])
      experiments[[name]]$accCritMet <- unique(dat$subject[!dat$exclAcc])
    }
  }
  #Filter based on other experiments excluded, and on global settings
  completedAll <- Reduce(intersect, lapply(experiments, FUN = function(x){x[['completed']]}))
  accCritMetAll <- Reduce(intersect, lapply(experiments, FUN = function(x){x[['accCritMet']]}))
  
  for(name in names(experiments)){
    #Filter out the participants that did not make the cut in other experiments 
    if(hardCutAcc & hardCutN){
      criterion <- intersect(completedAll, accCritMetAll)
      experiments[[name]]$data$excl[!(experiments[[name]]$data$subject %in% criterion)] <- T
    } else if(hardCutN){
      experiments[[name]]$data$excl[!(experiments[[name]]$data$subject %in% completedAll)] <- T
    } else if(hardCutAcc){
      experiments[[name]]$data$excl[!(experiments[[name]]$data$subject %in% accCritMetAll)] <- T
    }
    filteredData <- loadData(experiments[[name]]$data,
                             RTLimits = experiments[[name]]$RTLimits,
                             factors = experiments[[name]]$factors,
                             match = experiments[[name]]$match)
    #data to be fit
    preppedData <- filteredData$preppedData
    # Add other useful things as attributes

    if(experiments[[name]]$RL){
      preppedData <- prepareForFittingRL(filteredData$origData, preppedData)
      experiments[[name]]$transFunc <- list(func = transform.RL, earlyTransform = F)
    } 
    
    attr(preppedData, 'settings') <- list(constants = experiments[[name]]$constants, 
                                          match = experiments[[name]]$match, 
                                          n_v = length(unique(preppedData$R)),
                                          RL = experiments[[name]]$RL,
                                          transFunc = experiments[[name]]$transFunc)
    
    #Updated original data, that only includes the participants that made the cut
    experiments[[name]]$preppedData <- preppedData
    experiments[[name]]$data <- filteredData$origData
  }
  return(experiments)
}

# Utility function for data loading
loadData <- function(origData, RTLimits = c(0.15, 2), factors = NULL, exclude = TRUE, match) {
  origData$subject <- as.factor(as.integer(origData$subject))
  if(exclude) origData <- droplevels(origData[!origData$excl,])
  # remove fast, slow, and null responses
  origData <- origData[origData$rt>RTLimits[1] & origData$rt < RTLimits[2] & !is.na(origData$rt),]
  print(paste0("using RT limits [", RTLimits[1], ",", RTLimits[2], "]"))
  
  # make PMWG-style
  preppedData <- origData
  if(min(preppedData[,match[[1]][1]]) == 0) warning(paste0("Following DDM convention, accumulators start at 1 \n
                                                    Your ", match[[1]][1], ' column starts at 0, we added 1'))
  preppedData$R <- factor(preppedData[,match[[1]][1]])
  #Only keep the columns that we're interested in
  preppedData <- preppedData[,c('subject', 'R', 'RT', factors)]
  
  return(list(preppedData=preppedData, origData=origData))
}


loadRData <- function(name, dataRoot='./data', prefix = 'data_'){
  #loads an RData file, and returns it
  fileName <- file.path(dataRoot, paste0(prefix, name, '.RData'))
  load(fileName)
  get(ls()[ls() != fileName])
}

trialsPerSession <- function(dat, n_sessions){
  #Get trial counts
  tmp <- dat %>% group_by(subject) %>% tally()
  #Get most frequent occurence
  mostFreq <- data.frame(sort(table(tmp$n),decreasing=TRUE), stringsAsFactors = F)
  mostFreq <- as.numeric(levels(droplevels(mostFreq[,1][mostFreq[,2] > n_sessions])))
  #Get the n trials that are most often associated with a session. Assumes two sessions!
  n_trials <- intersect(mostFreq, round(mostFreq/n_sessions))[1]
  return(n_trials)
}

filterSubs <- function(dat, n_trials, accCriterion, slack){
  dat$exclAcc <- F
  dat$exclN <- F
  dat$accuracy[which(is.na(dat$accuracy))] <- F
  for (sub in unique(dat$subject)){
    tmp <- dat[dat$subject == sub,]
    #Filter based on accuracy and number of trials separately
    #Cut them some slack with regards to the number of necessary trials. 
    if(nrow(tmp) < n_trials*(1-slack) | nrow(tmp) > n_trials*(1+slack)) tmp$exclN <- T
    if(is.na(mean(mean(tmp$accuracy, na.rm = T))) | (mean(tmp$accuracy, na.rm = T) < accCriterion)) tmp$exclAcc <- T
    if(!tmp$exclN) tmp$subject[((n_trials/2)+1):nrow(tmp)] <- tmp$subject + 1000 
    dat[dat$subject == sub,] <- tmp
  }
  dat$excl <- dat$exclN | dat$exclAcc
  return(dat)
}

# Functions required for fitting RL-EAM taken from the Miletic et al. 2020
prepareForFittingRL <- function(origData, preppedData) {
  # some checks first
  if(!'reward' %in% colnames(origData)) stop('Column `reward` is missing')
  if(!'stimulus_set' %in% colnames(origData)) stop('Column `stimulus_set` is missing')
  if(!'rt' %in% colnames(origData)) stop('Column `rt` is missing (are you using `RT`?)')
  if(!'correct' %in% colnames(origData)) stop('Column `correct` is missing')
  if('choice' %in% colnames(origData)) warning('Warning! Column `choice` will be overwritten')
  
  cvs <- list()
  choiceIdx <- list()
  for(sub in unique(origData$subject)) {
    d <- getSubChoicesRL(origData[origData$subject == sub,])
    cvs[[sub]] <- d$outcomes
    choiceIdx[[sub]] <- d$VVchoiceIdx
  }
  attr(preppedData, 'cvs') <- cvs
  attr(preppedData, 'VVchoiceIdx') <- choiceIdx
  return(preppedData)
}


getSubChoicesRL <- function(data){
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
  
  return(list(VVchoiceIdx=VVchoiceIdx, outcomes=outcomes))
}

