prepDataJoint <- function(experiments, n_sessions = 2, slack = 0.1, hardCutAcc = TRUE, hardCutN = TRUE, splitSess = F){
  for(name in names(experiments)){
    dat <- loadRData(name)
    dat <- dat[dat$block!= 9,]
    #Make sure they all have the same name
    dat <- dat %>% rename(pp = experiments[[name]]$subNCol) #this is a dumb workaround
    dat$accuracy <- dat[,experiments[[name]]$accCol] 
    n_trialsPerSession <- trialsPerSession(dat)

    if(splitSess){
      #If we want to compare sessions, create two datasets and model + settings in the experiments list.
      n_sessions <- 1
      for (i in 1:2){
        if(i == 2){
          experiments[[paste0(name, "2")]] <- experiments[[name]]
          name <- paste0(name, "2")
        } 
        tmp <- dat[dat$sessionNr == i,]
        #Filter based on acc criterion and necessary number of trials
        tmp <- filterSubs(tmp, n_trialsPerSession*n_sessions, experiments[[name]]$accCriterion, slack) 

        experiments[[name]]$exclAccLocal <- length(unique(tmp$pp[tmp$exclAcc]))/length(unique(tmp$pp))
        experiments[[name]]$exclNLocal <- length(unique(tmp$pp[tmp$exclN]))/length(unique(tmp$pp))
        
        #Store the experiments and which subjects are excluded on what reasons
        experiments[[name]]$data <- tmp
        experiments[[name]]$completed <- unique(tmp$pp[!tmp$exclN])
        experiments[[name]]$accCritMet <- unique(tmp$pp[!tmp$exclAcc])
        if(i == 1) modelName <- experiments[[name]]$modelName
        experiments[[name]]$modelName <- paste0(modelName, "ses", i)
      }
    } else{
      #Filter based on acc criterion and necessary number of trials
      dat <- filterSubs(dat, n_trialsPerSession*n_sessions, experiments[[name]]$accCriterion, slack)
      
      experiments[[name]]$exclAccLocal <- length(unique(dat$pp[dat$exclAcc]))/length(unique(dat$pp))
      experiments[[name]]$exclNLocal <- length(unique(dat$pp[dat$exclN]))/length(unique(dat$pp))
      
      #Store the experiments and which subjects are excluded on what reasons
      experiments[[name]]$data <- dat
      experiments[[name]]$completed <- unique(dat$pp[!dat$exclN])
      experiments[[name]]$accCritMet <- unique(dat$pp[!dat$exclAcc])
    }
  }
  #Filter based on other experiments excluded, and on global settings
  completedAll <- Reduce(intersect, lapply(experiments, FUN = function(x){x[['completed']]}))
  accCritMetAll <- Reduce(intersect, lapply(experiments, FUN = function(x){x[['accCritMet']]}))
  
  for(name in names(experiments)){
    #Filter out the participants that did not make the cut in other experiments 
    if(hardCutAcc & hardCutN){
      criterion <- intersect(completedAll, accCritMetAll)
      dat <- experiments[[name]]$data
      #Seems to be high, not sure if calculated correctly
      experiments[[name]]$exclAccGlobal <- length(which(!unique(dat$pp[!dat$exclAcc]) %in% accCritMetAll))/
        length(unique(dat$pp[!dat$exclAcc]))
      experiments[[name]]$exclTrialNGlobal <- length(which(!unique(dat$pp[!dat$exclN]) %in% completedAll))/
        length(unique(dat$pp[!dat$exclN]))

      experiments[[name]]$data$excl[!(dat$pp %in% criterion)] <- T
    } else if(hardCutN){
      dat <- experiments[[name]]$data
      experiments[[name]]$exclTrialNGlobal <- length(which(!unique(dat$pp[!dat$exclN]) %in% completedAll))/
        length(unique(dat$pp[!dat$exclN]))
      experiments[[name]]$exclAccGlobal <- 0
      experiments[[name]]$data$excl[!(dat$pp %in% completedAll)] <- T
    } else if(hardCutAcc){
      dat <- experiments[[name]]$data
      experiments[[name]]$exclAccGlobal <- length(which(!unique(dat$pp[!dat$exclAcc]) %in% accCritMetAll))/
        length(unique(dat$pp[!dat$exclAcc]))
      experiments[[name]]$exclTrialNGlobal <- 0
      experiments[[name]]$data$excl[!(dat$pp %in% accCritMetAll)] <- T
    }
    filteredData <- loadData(experiments[[name]]$data,
                          RL = experiments[[name]]$RL,
                          RTLimits = experiments[[name]]$RTLimits,
                          factors = experiments[[name]]$factors,
                          match = experiments[[name]]$match,
                          constants = experiments[[name]]$constants,
                          transFunc = experiments[[name]]$transFunc)
    #data to be fit
    experiments[[name]]$preppedData <- filteredData$preppedData
    #Updated original data, that only includes the participants that made the cut
    experiments[[name]]$data <- filteredData$origData
  }
  return(experiments)
}

trialsPerSession <- function(dat){
  #Get trial counts
  tmp <- dat %>% group_by(pp) %>% tally()
  #Get most frequent occurence
  mostFreq <- data.frame(sort(table(tmp$n),decreasing=TRUE), stringsAsFactors = F)
  mostFreq <- as.numeric(levels(droplevels(mostFreq[,1][mostFreq[,2] > 2])))
  #Get the n trials that are most often associated with a session. Assumes two sessions!
  n_trials <- intersect(mostFreq, round(mostFreq/2))[1]
  return(n_trials)
}

filterSubs <- function(dat, n_trials, accCriterion, slack){
  dat$exclAcc <- F
  dat$exclN <- F
  dat$accuracy[which(is.na(dat$accuracy))] <- F
  for (sub in unique(dat$pp)){
    tmp <- dat[dat$pp == sub,]
    #Filter based on accuracy and number of trials separately
    #Cut them some slack with regards to the number of necessary trials. 
    if(nrow(tmp) < n_trials*(1-slack) | nrow(tmp) > n_trials*(1+slack)) tmp$exclN <- T
    if(is.na(mean(mean(tmp$accuracy, na.rm = T))) | (mean(tmp$accuracy, na.rm = T) < accCriterion)) tmp$exclAcc <- T
    dat[dat$pp == sub,] <- tmp
  }
  dat$excl <- dat$exclN | dat$exclAcc
  return(dat)
}