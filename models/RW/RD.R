
addConstants <- function(pars, names, values){
  #Add possible constants to the data
  namesOld <- names(pars)
  for (value in values){
    pars[length(pars)+1] <- value
  }
  names(pars) <- c(namesOld, names)
  return(pars)
}

driftsMSIT <- function(pars, data){
  for(i in 1:3){
    FlIdx <- data[,'flank'] == i
    SimIdx <- data[,'pos'] == i #Simon effect
    MatchIdx <- data[,'uniq'] == i
    pars[FlIdx,i,'v'] <- pars[FlIdx,i,'v'] + pars[1,1,'vFlank']
    pars[SimIdx,i,'v'] <- pars[SimIdx,i,'v'] + pars[1,1,'vSimon']
    #pars[SimIdx,i,'B'] <- pars[SimIdx,i,'B'] + pars[1,1,'BSimon']
    for(j in 1:3){
      # We are currently in accumulator i, and we want to give the targets in accumulator i an extra bump based on position
      # So we loop over the three positions, and find for which values, the position == j and give a bump based on that idx
      PosIdx <- data[,'pos'] == j
      pars[MatchIdx&PosIdx,i,'v'] <- pars[MatchIdx&PosIdx,i,'v'] + pars[1,1,'vMatch'] + pars[1,1,paste0('vPos.', j)]
    }
    
  }
  return(pars)
}

driftsMSIT2 <- function(pars, data){
  for(i in 1:3){
    FlIdx <- data[,'flank'] == i
    SimIdx <- data[,'pos'] == i #Simon effect
    pars[FlIdx,i,'v'] <- pars[FlIdx,i,'v'] + pars[1,1,'vFlank']
    pars[SimIdx,i,'v'] <- pars[SimIdx,i,'v'] + pars[1,1,'vSimon']
    pars[SimIdx,i,'B'] <- pars[SimIdx,i,'B'] + pars[1,1,'BSimon']
    for(j in 1:3){
      # We are currently in accumulator i, and we want to give the targets in accumulator i an extra bump based on position
      # So we loop over the three positions, and find for which values, the position == j and give a bump based on that idx
      PosIdx <- data[,'pos'] == j
      #The bump is scaled by the number of stimuli that are presented with (1 + nScale)^(nStims - 1)
      posBoost <- pars[1,1,paste0('vPos.', j)]
      #Every stim gets the boost, based on numerical distance with exponential decay: distScale^numDist
      # distScale <- (pars[1,1,'distScale']^(abs(data[PosIdx,'uniq'] - i)))
      pars[PosIdx,i,'v'] <- pars[PosIdx,i,'v'] + (pars[1,1,'vMatch'] + posBoost) #* distScale
    }
  }
  return(pars)
}

driftsMSIT3 <- function(pars, data){
  for(i in 1:3){
    FlIdx <- data[,'flank'] == i
    PosIdx <- data[,'pos'] == i #Simon effect
    MatchIdx <- data[,'uniq'] == i
    pars[FlIdx,i,'Qs'] <- pars[FlIdx,i,'Qs'] + pars[1,1,'QFlank']
    pars[PosIdx,i,'Qs'] <- pars[PosIdx,i,'Qs'] + pars[1,1,'QSimon']
    pars[MatchIdx, i, 'Qs'] <- pars[MatchIdx, i, 'Qs'] + pars[1,1,'QMatch']
    #distScale <- (pars[1,1,'distScale']^(abs(data[,'uniq']) - i))
    #pars[,i,'Qs'] <- pars[,i,'Qs'] + (pars[1,1,'QMatch']) * distScale
    #pars[PosIdx,,'wD'] <- pars[PosIdx,,'wD'] + pars[1,1,paste0('wDPos.', i)]*((1 + pars[1,1,'nScale'])^(data[PosIdx,'nstims'] - 1))
  }
  return(pars)
}


driftsRePar <- function(pars, ...){
  v1Pars <- pars[grep("V0_", names(pars), fixed = F)]
  v2Pars <- v1Pars
  diffPars <- pars[grep("diff_", names(pars), fixed = F)]
  if(length(diffPars) > 0){
    v1Pars <- v1Pars - diffPars/2
    v2Pars <- v2Pars + diffPars/2
  }
  names(v1Pars) <- paste("v_1_", gsub("^[^_]*_", "", names(v1Pars)), sep = "")
  names(v2Pars) <- paste("v_2_", gsub("^[^_]*_", "", names(v2Pars)), sep = "")
  pars <- pars[!grepl("V0_", names(pars)) & !grepl("diff_", names(pars))]
  pars <- c(pars, v1Pars, v2Pars)
}

prepPars <- function(pars, data){
  #This is where the magic happens
  
  settings <- attr(data, 'settings')
  transFunc <- settings$transFunc
  RL <- settings$RL
  constants <- settings$constants
  match <- settings$match
  matchedPars <- names(match) #get the pars that are different for e.g. response coded 
  n_v <- settings$n_v
  
  if(!is.null(transFunc)){
    if(transFunc$earlyTransform) pars <- transFunc$func(pars, data)
  }
  # Find all pars that are varying between conditions
  varyPars <- pars[grep("[_]", names(pars), fixed = F)]
  pars <- pars[!names(pars) %in% names(varyPars)]
  pars <- addConstants(pars, names(constants), sapply(constants, FUN = function(x) x))
  # Get the common parameters among the varying pars so that we can create an array using these and the standard parameters. 
  commonPars <- unique(gsub("*_.*", "", names(varyPars))) #e.g. V0.SPD & V0.ACC -> V0 as common parameter
  pars <- addConstants(pars, commonPars, rep(NA, length(commonPars)))

  # Create an array with dimensions c(Available responses, ntrials, npars)
  pars <- array(rep(pars, each = nrow(data)*n_v), dim = c(nrow(data),n_v, length(pars)), dimnames = list(NULL, NULL, names(pars)))
  
  # For all the parameters varying between condition, of which their base parameter was added to the array as NA,
  # Fill in the NAs with the actual parameter values
  # Varying pars should be denoted with a <basePar>_<factor>.<level>_<factor.<level> format.
  # Or if part of the coding <basePar>_<accumulatorNr>_<factor>.<level>_<factor.<level> format.
  for (par in names(varyPars)){
    parSplit <- strsplit(par, "_")[[1]] 
    basePar <- parSplit[1]
    #Check if we are modifying a parameter we have multiple accumulators per response or not.
    if(basePar %in% matchedPars){
      #remove base par and it's equivalent accumulator from the list
      parN <- as.integer(parSplit[2]) #e.g. response code 1 = wrong 2 = correct. 
      factorPars <- parSplit[-c(1,2)] #remaining are the factor(s) in format <factor>.<level>
      idx <- which(basePar == matchedPars) #We need to know which matched parameter we're talking about
    } else {
      idx <- 1
      #Not a parameter involved in match map
      factorPars <- parSplit[-1] #remaining are the factor(s) in format <factor>.<level>
    }
    
    
    criterion <- list()
    parCoded <- F #Keep score of whether the current parameter is spread over accumulators or row dependent.
    for(factorPar in factorPars){
      factorSplit <- strsplit(factorPar, "[.]")[[1]]
      if(factorSplit[1] == match[[idx]][1]){ #the current factor is spread over accumulators not on trial-dependent level. 
        currentAcc <- factorSplit[2] #select current accumulator
        parCoded <- T #Remember that this parameter is (partially) spread over accumulators. 
      } else{ #Find matching values based on trial(row)-dependent factors
        criterion[[factorPar]] <- which(data[,factorSplit[1]] == factorSplit[2])
        if(!parCoded) currentAcc <- 1:n_v #Not a accumulator dependent factor, select all. 
      }
    }
    criterion <- Reduce(intersect, criterion) #calculate for which rows we should replace considering all relevant factors
    if(basePar %in% matchedPars & !is.na(match[[idx]][2])){ #For this parameter the values are determined on the match map e.g. response correctness.
      for(currentAcc in 1:n_v){ #Loop over different accumulators and find when it matches the specified match mapping. 
        if(parN == 2) #If e.g v_2 this means drift rate for correct. v_1 is drift rate for incorrect (based on DDM convention) 
          respCrit <- which(currentAcc == data[,match[[idx]][2]]) #See for which values the currentAccumulator meets the predefined correct map. 
        else{
          respCrit <- which(currentAcc != data[,match[[idx]][2]])
        }
        #update criterion (whilst storing old criterion for other accumulators.)
        if(!is.null(criterion)) respCrit <- intersect(respCrit, criterion)
        pars[respCrit,currentAcc,basePar] <- varyPars[par]
      }
    } else{
      if(basePar %in% matchedPars & is.na(match[[idx]][2])){
        #In this case there's no response coding, just accuracy coding
        currentAcc <- parN  
        if(is.null(criterion)) criterion <- 1:nrow(pars) #no criteria, fill in whole column
      }
      pars[criterion,currentAcc,basePar] <- varyPars[par]
    }
  }
  #Do some transformations
  #pars[,,'t0'] <- pnorm(pars[,,'t0'])
  #pars[,,'t0'][pars[,,'t0'] < 1e-5] <- 0 #protect against floating point error
  if(!is.null(transFunc)){
    if(!transFunc$earlyTransform) pars <- transFunc$func(pars, data)
  }
  
  return(pars)
}


transform.RL <- function(pars, data) {
  sub <- data$subject[1]
  pars[,,'aV'] <- pnorm(pars[,,'aV'])
  cvs <- attr(data,"cvs")[[as.character(sub)]]
  choiceIdx <- attr(data, "VVchoiceIdx")[[as.character(sub)]]
  
  ### SM
  
  # Create start point vector
  startValues <- rep(pars[,1,'Qs'], each=2, times=ncol(cvs))
  
  # learning rates matrix
  learningRates <- matrix(rep(pars[,1,'aV'], each=ncol(cvs)), ncol=ncol(cvs), byrow=TRUE)
  
  # call C
  updated <- adapt.c.dmc(startValues = startValues, 
                         learningRates = learningRates, 
                         feedback = cvs,
                         learningRule='SARSA')
  
  # add back data to 'pars' array
  pars[,,'Qs'] <- matrix(updated$adaptedValues[choiceIdx], ncol=2, byrow=T)[,c(2,1)]
  #Note that this returns them in the format Qvalue 1 is correct response, Qvalue 2 is incorrect response. 
  #So we switch them around in the next line to match DDM convention (which we keep throughout the rest of the code)
  return(pars)
}



likelihood.RD <- function(pars,data,min.like=1e-10, sample = F)   
{
  
  facs <- colnames(data)[-which(colnames(data) %in% c("subject", "RT", "R"))]
  if(!sample){
    if(any(pars['t0'] > data$rt)) return(-1e5)
  }
  
  pars <- prepPars(pars, data)
  if (sample){
    #We're interested in posterior prediction, instead of fitting
    n = nrow(data)
    out <- rWaldRace(n=n, A=pars[,,'A'], v=pars[,,'v'], B=pars[,,'B'], t0=pars[,,'t0'], s=pars[,,'s'])
    out <- cbind(out, data[,facs])
    colnames(out)[-c(1,2)] <- facs
    return(out)
  } else{

    #Get the likelihoods for the sampler
    return(sum(log(pmax(dWaldRace(rt=data$RT, response=data$R,
                                  A=pars[,,'A'], v=pars[,,'v'], B=pars[,,'B'], t0=pars[,,'t0'], st0= 0, s=pars[,,'s'],
                                  silent=TRUE), min.like, na.rm=TRUE))))
  }
}

likelihood.ARD <- function(pars,data,min.like=1e-10, sample = F)   
{
  
  facs <- colnames(data)[-which(colnames(data) %in% c("subject", "RT", "R"))]
  if(!sample){
    #Protection against bad ranges/values + speed up
    if(any(pars['t0'] > data$rt)) return(-1e5)
    if(pars['aV'] < -3.5) return(-1e5)
    if(pars['aV'] > 3) return(-1e5)
    if(pars['wD'] < 0) return(-1e5)
    if(pars['wD'] > 20) return(-1e5)
    if(pars['wS'] < -10) return(-1e5)
    if(pars['wS'] > 20) return(-1e5)
  }
  pars <- prepPars(pars, data)
  if (sample){
    #We're interested in posterior prediction, instead of fitting
    n = nrow(data)
    out <- rARD(n=n, A=pars[,,'A'], V0=pars[,,'V0'], B=pars[,,'B'], t0=pars[,,'t0'], s=pars[,,'s'],
                       wS = pars[,,'wS'], wD = pars[,,'wD'], Qs = pars[,,'Qs'])
    
    out <- cbind(out, data[,facs])
    colnames(out)[-c(1,2)] <- facs
    return(out)
  } else{
    #Get the likelihoods for the sampler
    return(sum(log(pmax(dARD(rt=data$RT, response=data$R,
                                  A=pars[,,'A'], V0=pars[,,'V0'], B=pars[,,'B'], t0=pars[,,'t0'], st0= 0, s=pars[,,'s'],
                                  wS = pars[,,'wS'], wD = pars[,,'wD'], Qs = pars[,,'Qs'],
                                  silent=TRUE), min.like, na.rm=TRUE))))
  }
}
