
addConstants <- function(pars, names, values){
  #Add possible constants to the data
  namesOld <- names(pars)
  for (value in values){
    pars[length(pars)+1] <- value
  }
  names(pars) <- c(namesOld, names)
  return(pars)
}

generalTransf <- function(pars, data){
  for (name in names(transFunc)){
    transfPar <- transFunc[[name]]
    prePars <- pars[grep(transfPar[[2]], substr(names(pars),1, nchar(transfPar[[2]])), fixed = F)]
    newPars <- transfPar[[1]](pars[name], prePars)
    if(length(transfPar) > 2){
      if(transfPar[[3]]) names(newPars) <- gsub(transfPar[[2]], name, names(newPars))
    }
    pars <- c(pars[!names(pars) %in% c(names(newPars), name)], newPars)
  }
}


competingTransf <- function(pars, ...){
  sPars <- pars[grep("s_", names(pars), fixed = F)]
  piPars <- pars[grep("pi_", names(pars), fixed = F)]
  piPars <- pnorm(piPars)
  v1Pars <- sPars*(1-piPars)
  v2Pars <- sPars*piPars
  names(v1Pars) <- paste("v_1_", gsub("^[^_]*_", "", names(piPars)), sep = "")
  names(v2Pars) <- paste("v_2_", gsub("^[^_]*_", "", names(piPars)), sep = "")
  pars <- pars[!grepl("s_", names(pars)) & !grepl("pi_", names(pars))]
  pars <- c(pars, v1Pars, v2Pars)
  return(pars)
}

driftsRePar2 <- function(pars, ...){
  V0Pars <- pars[grep("V0_", names(pars), fixed = F)]
  diffPars <- pars[grep("diff_", names(pars), fixed = F)]
  diffPlusPars <- pars[grep("diffPlus_", names(pars), fixed = F)]
  if(!is.null(diffPars)){
    v1Pars <- V0Pars - diffPars/2
    v2Pars <- V0Pars + diffPars/2
  }
  if(!is.null(diffPlusPars)){
    v1ParsNew <- numeric()
    v2ParsNew <- numeric()
    i <- 0
    baseDiffPlusPars <- gsub("[.].*", "", names(diffPlusPars))
    for (basePar in baseDiffPlusPars){
      curPars <- pars[grep(basePar, names(pars), fixed = F)]
      for (par in curPars){
        i <- i + 1
        tmp1 <- v1Pars - par/2
        names(tmp1) <- paste(names(tmp1), "_", gsub("^[^_]*_", "", names(curPars)[i]), sep = "")
        tmp2 <- v2Pars + par/2 
        names(tmp2) <- paste(names(tmp2), "_", gsub("^[^_]*_", "", names(curPars)[i]), sep = "")
        v1ParsNew <- c(v1ParsNew, tmp1)
        v2ParsNew <- c(v2ParsNew, tmp2)
      }
    }
  }
}


driftsRePar <- function(pars, ...){
  V0Pars <- pars[grep("V0_", names(pars), fixed = F)]
  diffPars <- pars[grep("diff_", names(pars), fixed = F)]
  diffPlusPars <- pars[grep("diffPlus_", names(pars), fixed = F)]
  if(!is.null(diffPars)){
    v1Pars <- V0Pars - diffPars/2
    v2Pars <- V0Pars + diffPars/2
  }
  names(v1Pars) <- paste("v_1_", gsub("^[^_]*_", "", names(V0Pars)), sep = "")
  names(v2Pars) <- paste("v_2_", gsub("^[^_]*_", "", names(V0Pars)), sep = "")
  if(!is.null(diffPlusPars)){
    v1ParsNew <- numeric()
    v2ParsNew <- numeric()
    i <- 0
    for(diffPlusPar in diffPlusPars){
      i <- i + 1
      tmp1 <- v1Pars - diffPlusPar/2
      names(tmp1) <- paste(names(tmp1), "_", gsub("^[^_]*_", "", names(diffPlusPars)[i]), sep = "")
      tmp2 <- v2Pars + diffPlusPar/2 
      names(tmp2) <- paste(names(tmp2), "_", gsub("^[^_]*_", "", names(diffPlusPars)[i]), sep = "")
      v1ParsNew <- c(v1ParsNew, tmp1)
      v2ParsNew <- c(v2ParsNew, tmp2)
    }
    v1Pars <- v1ParsNew
    v2Pars <- v2ParsNew
  }

  pars <- pars[!grepl("V0_", names(pars)) & !grepl("diff_", names(pars)) & !grepl("diffPlus_", names(pars))]
  pars <- c(pars, v1Pars, v2Pars)
}


prepPars <- function(pars, data){
  #This is where the magic happens

  transFunc <- attr(data, "transFunc")
  RL <- attr(data, "RL")
  constants <- attr(data, "constants")
  match <- attr(data, "match")
  matchedPars <- names(match) #get the pars that are different for e.g. response coded 
  n_v <- attr(data, "n_v")
  
  if(RL) transFunc <- c(transform.RL, F)
  if(!is.null(transFunc)){
    if(transFunc[[2]]) pars <- transFunc[[1]](pars, data)
  }
  # Find all pars that are varying between conditions
  varyPars <- pars[grep("[_]", names(pars), fixed = F)]
  pars <- pars[!names(pars) %in% names(varyPars)]
  

  
  # Get the common parameters among the varying pars so that we can create an array using these and the standard parameters. 
  tmp <- unique(gsub("*_.*", "", names(varyPars))) #e.g. V0.SPD & V0.ACC -> V0 as common parameter
  pars <- addConstants(pars, tmp, rep(NA, length(tmp)))
  if(!"v" %in% names(pars)) pars <- addConstants(pars, "v", NA) #Done for the RL advantage framework, makes sure there's a drift rate to fill.
  
  pars <- addConstants(pars, names(constants), sapply(constants, FUN = function(x) x))
  # Create an array with dimensions c(Available responses, ntrials, npars)
  pars <- array(rep(pars, each = nrow(data)*n_v), dim = c(n_v, nrow(data), length(pars)), dimnames = list(NULL, NULL, names(pars)))
  
  # For all the parameters varying between condition, of which their base parameter was added to the array as NA,
  # Fill in the NAs with the actual parameter values
  # Varying pars should be denoted with a <basePar>_<factor>.<level>_<factor.<level> format.
  # Or if part of the coding <basePar>_<accumulatorNr>_<factor>.<level>_<factor.<level> format.
  for (par in names(varyPars)){
    idx <- 1
    tmp1 <- strsplit(par, "_")[[1]] 
    basePar <- tmp1[1]
    #Check if we are modifying a parameter we have multiple accumulators per response or not.
    if(basePar %in% matchedPars){
      #remove base par and it's equivalent accumulator from the list
      parN <- as.integer(tmp1[2]) #e.g. response code 1 = wrong 2 = correct. 
      factorPars <- tmp1[c(-1,-2)] #remaining are the factor(s) in format <factor>.<level>
      idx <- which(basePar == matchedPars) #We need to know which matched parameter we're talking about
    } else {
      #Not a parameter involved in match map
      factorPars <- tmp1[-1] #remaining are the factor(s) in format <factor>.<level>
    }


    criterion <- list()
    parCoded <- F #Keep score of whether the current parameter is spread over accumulators or row dependent.
    for(factorPar in factorPars){
      tmp2 <- strsplit(factorPar, "[.]")[[1]]
      if(tmp2[1] == match[[idx]][1]){ #the current factor is spread over accumulators not on trial-dependent level. 
        currentAcc <- as.integer(tmp2[2]) #select current accumulator
        parCoded <- T #Remember that this parameter is (partially) spread over accumulators. 
      } else{ #Find matching values based on trial(row)-dependent factors
        criterion[[factorPar]] <- which(data[,tmp2[1]] == tmp2[2])
        if(!parCoded) currentAcc <- 1:n_v #Not a accumulator dependent factor, select all. 
      }
    }
    criterion <- Reduce(intersect, criterion) #calculate for which rows we should replace considering all relevant factors
    if(basePar %in% matchedPars & !is.na(match[[idx]][2])){ #For this parameter the values are determined on the match map e.g. response correctness.
      for(currentAcc in 1:n_v){ #Loop over different accumulators and find when it matches the specified match mapping. 
        if(parN == 2) #If e.g v_2 this means drift rate for correct. v_1 is drift rate for incorrect. 
          respCrit <- which(currentAcc == data[,match[[idx]][2]]) #See for which values the currentAccumulator meets the predefined correct map. 
        else{
          respCrit <- which(currentAcc != data[,match[[idx]][2]])
        }
        #update criterion (whilst storing old criterion for other accumulators.)
        if(!is.null(criterion)) respCrit <- intersect(respCrit, criterion)
        pars[currentAcc,respCrit,basePar] <- varyPars[par]
      }
    } else{
      if(basePar %in% matchedPars & is.na(match[[idx]][2])){
        #In this case the 
        currentAcc <- parN  
        criterion <- 1:ncol(pars)
      }
      pars[currentAcc,criterion,basePar] <- varyPars[par]
    }
  }
  #Do some transformations
  pars[,,'t0'] <- pnorm(pars[,,'t0'])
  pars[,,'t0'][pars[,,'t0'] < -5] <- 100 #This is a bit cheating, but makes it so that t0s don't approach minus infinity get a low likelihood
  if(!is.null(transFunc)){
    if(!transFunc[[2]]) pars <- transFunc[[1]](pars, data)
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
  startValues <- rep(pars[1,,'SR'], each=2, times=ncol(cvs))
  
  # learning rates matrix
  learningRates <- matrix(rep(pars[1,,'aV'], each=ncol(cvs)), ncol=ncol(cvs), byrow=TRUE)
  
  # call C
  updated <- adapt.c.dmc(startValues = startValues, 
                         learningRates = learningRates, 
                         feedback = cvs,
                         learningRule='SARSA')
  
  # add back data to 'pars' array
  pars[,,'SR'] <- matrix(updated$adaptedValues[choiceIdx], ncol=2, byrow=FALSE)
  
  pars[,,'v'] <- t(cbind(pars[1,,"V0"] + pars[1,,"wD"]*(pars[2,,"SR"]-pars[1,,"SR"]) + pars[1,,'wS'] * (pars[2,,"SR"]+pars[1,,"SR"]),
                  pars[2,,"V0"] + pars[2,,"wD"]*(pars[1,,"SR"]-pars[2,,"SR"]) + pars[2,,'wS'] * (pars[2,,"SR"]+pars[1,,"SR"])))
  
  
  return(pars)
}


likelihood.RD <- function(pars,data,min.like=1e-10, sample = F)   
{
  
  facs <- colnames(data)[-which(colnames(data) %in% c("subject", "RT", "R"))]
  
  pars <- prepPars(pars, data)

  if (sample){
    #We're interested in posterior prediction, instead of fitting
    n = nrow(data)
    out <- rWaldRaceSM(n=n, 
                       A=lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'A'])),
                       v=lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'v'])),
                       B=lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'B'])),
                       t0=lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'t0'])),
                       st0= 0,
                       s=lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'s'])),
                       silent=TRUE, simPerTrial=FALSE)
    
    out <- cbind(out, data[,facs])
    colnames(out)[-c(1,2)] <- facs
    return(out)
  } else{
    #Get the likelihoods for the sampler
    return(sum(log(pmax(dWaldRace(rt=data$RT, response=data$R,
                           A=lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'A'])),
                           v=lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'v'])),
                           B=lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'B'])),
                           t0=lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'t0'])),
                           st0= 0,
                           s=lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'s'])),
                           silent=TRUE), min.like, na.rm=TRUE))))
  }
}


likelihood.TRD <- function(pars,data,min.like=1e-10, sample = F)
{
  
  facs <- colnames(data)[-which(colnames(data) %in% c("subject", "RT", "R"))]
  
  pars <- prepPars(pars, data)
  
  print(dimnames(pars))
  if (sample){
    #We're interested in posterior prediction, instead of fitting
    n = nrow(data)
    out <- rWaldRace_Timing_NS(n=n, 
                       A=lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'A'])),
                       v_E =lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'v'])),
                       B_E =lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'B'])),
                       t0E = lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'t0'])),
                       st0= 0,
                       s_E =lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'s'])),
                       v_T = lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'vT'])),
                       s_T = lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'sT'])),
                       B_T = lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'BT'])),
                       t0T=lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'t0T'])),
                       silent=TRUE, simPerTrial=FALSE)
    
    out <- cbind(out, data[,facs])
    colnames(out)[-c(1:3)] <- facs
    return(out)
  } else{
    #Get the likelihoods for the sampler
    return(sum(log(pmax(dWaldRaceTiming(rt=data$RT, response=data$R,
                                  A=lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'A'])),
                                  v=lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'v'])),
                                  B=lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'B'])),
                                  t0=lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'t0'])),
                                  st0= 0,
                                  s=lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'s'])),
                                  v_T = lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'vT'])),
                                  s_T = lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'sT'])),
                                  B_T = lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'BT'])),
                                  A_T = lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'AT'])),
                                  t0T=lapply(1:nrow(pars), FUN = function(x) return(pars[x,,'t0T'])),
                                  silent=TRUE), min.like, na.rm=TRUE))))
  }
}





