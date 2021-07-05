library(bayesplot)
library(ggpubr)
library(corrplot)
library(snowfall)

source("utils/utilsPP.R")

# Utility functions for the user ------------------------------------------
# Added by Reilly Innes/Niek Stevenson
pmwg_IC <- function(samples, preppedData){
  
  # Mean log-likelihood of the overall (samples-stage) model, for each subject
  mean_ll <- apply(samples$samples$subj_ll[, samples$samples$stage == "sample"], 1, mean)
  
  # Mean of each parameter across iterations.
  # Keep dimensions for parameters and subjects
  mean_pars <- t(apply(samples$samples$alpha[,, samples$samples$stage == "sample"], 1:2, mean))
  
  # Name 'mean_pars' so it can be used by the log_like function
  colnames(mean_pars) <- samples$par_names
  
  # log-likelihood for each subject using their mean parameter vector
  mean_pars_ll <- numeric(ncol(mean_pars))
  data <- preppedData
    
    
  for (sub in unique(data$subject)) {
    mean_pars_ll[sub] <- samples$ll_func(mean_pars[sub, ], data = data[data$subject == sub,])
  }
  
  # mean deviance(-2*ll of all data) 
  # effective number of parameters(-2*ll of all data - -2*ll with mean theta)
  pD <- sum(-2 * mean_ll + 2 * mean_pars_ll)
  
  # DIC = mean deviance + effective number of parameters
  DIC <- sum(-4 * mean_ll + 2 * mean_pars_ll)
  
  # BPIC = mean deviance + 2*effective number of parameters 
  # Note this is the "easy" BPIC, instead of the complex 2007 one
  BPIC <- sum(-6 * mean_ll + 4 * mean_pars_ll)
  return(c("DIC " = DIC, "BPIC" = BPIC, "Effective parameters" = pD))
  
}

pmwg_chainPlots <- function(samples, subjectParPlot = T, parameterPlot = T, subjectLLPlot = T, PDF = F, path = NULL){
  #Creates chain plots per parameter per subject, per parameter and likelihood
  if(PDF) pdf(paste0(path, "/chainPlots.pdf"))
  if (subjectParPlot){
    par(mfrow = c(2, 2))
    for (par in samples$par_names){
      matplot(t(samples$samples$alpha[par,,samples$samples$stage == "sample"]),type="l", main = par, xlab = "samples", ylab = "ParameterValue")
    } 
  }
  par(mfrow=c(1,1))
  if(parameterPlot) matplot(t(samples$samples$theta_mu[, samples$samples$stage == "sample"]), type="l", main = "Paramater chains", ylab = "Parameter Value", xlab = "samples")
  
  if(subjectLLPlot) matplot(t(samples$samples$subj_ll[, samples$samples$stage == "sample"]), type="l", main = "LogLikelihood chains per subject", ylab = "LogLikelihood", xlab = "samples")
  dev.off()
}


pmwg_parHist <- function(samples, PDF = F, path = NULL, stage = "sample"){
  #Creates a histogram of the posterior distribution per parameter
  chains <- as.array(as_mcmc(samples, filter = stage))
  plot <- mcmc_hist(chains)
  if(PDF) pdf(paste0(path, "/parHist.pdf"), onefile = F)
  print(plot)
  if(PDF) dev.off()
  
}

pmwg_parIntervals <- function(samples, PDF = F, path = NULL, stage = "sample"){
  #Creates a interval plot of the posterior distribution per parameter
  chains <- as.array(as_mcmc(samples, filter = stage))
  plot <- mcmc_intervals(chains)
  if(PDF) pdf(paste0(path, "/parIntervals.pdf"), onefile = F)
  print(plot)
  if(PDF)dev.off()
}

pmwg_median <- function(samples){
  #Returns the median of the posterior parameter estimates per subject
  finalSamples <- samples$samples$alpha[,,samples$samples$stage == "sample"]
  medians <- apply(finalSamples, 1:2, median)
  return(medians)
}

pmwg_SD <- function(samples){
  #Returns the SD of the posterior parameter estimates per subject
  finalSamples <- samples$samples$alpha[,,samples$samples$stage == "sample"]
  SDs <- apply(finalSamples, 1:2, sd)
  return(SDs)
}

pmwg_corPlot <- function(samples){
  covcor <- pmwg_covariance(samples)
  print(corrplot(covcor$correlation))
}

pmwg_particlesPerSub <- function(samples){
  #Unique particles per sub 
  x<-apply(samples$samples$alpha[1,,-1]!=samples$samples$alpha[1,,-(samples$samples$idx)],1,sum)
  x[order(x)]
}

pmwg_parValues <- function(samples){
  #Returns array with all parameter estimates, see pmwg_median for the median
  tmp <- samples$samples$alpha[,,samples$samples$idx]
  round(tmp,3)
}


pmwg_covariance <- function(samples, cor = T, plot = T){
  #### Covariance matrix
  cov<-apply(samples$samples$theta_sig[,,samples$samples$idx-1000:samples$samples$idx] ,1:2, mean)
  colnames(cov)<-samples$par_names
  rownames(cov)<-samples$par_names
  if (plot){
    diagonal<-apply(samples$samples$theta_sig,3,diag)
    matplot(log(t(diagonal)), type="l", main = "covariance plot", xlab = "samples", ylab = "covariance")
    #tbh not really sure how to interpret this plot
  }
  
  if (cor) {
    cor<-cov2cor(cov) #correlation matrix
    return(list(correlation = cor, covariance = cov))
  }
}


pmwg_InterSesCor <- function(sampled, sampled2 = NULL){
  if(is.null(sampled2)){
    medPars <- pmwg_median(sampled)
    pars1 <- medPars[1:(nrow(medPars)/2),]
    pars2 <- medPars[((nrow(medPars)/2)+1):nrow(medPars),]
  } else{
    pars1 <- pmwg_median(sampled)
    pars2 <- pmwg_median(sampled2)
  }
  
  par(mfrow = c(3,2))
  mainNames <- gsub(".*[|]", "", rownames(medPars))
  for (i in 1:nrow(pars1)){
    xlim <- ylim <- c(min(min(pars1[i,]), min(pars2[i,])), max(max(pars1[i,]), max(pars2[i,])))
    plot(pars1[i,], pars2[i,], main = mainNames[i], xlab = "Ses1", ylab = "Ses2", xlim = xlim, ylim = ylim)
    tmp <- legend('bottomright', c(" ", " "), bty='n', xjust=1, 
                  text.width = strwidth("RMSE = 0.03"))
    text(tmp$rect$left + tmp$rect$w, tmp$text$y,
         c(paste0('r = ', round(cor(pars1[i,], pars2[i,]), 2)), 
           paste0('RMSE = ', round(rmse(pars1[i,], pars2[i,]), 2))), pos = 2)
    abline(a=0, b=1)
  }
}

pmwg_Simulate <- function(samples, dat, ll_func, preppedData = NULL, n = 5, RL = F, combns = NULL, 
                          addColumns = NULL, match = NULL, jointSamples = NULL, inputPars = samples$par_names,
                          Rev = F){
  samples_stage <- length(samples$samples$stage[samples$samples$stage == "sample"])
  iterations <- round(seq(from = (samples$samples$idx - samples_stage),
                          to = samples$samples$idx,
                          length.out = n))
  
  if(is.null(preppedData)) preppedData <- samples$data
  preppedData <- split(preppedData, preppedData$subject, drop = T)
  
  if(!is.null(jointSamples)){
    parsJoint <- jointSamples$samples$alpha[inputPars,,iterations]
    rownames(parsJoint) <- gsub(".*[|]", "", inputPars)
    pp_joint <- sapply(preppedData, post.predict, parsJoint, ll_func, simplify = F, USE.NAMES = T)
    pp2_joint <- sapply(names(pp_joint), addStimSetInfo, input=pp_joint, orig_dat=dat, 
                  addColumns = addColumns, RL = RL, match = match, Rev = Rev, simplify = F, USE.NAMES = T)
    if(!sfIsRunning()) sfInit(parallel=TRUE, cpus =30); sfLibrary(moments)
    pp_joint <- sapply(pp2_joint, posteriorCalc, combns = combns, simplify = F, USE.NAMES = T)
    
  }
  
  pars <- samples$samples$alpha[,,iterations]
  pp <- sapply(preppedData, post.predict, pars, ll_func, simplify = F, USE.NAMES = T)
  n_bins <- 10 #Only necessary for RL
  pp2 <- sapply(names(pp), addStimSetInfo, input=pp, orig_dat=dat, addColumns = addColumns, RL = RL, match = match, Rev = Rev, simplify = F, USE.NAMES = T)
  data2 <- sapply((names(preppedData)), addStimSetInfo, input=preppedData, orig_dat=dat, addColumns = addColumns, RL = RL, Rev = Rev, match = match, simplify = F, USE.NAMES = T)
  if(!sfIsRunning()) sfInit(parallel=TRUE, cpus =30); sfLibrary(moments)
  pp <- sfSapply(pp2, posteriorCalc, combns = combns, simplify = F, USE.NAMES = T)
  preppedData <- sfSapply(data2, posteriorCalc, combns = combns,  simplify = F, USE.NAMES = T)
  
  if(!is.null(jointSamples)){
    return(list('pp' = pp, 'pp_joint' = pp_joint, 'data' = preppedData))
  }
  
  return(list('pp' = pp, 'data' = preppedData))
}

rmse <- function(x, y) {
  sqrt(mean((x-y)^2))
}
