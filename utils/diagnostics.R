library(bayesplot)
library(ggpubr)
library(corrplot)
library(snowfall)

# Utility functions for the user ------------------------------------------
# Added by Reilly Innes/Niek Stevenson
InformationCriterion <- function(samples, preppedData = NULL, name = NULL){
  
  # Mean log-likelihood of the overall (samples-stage) model, for each subject
  mean_ll <- apply(samples$samples$subj_ll[, samples$samples$stage == "sample"], 1, mean)
  # Mean of each parameter across iterations.
  # Keep dimensions for parameters and subjects
  mean_pars <- t(apply(samples$samples$alpha[,, samples$samples$stage == "sample"], 1:2, mean))
  
  # Name 'mean_pars' so it can be used by the log_like function
  colnames(mean_pars) <- samples$par_names
  
  # log-likelihood for each subject using their mean parameter vector
  if(is.null(preppedData)) preppedData <- samples$experiments[[name]]$preppedData
  
  data <- preppedData
  mean_pars_ll <- numeric() #sloppy but hey, still fast enough
  for (sub in unique(data$subject)) {
    mean_pars_ll[sub] <- likelihood.RD(mean_pars[sub, ], data = data[data$subject == sub,])
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

chainPlots <- function(samples, subjectParPlot = T, parameterPlot = T, subjectLLPlot = T, PDF = F, path = NULL){
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


parHist <- function(samples, PDF = F, path = NULL, stage = "sample"){
  #Creates a histogram of the posterior distribution per parameter
  samples$samples$alpha[grep("t0", rownames(samples$samples$alpha)),,] <- pnorm(samples$samples$alpha[grep("t0", rownames(samples$samples$alpha)),,])
  samples$samples$alpha[grep("aV", rownames(samples$samples$alpha)),,] <- pnorm(samples$samples$alpha[grep("aV", rownames(samples$samples$alpha)),,])
  
  chains <- as.array(as_mcmc(samples, filter = stage))
  plot <- mcmc_hist(chains)
  if(PDF) pdf(paste0(path, "/parHist.pdf"), onefile = F)
  print(plot)
  if(PDF) dev.off()
  
}

parIntervals <- function(samples, PDF = F, path = NULL, stage = "sample"){
  samples$samples$alpha[grep("t0", rownames(samples$samples$alpha)),,] <- pnorm(samples$samples$alpha[grep("t0", rownames(samples$samples$alpha)),,])
  samples$samples$alpha[grep("aV", rownames(samples$samples$alpha)),,] <- pnorm(samples$samples$alpha[grep("aV", rownames(samples$samples$alpha)),,])
  
  #Creates a interval plot of the posterior distribution per parameter
  chains <- as.array(as_mcmc(samples, filter = stage))
  plot <- mcmc_intervals(chains)
  if(PDF) pdf(paste0(path, "/parIntervals.pdf"), onefile = F)
  print(plot)
  if(PDF)dev.off()
}

parMedian <- function(samples){
  #Returns the median of the posterior parameter estimates per subject
  finalSamples <- samples$samples$alpha[,,samples$samples$stage == "sample"]
  medians <- apply(finalSamples, 1:2, median)
  return(medians)
}

parSD <- function(samples){
  #Returns the SD of the posterior parameter estimates per subject
  finalSamples <- samples$samples$alpha[,,samples$samples$stage == "sample"]
  SDs <- apply(finalSamples, 1:2, sd)
  return(SDs)
}

corPlot <- function(samples){
  covcor <- parCovariance(samples, plot = F)
  print(corrplot(covcor$correlation))
}

parPerSub <- function(samples){
  #Unique particles per sub 
  x<-apply(samples$samples$alpha[1,,-1]!=samples$samples$alpha[1,,-(samples$samples$idx)],1,sum)
  x[order(x)]
}

parValues <- function(samples){
  #Returns array with all parameter estimates, see parMedian for the median
  tmp <- samples$samples$alpha[,,samples$samples$idx]
  round(tmp,3)
}


parCovariance <- function(samples, cor = T, plot = T){
  #### Covariance matrix
  cov<-apply(samples$samples$theta_sig[,,samples$samples$idx-1000:samples$samples$idx] ,1:2, mean)
  
  
  colnames(cov)<-substr(samples$par_names, 1, 22)
  rownames(cov)<-substr(samples$par_names, 1, 22)
  
  
  colnames(cov) <- gsub("MSITprocsB2", "MSIT", colnames(cov))
  colnames(cov) <- gsub("RB4V02Diff2B", "RB", colnames(cov))
  colnames(cov) <- gsub("rlsat2V02B", "RLSAT", colnames(cov))
  colnames(cov) <- gsub("RevStd", "Rev", colnames(cov))
  
  rownames(cov) <- gsub("MSITprocsB2", "MSIT", rownames(cov))
  rownames(cov) <- gsub("RB4V02Diff2B", "RB", rownames(cov))
  rownames(cov) <- gsub("rlsat2V02B", "RLSAT", rownames(cov))
  rownames(cov) <- gsub("RevStd", "Rev", rownames(cov))
  

  
  if (plot){
    diagonal<-apply(samples$samples$theta_sig,3,diag)
    matplot(log(t(diagonal)), type="l", main = "covariance plot", xlab = "samples", ylab = "covariance")  
  }
  
  if (cor) {
    cor<-cov2cor(cov) #correlation matrix
    return(list(correlation = cor, covariance = cov))
  }
}


profilePlot <- function(sampled, n_values = 9, PDF = F, path = NULL){
  pars <- rowMeans(parMedian(sampled))
  ll_func <- sampled$ll_func
  
  test_data <- ll_func(pars, sampled$data, sample = T)
  attr(test_data, "settings") <- attr(sampled$data, "settings")
  sequence <- seq(from = -.22,
                  to = .22,
                  length.out = n_values)
  par(mfrow = c(2,2))
  if(PDF) pdf(paste0(path, "/profilePlots.pdf"))
  for(param in names(pars)){
    testvalues <- sequence + pars[param]
    tmp <- unlist(lapply(testvalues, function (i){ 
      #for each test value, it will apply the function where i = testvalue.
      pars[param] <- i
      ll_func(pars = pars, data = test_data)
    }))
    p <- plot(x = testvalues,
              y = tmp,
              type = "b",
              main = param,
              xlab = "par values",
              ylab = "log-likelihood")
    print(p)
    abline(v = pars[param], col = "blue")
  }
  if(PDF) while (!is.null(dev.list()))  dev.off()
}




interSesCor <- function(sampled, sampled2 = NULL){
  if(is.null(sampled2)){
    medPars <- parMedian(sampled)
    pars1 <- medPars[1:(nrow(medPars)/2),]
    pars2 <- medPars[((nrow(medPars)/2)+1):nrow(medPars),]
  } else{
    pars1 <- parMedian(sampled)
    pars2 <- parMedian(sampled2)
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

chainCertainty <- function(samples, cut = 100, path = NULL){
  jointModelName <- paste(lapply(samples$experiments, function(x){return(x$modelName)}), collapse = "_")
  
  if(is.null(path)) path <- file.path("./figures", "joint", jointModelName)
  dir.create(path, recursive = T)
  pdf(paste0(path, "/chainsHyperVSSubject.pdf"))
  alpha <- samples$samples$alpha
  mu <- samples$samples$theta_mu
  
  alpha[grep("t0", rownames(alpha)),,] <- pnorm(alpha[grep("t0", rownames(alpha)),,])
  alpha[grep("aV", rownames(alpha)),,] <- pnorm(alpha[grep("aV", rownames(alpha)),,])
  alpha <- aperm(alpha)
  
  mu[grep("t0", rownames(mu)),] <- pnorm(mu[grep("t0", rownames(mu)),])
  mu[grep("aV", rownames(mu)),] <- pnorm(mu[grep("aV", rownames(mu)),])
  
  for(i in 1:dim(alpha)[3]){
    tmp <- alpha[cut:nrow(alpha),,i]
    df <- data.frame(iteration = cut:nrow(alpha), median = apply(tmp, 1, median), upper = apply(tmp, 1, quantile, .95), lower = apply(tmp, 1, quantile, .05))
    df$hyper <- mu[i,cut:ncol(mu)]
    plot <- ggplot(df, aes(iteration)) + 
      geom_line(aes(y = median), color = "purple") + 
      geom_line(aes(y = hyper), color = "red") + 
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, fill = "blue") +
      ggtitle(rownames(mu)[i])
    print(plot)
  }
  dev.off()
}


profilePlot <- function(sampled, n_values = 9, PDF = F, path = NULL){
  pars <- rowMeans(parMedian(sampled))
  ll_func <- sampled$ll_func
  
  test_data <- ll_func(pars, sampled$data, sample = T)
  attr(test_data, "settings") <- attr(sampled$data, "settings")
  sequence <- seq(from = -.22,
                  to = .22,
                  length.out = n_values)
  par(mfrow = c(2,2))
  if(PDF) pdf(paste0(path, "/profilePlots.pdf"))
  for(param in names(pars)){
    testvalues <- sequence + pars[param]
    tmp <- unlist(lapply(testvalues, function (i){ 
      #for each test value, it will apply the function where i = testvalue.
      pars[param] <- i
      ll_func(pars = pars, data = test_data)
    }))
    p <- plot(x = testvalues,
              y = tmp,
              type = "b",
              main = param,
              xlab = "par values",
              ylab = "log-likelihood")
    print(p)
    abline(v = pars[param], col = "blue")
  }
  if(PDF) while (!is.null(dev.list()))  dev.off()
}

profilePlot <- function(sampled, n_values = 9, PDF = F, path = NULL){
  pars <- rowMeans(parMedian(sampled))
  ll_func <- sampled$ll_func
  
  test_data <- ll_func(pars, sampled$data, sample = T)
  attr(test_data, "settings") <- attr(sampled$data, "settings")
  sequence <- seq(from = -.22,
                  to = .22,
                  length.out = n_values)
  par(mfrow = c(2,2))
  if(PDF) pdf(paste0(path, "/profilePlots.pdf"))
  for(param in names(pars)){
    testvalues <- sequence + pars[param]
    tmp <- unlist(lapply(testvalues, function (i){ 
      #for each test value, it will apply the function where i = testvalue.
      pars[param] <- i
      ll_func(pars = pars, data = test_data)
    }))
    p <- plot(x = testvalues,
              y = tmp,
              type = "b",
              main = param,
              xlab = "par values",
              ylab = "log-likelihood")
    print(p)
    abline(v = pars[param], col = "blue")
  }
  if(PDF) while (!is.null(dev.list()))  dev.off()
}





rmse <- function(x, y) {
  sqrt(mean((x-y)^2))
}
