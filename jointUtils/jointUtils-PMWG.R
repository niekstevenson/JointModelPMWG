source("jointUtils/jointUtilsPP.R")

pmwg_JointIC <- function(samples, experiments){
  parPreFixs <- gsub("[|].*", "", samples$par_names)
  ICs <- list()
  i <- 0
  for(model in unique(parPreFixs)){
    i <- i + 1
    data <- experiments[[i]]$preppedData
    parIdx <- samples$par_names[which(parPreFixs == model)]
    # Identify number of subjects
    nsubj <- length(unique(samples$data$subject))
    
    # Mean log-likelihood of the overall (samples-stage) model, for each subject
    mean_ll <- apply(samples$samples$subj_ll[, samples$samples$stage == "sample"], 1, mean)
    
    # Mean of each parameter across iterations.
    # Keep dimensions for parameters and subjects
    mean_pars <- t(apply(samples$samples$alpha[parIdx,, samples$samples$stage == "sample"], 1:2, mean))
    
    # Name 'mean_pars' so it can be used by the log_like function
    colnames(mean_pars) <- gsub(".*[|]", "", samples$par_names[which(parPreFixs == model)])
    
    # log-likelihood for each subject using their mean parameter vector
    mean_pars_ll <- numeric(ncol(mean_pars))
    for (sub in unique(data$subject)) {
      mean_pars_ll[sub] <- experiments[[i]]$llFunc(mean_pars[sub, ], data = data[data$subject == sub,])
    }
    
    # mean deviance(-2*ll of all data) 
    # effective number of parameters(-2*ll of all data - -2*ll with mean theta)
    pD <- sum(-2 * mean_ll + 2 * mean_pars_ll)
    
    # DIC = mean deviance + effective number of parameters
    DIC <- sum(-4 * mean_ll + 2 * mean_pars_ll)
    
    # BPIC = mean deviance + 2*effective number of parameters 
    # Note this is the "easy" BPIC, instead of the complex 2007 one
    BPIC <- sum(-6 * mean_ll + 4 * mean_pars_ll)
    ICs[[model]] <- list("DIC " = DIC, "BPIC" = BPIC, "Effective parameters" = pD)
    
  }
  return(ICs)
}


pmwg_jointPostPredict <- function(samples, experiments){
  parPreFixs <- gsub("[|].*", "", samples$par_names)
  i <- 0
  jointModelName <- paste(lapply(experiments, function(x){return(x$modelName)}), collapse = "_")
  for (model in unique(parPreFixs)){
    Rev <- F
    i <- i + 1
    path <- file.path("./figures", "joint", jointModelName, experiments[[i]]$modelName)
    dir.create(path, recursive = T)


    origData <- experiments[[i]]$data
    preppedData <- experiments[[i]]$preppedData
    factors <- experiments[[i]]$factors
    combns <- experiments[[i]]$factors
    
    if(i == 3){
      combns <- lapply(as.list(combns), FUN = function(x) return(c("bin", x)))
    }
    if(i == 4){
      Rev = T
      combns <- "bin"
    }
    llFunc <- experiments[[i]]$llFunc
    currentPars <- samples$par_names[which(parPreFixs == model)]
    load(paste0("samples/", experiments[[i]]$modelName, ".RData"))
    singleSamples <- sampled
    pmwg_jointParInt(samples, singleSamples, pars = currentPars, PDF = T, path = path)
    
    tmp <- pmwg_Simulate(samples = singleSamples, dat = origData, ll_func = llFunc, preppedData = preppedData, RL = attr(preppedData, "RL"), 
                         combns = combns, n = 3, match = attr(preppedData, "match"), jointSamples = samples, inputPars = currentPars, Rev = Rev)
    ppInfo <- getDescriptives(tmp$pp, combns = combns)
    dataInfo <- getDescriptives(tmp$data, combns = combns)
    jointInfo <- getDescriptives(tmp$pp_joint, combns = combns)
    if(i == 1 | i == 2){
      pmwg_jointPpPlot(ppInfo, dataInfo, jointInfo, factors, PDF = T, path = path)
    }
    if(i == 3){
      debug(pmwg_jointRLSATPlot)
      pmwg_jointRLSATPlot(ppInfo, dataInfo, jointInfo, PDF = T, path = path)
    }
    if(i == 4){
      pmwg_jointRevPlot(ppInfo, dataInfo, jointInfo, PDF = T, path = path)
    }
  }
}



pmwg_jointPpPlot <- function(ppInfo, dataInfo, jointInfo, factors, PDF = F, path = NULL){
  OneDF <- function(df){
    newdf <- data.frame()
    for (j in c("Correct", "Error", "Acc")){
      if (j != "Acc"){
        for(rt in c("RT.10.", "RT.50.", "RT.90.")){
          newdf <- rbind(newdf, df[[paste0(rt, j, "By", Factor)]] %>% rename(column = rt) %>% mutate(Class = j, quant = rt))
        }
      } else{
        newdf <- rbind(newdf, df[[paste0(j, "By", Factor)]] %>% rename(column = "acc") %>% mutate(Class = j, quant = "Acc"))
      }
    }
    return(newdf)
  }
  
  #Plot utility function for non RL for factors of interest
  for(Factor in factors){
    ppdf <- OneDF(ppInfo)
    ppdf$source <- "Single"
    ppdfMeans <- aggregate(column ~get(Factor)*Class*source*quant, ppdf, mean)
    colnames(ppdfMeans)[1] <- Factor
    
    ppdf_joint <- OneDF(jointInfo)
    ppdf_joint$source <- "Joint"
    ppdf_jointMeans <- aggregate(column ~get(Factor)*Class*source*quant, ppdf_joint, mean)
    colnames(ppdf_jointMeans)[1] <- Factor
    
    datadf <- OneDF(dataInfo)
    datadf$source <- "data"
    
    df <- rbind(ppdf, datadf, ppdf_joint)
    
    datadf <- datadf %>% select(-reps)
    meandf <- rbind(ppdfMeans, datadf, ppdf_jointMeans)
    
    rtlim <- c(min(df$column[df$Class == "Correct" | df$Class == "Error" ] * 0.9), 
               max((df$column[df$Class == "Correct" | df$Class == "Error" ] * 1.1)))
    acclim <- c(min(df$column[df$Class == "Acc"] * 0.9), max((df$column[df$Class == "Acc"] * 1.1)))
    
    p1 <- ggplot(df %>% filter(Class == "Acc"), aes(x = as.factor(get(Factor)), y = column, group = source, color = source)) + 
      geom_point(position=position_dodge(width=0.4), alpha = .3)  +
      geom_point(data = meandf %>% filter(Class == "Acc"), size = 2.5, position=position_dodge(width=0.4))  +
      theme_bw() + 
      labs(x = Factor, y= "Acc") +
      coord_cartesian(ylim = acclim)
    p2 <- ggplot(df %>% filter(Class == "Correct"), aes(x = as.factor(get(Factor)), y = column, group = source, color = source)) + 
      geom_point(position=position_dodge(width=0.4), alpha = .3)  +
      geom_point(data = meandf %>% filter(Class == "Correct"), size = 2.5, position=position_dodge(width=0.4))  +
      theme_bw() + 
      labs(x = Factor, y= "RT Correct") +
      coord_cartesian(ylim =  rtlim)
    p3 <- ggplot(df %>% filter(Class == "Error"), aes(x = as.factor(get(Factor)), y = column, group = source, color = source)) + 
      geom_point(position=position_dodge(width=0.4), alpha = .3)  +
      geom_point(data = meandf %>% filter(Class == "Error"), size = 2.5, position=position_dodge(width=0.4))  +
      theme_bw() + 
      labs(x = Factor, y= "RT Error") +
      coord_cartesian(ylim = rtlim)
    if(PDF) png(paste0(path, '/PP_', Factor, '.png'), width=3, height=8, units = "in", res = 300)
    plot <- ggarrange(p1, p2, p3, nrow = 3, common.legend = T)
    print(plot)
    dev.off()
  }
}

pmwg_jointParInt <- function(joint, single, pars, PDF = F, path = NULL){
  single <- as.array(as_mcmc(single, filter = "sample"))
  joint <- as.array(as_mcmc(joint, filter = "sample"))
  joint <- joint[, pars]
  colnames(joint) <- colnames(single)
  
  # use bayesplot::mcmc_intervals_data() function to get intervals data in format easy to pass to ggplot
  combined <- rbind(mcmc_intervals_data(single), mcmc_intervals_data(joint))
  combined$model <- rep(c("Single", "Joint"), each = ncol(single))
  
  theme_set(bayesplot::theme_default())
  pos <- position_nudge(y = ifelse(combined$model == "Joint", 0, 0.15))
  intervals <- ggplot(combined, aes(x = m, y = parameter, color = model)) + 
    geom_linerange(aes(xmin = l, xmax = h), position = pos, size=2)+
    geom_linerange(aes(xmin = ll, xmax = hh), position = pos)+
    geom_point(position = pos)+
    coord_flip()+
    geom_vline(xintercept=0,linetype="dashed") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  if(PDF) png(paste0(path, '/jointParInt.png'), width = 5, height = 3.5, units = "in", res = 300)
  print(intervals)
  dev.off()
}

