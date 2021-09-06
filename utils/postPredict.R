source("utils/RLPlots.R")

postCheck <- function(experiments, PP = T, n = 25){
  #Written so that automatically plots are made after running, doesn't work for predictive plots RL yet
  for (exp in experiments){
    path <- file.path("figures", exp$modelName)
    dir.create(path)
    load(paste0("samples/", exp$modelName, ".RData"))
    InformationCriterion(sampled, preppedData = exp$preppedData)
    chainPlots(sampled, subjectParPlot = T, parameterPlot = T, subjectLLPlot = T, PDF = T, path = path)
    parIntervals(sampled, PDF = T, path = path)
    parHist(sampled, PDF = T, path = path)
    
    if(PP){
      dat <- exp$data
      factors = exp$factors
      both <- postSimulate(samples = sampled, dat = dat, ll_func = likelihood.RD,
                           RL = exp$RL, n = n, match = exp$match)
      postCDF(both$pp, both$data, factors, PDF = T, path = path)
      postHist(both$pp, both$data, factors, PDF = T, path = path)
      ppPlot(both$pp, both$data, factors, PDF = T, path = path)
    }
    
  }
}

postSimulate <- function(samples, dat, ll_func, preppedData = NULL, n = 5, RL = F, match = NULL, jointSamples = NULL, inputPars = samples$par_names,
                          Rev = F){
  samples_stageSingle <- length(samples$samples$stage[samples$samples$stage == "sample"])
  iterationsSingle <- round(seq(from = (samples$samples$idx - samples_stageSingle),
                                to = samples$samples$idx,
                                length.out = n))
  
  if(is.null(preppedData)) preppedData <- samples$data
  preppedData <- split(preppedData, preppedData$subject, drop = T)
  
  if(!is.null(jointSamples)){
    samples_stageJoint <- length(jointSamples$samples$stage[jointSamples$samples$stage == "sample"])
    iterationsJoint <- round(seq(from = (jointSamples$samples$idx - samples_stageJoint),
                                 to = jointSamples$samples$idx,
                                 length.out = n))
    parsJoint <- jointSamples$samples$alpha[inputPars,,iterationsJoint]
    rownames(parsJoint) <- gsub(".*[|]", "", inputPars)
    pp_joint <- sapply(preppedData, post.predict, parsJoint, ll_func, simplify = F, USE.NAMES = T)
    pp_joint <- sapply(names(pp_joint), addStimSetInfo, input=pp_joint, orig_dat=dat, RL = RL, match = match, Rev = Rev, simplify = F, USE.NAMES = T)
  }
  
  pars <- samples$samples$alpha[,,iterationsSingle]
  
  pp <- sapply(preppedData, post.predict, pars, ll_func, simplify = F, USE.NAMES = T)
  pp <- sapply(names(pp), addStimSetInfo, input=pp, orig_dat=dat, RL = RL, match = match, 
                Rev = Rev, simplify = F, USE.NAMES = T)
  preppedData <- sapply((names(preppedData)), addStimSetInfo, 
                  input=preppedData, orig_dat=dat, RL = RL, Rev = Rev, match = match, simplify = F, USE.NAMES = T)
  
  if(!is.null(jointSamples)){
    return(list('pp' = pp, 'pp_joint' = pp_joint, 'data' = preppedData))
  }
  
  return(list('pp' = pp, 'data' = preppedData))
}

post.predict <- function(data, pars, ll_func){
  sub <- data$subject[1]
  df <- data.frame()
  for(i in 1:dim(pars)[3]){
    CurPars <- pars[,sub,i]
    tmp <- ll_func(CurPars, data, sample = T)
    tmp$reps <- i
    df <- rbind(df, tmp)
  }
  df$subject <- sub
  return(df)
}

addStimSetInfo <- function(x, input, orig_dat, nBins=10, subNCol = "subjectNumber", RL = F, match, Rev = F) {
  df <- input[[x]]
  
  # make data compatible with pp format
  if(!'reps' %in% colnames(df)) {
    df$reps = 1
  }
  idx = df$reps == 1
  
  if(RL){
    df$acc <- as.integer(df$R == 2)
    # add trialN by stim set, and bin, for first rep
    df$stimSet <- orig_dat[orig_dat[subNCol]==x, 'stimulus_set']
    df$ease <- orig_dat[orig_dat[subNCol]==x, 'ease']
    df$trialNstimSet <- NA
    df$bin <- NA
    if(Rev){
      df$reversed <- orig_dat[orig_dat[subNCol]==x, 'reversed']
      for(stimSet in unique(df[idx, 'stimSet'])) {
        idx2 <- df$stimSet == stimSet
        df[idx & idx2, 'trialNstimSet'][!df[idx & idx2, 'reversed']] <- -(sum(!df[idx & idx2, 'reversed'])-1):0
        df[idx & idx2, 'trialNstimSet'][df[idx & idx2, 'reversed']] <- 1:sum(!df[idx & idx2, 'reversed'])
        df[idx & idx2, 'bin'] <- roundUp(df[idx & idx2, 'trialNstimSet'], 7)/7
      }
    } else{
      for(stimSet in unique(df[idx, 'stimSet'])) {
        idx2 <- df$stimSet == stimSet
        df[idx & idx2, 'trialNstimSet'] <- seq(1, sum(idx & idx2))
        df[idx & idx2, 'bin'] <- as.numeric(cut(seq(1, sum(idx & idx2)), nBins))
      }
    }
    #We use this to copy for other reps (and other factors of interest)
    if(max(df$reps) > 1) {
      for(colName in c('trialNstimSet', 'bin')) {
        df[, colName] <- rep(df[idx, colName], max(df$reps))
      }
    }
  }
  
  # copy&paste for all other reps

  if(!is.na(match[[1]][2])) df$acc <- df$R == df[,match[[1]][2]]
  return(df)
}

postSummary <- function(data, combns, source, quantiles = c(0.1, 0.5, 0.9)){
  df <- dplyr::bind_rows(data)
  dfRTs <-  df %>% group_by(acc, subject, across(combns), reps) %>% 
    summarise(value = quantile(RT, quantiles), quant = quantiles) %>% 
    ungroup() %>% group_by(acc, across(combns), quant, reps) %>% summarise_all(funs(mean))
  
  dfAcc <- df %>% group_by(subject, across(combns), reps) %>% 
    summarise(value = mean(acc), quant = 9) %>% 
    ungroup() %>% group_by(across(combns), quant, reps) %>% summarise_all(funs(mean)) %>% mutate(acc = 9)
  
  df <- rbind(dfAcc, dfRTs) %>% mutate(source = source) %>% select(-subject)
}

ppPlot <- function(pp, data, factors, PDF = F, path = NULL, joint = NULL){
  #Plot utility function for non RL for factors of interest
  for(Factor in factors){
    ppSummary <- postSummary(pp, Factor, 'single')
    dataSummary <- postSummary(data, Factor, 'data')
    df <- rbind(ppSummary, dataSummary) %>% as.data.frame()
    
    ppMeans <- ppSummary %>% group_by_(Factor, "acc", "quant", "source") %>% summarise(value = mean(value))
    means <- rbind(ppMeans, dataSummary %>% dplyr::select(-reps))
    dodgeWdth = .3
    
    if(!is.null(joint)){
      jointSummary <- postSummary(joint, Factor, 'joint')
      jointMeans <- jointSummary %>% group_by_(Factor, "acc", "quant", "source") %>% summarise(value = mean(value))
      df <- rbind(df, jointSummary %>% as.data.frame())
      means <- rbind(means, jointMeans)
      dodgeWdth = .4
    }
    
    nLevels <- length(unique((df[,Factor])))
    
    rtlim <- c(min(df$value[df$acc != 9 ] * 0.9), max((df$value[df$acc != 9] * 1.1)))
    acclim <- c(min(df$value[df$acc == 9] * 0.9), max((df$value[df$acc == 9] * 1.1)))
    p1 <- ggplot(df %>% filter(acc == 9), aes(x = as.factor(get(Factor)), y = value, group = source, color = source)) + 
      geom_point(position=position_dodge(width=dodgeWdth), alpha = .3)  +
      geom_point(data = means %>% filter(acc == 9), size = 2.5, position=position_dodge(width=dodgeWdth))  +
      theme_bw() + 
      labs(x = Factor, y= "Acc") +
      coord_cartesian(ylim = acclim)
    p2 <- ggplot(df %>% filter(acc == 1), aes(x = as.factor(get(Factor)), y = value, group = source, color = source)) + 
      geom_point(position=position_dodge(width=dodgeWdth), alpha = .3)  +
      geom_point(data = means %>% filter(acc == 1), size = 2.5, position=position_dodge(width=dodgeWdth))  +
      theme_bw() + 
      labs(x = Factor, y= "RT Correct") +
      coord_cartesian(ylim =  rtlim)
    p3 <- ggplot(df %>% filter(acc == 0), aes(x = as.factor(get(Factor)), y = value, group = source, color = source)) + 
      geom_point(position=position_dodge(width=dodgeWdth), alpha = .3)  +
      geom_point(data = means %>% filter(acc == 0), size = 2.5, position=position_dodge(width=dodgeWdth))  +
      theme_bw() + 
      labs(x = Factor, y= "RT Error") +
      coord_cartesian(ylim = rtlim)
    if(PDF) png(paste0(path, '/PP_', Factor, '.png'), width=max(3,1*nLevels), height=8, units = "in", res = 300)
    plot <- ggarrange(p1, p2, p3, nrow = 3, common.legend = T)
    print(plot)
    dev.off()
  }
}


postHist <- function(pp, data, factors, PDF = F, path = NULL, joint = NULL){
  for(Factor in factors){
    #For now based on median of the reps, because other ways are more complicated (but better)
    pp <- dplyr::bind_rows(pp) %>% filter(reps == floor(median(reps))) %>% mutate(source = "single")
    data <- dplyr::bind_rows(data) %>% mutate(source = "data")
    df <- rbind(pp, data)
    
    if(!is.null(joint)){
      joint <- dplyr::bind_rows(joint) %>% filter(reps == floor(median(reps))) %>% mutate(source = "joint")
      df <- rbind(df, joint)
    }
    
    df$RT[df$acc == 0] <- -df$RT[df$acc==0]
    
    entries <- list()
    nEntries <- length(unique(df[,Factor]))
    for(i in unique(df[,Factor])){
      plot <- ggplot(df[df[,Factor] == i,], aes(x=RT, fill=source, col = source)) +
        geom_density(alpha=.6) +
        theme_bw() +
        coord_cartesian(xlim = c(-quantile(data$RT, 0.99), quantile(data$RT, 0.99))) + 
        ggtitle(paste0(Factor, ": ", i))
      entries <- c(entries, list(plot))
    }
    if(PDF) png(paste0(path, '/Hist_', Factor, '.png'), height=3*min(nEntries, 4),  width=3*ceiling(nEntries/4), units = "in", res = 300)
    plot <- do.call(ggarrange, c(entries, nrow = min(nEntries, 4), ncol = ceiling(nEntries/4), common.legend = T))
    print(plot)
    dev.off()
  }
}

postCDF <- function(pp, data, factors, PDF = F, path = NULL, joint = NULL){
  quants <- seq(0.1, 0.9, 0.2)
  for(Factor in factors){
    ppSummary <- postSummary(pp, Factor, 'single', quantiles = quants)
    dataSummary <- postSummary(data, Factor, 'data', quantiles = quants)
    df <- rbind(ppSummary, dataSummary) %>% filter(acc != 9)
    
    ppMeanAcc <- dplyr::bind_rows(pp) %>% group_by_(Factor, "reps") %>% summarise(accMean = mean(acc)) %>% mutate(source = "single")
    dataMeanAcc <- dplyr::bind_rows(data) %>% group_by_(Factor, "reps") %>% summarise(accMean = mean(acc)) %>% mutate(source = "data")
    meanAccs <- rbind(ppMeanAcc, dataMeanAcc)
    
    if(!is.null(joint)){
      jointSummary <- postSummary(joint, Factor, 'joint', quantiles = quants)
      jointMeanAcc <- dplyr::bind_rows(joint) %>% group_by_(Factor, "reps") %>% summarise(accMean = mean(acc)) %>% mutate(source = "joint")
      df <- rbind(df, jointSummary)
      meanAccs <- rbind(meanAccs, jointMeanAcc)
    }
    
    df <- inner_join(meanAccs, df, by = c("source", "reps", Factor))
    dfMeans <- df %>% group_by_(Factor, "source", "quant", "acc") %>% summarise(value = mean(value), accMean = mean(accMean)) %>% as.data.frame()
    
    df <- df %>% filter(reps %in% floor(seq(1,max(reps), length.out = 3))) %>% as.data.frame()
    entries <- list()
    nEntries <- length(unique(df[,Factor]))
    for(i in unique(df[,Factor])){
      plot <- ggplot() +
        geom_point(data = df[df[,Factor] == i,] %>% filter(acc == 1), aes(x = value, y = quant*accMean, color = source), alpha = .6) +
        geom_point(data = df[df[,Factor] == i,] %>% filter(acc == 0), aes(x = value, y = quant*(1-accMean), color = source), alpha = .6) +
        geom_line(data = dfMeans[dfMeans[,Factor] == i,] %>% filter(acc == 1), aes(x = value, y = quant*accMean, color = source), size = 0.6) +
        geom_line(data = dfMeans[dfMeans[,Factor] == i,] %>% filter(acc == 0), aes(x = value, y = quant*(1-accMean), color = source), size = 0.6) +
        labs(x = "RT", y= "Defective quantile") +
        theme_bw() +
        ggtitle(paste0(Factor, ": ", i))
      entries <- c(entries, list(plot))
    }
    if(PDF) png(paste0(path, '/CDF_', Factor, '.png'), height=3*min(nEntries, 4),  width=3*ceiling(nEntries/4), units = "in", res = 300)
    plot <- do.call(ggarrange, c(entries, nrow = min(nEntries, 4), ncol = ceiling(nEntries/4), common.legend = T))
    print(plot)
    dev.off()
  }
}

roundUp <- function(x,to=10)
{
  to*(x%/%to + as.logical(x%%to))
}

