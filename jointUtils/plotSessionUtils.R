#So these are not necessarily well made, but just utility functions to plot some statistics over sessions/stimulus-sets (for RL)/bins


pmwg_binRevPlot <- function(rev){
  #This function makes a reversal plot for every stimulus set for the RL-Rev. This way we can check whether they learn from the reversals.
  roundUp <- function(x,to=10)
  {
    to*(x%/%to + as.logical(x%%to))
  }
  data <- rev$Rev$data
  data$correct <- data$correctTotal
  data$correct <- data$correct - 1

  newData <- data.frame()
  rtData <- data.frame()
  accData <- data.frame()
  for(sub in unique(data$pp)){
    inclSub <- T
    for(ses in c(1:2)){
      tmp <- data[data$pp == sub & data$sessionNr == ses,]
      for (setNr in unique(tmp$stimulus_set)){
        tmp2 <- tmp[tmp$stimulus_set == setNr,]
        if(!nrow(tmp2) < 1){
          if(!((mean(tmp2$reversed) < 0.2) | (mean(tmp2$reversed) > 0.8))){
            tmp2$trialNRev[!tmp2$reversed] <- -(nrow(tmp2[!tmp2$reversed,]) -1):1
            tmp2$trialNRev[tmp2$reversed] <- 1:nrow(tmp2[tmp2$reversed,])
            tmp2$trialNRev <- tmp2$trialNRev
            tmp2$trialNRev <- roundUp(tmp2$trialNRev, 5)
            tmp2$ses <- ses
            tmp2 <- tmp2[,c("pp", "rt", "correct", "trialNRev", "stimulus_set", "ses")]
            newData <- rbind(newData, tmp2)
          }
        } else{
          inclSub <- F
        }
      }
    }
    if(inclSub){
      #lazy programming
      tmp <- newData[newData$pp == sub,]
      tmpRT <- do.call(data.frame, aggregate(rt~correct*pp*trialNRev*stimulus_set*ses, tmp, mean))
      tmpAcc <- aggregate(correct~pp*trialNRev*stimulus_set*ses, tmp, mean)
      rtData <- rbind(rtData, tmpRT)
      accData <- rbind(accData, tmpAcc)
    }
    
  }
  rtData <- rtData[rtData$trialNRev < 31,]
  rtData <- rtData[rtData$trialNRev > - 31,]
  accData <- accData[accData$trialNRev < 31,]
  accData <- accData[accData$trialNRev > - 31,]
  rt <- do.call(data.frame, aggregate(rt ~ trialNRev*stimulus_set*ses, rtData, FUN = mean))
  rt$trialNRev <- rt$trialNRev + 100*rt$stimulus_set
  acc <- do.call(data.frame, aggregate(correct ~ trialNRev*stimulus_set*ses, accData, FUN = mean))
  acc$trialNRev <- acc$trialNRev + 100*acc$stimulus_set
  colnames(rt)[c(1,4)] <- c("var","mean")
  colnames(acc)[c(1,4)] <- c("var","mean")
  par(mfcol = c(1,2))
  accLim <- c(0.2,0.9)
  rtLim <- c(min(rt$mean)*0.9, max(rt$mean)*1.1)
  for (j in 0:1){
    tmp <- rt
    ylab <- "RT"
    ylim <- rtLim
    if(j == 0){tmp <- acc; ylab <- "Accuracy"; ylim <- accLim}
    plot(0,0, type = "n", ylim = ylim, xlim = c(min(rt$var), max(rt$var)), xlab = "set",xaxt="n", ylab = ylab, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    for(i in unique(acc$stimulus_set)){
      tmp2 <- tmp[tmp$stimulus_set == i,]
      lines(tmp2$var, tmp2$mean, lwd = 2)
      if(j == 1) abline(v = 100*tmp$stimulus_set, lwd = 1, lty = 2, col = "red")
    }
  }

}



pmwg_binPlot <- function(experiments, n_bins = 10){
  #This function makes binplots accross the two sessions, so we can see speed up/increase in accuracy between sessions
  for (exp in experiments){
    
    data <- exp$data
    if (exp$name == "rlsat" | exp$name == "Rev"){
      if(exp$name == "Rev") data$correct <- data$correctTotal
      data$correct <- data$correct - 1
    }
    data$bin <- NA
    newData <- data.frame()
    rtData <- data.frame()
    accData <- data.frame()
    for(sub in unique(data$pp)){
      inclSub <- T
      for(ses in c(1:2)){
        tmp <- data[data$pp == sub & data$sessionNr == ses,]
        if(!nrow(tmp) < 1){
          tmp$bin <- as.numeric(cut(1:nrow(tmp), n_bins)) + (ses-1)*n_bins
          tmp <- tmp[,c("pp", "rt", "correct", "bin")]
          newData <- rbind(newData, tmp)
        } else{
          inclSub <- F
        }
        
      }
      if(inclSub){
        #lazy programming
        tmp <- newData[newData$pp == sub,]
        tmpRT <- do.call(data.frame, aggregate(rt~correct*pp*bin, tmp, quantile, probs=seq(.1, .9, .4)))
        tmpAcc <- aggregate(correct~pp*bin, tmp, mean)
        rtData <- rbind(rtData, tmpRT)
        accData <- rbind(accData, tmpAcc)
      }
      
    }
    makeBinPlot(rtData, accData, "bin", xaxs = c(0, 20), xaxstick = seq(2, 20, 2))
  }
  
}

pmwg_sesComparePlot <- function(experiments){
  #This function plots mean acc/rt for each task depending on where they were presented in the task order, for both the 1st and 2nd session.
  for (exp in experiments){
    data <- exp$data
    if (exp$name == "rlsat" | exp$name == "Rev"){
      if(exp$name == "Rev") data$correct <- data$correctTotal
      data$correct <- data$correct - 1
    }
    newData <- data.frame()
    rtData <- data.frame()
    accData <- data.frame()
    for(sub in unique(data$pp)){
      inclSub <- T
      for(ses in c(1:2)){
        tmp <- data[data$pp == sub & data$sessionNr == ses,]
        if(!nrow(tmp) < 1){
          if(ses == 1){
            tmp$taskNr <- min(which(strsplit(tmp$sesOrder[1], "_")[[1]] == exp$name))
          } else{
            tmp$taskNr <- max(which(strsplit(tmp$sesOrder[1], "_")[[1]] == exp$name))
          }
          tmp <- tmp[,c("pp", "rt", "correct", "taskNr")]
          newData <- rbind(newData, tmp)
        } else{
          inclSub <- F
        }
      }
      if(inclSub){
        #lazy programming
        tmp <- newData[newData$pp == sub,]
        tmpRT <- do.call(data.frame, aggregate(rt~correct*pp*taskNr, tmp, quantile, probs=seq(.1, .9, .4)))
        tmpAcc <- aggregate(correct~pp*taskNr, tmp, mean)
        rtData <- rbind(rtData, tmpRT)
        accData <- rbind(accData, tmpAcc)
      }
    }
    makeBinPlot(rtData, accData, "taskNr", xaxs = c(0, 10), xaxstick = seq(2, 10, 2))
  }
}

pmwg_withinSesComparePlot <- function(experiments){
  #This function makes a plot that compares whether it mattered if the current task was presented as first or second in it's respective session e.g.
  #if session one = SSTMSIT then MSIT is 2nd. Third corresponds to first in the second session and fourth to second in the second session. 
  for (exp in experiments){
    data <- exp$data
    if (exp$name == "rlsat" | exp$name == "Rev"){
      if(exp$name == "Rev") {
        data$correct <- data$correctTotal; 
        coupledName <- "RB";
      }
      data$correct <- data$correct - 1
    }
    if(exp$name == "MSIT") coupledName <- "SST"
    if(exp$name == "RB") coupledName <- "Rev"
    newData <- data.frame()
    rtData <- data.frame()
    accData <- data.frame()
    for(sub in unique(data$pp)){
      inclSub <- T
      for(ses in c(1:2)){
        tmp <- data[data$pp == sub & data$sessionNr == ses,]
        if(!nrow(tmp) < 1){
          tasks <- strsplit(tmp$sesOrder[1], "_")[[1]]
          if(ses == 1){
            idx <- min(which(tasks == exp$name))
            tmp$withinSesNr <- 1
            if(idx > 1){
              if(tasks[idx -1] == coupledName) tmp$withinSesNr <- 2
            }
          } else{
            idx <- max(which(tasks == exp$name))
            tmp$withinSesNr <- 3
            if(tasks[idx -1] == coupledName){
              if(tasks[idx -2] == exp$name){
                if(tasks[idx -3] == coupledName){
                  tmp$withinSesNr <- 4
                }
              } else{
                tmp$withinSesNr <- 4
              }
            }
          }
          tmp <- tmp[,c("pp", "rt", "correct", "withinSesNr")]
          newData <- rbind(newData, tmp)
        } else{
          inclSub <- F
        }
      }
      if(inclSub){
        #lazy programming
        tmp <- newData[newData$pp == sub,]
        tmpRT <- do.call(data.frame, aggregate(rt~correct*pp*withinSesNr, tmp, quantile, probs=seq(.1, .9, .4)))
        tmpAcc <- aggregate(correct~pp*withinSesNr, tmp, mean)
        rtData <- rbind(rtData, tmpRT)
        accData <- rbind(accData, tmpAcc)
      }
    }
    makeBinPlot(rtData, accData, "withinSesNr", xaxs = c(0, 4), xaxstick = seq(1, 4, 1))
  }  
}

makeBinPlot <- function(rtData, accData, var, xaxs, xaxstick){
  #General utility function for plotting data over time frames on the x-axis. 
  std <- function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)))
  rt10 <- do.call(data.frame, aggregate(rt.10. ~ get(var)*correct, rtData, FUN = std))
  rt50 <- do.call(data.frame, aggregate(rt.50. ~ get(var)*correct, rtData, FUN = std))
  rt90 <- do.call(data.frame, aggregate(rt.90. ~ get(var)*correct, rtData, FUN = std))
  acc <- do.call(data.frame, aggregate(correct ~ get(var), accData, FUN = std))
  colnames(rt10)[c(1,3,4)] <- c("var","mean", "SE")
  colnames(rt50)[c(1,3,4)] <- c("var","mean", "SE")
  colnames(rt90)[c(1,3,4)] <- c("var","mean", "SE")
  colnames(acc)[1:3] <- c("var","mean", "SE")
  par(mfcol = c(1,3))
  accLim <- as.numeric(quantile(accData$correct, probs = c(0.05, 0.95)))
  rtLim <- c(min(rtData$rt.10.[rtData$correct == 1])*1.1, max(rtData$rt.90.[rtData$correct == 1])*0.8)
  plot(acc$var, acc$mean, type = "l", ylim = accLim, xlab = var, ylab = "Accuracy",
       cex.lab=1.3, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, xaxt = "n")
  arrows(x0=acc$var, y0=acc$mean-acc$SE, x1=acc$var, y1=acc$mean+acc$SE, code=3, angle=90, length=0.1)
  axis(1, at=xaxstick, lwd=1.5)
  for(i in c(1,0)){
    ylab <- "Error RT"
    if(i == 1) ylab <- "Correct RT"
    plot(0,0, type = "n", ylim = rtLim, xlim = xaxs, xlab = var, xaxt = "n", ylab = ylab, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    tmp10 <- rt10[rt10$correct == i,]
    tmp50 <- rt50[rt50$correct == i,]
    tmp90 <- rt90[rt50$correct == i,]
    lines(tmp10$var, tmp10$mean)
    arrows(x0=tmp10$var, y0=tmp10$mean-tmp10$SE, x1=tmp10$var, y1=tmp10$mean+tmp10$SE, code=3, angle=90, length=0.1)
    lines(tmp50$var, tmp50$mean)
    arrows(x0=tmp50$var, y0=tmp50$mean-tmp50$SE, x1=tmp50$var, y1=tmp50$mean+tmp50$SE, code=3, angle=90, length=0.1)
    lines(tmp90$var, tmp90$mean)
    arrows(x0=tmp90$var, y0=tmp90$mean-tmp90$SE, x1=tmp90$var, y1=tmp90$mean+tmp90$SE, code=3, angle=90, length=0.1)
    axis(1, at=xaxstick, lwd=1.5)
  }
}


