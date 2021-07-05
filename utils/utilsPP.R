
# Some functions used by the utility functions ----------------------------
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

addStimSetInfo <- function(x, input, orig_dat, addColumns=NULL, nBins=10, subNCol = "pp", RL = F, match, Rev = F) {
  df <- input[[x]]
  
  # make data compatible with pp format
  if(!'reps' %in% colnames(df)) {
    df$reps = 1
  }
  idx = df$reps == 1
  # add other columns from original data if wanted/needed
  if(!is.null(addColumns)) {
    for(colName in addColumns) {
      df[idx, colName] <- orig_dat[orig_dat[,subNCol]==x, colName]
    }
  }
  
  
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
    addColumns <- c(addColumns, 'trialNstimSet', 'bin')
  }
  
  # copy&paste for all other reps
  if(max(df$reps) > 1) {
    for(colName in addColumns) {
      df[, colName] <- rep(df[idx, colName], max(df$reps))
    }
  }
  if(!is.na(match[[1]][2])) df$acc <- df$R == df[,match[[1]][2]]
  return(df)
}

pmwg_ppPlot <- function(ppInfo, dataInfo, factors, PDF = F, path = NULL){
  #Plot utility function for non RL for factors of interest
  for(Factor in factors){
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
    
    ppdf <- OneDF(ppInfo)
    ppdf$source <- "model"
    ppdfMeans <- aggregate(column ~get(Factor)*Class*source*quant, ppdf, mean)
    colnames(ppdfMeans)[1] <- Factor
    datadf <- OneDF(dataInfo)
    datadf$source <- "data"
    
    df <- rbind(ppdf, datadf)
    
    datadf <- datadf %>% select(-reps)
    meandf <- rbind(ppdfMeans, datadf)
    
    rtlim <- c(min(df$column[df$Class == "Correct" | df$Class == "Error" ] * 0.9), 
               max((df$column[df$Class == "Correct" | df$Class == "Error" ] * 1.1)))
    acclim <- c(min(df$column[df$Class == "Acc"] * 0.9), max((df$column[df$Class == "Acc"] * 1.1)))
    
    p1 <- ggplot(df %>% filter(Class == "Acc"), aes(x = as.factor(get(Factor)), y = column, group = source, color = source)) + 
      geom_point(position=position_dodge(width=0.3), alpha = .3)  +
      geom_point(data = meandf %>% filter(Class == "Acc"), size = 2.5, position=position_dodge(width=0.3))  +
      theme_bw() + 
      labs(x = Factor, y= "Acc") +
      coord_cartesian(ylim = acclim)
    p2 <- ggplot(df %>% filter(Class == "Correct"), aes(x = as.factor(get(Factor)), y = column, group = source, color = source)) + 
      geom_point(position=position_dodge(width=0.3), alpha = .3)  +
      geom_point(data = meandf %>% filter(Class == "Correct"), size = 2.5, position=position_dodge(width=0.3))  +
      theme_bw() + 
      labs(x = Factor, y= "RT Correct") +
      coord_cartesian(ylim =  rtlim)
    p3 <- ggplot(df %>% filter(Class == "Error"), aes(x = as.factor(get(Factor)), y = column, group = source, color = source)) + 
      geom_point(position=position_dodge(width=0.3), alpha = .3)  +
      geom_point(data = meandf %>% filter(Class == "Error"), size = 2.5, position=position_dodge(width=0.3))  +
      theme_bw() + 
      labs(x = Factor, y= "RT Error") +
      coord_cartesian(ylim = rtlim)
    if(PDF) pdf(paste0(path, '/PP_', Factor, '.pdf'), width=3, height=8, onefile = F)
    plot <- ggarrange(p1, p2, p3, nrow = 3, common.legend = T)
    print(plot)
    dev.off()
  }
}


pmwg_RLSATPlot <- function(ppInfo, dataInfo, reversal = F, PDF = F, path = NULL){
  
  # Plot --------------------------------------------------------------------
  layoutM <- matrix(1:25, nrow=5, byrow=TRUE)
  layoutM[c(1, 5),1:2] <- 1
  layoutM[c(1, 5),4:5] <- 9
  layoutM[2:4,1:2] <- 2:7 #matrix(c(2:13), nrow=3, byrow=TRUE)
  layoutM[2:4,4:5] <- 10:15 #matrix(c(2:13), nrow=3, byrow=TRUE)
  layoutM[,3] <- 8
  layoutM
  if(PDF) pdf(paste0(path, '/BinPlotPerCue.pdf'), width=7, height=7/4*3)
  layout(layoutM, heights = c(0.01, .8, 1, 1, 0.01), widths=c(1,1,.1,1,1))
  par(oma=c(3,4,2,0), mar=c(0, 0, 1, 0.5) + 0.1, #mfcol=c(3,4), 
      mgp=c(2.75,.75,0), las=1, bty='l')
  i <- 0
  data.cex=1.5
  corrRTylim <- errRTylim <- c(.35,1.1)
  
  plot.new()
  if(i == 0) mtext('RL-ARD', side=3, cex=.66*1.2, font=2, line=1)
  for(cue in c('SPD', 'ACC')) {
    i <- i+1
    idxD <- dataInfo$AccBybinBycue$cue==cue
    idxM <- ppInfo$AccBybinBycue$cue==cue
    
    plotDataPPBins(data=dataInfo$AccBybinBycue[idxD,], pp=ppInfo$AccBybinBycue[idxM,],
                   xaxt='n', draw.legend = i==1, data.cex = data.cex,
                   dep.var='acc', ylab='', xlab = '', yaxt='n',
                   legend.pos='bottomright', ylim=c(0.5, 0.9), hline.by=0.1)
    axis(1, at=seq(2, 10, 2), labels=rep(NA, 5), lwd=2)
    if(i == 1) {
      mtext('Accuracy', side=2, cex=.66, line=3, las=0, font=1)
      axis(2, at=seq(.5, .9, .1), lwd=1.5)
    } else {
      axis(2, at=seq(.5, .9, .1), labels=rep(NA, 5), lwd=1.5)
    }
    if(i == 1) title('Speed')
    if(i == 2) title('Accuracy')
    if(i == 3) title('Speed')
    if(i == 4) title('Accuracy')
    
    ##
    plotDataPPBins(data=dataInfo$RT.10.CorrectBybinBycue[idxD,], pp=ppInfo$RT.10.CorrectBybinBycue[idxM,], dep.var='RT.10.', 
                   ylim=corrRTylim, xaxt='n', ylab='', yaxt='n', draw.legend = FALSE, data.cex = data.cex)
    plotDataPPBins(data=dataInfo$RT.50.CorrectBybinBycue[idxD,], pp=ppInfo$RT.50.CorrectBybinBycue[idxM,], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
    plotDataPPBins(data=dataInfo$RT.90.CorrectBybinBycue[idxD,], pp=ppInfo$RT.90.CorrectBybinBycue[idxM,], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
    axis(1, at=seq(2, 10, 2), labels=NA, lwd=1.5)
    if(i == 1) {
      mtext('Correct RTs (s)', side=2, cex=.66, line=3, las=0, font=1)
      axis(2, seq(.2, 1.6, .2), lwd=1.5)
    } else {
      axis(2, seq(.2, 1.6, .2), labels=NA, lwd=1.5)
    }
    
    ##
    plotDataPPBins(data=dataInfo$RT.10.ErrorBybinBycue[idxD,], pp=ppInfo$RT.10.ErrorBybinBycue[idxM,], dep.var='RT.10.', ylim=errRTylim, xaxt='n', ylab='', yaxt='n', draw.legend = FALSE, data.cex = data.cex)
    plotDataPPBins(data=dataInfo$RT.50.ErrorBybinBycue[idxD,], pp=ppInfo$RT.50.ErrorBybinBycue[idxM,], dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
    plotDataPPBins(data=dataInfo$RT.90.ErrorBybinBycue[idxD,], pp=ppInfo$RT.90.ErrorBybinBycue[idxM,], dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
    if(i == 1) {
      mtext('Error RTs (s)', side=2, cex=.66, line=3, las=0, font=1)
      axis(2, seq(.4, 1.2, .2), lwd=1.5)
    } else {
      axis(2, seq(.4, 1.2, .2), labels=NA, lwd=1.5)
    }
    axis(1, at=seq(2, 10, 2), lwd=1.5)
    mtext('Trial bin', side=1, cex=.66, line=2)
    
  }
  
  dev.off()
}


plotDataPPBins <- function(data, pp, joint_pp = NULL, dep.var='RT', xaxis='bin', colorM='blue', colorJ= 'green', colorD=1, 
                           draw.legend=TRUE, legend.pos='topright', plot.new=TRUE, draw.polygon=TRUE,
                           xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, data.cex=2, data.lwd=2,
                           hline.by=0.1, axhlines=NULL, plot.model.points=TRUE,
                           vline.by=2, axvlines=NULL, xaxt='s', yaxt='s') {
  #Utility function for RL plots taken from Steven, still working on some nicer automatic RL plots per factor.
  if(is.null(xlim))  xlim = range(data[,xaxis])+c(-.5, .5)
  if(is.null(ylim)) {
    ylimD = range(data[,dep.var])*c(.9, 1.1)
    ylimM = range(pp[,dep.var])*c(.9, 1.1)
    ylim = c(min(ylimD[1], ylimM[1]), max(ylimD[2], ylimM[2]))
  }
  if(is.null(xlab)) xlab <- xaxis
  if(is.null(ylab)) ylab <- dep.var
  
  # empty canvas
  if(plot.new) {
    plot(0,0, type='n', xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, xaxt=xaxt, yaxt=yaxt)
    if(!is.null(axhlines)) {
      abline(h=axhlines, col='lightgrey', lty=1)
    } else {
      abline(h=seq(0, 5, hline.by), col='lightgrey', lty=1)
    }
    if(!is.null(axvlines)) {
      abline(v=axvlines, col='lightgrey', lty=1)
    } else {
      abline(v=seq(-20, 20, vline.by), col='lightgrey', lty=1)
    }
  }
  
  if(draw.polygon) {
    lowerQ <- aggregate(as.formula(paste0(dep.var, '~bin')), pp, quantile, .025)
    upperQ <- aggregate(as.formula(paste0(dep.var, '~bin')), pp, quantile, .975)
    xs <- c(lowerQ[,xaxis], rev(lowerQ[,xaxis]))
    ys <- c(lowerQ[,dep.var], rev(upperQ[,dep.var]))
    polygon(xs, ys, col=rgb(col2rgb(colorM)[1]/255, col2rgb(colorM)[2]/255, col2rgb(colorM)[3]/255, alpha=.3), lty = NULL, border=NA)
  }
  
  if(!is.null(joint_pp)) {
    lowerQ <- aggregate(as.formula(paste0(dep.var, '~bin')), joint_pp, quantile, .025)
    upperQ <- aggregate(as.formula(paste0(dep.var, '~bin')), joint_pp, quantile, .975)
    xs <- c(lowerQ[,xaxis], rev(lowerQ[,xaxis]))
    ys <- c(lowerQ[,dep.var], rev(upperQ[,dep.var]))
    polygon(xs, ys, col=rgb(col2rgb(colorJ)[1]/255, col2rgb(colorJ)[2]/255, col2rgb(colorJ)[3]/255, alpha=.3), lty = NULL, border=NA)
  }
  
  # model
  if(plot.model.points) points(pp[,xaxis], pp[,dep.var], pch=20, col=colorM, cex=.01)
  
  # data
  points(data[,xaxis], data[,dep.var], col=colorD, pch=20, cex=data.cex)
  lines(data[,xaxis], data[,dep.var], col=colorD, lwd=data.lwd)
  if(!draw.polygon) {
    if(draw.legend) legend(legend.pos, c('Data', 'Model'), lty=c(1, NA), lwd=c(2, 2), col=c(colorD,colorM), pch=c(20, 20), bty='n')
  } else {
    if(draw.legend) {
      if(!is.null(joint_pp)){
        legend(legend.pos, c('Data', 'Model', 'Joint'), lty=c(1, NA, NA), 
               fill=c(NA, rgb(col2rgb(colorM)[1]/255, col2rgb(colorM)[2]/255, col2rgb(colorM)[3]/255, alpha=.3), 
                      rgb(col2rgb(colorJ)[1]/255, col2rgb(colorJ)[2]/255, col2rgb(colorJ)[3]/255, alpha=.3)),
               border=c(NA, 'black', 'black'),
               lwd=c(2, NA, NA), col=c(colorD, colorM, colorJ), pch=c(20, NA, NA), bty='n')
      }
      else{
      legend(legend.pos, c('Data', 'Model'), lty=c(1, NA), 
             fill=c(NA, rgb(col2rgb(colorM)[1]/255, col2rgb(colorM)[2]/255, col2rgb(colorM)[3]/255, alpha=.3)),
             border=c(NA, 'black'),
             lwd=c(2, NA), col=c(colorD,colorM), pch=c(20, NA), bty='n')
      }
    }
  }
}

posteriorCalc <- function(df, combns){
  if(!'acc' %in% colnames(df)) stop('Accuracy column is missing in data')
  for (comb in combns){
    #Hideous but hey, that's ok! Calculate accuracy, tmp90s, RT50s and RT90s for all factors per subject
    attr(df, paste0("qRTsBy", paste(comb, collapse = "By"))) <- do.call(data.frame, 
                                                                       aggregate(as.formula(paste0("RT~reps*", paste(comb, collapse = "*"))), 
                                                                                 df, quantile, probs=seq(.1, .9, .4)))
    attr(df, paste0("qRTsCorrectBy", paste(comb, collapse = "By"))) <- do.call(data.frame, 
                                                                              aggregate(as.formula(paste0("RT~reps*", paste(comb, collapse = "*"))), 
                                                                                        df[df$acc==1,], quantile, probs=seq(.1, .9, .4)))
    attr(df, paste0("qRTsErrorBy", paste(comb, collapse = "By"))) <- do.call(data.frame, 
                                                                            aggregate(as.formula(paste0("RT~reps*", paste(comb, collapse = "*"))), 
                                                                                      df[df$acc==0,], quantile, probs=seq(.1, .9, .4)))
    attr(df, paste0("AccBy", paste(comb, collapse = "By"))) <- do.call(data.frame, 
                                                                      aggregate(as.formula(paste0("acc~reps*", paste(comb, collapse = "*"))), 
                                                                                df, mean))
  }
  return(df)
}

calcDescriptives <- function(x, dep.var='RT', attr.name='RTsOverBins', id.var1='~reps*bin*s', id.var2='~reps*bin') {
  #Get the subject Acc, tmp90s, RT50s, RT90s and calculate group-level for this specific factor
  #First per subject
  allOverTime <- do.call(rbind, (lapply(names(x), function(y) {tmp <- attr(x[[y]], attr.name); tmp$s <- y; tmp})))
  if ("reps" %in% colnames(allOverTime)){
    id.var1 <- paste0("reps*", id.var1)
    if(!is.null(id.var2)){
      id.var2 <- paste0("reps*", id.var2)
    }
  } 
  
  form1 <- as.formula(paste0(dep.var, "~", id.var1))
  res <- aggregate(form1, allOverTime, mean)
  #Then group-level
  if(!is.null(id.var2)) {
    form2 <- as.formula(paste0(dep.var, "~", id.var2))
    res2 <- aggregate(form2, aggregate(form1, res, mean), mean)
  }
  return(res2)
}

getDescriptives <- function(df, combns){
  #Get the subject Acc, tmp90s, RT50s, RT90s and calculate group-level for all factors.
  descriptives <- list()
  for (comb in combns){
    for(descr in c("", "Correct", "Error")){
      #At this point send help
      for(quant in c("RT.10.", "RT.50.", "RT.90.")){
        descriptives[[paste0(quant, descr, "By", paste(comb, collapse = "By"))]] <- calcDescriptives(df, dep.var = quant, 
                                                                                                    attr.name = paste0("qRTs", descr, "By", paste(comb, collapse = "By")),
                                                                                                    id.var1 = paste(c(comb, "s"), collapse = "*"),
                                                                                                    id.var2 = paste(comb, collapse = "*"))     
      }
    }
    descriptives[[paste0("AccBy", paste(comb, collapse = "By"))]] <- calcDescriptives(df, dep.var = "acc", 
                                                                                     attr.name = paste0("AccBy", paste(comb, collapse = "By")),
                                                                                     id.var1 = paste(c(comb, "s"), collapse = "*"),
                                                                                     id.var2 = paste(comb, collapse = "*"))   
  }
  return(descriptives)
}

roundUp <- function(x,to=10)
{
  to*(x%/%to + as.logical(x%%to))
}
