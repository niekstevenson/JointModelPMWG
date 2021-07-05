pmwg_jointRLSATPlot <- function(ppInfo, dataInfo, jointInfo, reversal = F, PDF = F, path = NULL){
  
  # Plot --------------------------------------------------------------------
  layoutM <- matrix(1:25, nrow=5, byrow=TRUE)
  layoutM[c(1, 5),1:2] <- 1
  layoutM[c(1, 5),4:5] <- 9
  layoutM[2:4,1:2] <- 2:7 #matrix(c(2:13), nrow=3, byrow=TRUE)
  layoutM[2:4,4:5] <- 10:15 #matrix(c(2:13), nrow=3, byrow=TRUE)
  layoutM[,3] <- 8
  layoutM
  if(PDF) png(paste0(path, '/BinPlotPerCue.png'), width=7, height=7/4*3, units = "in", res = 300)
  layout(layoutM, heights = c(0.01, .8, 1, 1, 0.01), widths=c(1,1,.1,1,1))
  par(oma=c(3,4,2,0), mar=c(0, 0, 1, 0.5) + 0.1, mgp=c(2.75,.75,0), las=1, bty='l')
  i <- 0
  data.cex=1.5
  corrRTylim <- errRTylim <- c(.35,1.1)
  
  plot.new()
  if(i == 0) mtext('RL-ARD', side=3, cex=.66*1.2, font=2, line=1)
  for(cue in c('SPD', 'ACC')) {
    i <- i+1
    idxD <- dataInfo$AccBybinBycue$cue==cue
    idxM <- ppInfo$AccBybinBycue$cue==cue
    idxJ <- jointInfo$AccBybinBycue$cue==cue
    
    plotDataPPBins(data=dataInfo$AccBybinBycue[idxD,], pp=ppInfo$AccBybinBycue[idxM,], joint_pp = jointInfo$AccBybinBycue[idxM,],
                   xaxt='n', draw.legend = i==1, data.cex = data.cex,
                   dep.var='acc', ylab='', xlab = '', yaxt='n',
                   legend.pos='topleft', ylim=c(0.5, 0.9), hline.by=0.1)
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
    plotDataPPBins(data=dataInfo$RT.10.CorrectBybinBycue[idxD,], pp=ppInfo$RT.10.CorrectBybinBycue[idxM,], joint_pp=jointInfo$RT.10.CorrectBybinBycue[idxJ,],
                   dep.var='RT.10.', ylim=corrRTylim, xaxt='n', ylab='', yaxt='n', draw.legend = FALSE, data.cex = data.cex)
    plotDataPPBins(data=dataInfo$RT.50.CorrectBybinBycue[idxD,], pp=ppInfo$RT.50.CorrectBybinBycue[idxM,], joint_pp=jointInfo$RT.50.CorrectBybinBycue[idxJ,],
                   dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
    plotDataPPBins(data=dataInfo$RT.90.CorrectBybinBycue[idxD,], pp=ppInfo$RT.90.CorrectBybinBycue[idxM,], joint_pp=jointInfo$RT.90.CorrectBybinBycue[idxJ,],
                   dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
    axis(1, at=seq(2, 10, 2), labels=NA, lwd=1.5)
    if(i == 1) {
      mtext('Correct RTs (s)', side=2, cex=.66, line=3, las=0, font=1)
      axis(2, seq(.2, 1.6, .2), lwd=1.5)
    } else {
      axis(2, seq(.2, 1.6, .2), labels=NA, lwd=1.5)
    }
    
    ##
    plotDataPPBins(data=dataInfo$RT.10.ErrorBybinBycue[idxD,], pp=ppInfo$RT.10.ErrorBybinBycue[idxM,], joint_pp=jointInfo$RT.10.ErrorBybinBycue[idxJ,],
                   dep.var='RT.10.', ylim=errRTylim, xaxt='n', ylab='', yaxt='n', draw.legend = FALSE, data.cex = data.cex)
    plotDataPPBins(data=dataInfo$RT.50.ErrorBybinBycue[idxD,], pp=ppInfo$RT.50.ErrorBybinBycue[idxM,], joint_pp=jointInfo$RT.50.ErrorBybinBycue[idxJ,],
                   dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
    plotDataPPBins(data=dataInfo$RT.90.ErrorBybinBycue[idxD,], pp=ppInfo$RT.90.ErrorBybinBycue[idxM,], joint_pp=jointInfo$RT.90.ErrorBybinBycue[idxJ,],
                   dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
    if(i == 1) {
      mtext('Error RTs (s)', side=2, cex=.66, line=3, las=0, font=1)
      axis(2, seq(.4, 1.2, .2), lwd=1.5)
    } else {
      axis(2, seq(.4, 1.2, .2), labels=NA, lwd=1.5)
    }
    axis(1, at=seq(-10, 10, 2), lwd=1.5)
    mtext('Trial bin', side=1, cex=.66, line=2)
    
  }
  
  dev.off()
}

pmwg_jointRevPlot <- function(ppInfo, dataInfo, jointInfo, reversal = F, PDF = F, path = NULL){
  if(PDF) png(paste0(path, '/BinPlot.png'), width=3.5, height=7/4*3, units = "in", res = 300)
  par(oma=c(3,4,1,0), mar=c(0, 0, 1, 0.5) + 0.1, mfcol=c(2,1), mgp=c(2.75,.75,0), las=1, bty='l')
  data.cex=1.5
  corrRTylim <- errRTylim <- c(.25, 1.1)
  
  plotDataPPBins(data=dataInfo$AccBybin, pp=ppInfo$AccBybin, joint_pp = jointInfo$AccBybin,
                 xaxt='n', draw.legend = T, data.cex = data.cex,
                 dep.var='acc', ylab='', xlab = '', yaxt='n',
                 legend.pos='bottomleft', ylim=c(0.2, 0.85), hline.by=0.1)
  axis(1, at=seq(-10, 10, 2), labels=rep(NA, 5), lwd=2)
  mtext('proportion choice A', side=2, cex=.66, line=3, las=0, font=1)
  axis(2, at=seq(.1, .9, .2), lwd=1.5)
  title("RL-ARD")
  
  ##
  plotDataPPBins(data=dataInfo$RT.10.Bybin, pp=ppInfo$RT.10.Bybin, joint_pp=jointInfo$RT.10.Bybin, dep.var='RT.10.', 
                 ylim=corrRTylim, xaxt='n', ylab='', yaxt='n', draw.legend = FALSE, data.cex = data.cex)
  plotDataPPBins(data=dataInfo$RT.50.Bybin, pp=ppInfo$RT.50.Bybin, joint_pp=jointInfo$RT.50.Bybin, dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  plotDataPPBins(data=dataInfo$RT.90.Bybin, pp=ppInfo$RT.90.Bybin, joint_pp=jointInfo$RT.90.Bybin, dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
  axis(1, at=seq(-10, 10, 2), labels=NA, lwd=1.5)
  mtext('RTs (s)', side=2, cex=.66, line=3, las=0, font=1)
  axis(2, seq(.2, 1.2, .2), lwd=1.5)
  axis(1, at=seq(2, 10, 2), lwd=1.5)
  mtext('Trial bin', side=1, cex=.66, line=2)
  dev.off()
}


