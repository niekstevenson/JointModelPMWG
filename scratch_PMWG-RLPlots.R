# Plot posterior predictives
setwd("~/RL-PMWG/")
source("utils/utils-PMWG.R")
source("models/RW/RD.R")
load("samples/rlsat2V02B.RData")
dat <- experiments$rlsat$data
factors <- c("stimSet", "bin", "ease", "cue")
tmp <- pmwg_Simulate(sampled, dat, likelihood.RD, RL = T, factors = factors, n = 100)
ppInfo <- getDescriptives(tmp$pp, factors = factors)
dataInfo <- getDescriptives(tmp$data, factors = factors)


#Make some nice function of this
pdf('./figures/testrlsatFinal.pdf', width=3.5, height=7/4*3)
par(oma=c(3,4,1,0), mar=c(0, 0, 1, 0.5) + 0.1, mfcol=c(3,1), mgp=c(2.75,.75,0), las=1, bty='l')
data.cex=1.5
corrRTylim <- errRTylim <- c(.25, 1.1)

plotDataPPBins(data=dataInfo$AccBybin, pp=ppInfo$AccBybin,
               xaxt='n', draw.legend = T, data.cex = data.cex,
               dep.var='acc', ylab='', xlab = '', yaxt='n',
               legend.pos='topleft', ylim=c(0.2, 0.85), hline.by=0.1)
axis(1, at=seq(2, 10, 2), labels=rep(NA, 5), lwd=2)
mtext('Accuracy', side=2, cex=.66, line=3, las=0, font=1)
axis(2, at=seq(.5, .9, .1), lwd=1.5)
title("RL-EAM")

##
plotDataPPBins(data=dataInfo$RT.10.CorrectBybin, pp=ppInfo$RT.10.CorrectBybin, dep.var='RT.10.', 
               ylim=corrRTylim, xaxt='n', ylab='', yaxt='n', draw.legend = FALSE, data.cex = data.cex)
plotDataPPBins(data=dataInfo$RT.50.CorrectBybin, pp=ppInfo$RT.50.CorrectBybin, dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
plotDataPPBins(data=dataInfo$RT.90.CorrectBybin, pp=ppInfo$RT.90.CorrectBybin, dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
axis(1, at=seq(2, 10, 2), labels=NA, lwd=1.5)
mtext('Correct RTs (s)', side=2, cex=.66, line=3, las=0, font=1)
axis(2, seq(.2, 1.2, .2), lwd=1.5)

##
plotDataPPBins(data=dataInfo$RT.10.ErrorBybin, pp=ppInfo$RT.10.ErrorBybin, dep.var='RT.10.', ylim=errRTylim, xaxt='n', ylab='', yaxt='n', draw.legend = FALSE, data.cex = data.cex)
plotDataPPBins(data=dataInfo$RT.50.ErrorBybin, pp=ppInfo$RT.50.ErrorBybin, dep.var='RT.50.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
plotDataPPBins(data=dataInfo$RT.90.ErrorBybin, pp=ppInfo$RT.90.ErrorBybin, dep.var='RT.90.', plot.new = FALSE, draw.legend=FALSE, data.cex = data.cex)
mtext('Error RTs (s)', side=2, cex=.66, line=3, las=0, font=1)
axis(2, seq(.2, 1.2, .2), lwd=1.5)
axis(1, at=seq(2, 10, 2), lwd=1.5)
mtext('Trial bin', side=1, cex=.66, line=2)
dev.off()



# Plot --------------------------------------------------------------------
layoutM <- matrix(1:25, nrow=5, byrow=TRUE)
layoutM[c(1, 5),1:2] <- 1
layoutM[c(1, 5),4:5] <- 9
layoutM[2:4,1:2] <- 2:7 #matrix(c(2:13), nrow=3, byrow=TRUE)
layoutM[2:4,4:5] <- 10:15 #matrix(c(2:13), nrow=3, byrow=TRUE)
layoutM[,3] <- 8
layoutM

pdf('./figures/testrlsatFinal2.pdf', width=7, height=7/4*3)
layout(layoutM, heights = c(0.01, .8, 1, 1, 0.01), widths=c(1,1,.1,1,1))
par(oma=c(3,4,2,0), mar=c(0, 0, 1, 0.5) + 0.1, #mfcol=c(3,4), 
    mgp=c(2.75,.75,0), las=1, bty='l')
i <- 0
data.cex=1.5
corrRTylim <- errRTylim <- c(.35,1.1)

plot.new()
if(i == 0) mtext('RL-EAM', side=3, cex=.66*1.2, font=2, line=1)
if(i == 2) {plot.new(); mtext('RL-fARD', side=3, cex=.66*1.2, font=2, line=1)}
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

