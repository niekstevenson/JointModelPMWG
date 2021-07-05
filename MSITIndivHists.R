df <- experiments$MSIT$preppedData

df$acc <- df$resp == df$uniq
pdf('./figures/MSITrtsPerSub.pdf', width=7, height=7/4*3)
par(mfrow = c(1, 2))
for(sub in unique(df$subject)){

  dfSub <- df[df$subject == sub,]
  accR1 <- mean(dfSub[dfSub$resp == 1,]$resp == dfSub[dfSub$resp == 1,]$uniq)
  accR2 <- mean(dfSub[dfSub$resp == 2,]$resp == dfSub[dfSub$resp == 2,]$uniq)
  accR3 <- mean(dfSub[dfSub$resp == 3,]$resp == dfSub[dfSub$resp == 3,]$uniq)
  
  R1 <- mean(dfSub$R == 1)
  R2 <- mean(dfSub$R == 2)
  R3 <- mean(dfSub$R == 3)
  p <- ggplot(dfSub, aes(x=RT, fill=R)) +
    geom_histogram(aes(y=..density..), position = "dodge", alpha = 1, bins = 30) +
    geom_density(alpha=0.4)+
    ggtitle(paste0("Subject ", sub)) +
    theme_bw() +
    annotate("text", x = 0.95*max(dfSub$RT), y = c(3.3, 3.65, 4), label = c(paste0('Acc R1 = ', round(accR1, 2)),
                                                                                paste0('Acc R2 = ', round(accR2, 2)),
                                                                                 paste0('Acc R3 = ', round(accR3, 2)))) + 
    annotate("text", x = 0.75*max(dfSub$RT), y = c(3.3, 3.65, 4), label = c(paste0('% R1 = ', round(R1, 2)),
                                                                            paste0('% R2 = ', round(R2, 2)),
                                                                            paste0('% R3 = ', round(R3, 2))))
  
  
  print(p)
  # tmp <- legend('topright', c(" ", " ", " "), bty='n', xjust=1, 
  #               text.width = strwidth("Acc = 0.03"))
  # text(tmp$rect$left + tmp$rect$w, tmp$text$y,
  #      c(paste0('Acc R1 = ', round(accR1, 2)),
  #       paste0('Acc R2 = ', round(accR2, 2)),
  #       paste0('Acc R2 = ', round(accR3, 2))), pos = 2)
}
dev.off()
par(mfrow = c(1, 1))
newdf <- data.frame(sub = unique(df$subject), meanR2 = NA, accR2 = NA)
i <- 0
for (sub in unique(df$subject)){
  dfSub <- df[df$subject == sub,]
  i <- i + 1
  newdf$meanR2[i] <- mean(dfSub$R == 2)
  newdf$accR2[i] <- mean(dfSub$resp == dfSub$uniq)
}
plot(newdf$meanR2, newdf$accR2, xlab = "proportion R2", ylab = "percentage correct for R2")
cor(newdf$meanR2, newdf$accR2)
