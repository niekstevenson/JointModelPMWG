msit <- unique(experiments$MSIT$preppedData$subject) %>% levels()
rb <- unique(experiments$RB$preppedData$subject) %>% levels()
rlsat <- unique(experiments$rlsat$preppedData$subject) %>% levels()
rev <- unique(experiments$Rev$preppedData$subject) %>% levels()

allSubs <- c(msit, rb, rlsat, rev)
unqSubs <- unique(allSubs)
totaldf <- data.frame()

df <- data.frame(sub = unqSubs, experiment = "MSIT", done = unqSubs %in% msit)
totaldf <- rbind(totaldf, df)

df <- data.frame(sub = unqSubs, experiment = "RB", done = unqSubs %in% rb)
totaldf <- rbind(totaldf, df)

df <- data.frame(sub = unqSubs, experiment = "RLSAT", done = unqSubs %in% rlsat)
totaldf <- rbind(totaldf, df)

df <- data.frame(sub = unqSubs, experiment = "Rev", done = unqSubs %in% rev)
totaldf <- rbind(totaldf, df)

totaldf$sub <- as.character(totaldf$sub)


ggplot(totaldf, aes(x = reorder(sub, done), y = experiment, fill= done)) + 
  geom_tile(color = "white", size = 0.3) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))