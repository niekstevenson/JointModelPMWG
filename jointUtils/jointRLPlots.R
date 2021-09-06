jointRLSATPlot <- function(pp, data, joint, reversal = F, PDF = F, path = NULL){
  ppSummary <- postSummary(pp, c("cue", "bin"), "single")

  dataSummary <- postSummary(data, c("cue", "bin"), "data")
  df <- rbind(ppSummary, dataSummary)
  
  if(!is.null(joint)){
    jointSummary <- postSummary(joint, c("cue", "bin"), "joint")
    df <- rbind(df, jointSummary)
  }
  
  df <- df %>% mutate(bin = as.factor(bin)) %>% group_by(cue, bin, source, quant, acc) %>% 
    summarise(lower = quantile(value, 0.05), higher = quantile(value, 0.95), value = mean(value))
  
  rtlim <- c(min(df$value[df$acc != 9 ]), max((df$value[df$acc != 9])))
  acclim <- c(min(df$value[df$acc == 9]), max((df$value[df$acc == 9])))
  
  entries <- list()
  for(j in c(9, 1, 0)){
    for(i in c('SPD', 'ACC')) {
      if(j == 9){
        ylab <- "Accuracy"
        ylim <- acclim
      } else if(j == 1){
        ylab <- "Correct RTs"
        ylim <- rtlim
      } else{
        ylab <- "Error RTs"
        ylim <- rtlim
      }
      p <- ggplot(df %>% filter(cue == i, acc == j), aes(x = bin, y = value, group = interaction(source, quant), colour = source)) +
        geom_ribbon(aes(ymin = lower, ymax = higher, fill = source), alpha = .3, size = 0.35) +
        geom_line(data = df %>% filter(cue == i, acc == j, source == "data"), size = .9) + 
        geom_point(data = df %>% filter(cue == i, acc == j, source == "data"), size = 1.5) + 
        labs(x = "bin", y= ylab) +
        coord_cartesian(ylim = ylim) + 
        theme_bw()
      if(j == 9){
        p <- p + ggtitle(i)
      }
      entries <- c(entries, list(p))
    }
  }
  if(PDF) png(paste0(path, '/BinPlotPerCue.png'), width=6, height=8, units = "in", res = 300)
  plot <- do.call(ggarrange, c(entries, nrow = 3, ncol = 2, common.legend = T))
  print(plot)
  dev.off()
}

postSummaryRev <- function(data, combns, source){
  df <- dplyr::bind_rows(data)
  dfRTs <-  df %>% group_by(subject, across(combns), reps) %>% 
    summarise(value = quantile(RT, c(0.1, 0.5, 0.9)), quant = c(0.1, 0.5, 0.9)) %>% 
    ungroup() %>% group_by(across(combns), quant, reps) %>% summarise_all(funs(mean)) %>% mutate(acc = 0.5)
  
  dfAcc <- df %>% group_by(subject, across(combns), reps) %>% 
    summarise(value = mean(acc), quant = 9) %>% 
    ungroup() %>% group_by(across(combns), quant, reps) %>% summarise_all(funs(mean)) %>% mutate(acc = 9)
  
  df <- rbind(dfAcc, dfRTs) %>% mutate(source = source) %>% select(-subject)
}


jointRevPlot <- function(pp, data, joint, reversal = F, PDF = F, path = NULL){
  ppSummary <- postSummaryRev(pp, c("bin"), "single")
  dataSummary <- postSummaryRev(data, c("bin"), "data")
  df <- rbind(ppSummary, dataSummary) 
  
  if(!is.null(joint)){
    jointSummary <- postSummaryRev(joint, c("bin"), "joint")
    df <- rbind(df, jointSummary) 
  }
  
  df <- df %>% group_by(bin, source, quant, acc) %>% mutate(bin = as.factor(bin)) %>% 
    summarise(lower = quantile(value, 0.05), higher = quantile(value, 0.95), value = mean(value))
  
  rtlim <- c(min(df$value[df$acc != 9 ]), max((df$value[df$acc != 9])))
  acclim <- c(min(df$value[df$acc == 9]), max((df$value[df$acc == 9])))
  
  entries <- list()
  for(j in c(9, 0.5)){
    if(j == 9){
      ylab <- "% choice A"
      ylim <- acclim
    } else if(j == 0.5){
      ylab <- "RTs"
      ylim <- rtlim
    }
    p <- ggplot(df %>% filter(acc == j), aes(x = bin, y = value, group = interaction(source, quant), colour = source)) +
      geom_ribbon(aes(ymin = lower, ymax = higher, fill = source), alpha = .3, size = 0.35) +
      geom_line(data = df %>% filter(acc == j, source == "data"), size = .9) + 
      geom_point(data = df %>% filter(acc == j, source == "data"),  size = 1.5) + 
      labs(x = "bin", y= ylab) +
      coord_cartesian(ylim = ylim) + 
      theme_bw()
    entries <- c(entries, list(p))
  }
  if(PDF) png(paste0(path, '/BinPlotPerCue.png'), width=3, height=5.3, units = "in", res = 300)
  plot <- do.call(ggarrange, c(entries, nrow = 2, ncol = 1, common.legend = T))
  print(plot)
  dev.off()
}


