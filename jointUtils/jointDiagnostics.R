jointParInt <- function(joint, single, pars, PDF = F, path = NULL){
  joint$samples$alpha[grep("t0", rownames(joint$samples$alpha)),,] <- pnorm(joint$samples$alpha[grep("t0", rownames(joint$samples$alpha)),,])
  joint$samples$alpha[grep("aV", rownames(joint$samples$alpha)),,] <- pnorm(joint$samples$alpha[grep("aV", rownames(joint$samples$alpha)),,])
  
  single$samples$alpha[grep("t0", rownames(single$samples$alpha)),,] <- pnorm(single$samples$alpha[grep("t0", rownames(single$samples$alpha)),,])
  single$samples$alpha[grep("aV", rownames(single$samples$alpha)),,] <- pnorm(single$samples$alpha[grep("aV", rownames(single$samples$alpha)),,])
  
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


jointHyperSDPlot <- function(samples, experiments){
  samples$samples$theta_mu[grep("t0", rownames(samples$samples$theta_mu)),] <- pnorm(samples$samples$theta_mu[grep("t0", rownames(samples$samples$theta_mu)),])
  samples$samples$theta_mu[grep("aV", rownames(samples$samples$theta_mu)),] <- pnorm(samples$samples$theta_mu[grep("aV", rownames(samples$samples$theta_mu)),])
  jointParSDs <- apply(as_mcmc(samples, filter = "sample"), 2, sd) 
  singleParSDs <- numeric()
  for(i in 1:length(experiments)){
    load(paste0("samples/", experiments[[i]]$modelName, ".RData"))
    sampled$samples$theta_mu[grep("t0", rownames(sampled$samples$theta_mu)),] <- pnorm(sampled$samples$theta_mu[grep("t0", rownames(sampled$samples$theta_mu)),])
    sampled$samples$theta_mu[grep("aV", rownames(sampled$samples$theta_mu)),] <- pnorm(sampled$samples$theta_mu[grep("aV", rownames(sampled$samples$theta_mu)),])
    
    parSDs <- apply(as_mcmc(sampled, filter = "sample"), 2, sd) 
    singleParSDs <- c(singleParSDs, parSDs)
  }
  #parNames <- gsub(".*[|]", "", names(jointParSDs))
  jointParSDs <- c(jointParSDs, singleParSDs)
  df <- data.frame(SD = jointParSDs, group = c(rep("joint", length(jointParSDs)/2), rep("single", length(jointParSDs)/2)), 
                   name = rep(names(jointParSDs)[1:(length(jointParSDs)/2)], 2))
  
  ggplot(data=df, aes(x=name, y=SD, fill=group)) +
    geom_bar(stat="identity", color="black", position=position_dodge())+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
}

jointSubSDPlot <- function(samples, experiments){
  samples$samples$alpha[grep("t0", rownames(samples$samples$alpha)),,] <- pnorm(samples$samples$alpha[grep("t0", rownames(samples$samples$alpha)),,])
  samples$samples$alpha[grep("aV", rownames(samples$samples$alpha)),,] <- pnorm(samples$samples$alpha[grep("aV", rownames(samples$samples$alpha)),,])
  jointParSDs <- rowMeans(parSD(samples))  
  singleParSDs <- numeric()
  for(i in 1:length(experiments)){
    load(paste0("samples/", experiments[[i]]$modelName, ".RData"))
    sampled$samples$alpha[grep("t0", rownames(sampled$samples$alpha)),,] <- pnorm(sampled$samples$alpha[grep("t0", rownames(sampled$samples$alpha)),,])
    sampled$samples$alpha[grep("aV", rownames(sampled$samples$alpha)),,] <- pnorm(sampled$samples$alpha[grep("aV", rownames(sampled$samples$alpha)),,])
    parSDs <- rowMeans(parSD(sampled)) 
    singleParSDs <- c(singleParSDs, parSDs)
  }
  jointParSDs <- c(jointParSDs, singleParSDs)
  df <- data.frame(SD = jointParSDs, group = c(rep("joint", length(jointParSDs)/2), rep("single", length(jointParSDs)/2)), 
                   name = rep(names(jointParSDs)[1:(length(jointParSDs)/2)], 2))
  df$name <- gsub("MSITprocsB2", "MSIT", df$name)
  df$name <- gsub("RB4V02Diff2B", "RB", df$name)
  df$name <- gsub("rlsat2V02B", "RLSAT", df$name)
  df$name <- gsub("RevStd", "Rev", df$name)
  
  
  ggplot(data=df, aes(x=name, y=SD, fill=group)) +
    geom_bar(stat="identity", color="black", position=position_dodge())+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
}

jointChainPlots <- function(samples, PDF = F, path = NULL){
  if(PDF) pdf(paste0(path, "/chainPlots.pdf"))
  for (exp in samples$experiments){
    name <- exp$modelName
    idx <- grep(name, samples$par_names, fixed = F)
    for (par in samples$par_names[idx]){
      matplot(t(samples$samples$alpha[par,1:6,samples$samples$stage == "sample"]),type="l", main = par, xlab = "samples", ylab = "ParameterValue")
    } 
    matplot(t(samples$samples$theta_mu[idx, samples$samples$stage == "sample"]), type="l", main = paste0(name, ": Paramater chains"), ylab = "Parameter Value", xlab = "samples")
  }
  matplot(t(samples$samples$subj_ll[, samples$samples$stage == "sample"]), type="l", main = "LogLikelihood chains per subject", ylab = "LogLikelihood", xlab = "samples")
  dev.off()
  if(PDF) while (!is.null(dev.list()))  dev.off()
  
}

