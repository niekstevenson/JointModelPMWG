pmwg_jointParRec <- function(samples, experiments){
  parPreFixs <- gsub("[|].*", "", samples$par_names)
  i <- 0
  jointModelName <- paste(lapply(experiments, function(x){return(x$modelName)}), collapse = "_")
  medianPars <- pmwg_median(sampled)
  path <- file.path("./data", "ParRec", "joint", jointModelName)
  dir.create(path, recursive = T)
  epsilons <- c(0.3, 0.3, 0.5, 0.5)
  for (model in unique(parPreFixs)){
    i <- i + 1
    currentPars <- samples$par_names[which(parPreFixs == model)]
    currentPars <- medianPars[currentPars,]
    rownames(currentPars) <- experiments[[i]]$parNames
    
    origData <- experiments[[i]]$data
    data <- experiments[[i]]$preppedData
    RL <- attr(data, "RL")
    match <- attr(data, "match")
    data <- split(data, data$subject, drop = T)
    dat <- data.frame()
    for (name in names(data)){
      pp <- likelihood.RD(pars = currentPars[,name], data = data[[name]], sample = T)
      pp$subject <- name
      dat <- rbind(dat, pp)
    }
    
    factors <- experiments[[i]]$factors
    
    if(RL){
      factors <- c(factors, c("stimulus_set", "ease", "p_win_left", "p_win_right"))
      if(i == 4) factors <- c(factors, 'reversed')
    }

    dat[,factors] <- origData[,factors]
    if(RL){
      dat$p_win <- pmax(dat$p_win_left, dat$p_win_right)
      if(i == 4){
        dat$correctTotal[dat$reversed] <- ifelse(dat$R[dat$reversed] == 1, 2, 1)
        dat$correctTotal[!dat$reversed] <- dat$R
        dat$rewardTotal[dat$R == 2] <- rbinom(sum(dat$R == 2), 1, dat$p_win[dat$R == 2])
        dat$rewardTotal[dat$R == 1] <- rbinom(sum(dat$R == 1), 1, 1-dat$p_win[dat$R == 1])
        dat$reward <- ifelse(dat$reversed, ifelse(dat$rewardTotal == 1, 0, 1), dat$rewardTotal)
      } else{
        dat$reward[dat$R == 2] <- rbinom(sum(dat$R == 2), 1, dat$p_win[dat$R == 2])
        dat$reward[dat$R == 1] <- rbinom(sum(dat$R == 1), 1, 1-dat$p_win[dat$R == 1])
      }
    }
    dat[,match[[1]][1]] <- dat$R
    dat <- dat %>% select(-R)
    dat <- dat %>% rename(subjectNumber = subject, rt = RT)
    save(dat, file = paste0(path, "/", experiments[[i]]$modelName, ".RData"))
  }
}


pmwg_parRec <- function(experiments){
  for (exp in experiments){
    load(paste0("samples/", exp$modelName, ".RData"))
    medianPars <- pmwg_median(sampled)
    data <- sampled$data
    data <- split(data, data$subject, drop = T)
    df <- data.frame()
    for (name in names(data)){
      pp <- likelihood.RD(pars = medianPars[,name], data = data[[name]], sample = T)
      df <- rbind(df, pp)
    }
    
    attr(df, "RL") <- attr(sampled$data, "RL")
    attr(df, "constants") <- attr(sampled$data, "constants")
    attr(df, "match")<- attr(sampled$data, "match")
    attr(df, "n_v") <- attr(sampled$data, "n_v")
    
    priors <- list(
      theta_mu_mean = rep(0.1, length(exp$parNames)),
      theta_mu_var = diag(rep(9, length(exp$parNames)))
    )
    
    # Create the Particle Metropolis within Gibbs sampler object
    sampler <- pmwgs(
      data = df,
      pars = exp$parNames,
      ll_func = exp$llFunc,
      prior = priors
    )
    sampler = init(sampler, start_mu = exp$startPoints)
    burned <- run_stage(sampler, 
                        stage = "burn",
                        iter = 1000,
                        particles = 100,
                        n_cores = 20,
    )
    
    adapted <- run_stage(burned, 
                         stage = "adapt", 
                         iter = 10000, #Set up extremely high, terminates early anyway
                         particles = 100,
                         n_cores = 20,
    )
    
    sampled <- run_stage(adapted, 
                         stage = "sample",
                         iter = 1000, 
                         particles = 100,
                         n_cores = 20,
    )
    
    save(sampled, file = paste0("samples/ParRec/", exp$modelName, ".RData"))
  }
}

rmse <- function(x, y) {
  sqrt(mean((x-y)^2))
}

medParsRec <- pmwg_median(sampled)
medPars <- pmwg_median(sampled)

par(mfrow = c(3,2))
for (i in 1:nrow(medPars)){
  plot(medPars[i,], medParsRec[i,], main = rownames(medPars)[i], xlab = "Data generating", ylab = "Recovered")
  tmp <- legend('bottomright', c(" ", " "), bty='n', xjust=1, 
                text.width = strwidth("RMSE = 0.03"))
  text(tmp$rect$left + tmp$rect$w, tmp$text$y,
       c(paste0('r = ', round(cor(medPars[i,], medParsRec[i,]), 2)), 
         paste0('RMSE = ', round(rmse(medPars[i,], medParsRec[i,]), 2))), pos = 2)
  abline(a=0, b=1)
}