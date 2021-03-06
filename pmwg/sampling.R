### This is the pmwg sampler for standard puproses. It includes a multivariate normal with full covariance matrix on the group level. 

library(mvtnorm) ## For the multivariate normal.
library(MASS) ## For matrix inverse.
library(MCMCpack) #For inverse wishart
library(lme4)
library(parallel)
library(abind)
library(checkmate)
library(Rcpp)

source("pmwg/utils.R")
source("pmwg/messaging.R")

pmwgs <- function(data, pars, ll_func, prior = NULL) {
  ###Gets/sets priors, creates pmwgs object and stores essentials
  # Descriptives
  n_pars <- length(pars)
  subjects <- unique(data$subject)
  n_subjects <- length(subjects)

  # Hyperparameters
  v_half <- 2 # hyperparameter on Σ prior (Half-t degrees of freedom)
  A_half <- 1 # hyperparameter on Σ prior (Half-t scale) #nolint
  
  # Storage for the samples.
  samples <- sample_store(pars, subjects)
  # Checking and default priors
  if (is.null(prior)) {
    prior <- list(theta_mu_mean = rep(0, n_pars), theta_mu_var = diag(rep(1, n_pars)))
  }
  # Things I save rather than re-compute inside the loops.
  if(is.matrix(prior$theta_mu_var)){
    prior$theta_mu_invar <- MASS::ginv(prior$theta_mu_var) #Inverse of the matrix
  } else{
    prior$theta_mu_invar <- 1/prior$theta_mu_var #Inverse of the matrix
  }
  sampler <- list(
    data = data,
    par_names = pars,
    n_pars = n_pars,
    n_subjects = n_subjects,
    subjects = subjects,
    prior = prior,
    ll_func = ll_func,
    samples = samples,
    init = FALSE
  )
  #Hyper parameters
  attr(sampler, "v_half") <- v_half
  attr(sampler, "A_half") <- A_half
  class(sampler) <- "pmwgs"
  sampler
}

init <- function(pmwgs, start_mu = NULL, start_var = NULL,
         display_progress = TRUE, particles = 1000, n_cores = 1, epsilon = NULL, useC = T) {
  # Gets starting points for the mcmc process
  # If no starting point for group mean just use zeros
  if (is.null(start_mu)) start_mu <- stats::rnorm(pmwgs$n_pars, sd = 1)
  # If no starting point for group var just sample some
  if (is.null(start_var)) {
    #If prior on covariances is a vector, assume diagonal only, otherwise assume full cvs structure
    if(is.matrix(pmwgs$prior$theta_mu_var)){
      start_var <- MCMCpack::riwish(pmwgs$n_pars * 3,diag(pmwgs$n_pars))
    } else{
      start_var <- diag(1/rgamma(pmwgs$n_pars, 10, 5)) #But stupid maybe as startpoint
    }
  }
  
  if(useC){
    sourceCpp("pmwg/utilityFunctions.cpp")
    dmv <<- dmvnrm_arma_fast
  } else{
    dmv <<- mvtnorm::dmvnorm
  }
  
  # Sample the mixture variables' initial values.
  a_half <- 1 / stats::rgamma(n = pmwgs$n_pars, shape = 2, rate = 1)
  # Create and fill initial random effects for each subject
  likelihoods <- array(NA_real_, dim = c(pmwgs$n_subjects))
  if(n_cores > 1){
    proposals <- mclapply(X=1:pmwgs$n_subjects,FUN=start_proposals,start_mu = start_mu, 
                        start_var = start_var, n_particles = particles, pmwgs = pmwgs, mc.cores = n_cores)
  } else{
    proposals <- lapply(X=1:pmwgs$n_subjects,FUN=start_proposals,start_mu = start_mu, 
               start_var = start_var, n_particles = particles, pmwgs = pmwgs)
  }
  proposals <- simplify2array(proposals)
  pmwgs$init <- TRUE
  pmwgs$samples$theta_mu[, 1] <- start_mu
  pmwgs$samples$theta_var[, , 1] <- start_var
  pmwgs$samples$alpha[, , 1] <- do.call(cbind, proposals[1,])
  pmwgs$samples$last_theta_var_inv <- MASS::ginv(start_var)
  pmwgs$samples$subj_ll[, 1] <- unlist(proposals[2,])
  pmwgs$samples$a_half[, 1] <- a_half
  pmwgs$samples$idx <- 1
  pmwgs$samples$epsilon[,1] <- rep(set_epsilon(pmwgs$n_pars, epsilon), pmwgs$n_subjects)
  pmwgs$samples$origin[,1] <- rep(2, pmwgs$n_subjects)
  return(pmwgs)
}

start_proposals <- function(s, start_mu, start_var, n_particles, pmwgs){
  #Draw the first start point
  proposals <- particle_draws(n_particles, start_mu, start_var)
  colnames(proposals) <- rownames(pmwgs$samples$theta_mu) # preserve par names
  lw <- apply(proposals,1,pmwgs$ll_func,data = pmwgs$data[pmwgs$data$subject == pmwgs$subjects[s], ])
  weight <- exp(lw - max(lw))
  idx <- sample(x = n_particles, size = 1, prob = weight)
  return(list(proposal = proposals[idx,], ll = lw[idx]))
}

gibbs_step_diag<- function(sampler){
  # Gibbs step for diagonal only
  # Get single iter versions, tmu = theta_mu, tvar = theta_var
  last <- last_sample(sampler$samples)
  hyper <- attributes(sampler)
  prior <- sampler$prior
  
  prior$theta_mu_invar <- diag(prior$theta_mu_invar)
  last$tvinv <- diag(last$tvinv)
  
  #Mu
  var_mu = 1.0 / (sampler$n_subjects * last$tvinv + prior$theta_mu_invar)
  mean_mu = var_mu * ((apply(last$alpha, 1, sum) * last$tvinv + prior$theta_mu_mean * prior$theta_mu_invar))
  tmu <- rnorm(sampler$n_pars, mean_mu, sd = sqrt(var_mu))
  names(tmu) <- sampler$par_names
  
  #InvGamma alternative (probably inferior) prior
  # shape = hyper$shape + sampler$n_subjects / 2
  # rate = hyper$rate + rowSums( (last$alpha-tmu)^2 ) / 2
  # tvinv = rgamma(n=sampler$n_pars, shape=shape, rate=rate)
  # tvar = 1/tvinv
  tvinv = rgamma(n=sampler$n_pars, shape=hyper$v_half/2 + sampler$n_subjects/2, rate=hyper$v_half/last$a_half + 
                   rowSums( (last$alpha-tmu)^2 ) / 2)
  tvar = 1/tvinv
  #Contrary to standard pmwg I use shape, rate for IG()
  a_half <- 1 / stats::rgamma(n = sampler$n_pars, shape = (hyper$v_half + sampler$n_pars) / 2,
                              rate = hyper$v_half * tvinv + 1/hyper$A_half)
  return(list(tmu = tmu, tvar = diag(tvar), tvinv = diag(tvinv), a_half = a_half, alpha = last$alpha))
}

gibbs_step_full <- function(sampler){
  # Gibbs step for group means, with full covariance matrix estimation
  # tmu = theta_mu, tvar = theta_var
  last <- last_sample(sampler$samples)
  hyper <- attributes(sampler)
  prior <- sampler$prior
  
  # Here mu is group mean, so we are getting mean and variance
  var_mu <- MASS::ginv(sampler$n_subjects * last$tvinv + prior$theta_mu_invar)
  mean_mu <- as.vector(var_mu %*% (last$tvinv %*% apply(last$alpha, 1, sum) +
                                     prior$theta_mu_invar %*% prior$theta_mu_mean))
  chol_var_mu <- t(chol(var_mu)) # t() because I want lower triangle.
  # New sample for mu.
  tmu <- mvtnorm::rmvnorm(1, mean_mu, chol_var_mu %*% t(chol_var_mu))[1, ]
  names(tmu) <- sampler$par_names
  
  # New values for group var
  theta_temp <- last$alpha - tmu
  cov_temp <- (theta_temp) %*% (t(theta_temp))
  B_half <- 2 * hyper$v_half * diag(1 / last$a_half) + cov_temp # nolint
  tvar <- MCMCpack::riwish(hyper$v_half + sampler$n_pars - 1 + sampler$n_subjects, B_half) # New sample for group variance
  tvinv <- MASS::ginv(tvar)
  
  # Sample new mixing weights.
  a_half <- 1 / stats::rgamma(n = sampler$n_pars,shape = (hyper$v_half + sampler$n_pars) / 2,
                              rate = hyper$v_half * diag(tvinv) + hyper$A_half)
  return(list(tmu = tmu,tvar = tvar,tvinv = tvinv,a_half = a_half,alpha = last$alpha))
}

new_particle <- function (s, data, num_particles, parameters, eff_mu = NULL, 
                          eff_var = NULL, mix_proportion = c(0.5, 0.5, 0), 
                          likelihood_func = NULL, epsilon = NULL, subjects) 
{
  eff_mu <- eff_mu[, s]
  eff_var <- eff_var[, , s]
  mu <- parameters$tmu
  var <- parameters$tvar
  subj_mu <- parameters$alpha[, s]
  particle_numbers <- numbers_from_proportion(mix_proportion, num_particles)
  cumuNumbers <- cumsum(particle_numbers)
  pop_particles <- particle_draws(particle_numbers[1], mu, 
                                  var)
  ind_particles <- particle_draws(particle_numbers[2], subj_mu, 
                                  var * epsilon[s]^2)
  if(mix_proportion[3] == 0){
    eff_particles <- NULL
  } else{
    eff_particles <- particle_draws(particle_numbers[3], eff_mu, eff_var)
  }
  proposals <- rbind(pop_particles, ind_particles, eff_particles)
  colnames(proposals) <- names(mu)
  proposals[1, ] <- subj_mu
  lw <- apply(proposals, 1, likelihood_func, data = data[data$subject==subjects[s],])
  lp <- dmv(x = proposals, mean = mu, sigma = var, 
                         log = TRUE)
  prop_density <- dmv(x = proposals, mean = subj_mu, 
                                   sigma = var * (epsilon[s]^2))
  if (mix_proportion[3] == 0) {
    eff_density <- 0
  }
  else {
    eff_density <- dmv(x = proposals, mean = eff_mu, sigma = eff_var)
  }
  lm <- log(mix_proportion[1] * exp(lp) + (mix_proportion[2] * 
                                             prop_density) + (mix_proportion[3] * eff_density))
  l <- lw + lp - lm
  weights <- exp(l - max(l))
  idx <- sample(x = num_particles, size = 1, prob = weights)
  origin <- min(which(idx <= cumuNumbers))
  return(list(ll = lw[idx], origin = origin, proposal = proposals[idx, ]))
}

get_conditionals_diag <- function(s, samples, n_pars, iteration){
  pts2_unwound <- apply(samples$theta_var,3,diag)
  all_samples <- rbind(samples$alpha[, s,],samples$theta_mu,pts2_unwound)
  mu_tilde <- apply(all_samples, 1, mean)
  var_tilde <- stats::var(t(all_samples))
  condmvn <- condMVNorm::condMVN(mean = mu_tilde, sigma = var_tilde,
                                 dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                                 X.given = c(samples$theta_mu[,iteration], diag(samples$theta_var[,,iteration])))
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
}

get_conditionals_full <- function(s, samples, n_pars, iteration){
  pts2_unwound <- apply(samples$theta_var,3,unwind)
  all_samples <- rbind(samples$alpha[, s,],samples$theta_mu,pts2_unwound)
  mu_tilde <- apply(all_samples, 1, mean)
  var_tilde <- stats::var(t(all_samples))
  condmvn <- condMVNorm::condMVN(mean = mu_tilde, sigma = var_tilde,
                                 dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                                 X.given = c(samples$theta_mu[,iteration], unwind(samples$theta_var[,,iteration])))
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
}

run_stage <- function(pmwgs,
                      stage,
                      iter = 1000,
                      particles = 1000,
                      display_progress = TRUE,
                      n_cores = 1,
                      n_unique = ifelse(stage == "adapt", 20, NA),
                      epsilon = NULL,
                      pstar = NULL,
                      mix = NULL,
                      diag_only = F,
                      pdist_update_n = ifelse(stage == "sample", 500, NA),
                      useC = T
                      ) {
  # Set defaults for NULL values
  mix <- set_mix(stage, mix)
  # Set necessary local variables
  .n_unique <- n_unique
  # Set stable (fixed) new_sample argument for this run
  n_pars <- length(pmwgs$par_names)
  
  if(useC){
    sourceCpp("pmwg/utilityFunctions.cpp")
    dmv <<- dmvnrm_arma_fast
  } else{
    dmv <<- mvtnorm::dmvnorm
  }
  
  if(is.matrix(pmwgs$prior$theta_mu_var)){
    gibbs_step <- gibbs_step_full
    get_conditionals <- get_conditionals_full
  } else{
    gibbs_step <- gibbs_step_diag
    get_conditionals <- get_conditionals_diag
  }
  
  # Display stage to screen
  msgs <- list(
    burn = "Phase 1: Burn in\n",
    adapt = "Phase 2: Adaptation\n",
    sample = "Phase 3: Sampling\n"
  )
  cat(msgs[[stage]])
  
  alphaStar=-qnorm(pstar/2) #Idk about this one
  n0=round(5/(pstar*(1-pstar))) #Also not questioning this math for now
  if(is.null(epsilon)){
    epsilon <- pmwgs$samples$epsilon[,ncol(pmwgs$samples$epsilon)]
  }
  if(length(epsilon) == 1){
    epsilon <- rep(epsilon, pmwgs$n_subjects)
  }
  
  # Build new sample storage
  pmwgs <- extend_sampler(pmwgs, iter, stage)
  # create progress bar
  if (display_progress) {
    pb <- accept_progress_bar(min = 0, max = iter)
  }
  start_iter <- pmwgs$samples$idx
  
  eff_mu <- NULL
  eff_var <- NULL
  data <- pmwgs$data
  subjects <- pmwgs$subjects
  # Main iteration loop
  for (i in 1:iter) {
    accRate <- mean(accept_rate(pmwgs))
    if (display_progress) {
      update_progress_bar(pb, i, extra = accRate)
    }
    # Create/update efficient proposal distribution if we are in sampling phase.
    if(stage == "sample" & (i %% pdist_update_n == 0 || i == 1)){
      test_samples <- extract_samples(pmwgs, stage = c("adapt", "sample"))
      iteration <- dim(test_samples$theta_mu)[2]
      if(n_cores > 1){
        conditionals=parallel::mclapply(X = 1:pmwgs$n_subjects,FUN = get_conditionals,samples = test_samples, 
                              n_pars, iteration, mc.cores = n_cores)
      } else{
        conditionals=lapply(X = 1:pmwgs$n_subjects,FUN = get_conditionals,samples = test_samples, 
                            n_pars, iteration)
      }
      conditionals <- simplify2array(conditionals)
      eff_mu <- do.call(cbind, conditionals[1,])
      eff_var <- abind(conditionals[2,], along = 3)
    }
    
    pars <- gibbs_step(pmwgs)
    if(n_cores > 1){
      proposals=parallel::mclapply(X=1:pmwgs$n_subjects,FUN = new_particle, data, particles, pars, eff_mu, 
                   eff_var, mix, pmwgs$ll_func, epsilon, subjects, mc.cores =n_cores)
    } else{
      proposals=lapply(X=1:pmwgs$n_subjects, FUN = new_particle, data, particles, pars, eff_mu, 
                 eff_var, mix, pmwgs$ll_func, epsilon, subjects)
    }
    proposals <- array(unlist(proposals), dim = c(pmwgs$n_pars + 2, pmwgs$n_subjects))
    ll <- proposals[1,]
    origin <- proposals[2,]
    alpha <- proposals[3:(pmwgs$n_pars+2),]

    j <- start_iter + i
    
    pmwgs$samples$theta_mu[, j] <- pars$tmu
    pmwgs$samples$theta_var[, , j] <- pars$tvar
    pmwgs$samples$last_theta_var_inv <- pars$tvinv
    pmwgs$samples$alpha[, , j] <- alpha
    pmwgs$samples$idx <- j
    pmwgs$samples$subj_ll[, j] <- ll
    pmwgs$samples$origin[,j] <- origin
    pmwgs$samples$a_half[, j] <- pars$a_half
    
    if(!is.null(pstar)){
      if(j > n0){
        acc <- alpha[1,] != pmwgs$samples$alpha[1,,(j-1)]
        epsilon<-update.epsilon(epsilon^2, acc, pstar, j, pmwgs$n_pars, alphaStar)
      }
    }
    
    pmwgs$samples$epsilon[,j] <- epsilon
    
    if (stage == "adapt") {
      res <- test_sampler_adapted(pmwgs, n_unique, i, n_cores, get_conditionals)
      if (res == "success") {
        break
      } else if (res == "increase") {
        n_unique <- n_unique + .n_unique
      }
    }
  }
  if (display_progress) close(pb)
  if (stage == "adapt") {
    if (i == iter) {
      message(paste(
        "Particle Metropolis within Gibbs Sampler did not",
        "finish adaptation phase early (all", i, "iterations were",
        "run).\nYou should examine your samples and perhaps start",
        "a longer adaptation run."
      ))
    } else {
      pmwgs <- trim_na(pmwgs)
    }
  }
  return(pmwgs)
}

test_sampler_adapted <- function(pmwgs, n_unique, i, n_cores, conditionals_func) {
  n_pars <- length(pmwgs$par_names)
  if (i < n_unique) {
    return("continue")
  }
  test_samples <- extract_samples(pmwgs, stage = "adapt")
  # Only need to check uniqueness for one parameter
  first_par <- test_samples$alpha[1, , ]
  # Split the matrix into a list of vectors by subject
  # Needed for the case where every sample is unique for all subjects
  first_par_list <- split(first_par, seq(NROW(first_par)))
  # Get unique pars (new accepted particles) and check length for
  # all subjects is greater than unq_vals
  n_unique_sub <- lapply(lapply(first_par_list, unique), length)
  if (all(n_unique_sub > n_unique)) {
    message("Enough unique values detected: ", n_unique)
    message("Testing proposal distribution creation")
    attempt <- tryCatch({
      if(n_cores > 1){
        mclapply(X = 1:pmwgs$n_subjects,FUN = conditionals_func,samples = test_samples, n_pars, i, mc.cores = n_cores)
      } else{
        lapply(X = 1:pmwgs$n_subjects,FUN = conditionals_func,samples = test_samples, n_pars, i)
      }
    },error=function(e) e, warning=function(w) w)
    if (any(class(attempt) %in% "warning")) {
      warning("An problem was encountered creating proposal distribution")
      warning("Increasing required unique values and continuing adaptation")
      return("increase")
    }
    else {
      message("Successfully adapted after ", i, "iterations - stopping early")
      return("success")
    }
  }
  return("continue")
}


