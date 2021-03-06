extract_samples <- function(sampler, stage = c("adapt", "sample")) {
  samples <- sampler$samples
  stage_filter <- samples$stage %in% stage
  sampled_filter <- seq_along(samples$stage) <= samples$idx
  
  list(
    group_mu = samples$group_mu[,, stage_filter & sampled_filter],
    group_var = samples$group_var[,,, stage_filter & sampled_filter],
    alpha = samples$alpha[, , stage_filter & sampled_filter]
  )
}

particle_draws <- function(n, mu, covar) {
  if (n <= 0) {
    return(NULL)
  }
  mvtnorm::rmvnorm(n, mu, covar)
}

update.epsilon<- function(epsilon2, acc, p, i, d, alpha) {
  c=((1-1/d)*sqrt(2*pi)*exp(alpha^2/2)/(2*alpha) + 1/(d*p*(1-p)))
  Theta=log(sqrt(epsilon2))
  Theta=Theta+c*(acc-p)/max(200, i/d)
  return(exp(Theta))
}

unwind <- function(var_matrix, ...) {
  y <- t(chol(var_matrix))
  diag(y) <- log(diag(y))
  y[lower.tri(y, diag = TRUE)]
}

wind <- function(var_vector, ...) {
  n <- sqrt(2 * length(var_vector) + 0.25) - 0.5 ## Dim of matrix.
  if ((n * n + n) != (2 * length(var_vector))) stop("Wrong sizes in unwind.")
  out <- array(0, dim = c(n, n))
  out[lower.tri(out, diag = TRUE)] <- var_vector
  diag(out) <- exp(diag(out))
  out %*% t(out)
}


extend_sampler <- function(sampler, n_samples, stage) {
  old <- sampler$samples
  par_names <- sampler$par_names
  subject_ids <- sampler$subjects
  groups <- sampler$groups
  start <- old$idx + 1
  end <- old$idx + n_samples
  
  new_group_mu <- array(NA_real_,dim = dim(old$group_mu) + c(0, 0, n_samples), dimnames = list(par_names, groups, NULL))
  new_group_mu[,, - (start:end)] <- old$group_mu
  sampler$samples$group_mu <- new_group_mu
  
  new_group_var <- array(NA_real_,dim = dim(old$group_var) + c(0, 0, 0, n_samples),dimnames = list(par_names, par_names, groups, NULL))
  new_group_var[,,, - (start:end)] <- old$group_var
  sampler$samples$group_var <- new_group_var
  
  new_group_a_half <- array(NA_real_,dim = dim(old$group_a_half) + c(0,0, n_samples), dimnames = list(par_names, groups, NULL))
  new_group_a_half[,, - (start:end)] <- old$group_a_half
  sampler$samples$group_a_half <- new_group_a_half
  
  new_epsilon <- array(NA_real_,dim = dim(old$epsilon) + c(0, n_samples),dimnames = list(subject_ids, NULL))
  new_epsilon[, - (start:end)] <- old$epsilon
  sampler$samples$epsilon <- new_epsilon
  
  new_orig <- array(NA_real_,dim = dim(old$origin) + c(0, n_samples),dimnames = list(subject_ids, NULL))
  new_orig[, - (start:end)] <- old$origin
  sampler$samples$origin <- new_orig
  
  new_tmu <- array(NA_real_,dim = dim(old$theta_mu) + c(0, n_samples), dimnames = list(par_names, NULL))
  new_tmu[, - (start:end)] <- old$theta_mu
  sampler$samples$theta_mu <- new_tmu
  
  new_tvar <- array(NA_real_,dim = dim(old$theta_var) + c(0, 0, n_samples),dimnames = list(par_names, par_names, NULL))
  new_tvar[, , - (start:end)] <- old$theta_var
  sampler$samples$theta_var <- new_tvar
  
  new_alph <- array(NA_real_,dim = dim(old$alpha) + c(0, 0, n_samples),dimnames = list(par_names, subject_ids, NULL))
  new_alph[, , - (start:end)] <- old$alpha
  sampler$samples$alpha <- new_alph
  
  new_sll <- array(NA_real_,dim = dim(old$subj_ll) + c(0, n_samples),dimnames = list(subject_ids, NULL))
  new_sll[, - (start:end)] <- old$subj_ll
  sampler$samples$subj_ll <- new_sll
  
  new_ahalf <- array(NA_real_,dim = dim(old$a_half) + c(0, n_samples), dimnames = list(par_names, NULL))
  new_ahalf[, - (start:end)] <- old$a_half
  sampler$samples$a_half <- new_ahalf
  
  sampler$samples$stage <- c(old$stage, rep(stage, n_samples))
  sampler
}

trim_na <- function(sampler) {
  idx <- sampler$samples$idx
  sampler$samples$group_a_half <- sampler$samples$group_a_half[,,1:idx]
  sampler$samples$group_var <- sampler$samples$group_var[,,,1:idx]
  sampler$samples$group_mu <- sampler$samples$group_mu[,,1:idx]
  sampler$samples$theta_mu <- sampler$samples$theta_mu[, 1:idx]
  sampler$samples$theta_var <- sampler$samples$theta_var[, , 1:idx]
  sampler$samples$alpha <- sampler$samples$alpha[, , 1:idx]
  sampler$samples$subj_ll <- sampler$samples$subj_ll[, 1:idx]
  sampler$samples$a_half <- sampler$samples$a_half[, 1:idx]
  sampler$samples$stage <- sampler$samples$stage[1:idx]
  sampler$samples$epsilon <- sampler$samples$epsilon[,1:idx]
  sampler$samples$origin <- sampler$samples$origin[, 1:idx]
  sampler
}

relabel_samples <- function(sampler, indices, from="burn", to="adapt") {
  old_stage <- sampler$samples$stage
  if (!all(old_stage[indices] %in% from)) {
    stop(paste("Not all samples were from the", from, "stage"))
  }
  sampler$samples$stage[indices] <- to
  sampler
}


as_mcmc <- function(sampler, selection = "theta_mu", filter = stages) {
  stages <- c("burn", "adapt", "sample")
  if (all(filter %in% stages)) {
    filter <- which(sampler$samples$stage %in% filter)
  } else if (!all(filter %in% 1:sampler$samples$idx)) {
    stop("filter is not a vector of stage names, or integer vector of indices")
  }
  
  if (selection == "theta_mu") {
    return(coda::mcmc(t(sampler$samples$theta_mu[, filter])))
  } else if (selection == "theta_var") {
    tvar <- sampler$samples$theta_var[, , filter]
    return(stats::setNames(lapply(
      seq(dim(tvar)[1]),
      function(x) {
        coda::mcmc(t(tvar[x, , ]))
      }
    ), sampler$par_names))
  } else if (selection == "alpha") {
    alpha <- sampler$samples$alpha[, , filter]
    return(stats::setNames(lapply(
      seq(dim(alpha)[2]),
      function(x) {
        coda::mcmc(t(alpha[, x, ]))
      }
    ), sampler$subjects))
  }
  stop("Argument `selection` should be one of theta_mu, theta_var, alpha")
}

sample_store <- function(par_names, subject_ids, groups, iters = 1, stage = "init") {
  n_pars <- length(par_names)
  n_subjects <- length(subject_ids)
  n_groups <- length(groups)
  list(
    epsilon = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL)),
    origin = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL)),
    alpha = array(NA_real_,dim = c(n_pars, n_subjects, iters),dimnames = list(par_names, subject_ids, NULL)),
    group_mu = array(NA_real_,dim = c(n_pars, n_groups, iters), dimnames = list(par_names, groups, NULL)),
    group_var = array(NA_real_,dim = c(n_pars, n_pars, n_groups, iters), dimnames = list(par_names, par_names, groups, NULL)),
    theta_mu = array(NA_real_,dim = c(n_pars, iters), dimnames = list(par_names, NULL)),
    theta_var = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL)),
    stage = array(stage, iters),
    subj_ll = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL)),
    group_a_half = array(NA_real_,dim = c(n_pars, n_groups, iters),dimnames = list(par_names, groups, NULL)),
    a_half = array(NA_real_,dim = c(n_pars, iters),dimnames = list(par_names, NULL))
  )
}

unwind <- function(var_matrix, ...) {
  y <- t(chol(var_matrix))
  diag(y) <- log(diag(y))
  return(y[lower.tri(y, diag = TRUE)])
}


set_epsilon <- function(n_pars, epsilon) {
  if (checkmate::test_null(epsilon)) {
    if (n_pars > 15) {
      epsilon <- 0.1
    } else if (n_pars > 10) {
      epsilon <- 0.3
    } else {
      epsilon <- 0.5
    }
    message(
      sprintf(
        "Epsilon has been set to %.1f based on number of parameters",
        epsilon
      )
    )
  }
  epsilon
}

set_mix <- function(stage, mix) {
  if (checkmate::test_null(mix)) {
    if (stage == "sample") {
      mix <- c(0.1, 0.2, 0.7)
    } else {
      mix <- c(0.5, 0.5, 0.0)
    }
    message(
      sprintf(
        "mix has been set to c(%s) based on the stage being run",
        paste(mix, collapse = ", ")
      )
    )
  }
  mix
}

last_sample_group <- function(store, group) {
  list(
    tmu = store$group_mu[,group, store$idx],
    tvar = store$group_var[, , group,store$idx],
    alpha = store$alpha[, , store$idx],
    tvinv = store$last_group_var_inv[,,group],
    a_half = store$group_a_half[, group,store$idx],
    theta_mu_mean = store$theta_mu[, store$idx],
    theta_mu_invar = store$last_theta_var_inv
  )
}

last_sample_nested <- function(store) {
  list(
    tmu = store$theta_mu[, store$idx],
    tvar = store$theta_var[, , store$idx],
    group_mu = store$group_mu[, ,store$idx],
    tvinv = store$last_theta_var_inv,
    a_half = store$a_half[, store$idx]
  )
}

numbers_from_proportion <- function(mix_proportion, num_particles = 1000) {
  numbers <- stats::rbinom(n = 2, size = num_particles, prob = mix_proportion)
  if (mix_proportion[3] == 0) {
    numbers[3] <- 0
    numbers[2] <- num_particles - numbers[1]
  } else {
    numbers[3] <- num_particles - sum(numbers)
  }
  numbers
}


accept_rate <- function(pmwgs, window_size = 200) {
  n_samples <- pmwgs$samples$idx
  if (is.null(n_samples) || n_samples < 3) {
    return(array(0, dim(pmwgs$samples$alpha)[2]))
  }
  if (n_samples <= window_size) {
    start <- 1
    end <- n_samples
  } else {
    start <- n_samples - window_size + 1
    end <- n_samples
  }
  vals <- pmwgs$samples$alpha[1, , start:end]
  apply(
    apply(vals, 1, diff) != 0, # If diff != 0
    2,
    mean
  )
}