
# Setup -------------------------------------------------------------------
rm(list=ls())
library(parallel)
library(pmwg)
library(dplyr)
setwd("~/RL-PMWG/")

source("utils/utils-PMWG.R")
source("utils/DataPrep-PMWG.R")
source("jointUtils/DataPrepJoint-PMWG.R")
source("jointUtils/jointUtils-PWMG.R")
source("models/RW/RD.R")
load("samples/MSIT4v2B.RData")



true_x <- rowMeans(pmwg_median(sampled))

test_data <- likelihood.RD(true_x, sampled$data, sample = T)
attr(test_data, "RL") <- attr(sampled$data, "RL")
attr(test_data, "constants") <- attr(sampled$data, "constants")
attr(test_data, "match")<- attr(sampled$data, "match")
attr(test_data, "n_v") <- attr(sampled$data, "n_v")


test_ll <- function(x, 
                    n_values = 9,
                    synth_data,
                    dims,
                    ll_func,
                    server = FALSE,...) {
  pars <- names(x)
  xtmp <- x
  sequence <- seq(from = -.22,
                  to = .22,
                  length.out = n_values)
  op <- par(mfrow = dims)
  par_likelihoods <- lapply(setNames(pars, pars), function(p){
    testvalues <- sequence + true_x[[p]]
    tmp <- unlist(lapply(testvalues, function (i){ 
      #for each test vlaue, it will apply the function where i = testvalue.
      xtmp[[p]] <- i
      ll_func(pars = xtmp, data = synth_data)
    }))
    
    plot(x = testvalues,
         y = tmp,
         type = "b",
         main = p,
         xlab = "log par values",
         ylab = "log-likelihood",
         ...)
    abline(v = true_x[[p]], col = "red")
    return(tmp)})
  par(op)
  return(par_likelihoods)
}

ll_plots <- test_ll(x = true_x,
                    n_values = 9,
                    synth_data = test_data,
                    dims = c(2, 2),
                    ll_func = likelihood.RD
) 