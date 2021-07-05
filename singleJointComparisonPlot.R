
load("~/RL-PMWG/samples/Joint/Joint_MSIT4v2B_RB8v_rlsat2V02B_RevStd.RData")
joint <- as.array(as_mcmc(sampled, filter = "sample"))

load("~/RL-PMWG/samples/MSIT4v2B.RData")
msit <- as.array(as_mcmc(sampled, filter = "sample"))

load("~/RL-PMWG/samples/RB8v.RData")
rb <- as.array(as_mcmc(sampled, filter = "sample"))

load("~/RL-PMWG/samples/rlsat2V02B.RData")
rlsat <- as.array(as_mcmc(sampled, filter = "sample"))

load("~/RL-PMWG/samples/RevStd.RData")
rev <- as.array(as_mcmc(sampled, filter = "sample"))

single <- cbind(msit, rb, rlsat, rev)
colnames(single) <- colnames(joint)

# use bayesplot::mcmc_intervals_data() function to get intervals data in format easy to pass to ggplot
library(bayesplot)
combined <- rbind(mcmc_intervals_data(single), mcmc_intervals_data(joint))
combined$model <- rep(c("Single", "Joint"), each = ncol(single))

library(ggplot2)
theme_set(bayesplot::theme_default())
pos <- position_nudge(y = ifelse(combined$model == "Model 2", 0, 0.1))
ggplot(combined, aes(x = m, y = parameter, color = model)) + 
  geom_linerange(aes(xmin = l, xmax = h), position = pos, size=2)+
  geom_linerange(aes(xmin = ll, xmax = hh), position = pos)+
  geom_point(position = pos)+
  coord_flip()+
  geom_vline(xintercept=0,linetype="dashed")