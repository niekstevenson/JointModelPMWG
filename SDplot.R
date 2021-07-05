
load("~/RL-PMWG/samples/Joint/Joint_MSIT4v2B_RB8v_rlsat2V02B_RevStd.RData")
jointParSDs <- apply(as_mcmc(sampled, filter = "sample"), 2, sd) 

load("~/RL-PMWG/samples/MSIT4v2B.RData")
msitParSDs <- apply(as_mcmc(sampled, filter = "sample"), 2, sd) 

load("~/RL-PMWG/samples/RB8v.RData")
rbParSDs <- apply(as_mcmc(sampled, filter = "sample"), 2, sd) 

load("~/RL-PMWG/samples/rlsat2V02B.RData")
rlsatParSDs <- apply(as_mcmc(sampled, filter = "sample"), 2, sd) 

load("~/RL-PMWG/samples/RevStd.RData")
revParSDs <- apply(as_mcmc(sampled, filter = "sample"), 2, sd) 

jointParSDs <- c(jointParSDs, msitParSDs, rbParSDs, rlsatParSDs, revParSDs)
df <- data.frame(SD = jointParSDs, group = c(rep("joint", length(jointParSDs)/2), rep("single", length(jointParSDs)/2)), 
                                                name = rep(names(jointParSDs)[1:(length(jointParSDs)/2)], 2))

ggplot(data=df, aes(x=name, y=SD, fill=group)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



load("~/RL-PMWG/samples/Joint/Joint_MSIT4v2B_RB8v_rlsat2V02B_RevStd.RData")
jointParSDs <- rowMeans(pmwg_SD(sampled))

load("~/RL-PMWG/samples/MSIT4v2B.RData")
msitParSDs <- rowMeans(pmwg_SD(sampled))

load("~/RL-PMWG/samples/RB8v.RData")
rbParSDs <- rowMeans(pmwg_SD(sampled)) 

load("~/RL-PMWG/samples/rlsat2V02B.RData")
rlsatParSDs <- rowMeans(pmwg_SD(sampled))

load("~/RL-PMWG/samples/RevStd.RData")
revParSDs <- rowMeans(pmwg_SD(sampled))

jointParSDs <- c(jointParSDs, msitParSDs, rbParSDs, rlsatParSDs, revParSDs)
df <- data.frame(SD = jointParSDs, group = c(rep("joint", length(jointParSDs)/2), rep("single", length(jointParSDs)/2)), 
                 name = rep(names(jointParSDs)[1:(length(jointParSDs)/2)], 2))

ggplot(data=df, aes(x=name, y=SD, fill=group)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
