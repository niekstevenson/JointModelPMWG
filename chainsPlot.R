library(tidyr)
library(ggplot2)
df <- sampled$samples$theta_mu[,15000:sampled$samples$idx]
tmp<-as.data.frame(t(df))
tmp2 <- pivot_longer(tmp, cols = everything(),names_to = "Parameter", values_to = "value")
tmp2$Iteration = rep(seq(15000:sampled$samples$idx),each = sampled$n_pars)
#ggplot(tmp2, aes(x=Iteration, y=value))+geom_line()+facet_wrap(~Parameter, nrow = sampled$n_pars)+theme_bw()
ggplot(tmp2, aes(x=Iteration, y=value))+geom_line()+facet_wrap(~Parameter, scales = "free")+theme_bw()
#ggplot(tmp2, aes(x=Iteration, y=value, colour = Parameter))+geom_line()+theme_bw()