curExp <- "Rev"
curSub <- 2
exp <- experiments[[curExp]]
data <- exp$preppedData

data <- data[data$subject == curSub,]



# Predicted ---------------------------------------------------------------

parsPredict <- sampled$samples$alpha
rownames(parsPredict) <- gsub(".*[|]", "", rownames(parsPredict))
samples_stage <- length(sampled$samples$stage[sampled$samples$stage == "sample"])
iterations <- round(seq(from = (sampled$samples$idx - samples_stage),
                        to = sampled$samples$idx,
                        length.out = 100))


parsPredict <- parsPredict[exp$parNames, 1,iterations]

dfPredict <- data.frame()
for(i in 1:dim(parsPredict)[2]){
  tmp <- likelihood.RD(parsPredict[,i], data, sample = T)
  tmp$reps <- i
  dfPredict <- rbind(dfPredict, tmp)
}
dfPredict$subject <- curSub


pp_predict <- list()
pp_predict[[curSub]] <- dfPredict
if(RL = )
pp_predict <- addStimSetInfo(x =curSub, input=pp_predict, orig_dat=exp$data, 
                    addColumns = NULL, RL = RL, match = exp$match, Rev = F)
pp_predict <- posteriorCalc(pp_predict, combns = exp$factors)


predictInfo <- list()
predictInfo$'2' <- pp_predict
predictInfo <- getDescriptives(predictInfo, combns = exp$factors)


# data --------------------------------------------------------------------
pp_data <- list()
pp_data[[curSub]] <- data
pp_data <- addStimSetInfo(x =curSub, input=pp_data, orig_dat=exp$data, 
                     addColumns = NULL, RL = F, match = exp$match, Rev = F)
pp_data <- posteriorCalc(pp_data, combns = exp$factors)


dataInfo <- list()
dataInfo$'2' <- pp_data
dataInfo <- getDescriptives(dataInfo, combns = exp$factors)



# hyperparameters ---------------------------------------------------------


parsHyper <- sampled$samples$theta_mu
rownames(parsHyper) <- gsub(".*[|]", "", rownames(parsHyper))
parsHyper <- parsHyper[exp$parNames,iterations]


dfMed <- data.frame()
for(i in 1:dim(parsPredict)[2]){
  tmp <- likelihood.RD(parsHyper[,i], data, sample = T)
  tmp$reps <- i
  dfMed <- rbind(dfMed, tmp)
}
dfMed$subject <- curSub

pp_med <- list()
pp_med[[curSub]] <- dfMed
pp_med <- addStimSetInfo(x =curSub, input=pp_med, orig_dat=exp$data, 
                     addColumns = NULL, RL = F, match = exp$match, Rev = F)
pp_med <- posteriorCalc(pp_med, combns = exp$factors)


medInfo <- list()
medInfo$'2' <- pp_med
medInfo <- getDescriptives(medInfo, combns = exp$factors)

path <- paste0("./figures/predict/", curExp, "_ExclSub", curSub)
dir.create(path, recursive = T)

# pmwg_ppPlot(predictInfo, dataInfo, factors = exp$factors, PDF = T, path = "./figures/predictRB2")
pmwg_jointPpPlot(ppInfo = medInfo, dataInfo = dataInfo, jointInfo = predictInfo, factors = exp$factors, 
                 PDF = T, path = path)


# rownames(pars) <- gsub(".*[|]", "", rownames(pars))
# pars <- pars[exp$parNames, 1]
# 
# dat <- likelihood.RD(pars, data, sample = T)
# 
# attr(dat, "RL") <- attr(data, "RL")
# attr(dat, "constants") <- attr(data, "constants")
# attr(dat, "match")<- attr(data, "match")
# attr(dat, "n_v") <- attr(data, "n_v")
# attr(dat, "transFunc") <- attr(data, "transFunc")
# 
# likelihood.RD(pars, data)
# likelihood.RD(pars, dat)

medParsHyper <- apply(sampled$samples$theta_mu, 1, median)
names(medParsHyper) <- gsub(".*[|]", "", names(medParsHyper))
medParsHyper <- medParsHyper[exp$parNames]

medParsPredict <- apply(sampled$samples$alpha[,1,], 1, median)
names(medParsPredict) <- gsub(".*[|]", "", names(medParsPredict))
medParsPredict <- medParsPredict[exp$parNames]


likelihood.RD(medParsPredict, data)
likelihood.RD(medParsHyper, data)


