Read me Reilly:
Folders you'll need to make yourself: 
samples, figures and samples/Joint

Structure of running:
1. Set up experiments list, has a specific format that could do with some user friendliness
2. Send experiments to prepDataJoint (can also be a single experiment).
prepDataJoint can filter based on accuracy, rt and whether they met these criteria in the other data-sets you included (if multiple)
it includes the filtered data in each list entry of the experiments list. 

Then you can send it to run() or jointRun().
they perform similar roles, but jointRun integrates all data in one large dataset and makes the parameters identifiable as belonging to the ll function you want to apply to this single datasets. The jointLL makes sure these are all handled. Then for both:
pmwg runs the likelihood function you specified. The prepPars function works for every likelihood you specify. It creates a design matrix based on the parameter names you set in the experiments list. 

It saves a pmwg object with all samples information together with the experiments list that includes all sampling settings and parameter settings. This makes sure that posterior checks can be easily applied. 

I created a multitude of functions you can use to check your results. They are mostly oriented towards joint models or comparing joint models to results from single models but they are pretty modular and can be easily adjusted. For jointModels you can do jointPostCheck, for single models postCheck. There's also a parameterRecovery (for now only for joint models) and a profilePlot function, and I think more.

The goal was to be able to run and compare different joint models with only a few lines of code that differ between them. 

Some stuff is still a bit ugly, but that's mostly the RL stuff, which I needed. Also the RD.R (in models), contains a lot of transformation functions I tried for the MSIT, but other's won't use. 

