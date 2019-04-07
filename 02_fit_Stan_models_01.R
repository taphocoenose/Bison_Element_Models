# IMPORT ARCHAEOLOGICAL BISON ELEMENT DATASETS AND FIT STAN
# MODELS TO ESTIMATE THE EFFECTS OF VARIOUS DIETARY UTILITY,
# TAPHONOMIC, AND SITE TYPE COVARIATES ON ELEMENT FREQUENCIES.

# Load libraries
library(parallel)
library(rstan)
library(loo)

# Load data
load("01_Out.RData")

###############################################
# PREPARE DATA FRAME FOR PARALLEL CHAIN FITS #
#############################################

# Model names
model_names <- list.files(path="Stan models")
model_names <- model_names[grepl(".stan", model_names)]
model_names <- model_names[grepl("m1", model_names)]
model_names <- paste0("Stan models/", model_names)
# Set number of models to fit
n_models <- length(model_names)
# Set number of datasets X models to fit
n_datamodels <- n_models*length(ModelData)
# Set chains per model
n_chains <- 4
# Set random seeds
seeds <- sort(rep(sample(1:1000, size=n_datamodels), n_chains))
# Data frame for lapply of model fits
model_fit_pars <- data.frame(stan_model=rep(sort(rep(model_names, n_chains)),
                                        length(ModelData)),
                             seed=seeds, chain_ID=rep(1:n_chains, n_datamodels),
                             model_ID=rep(sort(rep(1:n_models, n_chains)),
                                      length(ModelData)),
                             DatasetID=sort(rep(1:length(ModelData), 
                                                n_chains*n_models)),
                             stringsAsFactors=FALSE)
model_fit_pars$datamodels <- paste0(model_fit_pars$model_ID, "_",
                                    model_fit_pars$DatasetID)
datamodels <- unique(model_fit_pars$datamodels)

# Create cluster and export RStan and LOO
cl <- makeCluster(detectCores()-1)
clusterEvalQ(cl, library(rstan))
clusterEvalQ(cl, library(loo))


###############################################
######## COMPILE AND FIT MODELS ##############
#############################################


# Export variables to cluster
clusterExport(cl, varlist=c("model_fit_pars", "ModelData"))

# Compile Stan models
comp_mods <- parLapply(cl, model_names,  function(x){
  # Allow Stan to save compiled model for use on each core
  rstan_options(auto_write=TRUE)
  # Compile model
  stan(x, data=ModelData[[1]], chains=1, iter=1)
})

# Export variables to cluster
clusterExport(cl, "comp_mods")


# Fit Stan models
sampled_chains <- parLapply(cl, 1:nrow(model_fit_pars), function(x){
  
  # Create sub directory for stan output
  subdir <- "Stan models/Out"
  if(!dir.exists(subdir)) dir.create(subdir)
  
  # Retrieve data
  md <- ModelData[[model_fit_pars$DatasetID[x]]]
  
  # Sink output to text file
  m_fname <- model_fit_pars$stan_model[x]
  m_fname <- paste0(substr(m_fname, 13, nchar(m_fname)-5), 
                    md$model_name[1])
  sink(paste0(subdir, "/out_", m_fname, ".txt"), append=TRUE)
  
  require(rstan)
  
  # Retrieve model x ID, chain ID, and unique seed for model x
  mID <- comp_mods[[model_fit_pars$model_ID[x]]]
  cID <- model_fit_pars$chain_ID[x]
  s <- model_fit_pars$seed[x]
  
  # Fit model x
  m <- stan(fit=mID, data=md, chains=1, chain_id=cID,
            seed=s, iter=1e4, warmup=7500, cores=1, 
            verbose=TRUE, save_warmup=FALSE,
            control=list(max_treedepth=13,
                         adapt_delta=0.95))
  
  # End sink
  sink()
  
  # Return model
  return(m)
  
})


# End previous cluster and create new one that uses
# fewer cores to avoid memory issues
stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(ceiling(detectCores()/4))

# Export variables and packages to cluster
clusterExport(cl, varlist=c("sampled_chains", "model_names",
                            "model_fit_pars", "ModelData"))
clusterEvalQ(cl, library(rstan))
clusterEvalQ(cl, library(loo))

# Aggregate chains into model fits
Model_List <- parLapply(cl, datamodels, function(i){
  
  chainInds <- which(model_fit_pars$datamodels==i)
  modelID <- model_fit_pars$model_ID[chainInds][1]
  dataID <-  model_fit_pars$DatasetID[chainInds][1]
  dataName <- ModelData[[dataID]]$model_name
  modelName <- model_fit_pars$stan_model[chainInds][1]
  modelName <- substr(modelName, 13, nchar(modelName)-5)
  
  # Aggregate chains for model i
  model_i <- sflist2stanfit(sampled_chains[chainInds])
  
  # LOO for model i
  loo_i <- loo(model_i)
  # Extract posterior samples from model i
  i_samples <- extract(model_i)
  
  # Summary of model i
  s_i <- summary(model_i)$summary
  
  # Creat object to store LOO model comparison
  LOO_comp <- data.frame(model=modelName,
                         data=dataName,
                         min_neff=min(s_i[,9], na.rm=TRUE),
                         med_neff=median(s_i[,9], na.rm=TRUE),
                         max_rhat=max(s_i[,10], na.rm=TRUE),
                         n_divergent=get_num_divergent(model_i),
                         max_treedepth=get_num_max_treedepth(model_i),
                         elpd=loo_i$elpd_loo,
                         elpd_SE=loo_i$se_elpd_loo,
                         e_pars=loo_i$p_loo,
                         e_pars_SE=loo_i$se_p_loo,
                         IC=loo_i$looic,
                         IC_SE=loo_i$se_looic,
                         stckng_wt=NA,
                         Par_k_d_negInf_0.5=length(which(loo_i$diagnostics$pareto_k<0.5)),
                         Par_k_d_0.5_0.7=length(which(loo_i$diagnostics$pareto_k>0.5 &
                                                      loo_i$diagnostics$pareto_k<0.7)),
                         Par_k_d_0.7_1=length(which(loo_i$diagnostics$pareto_k>0.7 &
                                                    loo_i$diagnostics$pareto_k<1)),
                         Par_k_d_1_Inf=length(which(loo_i$diagnostics$pareto_k>1)),
                         StringsAsFactors=FALSE)

  # Pointwise log likelihood
  lpd_pt_matrix <- loo_i$pointwise[,1]
  
  return(list(model_name=paste0(modelName, dataName), 
              fitted_model=model_i, 
              parameters_summary=s_i, 
              model_summary=LOO_comp,
              loo=loo_i,
              lpd=lpd_pt_matrix, 
              post_samples=i_samples))
  
})

# End cluster
stopCluster(cl)

# Extract data frame of LOO model comparison
LOO_comp <- lapply(Model_List, function(j) j$model_summary)
LOO_comp <- do.call("rbind", LOO_comp)
# Extract log likelihood matrix and generate stacking weights
lpd_pt_matrix <- sapply(Model_List, function(j) j$lpd)
if(nrow(LOO_comp)>1) LOO_comp$stckng_wt <- stacking_weights(lpd_pt_matrix)


# Save output of LOO
save(LOO_comp, file="02_m1_Out_LOO.RData")

# Save Models List in chunks
i1 <- seq(1, 101, 2)[seq(1, 101, 2) %in% 1:length(ModelData)]
i2 <- seq(2, 102, 2)[seq(2, 102, 2) %in% 1:length(ModelData)]
if(length(i1)>length(i2)) i2 <- c(i2, i1[length(i1)])
for(i in 1:(ceiling(length(ModelData)/2))){
  
  # Surbset models
  output_models <- Model_List[i1[i]:i2[i]]
  
  out_file <- paste0("02_m1_Out_Models_0", i, ".RData")
  
  save(output_models, file=out_file)
  
  rm(output_models, out_file)
  
}
