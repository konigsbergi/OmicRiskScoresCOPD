#load necessary R packages
library(caret)
library(glmnet)
library(glmnetUtils)
library(dplyr)
library(doParallel)
library(MASS) #for negative binomial modeling

#set seed for reproducibility
set.seed(11)

#register backend for parallelization
registerDoParallel(16)

#load phenotype data
pd = read.csv("/projects/ikonigsberg@xsede.org/ors_topmed/pd_main_overlap.csv", row.names = 1)

#Remove PRISM from GOLD definition
pd[!is.na(pd$finalGold_P2) & pd$finalGold_P2 == -1, ]$finalGold_P2 = NA
pd$GoldCategory_P2 = as.factor(pd$finalGold_P2)

#Define COPD by FEV1/FVC
pd$copd_P2 = NA
pd[!is.na(pd$FEV1_FVC_post_P2) & pd$FEV1_FVC_post_P2 < .7, ]$copd_P2 = 1
pd[!is.na(pd$FEV1_FVC_post_P2) & pd$FEV1_FVC_post_P2 > .7, ]$copd_P2 = 0

#Define additional change variables
pd$Change_P2_P3_lung_density_vnb = pd$lung_density_vnb_P3 - pd$lung_density_vnb_P2

#categorize variables to model
num_vars = c("Change_P2_P3_distwalked", "Change_P2_P3_FEV1_ml", "Change_P2_P3_FEV1_ml_yr", "Change_P2_P3_FEV1pp", "FEV1_post_P2", 
             "AWT_seg_Thirona_P2", "distwalked_P2", "FEV1_FVC_post_P2", "FEV1_post_P2", "FVC_post_P2", "lung_density_vnb_P2", 
             "Resting_SaO2_P2", "finalGold_P2", "Change_P2_P3_lung_density_vnb", "FEV6_post_P2", "Age_P2", "FEV1pp_post_P2", 
	     "FVCpp_post_P2", "PEF_post_P2", "FRC_TLC_ratio_Thirona_P2", "WallAreaPct_seg_Thirona_P2", "DLco_raw_P2", 
             "DLCO_GLI_pp_PbHb_adj_P2", "FEF2575_post_P2")
log_vars = c("pctEmph_Thirona_P2", "pctGasTrap_Thirona_P2", "BMI_P2", "ATS_PackYears_P2", "TLC_Thirona_P2", "FRC_Thirona_P2", "Pi10_Thirona_P2")
negbinom_vars = c("Exacerbation_Frequency_P2", "SGRQ_scoreTotal_P2", "BODE_P2")
bin_vars = c("Chronic_Bronchitis_P2", "Severe_Exacerbations_P2", "copd_P2", "smoking_status_P2") 
cat_vars = c("GoldCategory_P2", "MMRCDyspneaScor_P2") 

#select variables of interest
pd = pd %>% dplyr::select(sid, all_of(num_vars), all_of(log_vars), all_of(bin_vars), all_of(cat_vars), all_of(negbinom_vars))

#load processed omics data
prot = readRDS("/projects/ikonigsberg@xsede.org/ors_topmed/ProcessedProteins.rds")
met = readRDS("/projects/ikonigsberg@xsede.org/ors_topmed/ProcessedMetabolites.rds")
met = met %>% dplyr::select(!where(~all(.x %in% 0:1))) #remove binary features. 
cts = readRDS("/projects/ikonigsberg@xsede.org/ors_topmed/ProcessedNormalizedCounts.rds")

#log transform proteins and metabolites
prot = log(prot + 1)
met = log(met + 1)

rownames(pd) = pd$sid
pd = pd[order(rownames(pd)), ]

prot = prot[rownames(prot) %in% rownames(pd), ] 
met = met[rownames(met) %in% rownames(pd), ]
cts = cts[,colnames(cts) %in% rownames(pd)]

prot = prot[order(rownames(prot)), ]
met = met[order(rownames(met)), ]
cts = cts[, order(colnames(cts)) ]
cts = t(cts)

all(rownames(prot) == rownames(pd))
all(rownames(met) == rownames(pd))
all(rownames(cts) == rownames(pd))

#load training:testing ids
train_ids = readRDS("/projects/ikonigsberg@xsede.org/ors_topmed/train_ids_80_20.rds")
test_ids = readRDS("/projects/ikonigsberg@xsede.org/ors_topmed/test_ids_80_20.rds")

#protein models
lapply(num_vars, function(a){
  pd_sub = pd[pd$sid %in% train_ids, ]
  pd_sub = pd_sub[!is.na(eval(parse(text = paste0("pd_sub$", a)))), ]
  x = as.matrix(prot[rownames(prot) %in% rownames(pd_sub), ])
  y = eval(parse(text = paste0("pd_sub$", a)))
  set.seed(11)
  mod = cva.glmnet(x, y, type.measure = "mse", family = "gaussian", trace.it = T, na.action = na.omit, parallel = T, nlambda = 20)
  set.seed(11)
  saveRDS(mod, paste0("/scratch/alpine/ikonigsberg@xsede.org/ors_copd_final/model_protein_", a, ".rds"))
})

lapply(log_vars, function(a){
  pd_sub = pd[pd$sid %in% train_ids, ]
  pd_sub = pd_sub[!is.na(eval(parse(text = paste0("pd_sub$", a)))), ]
  x = as.matrix(prot[rownames(prot) %in% rownames(pd_sub), ])
  y = log(eval(parse(text = paste0("pd_sub$", a))))
  set.seed(11)
  mod = cva.glmnet(x, y, type.measure = "mse", family = "gaussian", trace.it = T, na.action = na.omit, parallel = T, nlambda = 20)
  set.seed(11)
  saveRDS(mod, paste0("/scratch/alpine/ikonigsberg@xsede.org/ors_copd_final/model_protein_", a, ".rds"))
})

lapply(negbinom_vars, function(a){
  pd_sub = pd[pd$sid %in% train_ids, ]
  pd_sub = pd_sub[!is.na(eval(parse(text = paste0("pd_sub$", a)))), ]
  x = as.matrix(prot[rownames(prot) %in% rownames(pd_sub), ])
  y = eval(parse(text = paste0("pd_sub$", a)))
  set.seed(11)
  mod = cva.glmnet(x, y, type.measure = "mse", family = negative.binomial(theta = 5), trace.it = T, na.action = na.omit, parallel = T, nlambda = 20)
  set.seed(11)
  saveRDS(mod, paste0("/scratch/alpine/ikonigsberg@xsede.org/ors_copd_final/model_protein_", a, ".rds"))
})

lapply(bin_vars, function(a){
  pd_sub = pd[pd$sid %in% train_ids, ]
  pd_sub = pd_sub[!is.na(eval(parse(text = paste0("pd_sub$", a)))), ]
  x = as.matrix(prot[rownames(prot) %in% rownames(pd_sub), ])
  y = eval(parse(text = paste0("pd_sub$", a)))
  y = as.factor(y)
  set.seed(11)
  mod = cva.glmnet(x, y, type.measure = "mse", family = "binomial", trace.it = T, na.action = na.omit, parallel = T, nlambda = 20)
  set.seed(11)
  saveRDS(mod, paste0("/scratch/alpine/ikonigsberg@xsede.org/ors_copd_final/model_protein_", a, ".rds"))
})

lapply(cat_vars, function(a){
  pd_sub = pd[pd$sid %in% train_ids, ]
  pd_sub = pd_sub[!is.na(eval(parse(text = paste0("pd_sub$", a)))), ]
  x = as.matrix(prot[rownames(prot) %in% rownames(pd_sub), ])
  y = eval(parse(text = paste0("pd_sub$", a)))
  y = as.factor(y)
  set.seed(11)
  mod = cva.glmnet(x, y, type.measure = "mse", family = "multinomial", trace.it = T, na.action = na.omit, parallel = T, nlambda = 20)
  set.seed(11)
  saveRDS(mod, paste0("/scratch/alpine/ikonigsberg@xsede.org/ors_copd_final/model_protein_", a, ".rds"))
})

###
#metabolite models
lapply(num_vars, function(a){
  pd_sub = pd[pd$sid %in% train_ids, ]
  pd_sub = pd_sub[!is.na(eval(parse(text = paste0("pd_sub$", a)))), ]
  x = as.matrix(met[rownames(met) %in% rownames(pd_sub), ])
  y = eval(parse(text = paste0("pd_sub$", a)))
  set.seed(11)
  mod = cva.glmnet(x, y, type.measure = "mse", family = "gaussian", trace.it = T, na.action = na.omit, parallel = T, nlambda = 20)
  set.seed(11)
  saveRDS(mod, paste0("/scratch/alpine/ikonigsberg@xsede.org/ors_copd_final/model_metabolite_", a, ".rds"))
})

lapply(log_vars, function(a){
  pd_sub = pd[pd$sid %in% train_ids, ]
  pd_sub = pd_sub[!is.na(eval(parse(text = paste0("pd_sub$", a)))), ]
  x = as.matrix(met[rownames(met) %in% rownames(pd_sub), ])
  y = log(eval(parse(text = paste0("pd_sub$", a))))
  set.seed(11)
  mod = cva.glmnet(x, y, type.measure = "mse", family = "gaussian", trace.it = T, na.action = na.omit, parallel = T, nlambda = 20)
  set.seed(11)
  saveRDS(mod, paste0("/scratch/alpine/ikonigsberg@xsede.org/ors_copd_final/model_metabolite_", a, ".rds"))
})

lapply(negbinom_vars, function(a){
  pd_sub = pd[pd$sid %in% train_ids, ]
  pd_sub = pd_sub[!is.na(eval(parse(text = paste0("pd_sub$", a)))), ]
  x = as.matrix(met[rownames(met) %in% rownames(pd_sub), ])
  y = eval(parse(text = paste0("pd_sub$", a)))
  set.seed(11)
  mod = cva.glmnet(x, y, type.measure = "mse", family = negative.binomial(theta = 5), trace.it = T, na.action = na.omit, parallel = T, nlambda = 20)
  set.seed(11)
  saveRDS(mod, paste0("/scratch/alpine/ikonigsberg@xsede.org/ors_copd_final/model_metabolite_", a, ".rds"))
})

lapply(bin_vars, function(a){
  print(a)
  pd_sub = pd[pd$sid %in% train_ids, ]
  pd_sub = pd_sub[!is.na(eval(parse(text = paste0("pd_sub$", a)))), ]
  x = as.matrix(met[rownames(met) %in% rownames(pd_sub), ])
  y = eval(parse(text = paste0("pd_sub$", a)))
  y = as.factor(y)
  set.seed(11)
  mod = cva.glmnet(x, y, type.measure = "mse", family = "binomial", trace.it = T, na.action = na.omit, parallel = T, nlambda = 20)
  set.seed(11)
  saveRDS(mod, paste0("/scratch/alpine/ikonigsberg@xsede.org/ors_copd_final/model_metabolite_", a, ".rds"))
})

lapply(cat_vars, function(a){
  pd_sub = pd[pd$sid %in% train_ids, ]
  pd_sub = pd_sub[!is.na(eval(parse(text = paste0("pd_sub$", a)))), ]
  x = as.matrix(met[rownames(met) %in% rownames(pd_sub), ])
  y = eval(parse(text = paste0("pd_sub$", a)))
  y = as.factor(y)
  set.seed(11)
  mod = cva.glmnet(x, y, type.measure = "mse", family = "multinomial", trace.it = T, na.action = na.omit, parallel = T, nlambda = 20)
  set.seed(11)
  saveRDS(mod, paste0("/scratch/alpine/ikonigsberg@xsede.org/ors_copd_final/model_metabolite_", a, ".rds"))
})

###
#transcript models
lapply(num_vars, function(a){
  pd_sub = pd[pd$sid %in% train_ids, ]
  pd_sub = pd_sub[!is.na(eval(parse(text = paste0("pd_sub$", a)))), ]
  x = as.matrix(cts[rownames(cts) %in% rownames(pd_sub), ])
  y = eval(parse(text = paste0("pd_sub$", a)))
  set.seed(11)
  mod = cva.glmnet(x, y, type.measure = "mse", family = "gaussian", trace.it = T, na.action = na.omit, parallel = T, nlambda = 20)
  set.seed(11)
  saveRDS(mod, paste0("/scratch/alpine/ikonigsberg@xsede.org/ors_copd_final/model_transcript_", a, ".rds"))
})

lapply(log_vars, function(a){
  pd_sub = pd[pd$sid %in% train_ids, ]
  pd_sub = pd_sub[!is.na(eval(parse(text = paste0("pd_sub$", a)))), ]
  x = as.matrix(cts[rownames(cts) %in% rownames(pd_sub), ])
  y = log(eval(parse(text = paste0("pd_sub$", a))))
  set.seed(11)
  mod = cva.glmnet(x, y, type.measure = "mse", family = "gaussian", trace.it = T, na.action = na.omit, parallel = T, nlambda = 20)
  set.seed(11)
  saveRDS(mod, paste0("/scratch/alpine/ikonigsberg@xsede.org/ors_copd_final/model_transcript_", a, ".rds"))
})

lapply(negbinom_vars, function(a){
  pd_sub = pd[pd$sid %in% train_ids, ]
  pd_sub = pd_sub[!is.na(eval(parse(text = paste0("pd_sub$", a)))), ]
  x = as.matrix(cts[rownames(cts) %in% rownames(pd_sub), ])
  y = eval(parse(text = paste0("pd_sub$", a)))
  set.seed(11)
  mod = cva.glmnet(x, y, type.measure = "mse", family = negative.binomial(theta = 5), trace.it = T, na.action = na.omit, parallel = T, nlambda = 20)
  set.seed(11)
  saveRDS(mod, paste0("/scratch/alpine/ikonigsberg@xsede.org/ors_copd_final/model_transcript_", a, ".rds"))
})


lapply(bin_vars, function(a){
  pd_sub = pd[pd$sid %in% train_ids, ]
  pd_sub = pd_sub[!is.na(eval(parse(text = paste0("pd_sub$", a)))), ]
  x = as.matrix(cts[rownames(cts) %in% rownames(pd_sub), ])
  y = eval(parse(text = paste0("pd_sub$", a)))
  y = as.factor(y)
  set.seed(11)
  mod = cva.glmnet(x, y, type.measure = "mse", family = "binomial", trace.it = T, na.action = na.omit, parallel = T, nlambda = 20)
  set.seed(11)
  saveRDS(mod, paste0("/scratch/alpine/ikonigsberg@xsede.org/ors_copd_final/model_transcript_", a, ".rds"))
})

lapply(cat_vars, function(a){
  pd_sub = pd[pd$sid %in% train_ids, ]
  pd_sub = pd_sub[!is.na(eval(parse(text = paste0("pd_sub$", a)))), ]
  x = as.matrix(cts[rownames(cts) %in% rownames(pd_sub), ])
  y = eval(parse(text = paste0("pd_sub$", a)))
  y = as.factor(y)
  set.seed(11)
  mod = cva.glmnet(x, y, type.measure = "mse", family = "multinomial", trace.it = T, na.action = na.omit, parallel = T, nlambda = 20)
  set.seed(11)
  saveRDS(mod, paste0("/scratch/alpine/ikonigsberg@xsede.org/ors_copd_final/model_transcript_", a, ".rds"))
})
