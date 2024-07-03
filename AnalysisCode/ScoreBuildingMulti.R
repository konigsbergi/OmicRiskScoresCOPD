#load necessary R packages
library(caret)
library(glmnet)
library(glmnetUtils)
library(dplyr)
library(doParallel)
library(multiview)

#set seed for reproducibility
set.seed(11)

#load phenotype data
pd = read.csv("/projects/ikonigsberg@xsede.org/ors_topmed/pd_main_overlap.csv", row.names = 1)

#categorize variables to model
num_vars = c("AWT_seg_Thirona_P2", "distwalked_P2", "FEV1_post_P2", "lung_density_vnb_P2")
log_vars = c("pctEmph_Thirona_P2")

#select variables of interest
pd = pd %>% dplyr::select(sid, all_of(num_vars), all_of(log_vars))

#load processed omics data
prot = readRDS("/projects/ikonigsberg@xsede.org/ors_topmed/ProcessedProteins.rds")
met = readRDS("/projects/ikonigsberg@xsede.org/ors_topmed/ProcessedMetabolites.rds")
met = met %>% dplyr::select(!where(~all(.x %in% 0:1))) #remove binary features. 
cts = readRDS("/projects/ikonigsberg@xsede.org/ors_topmed/ProcessedNormalizedCounts.rds")

rownames(pd) = pd$sid
pd = pd[order(rownames(pd)), ]

prot = prot[rownames(prot) %in% rownames(pd), ] 
met = met[rownames(met) %in% rownames(pd), ]
cts = cts[, colnames(cts) %in% rownames(pd), ]

#log transform proteins and metabolites
prot = log(prot + 1)
met = log(met + 1)

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

#define rhos to test
rhos = c(0, 0.5, 1)

#models for num_vars
out_num = sapply(num_vars, function(a){
  #a = num_vars[1]
  print(paste0("Starting on: ", a, "; ", Sys.time()))
  pd_sub = pd[pd$sid %in% train_ids, ]
  pd_sub = pd_sub[!is.na(eval(parse(text = paste0("pd_sub$", a)))), ]
  one = as.matrix(prot[rownames(prot) %in% rownames(pd_sub), ])
  two = as.matrix(cts[rownames(cts) %in% rownames(pd_sub), ])
  three = as.matrix(met[rownames(met) %in% rownames(pd_sub), ])
  pheno = eval(parse(text = paste0("pd_sub$", a)))
  cors = sapply(rhos, function(b){
    print(b)
    set.seed(11)
    mod = cv.multiview(list(one, two, three), pheno, nlambda = 20, type.measure = "mse", family = gaussian(), alpha = 1, rho = b)
    saveRDS(mod, paste0("/scratch/alpine/ikonigsberg@xsede.org/multi/rho_", b,"_model_", a, ".rds"))
    set.seed(11)
    pheno_test = pd[pd$sid %in% test_ids, ]
    one_test = as.matrix(prot[rownames(prot) %in% rownames(pheno_test), ])
    two_test = as.matrix(cts[rownames(cts) %in% rownames(pheno_test), ])
    three_test = as.matrix(met[rownames(met) %in% rownames(pheno_test), ])
    pheno_test = eval(parse(text = paste0("pheno_test$", a)))     
    pred = predict(mod, newx = list(one_test, two_test, three_test), s = mod$lambda.1se) 
    pred_df = as.data.frame(cbind(pred, pheno_test))
    saveRDS(pred_df, paste0("/scratch/alpine/ikonigsberg@xsede.org/multi/rho_", b,"_predictions_", a, ".rds"))
    cory = cor(pred_df$s1, pred_df$pheno_test, use = "na.or.complete")
    return(cory)
    })
      print(paste0(a, " is done! ", Sys.time()))
      return(cors)
})

out_num

#models for log vars
out_log = sapply(log_vars, function(a){
  #a = log_vars[1]
  print(paste0("Starting on: ", a, "; ", Sys.time()))
  pd_sub = pd[pd$sid %in% train_ids, ]
  pd_sub = pd_sub[!is.na(eval(parse(text = paste0("pd_sub$", a)))), ]
  one = as.matrix(prot[rownames(prot) %in% rownames(pd_sub), ])
  two = as.matrix(cts[rownames(cts) %in% rownames(pd_sub), ])
  three = as.matrix(met[rownames(met) %in% rownames(pd_sub), ])
  pheno = log(eval(parse(text = paste0("pd_sub$", a))))
  cors = sapply(rhos, function(b){
    print(b)
    set.seed(11)
    mod = cv.multiview(list(one, two, three), pheno, nlambda = 20, type.measure = "mse", family = gaussian(), alpha = 1, rho = b)
    saveRDS(mod, paste0("/scratch/alpine/ikonigsberg@xsede.org/multi/rho_", b,"_model_", a, ".rds"))
    set.seed(11)
    pheno_test = pd[pd$sid %in% test_ids, ]
    one_test = as.matrix(prot[rownames(prot) %in% rownames(pheno_test), ])
    two_test = as.matrix(cts[rownames(cts) %in% rownames(pheno_test), ])
    three_test = as.matrix(met[rownames(met) %in% rownames(pheno_test), ])
    pheno_test = log(eval(parse(text = paste0("pheno_test$", a))))  
    pred = predict(mod, newx = list(one_test, two_test, three_test), s = mod$lambda.1se) 
    pred_df = as.data.frame(cbind(pred, pheno_test))
    saveRDS(pred_df, paste0("/scratch/alpine/ikonigsberg@xsede.org/multi/rho_", b,"_predictions_", a, ".rds"))
    cory = cor(pred_df$s1, pred_df$pheno_test, use = "na.or.complete")
    return(cory)
    })
      print(paste0(a, " is done! ", Sys.time()))
      return(cors)
})

out_log


