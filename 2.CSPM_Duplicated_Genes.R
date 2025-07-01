rm(list = ls())
library(Seurat)
library(dplyr)
library(GSVA)
library(mclust)
library(ggplot2)
library(tibble)
library(caret)
library(tidyverse)
library(EasyML)
library(pROC)
library(ggprism)
library(PRROC)
source('script/gene.info.R')
cols <- c("#ED6355","#41A98E","#00B4D8","#EFA63A","#AB3282","#AA9A59",
          "#625D9E","#0073C2FF", "#868686FF", "#003C67FF","#8F7700FF", "#3B3B3BFF",
          'pink','firebrick3','black')


# This pipeline compares the performance of the Cellular Senescence Prediction Model (CSPM) across diverse duplicated gene control strategies.
# Example data is available at https://doi.org/10.7303/syn68545972.



# ---------------A: Remove duplicated genes (n>=2)-----------------------------

# 1. Data Preprocessing-----------------------------------------------------

# Load single-cell data and revised senescence gene sets
geneset <- readRDS('data/geneset_remove_n2.rds')
CAF.seurat <- readRDS("data/CAF.seurat.rds") 
CAF.seurat <- CAF.seurat %>% ScaleData(features = rownames(CAF.seurat))
exp <- as.data.frame(CAF.seurat@assays[["RNA"]]@data)
g <- gene.info(rownames(exp))
g <- g[g$gene_biotype=='protein_coding',]
exp <- exp[rownames(exp) %in% unique(g$symbol),] %>% as.matrix()

# GSVA analysis
ss <- gsva(exp, geneset, method = 'gsva', 
           kcdf = 'Gaussian',
           abs.ranking = F,
           mx.diff = T)
ss <- t(ss)

# GMM-based classification
Group <- data.frame(row.names = rownames(ss))
for (i in 1:ncol(ss)) {
  set.seed(123)
  fit <- Mclust(ss[,i],G = 2)
  Group[,colnames(ss)[i]] <- ifelse(fit$classification==2,'High','Low')
  pd <- data.frame(Cluster=Group[,colnames(ss)[i]],
                   Score=ss[,i])
  print(ggplot(pd, aes(x = Score, fill = Cluster)) +
          geom_histogram(aes(y = after_stat(density)), bins = 30, alpha = 0.6) +  
          geom_density(alpha = 0.7) +  
          theme_minimal() +
          labs(fill = "Cluster",title = colnames(ss)[i]))
}

Group$Type <- apply(Group,1,\(x){
  if (length(unique(x))==1) {
    return(unique(x))
  }else{
    return('Unclassfied')
  }
})

# DEGs selection
dea_cc <- Group[Group$Type!='Unclassfied',]
CAF.seurat <- AddMetaData(CAF.seurat,Group$Type,col.name = "Type")
Idents(CAF.seurat) <- CAF.seurat$Type
dea_genes <- FindMarkers(CAF.seurat,ident.1 = 'High',ident.2 = 'Low',logfc.threshold = 0.5,min.pct = 0.1,
                         min.diff.pct = 0.25,only.pos = T) %>% filter(p_val_adj<0.05,rownames(.) %in% unique(g$symbol)) %>% rownames() 

data_genes <- readRDS('data/data_genes.rds')
intergene <- dea_genes[dea_genes %in% data_genes] #intersecting genes with validation cohorts

# Load ground-truth validation datasets (filtered intersecting genes)
load('data/CSPM_validation_n2')
dim(exp_public);dim(anno_all) #validation cohort #1
dim(exp_self);dim(anno_self) #validation cohort #2

# 2. CSPM construction----------------------------------------------------------

# Select Pseudo ground-truth cells 
exp_scale <-  as.data.frame(CAF.seurat@assays[["RNA"]]@scale.data)
mexp <- exp_scale[intergene,rownames(Group)[Group$Type!='Unclassfied']] %>% t() %>% as.data.frame()

# Partition data into training and test sets
set.seed(123)
texp <- mexp[caret::createDataPartition(dea_cc$Type,p = 0.8)[[1]],]
tcc <- Group[rownames(texp),]

vexp <- mexp[!rownames(mexp)%in%rownames(texp),]
vcc <- Group[rownames(vexp),]

tcc$Type <- ifelse(tcc$Type=='High',"x1","x0")
vcc$Type <- ifelse(vcc$Type=='High',"x1","x0")

data.features = texp
data.targets = factor(tcc$Type)

# Select machine-learning models
mods <- c("rpart", "partDSA", "BstLm", "glmboost", "ENet", "rFerns", "evtree", "svmRadial", "RF", "svmLinear3","LDA", "gbm")

# Ensembel model training
fit_ensemble <- ensemble.train(data.features = data.features,data.targets = data.targets,methods = mods)
names(fit_ensemble) <- mods
#saveRDS(fit_ensemble,'data/fit_ensemble_n2.rds')


# Training set accuracy
res_train <- ensemble.predict(fit_ensemble,data.features,outcome = "x1")
res_train$Ensemble.res <- ifelse(res_train$Ensemble.ratio>0.5,"x1","x0")
sum(res_train$Ensemble.res==as.character(tcc$Type))/nrow(tcc)

# Test set accuracy
res_test <- ensemble.predict(fit_ensemble,data.features=vexp,outcome = "x1")
res_test$Ensemble.res <- ifelse(res_test$Ensemble.ratio>0.5,"x1","x0")
sum(res_test$Ensemble.res==as.character(vcc$Type))/nrow(vcc)


# 3. CSPM validation-----------------------------------------------------------

## 3.1 CSPM performance in validation cohort #1 (external cohort)####
exp <- exp_public #expression profile
dim(exp);dim(anno_all) #ground-truth labels

vexp1 <-as.data.frame(t((exp[intergene,]))) 
vcc1 <- anno_all

# Predict senescence using CSPM
res_test <- ensemble.predict(fit_ensemble,data.features=vexp1,outcome = "x1")
res_test$Ensemble.res <- ifelse(res_test$Ensemble.ratio>0.5,"x1","x0")
test <- mutate(res_test,Type=vcc1$Type,treat=vcc1$`treatment:ch1`)

# Calculate classification metrics
confusion_matrix <- table(test$Type,test$Ensemble.res,dnn=c("Actual", "Predicted"));confusion_matrix
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
recall <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
f1_score <- 2 * precision * recall / (precision + recall)
G_mean <- sqrt(confusion_matrix[2, 2] / sum(confusion_matrix[2, ]) * confusion_matrix[1, 1] / sum(confusion_matrix[1, ]))

predicted_probs <- test$Ensemble.ratio;actual <- test$Type
pr <- pr.curve(predicted_probs[actual == 'x1'], predicted_probs[actual == 'x0'], curve = TRUE)
Average_Precision <- pr$auc.integral

roc <- roc(test$Type,test$Ensemble.ratio)
roc_auc <- auc(roc)

cat(" Accuracy: ", accuracy, "\n",
    "Precision: ", precision, "\n",
    "Recall: ", recall, "\n",
    "F1 Score: ", f1_score, "\n",
    'G-mean:',G_mean,'\n',
    'AUC:',roc_auc,'\n',
    'Average Precision:',Average_Precision)

# Save Validation cohort #1 (external cohort) perfromace
metric <- data.frame(Metric=c("Accuracy", "Precision", "Recall", "F1 Score",'G-mean','Average Precision','AUC'),
                     Value = round(c(accuracy, precision, recall, f1_score, G_mean, Average_Precision,roc_auc),3))
metric$Metric <- factor(metric$Metric,levels = c("Accuracy", "Precision", "Recall", "F1 Score",'G-mean','Average Precision','AUC'))
#saveRDS(metric,'data/metric_val1_n2.rds')


## 3.2 CSPM performance in validation cohort #2 (in-house cohort)####
exp <- exp_self #exoressuib profile
dim(exp);dim(anno_self) #ground-truth labels

vexp1 <-as.data.frame(t((exp[intergene,]))) 
vcc1 <- anno_self

# Predict senescence using CSPM
res_test <- ensemble.predict(fit_ensemble,data.features=vexp1,outcome = "x1")
res_test$Ensemble.res <- ifelse(res_test$Ensemble.ratio>0.5,"x1","x0")
sum(res_test$Ensemble.res==as.character(vcc1$Type))/nrow(vcc1)
test <- mutate(res_test,Type=vcc1$Type,treat=vcc1$`treatment:ch1`)

# classifier performance
confusion_matrix <- table(test$Type,test$Ensemble.res,dnn=c("Actual", "Predicted"));confusion_matrix
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
recall <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
f1_score <- 2 * precision * recall / (precision + recall)
G_mean <- sqrt(confusion_matrix[2, 2] / sum(confusion_matrix[2, ]) * confusion_matrix[1, 1] / sum(confusion_matrix[1, ]))

predicted_probs <- test$Ensemble.ratio
actual <- test$Type
pr <- pr.curve(predicted_probs[actual == 'x1'], predicted_probs[actual == 'x0'], curve = TRUE)
Average_Precision <- pr$auc.integral

roc <- roc(test$Type,test$Ensemble.ratio)
roc_auc <- auc(roc);roc_auc

cat(" Accuracy: ", accuracy, "\n",
    "Precision: ", precision, "\n",
    "Recall: ", recall, "\n",
    "F1 Score: ", f1_score, "\n",
    'G-mean:',G_mean,'\n',
    'AUC:',roc_auc,'\n',
    'Average Precision:',Average_Precision)


# save classifier metrics
metric <- data.frame(Metric=c("Accuracy", "Precision", "Recall", "F1 Score",'G-mean','Average Precision','AUC'),
                     Value = round(c(accuracy, precision, recall, f1_score, G_mean, Average_Precision,roc_auc),3))
metric$Metric <- factor(metric$Metric,levels = c("Accuracy", "Precision", "Recall", "F1 Score",'G-mean','Average Precision','AUC'))
#saveRDS(metric,'data/metric_val2_n2.rds')


## 3.3 CSPM performance in pseudo ground-truth labels across single-cell datasets####

# Load CAF.seurat of each single-cell dataset
CAF_list <- readRDS("data/CAF_list.rds")

cmp <- data.frame(row.names = names(CAF_list))
d <- list()
for(j in names(CAF_list)){
  CAF.seurat <- CAF_list[[j]]
  CAF.seurat <- CAF.seurat %>% ScaleData(features = rownames(CAF.seurat))
  
  # Identify pseudo ground-truth labels
  exp <- as.data.frame(CAF.seurat@assays[["RNA"]]@data)
  g <- gene.info(rownames(exp))
  g <- g[g$gene_biotype=='protein_coding',]
  names <- unique(g$symbol)
  exp <- exp[rownames(exp)%in%names,]
  exp <- as.matrix(exp)

  ss <- gsva(exp, 
             geneset, 
             method = 'gsva', 
             kcdf = 'Gaussian',
             abs.ranking = F,
             mx.diff = T)
  ss <- t(ss)

  Group <- data.frame(row.names = rownames(ss))
  for (i in 1:ncol(ss)) {
    set.seed(123)
    fit <- Mclust(ss[,i],G = 2)
    Group[,colnames(ss)[i]] <- ifelse(fit$classification==2,'High','Low')
    pd <- data.frame(Cluster=Group[,colnames(ss)[i]],
                     Score=ss[,i])
    print(ggplot(pd, aes(x = Score, fill = Cluster)) +
            geom_histogram(aes(y = after_stat(density)), bins = 30, alpha = 0.6) +  
            geom_density(alpha = 0.7) + 
            theme_minimal() +
            labs(fill = "Cluster",title = colnames(ss)[i]))
  }
  
  Group$Type <- apply(Group,1,\(x){
    if (length(unique(x))==1) {
      return(unique(x))
    }else{
      return('Unclassfied')
    }
  })

  dea_cc <- Group[Group$Type!='Unclassfied',] # pseudo ground-truth labels

  # Predict senescence using CSPM
  exp_scale <-  as.data.frame(CAF.seurat@assays[["RNA"]]@scale.data)
  mexp <- exp_scale[intergene,rownames(Group)[Group$Type!='Unclassfied']] %>% t() %>% as.data.frame()
  
  if(j=='CNP0004138'){
    # Using only test set in discovery cohort
    set.seed(123)
    vexp <- mexp[!rownames(mexp)%in%rownames(texp),]
    vcc <- Group[rownames(vexp),]
    vcc$Type <- ifelse(vcc$Type=='High',"x1","x0")
    
    res_test <- ensemble.predict(fit_ensemble,data.features=vexp,outcome = "x1")
    res_test$Ensemble.res <- ifelse(res_test$Ensemble.ratio>0.5,"x1","x0")

    # classifier metrics
    confusion_matrix <- table(vcc$Type,res_test$Ensemble.res,dnn=c("Actual", "Predicted"))
    
    accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
    precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
    recall <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
    f1_score <- 2 * precision * recall / (precision + recall)
    G_mean <- sqrt(confusion_matrix[2, 2] / sum(confusion_matrix[2, ]) * confusion_matrix[1, 1] / sum(confusion_matrix[1, ]))
    
    predicted_probs <- res_test$Ensemble.ratio
    actual <- vcc$Type
    pr <- pr.curve(predicted_probs[actual == 'x1'], predicted_probs[actual == 'x0'], curve = TRUE)
    Average_Precision <- pr$auc.integral
    
    roc <- roc(vcc$Type,res_test$Ensemble.ratio)
    roc_auc <- auc(roc)
    
    d[[j]]$vcc <- vcc
    d[[j]]$res_test <- res_test
    
    cmp[j,'Accuracy'] <- accuracy
    cmp[j,'Precision'] <- precision
    cmp[j,'Recall'] <- recall
    cmp[j,'F1 Score'] <- f1_score
    cmp[j,'G-mean'] <- G_mean
    cmp[j,'Average Precision'] <- Average_Precision
    cmp[j,'AUC'] <- roc_auc
  }else{
    # Using all the pseudo ground-truth labels in other datasets
    vexp <- mexp
    vcc <- Group[rownames(vexp),]
    vcc$Type <- ifelse(vcc$Type=='High',"x1","x0")
    
    res_test <- ensemble.predict(fit_ensemble,data.features=vexp,outcome = "x1")
    res_test$Ensemble.res <- ifelse(res_test$Ensemble.ratio>0.5,"x1","x0")
    
    # classifier metrics
    confusion_matrix <- table(vcc$Type,res_test$Ensemble.res,dnn=c("Actual", "Predicted"))
    
    accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
    precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
    recall <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
    f1_score <- 2 * precision * recall / (precision + recall)
    G_mean <- sqrt(confusion_matrix[2, 2] / sum(confusion_matrix[2, ]) * confusion_matrix[1, 1] / sum(confusion_matrix[1, ]))
    
    predicted_probs <- res_test$Ensemble.ratio
    actual <- vcc$Type
    pr <- pr.curve(predicted_probs[actual == 'x1'], predicted_probs[actual == 'x0'], curve = TRUE)
    Average_Precision <- pr$auc.integral
    
    roc <- roc(vcc$Type,res_test$Ensemble.ratio)
    roc_auc <- auc(roc)
    
    d[[j]]$vcc <- vcc
    d[[j]]$res_test <- res_test
    
    cmp[j,'Accuracy'] <- accuracy
    cmp[j,'Precision'] <- precision
    cmp[j,'Recall'] <- recall
    cmp[j,'F1 Score'] <- f1_score
    cmp[j,'G-mean'] <- G_mean
    cmp[j,'Average Precision'] <- Average_Precision
    cmp[j,'AUC'] <- roc_auc
  }
}


# Merge pseudo ground-truth labels across datasets
vcc <- do.call(rbind, lapply(d, function(x) x$vcc))
res_test <- do.call(rbind, lapply(d, function(x) x$res_test))

confusion_matrix <- table(vcc$Type,res_test$Ensemble.res,dnn=c("Actual", "Predicted"))
confusion_matrix

accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
recall <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
f1_score <- 2 * precision * recall / (precision + recall)
G_mean <- sqrt(confusion_matrix[2, 2] / sum(confusion_matrix[2, ]) * confusion_matrix[1, 1] / sum(confusion_matrix[1, ]))

predicted_probs <- res_test$Ensemble.ratio
actual <- vcc$Type
pr <- pr.curve(predicted_probs[actual == 'x1'], predicted_probs[actual == 'x0'], curve = TRUE)
Average_Precision <- pr$auc.integral

roc <- roc(vcc$Type,res_test$Ensemble.ratio)
roc_auc <- auc(roc)

cmp['Merge','Accuracy'] <- accuracy
cmp['Merge','Precision'] <- precision
cmp['Merge','Recall'] <- recall
cmp['Merge','F1 Score'] <- f1_score
cmp['Merge','G-mean'] <- G_mean
cmp['Merge','Average Precision'] <- Average_Precision
cmp['Merge','AUC'] <- roc_auc

# save CSPM performance (remove dupliated genes n>=2)
#saveRDS(cmp,'data/cmp_n2.rds')



# ------------------B: Remove duplicated genes (n>=3)--------------------------

# 1. Data Preprocessing-----------------------------------------------------

# Load single-cell data and revised senescence gene sets
geneset <- readRDS('data/geneset_remove_n3.rds')
CAF.seurat <- readRDS("data/CAF.seurat.rds") 
CAF.seurat <- CAF.seurat %>% ScaleData(features = rownames(CAF.seurat))
exp <- as.data.frame(CAF.seurat@assays[["RNA"]]@data)
g <- gene.info(rownames(exp))
g <- g[g$gene_biotype=='protein_coding',]
exp <- exp[rownames(exp) %in% unique(g$symbol),] %>% as.matrix()

# GSVA analysis
ss <- gsva(exp, 
           geneset, 
           method = 'gsva', 
           kcdf = 'Gaussian',
           abs.ranking = F,
           mx.diff = T)
ss <- t(ss)

# GMM-based classification
Group <- data.frame(row.names = rownames(ss))
for (i in 1:ncol(ss)) {
  set.seed(123)
  fit <- Mclust(ss[,i],G = 2)
  Group[,colnames(ss)[i]] <- ifelse(fit$classification==2,'High','Low')
  pd <- data.frame(Cluster=Group[,colnames(ss)[i]],
                   Score=ss[,i])
  print(ggplot(pd, aes(x = Score, fill = Cluster)) +
          geom_histogram(aes(y = after_stat(density)), bins = 30, alpha = 0.6) +  
          geom_density(alpha = 0.7) +  
          theme_minimal() +
          labs(fill = "Cluster",title = colnames(ss)[i]))
}

Group$Type <- apply(Group,1,\(x){
  if (length(unique(x))==1) {
    return(unique(x))
  }else{
    return('Unclassfied')
  }
})

# DEGs selection
dea_cc <- Group[Group$Type!='Unclassfied',]
CAF.seurat <- AddMetaData(CAF.seurat,Group$Type,col.name = "Type")
Idents(CAF.seurat) <- CAF.seurat$Type
dea_genes <- FindMarkers(CAF.seurat,ident.1 = 'High',ident.2 = 'Low',logfc.threshold = 0.5,min.pct = 0.1,
                         min.diff.pct = 0.25,only.pos = T) %>% filter(p_val_adj<0.05,rownames(.) %in% unique(g$symbol)) %>% rownames() 

data_genes <- readRDS('data/data_genes.rds')
intergene <- dea_genes[dea_genes %in% data_genes] #intersecting genes with validation cohorts

# Load ground-truth validation datasets (filtered intersecting genes)
load('data/CSPM_validation_n3.rda')
dim(exp_public);dim(anno_all) #validation cohort #1
dim(exp_self);dim(anno_self) #validation cohort #2


# 2. CSPM construction----------------------------------------------------------

# Select Pseudo ground-truth cells 
exp_scale <-  as.data.frame(CAF.seurat@assays[["RNA"]]@scale.data)
mexp <- exp_scale[intergene,rownames(Group)[Group$Type!='Unclassfied']] %>% t() %>% as.data.frame()

# Partition data into training and test sets
set.seed(123)
texp <- mexp[caret::createDataPartition(dea_cc$Type,p = 0.8)[[1]],]
tcc <- Group[rownames(texp),]

vexp <- mexp[!rownames(mexp)%in%rownames(texp),]
vcc <- Group[rownames(vexp),]

tcc$Type <- ifelse(tcc$Type=='High',"x1","x0")
vcc$Type <- ifelse(vcc$Type=='High',"x1","x0")

data.features = texp
data.targets = factor(tcc$Type)

# Select machine-learning models
mods <- c("rpart", "partDSA", "BstLm", "glmboost", "ENet", "rFerns", "evtree", "svmRadial", "RF", "svmLinear3","LDA", "gbm")

# Ensembel model training
fit_ensemble <- ensemble.train(data.features = data.features,data.targets = data.targets,methods = mods)
names(fit_ensemble) <- mods
#saveRDS(fit_ensemble,'data/fit_ensemble_n3.rds')

# Training set accuracy
res_train <- ensemble.predict(fit_ensemble,data.features,outcome = "x1")
res_train$Ensemble.res <- ifelse(res_train$Ensemble.ratio>0.5,"x1","x0")
sum(res_train$Ensemble.res==as.character(tcc$Type))/nrow(tcc)

# Test set accuracy
res_test <- ensemble.predict(fit_ensemble,data.features=vexp,outcome = "x1")
res_test$Ensemble.res <- ifelse(res_test$Ensemble.ratio>0.5,"x1","x0")
sum(res_test$Ensemble.res==as.character(vcc$Type))/nrow(vcc)


# 3. CSPM validation-----------------------------------------------------------

## 3.1 CSPM performance in validation cohort #1 (external cohort)####
exp <- exp_public #expression profile
dim(exp);dim(anno_all) #ground-truth labels

vexp1 <-as.data.frame(t((exp[intergene,]))) 
vcc1 <- anno_all

# Predict senescence using CSPM
res_test <- ensemble.predict(fit_ensemble,data.features=vexp1,outcome = "x1")
res_test$Ensemble.res <- ifelse(res_test$Ensemble.ratio>0.5,"x1","x0")
test <- mutate(res_test,Type=vcc1$Type,treat=vcc1$`treatment:ch1`)

# Calculate classification metrics
confusion_matrix <- table(test$Type,test$Ensemble.res,dnn=c("Actual", "Predicted"));confusion_matrix
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
recall <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
f1_score <- 2 * precision * recall / (precision + recall)
G_mean <- sqrt(confusion_matrix[2, 2] / sum(confusion_matrix[2, ]) * confusion_matrix[1, 1] / sum(confusion_matrix[1, ]))

predicted_probs <- test$Ensemble.ratio;actual <- test$Type
pr <- pr.curve(predicted_probs[actual == 'x1'], predicted_probs[actual == 'x0'], curve = TRUE)
Average_Precision <- pr$auc.integral

roc <- roc(test$Type,test$Ensemble.ratio)
roc_auc <- auc(roc)

cat(" Accuracy: ", accuracy, "\n",
    "Precision: ", precision, "\n",
    "Recall: ", recall, "\n",
    "F1 Score: ", f1_score, "\n",
    'G-mean:',G_mean,'\n',
    'AUC:',roc_auc,'\n',
    'Average Precision:',Average_Precision)

# Save Validation cohort #1 (external cohort) perfromace
metric <- data.frame(Metric=c("Accuracy", "Precision", "Recall", "F1 Score",'G-mean','Average Precision','AUC'),
                     Value = round(c(accuracy, precision, recall, f1_score, G_mean, Average_Precision,roc_auc),3))
metric$Metric <- factor(metric$Metric,levels = c("Accuracy", "Precision", "Recall", "F1 Score",'G-mean','Average Precision','AUC'))
#saveRDS(metric,'data/metric_val1_n3.rds')


## 3.2 CSPM performance in validation cohort #2 (in-house cohort)####
exp <- exp_self #exoressuib profile
dim(exp);dim(anno_self) #ground-truth labels

vexp1 <-as.data.frame(t((exp[intergene,]))) 
vcc1 <- anno_self

# Predict senescence using CSPM
res_test <- ensemble.predict(fit_ensemble,data.features=vexp1,outcome = "x1")
res_test$Ensemble.res <- ifelse(res_test$Ensemble.ratio>0.5,"x1","x0")
sum(res_test$Ensemble.res==as.character(vcc1$Type))/nrow(vcc1)
test <- mutate(res_test,Type=vcc1$Type,treat=vcc1$`treatment:ch1`)

# classifier performance
confusion_matrix <- table(test$Type,test$Ensemble.res,dnn=c("Actual", "Predicted"));confusion_matrix
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
recall <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
f1_score <- 2 * precision * recall / (precision + recall)
G_mean <- sqrt(confusion_matrix[2, 2] / sum(confusion_matrix[2, ]) * confusion_matrix[1, 1] / sum(confusion_matrix[1, ]))

predicted_probs <- test$Ensemble.ratio
actual <- test$Type
pr <- pr.curve(predicted_probs[actual == 'x1'], predicted_probs[actual == 'x0'], curve = TRUE)
Average_Precision <- pr$auc.integral

roc <- roc(test$Type,test$Ensemble.ratio)
roc_auc <- auc(roc);roc_auc

cat(" Accuracy: ", accuracy, "\n",
    "Precision: ", precision, "\n",
    "Recall: ", recall, "\n",
    "F1 Score: ", f1_score, "\n",
    'G-mean:',G_mean,'\n',
    'AUC:',roc_auc,'\n',
    'Average Precision:',Average_Precision)


# save classifier metrics
metric <- data.frame(Metric=c("Accuracy", "Precision", "Recall", "F1 Score",'G-mean','Average Precision','AUC'),
                     Value = round(c(accuracy, precision, recall, f1_score, G_mean, Average_Precision,roc_auc),3))
metric$Metric <- factor(metric$Metric,levels = c("Accuracy", "Precision", "Recall", "F1 Score",'G-mean','Average Precision','AUC'))
#saveRDS(metric,'data/metric_val2_n3.rds')


## 3.3 CSPM performance across single-cell datasets####

# Load CAF.seurat of each single-cell dataset
CAF_list <- readRDS("data/CAF_list.rds")

cmp <- data.frame(row.names = names(CAF_list))
d <- list()
for(j in names(CAF_list)){
  CAF.seurat <- CAF_list[[j]]
  CAF.seurat <- CAF.seurat %>% ScaleData(features = rownames(CAF.seurat))
  
  # Identify pseudo ground-truth labels
  exp <- as.data.frame(CAF.seurat@assays[["RNA"]]@data)
  g <- gene.info(rownames(exp))
  g <- g[g$gene_biotype=='protein_coding',]
  names <- unique(g$symbol)
  exp <- exp[rownames(exp)%in%names,]
  exp <- as.matrix(exp)
  
  ss <- gsva(exp, 
             geneset, 
             method = 'gsva', 
             kcdf = 'Gaussian',
             abs.ranking = F,
             mx.diff = T)
  ss <- t(ss)
  
  Group <- data.frame(row.names = rownames(ss))
  for (i in 1:ncol(ss)) {
    set.seed(123)
    fit <- Mclust(ss[,i],G = 2)
    Group[,colnames(ss)[i]] <- ifelse(fit$classification==2,'High','Low')
    pd <- data.frame(Cluster=Group[,colnames(ss)[i]],
                     Score=ss[,i])
    print(ggplot(pd, aes(x = Score, fill = Cluster)) +
            geom_histogram(aes(y = after_stat(density)), bins = 30, alpha = 0.6) +  
            geom_density(alpha = 0.7) + 
            theme_minimal() +
            labs(fill = "Cluster",title = colnames(ss)[i]))
  }
  
  Group$Type <- apply(Group,1,\(x){
    if (length(unique(x))==1) {
      return(unique(x))
    }else{
      return('Unclassfied')
    }
  })
  
  dea_cc <- Group[Group$Type!='Unclassfied',] # pseudo ground-truth labels
  
  # Predict senescence using CSPM
  exp_scale <-  as.data.frame(CAF.seurat@assays[["RNA"]]@scale.data)
  mexp <- exp_scale[intergene,rownames(Group)[Group$Type!='Unclassfied']] %>% t() %>% as.data.frame()
  
  if(j=='CNP0004138'){
    # Using only test set in discovery cohort
    set.seed(123)
    vexp <- mexp[!rownames(mexp)%in%rownames(texp),]
    vcc <- Group[rownames(vexp),]
    vcc$Type <- ifelse(vcc$Type=='High',"x1","x0")
    
    res_test <- ensemble.predict(fit_ensemble,data.features=vexp,outcome = "x1")
    res_test$Ensemble.res <- ifelse(res_test$Ensemble.ratio>0.5,"x1","x0")
    
    # classifier metrics
    confusion_matrix <- table(vcc$Type,res_test$Ensemble.res,dnn=c("Actual", "Predicted"))
    
    accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
    precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
    recall <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
    f1_score <- 2 * precision * recall / (precision + recall)
    G_mean <- sqrt(confusion_matrix[2, 2] / sum(confusion_matrix[2, ]) * confusion_matrix[1, 1] / sum(confusion_matrix[1, ]))
    
    predicted_probs <- res_test$Ensemble.ratio
    actual <- vcc$Type
    pr <- pr.curve(predicted_probs[actual == 'x1'], predicted_probs[actual == 'x0'], curve = TRUE)
    Average_Precision <- pr$auc.integral
    
    roc <- roc(vcc$Type,res_test$Ensemble.ratio)
    roc_auc <- auc(roc)
    
    d[[j]]$vcc <- vcc
    d[[j]]$res_test <- res_test
    
    cmp[j,'Accuracy'] <- accuracy
    cmp[j,'Precision'] <- precision
    cmp[j,'Recall'] <- recall
    cmp[j,'F1 Score'] <- f1_score
    cmp[j,'G-mean'] <- G_mean
    cmp[j,'Average Precision'] <- Average_Precision
    cmp[j,'AUC'] <- roc_auc
  }else{
    # Using all the pseudo ground-truth labels in other datasets
    vexp <- mexp
    vcc <- Group[rownames(vexp),]
    vcc$Type <- ifelse(vcc$Type=='High',"x1","x0")
    
    res_test <- ensemble.predict(fit_ensemble,data.features=vexp,outcome = "x1")
    res_test$Ensemble.res <- ifelse(res_test$Ensemble.ratio>0.5,"x1","x0")
    
    # classifier metrics
    confusion_matrix <- table(vcc$Type,res_test$Ensemble.res,dnn=c("Actual", "Predicted"))
    
    accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
    precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
    recall <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
    f1_score <- 2 * precision * recall / (precision + recall)
    G_mean <- sqrt(confusion_matrix[2, 2] / sum(confusion_matrix[2, ]) * confusion_matrix[1, 1] / sum(confusion_matrix[1, ]))
    
    predicted_probs <- res_test$Ensemble.ratio
    actual <- vcc$Type
    pr <- pr.curve(predicted_probs[actual == 'x1'], predicted_probs[actual == 'x0'], curve = TRUE)
    Average_Precision <- pr$auc.integral
    
    roc <- roc(vcc$Type,res_test$Ensemble.ratio)
    roc_auc <- auc(roc)
    
    d[[j]]$vcc <- vcc
    d[[j]]$res_test <- res_test
    
    cmp[j,'Accuracy'] <- accuracy
    cmp[j,'Precision'] <- precision
    cmp[j,'Recall'] <- recall
    cmp[j,'F1 Score'] <- f1_score
    cmp[j,'G-mean'] <- G_mean
    cmp[j,'Average Precision'] <- Average_Precision
    cmp[j,'AUC'] <- roc_auc
  }
}


# Merge labels across datasets
vcc <- do.call(rbind, lapply(d, function(x) x$vcc))
res_test <- do.call(rbind, lapply(d, function(x) x$res_test))

confusion_matrix <- table(vcc$Type,res_test$Ensemble.res,dnn=c("Actual", "Predicted"))
confusion_matrix

accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
recall <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
f1_score <- 2 * precision * recall / (precision + recall)
G_mean <- sqrt(confusion_matrix[2, 2] / sum(confusion_matrix[2, ]) * confusion_matrix[1, 1] / sum(confusion_matrix[1, ]))

predicted_probs <- res_test$Ensemble.ratio
actual <- vcc$Type
pr <- pr.curve(predicted_probs[actual == 'x1'], predicted_probs[actual == 'x0'], curve = TRUE)
Average_Precision <- pr$auc.integral

roc <- roc(vcc$Type,res_test$Ensemble.ratio)
roc_auc <- auc(roc)

cmp['Merge','Accuracy'] <- accuracy
cmp['Merge','Precision'] <- precision
cmp['Merge','Recall'] <- recall
cmp['Merge','F1 Score'] <- f1_score
cmp['Merge','G-mean'] <- G_mean
cmp['Merge','Average Precision'] <- Average_Precision
cmp['Merge','AUC'] <- roc_auc

# save CSPM performance (remove duplicated genes n>=3)
# saveRDS(cmp,'data/cmp_n3.rds')


# ------------------------------C: Performance comparison-----------------------

# 1.Comparison in pseudo ground-truth labels--------------------------
cmp1 <- readRDS("data/cmp.rds") %>% mutate(version='retain') %>% rownames_to_column('cohort') # retain duplicated genes
cmp2 <- readRDS("data/cmp_n2.rds") %>% mutate(version='remove (n>=2)')%>% rownames_to_column('cohort') # remove duplicated genes n>=2
cmp3 <- readRDS("data/cmp_n3.rds") %>% mutate(version='remove (n>=3)')%>% rownames_to_column('cohort') # remove duplicated genes n>=3

cmp <- rbind(cmp1,cmp2,cmp3) %>% pivot_longer(cols = 2:8,  names_to = "Metric",values_to = "Value") %>% filter(cohort!='Merge')
cmp$version <- factor(cmp$version,levels = c('retain','remove (n>=2)','remove (n>=3)'))

# classifier metrics (Ref Figure 2A-C)
for(i in c('Accuracy','G-mean','AUC')){
  data <- cmp %>% filter(Metric %in% i)
ggplot(data, aes(x = cohort, y = Value, color = version, group = version)) +
    geom_line(size = 1) +  
    geom_point(size = 2) + 
    labs(title = "Classifier Performance in Pseudo Ground-Truth Labels Across Cohorts",x = "",y = i, color = 'Duplicated genes') +
    scale_color_manual(values = c('retain'='#ef7166',
                                  'remove (n>=2)'='#336385',
                                  'remove (n>=3)'='#efa63a')) +
    theme_bw(base_rect_size = 1) +
    theme(
      panel.grid.major = element_line(color = "grey80", size = 0.5),
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.position = "right",  
      legend.text = element_text(size = 12, colour = 'black'),
      legend.title = element_text(face = "bold", size = 12),
      axis.text.y= element_text(size = 12, colour = 'black'),
      axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1,size = 14,color = 'black'),  
      axis.title = element_text(face = "bold", colour = "black", size = 15),
      plot.title = element_text(face = "bold", colour = "black", size = 15, hjust = 0.5))
}

# classifier metrics (Ref Figure 2D)
cmp <- rbind(cmp1,cmp2,cmp3) %>% filter(cohort=='Merge')
data <- cmp %>% pivot_longer(cols = 2:8,names_to = 'Metric',values_to = 'Value')
data$version <- factor(data$version,levels = c('retain','remove (n>=2)','remove (n>=3)'))

ggplot(data, aes(x = Metric, y = Value, color = version, group = version)) +
  geom_line(size = 1) +  
  geom_point(size = 2) + 
  labs(title = "Classifier Performance in All Pseudo Ground-Truth Labels",x = "",y = 'Value', color = 'Duplicate genes') +
  scale_color_manual(values = c('retain'='#ef7166',
                                'remove (n>=2)'='#336385',
                                'remove (n>=3)'='#efa63a')) +
  theme_bw(base_rect_size = 1) +
  theme(
    panel.grid.major = element_line(color = "grey80", size = 0.5),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.position = "right",  
    legend.text = element_text(size = 12, colour = 'black'),
    legend.title = element_text(face = "bold", size = 12),
    axis.text.y= element_text(size = 12, colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1,size = 14,color = 'black'),  # x轴标签倾斜45度
    axis.title = element_text(face = "bold", colour = "black", size = 15),
    plot.title = element_text(face = "bold", colour = "black", size = 15, hjust = 0.5))


# 2. Comparison in Ground-truth labels---------------------------------

metric_val1 <- readRDS("data/metric_val1.rds") %>% mutate(version='retain')
metric_val2 <- readRDS("data/metric_val2.rds")%>% mutate(version='retain')

metric_val1_n2 <- readRDS("data/metric_val1_n2.rds") %>% mutate(version='remove (n>=2)')
metric_val2_n2 <- readRDS("data/metric_val2_n2.rds") %>% mutate(version='remove (n>=2)')

metric_val1_n3 <- readRDS("data/metric_val1_n3.rds") %>% mutate(version='remove (n>=3)')
metric_val2_n3 <- readRDS("data/metric_val2_n3.rds") %>% mutate(version='remove (n>=3)')


# Validation cohort #1
data <- rbind(metric_val1,metric_val1_n2,metric_val1_n3)
data$version <- factor(data$version,levels = c('retain','remove (n>=2)','remove (n>=3)'))

ggplot(data, aes(x = Metric, y = Value, color = version, group = version)) +
  geom_line(size = 1) +  
  geom_point(size = 2) + 
  labs(title = "Classifier Performance in Validation Cohort #1 (external cohort)",x = "",y = 'Value', color = 'Duplicated genes') +
  scale_color_manual(values = c('retain'='#ef7166',
                                'remove (n>=2)'='#336385',
                                'remove (n>=3)'='#efa63a')) +
  theme_bw(base_rect_size = 1) +
  theme(
    panel.grid.major = element_line(color = "grey80", size = 0.5),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.position = "right",  
    legend.text = element_text(size = 12, colour = 'black'),
    legend.title = element_text(face = "bold", size = 12),
    axis.text.y= element_text(size = 12, colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1,size = 14,color = 'black'),  # x轴标签倾斜45度
    axis.title = element_text(face = "bold", colour = "black", size = 15),
    plot.title = element_text(face = "bold", colour = "black", size = 15, hjust = 0.5))
ggsave('Figuredir/RefFigure1E.pdf')


# Validation cohort #2
data <- rbind(metric_val2,metric_val2_n2,metric_val2_n3)
data$version <- factor(data$version,levels = c('retain','remove (n>=2)','remove (n>=3)'))
ggplot(data, aes(x = Metric, y = Value, color = version, group = version)) +
  geom_line(size = 1) +  
  geom_point(size = 2) + 
  labs(title = "Classifier Performance in Validation Cohort #2 (in-house)",x = "",y = 'Value', color = 'Duplicated genes') +
  scale_color_manual(values = c('retain'='#ef7166',
                                'remove (n>=2)'='#336385',
                                'remove (n>=3)'='#efa63a')) +
  theme_bw(base_rect_size = 1) +
  theme(
    panel.grid.major = element_line(color = "grey80", size = 0.5),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.position = "right",  
    legend.text = element_text(size = 12, colour = 'black'),
    legend.title = element_text(face = "bold", size = 12),
    axis.text.y= element_text(size = 12, colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1,size = 14,color = 'black'),  # x轴标签倾斜45度
    axis.title = element_text(face = "bold", colour = "black", size = 15),
    plot.title = element_text(face = "bold", colour = "black", size = 15, hjust = 0.5))
ggsave('Figuredir/RefFigure1F.pdf')


