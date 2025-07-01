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


# This pipeline shows the development and validation processes of the Cellular Senescence Prediction Model (CSPM).
# Example data is available at https://doi.org/10.7303/syn68545972.


# 1. Data Preprocessing-----------------------------------------------------

# Load single-cell data and senescence gene sets
geneset <- readRDS('data/geneset.rds')
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
load('data/CSPM_validation.rda')
dim(exp_public);dim(anno_all) #validation cohort #1
dim(exp_self);dim(anno_self) #validation cohort #2

# 2. CSPM construction----------------------------------------------------------

## 2.1 Model training####

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
#saveRDS(fit_ensemble,'data/fit_ensemble.rds')


# Training set accuracy
res_train <- ensemble.predict(fit_ensemble,data.features,outcome = "x1")
res_train$Ensemble.res <- ifelse(res_train$Ensemble.ratio>0.5,"x1","x0")
sum(res_train$Ensemble.res==as.character(tcc$Type))/nrow(tcc)

# Test set accuracy
res_test <- ensemble.predict(fit_ensemble,data.features=vexp,outcome = "x1")
res_test$Ensemble.res <- ifelse(res_test$Ensemble.ratio>0.5,"x1","x0")
sum(res_test$Ensemble.res==as.character(vcc$Type))/nrow(vcc)

## 2.2 Identify senescent CAFs in discovery cohort####

# Predict senescence state of unclassified cells
iexp <- exp_scale[intergene,rownames(Group)[Group$Type=='Unclassfied']] %>% t() %>% as.data.frame()
res_ident <- ensemble.predict(fit_ensemble,data.features=iexp,outcome = "x1")
res_ident$Ensemble.res <- ifelse(res_ident$Ensemble.ratio>0.5,"x1","x0")

data <- CAF.seurat@meta.data %>% mutate(Type_ident=Group$Type)
n=match(rownames(Group)[Group$Type=='Unclassfied'],rownames(data))
data$Type_ident[n] <- res_ident$Ensemble.res
data$Type_ident[data$Type_ident=="x1"] <- "High"
data$Type_ident[data$Type_ident=="x0"] <- "Low"
data$Type_ident <- ifelse(data$Type_ident=="High","sCAF","nsCAF")
data$Type_ident <- factor(data$Type_ident)
CAF.seurat@meta.data <- data

table(CAF.seurat$Type_ident)
#saveRDS(CAF.seurat,'data/CAF.seurat_ident.rds')

# Define sCAF signature
Idents(CAF.seurat) <- CAF.seurat$Type_ident
markers <- FindAllMarkers(CAF.seurat, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25) %>% filter(p_val_adj<0.05)

split_markers <- split(markers,markers$cluster)
split_markers <- lapply(split_markers, \(x) 
                        filter(x,avg_log2FC>0.25,
                               pct.1>0.4,pct.2<0.5,
                               pct.1-pct.2>0.2))

sCAF_signature <- split_markers$sCAF$gene

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

# ROC curve
ggroc(roc, legacy.axes = TRUE,alpha = 1, colour = '#327eb9',  linetype = 1, size = 1)+ 
  geom_segment(aes(x = -0.005, y = -0.005, xend = 1.005, yend = 1.005), linetype = "dashed", color = "#868686FF") +
  geom_text(aes(x = 0.7, y = 0.34, label = paste("AUC =", round(roc_auc, 2))), size = 6, color = "black") +
  labs(x = "1 - Specificity", y = "Sensitivity",title = "Validation Cohort #1") +
  theme_bw(base_rect_size = 1.8)+
  theme(axis.title = element_text(size = 18,face = 'bold'),
        axis.text.x = element_text(colour = 'black', size = 12),
        axis.text.y = element_text(colour = 'black', size = 12),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        panel.grid = element_blank(),
        panel.grid.major =element_line(color = "#cacfd2", linetype = "dashed"),
  )+
  scale_x_continuous(expand = c(0.02, 0),
                     breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(expand = c(0.02,0),
                     breaks = seq(0, 1, by = 0.2))
ggsave('Figure_dir/Figure1D.pdf')

# Classifier metrics
data <- data.frame(d=c("Accuracy", "Precision", "Recall", "F1 Score",'G-mean','Average Precision'),
                   Value = round(c(accuracy, precision, recall, f1_score, G_mean, Average_Precision),4))
data$d <- factor(data$d,levels = c("Accuracy", "Precision", "Recall", "F1 Score",'G-mean','Average Precision'))

ggplot(data, aes(x = Value, y = d)) +  
  geom_segment(aes(x = 0, xend = Value, y = d, yend = d), color='black',linewidth = 0.85)+
  geom_point(aes(color=d,fill = d),size = 10, shape = 21) +  
  geom_text(aes(label = Value), hjust = 0.5, size = 5)+
  scale_fill_manual(values = c('#ed6355','#41a98e','#00b4d8','#efa63a','#ab3282','#f4ee72'))+
  scale_color_manual(values = c('#ed6355','#41a98e','#00b4d8','#efa63a','#ab3282','#f4ee72'))+
  theme_classic()+
  labs(x="Value",y="",fill="",title = 'Classifier Metrics')+
  theme(axis.text.y = element_text(size = 16, colour = 'black',vjust = 0.5),
        axis.text.x = element_text(size = 12,colour = 'black',vjust = 0.5),
        axis.title.x = element_text(size = 18,colour = 'black',face = 'bold'),
        plot.title = element_text(hjust = 0.5,vjust = 1, size = 18, face = "bold"),
        axis.line = element_line(linewidth = 0.8),
        legend.position = 'none',
        plot.margin = margin(0,0.5,0,0,'cm'))
ggsave('Figure_dir/Figure1E.pdf')

# Save Validation cohort #1 (external cohort) perfromace
metric <- data.frame(Metric=c("Accuracy", "Precision", "Recall", "F1 Score",'G-mean','Average Precision','AUC'),
                     Value = round(c(accuracy, precision, recall, f1_score, G_mean, Average_Precision,roc_auc),3))
metric$Metric <- factor(metric$Metric,levels = c("Accuracy", "Precision", "Recall", "F1 Score",'G-mean','Average Precision','AUC'))
#saveRDS(metric,'data/metric_val1.rds')


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

# ROC curve
ggroc(roc, legacy.axes = TRUE,
      alpha = 1, colour = '#327eb9',  linetype = 1, size = 1)+ 
  geom_segment(aes(x = -0.005, y = -0.005, xend = 1.005, yend = 1.005), linetype = "dashed", color = "#868686FF") +
  geom_text(aes(x = 0.7, y = 0.34, label = paste("AUC =", round(roc_auc, 2))), size = 6, color = "black") +
  labs(x = "1 - Specificity", y = "Sensitivity",title = "Validation Cohort #2 (in-house)") +
  theme_bw(base_rect_size = 1.8)+
  theme(axis.title = element_text(size = 18,face = 'bold'),
        axis.text.x = element_text(colour = 'black', size = 12),
        axis.text.y = element_text(colour = 'black', size = 12),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        panel.grid = element_blank(),
        panel.grid.major =element_line(color = "#cacfd2", linetype = "dashed"),
  )+
  scale_x_continuous(expand = c(0.02, 0),
                     breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(expand = c(0.02,0),
                     breaks = seq(0, 1, by = 0.2))
ggsave('Figure_dir/FigureS1C.pdf')

# save classifier metrics
metric <- data.frame(Metric=c("Accuracy", "Precision", "Recall", "F1 Score",'G-mean','Average Precision','AUC'),
                     Value = round(c(accuracy, precision, recall, f1_score, G_mean, Average_Precision,roc_auc),3))
metric$Metric <- factor(metric$Metric,levels = c("Accuracy", "Precision", "Recall", "F1 Score",'G-mean','Average Precision','AUC'))
#saveRDS(metric,'data/metric_val2.rds')


## 3.3 CSPM performance across single-cell datasets####

# Load CAF.seurat of each single-cell dataset
CAF_list <- readRDS("data/CAF_list.rds")

### 3.3.1 Validation in pseudo ground-truth labels across datasets####
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

# Classifier accuracy
data <- data.frame(d=c(rownames(cmp)[-9],'CNP0004138-Test Set'),Value = round(cmp$Accuracy,3))
data$d <- factor(data$d,levels = data$d)
ggplot(data, aes(x = Value, y = d)) +  
  geom_segment(aes(x = 0, xend = Value, y = d, yend = d), color='black',linewidth = 0.7)+
  geom_point(aes(color=d,fill = d),size = 10, shape = 21) +  
  geom_text(aes(label = Value), hjust = 0.55, size = 5)+
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols)+
  theme_classic()+
  labs(x="Accuracy",y="",fill="",title = 'Classifier Performance\nin Pseudo Ground Truth Labels')+
  theme(axis.text.y = element_text(size = 14, colour = 'black',vjust = 0.5),
        axis.text.x = element_text(size = 12,colour = 'black',vjust = 0.5),
        axis.title.x = element_text(size = 18,colour = 'black',face = 'bold'),
        plot.title = element_text(hjust = 0.5,vjust = 1, size = 18, face = "bold"),
        axis.line = element_line(linewidth = 0.7),
        legend.position = 'none',
        plot.margin = margin(0,0.5,0,0,'cm'))
ggsave('Figure_dir/FigureS2A.pdf')

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

# save CSPM performance 
# saveRDS(cmp,'data/cmp.rds')

### 3.3.2 Identify senescent CAFs across datasets####
CAF_list_ident <- lapply(names(CAF_list),function(x){
  CAF.seurat <- CAF_list[[x]]
  if(x=='CNP0004138'){
    # Discovery cohort
    CAF.seurat <- readRDS("data/CAF.seurat_ident.rds")
    CAF.seurat$Major.type2 <- CAF.seurat$Type_ident
  }else{
    # Other independent datasets
    CAF.seurat <- CAF.seurat %>% ScaleData(features = rownames(CAF.seurat))
    
    exp_scale <-  as.data.frame(CAF.seurat@assays[["RNA"]]@scale.data)
    iexp <- exp_scale[intergene,] %>% t() %>% as.data.frame()
    
    res_ident <- ensemble.predict(fit_ensemble,data.features=iexp,outcome = "x1")
    res_ident$Ensemble.res <- ifelse(res_ident$Ensemble.ratio>0.5,"x1","x0")
    table(res_ident$Ensemble.res)
    
    CAF.seurat$Type_ident <- res_ident$Ensemble.res
    CAF.seurat$Type_ident <- ifelse(CAF.seurat$Type_ident=="x1","sCAF","nsCAF")
    CAF.seurat$Major.type2 <- factor(CAF.seurat$Type_ident)
  }
  return(CAF.seurat)
})
names(CAF_list_ident) <- names(CAF_list)

lapply(CAF_list_ident,function(x){
  table(x$Major.type2)
})

#saveRDS(CAF_list_ident,'data/CAF_list_ident.rds')


