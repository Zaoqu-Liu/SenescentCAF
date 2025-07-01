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

# This pipeline compares the performance of the Cellular Senescence Prediction Model (CSPM) with GSVA and AUCell in two ground-truth validation cohorts.
# Example data is available at https://doi.org/10.7303/syn68545972.


# 1. Data preprocessiiong----------------------------------------

# CSPM performance 
metric_val1 <- readRDS("data/metric_val1.rds") %>% 
  pivot_wider(names_from = Metric,values_from = Value) %>% 
  mutate(category='CSPM')
metric_val2 <- readRDS("data/metric_val2.rds") %>% 
  pivot_wider(names_from = Metric,values_from = Value) %>% 
  mutate(category='CSPM')

readRDS('data/geneset.rds') # senescence gene sets

# Load validation datasets (all available genes)
load('data/validation_data.rda')
dim(exp1);dim(anno1) #validation cohort #1
dim(exp2);dim(anno2) #validation cohort #2


# 2. CSPM VS GSVA------------------------------------------------

# validation cohort #1
ss <- gsva(exp1, 
           geneset, 
           method = 'gsva', 
           kcdf = 'Gaussian',
           abs.ranking = F,
           mx.diff = T)
ss <- t(ss)

pl <- list()
for(probs in c(0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8)){ 
  df <- data.frame()
  for(i in names(geneset)){
    d <- ss[, i, drop = FALSE] %>% as.data.frame()
    threshold <- quantile(d[[i]], probs = probs, na.rm = TRUE) #select threshold
    d$Type <- ifelse(d[[i]] >= threshold, "x1", "x0")
    
    # GSVA classifier performance
    confusion_matrix <- table(anno1$Type,d$Type,dnn=c("Actual", "Predicted"))
    
    accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
    precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
    recall <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
    f1_score <- 2 * precision * recall / (precision + recall)
    G_mean <- sqrt(confusion_matrix[2, 2] / sum(confusion_matrix[2, ]) * confusion_matrix[1, 1] / sum(confusion_matrix[1, ]))
    
    predicted_probs <- d[[i]]
    actual <- anno1$Type
    pr <- pr.curve(predicted_probs[actual == 'x1'], predicted_probs[actual == 'x0'], curve = TRUE)
    Average_Precision <- pr$auc.integral
    
    roc <- roc(anno1$Type,d[[i]])
    roc_auc <- auc(roc)
    
    dd <- data.frame(d=c("Accuracy", "Precision", "Recall", "F1 Score",'G-mean','Average Precision','AUC'),
                     Value = c(accuracy, precision, recall, f1_score, G_mean, Average_Precision,roc_auc)) %>% 
      pivot_wider(names_from = d,values_from = Value) %>% 
      mutate(category=i)
    df <- rbind(df,dd)
  }
  data <- rbind(metric_val1,df) %>% pivot_longer(cols = 1:7,names_to = 'Metric',values_to = 'Value') %>%
    mutate(category = factor(category, levels = rev(c('CSPM',names(geneset)))))
  
  top_percent <- (1 - probs) * 100
  threshold_label <- sprintf("Threshold: Top%.0f%%", top_percent)
  
  pl[[as.character(probs)]] <-   ggplot(data, aes(x = Metric, y = Value, color = category, group = category)) +
    geom_line(size = 1) +  
    geom_point(size = 2) + 
    labs(title = threshold_label,x = "",y = "Metric Value", color = "Method") +
    scale_color_manual(values = rev(cols[1:9]),guide = guide_legend(reverse = T)) +
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
      plot.title = element_text(face = "bold", colour = "black", size = 15, hjust = 0.5)
    )
}
wrap_plots(pl,nrow = 3)+plot_layout(guides = "collect")+
  plot_annotation(title = 'Classifier Performance vs GSVA in Validation Cohort #1')&
  theme(plot.title = element_text(size = 20,colour = 'black',face = 'bold',hjust = 0.3))
ggsave('Figuredir/RefFigure4.pdf')


# validation cohort #2
ss <- gsva(exp2, 
           geneset, 
           method = 'gsva', 
           kcdf = 'Gaussian',
           abs.ranking = F,
           mx.diff = T)
ss <- t(ss)

pl <- list()
for(probs in c(0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8)){
  df <- data.frame()
  for(i in names(geneset)){
    d <- ss[, i, drop = FALSE] %>% as.data.frame()
    threshold <- quantile(d[[i]], probs = probs, na.rm = TRUE) #select threshold
    d$Type <- ifelse(d[[i]] >= threshold, "x1", "x0")
    
    # GSVA classifier formance
    confusion_matrix <- table(anno2$Type,d$Type,dnn=c("Actual", "Predicted"))
    
    accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
    precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
    recall <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
    f1_score <- 2 * precision * recall / (precision + recall)
    G_mean <- sqrt(confusion_matrix[2, 2] / sum(confusion_matrix[2, ]) * confusion_matrix[1, 1] / sum(confusion_matrix[1, ]))
    
    predicted_probs <- d[[i]]
    actual <- anno2$Type
    pr <- pr.curve(predicted_probs[actual == 'x1'], predicted_probs[actual == 'x0'], curve = TRUE)
    Average_Precision <- pr$auc.integral
    
    roc <- roc(anno2$Type,d[[i]])
    roc_auc <- auc(roc)
    
    dd <- data.frame(d=c("Accuracy", "Precision", "Recall", "F1 Score",'G-mean','Average Precision','AUC'),
                     Value = c(accuracy, precision, recall, f1_score, G_mean, Average_Precision,roc_auc)) %>% 
      pivot_wider(names_from = d,values_from = Value) %>% 
      mutate(category=i)
    df <- rbind(df,dd)
  }
  data <- rbind(metric_val2,df) %>% pivot_longer(cols = 1:7,names_to = 'Metric',values_to = 'Value') %>%
    mutate(category = factor(category, levels = rev(c('CSPM',names(geneset)))))
  
  top_percent <- (1 - probs) * 100
  threshold_label <- sprintf("Threshold: Top%.0f%%", top_percent)
  
  pl[[as.character(probs)]] <-   ggplot(data, aes(x = Metric, y = Value, color = category, group = category)) +
    geom_line(size = 1) +  
    geom_point(size = 2) + 
    labs(title = threshold_label,x = "",y = "Metric Value", color = "Method") +
    scale_color_manual(values = rev(cols[1:9]),guide = guide_legend(reverse = T)) +
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
      plot.title = element_text(face = "bold", colour = "black", size = 15, hjust = 0.5)
    )
}
wrap_plots(pl,nrow = 3)+plot_layout(guides = "collect")+
  plot_annotation(title = 'Classifier Performance vs GSVA in Validation Cohort #2 (in-house cohort)')&
  theme(plot.title = element_text(size = 20,colour = 'black',face = 'bold',hjust = 0.3))
ggsave('Figuredir/RefFigure5.pdf')


# 3. CSPM VS AUCell--------------------------------------------------

# validation cohort #1
cells_rankings <- AUCell_buildRankings(exp1, nCores=10)  
cells_AUC <- AUCell_calcAUC(geneset, cells_rankings,nCores = 10, 
                            aucMaxRank = nrow(cells_rankings)*0.1)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, assign=F, nCores=10)
AUC <- getAUC(cells_AUC) %>% t() %>% as.data.frame()

df <- data.frame()
for(i in names(geneset)){
  threshold <- cells_assignment[[i]]$aucThr$thresholds['Global_k1','threshold'] #select threshold
  d <- AUC[,i,drop=F] %>% mutate(Type=ifelse(.data[[i]] >= threshold, "x1", "x0"))
  
  # AUCell classifier performance
  confusion_matrix <- table(anno1$Type,d$Type,dnn=c("Actual", "Predicted"))
  
  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
  precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
  recall <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
  f1_score <- 2 * precision * recall / (precision + recall)
  G_mean <- sqrt(confusion_matrix[2, 2] / sum(confusion_matrix[2, ]) * confusion_matrix[1, 1] / sum(confusion_matrix[1, ]))
  
  predicted_probs <- d[[i]]
  actual <- anno1$Type
  pr <- pr.curve(predicted_probs[actual == 'x1'], predicted_probs[actual == 'x0'], curve = TRUE)
  Average_Precision <- pr$auc.integral
  
  roc <- roc(anno1$Type,d[[i]])
  roc_auc <- auc(roc)
  
  dd <- data.frame(d=c("Accuracy", "Precision", "Recall", "F1 Score",'G-mean','Average Precision','AUC'),
                   Value = c(accuracy, precision, recall, f1_score, G_mean, Average_Precision,roc_auc)) %>% 
    pivot_wider(names_from = d,values_from = Value) %>% 
    mutate(category=i)
  df <- rbind(df,dd)
}

data <- rbind(metric_val1,df) %>% pivot_longer(cols = 1:7,names_to = 'Metric',values_to = 'Value') %>%
  mutate(category = factor(category, levels = rev(c('CSPM',names(geneset)))))

p1 <- ggplot(data, aes(x = Metric, y = Value, color = category, group = category)) +
  geom_line(size = 1) +  
  geom_point(size = 2) + 
  labs(title = 'Classifier Performance vs AUCell in Validation Cohort #1',x = "",y = "Metric Value", color = "Method") +
  scale_color_manual(values = rev(cols[1:9]),guide = guide_legend(reverse = T)) +
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
    plot.title = element_text(face = "bold", colour = "black", size = 18, hjust = 0.5)
  )

# validation cohort #2
cells_rankings <- AUCell_buildRankings(exp2, nCores=10)  
cells_AUC <- AUCell_calcAUC(geneset, cells_rankings,nCores = 10, 
                            aucMaxRank = nrow(cells_rankings)*0.1)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, assign=F, nCores=10)
AUC <- getAUC(cells_AUC) %>% t() %>% as.data.frame()

df <- data.frame()
for(i in names(geneset)){
  threshold <- cells_assignment[[i]]$aucThr$thresholds['Global_k1','threshold'] #sekect threshold
  d <- AUC[,i,drop=F] %>% mutate(Type=ifelse(.data[[i]] >= threshold, "x1", "x0"))
  
  # AUCell performance
  confusion_matrix <- table(anno2$Type,d$Type,dnn=c("Actual", "Predicted"))
  
  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
  precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
  recall <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
  f1_score <- 2 * precision * recall / (precision + recall)
  G_mean <- sqrt(confusion_matrix[2, 2] / sum(confusion_matrix[2, ]) * confusion_matrix[1, 1] / sum(confusion_matrix[1, ]))
  
  predicted_probs <- d[[i]]
  actual <- anno2$Type
  pr <- pr.curve(predicted_probs[actual == 'x1'], predicted_probs[actual == 'x0'], curve = TRUE)
  Average_Precision <- pr$auc.integral
  
  roc <- roc(anno2$Type,d[[i]])
  roc_auc <- auc(roc)
  
  dd <- data.frame(d=c("Accuracy", "Precision", "Recall", "F1 Score",'G-mean','Average Precision','AUC'),
                   Value = c(accuracy, precision, recall, f1_score, G_mean, Average_Precision,roc_auc)) %>% 
    pivot_wider(names_from = d,values_from = Value) %>% 
    mutate(category=i)
  df <- rbind(df,dd)
}

data <- rbind(metric_val2,df) %>% pivot_longer(cols = 1:7,names_to = 'Metric',values_to = 'Value') %>%
  mutate(category = factor(category, levels = rev(c('CSPM',names(geneset)))))

p2 <- ggplot(data, aes(x = Metric, y = Value, color = category, group = category)) +
  geom_line(size = 1) +  
  geom_point(size = 2) + 
  labs(title = 'Classifier Performance vs AUCell in Validation Cohort #2 (in-house cohort)',x = "",y = "Metric Value", color = "Method") +
  scale_color_manual(values = rev(cols[1:9]),guide = guide_legend(reverse = T)) +
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
    plot.title = element_text(face = "bold", colour = "black", size = 18, hjust = 0.5)
  )


p1+p2+plot_layout(guides = 'collect')
ggsave('Figuredir/RefFigure6.pdf')