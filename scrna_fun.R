###分析主流程###
library(Seurat)
library(harmony)
library(parallel)
library(dplyr)
library(reshape2)
library(data.table)
library(HGNChelper)
library(readxl)
library(reshape2)
library(stringr)
library(tibble)
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
library(enrichplot)
library(tidyverse)
options(bitmapType='cairo')

##使用SCTransform + harmony 
#细胞过滤函数
filter <- function(sce){
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^M(?i:T)-")
  sce[["percent.HB"]] <- PercentageFeatureSet(sce, pattern = "H(?i:B)")   ###(?i:)表示该字符不区分大小写
  ##设定阈值标准，主要是nFeature_RNA，设定90%阈值,
  #smart-seq是不是不用过滤
  maxnF <- max(sce@meta.data$nFeature_RNA)
  nF_cutoff <- round(0.9*maxnF)
  #VlnPlot(sce, features = c("nFeature_RNA", "percent.mt", "percent.HB","percent.RP"), ncol = 4)
  sce <- subset(sce, subset = nFeature_RNA > 200 & nFeature_RNA < nF_cutoff & percent.mt < 10)
}
#聚类函数，分为两个，一个是一个GSE有多个样本需要整合分析（cluster_multi_SRR），另外是一个GSE只有一个样本，不需要整合分析（cluster_one_SRR）
cluster_multi_SRR <- function(sce_list,GSE){
  ye <- sce_list[[2]]
  if(length(sce_list)>2){
    for(k in 3:length(sce_list)){
      ye <- append(ye,sce_list[[k]])
    }
  }else{
    ye <- ye
  }
  scRNA <- merge(sce_list[[1]],y = ye,project = GSE)
  ##SCT标准化数据，代替NormalizeData,FindVariableFeatures,ScaleData三个函数
  scRNA <- SCTransform(scRNA,vars.to.regress = "percent.mt",verbose = FALSE)
  scRNA <- RunPCA(scRNA)
  scRNA <- RunHarmony(scRNA,group.by.vars = "orig.ident",assay.use = "SCT",max.iter.harmony = 10)
  # group.by.vars参数是设置按哪个分组来整合
  # max.iter.harmony设置迭代次数，默认是10。运行RunHarmony结果会提示在迭代多少次后完成了收敛。
  # RunHarmony函数中有个lambda参数，默认值是1，决定了Harmony整合的力度。lambda值调小，整合力度变大，反之。（只有这个参数影响整合力度，调整范围一般在0.5-2之间）
  
  ###确认最佳PC
  # Determine percent of variation associated with each PC
  pct <- scRNA [["pca"]]@stdev / sum( scRNA [["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  #主成分累积贡献大于90%
  #PC本身对方差贡献小于5%
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC  两个连续PCs之间差异小于0.1%
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  # last point where change of % of variation is more than 0.1%.
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  pc.num <- 1:pcs
  res <-  c(0.2,0.4,0.6,0.8,1.0) #改为5个resolution 20220616
  
  scRNA <- RunTSNE(scRNA, reduction="harmony", dims=pc.num) %>% RunUMAP(reduction="harmony", dims=pc.num) %>%
    FindNeighbors(reduction="harmony", dims=pc.num) %>% FindClusters(resolution=res)
  return(scRNA)
}
cluster_one_SRR <- function(sce_list,GSE){
  scRNA <- SCTransform(sce_list[[1]],vars.to.regress = "percent.mt",verbose = FALSE)
  scRNA <- RunPCA(scRNA)
  pct <- scRNA [["pca"]]@stdev / sum( scRNA [["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  pcs <- min(co1, co2)
  pc.num <- 1:pcs
  res <-  c(0.2,0.4,0.6,0.8,1.0)
  scRNA <- RunTSNE(scRNA, reduction="pca", dims=pc.num) %>% RunUMAP(reduction="pca", dims=pc.num) %>%
    FindNeighbors(reduction="pca", dims=pc.num) %>% FindClusters(resolution=res)
  return(scRNA)
}

#添加困惑度参数   2022.4.4
cluster_multi_SRR_smart <- function(sce_list,GSE){
  ye <- sce_list[[2]]
  if(length(sce_list)>2){
    for(k in 3:length(sce_list)){
      ye <- append(ye,sce_list[[k]])
    }
  }else{
    ye <- ye
  }
  scRNA <- merge(sce_list[[1]],y = ye,project = GSE)
  ##SCT标准化数据，代替NormalizeData,FindVariableFeatures,ScaleData三个函数
  ## 报错 cell attribute "log_umi" contains NA, NaN, or infinite value，怀疑是测序深度过低导致umi太小, 解决办法，去掉所有测序reads小于2的样本
  scRNA <- SCTransform(scRNA,vars.to.regress = "percent.mt",verbose = FALSE)
  
  scRNA <- RunPCA(scRNA)
  scRNA <- RunHarmony(scRNA,group.by.vars = "orig.ident",assay.use = "SCT",max.iter.harmony = 10)
  # group.by.vars参数是设置按哪个分组来整合
  # max.iter.harmony设置迭代次数，默认是10。运行RunHarmony结果会提示在迭代多少次后完成了收敛。
  # RunHarmony函数中有个lambda参数，默认值是1，决定了Harmony整合的力度。lambda值调小，整合力度变大，反之。（只有这个参数影响整合力度，调整范围一般在0.5-2之间）
  
  ###确认最佳PC
  # Determine percent of variation associated with each PC
  pct <- scRNA [["pca"]]@stdev / sum( scRNA [["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  #主成分累积贡献大于90%
  #PC本身对方差贡献小于5%
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC  两个连续PCs之间差异小于0.1%
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  # last point where change of % of variation is more than 0.1%.
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  pc.num <- 1:pcs
  res <-  c(0.2,0.4,0.6,0.8,1.0)
  perplexity <- floor((dim(scRNA)[2]-1)/3 - 1)
  scRNA <- RunTSNE(scRNA, reduction="harmony", dims=pc.num, perplexity = perplexity) %>% RunUMAP(reduction="harmony", dims=pc.num) %>%
    FindNeighbors(reduction="harmony", dims=pc.num) %>% FindClusters(resolution=res)
  return(scRNA)
}
cluster_one_SRR_smart <- function(sce_list,GSE){
  scRNA <- SCTransform(sce_list[[1]],vars.to.regress = "percent.mt",verbose = FALSE)
  scRNA <- RunPCA(scRNA)
  pct <- scRNA [["pca"]]@stdev / sum( scRNA [["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  pcs <- min(co1, co2)
  pc.num <- 1:pcs
  res <-  c(0.2,0.4,0.6,0.8,1.0)
  perplexity <- floor((dim(scRNA)[2]-1)/3 - 1)
  scRNA <- RunTSNE(scRNA, reduction="pca", dims=pc.num, perplexity = perplexity) %>% RunUMAP(reduction="pca", dims=pc.num) %>%
    FindNeighbors(reduction="pca", dims=pc.num) %>% FindClusters(resolution=res)
  return(scRNA)
}

#输出2D坐标点信息（包括x，y轴位置（2位小数），cluster信息，细胞定义结果，后续加样本基本信息）,虽然是tsne命名，但是也包含了UMAP数据
output_2D <- function(scRNA,rslt,studyid){
  tsne_2D_list <- list()
  for(t in unique(scRNA$tissue)){
    tissue <- t
    scRNA_subtissue <- subset(scRNA, tissue == t)
    
    org <- unique(scRNA_subtissue$organism)
    res <-  c(0.2,0.4,0.6,0.8,1.0)
    cluster_info <- scRNA_subtissue@meta.data[,paste0("SCT_snn_res.",res)]
    #输出2D坐标
    tsne_2D <- round(as.data.frame(scRNA_subtissue[['tsne']]@cell.embeddings),2) %>% cbind(round(as.data.frame(scRNA_subtissue[['umap']]@cell.embeddings),2)) %>% cbind(cluster_info)
    cell_id <- rownames(tsne_2D)
    
    # prepare gene sets
    sctype = openxlsx::read.xlsx(db_)
    if(tissue %in% unique(sctype$tissueType)){
      gs_list = gene_sets_prepare(db_, tissue)
      es.max = sctype_score(scRNAseqData = scRNA_subtissue[["SCT"]]@scale.data, scaled = TRUE, 
                            gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
      for(res in c('SCT_snn_res.0.2','SCT_snn_res.0.4','SCT_snn_res.0.6','SCT_snn_res.0.8','SCT_snn_res.1')){
        cL_resutls = switch(res,
                            # 'SCT_snn_res.0.1' = do.call("rbind", lapply(unique(scRNA_subtissue@meta.data$SCT_snn_res.0.1), function(cl){
                            #   es.max.cl = sort(rowSums(es.max[ ,rownames(scRNA_subtissue@meta.data[scRNA_subtissue@meta.data$SCT_snn_res.0.1==cl, ])]), decreasing = !0)
                            #   head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNA_subtissue@meta.data$SCT_snn_res.0.1==cl)), 10)
                            # })),
                            'SCT_snn_res.0.2' = do.call("rbind", lapply(unique(scRNA_subtissue@meta.data$SCT_snn_res.0.2), function(cl){
                              es.max.cl = sort(rowSums(es.max[ ,rownames(scRNA_subtissue@meta.data[scRNA_subtissue@meta.data$SCT_snn_res.0.2==cl, ])]), decreasing = !0)
                              head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNA_subtissue@meta.data$SCT_snn_res.0.2==cl)), 10)
                            })),
                            # 'SCT_snn_res.0.3' = do.call("rbind", lapply(unique(scRNA_subtissue@meta.data$SCT_snn_res.0.3), function(cl){
                            #   es.max.cl = sort(rowSums(es.max[ ,rownames(scRNA_subtissue@meta.data[scRNA_subtissue@meta.data$SCT_snn_res.0.3==cl, ])]), decreasing = !0)
                            #   head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNA_subtissue@meta.data$SCT_snn_res.0.3==cl)), 10)
                            # })),
                            'SCT_snn_res.0.4' = do.call("rbind", lapply(unique(scRNA_subtissue@meta.data$SCT_snn_res.0.4), function(cl){
                              es.max.cl = sort(rowSums(es.max[ ,rownames(scRNA_subtissue@meta.data[scRNA_subtissue@meta.data$SCT_snn_res.0.4==cl, ])]), decreasing = !0)
                              head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNA_subtissue@meta.data$SCT_snn_res.0.4==cl)), 10)
                            })),
                            'SCT_snn_res.0.6' = do.call("rbind", lapply(unique(scRNA_subtissue@meta.data$SCT_snn_res.0.6), function(cl){
                              es.max.cl = sort(rowSums(es.max[ ,rownames(scRNA_subtissue@meta.data[scRNA_subtissue@meta.data$SCT_snn_res.0.6==cl, ])]), decreasing = !0)
                              head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNA_subtissue@meta.data$SCT_snn_res.0.6==cl)), 10)
                            })),
                            'SCT_snn_res.0.8' = do.call("rbind", lapply(unique(scRNA_subtissue@meta.data$SCT_snn_res.0.8), function(cl){
                              es.max.cl = sort(rowSums(es.max[ ,rownames(scRNA_subtissue@meta.data[scRNA_subtissue@meta.data$SCT_snn_res.0.8==cl, ])]), decreasing = !0)
                              head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNA_subtissue@meta.data$SCT_snn_res.0.8==cl)), 10)
                            })),
                            'SCT_snn_res.1' = do.call("rbind", lapply(unique(scRNA_subtissue@meta.data$SCT_snn_res.1), function(cl){
                              es.max.cl = sort(rowSums(es.max[ ,rownames(scRNA_subtissue@meta.data[scRNA_subtissue@meta.data$SCT_snn_res.1==cl, ])]), decreasing = !0)
                              head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNA_subtissue@meta.data$SCT_snn_res.1==cl)), 10)
                            }))
        )
        sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
        sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
        celltype_sub <- sctype_scores[,c(1,2)]
        colnames(celltype_sub) <- c(res,paste0(res,'_celltype'))
        celltype_sub=celltype_sub[!duplicated(celltype_sub),]
        celltype_sub=celltype_sub[!duplicated(celltype_sub[,1]),] ###20220705 fixed 根据第一列去重
        tsne_2D <- left_join(tsne_2D,celltype_sub,by=res)
      }
      rownames(tsne_2D) <- cell_id
    }else{
      tsne_2D <- tsne_2D
    }
    tsne_2D_list[[t]] <- tsne_2D
  }
  tsne_2D <- data.frame()
  for(i in 1:length(tsne_2D_list)){
    tsne_2D <- rbind(tsne_2D,tsne_2D_list[[i]])
  }
  write.csv(tsne_2D,file = paste0(rslt,'/',studyid,'_tsne_2D.csv'),row.names = TRUE)
}
output_3D <- function(scRNA,rslt,studyid){
  pct <- scRNA [["pca"]]@stdev / sum( scRNA [["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  pcs <- min(co1, co2)
  pc.num <- 1:pcs
  scRNA <- RunTSNE(scRNA,dims = pc.num,dim.embed = 3) %>% RunUMAP(dims = pc.num,n.components = 3L)
  tsne_3D_list <- list()
  for(t in unique(scRNA$tissue)){    
    tissue <- t
    scRNA_subtissue <- subset(scRNA, tissue == t)
    org <- unique(scRNA_subtissue$organism)
    res <-  c(0.2,0.4,0.6,0.8,1.0)
    cluster_info <- scRNA_subtissue@meta.data[,paste0("SCT_snn_res.",res)]
    #输出3D坐标
    tsne_3D <- round(as.data.frame(scRNA_subtissue[['tsne']]@cell.embeddings),2) %>% cbind(round(as.data.frame(scRNA_subtissue[['umap']]@cell.embeddings),2)) %>% cbind(cluster_info)
    cell_id <- rownames(tsne_3D)
    # prepare gene sets
    sctype = openxlsx::read.xlsx(db_)
    if(tissue %in% unique(sctype$tissueType)){
      gs_list = gene_sets_prepare(db_, tissue)
      es.max = sctype_score(scRNAseqData = scRNA_subtissue[["SCT"]]@scale.data, scaled = TRUE, 
                            gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
      for(res in c('SCT_snn_res.0.2','SCT_snn_res.0.4','SCT_snn_res.0.6','SCT_snn_res.0.8','SCT_snn_res.1')){
        cL_resutls = switch(res,
                            # 'SCT_snn_res.0.1' = do.call("rbind", lapply(unique(scRNA_subtissue@meta.data$SCT_snn_res.0.1), function(cl){
                            #   es.max.cl = sort(rowSums(es.max[ ,rownames(scRNA_subtissue@meta.data[scRNA_subtissue@meta.data$SCT_snn_res.0.1==cl, ])]), decreasing = !0)
                            #   head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNA_subtissue@meta.data$SCT_snn_res.0.1==cl)), 10)
                            # })),
                            'SCT_snn_res.0.2' = do.call("rbind", lapply(unique(scRNA_subtissue@meta.data$SCT_snn_res.0.2), function(cl){
                              es.max.cl = sort(rowSums(es.max[ ,rownames(scRNA_subtissue@meta.data[scRNA_subtissue@meta.data$SCT_snn_res.0.2==cl, ])]), decreasing = !0)
                              head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNA_subtissue@meta.data$SCT_snn_res.0.2==cl)), 10)
                            })),
                            # 'SCT_snn_res.0.3' = do.call("rbind", lapply(unique(scRNA_subtissue@meta.data$SCT_snn_res.0.3), function(cl){
                            #   es.max.cl = sort(rowSums(es.max[ ,rownames(scRNA_subtissue@meta.data[scRNA_subtissue@meta.data$SCT_snn_res.0.3==cl, ])]), decreasing = !0)
                            #   head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNA_subtissue@meta.data$SCT_snn_res.0.3==cl)), 10)
                            # })),
                            'SCT_snn_res.0.4' = do.call("rbind", lapply(unique(scRNA_subtissue@meta.data$SCT_snn_res.0.4), function(cl){
                              es.max.cl = sort(rowSums(es.max[ ,rownames(scRNA_subtissue@meta.data[scRNA_subtissue@meta.data$SCT_snn_res.0.4==cl, ])]), decreasing = !0)
                              head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNA_subtissue@meta.data$SCT_snn_res.0.4==cl)), 10)
                            })),
                            'SCT_snn_res.0.6' = do.call("rbind", lapply(unique(scRNA_subtissue@meta.data$SCT_snn_res.0.6), function(cl){
                              es.max.cl = sort(rowSums(es.max[ ,rownames(scRNA_subtissue@meta.data[scRNA_subtissue@meta.data$SCT_snn_res.0.6==cl, ])]), decreasing = !0)
                              head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNA_subtissue@meta.data$SCT_snn_res.0.6==cl)), 10)
                            })),
                            'SCT_snn_res.0.8' = do.call("rbind", lapply(unique(scRNA_subtissue@meta.data$SCT_snn_res.0.8), function(cl){
                              es.max.cl = sort(rowSums(es.max[ ,rownames(scRNA_subtissue@meta.data[scRNA_subtissue@meta.data$SCT_snn_res.0.8==cl, ])]), decreasing = !0)
                              head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNA_subtissue@meta.data$SCT_snn_res.0.8==cl)), 10)
                            })),
                            'SCT_snn_res.1' = do.call("rbind", lapply(unique(scRNA_subtissue@meta.data$SCT_snn_res.1), function(cl){
                              es.max.cl = sort(rowSums(es.max[ ,rownames(scRNA_subtissue@meta.data[scRNA_subtissue@meta.data$SCT_snn_res.1==cl, ])]), decreasing = !0)
                              head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNA_subtissue@meta.data$SCT_snn_res.1==cl)), 10)
                            }))
        )
        sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
        sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
        celltype_sub <- sctype_scores[,c(1,2)]
        colnames(celltype_sub) <- c(res,paste0(res,'_celltype'))
        celltype_sub=celltype_sub[!duplicated(celltype_sub),]
        celltype_sub=celltype_sub[!duplicated(celltype_sub[,1]),] ###20220705 fixed 根据第一列去重
        tsne_3D <- left_join(tsne_3D,celltype_sub,by=res)
      }
      rownames(tsne_3D) <- cell_id
    }else{
      tsne_3D <- tsne_3D
    }
    tsne_3D_list[[t]] <- tsne_3D
  }
  tsne_3D <- data.frame()
  for(i in 1:length(tsne_3D_list)){
    tsne_3D <- rbind(tsne_3D,tsne_3D_list[[i]])
  }
  write.csv(tsne_3D,file = paste0(rslt,'/',studyid,'_tsne_3D.csv'),row.names = TRUE)
}
output_3D_smart <- function(scRNA,rslt,studyid){
  pct <- scRNA [["pca"]]@stdev / sum( scRNA [["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  pcs <- min(co1, co2)
  pc.num <- 1:pcs
  perplexity <- floor((dim(scRNA)[2]-1)/3 - 1)
  scRNA <- RunTSNE(scRNA,dims = pc.num,dim.embed = 3,perplexity = perplexity,check_duplicates=FALSE) %>% RunUMAP(dims = pc.num,n.components = 3L)
  tsne_3D_list <- list()
  for(t in unique(scRNA$tissue)){
    # tissue <- t
    scRNA_subtissue <- subset(scRNA, tissue == t)
    org <- unique(scRNA_subtissue$organism)
    res <-  c(0.2,0.4,0.6,0.8,1.0)
    cluster_info <- scRNA_subtissue@meta.data[,paste0("SCT_snn_res.",res)]
    #输出3D坐标
    tsne_3D <- round(as.data.frame(scRNA_subtissue[['tsne']]@cell.embeddings),2) %>% cbind(round(as.data.frame(scRNA_subtissue[['umap']]@cell.embeddings),2)) %>% cbind(cluster_info)
    cell_id <- rownames(tsne_3D)
    # prepare gene sets
    sctype = openxlsx::read.xlsx(db_)
    if(tissue %in% unique(sctype$tissueType)){
      gs_list = gene_sets_prepare(db_, tissue)
      es.max = sctype_score(scRNAseqData = scRNA_subtissue[["SCT"]]@scale.data, scaled = TRUE, 
                            gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
      for(res in c('SCT_snn_res.0.2','SCT_snn_res.0.4','SCT_snn_res.0.6','SCT_snn_res.0.8','SCT_snn_res.1')){
        cL_resutls = switch(res,
                            # 'SCT_snn_res.0.1' = do.call("rbind", lapply(unique(scRNA_subtissue@meta.data$SCT_snn_res.0.1), function(cl){
                            #   es.max.cl = sort(rowSums(es.max[ ,rownames(scRNA_subtissue@meta.data[scRNA_subtissue@meta.data$SCT_snn_res.0.1==cl, ])]), decreasing = !0)
                            #   head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNA_subtissue@meta.data$SCT_snn_res.0.1==cl)), 10)
                            # })),
                            'SCT_snn_res.0.2' = do.call("rbind", lapply(unique(scRNA_subtissue@meta.data$SCT_snn_res.0.2), function(cl){
                              es.max.cl = sort(rowSums(es.max[ ,rownames(scRNA_subtissue@meta.data[scRNA_subtissue@meta.data$SCT_snn_res.0.2==cl, ])]), decreasing = !0)
                              head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNA_subtissue@meta.data$SCT_snn_res.0.2==cl)), 10)
                            })),
                            # 'SCT_snn_res.0.3' = do.call("rbind", lapply(unique(scRNA_subtissue@meta.data$SCT_snn_res.0.3), function(cl){
                            #   es.max.cl = sort(rowSums(es.max[ ,rownames(scRNA_subtissue@meta.data[scRNA_subtissue@meta.data$SCT_snn_res.0.3==cl, ])]), decreasing = !0)
                            #   head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNA_subtissue@meta.data$SCT_snn_res.0.3==cl)), 10)
                            # })),
                            'SCT_snn_res.0.4' = do.call("rbind", lapply(unique(scRNA_subtissue@meta.data$SCT_snn_res.0.4), function(cl){
                              es.max.cl = sort(rowSums(es.max[ ,rownames(scRNA_subtissue@meta.data[scRNA_subtissue@meta.data$SCT_snn_res.0.4==cl, ])]), decreasing = !0)
                              head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNA_subtissue@meta.data$SCT_snn_res.0.4==cl)), 10)
                            })),
                            'SCT_snn_res.0.6' = do.call("rbind", lapply(unique(scRNA_subtissue@meta.data$SCT_snn_res.0.6), function(cl){
                              es.max.cl = sort(rowSums(es.max[ ,rownames(scRNA_subtissue@meta.data[scRNA_subtissue@meta.data$SCT_snn_res.0.6==cl, ])]), decreasing = !0)
                              head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNA_subtissue@meta.data$SCT_snn_res.0.6==cl)), 10)
                            })),
                            'SCT_snn_res.0.8' = do.call("rbind", lapply(unique(scRNA_subtissue@meta.data$SCT_snn_res.0.8), function(cl){
                              es.max.cl = sort(rowSums(es.max[ ,rownames(scRNA_subtissue@meta.data[scRNA_subtissue@meta.data$SCT_snn_res.0.8==cl, ])]), decreasing = !0)
                              head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNA_subtissue@meta.data$SCT_snn_res.0.8==cl)), 10)
                            })),
                            'SCT_snn_res.1' = do.call("rbind", lapply(unique(scRNA_subtissue@meta.data$SCT_snn_res.1), function(cl){
                              es.max.cl = sort(rowSums(es.max[ ,rownames(scRNA_subtissue@meta.data[scRNA_subtissue@meta.data$SCT_snn_res.1==cl, ])]), decreasing = !0)
                              head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNA_subtissue@meta.data$SCT_snn_res.1==cl)), 10)
                            }))
        )
        sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
        sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
        celltype_sub <- sctype_scores[,c(1,2)]
        colnames(celltype_sub) <- c(res,paste0(res,'_celltype'))
        celltype_sub=celltype_sub[!duplicated(celltype_sub),]
        tsne_3D <- left_join(tsne_3D,celltype_sub,by=res)
      }
      rownames(tsne_3D) <- cell_id
    }else{
      tsne_3D <- tsne_3D
    }
    tsne_3D_list[[t]] <- tsne_3D
  }
  tsne_3D <- data.frame()
  for(i in 1:length(tsne_3D_list)){
    tsne_3D <- rbind(tsne_3D,tsne_3D_list[[i]])
  }
  write.csv(tsne_3D,file = paste0(rslt,'/',studyid,'_tsne_3D.csv'),row.names = TRUE)
}

###细胞占比输出
abundance <- function(scRNA,rslt,studyid){
  dir.create(paste0(rslt,'/cluster'))
  org <- unique(scRNA$organism)
  res <- c(0.2,0.4,0.6,0.8,1.0)
  cluster_list <- list()
  # cluster_list[['res_0.1']] <- scRNA$SCT_snn_res.0.1
  cluster_list[['res_0.2']] <- scRNA$SCT_snn_res.0.2
  # cluster_list[['res_0.3']] <- scRNA$SCT_snn_res.0.3
  cluster_list[['res_0.4']] <- scRNA$SCT_snn_res.0.4
  cluster_list[['res_0.6']] <- scRNA$SCT_snn_res.0.6
  cluster_list[['res_0.8']] <- scRNA$SCT_snn_res.0.8
  cluster_list[['res_1']] <- scRNA$SCT_snn_res.1
  for(clu in c(0.2,0.4,0.6,0.8,1.0)){
    table <- table(cluster_list[[paste0('res_',clu)]],scRNA$orig.ident)  ####fixed on 20220409
    abundance <- cbind(row.names(table),table)
    colnames(abundance)[1]="Cluster"
    write.csv(abundance, paste0(rslt,'/cluster/',studyid,"_SCT_snn_res.",clu,'_cluster_abundance.csv'), row.names = FALSE)
    cluster_present <- prop.table(as.matrix(table),1) #按行计算每个样本在不同celltype中的细胞比例
    cluster_present <- as.data.frame(cbind(Cluster=0:(length(levels(cluster_list[[paste0('res_',clu)]]))-1),
                                           cluster_present))
    write.csv(cluster_present, paste0(rslt,'/cluster/',studyid,"_SCT_snn_res.",clu,'_cluster_persent.csv'), row.names = FALSE)
  }
}
#按celltype输出 add on 20220615
abundance_ct <- function(scRNA,rslt,studyid){
  dir.create(paste0(rslt,'/celltype'))
  org <- unique(scRNA$organism)
  res <- c(0.2,0.4,0.6,0.8,1.0)
  celltypeinfo <- read_csv(paste0(rslt,'/',studyid,'_tsne_2D.csv')) %>% as.data.frame() %>% .[,c(1,11:15)] ##所有是从13:19
  celltype_list <- list()
  # celltype_list[['res_0.1']] <- as.factor(celltypeinfo[,2])
  celltype_list[['res_0.2']] <- as.factor(celltypeinfo[,2])
  # celltype_list[['res_0.3']] <- as.factor(celltypeinfo[,4])
  celltype_list[['res_0.4']] <- as.factor(celltypeinfo[,3])
  celltype_list[['res_0.6']] <- as.factor(celltypeinfo[,4])
  celltype_list[['res_0.8']] <- as.factor(celltypeinfo[,5])
  celltype_list[['res_1']] <- as.factor(celltypeinfo[,6])
  # names(celltype_list[['res_0.1']]) <- celltypeinfo[,1]
  names(celltype_list[['res_0.2']]) <- celltypeinfo[,1]
  # names(celltype_list[['res_0.3']]) <- celltypeinfo[,1]
  names(celltype_list[['res_0.4']]) <- celltypeinfo[,1]
  names(celltype_list[['res_0.6']]) <- celltypeinfo[,1]
  names(celltype_list[['res_0.8']]) <- celltypeinfo[,1]
  names(celltype_list[['res_1']]) <- celltypeinfo[,1]
  for(clu in c(0.2,0.4,0.6,0.8,1.0)){
    table <- table(celltype_list[[paste0('res_',clu)]],scRNA$orig.ident)  ####fixed on 20220409
    abundance <- cbind(row.names(table),table)
    colnames(abundance)[1]="Celltype"
    write.csv(abundance, paste0(rslt,'/celltype/',studyid,"_SCT_snn_res.",clu,'_celltype_abundance.csv'), row.names = FALSE)
    cluster_present <- prop.table(as.matrix(table),1) #按行计算每个样本在不同celltype中的细胞比例
    cluster_present <- as.data.frame(cbind(Celltype=levels(celltype_list[[paste0('res_',clu)]]),cluster_present))
    write.csv(cluster_present, paste0(rslt,'/celltype/',studyid,"_SCT_snn_res.",clu,'_celltype_persent.csv'), row.names = FALSE)
  }
}

###基因表达数据输出
output_exp <- function(scRNA,rslt,studyid){
  org <- unique(scRNA$organism)
  exp <- scRNA@assays$SCT@data %>% as.data.frame()  ###这里提供归一化处理的表达谱counts值（SCT） 20220622改为使用data矩阵
  exp$gene <- rownames(exp)
  exp_melt <- as.data.frame(reshape2::melt(exp,id.vars = "gene",variable.name="cell_name"))
  exp_melt_no0 <- exp_melt[exp_melt$value > 0,] ###去除0值
  exp_melt_no0 <- exp_melt_no0[complete.cases(exp_melt_no0[,3]),] ###去除空值
  write.csv(exp_melt_no0,file = paste0(rslt,'/',studyid,'_exp_melt_filter.csv'),row.names = F)
}
###marker基因输出
marker <- function(scRNA,rslt,studyid){
  org <- unique(scRNA$organism)
  res <- c(0.2,0.4,0.6,0.8,1.0)
  for(clu in c(0.2,0.4,0.6,0.8,1.0)){
    Idents(scRNA) <- paste0("SCT_snn_res.",clu)
    cluster_marker <- FindAllMarkers(scRNA, assay = 'SCT',#设置assay为RNA#
                                     slot = 'data', only.pos = TRUE, verbose = TRUE)
    write.csv(cluster_marker,file = paste0(rslt,'/cluster/',studyid,"_SCT_snn_res.",clu,"_marker.csv"),row.names = F)
  }
}
##按celltype输出
marker_ct <- function(scRNA,rslt,studyid){
  org <- unique(scRNA$organism)
  res <- c(0.2,0.4,0.6,0.8,1.0)
  for(clu in c(0.2,0.4,0.6,0.8,1.0)){
    Idents(scRNA) <- paste0("SCT_snn_res.",clu,"_celltype")
    celltype_marker <- FindAllMarkers(scRNA, assay = 'SCT',#设置assay为RNA#
                                      slot = 'data', only.pos = TRUE, verbose = TRUE)
    write.csv(celltype_marker,file = paste0(rslt,'/celltype/',studyid,"_SCT_snn_res.",clu,"_marker.csv"),row.names = F)
  }
}
####不同组（如果有的话）进行差异分析
diff <- function(scRNA,rslt,res,control,case,studyid){ #rslt <- "/data/panlab/maqinfeng/organoid/sc_outdir/celseq/step3_result/GSE113561"
  future::plan('multisession',workers=8) ###调用多核，4核
  scRNA_bak <- scRNA
  Idents(scRNA_bak) <- res
  all_diff <- data.frame()
  for(clu in 0:(length(table(Idents(scRNA_bak)))-1)){
    sub_object <- subset(scRNA_bak,idents = clu)
    Idents(sub_object) <- 'group'
    tryCatch({
      diff_df <- FindMarkers(sub_object,ident.1 = case, ident.2 = control,group.by = 'group')
      diff_df$clusters <- clu
      diff_df$resolution <- res
      diff_df = diff_df %>% rownames_to_column('gene')
      all_diff <- rbind(all_diff,diff_df)
    },error=function(e){
      # sink('D:/disease project/scrna_file/analysis/scrna_analysis/diff_errlog.txt', append = T)
      sink(paste0(errorpath,'diff_errlog.txt'),append = T)
      print(paste0(studyid,'_',res,'_',clu))
      sink()
    })
  }
  future::plan('multisession',workers=1) ###调用多核，4核
  write.csv(all_diff,paste0(rslt,'/cluster/',studyid,'_',res,'_diff.csv'),row.names = F)
}
####不同组（如果有的话）进行差异分析,按celltype
diff_ct <- function(scRNA,rslt,res,control,case,studyid){ #rslt <- "/data/panlab/maqinfeng/organoid/sc_outdir/celseq/step3_result/GSE113561"
  future::plan('multisession',workers=8) ###调用多核，4核
  scRNA_bak <- scRNA
  Idents(scRNA_bak) <- res
  all_diff <- data.frame()
  for(clu in names(table(Idents(scRNA_bak)))){
    sub_object <- subset(scRNA_bak,idents = clu)
    Idents(sub_object) <- 'group'
    tryCatch({
      diff_df <- FindMarkers(sub_object,ident.1 = case, ident.2 = control,group.by = 'group')
      diff_df$clusters <- clu
      diff_df$resolution <- res
      diff_df = diff_df %>% rownames_to_column('gene')
      all_diff <- rbind(all_diff,diff_df)
    },error=function(e){
      # sink('/public/home/panjianbo/mqf/scRNA_MQF/Seurat_pipline/diff_errlog.txt', append = T)
      sink(paste0(errorpath,'diff_errlog.txt'),append = T)
      print(paste0(studyid,'_',res,'_celltype',clu))
      sink()
    })
    future::plan('multisession',workers=1) ###调用多核，4核
  }
  write.csv(all_diff,paste0(rslt,'/celltype/',studyid,'_',res,'_diff.csv'),row.names = F)
}


###不同组，整体差异分析
diff_whole <- function(scRNA,rslt,control,case,studyid){
  future::plan('multisession',workers=8) ###调用多核，4核
  scRNA_bak <- scRNA
  Idents(scRNA_bak) <- 'group'
  whole_diff <- FindMarkers(scRNA_bak,ident.1 = case, ident.2 = control,group.by = 'group',logfc.threshold = -1, min.pct = 0)
  whole_diff = whole_diff %>% rownames_to_column('gene')
  future::plan('multisession',workers=1) ###调用多核，4核
  write.csv(whole_diff,paste0(rslt,'/',studyid,'_whole_diff.csv'),row.names = F)
}



organism_orgDb = data.frame(organism = c("Mus musculus","Homo sapiens"),
                            OrgDbs = c('org.Mm.eg.db','org.Hs.eg.db'))

output_files = function(rslt,studyid){
  allDiff = read_delim(paste0(rslt,'/',studyid,'_whole_diff.csv'))
  studyinfo <- dataset[dataset$DSAid == studyid,]
  # organism = studyinfo %>% dplyr::filter(Pertorgid==studyid) %>% .$Organism %>% gsub('_',' ',.)
  organism = studyinfo$Organism
  orgDb = organism_orgDb[match(organism,organism_orgDb$organism),2]
  file_trans = bitr(unlist(allDiff[,1]),"SYMBOL","ENTREZID",orgDb)
  colnames(allDiff)[1]=colnames(file_trans)[1]
  allDiff=inner_join(file_trans,allDiff,by='SYMBOL')
  allDiff=allDiff[-1]
  write.table(allDiff,paste0(rslt,'/',studyid,'_whole_diff.csv'),row.names = F,sep = "\t",quote = F)
  # output the most diff genes(abs(log2FC) top 500)
  mostDiff = allDiff %>% dplyr::filter(p_val < 0.05 & abs(avg_log2FC) > 0) %>% arrange(desc(abs(avg_log2FC))) %>% head(2000)
  # the most up
  upmostdiff = mostDiff %>% dplyr::filter(avg_log2FC > 0)
  write.table(upmostdiff,paste0(rslt,"/common_analysis/",studyid,"_upmostdiff.txt"),row.names = F,sep = "\t",quote = F)
  # the most down
  downmostdiff = mostDiff %>% dplyr::filter(avg_log2FC < 0)
  write.table(downmostdiff,paste0(rslt,"/common_analysis/",studyid,"_downmostdiff.txt"),row.names = F,sep = "\t",quote = F)
}

localdbpath = src  #enrich path
enrich_fun = function(rslt,studyid){
  if(file.exists(paste0(rslt,'/',studyid,'_whole_diff.csv'))){
    organism = dataset[dataset$DSAid == studyid,]$Organism
    orgDb = organism_orgDb[match(organism,organism_orgDb$organism),2]
    if(organism == 'Homo sapiens'){
      load(paste0(localdbpath,'GO_KEGG_Hs_ENTREZID_2022May18.RData'))
    }else{
      load(paste0(localdbpath,'GO_KEGG_Mm_ENTREZID_2022May18.RData'))
    }
    up = read.table(paste0(rslt,'/common_analysis/',studyid,'_upmostdiff.txt'),sep = '\t',header = T) %>% .$ENTREZID
    down = read.table(paste0(rslt,'/common_analysis/',studyid,'_downmostdiff.txt'),sep = '\t',header = T) %>% .$ENTREZID
    if(length(up) != 0 | length(down) != 0){
      DEGs = c(up,down)
      tryCatch({
        enrich_GO_plot(DEGs,studyid,GO_TERM2GENE,GO_TERM2NAME,GO_TERM2ONT,orgDb,rslt)
        GO = read_tsv(paste0(rslt,'/common_analysis/',studyid,'_GO.tsv'))
        BP = GO %>% dplyr::filter(ONT == "BP");CC = GO %>% dplyr::filter(ONT == "CC");MF = GO %>% dplyr::filter(ONT == "MF");
        p = BUBBLE(BP[1:10,],"GO")
        ggsave(paste0(rslt,"/common_analysis/pics/",studyid,"_BP.png"),plot = p,width = 12,height = 8)
        p = BUBBLE(CC[1:10,],"GO")
        ggsave(paste0(rslt,"/common_analysis/pics/",studyid,"_CC.png"),plot = p,width = 12,height = 8)
        p = BUBBLE(MF[1:10,],"GO")
        ggsave(paste0(rslt,"/common_analysis/pics/",studyid,"_MF.png"),plot = p,width = 12,height = 8)
      },error=function(e){
        sink(paste0(rslt,'/common_analysis/',studyid,'_GO.fail'))
        sink()
      })
      tryCatch({
        enrich_KEGG_plot(DEGs,studyid,KEGG_TERM2GENE,KEGG_TERM2NAME,orgDb,rslt)
        KEGG = read_tsv(paste0(rslt,"/common_analysis/",studyid,'_KEGG.tsv'))
        p = BUBBLE(KEGG[1:10,],"KEGG")
        ggsave(paste0(rslt,"/common_analysis/pics/",studyid,"_KEGG.png"),plot =p,width = 12,height = 8)
      },error=function(e){
        sink(paste0(rslt,"/common_analysis/",studyid,'_KEGG.fail'))
        sink()
      })
      tryCatch({
        enrich_marker(up,down,studyid,orgDb,organism,rslt)
        marker = read_tsv(paste0(rslt,"/common_analysis/",studyid,'_marker.tsv'))
      },error=function(e){
        sink(paste0(rslt,"/common_analysis/",studyid,'_marker.fail'))
        sink()
      })
      
      # p = BUBBLE(marker[1:10,],"KEGG")
      # ggsave(paste0('./result/',studyid,"/pics/",studyid,"_marker.png"),plot = p,width = 12,height = 8)
    }else{
      print('enrich fail: the most diff fail!')
    }
  }else{
    print('differential files not exist!')
  }
  ### GSEA ###
  rank = read.table(paste0(rslt,"/",studyid,"_whole_diff.csv"),header = T,sep = "\t")
  rankgenes = rank$avg_log2FC
  names(rankgenes) = rank$ENTREZID
  tryCatch({
    GSEA = clusterProfiler::GSEA(rankgenes %>% sort(.,decreasing = T),pvalueCutoff = 1,TERM2GENE = KEGG_TERM2GENE) %>% clusterProfiler::setReadable(.,orgDb,keyType="ENTREZID")
    GSEAresult = GSEA %>% .@result %>% left_join(KEGG_TERM2NAME,by = c('ID'='KEGGPATHID')) %>% dplyr::select(-Description) %>% dplyr::select(ID,NAME,everything()) %>% dplyr::rename(Description = 'NAME')
    GSEA@result = GSEAresult
    save(GSEA,file = paste0(rslt,'/common_analysis/',"GSEA.RData"))
    write_delim(GSEAresult,paste0(rslt,"/common_analysis/",studyid,'_GSEA_terms.txt'),delim = '\t')
  },error=function(e){
    sink(paste0(rslt,'/common_analysis/',studyid,'_GSEA.fail'))
    sink()
  })
}

enrich_GO_plot = function(DEGs,studyid,GO_TERM2GENE,GO_TERM2NAME,GO_TERM2ONT,orgDb,rslt){
  x <- clusterProfiler::enricher(DEGs,TERM2GENE = GO_TERM2GENE,pvalueCutoff = 1,qvalueCutoff = 1)
  x = clusterProfiler::setReadable(x,orgDb,keyType="ENTREZID") %>% dplyr::select(-Description)
  GOenrich_result = x@result %>% left_join(GO_TERM2NAME,by = c('ID'='PATHID')) %>% 
    left_join(GO_TERM2ONT,by = c('ID'='PATHID')) %>% dplyr::select(ID,NAME,ONT,everything()) %>% dplyr::rename(Description = 'NAME')
  write_tsv(GOenrich_result,paste0(rslt,"/common_analysis/",studyid,"_GO.tsv"))
  # GObubble(GOenrich_result_batch2,regulate,DEG_length)
}
enrich_KEGG_plot = function(DEGs,studyid,KEGG_TERM2GENE,KEGG_TERM2NAME,orgDb,rslt){
  x <- clusterProfiler::enricher(DEGs,TERM2GENE = KEGG_TERM2GENE,pvalueCutoff =1,qvalueCutoff = 1)
  x = clusterProfiler::setReadable(x,orgDb,keyType="ENTREZID") %>% dplyr::select(-Description)
  KEGGenrich_result = x@result %>% left_join(KEGG_TERM2NAME,by = c('ID'='KEGGPATHID')) %>% dplyr::select(ID,NAME,everything()) %>% dplyr::rename(Description = 'NAME')
  write_tsv(KEGGenrich_result,paste0(rslt,"/common_analysis/",studyid,"_KEGG.tsv"))
  # KEGGbubble(KEGGenrich_result_batch2,regulate)
}

enrich_marker = function(up,down,studyid,orgDb,organism,rslt){
  if(organism == 'Homo sapiens'){
    load(paste0(localdbpath,'PanglaoDBCellMarker_TERM2GENE_Hs.RData'))
  }else{
    load(paste0(localdbpath,'PanglaoDBCellMarker_TERM2GENE_Mm.RData'))
  }
  x_up = tryCatch({
    clusterProfiler::enricher(up,TERM2GENE = PanglaoDBCellMarker_TERM2GENE,minGSSize = 1,pvalueCutoff =1,qvalueCutoff = 1) %>% 
      clusterProfiler::setReadable(.,orgDb,keyType="ENTREZID") %>% .@result %>% separate(Description,c('Source','Description'),'_') %>% 
      add_column(Regulation = 'up',.before = 'GeneRatio')
  },error=function(e){
    data.frame()
  })
  x_down = tryCatch({
    clusterProfiler::enricher(down,TERM2GENE = PanglaoDBCellMarker_TERM2GENE,minGSSize = 1,pvalueCutoff =1,qvalueCutoff = 1) %>% 
      clusterProfiler::setReadable(.,orgDb,keyType="ENTREZID") %>% .@result %>% separate(Description,c('Source','Description'),'_') %>% 
      add_column(Regulation = 'down',.before = 'GeneRatio')
  },error=function(e){
    data.frame()
  })
  x = bind_rows(x_up,x_down)
  write_tsv(x,paste0(rslt,'/common_analysis/',studyid,"_marker.tsv"))
}

BUBBLE = function(enrich,type){
  if(toupper(type) == "GO" | toupper(type) == "KEGG"){
    enrich = enrich %>% na.omit()
    # enrich = enrich %>% separate(GeneRatio,c("enrichcount","genecount"),"/",2)
    enrich = enrich %>% separate(GeneRatio,c("enrichcount","genecount"),"/")
    DEG_length = enrich$genecount[1] %>% as.numeric()
    enrich$GeneRatio = apply(enrich %>% dplyr::select(enrichcount,genecount),1,function(x)as.numeric(x[1])/as.numeric(x[2]))
    enrich = enrich %>% arrange(desc(GeneRatio))
    x=enrich$GeneRatio
    y<-factor(enrich$Description,levels = rev(enrich$Description))
    ggplot(enrich,aes(x,y))+
      geom_point(aes(size=Count,color=p.adjust))+ 
      guides(shape = guide_legend(order=1,override.aes=list(size=3.5)),color = guide_colourbar(order=2),size = guide_legend(order=3))+
      scale_color_gradient(low = "red", high = "blue")+ 
      labs(color=expression(p.adjust),size="Count",x="GeneRatio",y="",title="")+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5),
            #去网格去背景色
            axis.line = element_line(colour = "black"), axis.text = element_text(color = "black",size = 14),
            #刻度字体大小
            legend.text = element_text(size = 14),legend.title=element_text(size=14),
            axis.title.x = element_text(size = 14))+scale_size_continuous(range=c(4,8))
  }else{
    print("Invalid data type input!")
  }
}

plot_GSEA = function(rslt,studyid){
  dir.create(paste0(rslt,"/common_analysis/",'pics/GSEA'))
  load(paste0(rslt,"/common_analysis/","GSEA.RData"))
  GSEA@result = GSEA@result %>% dplyr::filter(pvalue < 0.05)
  tryCatch({
    for (i in 1:dim(GSEA@result)[1]) {
      pathID = GSEA@result$ID[i]
      theme_update(plot.subtitle = element_text(hjust = 0.5))
      p = gseaplot2(GSEA,geneSetID = i,title = '',color = "green")
      ggsave(paste0(rslt,"/common_analysis/",'pics/GSEA/',pathID,'.png'))
    }
  },error=function(e){
    sink(paste0(rslt,'/common_analysis/',studyid,'_GSEA_plot.fail'))
    sink()
  })
}

