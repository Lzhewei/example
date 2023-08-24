


organism_orgDb = data.frame(organism = c("Homo sapiens","Mus musculus"),
                            OrgDbs = c('org.Hs.eg.db','org.Mm.eg.db'))
##### down stream files generate #####
### log2 判断 (by GEO2R)
iflog2 = function(exprSet){
  ex <- exprSet
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) {
    ex[ex <= 0] <- NaN
    exprSet <- log2(ex)
    print('log2 transform finished')
  }else{print('log2 transform not needed')}
  return(exprSet)
}

# ### main function
# main = function(studyid){
#   study = group_list %>% filter(.data[['DSAid']] == studyid)
#   GSE = study$Accession
#   GPL = study$Platform
#   if(file.exists(paste0(readdir,GSE,"_",GPL,".txt"))){
#     control = unlist(str_split(study$Control,";"))
#     case = unlist(str_split(study$Case,";"))
#     if(length(control) == 1 & length(case) == 1 ){
#       # write_delim(study,'./study_singlesample0531.txt',append = T,delim = '\t')
#       diff_fun_nr(studyid,GSE,control,case,resultdir,GPL)
#     }else{
#       diff_fun(studyid,GSE,control,case,resultdir,GPL)
#     }
#   }else{
#     write_delim(study,paste0(errpath,'study_notexist.txt'),append = T,delim = '\t')
#   }
# }

### diff & generate raw files
diff_fun = function(studyid,GSE,control,case,resultdir,GPL){
  class<-factor(c(rep("control",length(control)),rep("case",length(case))))
  exprSet = read_delim(paste0(readdir,GSE,'_',GPL,".txt"),delim = "\t")
  exprSet[exprSet == 'null'] <- NA
  exprSet[exprSet == 'NULL'] <- NA
  exprSet = exprSet %>% dplyr::select(ENTREZID,all_of(control),all_of(case)) %>% 
    na.omit() %>% dplyr::mutate_all(.,as.numeric) %>% 
    dplyr::group_by(ENTREZID) %>%   
    dplyr::summarise(across(everything(),mean)) %>%
    filter(!is.na(ENTREZID)) %>%
    column_to_rownames(var = 'ENTREZID')
  
  if(dim(exprSet)[1] == 0){
  }
  # add average/median expression
  addcol = exprSet %>% dplyr::select(all_of(control)) %>% rowMeans(na.rm = T) %>% as.data.frame() %>% rownames_to_column('GENEID') %>% dplyr::rename(AveExpr_Control='.') %>% 
    full_join(exprSet %>% dplyr::select(all_of(case)) %>% rowMeans(na.rm = T) %>% as.data.frame() %>% rownames_to_column('GENEID') %>% dplyr::rename(AveExpr_Case='.')) %>% 
    full_join(exprSet %>% dplyr::select(all_of(control)) %>% as.matrix() %>% Biobase::rowMedians(na.rm = T) %>% as.data.frame() %>% add_column(GENEID = rownames(exprSet)) %>% dplyr::rename(MidExpr_Control='.')) %>% 
    full_join(exprSet %>% dplyr::select(all_of(case)) %>% as.matrix() %>% Biobase::rowMedians(na.rm = T) %>% as.data.frame() %>% add_column(GENEID = rownames(exprSet)) %>% dplyr::rename(MidExpr_Case='.'))
  log2_exprSet = iflog2(exprSet)
  ### limma
  design=model.matrix(~0+class)
  dimnames(design) = list(c(control,case),c("case","control"))
  contrast.matrix <- makeContrasts('case-control',levels = design)
  fit=lmFit(log2_exprSet,design) %>% contrasts.fit(., contrast.matrix) %>% eBayes(.)
  allDiff=topTable(fit,coef=1,number=Inf) %>% rownames_to_column("GENEID")
  # bind addcol
  allDiff = allDiff %>% full_join(addcol) %>% dplyr::select(GENEID,starts_with('Ave'),starts_with('Mid'),everything())
  dir.create(paste0(resultdir,studyid))
  dir.create(paste0(resultdir,studyid,"/pics"))
  # add colume regulation(threhold pvalue 0.05 log2FC top1000)
  regulate = allDiff %>% filter(P.Value < 0.05 & abs(logFC) > 0) %>% arrange(desc(abs(logFC))) %>% head(1000)
  active_genes = regulate %>% filter(logFC > 0) %>% .$GENEID
  negtive_genes = regulate %>% filter(logFC < 0) %>% .$GENEID
  allDiff$regulation = 'nosig'
  if(length(active_genes) != 0){
    allDiff[allDiff$GENEID %in% active_genes,]$regulation = 'up'
  }
  if(length(negtive_genes) != 0){
    allDiff[allDiff$GENEID %in% negtive_genes,]$regulation = 'down'
  }
  # add colume neglog10p
  allDiff$neglog10p = -log10(allDiff$P.Value)
  # delete columes
  allDiff = allDiff %>% dplyr::select(-c('t','B','AveExpr'))
  # sort logFC for GSEA
  allDiff = allDiff %>% arrange(desc(logFC))
  # unite colnames
  allDiff = allDiff %>% dplyr::rename(GeneID = 'GENEID',Log2FC = 'logFC',PValue = 'P.Value',AdjPValue = 'adj.P.Val')
  # add fold change
  # output alldiff file
  write.table(allDiff,paste0(resultdir,studyid,"/",studyid,"_alldiff.txt"),row.names = F,sep = "\t",quote = F)
  log2_exprSet = log2_exprSet %>% 
    unite(col = !!paste0(control,collapse = '.'),all_of(control),sep = ';', remove = TRUE) %>% 
    unite(col = !!paste0(case,collapse = '.'),all_of(case),sep = ';', remove = TRUE)
  log2_exprSet = log2_exprSet %>% rownames_to_column("GeneID")
  write.table(log2_exprSet,paste0(resultdir,studyid,"/",studyid,"_exp.txt"),row.names = F,sep = "\t",quote = F)
}

# no replicates
diff_fun_nr = function(studyid,GSE,control,case,resultdir,GPL){
  class<-factor(c(rep("control",length(control)),rep("case",length(case))))
  ### 选择表达谱
  exprSet = read.table(paste0(readdir,GSE,'_',GPL,".txt"),sep = "\t",header = T) %>% 
    dplyr::select(ENTREZID,all_of(control),all_of(case)) %>% 
    group_by(ENTREZID) %>% 
    summarise(across(everything(),mean)) %>% 
    column_to_rownames(var = 'ENTREZID')
  if(dim(exprSet)[1] == 0){
    stop(paste0(studyid,' expression profile is empty'))
  }
  log2_exprSet = iflog2(exprSet)
  # edgeR test
  exprSet = exprSet %>% dplyr::select(all_of(control),all_of(case))
  #构建DEGList
  groups <- factor(c(rep("control",length(control)),rep("case",length(case))),levels = c("control","case"))
  # drop negative
  exprSet = exprSet[exprSet[,1] > 0 & exprSet[,2] > 0,]
  if(dim(exprSet)[1] != 0){
    y <- DGEList(counts = exprSet, group = groups)
    #数据过滤，保留至少在一个样本里有表达的基因（CPM > 1）
    keep <- rowSums(cpm(y)>1) >= 1
    y <- y[keep, ,keep.lib.sizes = FALSE]
    #标准化
    y <- calcNormFactors(y)
    #差异分析，设置bcv，人使用0.4，小鼠设为0.1
    bcv = 0.1
    # if(species == 'homo_sapiens'){bcv = 0.4}else{bcv = 0.1}
    et <- exactTest(y, dispersion = bcv ^ 2)
    etab <- as.data.frame(et$table)
    etab$p.adjust <- p.adjust(etab$PValue, method="fdr")
    etab <- etab[,c('logFC','PValue','p.adjust')]
    colnames(etab) <- c('Log2FC','PValue','AdjPValue')
    etab = etab %>% rownames_to_column('GeneID')
    # add average/median expression
    addcol = exprSet %>% dplyr::select(all_of(case)) %>% rowMeans(na.rm = T) %>% as.data.frame() %>% rownames_to_column('GENEID') %>% dplyr::rename(AveExpr_Case='.') %>% 
      full_join(exprSet %>% dplyr::select(all_of(control)) %>% rowMeans(na.rm = T) %>% as.data.frame() %>% rownames_to_column('GENEID') %>% dplyr::rename(AveExpr_Control='.')) %>% 
      full_join(exprSet %>% dplyr::select(all_of(case)) %>% as.matrix() %>% Biobase::rowMedians(na.rm = T) %>% as.data.frame() %>% add_column(GENEID = rownames(exprSet)) %>% dplyr::rename(MidExpr_Case='.')) %>% 
      full_join(exprSet %>% dplyr::select(all_of(control)) %>% as.matrix() %>% Biobase::rowMedians(na.rm = T) %>% as.data.frame() %>% add_column(GENEID = rownames(exprSet)) %>% dplyr::rename(MidExpr_Control='.'))
    allDiff = etab %>% full_join(addcol,by = c('GeneID'='GENEID')) %>% dplyr::select(GeneID,starts_with('Ave'),starts_with('Mid'),everything())
    # if(allDiff %>% filter(PValue < 0.05) %>% dim() %>% .[1] != 0){
    dir.create(paste0(resultdir,studyid))
    dir.create(paste0(resultdir,studyid,"/pics"))
    # add colume regulation(threhold pvalue 0.05 log2FC top1000)
    regulate = allDiff %>% filter(PValue < 0.05 & abs(Log2FC) > 0) %>% arrange(desc(abs(Log2FC))) %>% head(1000)
    active_genes = regulate %>% filter(Log2FC > 0) %>% .$GeneID
    negtive_genes = regulate %>% filter(Log2FC < 0) %>% .$GeneID
    allDiff$regulation = 'nosig'
    if(length(active_genes) != 0){
      allDiff[allDiff$GeneID %in% active_genes,]$regulation = 'up'
    }
    if(length(negtive_genes) != 0){
      allDiff[allDiff$GeneID %in% negtive_genes,]$regulation = 'down'
    }
    allDiff = allDiff %>% arrange(desc(Log2FC))
    write.table(allDiff,paste0(resultdir,studyid,"/",studyid,"_alldiff.txt"),row.names = F,sep = "\t",quote = F)
    log2_exprSet = log2_exprSet %>% dplyr::select(all_of(control),all_of(case)) %>% rownames_to_column("GeneID")
    write.table(log2_exprSet,paste0(resultdir,studyid,"/",studyid,"_exp.txt"),row.names = F,sep = "\t",quote = F)
  }else{
    print('sth wrong with profile!')
    write_delim(c(studyid) %>% as.data.frame(),paste0(errpath,'study_profilefail.txt'),append = T,delim = '\t')
  }
  
}

