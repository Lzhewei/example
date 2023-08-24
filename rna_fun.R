
###有生物学重复
diff_fun <- function(control,case,organism,studyid,profiles){
  library(DESeq2)
  group = org_SRR_GSE_GSM[org_SRR_GSE_GSM$GSM %in% c(control,case),]
  #%>% distinct(.,SRR,.keep_all = T)
  if (length(setdiff(group$SRR,colnames(profiles[[organism]]))) ==0) {
    profile <- profiles[[organism]] %>% dplyr::select(c('GeneID','Length',group$SRR)) #注意列名对不对得上！！！
    if(profile[,3] %>% table() %>% as_tibble() %>% dim() %>% .[1] != 1){
      colnames(profile)[3:(dim(group)[1]+2)] = profile %>% dplyr::select(starts_with('SRR')) %>% colnames() %>% as.data.frame() %>% dplyr::rename(SRR = '.') %>% left_join(group) %>% .$GSM
      # colnames(profile)[3:(dim(group)[1]+2)] = profile %>% dplyr::select(starts_with('ERR')) %>% colnames() %>% as.data.frame() %>% dplyr::rename(SRR = '.') %>% left_join(group) %>% .$GSM
      genelength = profile %>% dplyr::select(c('GeneID','Length'))
      profile = profile %>% na.omit() %>% dplyr::select(c('GeneID','Length',control,case)) %>% column_to_rownames('GeneID') %>% dplyr::select(-Length)
      #profile = profile %>% dplyr::select(c('GeneID','Length',control,case)) %>% column_to_rownames('GeneID') %>% dplyr::select(-Length)
      ###去除所有样本都为0（没有检测到表达的基因）的基因
      profile <- profile[rowSums(profile)>0,]
      #更新genelength列表
      genelength <- subset(genelength,genelength$GeneID %in% rownames(profile))
      ###这一步很关键，要明白condition这里是因子，不是样本名称
      condition <- factor(c(rep("control",length(control)),rep("case",length(case))),levels = c("control","case"))
      colData <- data.frame(row.names = colnames(profile),condition)
      ##构建dds对象，开始DESeq流程，dds = DESeqDataSet Object
      dds <- DESeqDataSetFromMatrix(profile, colData, design = ~condition)
      dds <- DESeq(dds,minReplicatesForReplace = 40)
      ###查看结果
      res <- results(dds) %>% as.data.frame() %>% rownames_to_column('GeneID') %>% mutate(GeneID = as.numeric(GeneID)) %>% left_join(genelength[!duplicated(genelength),])
      ##计算TPM矩阵
      len <- genelength$Length
      kb <- len / 1000
      RPKM <- profile / kb
      TPM <- t(t(RPKM)/colSums(RPKM) * 1000000) %>% as.data.frame()
      TPM$GeneID <- rownames(TPM)
      TPM <- TPM[,c(ncol(TPM),1:(ncol(TPM)-1))]
      #计算中位数和平均值
      controlmeans <- TPM %>% dplyr::select(all_of(control)) %>% rowMeans()
      casemeans <- TPM %>% dplyr::select(all_of(case)) %>% rowMeans()
      controlmedian <- TPM %>% dplyr::select(all_of(control)) %>% as.matrix() %>% rowMedians()
      casemedian <- TPM %>% dplyr::select(all_of(case)) %>% as.matrix() %>% rowMedians()
      res = res %>% add_column(AveExpr_Control = controlmeans,AveExpr_Case = casemeans,MidExpr_Control = controlmedian,MidExpr_Case = casemedian,.after = 'GeneID')
      ###将TPM矩阵中每一列的表达值进行转化，只分成control和case两列，一组多个样本用“;” 隔开
      TPM = TPM %>% 
        unite(col = !!paste0(control,collapse = ';'),all_of(control),sep = ';', remove = TRUE) %>% 
        unite(col = !!paste0(case,collapse = ';'),all_of(case),sep = ';', remove = TRUE)
      pvalue=0.05
      allDiff = res
      regulate = allDiff %>% filter(pvalue < 0.05 & abs(log2FoldChange) > 0) %>% arrange(desc(abs(log2FoldChange))) %>% head(1000)
      active_genes = regulate %>% filter(log2FoldChange > 0) %>% .$GeneID
      negtive_genes = regulate %>% filter(log2FoldChange < 0) %>% .$GeneID
      allDiff$regulation = 'nosig'
      if(length(active_genes) != 0){
        allDiff[allDiff$GeneID %in% active_genes,]$regulation = 'up'
      }
      if(length(negtive_genes) != 0){
        allDiff[allDiff$GeneID %in% negtive_genes,]$regulation = 'down'
      }
      allDiff <- allDiff %>% dplyr::select(GeneID,everything(),-baseMean,-lfcSE,-stat,-Length)
      # add colume neglog10p
      allDiff$neglog10p = -log10(allDiff$pvalue)
      # unite colnames
      allDiff = allDiff %>% dplyr::rename(Log2FC = 'log2FoldChange',PValue = 'pvalue',AdjPValue = 'padj')
      
      dir.create(paste0(resultdir,studyid))
      dir.create(paste0(resultdir,studyid,"/pics"))
      # sort logFC for GSEA
      allDiff = allDiff %>% arrange(desc(Log2FC))
      # output alldiff file
      write.table(allDiff,paste0(resultdir,studyid,"/",studyid,"_alldiff.txt"),row.names = F,sep = "\t",quote = F)
      write.table(TPM,paste0(resultdir,studyid,"/",studyid,"_exp.txt"),row.names = F,sep = "\t",quote = F)
    }else{
      write_delim(c(studyid) %>% as.data.frame(),paste0(errpath,'/study_noexpression_failnewrun.txt'),append = T,delim = '\t')
    }
  }else{
    fail_err=setdiff(group$SRR,colnames(profiles[[organism]]))%>% paste0(.,collapse = ';')
    write_delim(cbind(studyid,study$Accession,fail_err) %>% as.data.frame(),paste0(errpath,'/study_err_notexist_failnewrun.txt'),append = T,delim = '\t')
  }
  
}

###无生物学重复
diff_fun_U <- function(control,case,organism,studyid,profiles){
  library(edgeR)
  group = org_SRR_GSE_GSM[org_SRR_GSE_GSM$GSM %in% c(control,case),]
  if (length(setdiff(group$SRR,colnames(profiles[[organism]]))) ==0) {
    profile <- profiles[[organism]] %>% dplyr::select(c('GeneID','Length',group$SRR))
    # some profile may got zero expression each gene, which cannot used to analysis
    if(profile[,3] %>% table() %>% as_tibble() %>% dim() %>% .[1] != 1){
      colnames(profile)[3:(dim(group)[1]+2)] = profile %>% dplyr::select(starts_with('SRR')) %>% colnames() %>% as.data.frame() %>% dplyr::rename(SRR = '.') %>% left_join(group) %>% .$GSM
      # colnames(profile)[3:(dim(group)[1]+2)] = profile %>% dplyr::select(starts_with('ERR')) %>% colnames() %>% as.data.frame() %>% dplyr::rename(SRR = '.') %>% left_join(group) %>% .$GSM
      genelength = profile %>% dplyr::select(c('GeneID','Length'))
      profile = profile %>% na.omit() %>% dplyr::select(c('GeneID','Length',control,case)) %>% column_to_rownames('GeneID') %>% dplyr::select(-Length)
      #profile = profile %>% dplyr::select(c('GeneID','Length',control,case)) %>% column_to_rownames('GeneID') %>% dplyr::select(-Length)
      
      #构建DEGList
      groups <- factor(c(rep("control",length(control)),rep("case",length(case))),levels = c("control","case"))
      y <- DGEList(counts = profile, group = groups)
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
      
      ##计算TPM矩阵
      len <- genelength$Length
      kb <- len / 1000
      RPKM <- profile / kb
      TPM <- t(t(RPKM)/colSums(RPKM) * 1000000) %>% as.data.frame()
      
      #计算中位数和平均值
      controlmeans <- TPM %>% dplyr::select(all_of(control)) %>% rowMeans()
      casemeans <- TPM %>% dplyr::select(all_of(case)) %>% rowMeans()
      controlmedian <- TPM %>% dplyr::select(all_of(control)) %>% as.matrix() %>% rowMedians()
      casemedian <- TPM %>% dplyr::select(all_of(case)) %>% as.matrix() %>% rowMedians()
      TPM = TPM %>% rownames_to_column('GeneID')
      res = right_join(tibble(GeneID = TPM$GeneID,AveExpr_Control = controlmeans,AveExpr_Case = casemeans,MidExpr_Control = controlmedian,MidExpr_Case = casemedian),etab)
      
      ###将TPM矩阵中每一列的表达值进行转化，只分成control和case两列，一组多个样本用“;” 隔开
      TPM = TPM %>% 
        unite(col = !!paste0(control,collapse = ';'),all_of(control),sep = ';', remove = TRUE) %>% 
        unite(col = !!paste0(case,collapse = ';'),all_of(case),sep = ';', remove = TRUE)
      pvalue=0.05
      #
      allDiff = res
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
      # add colume neglog10p
      allDiff$neglog10p = -log10(allDiff$PValue)
      dir.create(paste0(resultdir,studyid))
      dir.create(paste0(resultdir,studyid,"/pics"))
      # sort logFC for GSEA
      allDiff = allDiff %>% arrange(desc(Log2FC))
      # output alldiff file
      write.table(allDiff,paste0(resultdir,studyid,"/",studyid,"_alldiff.txt"),row.names = F,sep = "\t",quote = F)
      write.table(TPM,paste0(resultdir,studyid,"/",studyid,"_exp.txt"),row.names = F,sep = "\t",quote = F)
    }else{
      write_delim(c(studyid) %>% as.data.frame(),paste0(errpath,'/study_noexpression_failnewrun.txt'),append = T,delim = '\t')
    }
  }else{
    fail_err=setdiff(group$SRR,colnames(profiles[[organism]]))%>% paste0(.,collapse = ';')
    write_delim(cbind(studyid,study$Accession,fail_err) %>% as.data.frame(),'./study_err_notexist_failnewrun.txt',append = T,delim = '\t')
    write_delim(cbind(studyid,study$Accession,fail_err) %>% as.data.frame(),paste0(errpath,'/study_err_notexist_failnewrun.txt'),append = T,delim = '\t')
  }
  
}
