output_files = function(studyid,control,case){
  allDiff = read_delim(paste0(resultdir,studyid,'/',studyid,'_alldiff.txt'))
  # output the most diff genes(abs(log2FC) top 1000)
  mostDiff = allDiff %>% filter(PValue < 0.05 & abs(Log2FC) > 0) %>% arrange(desc(abs(Log2FC))) %>% head(1000)
  # the most up
  upmostdiff = mostDiff %>% filter(Log2FC > 0)
  write.table(upmostdiff,paste0(resultdir,studyid,"/",studyid,"_upmostdiff.txt"),row.names = F,sep = "\t",quote = F)
  # the most down
  downmostdiff = mostDiff %>% filter(Log2FC < 0)
  write.table(downmostdiff,paste0(resultdir,studyid,"/",studyid,"_downmostdiff.txt"),row.names = F,sep = "\t",quote = F)
}

# enrich_GO_plot = function(DEGs,studyid,GO_TERM2GENE,GO_TERM2NAME,GO_TERM2ONT,orgDb){
#   x <- clusterProfiler::enricher(DEGs,TERM2GENE = GO_TERM2GENE,pvalueCutoff = 1,qvalueCutoff = 1)
#   x = clusterProfiler::setReadable(x,orgDb,keyType="ENTREZID") %>% dplyr::select(-Description)
#   GOenrich_result = x@result %>% left_join(GO_TERM2NAME,by = c('ID'='PATHID')) %>% 
#     left_join(GO_TERM2ONT,by = c('ID'='PATHID')) %>% dplyr::select(ID,NAME,ONT,everything()) %>% dplyr::rename(Description = 'NAME')
#   write_tsv(GOenrich_result,paste0(resultdir,studyid,'/',studyid,"_GO.tsv"))
# }
# enrich_KEGG_plot = function(DEGs,studyid,KEGG_TERM2GENE,KEGG_TERM2NAME,orgDb){
#   x <- clusterProfiler::enricher(DEGs,TERM2GENE = KEGG_TERM2GENE,pvalueCutoff =1,qvalueCutoff = 1)
#   x = clusterProfiler::setReadable(x,orgDb,keyType="ENTREZID") %>% dplyr::select(-Description)
#   KEGGenrich_result = x@result %>% left_join(KEGG_TERM2NAME,by = c('ID'='KEGGPATHID')) %>% dplyr::select(ID,NAME,everything()) %>% dplyr::rename(Description = 'NAME')
#   write_tsv(KEGGenrich_result,paste0(resultdir,studyid,'/',studyid,"_KEGG.tsv"))
# }

enrich_marker = function(up,down,studyid,orgDb){
  if(orgDb == 'org.Hs.eg.db'){
    load(paste0(localdbpath,'PanglaoDBCellMarker_TERM2GENE_Hs.RData'))
  }else{
    load(paste0(localdbpath,'PanglaoDBCellMarker_TERM2GENE_Mm.RData'))
  }
  x_up = clusterProfiler::enricher(up,TERM2GENE = PanglaoDBCellMarker_TERM2GENE,minGSSize = 1,pvalueCutoff =1,qvalueCutoff = 1) %>% 
    clusterProfiler::setReadable(.,orgDb,keyType="ENTREZID") %>% .@result %>% separate(col=Description,into=c('Source','Description'),sep='_') %>% 
    add_column(Regulation = 'up',.before = 'GeneRatio')
  x_down = clusterProfiler::enricher(down,TERM2GENE = PanglaoDBCellMarker_TERM2GENE,minGSSize = 1,pvalueCutoff =1,qvalueCutoff = 1) %>% 
    clusterProfiler::setReadable(.,orgDb,keyType="ENTREZID") %>% .@result %>% separate(col=Description,into=c('Source','Description'),sep='_') %>% 
    add_column(Regulation = 'down',.before = 'GeneRatio')
  x = bind_rows(x_up,x_down)
  write_tsv(x,paste0(resultdir,studyid,'/',studyid,"_marker.tsv"))
}

enrich_tf = function(up,down,studyid,orgDb){
  if(orgDb == 'org.Hs.eg.db'){
    load(paste0(localdbpath,'TRRUSTv2_TERM2GENE_Hs.RData'))
    
  }else{
    load(paste0(localdbpath,'TRRUSTv2_TERM2GENE_Mm.RData'))
  }
  x_up = clusterProfiler::enricher(up,TERM2GENE = tfenrich_TERM2GENE,minGSSize = 1,pvalueCutoff =1,qvalueCutoff = 1) %>% 
    clusterProfiler::setReadable(.,orgDb,keyType="ENTREZID") %>% .@result %>% separate(col=Description,into=c('Source','Description'),sep='_') %>% 
    add_column(Regulation = 'up',.before = 'GeneRatio')
  x_down = clusterProfiler::enricher(down,TERM2GENE = tfenrich_TERM2GENE,minGSSize = 1,pvalueCutoff =1,qvalueCutoff = 1) %>% 
    clusterProfiler::setReadable(.,orgDb,keyType="ENTREZID") %>% .@result %>% separate(col=Description,into=c('Source','Description'),sep='_') %>% 
    add_column(Regulation = 'down',.before = 'GeneRatio')
  x = bind_rows(x_up,x_down)
  write_tsv(x,paste0(resultdir,studyid,'/',studyid,"_tf.tsv"))
}

BUBBLE = function(enrich,type){
  if(toupper(type) == "GO" | toupper(type) == "KEGG"){
    enrich = enrich %>% na.omit()
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

### GO ### KEGG ###
enrich_fun = function(studyid){
  organism = group_list %>% filter(.data[['DSAid']]==studyid) %>% .$Organism
  organism=gsub('_',' ',organism)
  orgDb = organism_orgDb[match(organism,organism_orgDb$organism),2]
  # load(paste0(localdbpath,list.files(localdbpath) %>% grep(str_split(orgDb,'\\.',simplify = T)[2],.,value = T)))
  load(paste0(localdbpath,list.files(localdbpath) %>% grep('GO_KEGG',.,value = T) %>% grep(str_split(orgDb,'\\.',simplify = T)[2],.,value = T)))
  # up = read.table(paste0(resultdir,studyid,'/',studyid,'_upmostdiff.txt'),sep = '\t',header = T) %>% .$GeneID
  # down = read.table(paste0(resultdir,studyid,'/',studyid,'_downmostdiff.txt'),sep = '\t',header = T) %>% .$GeneID
  up = read_delim(paste0(resultdir,studyid,'/',studyid,'_upmostdiff.txt'),delim = '\t') %>% .$GeneID
  down = read_delim(paste0(resultdir,studyid,'/',studyid,'_downmostdiff.txt'),delim = '\t') %>% .$GeneID
  if(length(up) != 0 | length(down) != 0){
    DEGs = c(up,down)
    tryCatch({
      enrich_GO_plot(DEGs,studyid,GO_TERM2GENE,GO_TERM2NAME,GO_TERM2ONT,orgDb)
      GO = read_tsv(paste0(resultdir,studyid,'/',studyid,'_GO.tsv'))
      BP = GO %>% filter(ONT == "BP");CC = GO %>% filter(ONT == "CC");MF = GO %>% filter(ONT == "MF");
      p = BUBBLE(BP[1:10,],"GO")
      ggsave(paste0(resultdir,studyid,"/pics/",studyid,"_BP.png"),plot = p,width = 12,height = 8)
      p = BUBBLE(CC[1:10,],"GO")
      ggsave(paste0(resultdir,studyid,"/pics/",studyid,"_CC.png"),plot = p,width = 12,height = 8)
      p = BUBBLE(MF[1:10,],"GO")
      ggsave(paste0(resultdir,studyid,"/pics/",studyid,"_MF.png"),plot = p,width = 12,height = 8)
    },error=function(e){
      sink(paste0(resultdir,studyid,'/',studyid,'_GO.fail'))
      sink()
    })
    tryCatch({
      enrich_KEGG_plot(DEGs,studyid,KEGG_TERM2GENE,KEGG_TERM2NAME,orgDb)
      KEGG = read_tsv(paste0(resultdir,studyid,'/',studyid,'_KEGG.tsv'))
      p = BUBBLE(KEGG[1:10,],"KEGG")
      ggsave(paste0(resultdir,studyid,"/pics/",studyid,"_KEGG.png"),plot =p,width = 12,height = 8)
    },error=function(e){
      sink(paste0(resultdir,studyid,'/',studyid,'_KEGG.fail'))
      sink()
    })
    tryCatch({
      enrich_marker(up,down,studyid,orgDb)
    },error=function(e){
      sink(paste0(resultdir,studyid,'/',studyid,'_marker.fail'))
      sink()
    })
    tryCatch({
      enrich_tf(up,down,studyid,orgDb)
      # marker = read_tsv(paste0(resultdir,studyid,'/',studyid,'_marker.tsv'))
    },error=function(e){
      sink(paste0(resultdir,studyid,'/',studyid,'_tf.fail'))
      sink()
    })
    ### GSEA ###
  }else{
    print('enrich fail: the most diff fail!')
  }
  rank = read_delim(paste0(resultdir,studyid,"/",studyid,"_alldiff.txt"),delim = '\t')
  rankgenes = rank$Log2FC
  names(rankgenes) = rank$GeneID
  tryCatch({
    GSEA = clusterProfiler::GSEA(rankgenes %>% sort(.,decreasing = T),pvalueCutoff = 1,TERM2GENE = KEGG_TERM2GENE) %>% clusterProfiler::setReadable(.,orgDb,keyType="ENTREZID")
    GSEAresult = GSEA %>% .@result %>% left_join(KEGG_TERM2NAME,by = c('ID'='KEGGPATHID')) %>% dplyr::select(-Description) %>% dplyr::select(ID,NAME,everything()) %>% dplyr::rename(Description = 'NAME')
    GSEA@result = GSEAresult
    save(GSEA,file = paste0(resultdir,studyid,"/GSEA.RData"))
    write_delim(GSEAresult,paste0(resultdir,studyid,"/",studyid,'_GSEA_terms.txt'),delim = '\t')
  },error=function(e){
    sink(paste0(resultdir,studyid,'/',studyid,'_GSEA.fail'))
    sink()
  })
}

plot_GSEA = function(studyid){
  dir.create(paste0(resultdir,studyid,'/pics/GSEA'))
  load(paste0(resultdir,studyid,'/GSEA.RData'))
  GSEA@result = GSEA@result %>% filter(pvalue < 0.05)
  for (i in 1:dim(GSEA@result)[1]) {
    pathID = GSEA@result$ID[i]
    p = enrichplot::gseaplot2(GSEA,geneSetID = i,title = '',color = "green")
    png(paste0(resultdir,studyid,'/pics/GSEA/',pathID,'.png'),width = 500,height = 500)
    print(p)
    dev.off()
  }
}

