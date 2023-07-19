### set workspace
setwd("~/SARC_proj/")

##load 45 immune-related immune signatures
gene <- read.csv("ImmuneSignatures.csv")[, 1:2]
genelist <- read.table("ImmuneRelatedGeneList.txt",header = T,sep = "\t")
colnames(genelist)
genelist <- genelist[,c("Symbol","Category")]
colnames(genelist) <- colnames(gene)
genes <- rbind(gene,genelist)
ImmuneRelated_list<- split(as.matrix(genes)[,1], genes[,2])
    
##load gene expression data
gene_exp.mat = read.table("TCGA_SARC_geneExp.txt",header = T,sep = "\t")
gene_exp.df = gene_exp.mat[!duplicated(gene_exp.mat$gene),]
row.names(gene_exp.df) = gene_exp.df$gene
gene_exp.df = gene_exp.df[,-c(1:2)]
#retrieve the tumor samples
gene_exp.tumor.df = gene_exp.df[,which(substring(colnames(gene_exp.df),14,15) == "01")]


#perform ssGSEA for 45 immune-related signatures among sarcoma patients
if(T){
  library(GSVA)
  gsva<- gsva(as.matrix(gene_exp.tumor.df), ImmuneRelated_list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
  gsva1<- t(scale(t(gsva)))
}
    
# 分类
if(T){
  library(cluster)
  # install.packages("dendextend"); 
  library(dendextend)
  
  #clust and heatmap
  col_clust <- hclust(dist(t(gsva1), method = 'euclidean'), method = 'complete')
  plot(col_clust,labels=FALSE)
  
  plot(as.dendrogram(col_clust), cex = 0.6, leaflab = "none")
  
  library(ggplot2)
  dend <- color_branches(col_clust, k = 3,col=c("green","blue","red"))
  hclust.p = ggplot(as.ggdend(dend),labels = FALSE,theme = theme_minimal()) +
    labs(x="Sample", y="Height")+
    theme_bw()+
    scale_y_continuous(expand = c(0, 0))+
    scale_x_continuous(expand = c(0, 1))+
    coord_cartesian(ylim=c(0,27))+
    theme(axis.text.x = element_text(color="Black", size=8,hjust=0.5,vjust = 0),
          axis.text.y = element_text(color ="Black", size = 8, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.title.x = element_text(face="bold", color ="Black", size = 10, angle = 0, hjust = .5, vjust = 0.5),
          axis.title.y = element_text(face="bold",color ="Black", size = 10, angle = 90, hjust = 0.5, vjust = 1.5),
          legend.position=c(0.88,0.75), legend.key = element_blank(), legend.background = element_blank(),
          legend.key.size = unit(6, "pt"),
          legend.text = element_text(color ="Black", size = 9, angle = 0, hjust = 0, vjust = 0.5),
          legend.title = element_text(face="bold",color ="Black", size = 10, angle = 0, hjust = .5, vjust = .5),
          strip.text = element_text(face="bold", color ="Black", size = 12)
    )
  
  
  dend <- color_branches(col_clust, k = 3,col=c("green","blue","red") ,groupLabels=c("C1(Low)","C3(Moderate)","C2(High)")) 
  plot(dend)
  
  library(gplots)
  z <- as.matrix(gsva1)
  heatmap.2( z,
             col=bluered,
             trace="none",key=T,keysize=1,scale = "none",
             Colv=dend,margins=c(8,16),cexRow = 0.8,labCol = NA,
             dendrogram="both",Rowv=T)
  
  #get clusters
  result <- cutree(col_clust,k=3)
  table(result)
  result.mat <- data.frame(result)
  result.mat <- data.frame("sample" = row.names(result.mat),result.mat)
  result.mat$result <- ifelse(result.mat$result== 1,"cluster1",ifelse(result.mat$result== 2,"cluster2","cluster3"))
  
  
}
    
    
## validate the result using consensus cluster analysis
if(T){
  # BiocManager::install("ConsensusClusterPlus")
  library(ConsensusClusterPlus)
  
  res <- ConsensusClusterPlus(as.matrix(gsva1), maxK = 8, reps = 1000, pItem = 0.8, pFeature = 1, distance = "euclidean",
                              clusterAlg = "km", corUse = "everything", seed=123456, innerLinkage = "complete", finalLinkage = "complete",
                              plot="pdf", writeTable=T,title = "03.untitled_consensus_cluster")
  icl <- calcICL(res, title = "03.untitled_consensus_cluster",plot = 'pdf')
  col_clust2 <- res[[3]][['consensusTree']]
  
  plot(col_clust2)
  
  dend2 <- color_branches(as.dendrogram(col_clust2), k = 3,col=c("blue","red","green") ,groupLabels=c("C3","C2","C1")) 
  plot(dend2)
  
  library(pheatmap)
  z <- as.matrix(gsva1)
  heatmap.2( z,
             #col=bluered(length(seq(-5,5,0.01))-1),breaks = seq(-5,5,0.01),
             col=bluered,
             trace="none",key=T,keysize=1,scale = "none",
             Colv=as.dendrogram(col_clust2),margins=c(8,16),cexRow = 0.8,labCol = NA,
             dendrogram="both",Rowv=T)
  
  result2 <- data.frame(res[[3]][['consensusClass']])
  result.mat2 <- data.frame("sample" = row.names(result2),result2)
  colnames(result.mat2)[2] <- "result"
  result.mat2$result <- ifelse(result.mat2$result== 1,"cluster1",ifelse(result.mat2$result== 2,"cluster2","cluster3"))
}

## limma was used to identify the differentially infiltrated signatures
if(T){
  library(limma)
  #cancer=1, normal=0. 
  condition_table <- result.mat
  condition_table_for_limma <- model.matrix(~0+ condition_table$result)
  colnames(condition_table_for_limma) <- gsub("condition_table\\$result","", colnames(condition_table_for_limma))
  lmfit <- lmFit(gsva1, condition_table_for_limma)
  
  contrast.matrix<-makeContrasts(cluster1-cluster2,levels=condition_table_for_limma)
  fit2<-contrasts.fit(lmfit,contrast.matrix)
  fit2<-eBayes(fit2)
  fit2$coefficients
  tT=topTable(fit2,adjust='fdr',coef=1,number=Inf,p.value=1)
  sig_cluster1TO2 <- tT[which(abs(tT$logFC) > 1 & tT$adj.P.Val < 0.05),]
  dim(sig_cluster1TO2)
  # [1] 43  6
  
  contrast.matrix<-makeContrasts(cluster1-cluster3,levels=condition_table_for_limma)
  fit2<-contrasts.fit(lmfit,contrast.matrix)
  fit2<-eBayes(fit2)
  fit2$coefficients
  tT=topTable(fit2,adjust='fdr',coef=1,number=Inf,p.value=1)
  sig_cluster1TO3 <- tT[which(abs(tT$logFC) > 1 & tT$adj.P.Val < 0.05),]
  dim(sig_cluster1TO3)
  # [1] 29  6
  
  contrast.matrix<-makeContrasts(cluster2-cluster3,levels=condition_table_for_limma)
  fit2<-contrasts.fit(lmfit,contrast.matrix)
  fit2<-eBayes(fit2)
  fit2$coefficients
  tT=topTable(fit2,adjust='fdr',coef=1,number=Inf,p.value=1)
  sig_cluster2TO3 <- tT[which(abs(tT$logFC) > 1 & tT$adj.P.Val < 0.05),]
  dim(sig_cluster2TO3)
  # [1] 23  6
  
  sig_gsva <- gsva1[unique(c(row.names(sig_cluster1TO2),row.names(sig_cluster1TO3),row.names(sig_cluster2TO3))),]
  z2 <- as.matrix(sig_gsva)
  heatmap.2( z2,
			 col=bluered(length(seq(-3,3,0.01))-1),breaks = seq(-3,3,0.01),
			 # col=bluered,
			 trace="none",key=T,keysize=1,scale = "none",
			 Colv=dend,labCol = NA,margins=c(8,16),cexRow = 0.8,
			 dendrogram="both",Rowv=T)
			 
}

save.image("classification.RData")    