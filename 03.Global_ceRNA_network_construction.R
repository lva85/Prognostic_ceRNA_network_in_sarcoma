### set workspace
setwd("~/SARC_proj/")

load("classification.RData")

##load gene expression data
if(T){
	gene_exp.mat = read.table("TCGA_SARC_geneExp.txt",header = T,sep = "\t")
	gene_exp.df = gene_exp.mat[!duplicated(gene_exp.mat$gene),]
	row.names(gene_exp.df) = gene_exp.df$gene
	gene_exp.df = gene_exp.df[,-c(1:2)]
	#retrieve the tumor samples
	gene_exp.tumor.df = gene_exp.df[,which(substring(colnames(gene_exp.df),14,15) == "01")]
	#normalization
	gene_exp.tumor.df = gene_exp.tumor.df[which(apply(gene_exp.tumor.df,1,sum)>0),]
	gene_exp.tumor.df.normalized <- t(scale(t(gene_exp.tumor.df)))

}

##construct the global ceRNA network based on the classification
  if(T){
    ## DEGs between C1 and C2
    if(T){
      library(limma) 
      condition_table <- result.mat
      condition_table_for_limma <- model.matrix(~0+ condition_table$result)
      colnames(condition_table_for_limma) <- gsub("condition_table\\$result","", colnames(condition_table_for_limma))
      lmfit <- lmFit(gene_exp.tumor.df.normalized, condition_table_for_limma)
      
      contrast.matrix<-makeContrasts(cluster2-cluster1,levels=condition_table_for_limma)
      fit2<-contrasts.fit(lmfit,contrast.matrix)
      fit2<-eBayes(fit2)
      tT=topTable(fit2,adjust='fdr',coef=1,number=Inf,p.value=1)
      DEGs_cluster_HighTOLow <- tT[which(abs(tT$logFC) > 1 & tT$adj.P.Val < 0.05),]
      dim(DEGs_cluster_HighTOLow)
      
    }
    
	###GO/KEGG for DEGs
	if(T){
	library(ggplot2)
	library(stringr)
	library(clusterProfiler)

	DEGs_cluster_HighTOLow.df <- bitr(row.names(DEGs_cluster_HighTOLow), fromType="SYMBOL",
							  toType="ENTREZID", 
							  OrgDb = "org.Hs.eg.db")

	##GO
	DEGs_cluster_HighTOLow.go <- enrichGO(gene = DEGs_cluster_HighTOLow.df$ENTREZID, 
								  keyType = "ENTREZID",
								  OrgDb = "org.Hs.eg.db",
								  ont="BP",
								  pAdjustMethod = "BH",
								  pvalueCutoff  = 0.05,
								  qvalueCutoff  = 0.05,
								  readable      = TRUE)

	write.table(DEGs_cluster_HighTOLow.go[DEGs_cluster_HighTOLow.go$pvalue < 0.05 & DEGs_cluster_HighTOLow.go$qvalue < 0.05],
				paste0('GO_for_DEGs_cluster_HighTOLow.ClusterProfile.txt'),row.names = F,sep = "\t",quote = F)

	##KEGG
	library(DO.db)
	DEGs_cluster_HighTOLow.ego <- enrichKEGG(
	  gene = DEGs_cluster_HighTOLow.df$ENTREZID,
	  keyType = "kegg",
	  organism  = 'hsa',
	  pvalueCutoff  = 0.05,
	  pAdjustMethod  = "BH",
	  qvalueCutoff  = 0.05
	)

	DEGs_cluster_HighTOLow.ego.2=DOSE::setReadable(DEGs_cluster_HighTOLow.ego, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
	
	write.table(DEGs_cluster_HighTOLow.ego.2[DEGs_cluster_HighTOLow.ego.2$pvalue < 0.05 & DEGs_cluster_HighTOLow.ego.2$qvalue < 0.05],
				"KEGG_for_DEGs_cluster_HighTOLow.ClusterProfile.txt",row.names = F,sep = "\t",quote = F)

	sig.n = dim(DEGs_cluster_HighTOLow.ego.2[DEGs_cluster_HighTOLow.ego.2$pvalue < 0.05 & DEGs_cluster_HighTOLow.ego.2$qvalue < 0.05])[1]

	
	ggplot(DEGs_cluster_HighTOLow.ego.2@result[1:20,],aes(x = GeneRatio,y = Description))+
	  geom_point(aes(color = p.adjust,
					 size = Count))+
	  scale_size(range=c(2, 14)) + 
	  scale_color_gradient(low = "red", high = "blue")+
	  xlab("GeneRatio")+
	  ylab("") + 
	  theme_bw()+
	  guides(
		color = guide_colorbar(reverse = TRUE))
	}

	##annotate the mRNA/lncRNA in the DE result
	if(T){
	  #load gtf annotation file（download from GENCODE）
      # BiocManager::install("rtracklayer")
      # BiocManager::install("SummarizedExperiment")
      library(rtracklayer)
      library(SummarizedExperiment)
      
      gtf <- rtracklayer::import('gencode.v38.annotation.gtf')
      gtf_df <- as.data.frame(gtf)
      
      gtf_df.gene = gtf_df[which(gtf_df$type == "gene"),]
      gtf_df.gene = gtf_df.gene[!duplicated(gtf_df.gene$gene_name),]
      
	  DEGs_cluster_HighTOLow.type = data.frame("gene"=row.names(DEGs_cluster_HighTOLow),
                                               "id"=gtf_df.gene$gene_id[match(row.names(DEGs_cluster_HighTOLow),gtf_df.gene$gene_name)],
                                               "type"=gtf_df.gene$gene_type[match(row.names(DEGs_cluster_HighTOLow),gtf_df.gene$gene_name)])
      DElncRNAs = DEGs_cluster_HighTOLow.type[which(DEGs_cluster_HighTOLow.type$type == "lncRNA"),]
      DEgenes = DEGs_cluster_HighTOLow.type[which(DEGs_cluster_HighTOLow.type$type == "protein_coding"),]
      write.table(DElncRNAs,
                  "04.ceRNA_interaction_network/01.DElncRNA.txt",row.names=F,col.names = F, sep = "\t",quote = F)
      
      DEGs_cluster_HighTOLow.type.ceRNA = DEGs_cluster_HighTOLow.type[which(DEGs_cluster_HighTOLow.type$type %in% c("lncRNA","protein_coding")),]
      
      table(DEGs_cluster_HighTOLow.type.ceRNA$type)
      # lncRNA protein_coding 
      # 127           1645 
	  
	}
	
    ##DE miRNA
    if(T){
      miRNA_exp = read.table("SARC_dataPrep.csv",header = T,sep = ",")
      miRNA_exp.df = miRNA_exp[,grep("reads_per_million_miRNA_mapped_TCGA",colnames(miRNA_exp))]
      colnames(miRNA_exp.df) = gsub("reads_per_million_miRNA_mapped_","",colnames(miRNA_exp.df))
      row.names(miRNA_exp.df) = miRNA_exp$miRNA_ID
      dim(miRNA_exp.df)
      # [1] 1881  259
      
      ##filtration, normalization
      miRNA_exp.df.normalized = miRNA_exp.df[which(apply(miRNA_exp.df,1,sum)>0),]
      dim(miRNA_exp.df.normalized)
      # [1] 1469  259
      miRNA_exp.df.normalized <- t(scale(t(miRNA_exp.df.normalized)))
      
      ## retrive the samples in the classification result
      miRNA_exp.df.normalized = miRNA_exp.df.normalized[,which(substring(colnames(miRNA_exp.df.normalized),1,12) %in% survival.info$PATIENT_ID)]
      
      library(limma)
      condition_table <- data.frame("sample" = colnames(miRNA_exp.df.normalized),
                                    "result" = survival.info$cluster[match(substring(colnames(miRNA_exp.df.normalized),1,12),survival.info$PATIENT_ID)])
      condition_table_for_limma <- model.matrix(~0+ condition_table$result)
      colnames(condition_table_for_limma) <- gsub("condition_table\\$result","", colnames(condition_table_for_limma))
      lmfit <- lmFit(miRNA_exp.df.normalized, condition_table_for_limma)
      
      contrast.matrix<-makeContrasts(cluster2-cluster1,levels=condition_table_for_limma)
      fit2<-contrasts.fit(lmfit,contrast.matrix)
      fit2<-eBayes(fit2)
      tT=topTable(fit2,adjust='fdr',coef=1,number=Inf,p.value=1)
      DE_miRNA_cluster_HighTOLow <- tT[which(tT$adj.P.Val < 0.05),]
      dim(DE_miRNA_cluster_HighTOLow)
      # [1] 171   6
      
       
    }
    
    ## DElncRNAs-DEmiRNA interaction 
    if(T){
      ## because the miRNA in the miRCode database is miRNA family, so we need to link the miRNA family with mature miRNA， miRBaseConverter package was used to obtain the miRNA family information
      # BiocManager::install("miRBaseConverter")
      library(miRBaseConverter)
	  
	  ## load hsa_hairpin2mature (this file is the miRNA mature and haipin pair information obtained from miRBase database)
      h2m = read.table("hsa_hairpin2mature.txt",header = F,sep = "\t")
	  
      miRNANames=h2m$V2
      version=checkMiRNAVersion(miRNANames,verbose = FALSE)
      result=miRNA_NameToAccession(miRNANames,version=version)
      Accessions=result$Accession
      Family_Info=checkMiRNAFamily(Accessions)
      Family_Info2 = Family_Info[which(!is.na(Family_Info$Family)),]
      Family_Info2 = Family_Info2[!duplicated(Family_Info2),]
      
      ##because the miRNA family information in the miRCode database is collapsed, so we used home script to uncollapse it for further match
      DElncRNA_miRNA_interaction.uniq.miRNAfamily = read.table("DElncRNA_miRNA_interaction.uniq.miRNAfamily.uniq.txt",header = F)
      
      ### obtain the pre-miRNA through the miRNA family information
      DElncRNA_miRNA_interaction.final = merge(DElncRNA_miRNA_interaction.uniq.miRNAfamily,
                                               Family_Info2,by.x="V5",by.y="Family",all.x=T)
      
      
      DElncRNA_miRNA_interaction.final = DElncRNA_miRNA_interaction.final[which(!is.na(DElncRNA_miRNA_interaction.final$miRNAName_v21)),]
      dim(DElncRNA_miRNA_interaction.final)
      # [1] 6718    8
      length(unique(DElncRNA_miRNA_interaction.final$V4))
      # [1] 95
      
      DElncRNAs$id2 = substring(DElncRNAs$id,1,15)
      
      DElncRNA_miRNA_interaction.final2 = merge(DElncRNA_miRNA_interaction.final,DElncRNAs,by.x="V4",by.y="id2",all.x=T)
      
      
      ###obtain the interaction pairs for DElncRNAs
      DElncRNA_miRNA_interaction.final3 = DElncRNA_miRNA_interaction.final2[,c("gene","miRNAName_v21")]
      colnames(DElncRNA_miRNA_interaction.final3) = c("DElncRNA","miRNA")
      
      ###obtain the interaction pairs for DElncRNAs and DEmiRNAs
      DElncRNA_DEmiRNA_interaction = DElncRNA_miRNA_interaction.final3[which(DElncRNA_miRNA_interaction.final3$miRNA %in%
                                                                               row.names(DE_miRNA_cluster_HighTOLow) ),]
      
	  ##obtain the mature miRNAs for the pre-miRNAs in the DElncRNA-DEmiRNA interaction pairs for further matching the DEmRNA-DEmiRNAs interaction pairs
      DEpremiRNA.left = unique(DElncRNA_DEmiRNA_interaction$DEmiRNA)
      DEmaturemiRNA.left = h2m[which(h2m$V2 %in% DEpremiRNA.left),] #[1] 99
      
    }
    
    ## DEmiRNA-DEmRNA interaction
    if(T){
      ## obtain the shared DEmRNA-DEmiRNA interaction pairs from miRDB, TargetScan, miRWalk using home-script
      DEmaturemiRNA.left.target.common = read.table("04.ceRNA_interaction_network/DEmaturemiRNA.left.target.common.txt",header = F,sep = "\t")
      
      ##obtain the pre-miRNA information
      DEmaturemiRNA.left.target.common.matchPre = merge(DEmaturemiRNA.left.target.common,
                                                        h2m,by.x="V1",by.y="V1",all.x=T)
      
      table(unique(DEmaturemiRNA.left.target.common.matchPre$V2.y) %in% unique(DEmaturemiRNA.left$V2))
      #54
      
      ##further filter out the DEmRNA-DEmiRNA interaction pairs
      DEmiRNA_DEmRNA_interaction = DEmaturemiRNA.left.target.common.matchPre[which(DEmaturemiRNA.left.target.common.matchPre$V2.x %in% 
                                                                                     DEgenes$gene ),]
      colnames(DEmiRNA_DEmRNA_interaction) = c("DEmature_miRNA","DEmRNA","DEpre_miRNA")
      
      
    }
    
    ##organize the DElncRNA-DEmiRNA interaction pairs, DEmiRNA-DEmRNA interaction pairs
    if(T){
      ##one miRNA in the DElncRNA-DEmiRNA interaction pairs doesn't have the DEmRNA target, will remove it from DElncRNA-DEmiRNA interaction pairs
      unique(DElncRNA_DEmiRNA_interaction$DEmiRNA) %in% unique(DEmiRNA_DEmRNA_interaction$DEpre_miRNA)
      # hsa-mir-320a
      
      DElncRNA_DEmiRNA_interaction2 = DElncRNA_DEmiRNA_interaction[which(DElncRNA_DEmiRNA_interaction$DEmiRNA %in%
                                                                           DEmiRNA_DEmRNA_interaction$DEpre_miRNA     ),]
      
    }
    
    ## merge DElncRNA-DEmiRNA-DEmRNA interaction pairs for Cytoscape visulization
    if(T){
      DElncRNA_DEmiRNA_DEmRNA_interaction =data.frame(
        "node1" = c(DElncRNA_DEmiRNA_interaction2$DElncRNA,DEmiRNA_DEmRNA_interaction$DEmRNA),
        "node2" = c(DElncRNA_DEmiRNA_interaction2$DEmiRNA,DEmiRNA_DEmRNA_interaction$DEpre_miRNA)
      )
      
      node_type = data.frame("node" = c(unique(DElncRNA_DEmiRNA_interaction2$DElncRNA),
                                        unique(DEmiRNA_DEmRNA_interaction$DEmRNA),
                                        unique(c(DEmiRNA_DEmRNA_interaction$DEpre_miRNA,
                                                 DElncRNA_DEmiRNA_interaction2$DEmiRNA)) ),
                             "type" = c(rep("lncRNA",length(unique(DElncRNA_DEmiRNA_interaction2$DElncRNA)) ),
                                        rep("mRNA",length(unique(DEmiRNA_DEmRNA_interaction$DEmRNA)) ),
                                        rep("miRNA",length(unique(c(DEmiRNA_DEmRNA_interaction$DEpre_miRNA,
                                                                    DElncRNA_DEmiRNA_interaction2$DEmiRNA))) ) )
      )
      
      node_type$reg = apply(node_type,1,function(x){
        if(x[1] %in% row.names(DEGs_cluster_HighTOLow)){
          ifelse(DEGs_cluster_HighTOLow[x[1],"logFC"] >0,"up","down" )
        }else if(x[1] %in% row.names(DE_miRNA_cluster_HighTOLow)){
          ifelse(DE_miRNA_cluster_HighTOLow[x[1],"logFC"] >0,"up","down" )
        }
      }) 
      
      
      write.table(DElncRNA_DEmiRNA_DEmRNA_interaction,"DElncRNA_DEmiRNA_DEmRNA_interaction.txt",row.names = F,sep = "\t",quote = F)
      write.table(node_type,"DElncRNA_DEmiRNA_DEmRNA_nodeType.txt",row.names = F,sep = "\t",quote = F)
      
      table(node_type$type)
      # lncRNA  miRNA   mRNA 
      # 94     54    993 
    }
    
  }
  
  
  
