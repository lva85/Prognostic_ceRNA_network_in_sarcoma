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

#DEGs_diff.df
cemiRNA.exp = miRNA_exp.df.normalized[row.names(miRNA_exp.df.normalized) %in% node_type$node[which(node_type$type == "miRNA")],]

cemiRNA.exp.t = cbind("patient_id" = substring(colnames(cemiRNA.exp),1,12),
					"sample" = colnames(cemiRNA.exp),t(cemiRNA.exp))
cemiRNA.exp.add.SurInfo = merge(survival.info,cemiRNA.exp.t,
							  by.x = "PATIENT_ID",by.y="patient_id",all.y=T)

colnames(cemiRNA.exp.add.SurInfo)[70:ncol(cemiRNA.exp.add.SurInfo)] = gsub("-","_",colnames(cemiRNA.exp.add.SurInfo)[70:ncol(cemiRNA.exp.add.SurInfo)])

##transform the expression matrix into numeric
cemiRNA.exp.add.SurInfo[,70:ncol(cemiRNA.exp.add.SurInfo)] = lapply(cemiRNA.exp.add.SurInfo[,70:ncol(cemiRNA.exp.add.SurInfo)],as.numeric)

cemiRNA.exp.add.SurInfo = cemiRNA.exp.add.SurInfo[which(cemiRNA.exp.add.SurInfo$OS_MONTHS<120),]

####loop to group the expression of miRNAs into high and low, further to analyze the OS assoication
cemiRNA_exp.L2H.survival = list()
cemiRNA_exp.L2H.survival.sig = list()
cemiRNA_exp.L2H.survival.sig.id = list()
cemiRNA_exp.L2H.survival.sig.p = list()
j=1
for(i in 1:length(colnames(cemiRNA.exp.add.SurInfo)[70:ncol(cemiRNA.exp.add.SurInfo)])){ # i=9
id = colnames(cemiRNA.exp.add.SurInfo)[70:ncol(cemiRNA.exp.add.SurInfo)][i]
iexp = cemiRNA.exp.add.SurInfo[,id]

imedian = median(iexp)
idf.mat = cemiRNA.exp.add.SurInfo
idf.mat$exp = ifelse(idf.mat[,id] <= imedian,0,1)

library("survival")
library("survminer")

ifit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ exp, data=idf.mat)

survdiff.h.l = survdiff(Surv(OS_MONTHS, OS_STATUS) ~ exp, data = idf.mat)
pvalue.h.l = 1 - pchisq(survdiff.h.l$chisq, length(survdiff.h.l$n) -1)

ip <- ggsurvplot(ifit, data = idf.mat,
				 pval = F,
				 pval.coord = c(0, 0.03),
				 size = 0.2,
				 palette=c("cornflowerblue","brown3"),  
				 legend = c(0.8,0.75),
				 legend.labs=c(paste0("Low (n=",nrow(idf.mat[idf.mat$exp == 0,]),")"),
							   paste0("High (n=",nrow(idf.mat[idf.mat$exp == 1,]),")") ), 
				 legend.title= "",
				 title="", 
				 ylab="Overall survival (%)",
				 xlab = " Time (Months)", 
				 censor.shape = 124,censor.size = 1,conf.int = FALSE, 
				 font.legend = c(3,"bold", "black"),
				 font.main = c(4, "bold", "black"),
				 font.x = c(3, "bold", "black"),
				 font.y = c(3, "bold", "black"),
				 break.x.by = 40,
				 font.tickslab = c(3,"black")
)

ip$plot = ip$plot + ggplot2::annotate("text",x = 55, y = 1,color="black",size = 1,fontface="bold",
									  label = paste0(gsub("_","-",id)," ( p=",format(pvalue.h.l, scientific = T,digits = 3)," )"))+
  theme(
	axis.ticks.x=element_line(size=0.15),
	axis.ticks.y=element_line(size=0.15),
		axis.line.x=element_line(size=0.2),
		axis.line.y=element_line(size=0.2)
  )

cemiRNA_exp.L2H.survival[[i]] = ip

if(pvalue.h.l < 0.05 ){
  cemiRNA_exp.L2H.survival.sig[[j]] = ip
  cemiRNA_exp.L2H.survival.sig.id[j] = id
  cemiRNA_exp.L2H.survival.sig.p[j] = format(pvalue.h.l, scientific = T,digits = 3)
  j = j + 1
}

}

arrange_ggsurvplots(cemiRNA_exp.L2H.survival.sig,  ncol = 3, nrow = 3,newpage = FALSE)

cemiRNA_exp.L2H.survival.sig.res = data.frame("miRNA" = gsub("_","-",unlist(cemiRNA_exp.L2H.survival.sig.id)),
                                                    "pvalue" = unlist(cemiRNA_exp.L2H.survival.sig.p))


### further construction prognosis-related ceRNA network based on the key miRNAs
DElncRNA_DEmiRNA_interaction2.survivalRelated = DElncRNA_DEmiRNA_interaction2[DElncRNA_DEmiRNA_interaction2$DEmiRNA %in% cemiRNA_exp.L2H.survival.sig.res$miRNA,]
      DEmiRNA_DEmRNA_interaction.survialRelated = DEmiRNA_DEmRNA_interaction[DEmiRNA_DEmRNA_interaction$DEpre_miRNA %in% cemiRNA_exp.L2H.survival.sig.res$miRNA,]
      
      survivalRelated.edge =data.frame(
        "node1" = c(DElncRNA_DEmiRNA_interaction2.survivalRelated$DElncRNA,DEmiRNA_DEmRNA_interaction.survialRelated$DEmRNA),
        "node2" = c(DElncRNA_DEmiRNA_interaction2.survivalRelated$DEmiRNA,DEmiRNA_DEmRNA_interaction.survialRelated$DEpre_miRNA)
      )
      
      survivalRelated.node = data.frame("node" = c(unique(DElncRNA_DEmiRNA_interaction2.survivalRelated$DElncRNA),
                                                   unique(DEmiRNA_DEmRNA_interaction.survialRelated$DEmRNA),
                                                   unique(c(DEmiRNA_DEmRNA_interaction.survialRelated$DEpre_miRNA,
                                                            DElncRNA_DEmiRNA_interaction2.survivalRelated$DEmiRNA)) ),
                                        "type" = c(rep("lncRNA",length(unique(DElncRNA_DEmiRNA_interaction2.survivalRelated$DElncRNA)) ),
                                                   rep("mRNA",length(unique(DEmiRNA_DEmRNA_interaction.survialRelated$DEmRNA)) ),
                                                   rep("miRNA",length(unique(c(DEmiRNA_DEmRNA_interaction.survialRelated$DEpre_miRNA,
                                                                               DElncRNA_DEmiRNA_interaction2.survivalRelated$DEmiRNA))) ) )
      )
      
      survivalRelated.node$reg = apply(survivalRelated.node,1,function(x){
        if(x[1] %in% row.names(DEGs_cluster_HighTOLow)){
          ifelse(DEGs_cluster_HighTOLow[x[1],"logFC"] >0,"up","down" )
        }else if(x[1] %in% row.names(DE_miRNA_cluster_HighTOLow)){
          ifelse(DE_miRNA_cluster_HighTOLow[x[1],"logFC"] >0,"up","down" )
        }
      }) 
	  
	  
### PPI construction and analysis for the genes in the prognosis-related ceRNA newtork
##提取survival related cemRNA
if(T){
  survivalRelated.node.mRNAs = survivalRelated.node[which(survivalRelated.node$type == "mRNA"),]
  dim(survivalRelated.node.mRNAs)
  # [1] 493   3
  
}

##PPI interaction netwrok
if(T){
  PPI.edge = read.table("08.PPI_and_Hub_for_cemRNA/string_interactions_short.tsv",header = T,sep = "\t")
  length(unique(c(PPI.edge$node1,PPI.edge$node2)))
  # [1] 400
  table(unique(c(PPI.edge$node1,PPI.edge$node2)) %in% survivalRelated.node.mRNAs$node )
  # FALSE  TRUE 
  # 6   394
  
  ###filter out the nodes only included in prognosis-related ceRNA network
  PPI.edge.filter = PPI.edge[which( (PPI.edge$node1 %in% survivalRelated.node.mRNAs$node ) &
									  (PPI.edge$node2 %in% survivalRelated.node.mRNAs$node )  ),]
  length(unique(c(PPI.edge.filter$node1,PPI.edge.filter$node2)))
  # [1] 394
  table( unique(c(PPI.edge.filter$node1,PPI.edge.filter$node2)) %in% survivalRelated.node.mRNAs$node )
  
  ###calculate the degree of nodes
  View(table(c(PPI.edge.filter$node1,PPI.edge.filter$node2)))
  PPI.edge.filter.node = data.frame(table(c(PPI.edge.filter$node1,PPI.edge.filter$node2)))
  colnames(PPI.edge.filter.node) = c("node","degree")
  
}

##hub genes in PPI interaction network
if(T){
  PPI.hubgene.Top30 = read.table("08.PPI_and_Hub_for_cemRNA/07.survivalRelated.cemRNA.PPI.edge.txt_Degree_top30 default node.csv",
								 header = T,sep = ",")
  
  PPI.hubgene.Top30 = PPI.hubgene.Top30[,c(2,1)]
  PPI.hubgene.Top30$Regulation = apply(PPI.hubgene.Top30,1,function(x){
	survivalRelated.node$reg[which(survivalRelated.node$node == x[1])]
  })
  
  
}

##KM analysis of TOP30 hub genes
if(T){
  #DEGs_diff.df
  Top30.hubGene.exp = gene_exp.tumor.df.normalized[row.names(gene_exp.tumor.df.normalized) %in% PPI.hubgene.Top30$name,]
  
  Top30.hubGene.exp.t = cbind("patient_id" = substring(colnames(Top30.hubGene.exp),1,12),
							  "sample" = colnames(Top30.hubGene.exp),t(Top30.hubGene.exp))
  Top30.hubGene.exp.add.SurInfo = merge(survival.info,Top30.hubGene.exp.t,
										by.x = "PATIENT_ID",by.y="patient_id",all.y=T)
  
  colnames(Top30.hubGene.exp.add.SurInfo)[70:ncol(Top30.hubGene.exp.add.SurInfo)] = gsub("-","_",colnames(Top30.hubGene.exp.add.SurInfo)[70:ncol(Top30.hubGene.exp.add.SurInfo)])
  
  ##transformation
  Top30.hubGene.exp.add.SurInfo[,70:ncol(Top30.hubGene.exp.add.SurInfo)] = lapply(Top30.hubGene.exp.add.SurInfo[,70:ncol(Top30.hubGene.exp.add.SurInfo)],as.numeric)
  
  Top30.hubGene.exp.add.SurInfo = Top30.hubGene.exp.add.SurInfo[which(Top30.hubGene.exp.add.SurInfo$OS_MONTHS<120),]
  
  ####group
  Top30.hubGene.exp.L2H.survival = list()
  Top30.hubGene.exp.L2H.survival.sig = list()
  Top30.hubGene.exp.L2H.survival.sig.id = list()
  Top30.hubGene.exp.L2H.survival.sig.p = list()
  j=1
  for(i in 1:length(colnames(Top30.hubGene.exp.add.SurInfo)[70:ncol(Top30.hubGene.exp.add.SurInfo)])){ # i=1
	id = colnames(Top30.hubGene.exp.add.SurInfo)[70:ncol(Top30.hubGene.exp.add.SurInfo)][i]
	iexp = Top30.hubGene.exp.add.SurInfo[,id]
	
	imedian = median(iexp)
	idf.mat = Top30.hubGene.exp.add.SurInfo
	idf.mat$exp = ifelse(idf.mat[,id] <= imedian,0,1)
	
	library("survival")
	library("survminer")
	
	ifit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ exp, data=idf.mat)
	
	survdiff.h.l = survdiff(Surv(OS_MONTHS, OS_STATUS) ~ exp, data = idf.mat)
	pvalue.h.l = 1 - pchisq(survdiff.h.l$chisq, length(survdiff.h.l$n) -1)
	
	ip <- ggsurvplot(ifit, data = idf.mat,
					 pval = F,
					 pval.coord = c(0, 0.03),
					 size = 0.6,
					 palette=c("Set1"),  
					 legend = c(0.8,0.75),
					 legend.title= "",
					 title="", 
					 ylab="Overall survival (%)",
					 xlab = " Time (Months)", 
					 censor.shape = 124,censor.size = 2,conf.int = FALSE, 
					 font.legend = c(8,"bold", "black"),
					 font.main = c(10, "bold", "black"),
					 font.x = c(8, "bold", "black"),
					 font.y = c(8, "bold", "black"),
					 break.x.by = 40,
					 font.tickslab = c(8,"black")
	)
	
	ip$plot = ip$plot + ggplot2::annotate("text",x = 55, y = 1,color="black",size = 4,fontface="bold",
										  label = paste0(gsub("_","-",id)," ( p=",format(pvalue.h.l, scientific = T,digits = 3)," )"))
	
	Top30.hubGene.exp.L2H.survival[[i]] = ip
	
	if(pvalue.h.l < 0.05 ){
	  Top30.hubGene.exp.L2H.survival.sig[[j]] = ip
	  Top30.hubGene.exp.L2H.survival.sig.id[j] = id
	  Top30.hubGene.exp.L2H.survival.sig.p[j] = format(pvalue.h.l, scientific = T,digits = 3)
	  j = j + 1
	}
	
  }
  
  arrange_ggsurvplots(Top30.hubGene.exp.L2H.survival.sig,  ncol = 3, nrow = 4,newpage = FALSE)
  
  Top30.hubGene.exp.L2H.survival.sig.res = data.frame("Gene" = gsub("_","-",unlist(Top30.hubGene.exp.L2H.survival.sig.id)),
													  "pvalue" = unlist(Top30.hubGene.exp.L2H.survival.sig.p))

}

##Top module
if(T){
  TopModule.gene = read.table("08.PPI_and_Hub_for_cemRNA/14.survivalRelated.cemRNA.PPI.edge.txt.TopModule.geneList.txt",header = F)
  
  library(ggplot2)
  library(stringr)
  library(clusterProfiler)
  
  TopModule.gene.df <- bitr(TopModule.gene$V1, fromType="SYMBOL",
							toType="ENTREZID", 
							OrgDb = "org.Hs.eg.db")
  
  ##GO
  TopModule.gene.go <- enrichGO(gene = TopModule.gene.df$ENTREZID, 
								keyType = "ENTREZID",
								OrgDb = "org.Hs.eg.db",
								ont="BP",
								pAdjustMethod = "BH",
								pvalueCutoff  = 0.05,
								qvalueCutoff  = 0.05,
								readable      = TRUE)
  
  dotplot(TopModule.gene.go, showCategory = 20) 
  
  ##KEGG
  library(DO.db)
  TopModule.gene.ego <- enrichKEGG(
	gene = TopModule.gene.df$ENTREZID,
	keyType = "kegg",
	organism  = 'hsa',
	pvalueCutoff  = 0.05,
	pAdjustMethod  = "BH",
	qvalueCutoff  = 0.05
  )
  
  TopModule.gene.ego.2=DOSE::setReadable(TopModule.gene.ego, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  dotplot(TopModule.gene.ego, showCategory = 20) 
  
}


##KM analysis for the genes in top genes and top module
if(T){
##process
if(T){
  #DEGs_diff.df
  key_genes = unique(c(PPI.hubgene.Top30$name,TopModule.gene$V1))
  key_genes.exp = gene_exp.tumor.df.normalized[row.names(gene_exp.tumor.df.normalized) %in% key_genes,]
  
  key_genes.exp.t = cbind("patient_id" = substring(colnames(key_genes.exp),1,12),
						"sample" = colnames(key_genes.exp),t(key_genes.exp))
  key_genes.exp.add.SurInfo = merge(survival.info,key_genes.exp.t,
								  by.x = "PATIENT_ID",by.y="patient_id",all.y=T)
  
  colnames(key_genes.exp.add.SurInfo)[70:ncol(key_genes.exp.add.SurInfo)] = gsub("-","_",colnames(key_genes.exp.add.SurInfo)[70:ncol(key_genes.exp.add.SurInfo)])
  
  ##tramsformation
  key_genes.exp.add.SurInfo[,70:ncol(key_genes.exp.add.SurInfo)] = lapply(key_genes.exp.add.SurInfo[,70:ncol(key_genes.exp.add.SurInfo)],as.numeric)
  
  key_genes.exp.add.SurInfo = key_genes.exp.add.SurInfo[which(key_genes.exp.add.SurInfo$OS_MONTHS<120),]
  
  ####group
  key_genes.exp.L2H.survival = list()
  key_genes.exp.L2H.survival.sig = list()
  key_genes.exp.L2H.survival.sig.id = list()
  key_genes.exp.L2H.survival.sig.p = list()
  j=1
  for(i in 1:length(colnames(key_genes.exp.add.SurInfo)[70:ncol(key_genes.exp.add.SurInfo)])){ # i=1
	id = colnames(key_genes.exp.add.SurInfo)[70:ncol(key_genes.exp.add.SurInfo)][i]
	iexp = key_genes.exp.add.SurInfo[,id]
	
	imedian = median(iexp)
	idf.mat = key_genes.exp.add.SurInfo
	idf.mat$exp = ifelse(idf.mat[,id] <= imedian,0,1)
	
	library("survival")
	library("survminer")
	
	ifit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ exp, data=idf.mat)
	
	survdiff.h.l = survdiff(Surv(OS_MONTHS, OS_STATUS) ~ exp, data = idf.mat)
	pvalue.h.l = 1 - pchisq(survdiff.h.l$chisq, length(survdiff.h.l$n) -1)
	
	ip <- ggsurvplot(ifit, data = idf.mat,
					 pval = F,
					 pval.coord = c(0, 0.03),
					 size = 0.2,
					 palette=c("cornflowerblue","brown3"),  
					 legend = c(0.8,0.75),
					 legend.labs=c(paste0("Low (n=",nrow(idf.mat[idf.mat$exp == 0,]),")"),
								   paste0("High (n=",nrow(idf.mat[idf.mat$exp == 1,]),")") ), 
					 legend.title= "",
					 title="", 
					 ylab="Overall survival (%)",
					 xlab = " Time (Months)", 
					 censor.shape = 124,censor.size = 1,conf.int = FALSE, 
					 font.legend = c(3,"bold", "black"),
					 font.main = c(4, "bold", "black"),
					 font.x = c(3, "bold", "black"),
					 font.y = c(3, "bold", "black"),
					 break.x.by = 40,
					 font.tickslab = c(3,"black")
	)
	
	ip$plot = ip$plot + ggplot2::annotate("text",x = 55, y = 1,color="black",size = 1,fontface="bold",
										  label = paste0(gsub("_","-",id)," ( p=",format(pvalue.h.l, scientific = T,digits = 3)," )"))+
	  theme(
		axis.ticks.x=element_line(size=0.15),
		axis.ticks.y=element_line(size=0.15),
		axis.line.x=element_line(size=0.2),
		axis.line.y=element_line(size=0.2)
	  )
	
	key_genes.exp.L2H.survival[[i]] = ip
	
	if(pvalue.h.l < 0.05 ){
	  key_genes.exp.L2H.survival.sig[[j]] = ip
	  key_genes.exp.L2H.survival.sig.id[j] = id
	  key_genes.exp.L2H.survival.sig.p[j] = format(pvalue.h.l, scientific = T,digits = 3)
	  j = j + 1
	}
	
  }
  
  arrange_ggsurvplots(key_genes.exp.L2H.survival.sig,  ncol = 3, nrow = 3,newpage = FALSE)
  
  key_genes.exp.L2H.survival.sig.res = data.frame("gene" = gsub("_","-",unlist(key_genes.exp.L2H.survival.sig.id)),
												"pvalue" = unlist(key_genes.exp.L2H.survival.sig.p))
  
}


}

##KM analysis for target kncRNAs of key miRNA
if(T){
##process
if(T){
  target_lncRNA.exp = gene_exp.tumor.df.normalized[row.names(gene_exp.tumor.df.normalized) %in% 
													 survivalRelated.node$node[survivalRelated.node$type=="lncRNA"],]
  
  target_lncRNA.exp.t = cbind("patient_id" = substring(colnames(target_lncRNA.exp),1,12),
						  "sample" = colnames(target_lncRNA.exp),t(target_lncRNA.exp))
  target_lncRNA.exp.add.SurInfo = merge(survival.info,target_lncRNA.exp.t,
									by.x = "PATIENT_ID",by.y="patient_id",all.y=T)
  
  colnames(target_lncRNA.exp.add.SurInfo)[70:ncol(target_lncRNA.exp.add.SurInfo)] = gsub("-","_",colnames(target_lncRNA.exp.add.SurInfo)[70:ncol(target_lncRNA.exp.add.SurInfo)])
  
  ##tramsformation
  target_lncRNA.exp.add.SurInfo[,70:ncol(target_lncRNA.exp.add.SurInfo)] = lapply(target_lncRNA.exp.add.SurInfo[,70:ncol(target_lncRNA.exp.add.SurInfo)],as.numeric)
  
  target_lncRNA.exp.add.SurInfo = target_lncRNA.exp.add.SurInfo[which(target_lncRNA.exp.add.SurInfo$OS_MONTHS<120),]
  
  ####group
  target_lncRNA.exp.L2H.survival = list()
  target_lncRNA.exp.L2H.survival.sig = list()
  target_lncRNA.exp.L2H.survival.sig.id = list()
  target_lncRNA.exp.L2H.survival.sig.p = list()
  j=1
  for(i in 1:length(colnames(target_lncRNA.exp.add.SurInfo)[70:ncol(target_lncRNA.exp.add.SurInfo)])){ # i=35
	id = colnames(target_lncRNA.exp.add.SurInfo)[70:ncol(target_lncRNA.exp.add.SurInfo)][i]
	iexp = target_lncRNA.exp.add.SurInfo[,id]
	
	imedian = median(iexp)
	idf.mat = target_lncRNA.exp.add.SurInfo
	idf.mat$exp = ifelse(idf.mat[,id] <= imedian,0,1)
	
	library("survival")
	library("survminer")
	
	ifit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ exp, data=idf.mat)
	
	survdiff.h.l = survdiff(Surv(OS_MONTHS, OS_STATUS) ~ exp, data = idf.mat)
	pvalue.h.l = 1 - pchisq(survdiff.h.l$chisq, length(survdiff.h.l$n) -1)
	
	ip <- ggsurvplot(ifit, data = idf.mat,
					 pval = F,
					 pval.coord = c(0, 0.03),
					 size = 0.2,
					 palette=c("cornflowerblue","brown3"), 
					 legend = c(0.8,0.75),
					 legend.labs=c(paste0("Low (n=",nrow(idf.mat[idf.mat$exp == 0,]),")"),
								   paste0("High (n=",nrow(idf.mat[idf.mat$exp == 1,]),")") ), 
					 legend.title= "",
					 title="", 
					 ylab="Overall survival (%)",
					 xlab = " Time (Months)", 
					 censor.shape = 124,censor.size = 1,conf.int = FALSE, 
					 font.legend = c(3,"bold", "black"),
					 font.main = c(4, "bold", "black"),
					 font.x = c(3, "bold", "black"),
					 font.y = c(3, "bold", "black"),
					 break.x.by = 40,
					 font.tickslab = c(3,"black")
	)
	
	
	ip$plot = ip$plot + ggplot2::annotate("text",x = 55, y = 1,color="black",size = 1,fontface="bold",
										  label = paste0(gsub("_","-",id)," ( p=",format(pvalue.h.l, scientific = T,digits = 3)," )"))+
	  theme(
		axis.ticks.x=element_line(size=0.15),
		axis.ticks.y=element_line(size=0.15),
		axis.line.x=element_line(size=0.2),
		axis.line.y=element_line(size=0.2)
	  )
	
	
	
	
	
	
	
	target_lncRNA.exp.L2H.survival[[i]] = ip
	
	if(pvalue.h.l < 0.05 ){
	  target_lncRNA.exp.L2H.survival.sig[[j]] = ip
	  target_lncRNA.exp.L2H.survival.sig.id[j] = id
	  target_lncRNA.exp.L2H.survival.sig.p[j] = format(pvalue.h.l, scientific = T,digits = 3)
	  j = j + 1
	}
	
  }
  
  arrange_ggsurvplots(target_lncRNA.exp.L2H.survival.sig,  ncol = 3, nrow = 3,newpage = FALSE)
  
  target_lncRNA.exp.L2H.survival.sig.res = data.frame("gene" = gsub("_","-",unlist(target_lncRNA.exp.L2H.survival.sig.id)),
												  "pvalue" = unlist(target_lncRNA.exp.L2H.survival.sig.p))
  
}

### co-expression analysis for the key miRNAs' targeting mRNAs and targeting lncRNAs in the prognosis-related ceRNA network
#cemRNA/celncRNA
survivalRelated_cemRNAs.celncRNA2 = gene_exp.tumor.df.normalized[which( row.names(gene_exp.tumor.df.normalized) %in% 
																		 c(key_genes.exp.L2H.survival.sig.res$gene, 
																		   target_lncRNA.exp.L2H.survival.sig.res$gene )),]
survivalRelated_cemRNAs.celncRNA.t2 = t(survivalRelated_cemRNAs.celncRNA2)
colnames(survivalRelated_cemRNAs.celncRNA.t2) = gsub("-","_",colnames(survivalRelated_cemRNAs.celncRNA.t2))
survivalRelated_cemRNAs.celncRNA.t2 = data.frame(survivalRelated_cemRNAs.celncRNA.t2)

##PCC
kmiRNA_corr2 <- round(cor(survivalRelated_cemRNAs.celncRNA.t2), 2)
library(ggcorrplot)
kmiRNA_p.mat2 <- format(cor_pmat(survivalRelated_cemRNAs.celncRNA.t2), scientific = T,digits = 3)

##transformation
kmiRNA_corr3 = kmiRNA_corr2[gsub("-","_",row.names(kmiRNA_corr2)) %in% key_genes.exp.L2H.survival.sig.res$gene ,
							gsub("_","-",colnames(kmiRNA_corr2)) %in% target_lncRNA.exp.L2H.survival.sig.res$gene]

kmiRNA_p.mat3 = kmiRNA_p.mat2[gsub("-","_",row.names(kmiRNA_p.mat2)) %in% key_genes.exp.L2H.survival.sig.res$gene ,
							 gsub("_","-",colnames(kmiRNA_p.mat2)) %in% target_lncRNA.exp.L2H.survival.sig.res$gene]

##reshape
library(reshape2)
kmiRNA_corr4 = cbind("gene" = row.names(kmiRNA_corr3),kmiRNA_corr3)
kmiRNA_corr4 = melt(kmiRNA_corr3,value.name = "gene")
View(kmiRNA_corr4)

kmiRNA_p.mat4 = cbind("gene" = row.names(kmiRNA_p.mat3),kmiRNA_p.mat3)
kmiRNA_p.mat4 = melt(kmiRNA_p.mat4, id.vars = "gene")
kmiRNA_p.mat4 = kmiRNA_p.mat4[(nrow(kmiRNA_p.mat3)+1):nrow(kmiRNA_p.mat4),]
View(kmiRNA_p.mat4)

##filter
candidate_pair = merge(kmiRNA_corr4, kmiRNA_p.mat4, by=c("Var1","Var2"))
colnames(candidate_pair) = c("gene", "lncRNA", "cor", "pvalue")
candidate_pair2 = candidate_pair[which(abs(candidate_pair$cor >= 0.4) & candidate_pair$pvalue < 0.05),]
View(candidate_pair2)

##filter
table(survivalRelated.edge$node2)
View(survivalRelated.edge)

m = 1
n = 1
node1 = list()
node2 = list()
for(i in 1:nrow(candidate_pair2)){ # i=1   
  igene = as.character( candidate_pair2[i,1] )
  ilncRNA = as.character( candidate_pair2[i,2] )
  for(j in 1:nrow(cemiRNA_exp.L2H.survival.sig.res)){ # j=1
	jmiRNA = cemiRNA_exp.L2H.survival.sig.res[j,1]
	if(nrow(survivalRelated.edge[which(survivalRelated.edge$node1 == igene & survivalRelated.edge$node2 == jmiRNA),]) == 1 &
	   nrow(survivalRelated.edge[which(survivalRelated.edge$node1 == ilncRNA & survivalRelated.edge$node2 == jmiRNA),]) == 1 ){
	  node1[m] = jmiRNA
	  node1[m+1] = jmiRNA
	  node2[n] = igene
	  node2[n+1] = ilncRNA
	  m = m + 2 
	  n = n + 2
	}
			  
  }
}

survivalRelated.edge.candidatePair = data.frame("node1" = unlist(node2), "node2" = unlist(node1))
survivalRelated.edge.candidatePair2 = survivalRelated.edge.candidatePair[!duplicated(survivalRelated.edge.candidatePair),]
View(survivalRelated.edge.candidatePair2)

##plot
if(T){
  
  survival_related.ceRNA.cor.plot2 = list()
  k=1
  
  for(i in 1:length(unique(survivalRelated.edge.candidatePair2$node2))){ 
	ikeymiRNA = unique(survivalRelated.edge.candidatePair2$node2)[i]
	itargets = survivalRelated.edge.candidatePair2$node1[survivalRelated.edge.candidatePair2$node2==ikeymiRNA]
	ipairs = candidate_pair2[which(candidate_pair2$gene %in% itargets & candidate_pair2$lncRNA %in% itargets),]
	
	for(j in 1:nrow(ipairs)){ 
	  jlncRNA = as.character( ipairs$lncRNA[j] )
	  jmRNA = as.character( ipairs$gene[j] )
	  jp <- ggscatter(survivalRelated_cemRNAs.celncRNA.t, x = jmRNA, y = jlncRNA,
					  title = paste0("Target of ",ikeymiRNA),
					  color = "#426671", size =0.05, 
					  add = "reg.line",  
					  add.params = list(color = "#764C29", fill = "#E7E1D7",size=0.2), 
					  conf.int = TRUE, 
					  xlab = gsub("_","-",jmRNA), ylab = gsub("_","-",jlncRNA) ) +
		stat_cor(method = "pearson", label.x = 3, label.sep = "\n", p.accuracy = 0.001, r.accuracy = 0.001, size = 1.5)+
		theme(title=element_text(color="Black", face="bold",size = 4),
			  axis.text.x = element_text(color="Black", size=3,hjust=0.5,vjust = 0),
			  axis.text.y = element_text(color ="Black", size = 3, angle = 0, hjust = 1, vjust = 0.5),
			  axis.title.x = element_text(face="bold", color ="Black", size = 4, angle = 0, hjust = .5, vjust = 0.5),
			  axis.title.y = element_text(face="bold",color ="Black", size = 4, angle = 90, hjust = 0.5, vjust = 1.5),
			  legend.position="top",
			  legend.text = element_text(color ="Black", size = 4, angle = 0, hjust = 0, vjust = 0.5),
			  legend.title = element_text(face="bold",color ="Black", size = 4, angle = 0, hjust = .5, vjust = .5),
			  strip.text = element_text(face="bold", color ="Black", size = 4),
			  axis.ticks.x=element_line(size=0.15),
			  axis.ticks.y=element_line(size=0.15),
			  axis.line.x=element_line(size=0.2),
			  axis.line.y=element_line(size=0.2)
			  )
	  
	  survival_related.ceRNA.cor.plot2[[k]] = jp
	  k = k + 1
	  
	}
	
  }
  
   ggarrange(plotlist = survival_related.ceRNA.cor.plot2, ncol = 5, nrow = 5)
  
}

survivalRelated.node.candidateNode = survivalRelated.node[survivalRelated.node$node %in% 
															unique(c(survivalRelated.edge.candidatePair2$node1,
																	 survivalRelated.edge.candidatePair2$node2)),]
